from copy import copy

from compas.datastructures import Network
from compas.geometry import transform_points
from compas.utilities import geometric_key
from compas_rhino.geometry import RhinoLine
from overalls.space_frame._mixin import SpaceFrameMixin
from overalls.utils import is_in_domain
from overalls.utils import to_list


class SpaceFrameNetwork(Network, SpaceFrameMixin):

    FROM_MESH = 0
    FROM_CONN = 1
    FROM_LINES = 2

    def __init__(self):
        super(SpaceFrameNetwork, self).__init__()
        self.update_default_node_attributes(
            parent_fkey=None, created_from=None, part_id=None
        )
        self.update_default_edge_attributes(
            parent_fkey=None,
            created_from=None,
            structural_data=None,
            part_id=None,
            ignore_edge=False,
        )

    def nv_degree(self, key, *args, **kwargs):
        return self.degree(key, *args, **kwargs)

    def nvs(self, *args, **kwargs):
        return self.nodes(*args, **kwargs)

    def nv_neighbors(self, key, *args, **kwargs):
        return self.neighbors(key, *args, **kwargs)

    def nv_coordinates(self, key, *args, **kwargs):
        return self.node_coordinates(key, *args, **kwargs)

    def nv_attribute(self, key, *args, **kwargs):
        return self.node_attribute(key, *args, **kwargs)

    def faces(self, *args, **kwargs):
        return iter(())

    def transform_network(self, transformation):
        nodes = [self.node_coordinates(key) for key in self.nodes()]
        xyz = transform_points(nodes, transformation)
        for index, (key, attr) in enumerate(self.nodes(True)):
            attr["x"] = xyz[index][0]
            attr["y"] = xyz[index][1]
            attr["z"] = xyz[index][2]

        self.build_dist_dict()

    def find_closest_node(
        self,
        nkey_origin,
        nkeys_search,
        length_domain=None,
        prefer_distant=False,
        **kwargs
    ):
        self.ensure_dist_dict()
        dists = copy(self.dist_dict[nkey_origin])

        for nkey in self.dist_dict.keys():
            if nkey == nkey_origin:
                continue
            if nkey not in nkeys_search:
                # print("Not in search")
                dists.pop(nkey)
                continue
            if nkey in self.neighbors(nkey_origin):
                # print("in nbors")
                dists.pop(nkey)
                continue
            if not is_in_domain(dists[nkey], length_domain):
                # print("not in domain")
                dists.pop(nkey)
                continue

        if len(dists) < 1:  # No matches
            return None

        dists = dists.items()
        dists.sort(key=lambda x: x[1], reverse=prefer_distant)
        keys, _ = zip(*dists)
        return keys[0]

    def ok_conn(self, u, v, degree_domain, min_angle):
        # print("Testing connection {}-{}".format(u, v))
        if self.has_edge(u, v, directed=False):
            # print("Not ok_conn because has_edge")
            return False

        if degree_domain:
            if not self.ok_degrees([u, v], degree_domain):
                # print("Not ok_conn because vertex_degree: {}-{}".format(u, v))
                return False

        if min_angle:
            if not self.ok_edges_angles([u, v], min_angle, additional_edge=(u, v)):
                # print("Not ok_conn because vertex_angles: {}-{}".format(u, v))
                return False

        # print("{}-{} passed ok_conn".format(u, v))
        return True

    def find_closest_nhood(
        self,
        nkey_origin,
        nkeys_search,
        n_results=3,
        length_domain=None,
        include_start_node=True,
        **kwargs
    ):
        closest_nkey = self.find_closest_node(
            nkey_origin, nkeys_search, length_domain=length_domain, **kwargs
        )
        nbors = []

        if include_start_node:
            nbors.append(closest_nkey)

        rings = 1
        nhood = []
        while nhood < n_results:
            nhood = self.neighborhood(closest_nkey, rings=rings)
            rings += 1

        nhood = nhood[:n_results]

        for nkey in nhood:
            if is_in_domain(self.dist_dict[nkey_origin][nkey], length_domain):
                nbors.append(nkey)

        return nbors

    def connect_nodes(
        self,
        start_nkey,
        end_nkeys,
        degree_domain=None,
        min_angle=None,
        max_n_conn=None,
    ):
        """Create a line node to node."""
        u = start_nkey
        end_nkeys = to_list(end_nkeys)

        n_conns = len(list(self.nodes_where({"created_from": self.FROM_CONN})))
        for v in end_nkeys:
            if v is None:
                continue

            if max_n_conn:
                if n_conns > max_n_conn:
                    break

            # max_degree-1 because we will add one
            min_degree, max_degree = degree_domain
            min_degree = min_degree + 1 if min_degree else None
            max_degree = max_degree - 1 if max_degree else None

            if self.ok_conn(u, v, (min_degree, max_degree), min_angle):
                self.add_edge(
                    u,
                    v,
                    created_from=self.FROM_CONN,
                    # parent_fkey=(start_fkey, end_fkey),
                )
                n_conns += 1

    @classmethod
    def from_space_frame_mesh(cls, mesh, scale="mm", check_ignore_edge_attr=True):
        # TODO: Bugfix, creates super long lines.
        if scale != "mm":
            raise NotImplementedError("Scaling needs fixing.")
        if scale not in ("m", "mm"):
            raise Exception("Scale needs to be mm or m.")

        network = cls()

        if scale == "m":
            T = mesh.get_scale_cgxform(1000)
            mesh.transform(T)

        for vkey in mesh.vertices():
            data = mesh.vertex_attributes(vkey)
            data.update({"created_from": cls.FROM_MESH})
            network.add_node(key=vkey, **data)

        for u, v in mesh.edges():
            data = mesh.edge_attributes((u, v))
            if check_ignore_edge_attr and data.get("ignore_edge"):
                continue
            data = {"created_from": cls.FROM_MESH}
            network.add_edge(u, v, **data)

        return network

    @staticmethod
    def _add_node_if_needed(network, xyz, node_gkeys, data={}):
        gkey = geometric_key(xyz)
        nkey = node_gkeys.get(gkey)

        if nkey is not None:
            return nkey

        x, y, z = xyz
        nkey = network.add_node(x=x, y=y, z=z, attr_dict=data)
        node_gkeys[gkey] = nkey

        return nkey

    @classmethod
    def from_rglines(cls, lines, scale="mm", attr_dicts=[]):
        if scale != "mm":
            raise NotImplementedError("Scaling needs fixing.")
        if scale not in ("m", "mm"):
            raise Exception("Scale needs to be mm or m.")

        network = cls()

        node_gkeys = {}

        for i, rgline in enumerate(lines):
            line = RhinoLine.from_geometry(rgline).to_compas()
            data = attr_dicts[i] if i < len(attr_dicts) else {}

            data.update({"created_from": cls.FROM_LINES})

            u = cls._add_node_if_needed(network, list(line.start), node_gkeys, data)
            v = cls._add_node_if_needed(network, list(line.end), node_gkeys, data)

            network.add_edge(u, v, attr_dict=data)

        if scale == "m":
            T = network.get_scale_cgxform(1000)
            network.transform_network(T)

        return network
