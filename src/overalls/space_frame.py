import math
from copy import copy
from itertools import combinations
from itertools import cycle

import Rhino.Geometry as rg
from compas.datastructures import Mesh
from compas.datastructures import Network
from compas.geometry import Point
from compas.geometry import Vector
from compas.utilities import geometric_key
from compas_rhino.geometry import RhinoLine
from overalls.utils import flatten
from overalls.utils import is_in_domain
from overalls.utils import to_list


class SpaceFrameMixin(object):
    def connected_edges(self, key):
        u = key

        if isinstance(self, Network):
            nbors = self.neighbors(u)
        else:
            nbors = self.vertex_neighbors(u)

        for v in nbors:
            yield (u, v)

    def ok_degree(self, keys, max_degree):
        if isinstance(self, Network):
            degree_func = self.degree
        else:
            degree_func = self.vertex_degree

        keys = to_list(keys)
        for key in keys:
            if degree_func(key) > max_degree:
                return False
        return True

    def ok_edges_angles(self, keys, min_angle, **kwargs):
        vkeys = to_list(keys)
        for vkey in vkeys:
            if not self.ok_edge_angles(vkey, min_angle, **kwargs):
                return False
        return True

    def ok_edge_angles(self, key, min_angle, additional_edge=None):
        edges = list(self.connected_edges(key))

        if additional_edge:
            edges.append(additional_edge)

        vectors = []
        for u, v in edges:
            if u != key:
                u, v = v, u
            vec = rg.Vector3d(*self.edge_vector(u, v))
            vec.Unitize()
            vectors.append(vec)

        to_compare = combinations(range(len(vectors)), 2)
        for v, v1 in list(to_compare):
            v, v1 = vectors[v], vectors[v1]
            if rg.Vector3d.VectorAngle(v, v1) < min_angle:
                return False

        return True

    def ok_edge_length(self, ekeys, length_domain):
        for u, v in ekeys:
            if not is_in_domain(self.edge_length(u, v), length_domain):
                return False
        return True

    def edge_to_rgline(self, u, v):
        pt_a, pt_b = self.edge_coordinates(u, v)
        pt_a = rg.Point3d(*pt_a)
        pt_b = rg.Point3d(*pt_b)
        return rg.Line(pt_a, pt_b)

    def to_rglines(self):
        lines = []

        for u, v in self.edges():
            lines.append(self.edge_to_rgline(u, v))

        return lines


class SpaceFrameNetwork(SpaceFrameMixin, Network):

    FROM_MESH = 0
    FROM_CONN = 1

    def __init__(self):
        super(SpaceFrameNetwork, self).__init__()
        self.dist_dict = {}

    def ensure_dist_dict(self):
        if set(self.nodes()) != set(self.dist_dict.keys()):
            self.build_dist_dict()

    def build_dist_dict(self):
        for nkey in self.nodes():
            self.dist_dict[nkey] = {}

        key_combinations = combinations(self.nodes(), 2)

        for u, v in key_combinations:
            x, y, z = self.node_coordinates(u)
            x1, y1, z1 = self.node_coordinates(v)

            dist = math.sqrt((x - x1) ** 2 + (y - y1) ** 2 + (z - z1) ** 2)

            self.dist_dict[u][v] = dist
            self.dist_dict[v][u] = dist

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

    def ok_conn(self, u, v, max_degree, min_angle):
        # print("Testing connection {}-{}".format(u, v))
        if self.has_edge(u, v, directed=False):
            print("Not ok_conn because has_edge")
            return False

        if max_degree:
            if not self.ok_degree([u, v], max_degree):
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
        self, start_nkey, end_nkeys, max_degree=None, min_angle=None, max_n_conn=None,
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
            if self.ok_conn(u, v, max_degree - 1, min_angle):
                self.add_edge(
                    u,
                    v,
                    created_from=self.FROM_CONN,
                    # parent_fkey=(start_fkey, end_fkey),
                )
                n_conns += 1

    @classmethod
    def from_space_frame_mesh(cls, mesh):
        network = cls()

        for vkey in mesh.vertices():
            x, y, z = mesh.vertex_coordinates(vkey)
            # data = mesh.vertex_attributes(vkey)
            # data.update({"created_from": cls.FROM_MESH})
            network.add_node(key=vkey, x=x, y=y, z=z)

        for u, v in mesh.edges():
            # data = mesh.edge_attributes((u, v))
            # data.update({"created_from": cls.FROM_MESH})
            network.add_edge(u, v)

        return network

    @classmethod
    def from_rglines(cls, lines, attr_dicts=None):
        def add_node_if_needed(xyz, node_gkeys):
            gkey = geometric_key(xyz)
            nkey = node_gkeys.get(gkey)

            if nkey is not None:
                return nkey

            x, y, z = xyz
            nkey = cls.add_node(x=x, y=y, z=z, created_from=cls.FROM_LINES)
            node_gkeys[gkey] = nkey

            return nkey

        network = cls()

        node_gkeys = {}

        for i, rgline in enumerate(lines):
            line = RhinoLine.from_geometry(rgline).to_compas()
            data = attr_dicts[i] if i < len(attr_dicts) else {}

            data.update({"created_from": cls.FROM_MESH})

            u = add_node_if_needed(list(line.start), node_gkeys)
            v = add_node_if_needed(list(line.end), node_gkeys)

            network.add_edge(u, v, attr_dict=data)

        return network


class SpaceFrameMesh(SpaceFrameMixin, Mesh):
    FROM_SUBD_TRI = 0
    FROM_SUBD_FRAME = 1

    def __init__(self):
        super(SpaceFrameMesh, self).__init__()
        self.update_default_face_attributes(parent_fkey=None, created_from=None)
        self.update_default_vertex_attributes(parent_fkey=None, created_from=None)
        self.update_default_edge_attributes(parent_fkey=None, created_from=None)

    def ok_subd(
        self,
        old_vkeys,
        new_vkeys,
        max_degree=None,
        min_angle=None,
        edge_length_domain=None,
        **kwargs
    ):
        new_vkeys = to_list(new_vkeys)

        if max_degree:
            if not self.ok_degree(old_vkeys, max_degree):
                return False

        if min_angle:
            if not self.ok_edges_angles(flatten([old_vkeys + new_vkeys]), min_angle):
                return False

        if edge_length_domain:
            ekeys = [self.connected_edges(vkey) for vkey in new_vkeys]
            ekeys = flatten(ekeys)
            uv_sets = set(frozenset((u, v)) for u, v in ekeys)  # unique edges
            if not self.ok_edge_length(uv_sets, edge_length_domain):
                return False

        return True

    def subdiv_faces(
        self, fkeys, n_iters, scheme=[0], rel_pos=0.5, rel_pos_z=0.0, **kwargs
    ):
        schemes = [
            self.face_subdiv_retriangulate,
            self.face_subdiv_mid_cent,
            self.face_subdiv_mids,
            self.face_subdiv_frame,
        ]

        repeating_n_iters = cycle(n_iters)
        repeating_schemes = cycle(scheme)
        repeating_rel_pos = cycle(to_list(rel_pos))
        repeating_rel_pos_z = cycle(to_list(rel_pos_z))

        new_subd_fkeys = []
        for fkey in fkeys:
            n_iters = next(repeating_n_iters)
            kwargs.update({"rel_pos": next(repeating_rel_pos)})
            kwargs.update({"rel_pos_z": next(repeating_rel_pos_z)})

            subdiv_func = schemes[next(repeating_schemes)]

            i = 0
            next_to_subd = [fkey]
            while i < n_iters:

                to_subd = next_to_subd
                next_to_subd = []

                for parent_fkey in to_subd:
                    parent_face_verts = self.face_vertices(parent_fkey)
                    parent_attrs = self.face_attributes(parent_fkey)

                    new_vkeys, new_fkeys = subdiv_func(parent_fkey, **kwargs)

                    if not self.ok_subd(parent_face_verts, new_vkeys, **kwargs):
                        self.undo_face_subd(
                            new_fkeys, parent_fkey, parent_face_verts, parent_attrs
                        )
                        continue

                    next_to_subd += new_fkeys
                    new_subd_fkeys += new_fkeys

                i += 1

        self.cull_vertices()

        return [fkey for fkey in new_fkeys if self.has_face(fkey)]

    def face_subdiv_retriangulate(self, fkey, rel_pos=None, rel_pos_z=0.0, **kwattr):
        xyz = Point(*self.face_center(fkey))

        vecs = []
        if rel_pos or rel_pos_z:
            for v_xyz in self.face_coordinates(fkey):
                vecs.append(xyz - Point(*v_xyz))

        if rel_pos:
            rel_pos = cycle(to_list(rel_pos))
            for v in vecs:
                factor = next(rel_pos)
                xyz -= v * factor

        if rel_pos_z:
            # set up max dist
            dists = [v.length for v in vecs]
            dists.sort()

            normal = Vector(*self.face_normal(fkey, unitized=True))
            z_vec = normal * dists[0] * rel_pos_z
            xyz += z_vec

        vkey, fkeys = self.insert_vertex(fkey, xyz=list(xyz), return_fkeys=True)

        ekeys = [(vkey, v) for v in self.vertex_neighbors(vkey)]  # noqa: F841

        data = kwattr
        data.update({"parent_fkey": fkey})

        # self.components_attributes(fkeys, [vkey], ekeys, data)

        return vkey, fkeys

    def face_subdiv_mid_cent(self, fkey, rel_pos=None, rel_pos_z=0.0, **kwattr):
        x, y, z = Point(*self.face_center(fkey))

        new_vkey_center = self.add_vertex(x=x, y=y, z=z)

        new_vkeys = []
        for u in self.face_vertices(fkey):
            v = self.face_vertex_descendant(fkey, u)
            x, y, z = self.edge_midpoint(u, v)
            new_vkeys.append(self.add_vertex(x=x, y=y, z=z))

        new_fkeys = []
        for corner, mid in zip(self.face_vertices(fkey), new_vkeys):
            vkeys = [corner, mid, new_vkey_center]
            new_fkeys.append(self.add_face(vkeys))

            next_corner = self.face_vertex_descendant(fkey, corner)
            vkeys = [mid, next_corner, new_vkey_center]
            new_fkeys.append(self.add_face(vkeys))

        self.delete_face(fkey)

        new_vkeys.append(new_vkey_center)

        return new_vkeys, new_fkeys

    def face_subdiv_mids(self, fkey, rel_pos=None, rel_poz_z=0.0, **kwattr):
        new_vkeys = []
        for u in self.face_vertices(fkey):
            v = self.face_vertex_descendant(fkey, u)
            x, y, z = self.edge_midpoint(u, v)
            new_vkeys.append(self.add_vertex(x=x, y=y, z=z))

        new_fkeys = []
        for i, corner in enumerate(self.face_vertices(fkey)):
            prev_mid = new_vkeys[(i - 1) % len(new_vkeys)]
            mid = new_vkeys[i]
            vkeys = [prev_mid, corner, mid]
            new_fkeys.append(self.add_face(vkeys))

        # middle face
        new_fkeys.append(self.add_face(new_vkeys))

        self.delete_face(fkey)

        return new_vkeys, new_fkeys

    def face_subdiv_frame(self, fkey, move_z=False, rel_dist=0.5, **kwattr):
        if rel_dist == 1:
            return self.face_subdiv_retriangulate(fkey, move_z=move_z ** kwattr)
        if rel_dist == 0:
            raise Exception("rel_dist can't be 0")

        face_center_pt = Point(*self.face_center(fkey))

        new_vkeys = []
        for x, y, z in self.face_coordinates(fkey):
            pt = Point(x, y, z)

            v = face_center_pt - pt
            pt += v * rel_dist
            new_vkeys.append(self.add_vertex(x=pt.x, y=pt.y, z=pt.z))

        face_verts = self.face_vertices(fkey)
        new_fkeys = []
        for j in range(len(face_verts)):
            vkeys = []
            vkeys.append(face_verts[j])
            vkeys.append(face_verts[(j + 1) % len(face_verts)])
            vkeys.append(new_vkeys[(j + 1) % len(new_vkeys)])
            vkeys.append(new_vkeys[j])

            new_fkeys.append(self.add_face(vkeys))

        # add new center face
        new_fkeys.append(self.add_face(new_vkeys))

        self.delete_face(fkey)

        return new_vkeys, new_fkeys

    def undo_face_subd(self, new_fkeys, old_fkey, old_face_verts, old_attrs):
        for fkey in new_fkeys:
            try:
                self.delete_face(fkey)
            except KeyError:
                pass

        self.add_face(old_face_verts, fkey=old_fkey, attr_dict=old_attrs)
