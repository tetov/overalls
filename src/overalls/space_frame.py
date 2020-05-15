import math
from copy import copy
from itertools import combinations
from itertools import cycle

import Rhino.Geometry as rg
from compas.datastructures import Mesh
from compas.datastructures import Network
from compas.geometry import Point
from compas.geometry import Scale
from compas.geometry import Translation
from compas.geometry import Vector
from compas.geometry import transform_points
from compas.utilities import geometric_key
from compas_ghpython.artists import MeshArtist
from compas_rhino.geometry import RhinoLine
from ghpythonlib.treehelpers import list_to_tree
from overalls.utils import flatten
from overalls.utils import is_in_domain
from overalls.utils import to_list


class SpaceFrameMixin(object):
    def get_bbox_origin(self):
        if hasattr(self, "nodes"):
            pts = [self.node_coordinates(key) for key in self.nodes()]
        else:
            pts = [self.vertex_coordinates(key) for key in self.vertices()]

        rg_pts = [rg.Point3d(*pt) for pt in pts]

        bbox = rg.BoundingBox(rg_pts)
        return bbox.Corner(True, True, True)

    def get_scale_rgxform(self, factor):
        scale_origin = self.get_bbox_origin()
        return rg.Transform.Scale(scale_origin, factor)

    def get_scale_cgxform(self, factor):
        scale_origin = self.get_bbox_origin()
        scale_origin = [scale_origin.X, scale_origin.Y, scale_origin.Z]

        return Translation(scale_origin) * Scale([factor] * 3)

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

    def to_rglines(self, ekeys=None):
        if ekeys is None:
            ekeys = self.edges()

        lines = []

        for u, v in ekeys:
            lines.append(self.edge_to_rgline(u, v))

        return list_to_tree(lines)

    def to_rglines_meters(self, ekeys=None):
        # TODO: Bugfix, creates super long lines.
        raise NotImplementedError("Scaling needs fixing.")
        T = self.get_scale_rgxform(0.001)
        lines = self.to_rglines(ekeys=ekeys)
        for line in lines:
            line.Transform(T)

        return lines


class SpaceFrameNetwork(SpaceFrameMixin, Network):

    FROM_MESH = 0
    FROM_CONN = 1
    FROM_LINES = 2

    def __init__(self):
        super(SpaceFrameNetwork, self).__init__()
        self.dist_dict = {}
        self.update_default_edge_attributes(created_from=None)
        self.update_default_node_attributes(created_from=None)

    def transform_network(self, transformation):
        nodes = [self.node_coordinates(key) for key in self.nodes()]
        xyz = transform_points(nodes, transformation)
        for index, (key, attr) in enumerate(self.nodes(True)):
            attr["x"] = xyz[index][0]
            attr["y"] = xyz[index][1]
            attr["z"] = xyz[index][2]

        self.build_dist_dict()

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
            # print("Not ok_conn because has_edge")
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
    def from_space_frame_mesh(cls, mesh, scale="mm"):
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
            x, y, z = mesh.vertex_coordinates(vkey)
            # data = mesh.vertex_attributes(vkey)
            # data.update({"created_from": cls.FROM_MESH})
            network.add_node(key=vkey, x=x, y=y, z=z)

        for u, v in mesh.edges():
            # data = mesh.edge_attributes((u, v))
            # data.update({"created_from": cls.FROM_MESH})
            network.add_edge(u, v)

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


class SpaceFrameMesh(SpaceFrameMixin, Mesh):
    FROM_SUBD_TRI = 0
    FROM_SUBD_FRAME = 1

    def __init__(self):
        super(SpaceFrameMesh, self).__init__()
        self.update_default_face_attributes(parent_fkey=None, created_from=None)
        self.update_default_vertex_attributes(parent_fkey=None, created_from=None)
        self.update_default_edge_attributes(parent_fkey=None, created_from=None)

    def delete_face(self, fkey):
        try:
            super(SpaceFrameMesh, self).delete_face(fkey)
        except KeyError:
            pass

    """
    def collapse_edges(self, min_edge_length):
        edges = set(list(self.edges()))
        # print(len(edges))
        edges -= set(self.edges_on_boundary(oriented=True))
        edges = [(u, v) for u, v in edges if self.edge_length(u, v) < min_edge_length]

        for u, v in edges:
            if self.has_edge((u, v)) and self.has_vertex(u) and self.has_vertex(v):
                # t = 0 if self.vertex_degree(u) < self.vertex_degree(v) else 1
                t = 0.5
                self.collapse_edge(u, v, t=t)
        """

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
                print("Not ok_subd due to vertex degrees")
                return False

        if min_angle:
            if not self.ok_edges_angles(flatten([old_vkeys + new_vkeys]), min_angle):
                print("Not ok_subd due to edge angles.")
                return False

        if edge_length_domain:
            ekeys = [self.connected_edges(vkey) for vkey in new_vkeys]
            ekeys = flatten(ekeys)
            uv_sets = set(frozenset((u, v)) for u, v in ekeys)  # unique edges
            if not self.ok_edge_length(uv_sets, edge_length_domain):
                print("Not ok_subd due to edge length constraints.")
                return False

        return True

    def _move_center_pt(self, pt, fkey, rel_pos, move_z):
        if not rel_pos and not move_z:
            return pt

        if rel_pos:
            rel_pos = cycle(to_list(rel_pos))
            for v_xyz in self.face_coordinates(fkey):
                v = pt - Point(*v_xyz)

                factor = next(rel_pos)
                pt += v * factor

        if move_z:
            # set up max dist
            normal = Vector(*self.face_normal(fkey, unitized=True))
            z_vec = normal * move_z
            pt += z_vec

        return pt

    def subdiv_faces(
        self, fkeys, n_iters, scheme=[0], rel_pos=None, move_z=None, **kwargs
    ):
        subdiv_funcs = [
            self.face_subdiv_pyramid,
            self.face_subdiv_mid_cent,
            self.face_subdiv_mids,
            self.face_subdiv_split,
            self.face_subdiv_frame,
        ]
        repeating_n_iters = cycle(n_iters)
        subdiv_cycler = cycle(scheme)
        repeating_rel_pos = cycle(to_list(rel_pos)) if rel_pos else None
        repeating_move_z = cycle(to_list(move_z)) if move_z else None

        fkeys = set(fkeys)
        while len(fkeys) > 0:
            fkey = fkeys.pop()
            n_iters = next(repeating_n_iters)
            if repeating_rel_pos:
                kwargs.update({"rel_pos": next(repeating_rel_pos)})
            if repeating_move_z:
                kwargs.update({"move_z": next(repeating_move_z)})

            subdiv_parent_face = next(subdiv_cycler)
            if not hasattr(subdiv_parent_face, "__next__"):
                subdiv_parent_face = cycle(to_list(subdiv_parent_face))

            i = 0
            next_to_subd = set([fkey])
            while i < n_iters:

                to_subd = next_to_subd
                next_to_subd = set()

                subdiv_iteration = next(subdiv_parent_face)
                if not hasattr(subdiv_iteration, "__next__"):
                    subdiv_iteration = cycle(to_list(subdiv_iteration))

                while len(to_subd) > 0:
                    parent_fkey = to_subd.pop()
                    parent_face_verts = self.face_vertices(parent_fkey)
                    # parent_attrs = self.face_attributes(parent_fkey)

                    subdiv_child_face = next(subdiv_iteration)

                    subdiv_func = subdiv_funcs[subdiv_child_face]

                    new_vkeys, new_fkeys, deleted_faces = subdiv_func(
                        parent_fkey, **kwargs
                    )

                    if not self.ok_subd(parent_face_verts, new_vkeys, **kwargs):
                        self.undo_face_subd(
                            new_fkeys, new_vkeys, deleted_faces,
                        )
                        continue

                    deleted_fkeys = set([fkey for fkey, _ in deleted_faces])

                    # add new fkeys to set for next iteration
                    next_to_subd |= set(new_fkeys)
                    next_to_subd -= deleted_fkeys

                    # remove deleted faces from input list of fkeys
                    fkeys -= deleted_fkeys

                    # remove deleted faces from list for this iteration
                    to_subd -= deleted_fkeys

                i += 1

        self.cull_vertices()

    def face_subdiv_pyramid(self, fkey, rel_pos=None, move_z=None, **kwattr):
        xyz = Point(*self.face_center(fkey))

        x, y, z = self._move_center_pt(xyz, fkey, rel_pos, move_z)

        face_halfedges = self.face_halfedges(fkey)

        deleted_faces = [(fkey, self.face_vertices(fkey))]
        self.delete_face(fkey)

        center_vkey = self.add_vertex(x=x, y=y, z=z)

        new_fkeys = []
        for u, v in face_halfedges:
            vkeys = [u, v, center_vkey]
            new_fkeys.append(self.add_face(vkeys))

        return [center_vkey], new_fkeys, deleted_faces

    def _split_edges(self, fkey, halfedges):
        new_fkeys = []
        new_vkeys = []
        deleted_faces = []

        for u, v in halfedges:
            adj_faces = [key for key in self.edge_faces(u, v) if key != fkey]
            adj_face = adj_faces.pop()
            adj_face_verts = self.face_vertices(adj_face) if adj_face else None

            new_vkeys.append(self.split_edge(u, v, allow_boundary=True))

            if adj_face:
                deleted_faces.append((adj_face, adj_face_verts))

                split_from = new_vkeys[-1]
                split_to = self.face_vertex_descendant(adj_face, split_from, n=2)

                self.split_face(adj_face, split_from, split_to)

                new_fkeys += self.edge_faces(split_from, split_to)
                self.delete_face(adj_face)

        return new_vkeys, new_fkeys, deleted_faces

    def face_subdiv_mid_cent(self, fkey, rel_pos=None, move_z=None, **kwattr):
        xyz = Point(*self.face_center(fkey))

        x, y, z = self._move_center_pt(xyz, fkey, rel_pos, move_z)

        face_halfedges = self.face_halfedges(fkey)

        new_vkeys, new_fkeys, deleted_faces = self._split_edges(fkey, face_halfedges)

        deleted_faces.append((fkey, self.face_vertices(fkey)))
        self.delete_face(fkey)

        center_vkey = self.add_vertex(x=x, y=y, z=z)

        for i, face_halfedge in enumerate(face_halfedges):
            face_vertex, _ = face_halfedge
            edge_mid = new_vkeys[i]
            prev_edge_mid = new_vkeys[(i - 1) % len(new_vkeys)]

            vkeys = [face_vertex, edge_mid, center_vkey, prev_edge_mid]
            new_fkeys.append(self.add_face(vkeys))

        new_vkeys.append(center_vkey)

        return new_vkeys, new_fkeys, deleted_faces

    def face_subdiv_mids(self, fkey, **kwattr):
        # remove rel_pos and rel_pos_z since they don't apply
        kwattr.pop("rel_pos", None)
        kwattr.pop("rel_pos_z", None)

        face_halfedges = self.face_halfedges(fkey)

        new_vkeys, new_fkeys, deleted_faces = self._split_edges(fkey, face_halfedges)

        deleted_faces.append((fkey, self.face_vertices(fkey)))
        self.delete_face(fkey)

        for i, halfedges in enumerate(face_halfedges):
            face_vertex, _ = halfedges
            prev_edge_mid = new_vkeys[(i - 1) % len(new_vkeys)]
            edge_mid = new_vkeys[i]
            vkeys = [prev_edge_mid, face_vertex, edge_mid]
            new_fkeys.append(self.add_face(vkeys))

        # middle face
        new_fkeys.append(self.add_face(new_vkeys))

        return new_vkeys, new_fkeys, deleted_faces

    def face_subdiv_split(
        self, fkey, tri_split_equal=True, shift_split=False, **kwattr
    ):
        # remove rel_pos and rel_pos_z since they don't apply
        kwattr.pop("rel_pos", None)
        kwattr.pop("rel_pos_z", None)

        face_halfedges = self.face_halfedges(fkey)
        # print("fkeys: {}, halfedges: {}".format(face_vertices, face_halfedges))

        u1, v1 = face_halfedges[0 + shift_split]

        new_vkeys, new_fkeys, deleted_faces = self._split_edges(fkey, [(u1, v1)])
        deleted_faces.append((fkey, self.face_vertices(fkey)))

        e_mid1 = new_vkeys[-1]

        if len(face_halfedges) == 3:
            u2, v2 = face_halfedges[1 + shift_split]
            if tri_split_equal:
                # tri split in half

                new_face_verts1 = [u1, e_mid1, v2]
                new_face_verts2 = [e_mid1, v1, v2]
            else:
                # tri face split into one quad and one tri

                _new_vkeys, _new_fkeys, _deleted_faces = self._split_edges(
                    fkey, [(u2, v2)]
                )
                new_vkeys += _new_vkeys
                new_fkeys += _new_fkeys
                deleted_faces += _deleted_faces

                e_mid2 = new_vkeys[-1]

                new_face_verts1 = [u1, e_mid1, e_mid2, v2]
                new_face_verts2 = [e_mid1, u2, e_mid2]

        else:
            # Quad face into two smaller quad faces

            # Get edge opposite (u1, v1)
            u2, v2 = face_halfedges[2 + shift_split]

            _new_vkeys, _new_fkeys, _deleted_faces = self._split_edges(fkey, [(u2, v2)])
            new_vkeys += _new_vkeys
            new_fkeys += _new_fkeys
            deleted_faces += _deleted_faces

            e_mid2 = new_vkeys[-1]

            new_face_verts1 = [u1, e_mid1, e_mid2, v2]
            new_face_verts2 = [e_mid1, v1, u2, e_mid2]

        self.delete_face(fkey)

        new_fkeys.append(self.add_face(new_face_verts1))
        new_fkeys.append(self.add_face(new_face_verts2))

        return new_vkeys, new_fkeys, deleted_faces

    def face_subdiv_frame(self, fkey, rel_pos=None, move_z=None, **kwattr):
        if rel_pos == 1 or rel_pos == [1, 1, 1]:
            return self.face_subdiv_pyramid(fkey, move_z=move_z, **kwattr)

        face_center_pt = Point(*self.face_center(fkey))
        face_normal = Vector(*self.face_normal(fkey))
        face_coordinates = self.face_coordinates(fkey)
        face_halfedges = self.face_halfedges(fkey)

        deleted_faces = [(fkey, self.face_vertices(fkey))]
        self.delete_face(fkey)

        if move_z is None:
            move_z = cycle(to_list(0))

        if not isinstance(move_z, cycle):
            move_z = cycle(to_list(move_z))

        if rel_pos is None:
            rel_pos = cycle(to_list(0.5))

        if not isinstance(rel_pos, cycle):
            rel_pos = cycle(to_list(rel_pos))

        new_vkeys = []
        for x, y, z in face_coordinates:
            pt = Point(x, y, z)
            factor = next(rel_pos)

            v = face_center_pt - pt
            pt += v * factor
            if move_z:
                z_factor = next(move_z)
                pt += face_normal * z_factor
            new_vkeys.append(self.add_vertex(x=pt.x, y=pt.y, z=pt.z))

        new_fkeys = []
        for i, uv in enumerate(face_halfedges):
            u, v = uv
            vkeys = []
            vkeys.append(u)
            vkeys.append(v)
            vkeys.append(new_vkeys[(i + 1) % len(new_vkeys)])
            vkeys.append(new_vkeys[i])

            new_fkeys.append(self.add_face(vkeys))

        # add new center face
        new_fkeys.append(self.add_face(new_vkeys))

        return new_vkeys, new_fkeys, deleted_faces

    def undo_face_subd(self, new_fkeys, new_vkeys, deleted_faces):
        # TODO: Figure out how a new vkey can also be part of deleted face
        deleted_fkeys, deleted_face_verts = zip(*deleted_faces)
        dups_removed = []
        for face_verts in deleted_face_verts:
            face_verts = [key for key in face_verts if key not in new_vkeys]
            dups_removed.append(face_verts)
        deleted_faces = zip(deleted_fkeys, dups_removed)

        for fkey in new_fkeys:
            self.delete_face(fkey)

        for vkey in new_vkeys:
            self.delete_vertex(vkey)

        for fkey, face_vertices in deleted_faces:
            try:
                self.add_face(face_vertices, fkey=fkey)
            except KeyError:
                print(fkey, face_vertices)
                print(
                    "Face key already present: {}. Vertex exists in mesh: {}".format(
                        self.has_face(fkey),
                        [self.has_vertex(key) for key in face_vertices],
                    )
                )
                raise

    def to_rgmesh(self):
        artist = MeshArtist(self)
        return artist.draw_mesh()

    def to_rgmesh_meters(self):
        raise NotImplementedError("Scaling needs fixing.")
        rg_mesh = self.to_rgmesh()
        T = self.get_scale_rgxform(0.001)
        rg_mesh.Transform(T)

        return rg_mesh
