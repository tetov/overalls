from itertools import chain
from itertools import cycle

import Rhino.Geometry as rg
from compas.datastructures import Mesh
from compas.datastructures import Network
from compas.geometry import Point
from compas.geometry import Vector
from compas.geometry import angle_vectors
from overalls.utils import flatten
from overalls.utils import is_in_domain
from overalls.utils import to_list


class ExtraEdge(object):
    def __init__(self, u, v, **kwattr):
        if u == v:
            raise Exception("Need two different vertices to create an edge.")

        self.u = u
        self.v = v

        self.attr = kwattr

    @property
    def conn_dict(self):
        return {self.u: self.v, self.v: self.u}

    def __eq__(self, other):
        return self.conn_dict == other.conn_dict


class SpaceFrame(Mesh):
    FROM_SUBD_TRI = 0
    FROM_SUBD_FRAME = 1
    FROM_CONN_CENTERS = 2
    FROM_CONN_VERTICES = 3

    def __init__(self):
        super(SpaceFrame, self).__init__()
        self.update_default_face_attributes(parent_fkey=None, created_from=None)
        self.update_default_vertex_attributes(parent_fkey=None, created_from=None)
        self.update_default_edge_attributes(parent_fkey=None, created_from=None)
        self._extra_edges = []

    def extra_edges(self):
        for edge in self._extra_edges:
            yield (edge.u, edge.v)

    def mesh_edges(self):
        return super(SpaceFrame, self).edges()

    def edges(self):
        return chain(self.mesh_edges(), self.extra_edges())

    def add_edge(self, u, v, **kwattr):
        edge = ExtraEdge(u, v, **kwattr)
        self._extra_edges.append(edge)

    def edge_attributes(self, ekey, keys, values):
        u, v = ekey
        for edge in self.extra_edges():
            tmp_obj = ExtraEdge(u, v)

            if tmp_obj == edge:
                new_attrs = {key: value for (key, value) in zip(keys, values)}
                data = edge.attr
                data.update(new_attrs)
                edge.attr = data
                break
        super(SpaceFrame, self).edge_attributes((u, v), keys, values)

    def edge_attribute(self, ekey, key, value):
        self.edge_attributes(ekey, to_list(key), to_list(value))

    def components_attributes(self, fkeys, vkeys, ekeys, attr):
        for vkey in vkeys:
            self.vertex_attributes(vkey, attr.keys(), attr.values())
        for u, v in ekeys:
            self.edge_attributes((u, v), attr.keys(), attr.values())
        for fkey in fkeys:
            self.face_attributes(fkey, attr.keys(), attr.values())

    def vertex_edges(self, u):
        nbors = self.vertex_neighbors(u)
        nbors += [
            edge.conn_dict[u]
            for edge in self.extra_edges()
            if u in edge.conn_dict.keys()
        ]

        return [(u, v) for v in nbors]

    def ok_vertex_degree(self, vkeys, max_degree):
        vkeys = to_list(vkeys)
        for vkey in vkeys:
            if self.vertex_degree(vkey) > max_degree:
                return False
        return True

    def ok_vertex_angles(self, vkeys, min_angle):
        vkeys = flatten(to_list(vkeys))
        for vkey in vkeys:
            edges = self.vertex_edges(vkey)
            vectors = [self.edge_vector(u, v) for u, v in edges]

            for v in vectors:
                to_check = [v1 for v1 in vectors if v1 != v]
                for v1 in to_check:
                    if abs(angle_vectors(v, v1, deg=True)) < min_angle:
                        return False
        return True

    def ok_edge_length(self, ekeys, length_domain):
        for u, v in ekeys:
            if not is_in_domain(self.edge_length(u, v), length_domain):
                return False
        return True

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
            if not self.ok_vertex_degree(old_vkeys, max_degree):
                return False

        if min_angle:
            if not self.ok_vertex_angles([old_vkeys + new_vkeys], min_angle):
                return False

        if edge_length_domain:
            ekeys = [self.vertex_edges(vkey) for vkey in new_vkeys]
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

    def find_closest_faces(
        self,
        fkeys_origin,
        fkeys_search,
        n_face_connections=2,
        dist_domain=(None, None),
        prefer_long=False,
    ):
        pt_cloud_dict = {}

        for fkey in fkeys_search:
            x, y, z = self.face_center(fkey)
            pt = rg.Point3d(x, y, z)
            pt_cloud_dict[pt] = fkey

        partners = []

        for fkey in fkeys_origin:
            dists = []
            x, y, z = self.face_center(fkey)
            pt1 = rg.Point3d(x, y, z)

            for pt2, key in pt_cloud_dict.iteritems():
                if key == fkey or key in self.face_neighbors(fkey):
                    continue
                dist = pt1.DistanceTo(pt2)
                if is_in_domain(dist, dist_domain):
                    dists.append((key, dist))

            if len(dists) < 1:  # No matches
                continue

            dists.sort(key=lambda x: x[1], reverse=prefer_long)
            keys, _ = zip(*dists)
            partners.append((fkey, keys[:n_face_connections]))
        return partners

    def find_closest_face_nhood(self, fkeys_origin, fkeys_search, *args, **kwargs):
        original_partners = self.find_closest_faces(
            fkeys_origin, fkeys_search, *args, **kwargs
        )

        partners = []
        for partner_a, partners_b in original_partners:
            nbors = set()
            for partner in partners_b:
                nbors.update(self.face_neighbors(partner))
            partners.append((partner_a, list(nbors)))
        return partners

    def connect_faces_center_centers(
        self, start_fkey, end_fkeys, max_degree=None,
    ):

        # connect center to center
        x, y, z = self.face_center(start_fkey)
        start_vkey, _ = self.face_subdiv_retriangulate(
            start_fkey, created_from=self.FROM_CONN_CENTERS
        )

        for fkey in end_fkeys:
            vkey, _ = self.face_subdiv_retriangulate(
                fkey, return_vkey=True, created_from=self.FROM_CONN_CENTERS
            )
            self.add_edge(
                start_vkey,
                vkey,
                created_from=self.FROM_CONN_CENTERS,
                parent_fkey=(start_fkey, fkey),
            )

    def connect_faces_vert_vert(
        self,
        start_fkeys,
        end_fkeys,
        rel_start_vkeys=None,
        rel_end_vkeys=None,
        max_degree=None,
    ):
        """Create a line between a vertex on one mesh to a vertex on another."""
        start_fkeys = to_list(start_fkeys)
        end_fkeys = to_list(end_fkeys)
        if rel_start_vkeys:
            rel_start_vkeys = to_list(rel_start_vkeys)
        if rel_end_vkeys:
            rel_end_vkeys = to_list(rel_end_vkeys)

        max_iters = max(len(start_fkeys), len(end_fkeys))

        for j in range(max_iters):
            start_fkey = start_fkeys[j % len(start_fkeys)]
            end_fkey = end_fkeys[j % len(end_fkeys)]

            start_vkeys = self.face_vertices(start_fkey)
            end_vkeys = self.face_vertices(end_fkey)

            if rel_start_vkeys:
                start_key_idx = rel_start_vkeys[j % len(rel_start_vkeys)]
            else:
                start_vkeys.sort(key=lambda x: self.vertex_degree(x))
                start_key_idx = 0
            if rel_end_vkeys:
                end_key_idx = rel_end_vkeys[j % len(rel_end_vkeys)]
            else:
                end_vkeys.sort(key=lambda x: self.vertex_degree(x))
                end_key_idx = 0

            u = start_vkeys[start_key_idx]
            v = end_vkeys[end_key_idx]

            if max_degree:
                u_ok = self.vertex_degree(u) <= max_degree
                v_ok = self.vertex_degree(v) <= max_degree
                if not u_ok or not v_ok:
                    continue

            self.add_edge(
                u,
                v,
                created_from=self.FROM_CONN_VERTICES,
                parent_fkey=(start_fkey, end_fkey),
            )

    def edge_to_rgline(self, u, v):
        pt_a, pt_b = self.edge_coordinates(u, v)
        pt_a = rg.Point3d(*pt_a)
        pt_b = rg.Point3d(*pt_b)
        return rg.Line(pt_a, pt_b)

    def to_rglines(self, separate_conn_lines=False):
        lines = []

        if separate_conn_lines:
            conn_lines = []
            for u, v in self.extra_edges():
                conn_lines.append(self.edge_to_rgline(u, v))
            edges = self.mesh_edges()
        else:
            edges = self.edges()

        for u, v in edges:
            lines.append(self.edge_to_rgline(u, v))

        if separate_conn_lines:
            return lines, conn_lines

        return lines

    def to_network(self):
        network = Network()

        for vkey in self.vertices():
            x, y, z = self.vertex_coordinates(vkey)
            network.add_node(x=x, y=y, z=z)
        for u, v in self.edges():
            network.add_edge(u, v)

        return network
