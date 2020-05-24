from collections import Sequence
from itertools import cycle
from itertools import ifilter

import Rhino.Geometry as rg
from compas.datastructures import Mesh
from compas.geometry import Point
from compas.geometry import Vector
from compas.utilities import geometric_key
from compas.utilities import pairwise
from compas_ghpython.artists import MeshArtist
from overalls.space_frame._mixin import SpaceFrameMixin
from overalls.utils import flatten
from overalls.utils import to_list


class SpaceFrameMesh(Mesh, SpaceFrameMixin):
    FROM_SUBD_TRI = 0
    FROM_SUBD_FRAME = 1

    def __init__(self):
        super(SpaceFrameMesh, self).__init__()

        self.update_default_face_attributes(
            parent_fkey=None, created_from=None, n_iters=0, part_id=None
        )
        self.update_default_vertex_attributes(
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
        return self.vertex_degree(key, *args, **kwargs)

    def nvs(self, *args, **kwargs):
        return self.vertices(*args, **kwargs)

    def nv_neighbors(self, key, *args, **kwargs):
        return self.vertex_neighbors(key, *args, **kwargs)

    def nv_coordinates(self, key, *args, **kwargs):
        return self.vertex_coordinates(key, *args, **kwargs)

    def nv_attribute(self, u, v, key, *args, **kwargs):
        return self.vertex_attribute(u, v, key, *args, **kwargs)

    def update_default_nv_attributes(self, **kwattr):
        return self.update_default_vertex_attributes(**kwattr)

    def delete_face(self, fkey):
        try:
            super(SpaceFrameMesh, self).delete_face(fkey)
        except KeyError:
            pass

    def weld(self, precision=None):
        """Adapted from compas.datastructures.mesh.join"""
        geo = geometric_key
        key_xyz = {key: self.vertex_coordinates(key) for key in self.vertices()}
        gkey_key = {
            geo(xyz, precision): (key, self.vertex_attributes(key))
            for key, xyz in key_xyz.items()
        }
        gkey_index = {gkey: index for index, gkey in enumerate(gkey_key)}

        vertices = [
            (key_xyz[key_data[0]], key_data[1]) for gkey, key_data in gkey_key.items()
        ]

        faces = []
        for fkey in self.faces():
            fdata = self.face_attributes(fkey)
            new_vkeys = []
            for key in self.face_vertices(fkey):
                coords = key_xyz[key]
                gkey = geo(coords)
                idx = gkey_index[gkey]
                new_vkeys.append(idx)
            faces.append((new_vkeys, fdata))

        faces[:] = [
            ([u for u, v in pairwise(face + face[:1]) if u != v], fdata)
            for face, fdata in faces
        ]

        original_vertices = list(self.vertices())
        for vkey in original_vertices:
            self.delete_vertex(vkey)
        self.cull_vertices()

        for fkey in self.faces():
            print("{} still exists".format(fkey))

        for xyz, data in vertices:
            data["x"], data["y"], data["z"] = xyz
            self.add_vertex(key=gkey_index[geo(xyz, precision)], **data)

        for vkeys, data in faces:
            fdata = {}
            for name in data:
                fdata[name] = data.get(name)
            self.add_face(vkeys, **fdata)

        self.cull_vertices()
        self.unify_cycles()

    def ok_subd(
        self,
        old_vkeys,
        new_vkeys,
        new_fkeys,
        degree_domain=None,
        min_angle=None,
        edge_length_domain=None,
        max_edge_ratio_diff=None,
        **kwargs
    ):
        new_vkeys = to_list(new_vkeys)

        if degree_domain:
            if not self.ok_degrees(old_vkeys, degree_domain):
                # print("Not ok_subd due to vertex degrees")
                return False

        if min_angle:
            if not self.ok_edges_angles(flatten([old_vkeys + new_vkeys]), min_angle):
                # print("Not ok_subd due to edge angles.")
                return False

        if edge_length_domain:
            ekeys = [self.connected_edges(vkey) for vkey in new_vkeys]
            ekeys = flatten(ekeys)
            uv_sets = set(frozenset((u, v)) for u, v in ekeys)  # unique edges
            if not self.ok_edge_lengths(uv_sets, edge_length_domain):
                # print("Not ok_subd due to edge length constraints.")
                return False

        if max_edge_ratio_diff:
            for fkey in new_fkeys:
                if not self.ok_edge_length_ratios(fkey, max_edge_ratio_diff):
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

    def estimate_subd_face_edge_angles(self, fkey):
        est_vec_angles = []
        for vkey in self.face_vertices(fkey):
            prev_vkey = self.face_vertex_ancestor(fkey, vkey)
            next_vkey = self.face_vertex_descendant(fkey, vkey)

            origin_xyz = self.vertex_coordinates(vkey)
            prev_xyz = self.vertex_coordinates(prev_vkey)
            next_xyz = self.vertex_coordinates(next_vkey)

            vec_a = rg.Vector3d(*prev_xyz) - rg.Vector3d(*origin_xyz)
            vec_b = rg.Vector3d(*next_xyz) - rg.Vector3d(*origin_xyz)

            vec_angle = rg.Vector3d.VectorAngle(vec_a, vec_b)
            est_vec_angles.append(vec_angle / 2)
        return est_vec_angles

    def estimate_subd_edge_length(self, fkey):
        cur_lengths = [self.edge_length(u, v) for u, v in self.face_halfedges(fkey)]

        avg_cur_length = sum(cur_lengths) / len(cur_lengths)

        return avg_cur_length / 2

    @staticmethod
    def _get_next_cycle(var):
        var = next(var)
        if hasattr(var, "__next__"):
            return var
        if not isinstance(var, Sequence):
            var = to_list(var)
        return cycle(var)

    def subdiv_faces_static(
        self,
        fkeys,
        n_iters=None,
        scheme=[0],
        rel_pos=None,
        max_iters=None,
        move_z=None,
        degree_domain=None,
        min_angle=None,
        edge_length_domain=None,
        **kwargs
    ):
        # repeating_n_iters = cycle(n_iters) if n_iters else None
        subdiv_cycler = cycle(scheme)
        # repeating_rel_pos = cycle(to_list(rel_pos)) if rel_pos else None
        # repeating_move_z = cycle(to_list(move_z)) if move_z else None

        subdiv_funcs = [
            self.face_subdiv_pyramid,
            self.face_subdiv_mid2center,
            self.face_subdiv_mids2mids,
            self.face_subdiv_split_in_half,
            self.face_subdiv_split_mid2mid,
            self.face_subdiv_frame,
        ]

        if not n_iters:
            raise NotImplementedError

        fkeys = set(fkeys)

        next_iter = fkeys
        deleted_fkeys = set()
        i = 0
        while i < n_iters:
            i += 1

            subdiv_func_per_iter = self._get_next_cycle(subdiv_cycler)
            to_subd = [fkey for fkey in next_iter if fkey not in deleted_fkeys]
            next_iter = []
            filtered_to_subd = []
            for fkey in to_subd:
                face_verts = self.face_vertices(fkey)
                halfedges = self.face_halfedges(fkey)

                # check current state versus requirements
                degree_min, degree_max = degree_domain
                degree_min = degree_min - 1 if degree_min else None
                degree_max = degree_max + 1 if degree_max else None
                if not self.ok_degrees(face_verts, (degree_min, degree_max)):
                    continue

                if not self.ok_edges_angles(face_verts, min_angle):
                    continue

                if not self.ok_edge_lengths(halfedges, (edge_length_domain[0], None)):
                    continue

                # check estimated subd state versus requirements
                new_edge_angles = self.estimate_subd_face_edge_angles(fkey)
                new_edge_angles.sort()

                if new_edge_angles[0] < min_angle:
                    continue

                new_edge_length = self.estimate_subd_edge_length(fkey)

                if new_edge_length < edge_length_domain[0]:
                    continue

                filtered_to_subd.append(fkey)

            for fkey in filtered_to_subd:
                # if not self.has_face(fkey):  # might have been deleted in prev subd
                # continue
                if fkey in deleted_fkeys:
                    continue

                subdiv_func_idx = next(subdiv_func_per_iter)
                subdiv_func = subdiv_funcs[subdiv_func_idx]

                new_vkeys, new_fkeys, deleted_faces = subdiv_func(fkey, **kwargs)

                newly_deleted_fkeys, _ = zip(*deleted_faces)
                deleted_fkeys |= set(newly_deleted_fkeys)

                next_iter += new_fkeys

    def subdiv_faces(
        self,
        fkeys,
        n_iters=None,
        scheme=[0],
        rel_pos=None,
        move_z=None,
        max_iters=None,
        **kwargs
    ):
        subdiv_funcs = [
            self.face_subdiv_pyramid,
            self.face_subdiv_mid2center,
            self.face_subdiv_mids2mids,
            self.face_subdiv_split_in_half,
            self.face_subdiv_split_mid2mid,
            self.face_subdiv_frame,
        ]
        repeating_n_iters = cycle(n_iters) if n_iters else None
        subdiv_cycler = cycle(scheme)
        repeating_rel_pos = cycle(to_list(rel_pos)) if rel_pos else None
        repeating_move_z = cycle(to_list(move_z)) if move_z else None

        fkeys = set(fkeys)

        while len(fkeys) > 0:
            fkey = fkeys.pop()

            if n_iters:
                n_iters = next(repeating_n_iters)
            elif self.face_attribute(fkey, "n_iters"):
                n_iters = self.face_attribute(fkey, "n_iters")
            else:
                n_iters = 0

            if max_iters:
                n_iters = n_iters if n_iters < max_iters else max_iters

            if repeating_rel_pos:
                kwargs.update({"rel_pos": next(repeating_rel_pos)})
            if repeating_move_z:
                kwargs.update({"move_z": next(repeating_move_z)})

            subdiv_func_per_parent_face = self._get_next_cycle(subdiv_cycler)

            i = 0
            next_to_subd = set([fkey])
            while i < n_iters:

                to_subd = next_to_subd
                next_to_subd = set()

                while len(to_subd) > 0:
                    parent_fkey = to_subd.pop()
                    parent_face_verts = self.face_vertices(parent_fkey)
                    part = self.face_attribute(parent_fkey, "part")
                    # parent_attrs = self.face_attributes(parent_fkey)

                    subd_func_idx = next(subdiv_func_per_parent_face)

                    subdiv_func = subdiv_funcs[subd_func_idx]

                    new_vkeys, new_fkeys, deleted_faces = subdiv_func(
                        parent_fkey, **kwargs
                    )

                    for vkey in new_vkeys:
                        self.vertex_attribute(vkey, "part", part)

                    if not self.ok_subd(
                        parent_face_verts, new_vkeys, new_fkeys, **kwargs
                    ):
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

                new_fkeys += self.split_face(adj_face, split_from, split_to)

                self.edge_attribute((split_from, split_to), "ignore_edge", True)

                self.delete_face(adj_face)

        return new_vkeys, new_fkeys, deleted_faces

    def face_subdiv_mid2center(self, fkey, rel_pos=None, move_z=None, **kwattr):
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

    def face_subdiv_mids2mids(self, fkey, **kwattr):
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

    def face_subdiv_split_in_half(self, fkey, shift_split=False, **kwattr):
        # remove rel_pos and rel_pos_z since they don't apply
        kwattr.pop("rel_pos", None)
        kwattr.pop("rel_pos_z", None)

        new_vkeys = []
        new_fkeys = []
        deleted_faces = [(fkey, self.face_vertices(fkey))]

        face_halfedges = self.face_halfedges(fkey)

        u1, v1 = face_halfedges[0 + shift_split]

        if len(face_halfedges) < 4:
            _new_vkeys, _new_fkeys, _deleted_faces = self._split_edges(fkey, [(u1, v1)])
            new_vkeys += _new_vkeys
            new_fkeys += _new_fkeys
            deleted_faces += _deleted_faces

            e_mid1 = new_vkeys[-1]

            u2, v2 = face_halfedges[1 + shift_split]

            new_face_verts1 = [u1, e_mid1, v2]
            new_face_verts2 = [e_mid1, v1, v2]

            new_fkeys.append(self.add_face(new_face_verts1))
            new_fkeys.append(self.add_face(new_face_verts2))

        else:
            split_from = u1
            split_to_idx = len(face_halfedges) // 2  # support ngon for fun
            split_to = self.face_vertex_ancestor(fkey, u1, n=split_to_idx)

            new_fkeys += self.split_face(fkey, split_from, split_to)

        self.delete_face(fkey)

        return new_vkeys, new_fkeys, deleted_faces

    def face_subdiv_split_mid2mid(self, fkey, shift_split=False, **kwattr):
        # remove rel_pos and rel_pos_z since they don't apply
        kwattr.pop("rel_pos", None)
        kwattr.pop("rel_pos_z", None)

        face_halfedges = self.face_halfedges(fkey)
        # print("fkeys: {}, halfedges: {}".format(face_vertices, face_halfedges))

        u1, v1 = face_halfedges[0 + shift_split]

        new_vkeys, new_fkeys, deleted_faces = self._split_edges(fkey, [(u1, v1)])
        deleted_faces.append((fkey, self.face_vertices(fkey)))

        e_mid1 = new_vkeys[-1]

        u2, v2 = face_halfedges[2 + shift_split]

        _new_vkeys, _new_fkeys, _deleted_faces = self._split_edges(fkey, [(u2, v2)])
        new_vkeys += _new_vkeys
        new_fkeys += _new_fkeys
        deleted_faces += _deleted_faces

        e_mid2 = new_vkeys[-1]

        if len(face_halfedges) == 3:

            new_face_verts1 = [u1, e_mid1, e_mid2, v2]
            new_face_verts2 = [e_mid1, u2, e_mid2]

        else:
            # Quad face into two smaller quad faces

            new_face_verts1 = [u1, e_mid1, e_mid2, v2]
            new_face_verts2 = [e_mid1, v1, u2, e_mid2]

        self.delete_face(fkey)

        new_fkeys.append(self.add_face(new_face_verts1))
        new_fkeys.append(self.add_face(new_face_verts2))

        return new_vkeys, new_fkeys, deleted_faces

    def edge_split_subd(self, u, v, to_trifaces=True):
        """Split a edge and adjacent faces.

        Parameters
        ----------
        u, v
            Edge key
        to_trifaces : :obj:`bool`, optional
            If subdivision should continue until all resulting faces are
            triangular.

        Returns
        -------
        :obj:`list` of :obj:`int`
            Identifier for the new vertex.
        :obj:`list` of :obj:`int`
            Identifiers for new faces
        :obj:`list`
            Placeholder for `deleted_faces`
        """
        fkeys = self.edge_faces(u, v)

        # remove None values
        fkeys = list(ifilter(None, fkeys))

        new_vkey = self.split_edge(u, v, allow_boundary=True)
        new_fkeys = []

        for fkey in fkeys:
            split_from = new_vkey
            split_to = self.face_vertex_descendant(fkey, split_from, n=2)

            new_fkeys += self.split_face(fkey, split_from, split_to)

            if to_trifaces:
                for fkey in new_fkeys[-2:]:
                    split_from = new_vkey
                    split_to = self.face_vertex_descendant(fkey, split_from, n=2)

                    try:
                        new_fkeys += self.split_face(fkey, split_from, split_to)
                    except ValueError:  # face is already triface
                        pass

        # Imitating the return values of face subds
        return [new_vkey], list(ifilter(None, new_fkeys)), []

    def face_subdiv_frame(self, fkey, rel_pos=None, move_z=None, **kwattr):
        """Subdivide a face by offsetting its edges inwards

        Creates ``verts+1`` new faces.

        Parameters
        ----------
        fkey : :obj:`int`
        rel_pos: :obj:`list` of :obj:`float`
        rel_pos: :obj:`float`

        Returns
        -------
        :obj:`list` of :obj:`int`
            Keys of newly created vertices.
        :obj:`list` of :obj:`int`
            Keys of newly created faces.
        :obj:`tuple` of :obj:`int` and :obj:`tuple` of :obj:`int`
            Face keys of removed faces and their vertex keys.
        """

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
        """Reverse a face subdivision.

        Parameters
        ----------
        new_fkeys : :obj:`list` of :obj:`int`
            Keys of faces created in last subdivision.
        new_vkeys : :obj:`list` of :obj:`int`
            Keys of vertices created in last subdivision.
        deleted_faces : :obj:`list` of :obj:`tuple`
            Keys of faces that were removed in last subdivision and their
            vertex key. E.g. ``(fkey, (vkey1, vkey2, vkey3))``.

        Raises
        ------
        :exc:`KeyError`
            If vertex key from deleted_faces is not found in mesh.
        """
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

    def n_iters_from_edge_data(
        self, face_label="n_iters", edge_label="structural_data"
    ):
        self.update_default_edge_attributes({edge_label: None})
        self.update_default_face_attributes({face_label: 0})
        for u, v in self.edges():
            data = self.edge_attribute((u, v), edge_label)
            if not data:
                continue

            adj_fkeys = self.edge_faces(u, v)
            adj_fkeys = [x for x in adj_fkeys if x]

            if len(adj_fkeys) < 2:
                fkey = adj_fkeys.pop()
                face_data = self.face_attribute(fkey, face_label)
                face_data += data
                self.face_attribute(fkey, face_label, face_data)
            else:
                fkey_a, fkey_b = adj_fkeys
                face_data_a = self.face_attribute(fkey_a, face_label)
                face_data_b = self.face_attribute(fkey_b, face_label)

                face_data_a += int(data // 2)
                face_data_b += int(data // 2 + data % 2)

                self.face_attribute(fkey_a, face_label, face_data_a)
                self.face_attribute(fkey_b, face_label, face_data_b)

    def to_rgmesh(self):
        """Convert to :class:`Rhino.Geometry.Mesh`

        Returns
        -------
        :class:`Rhino.Geometry.Mesh`
        """
        artist = MeshArtist(self)
        return artist.draw_mesh()

    def to_rgmesh_meters(self):
        raise NotImplementedError("Scaling needs fixing.")
        rg_mesh = self.to_rgmesh()
        T = self.get_scale_rgxform(0.001)
        rg_mesh.Transform(T)

        return rg_mesh
