from itertools import cycle

import ghpythonlib.components as ghcomp
import Rhino.Geometry as rg
from compas.datastructures import Mesh
from compas.geometry import Point
from compas.geometry import Vector
from overalls.utils import is_in_domain
from overalls.utils import to_list


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

    def is_valid(self):
        """This is sort of a test to see if mesh is weaved."""
        for u, v in self.edges():
            if self.edge_faces(u, v) > 2:
                return False
        return True

    def components_attributes(self, fkeys, vkeys, ekeys, attrs):
        for key, value in attrs.items():
            for vkey in vkeys:
                self.vertex_attribute(vkey, key, value)
            for u, v in ekeys:
                self.edge_attribute((u, v), key, value)
            for fkey in fkeys:
                self.face_attribute(fkey, key, value)

    def subdiv_faces(self, fkeys, n_iters, scheme=[0], **kwargs):
        schemes = [self.face_subdiv_retriangulate, self.face_subdiv_frame]

        repeating_n_iters = cycle(n_iters)
        repeating_schemes = cycle(scheme)

        new_fkeys = []
        for fkey in fkeys:
            n_iters = next(repeating_n_iters)
            subdiv_func = schemes[next(repeating_schemes)]

            i = 0
            next_to_subd = [fkey]
            while i < n_iters:
                to_subd = next_to_subd
                next_to_subd = []

                for fkey_prim in to_subd:
                    _, returned_fkeys = subdiv_func(fkey_prim, **kwargs)
                    next_to_subd += returned_fkeys
                    new_fkeys += returned_fkeys
                i += 1

        return new_fkeys

    def face_subdiv_retriangulate(self, fkey, move_z=0.0, rel_pos=None, **attrs):
        xyz = Point(*self.face_center(fkey))

        if rel_pos:
            if len(rel_pos) != len(list(self.face_vertices(fkey))):
                raise Exception(
                    "Length of rel_pos not same as amount of face vertices."
                )
            for factor, vkey in zip(rel_pos, self.face_vertices(fkey)):
                v_xyz = Point(*self.vertex_coordinates(vkey))
                v = (xyz - v_xyz) * 2
                xyz += v * factor

        if move_z:
            v = Vector(*self.face_normal(fkey))
            xyz += v * move_z

        vkey, fkeys = self.insert_vertex(fkey, xyz=list(xyz), return_fkeys=True)

        ekeys = [(vkey, v) for v in self.vertex_neighbors(vkey)]

        attrs.update({"parent_fkey": fkey})

        self.components_attributes(fkeys, [vkey], ekeys, attrs)

        return vkey, fkeys

    def face_subdiv_frame(self, fkey, move_z=False, rel_dist=0.5, **attrs):
        if rel_dist == 1:
            return self.face_subdiv_retriangulate(fkey, move_z=move_z, **attrs)
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

    def find_closest_faces(
        self, fkeys_origin, fkeys_search, n_face_connections=2, dist_domain=(None, None)
    ):
        pt_cloud_dict = {}
        pt_cloud_list = []

        for fkey in fkeys_search:
            x, y, z = self.face_center(fkey)
            pt = rg.Point3d(x, y, z)
            pt_cloud_dict[pt] = fkey
            pt_cloud_list.append(pt)

        partners = []
        for fkey in fkeys_origin:
            partner_pts = []

            x, y, z = self.face_center(fkey)
            pt = rg.Point3d(x, y, z)

            closest_pts, _, dist = ghcomp.ClosestPoints(
                pt, pt_cloud_list, n_face_connections
            )

            if n_face_connections < 2:
                closest_pts = [closest_pts]

            for cp in closest_pts:
                if is_in_domain(pt.DistanceTo(cp), dist_domain):
                    partner_pts.append(pt_cloud_dict[cp])

            if len(partner_pts) > 0:
                partners.append((fkey, partner_pts))
        return partners

    def find_closest_faces_neighbors(self, fkeys_origin, fkeys_search, *args, **kwargs):
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

    def find_closest_face_same_mesh(
        self, fkeys, n_face_connections=2, dist_domain=(None, None), prefer_long=False
    ):
        fkeys = list(fkeys)
        pt_cloud_dict = {}
        for fkey in fkeys:
            x, y, z = self.face_center(fkey)
            pt = rg.Point3d(x, y, z)
            pt_cloud_dict[pt] = fkey

        partners = []

        for fkey in fkeys:
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

    # CREATE CONNECTING LINES
    def connect_faces_center_centers(
        self, start_fkey, end_fkeys, max_degrees=None,
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
        max_degrees=None,
    ):
        """Create a line between a vertex on one mesh to a vertex on another.

        Parameters
        ----------
        meshes : :obj:`list` of :class:`compas.datastructures.Mesh`
            Mesh or meshes to connect.
        fkeys_a : list of int or int
            Face keys for faces on first mesh to connect. If length between fkeys_a
            and fkeys_b differ the shorter list will be wrapped.
        fkeys_a : list of int or int
            Face keys for faces on second mesh to connect. If length between fkeys_a
            and fkeys_b differ the shorter list will be wrapped.
        relative_vkeys_a : list of int or int, optional
            The face vertices on faces on mesh_a that will be connected.
        relative_vkeys_b : list of int or int, optional
            The face vertices on faces on mesh_b that will be connected.

        Returns
        -------
        list of Rhino.Geometry.Line
        """
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

            if max_degrees:
                u_ok = self.vertex_degree(u) <= max_degrees
                v_ok = self.vertex_degree(v) <= max_degrees
                if not u_ok or not v_ok:
                    continue

            self.add_edge(
                u,
                v,
                created_from=self.FROM_CONN_VERTICES,
                parent_fkey=(start_fkey, end_fkey),
            )
