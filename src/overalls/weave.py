from collections import Sequence

import ghpythonlib.components as ghcomp
import Rhino.Geometry as rg
from overalls.utils import is_in_domain


def find_closest_faces(
    meshes, face_keys, n_face_connections=2, dist_domain=(None, None)
):
    mesh_a, mesh_b = meshes
    face_keys_a, face_keys_b = face_keys

    pt_cloud_dict = {}
    pt_cloud_list = []
    for fkey in face_keys_b:
        x, y, z = mesh_b.face_center(fkey)
        pt = rg.Point3d(x, y, z)
        pt_cloud_dict[pt] = fkey
        pt_cloud_list.append(pt)

    partners = []
    for fkey_a in face_keys_a:
        partner_pts = []

        x, y, z = mesh_a.face_center(fkey_a)
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
            partners.append((fkey_a, partner_pts))
    return partners


def find_closest_faces_neighbors(meshes, *args, **kwargs):
    original_partners = find_closest_faces(meshes, *args, **kwargs)

    partners = []
    for partner_a, partners_b in original_partners:
        nbors = set()
        for partner in partners_b:
            nbors.update(meshes[1].face_neighbors(partner))
        partners.append((partner_a, list(nbors)))
    return partners


def find_closest_face_same_mesh(
    mesh, fkeys, n_face_connections=2, dist_domain=(None, None), prefer_long=False
):
    fkeys = list(fkeys)
    pt_cloud_dict = {}
    for fkey in fkeys:
        x, y, z = mesh.face_center(fkey)
        pt = rg.Point3d(x, y, z)
        pt_cloud_dict[pt] = fkey

    partners = []

    for fkey in fkeys:
        dists = []
        x, y, z = mesh.face_center(fkey)
        pt1 = rg.Point3d(x, y, z)
        for pt2, key in pt_cloud_dict.iteritems():

            if key == fkey or key in mesh.face_neighbors(fkey):
                continue
            dist = pt1.DistanceTo(pt2)
            if is_in_domain(dist, dist_domain):
                dists.append((key, dist))

        if len(dists) < 1:  # No matches
            continue

        dists.sort(key=lambda x: x[1], reverse=prefer_long)
        keys, _ = zip(*dists)
        print(n_face_connections)
        print(keys[:n_face_connections])
        partners.append((fkey, keys[:n_face_connections]))
    return partners


# CREATE CONNECTING LINES
def line_center_center(meshes, fkeys, network, max_degrees=None, **kwargs):
    def _edges_verts_center(mesh, network, fkey, nkey, mesh_id):
        for vkey in mesh.face_vertices(fkey):
            nkey_2 = next(network.nodes_where({"mesh": mesh_id, "vkey": vkey}))
            network.add_edge(nkey, nkey_2, mesh=mesh_id, fkey=fkey)

    if not isinstance(meshes, Sequence):
        meshes = [meshes] * 2

    mesh_a, mesh_b = meshes
    if mesh_a == mesh_b:
        mesh_a_id = 0
        mesh_b_id = mesh_a_id
    else:
        mesh_a_id = 0
        mesh_b_id = 1

    fkey_a, fkeys_b = fkeys

    # connect center to center
    x, y, z = mesh_a.face_center(fkey_a)
    center_node_a = network.add_node(x=x, y=y, z=z, mesh=mesh_a_id, fkey=fkey_a)

    nkeys_b = []
    for fkey in fkeys_b:
        try:
            nkey = next(network.nodes_where({"mesh": mesh_b_id, "fkey": fkey}))
            if network.degrees(nkey) > max_degrees:
                continue
        except StopIteration:
            x, y, z = mesh_b.face_center(fkey)
            nkey = network.add_node(x=x, y=y, z=z, mesh=mesh_b_id, fkey=fkey)
        nkeys_b.append(nkey)

    if len(nkeys_b) == 0:  # skip rest if no matches where found
        return

    # subdivide a face so there's a node to connect to at center of face
    _edges_verts_center(mesh_a, network, fkey_a, center_node_a, mesh_a_id)

    for nkey in nkeys_b:
        # subdivide face so there's a node to connect to at center of faces
        _edges_verts_center(mesh_b, network, fkey, nkey, mesh_b_id)

        # connect a to b
        network.add_edge(center_node_a, nkey)


def line_vert_vert(
    meshes, fkeys, network, rel_vkeys=[None, None], max_degrees=None, **kwargs
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
    if not isinstance(meshes, Sequence):
        meshes = [meshes] * 2

    mesh_a, mesh_b = meshes
    if mesh_a == mesh_b:
        mesh_a_id = 0
        mesh_b_id = mesh_a_id
    else:
        mesh_a_id = 0
        mesh_b_id = 1

    fkeys_a, fkeys_b = fkeys
    rel_vkeys_a, rel_vkeys_b = rel_vkeys

    len_fkeys_a = len(fkeys_a) if isinstance(fkeys_a, Sequence) else 1
    len_fkeys_b = len(fkeys_b) if isinstance(fkeys_b, Sequence) else 1
    max_iters = max(len_fkeys_a, len_fkeys_b)

    params = []
    for param in (fkeys_a, fkeys_b, rel_vkeys_a, rel_vkeys_b):
        if param is not None and not isinstance(param, Sequence):
            param = [param]
        params.append(param)

    fkeys_a, fkeys_b, rel_vkeys_a, rel_vkeys_b = params

    if fkeys_a is None or fkeys_b is None:
        raise ValueError("Neither fkeys_a nor fkeys_b can be None.")

    for j in range(max_iters):
        fkey_a = fkeys_a[j % len(fkeys_a)]
        fkey_b = fkeys_b[j % len(fkeys_b)]

        vkeys_a = mesh_a.face_vertices(fkey_a)
        vkeys_b = mesh_b.face_vertices(fkey_b)

        nkeys_a = list(
            [
                next(network.nodes_where({"mesh": mesh_a_id, "vkey": vkey}))
                for vkey in vkeys_a
            ]
        )
        nkeys_b = list(
            [
                next(network.nodes_where({"mesh": mesh_b_id, "vkey": vkey}))
                for vkey in vkeys_b
            ]
        )

        if rel_vkeys_a:
            key_idx_a = rel_vkeys_a[j % len(rel_vkeys_a)]
        else:
            nkeys_a.sort(key=lambda x: network.degree(x))
            key_idx_a = 0
        if rel_vkeys_b:
            key_idx_b = rel_vkeys_b[j % len(rel_vkeys_b)]
        else:
            nkeys_b.sort(key=lambda x: network.degree(x))
            key_idx_b = 0

        u = nkeys_a[key_idx_a]
        v = nkeys_b[key_idx_b]

        if max_degrees:
            u_ok = network.degrees(u) <= max_degrees
            v_ok = network.degrees(v) <= max_degrees
            if not u_ok or not v_ok:
                continue

        network.add_edge(u, v)
