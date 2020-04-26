from itertools import cycle

import Rhino.Geometry as rg


def subdiv(subd, fkeys, n_iters, scheme=[0], **kwargs):
    schemes = [face_subdiv_retriangulate, face_subdiv_frame]

    subd = subd.copy()

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
                # subd = subd.copy()

                subd, returned_fkeys = subdiv_func(subd, fkey_prim, **kwargs)
                next_to_subd += returned_fkeys
                new_fkeys += returned_fkeys
            i += 1

    return subd, new_fkeys


def face_subdiv_retriangulate(subd, fkey, **kwargs):
    x, y, z = subd.face_center(fkey)
    fkeys = []
    w = subd.add_vertex(x=x, y=y, z=z)
    for u, v in subd.face_halfedges(fkey):
        fkeys.append(subd.add_face([u, v, w]))
    del subd.face[fkey]
    # fix bug raised in https://github.com/compas-dev/compas/issues/522
    sets_verts = [set(subd.face_vertices(fkey_)) for fkey_ in fkeys]
    for i in range(len(sets_verts)):
        if len(sets_verts[i]) == 2:
            del subd.face[fkeys.pop(i)]
    return subd, fkeys


def face_subdiv_frame(subd, fkey, rel_dist=0.5):
    x, y, z = subd.face_center(fkey)
    face_center_pt = rg.Point3d(x, y, z)

    new_vkeys = []
    for x, y, z in subd.face_coordinates(fkey):
        pt = rg.Point3d(x, y, z)

        v = rg.Vector3d(face_center_pt) - rg.Vector3d(pt)
        pt += v * rel_dist
        new_vkeys.append(subd.add_vertex(x=pt.X, y=pt.Y, z=pt.Z))

    face_verts = subd.face_vertices(fkey)
    new_fkeys = []
    for j in range(len(face_verts)):
        vkeys = []
        vkeys.append(face_verts[j])
        vkeys.append(face_verts[(j + 1) % len(face_verts)])
        vkeys.append(new_vkeys[(j + 1) % len(new_vkeys)])
        vkeys.append(new_vkeys[j])

        new_fkeys.append(subd.add_face(vkeys))

    # add new center face
    new_fkeys.append(subd.add_face(new_vkeys))

    subd.delete_face(fkey)

    return subd, new_fkeys
