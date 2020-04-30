from itertools import cycle

from compas.datastructures import Mesh
from compas.geometry import Point


class OverallMesh(Mesh):
    def __init__(self):
        super(OverallMesh, self).__init__()

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
                    returned_fkeys = subdiv_func(fkey_prim, **kwargs)
                    next_to_subd += returned_fkeys
                    new_fkeys += returned_fkeys
                i += 1

        return new_fkeys

    def face_subdiv_retriangulate(self, fkey, **kwargs):
        _, fkeys = self.insert_vertex(fkey, return_fkeys=True)
        return fkeys

    def face_subdiv_frame(self, fkey, rel_dist=0.5):
        x, y, z = self.face_center(fkey)
        face_center_pt = Point(x, y, z)

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

        return new_fkeys
