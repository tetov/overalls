import Rhino.Geometry as rg
from compas.datastructures import Mesh
from compas.datastructures import Network


def rglines_from_edges(mesh, nkeys=None):
    if nkeys is None:
        nkeys = mesh.edges()
    lines = []
    for u, v in nkeys:
        pt_a, pt_b = mesh.edge_coordinates(u, v)
        pt_a = rg.Point3d(*pt_a)
        pt_b = rg.Point3d(*pt_b)
        line = rg.Line(pt_a, pt_b)
        lines.append(line)
    return lines


def to_compas_bugfix(mesh, cls=None):
    if not cls:
        cls = Mesh

    faces = []
    for face in mesh.faces:
        if face[0] == face[-1]:
            faces.append(face[:-1])
        elif face[-2] == face[-1]:
            faces.append(face[:-1])
        else:
            faces.append(face)
    return cls.from_vertices_and_faces(mesh.vertices, faces)


def meshes_to_network(meshes):
    network = Network()

    network.update_default_node_attributes(mesh=None, vkey=None, fkey=None)
    network.update_default_edge_attributes(mesh=None, fkey=None)

    for i, mesh in enumerate(meshes):
        for vkey in mesh.vertices():
            x, y, z = mesh.vertex_coordinates(vkey)
            network.add_node(x=x, y=y, z=z, mesh=i, vkey=vkey)
        for u, v in mesh.edges():
            u1 = next(network.nodes_where({"vkey": u, "mesh": i}))
            v1 = next(network.nodes_where({"vkey": v, "mesh": i}))
            network.add_edge(u1, v1, mesh=i)

    return network


def is_in_domain(value, domain):
    min, max = domain

    min_ok = value > min or min is None
    max_ok = value < max or max is None

    return min_ok and max_ok


def flatten(list_):
    return [item for sublist in list_ for item in sublist]
