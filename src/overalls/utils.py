import Rhino.Geometry as rg
from compas.datastructures import Network


def plane_normal_to_rgplane(plane_normal):
    pt, normal = plane_normal

    return rg.Plane(rg.Point3d(*pt), rg.Vector3d(*normal))


def rgplane_to_plane_normal(rgplane):
    pt = rgplane.Origin
    normal = rgplane.Normal
    return ((pt.X, pt.Y, pt.Z), (normal.X, normal.Y, normal.Z))


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
    print("val: {}, domain: {}".format(value, domain))
    min_, max_ = domain

    min_ok = value > min_ or min_ is None
    max_ok = value < max_ or max_ is None

    return min_ok and max_ok


def flatten(list_):
    return [item for sublist in list_ for item in sublist]


def is_str(obj):
    try:
        basestring
    except NameError:
        basestring = (str, bytes)
    return isinstance(obj, basestring)


def to_iter(obj):
    if not is_str(obj):
        try:
            return iter(obj)
        except TypeError:
            pass
    return iter([obj])


def to_list(obj):
    if not is_str(obj):
        try:
            return list(obj)
        except TypeError:
            pass
    return [obj]
