import math
from itertools import combinations
from itertools import ifilter

import Rhino.Geometry as rg
from compas.geometry import Scale
from compas.geometry import Translation
from compas.utilities import geometric_key
from ghpythonlib.treehelpers import list_to_tree
from overalls.utils import is_in_domain
from overalls.utils import to_list
from System.Collections.Generic import KeyNotFoundException


class SpaceFrameMixin(object):
    def __init__(self):
        self.dist_dict = {}

    def nv_degree(self, key, *args, **kwargs):
        raise NotImplementedError

    def nvs(self, *args, **kwargs):
        raise NotImplementedError

    def nv_neighbors(self, key, *args, **kwargs):
        raise NotImplementedError

    def nv_coordinates(self, key, *args, **kwargs):
        raise NotImplementedError

    def nv_attribute(self, u, v, key, *args, **kwargs):
        raise NotImplementedError

    def faces(self, *args, **kwargs):
        raise NotImplementedError

    def derive_edge_attr(self, u, v, key):
        u_val = self.nv_attribute(u, key)
        v_val = self.nv_attribute(v, key)

        # remove None values
        vals = list(ifilter(None, [u_val, v_val]))

        if len(vals) == 0:
            return
        elif len(vals) == 2 and vals[0] != vals[1]:
            return

        val = vals[0]

        self.edge_attribute((u, v), key, val)
        return val

    def derive_edges_attr(self, key, ekeys=None):
        if not ekeys:
            ekeys = self.edges()

        for u, v in ekeys:
            self.derive_edge_attr(u, v, key)

    def get_bbox_origin(self):
        pts = [self.nv_coordinates(key) for key in self.nvs()]

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

    def ensure_dist_dict(self):
        if set(self.nvs()) != set(self.dist_dict.keys()):
            self.build_dist_dict()

    def build_dist_dict(self):
        for nkey in self.nvs():
            self.dist_dict[nkey] = {}

        key_combinations = combinations(self.nvs(), 2)

        for u, v in key_combinations:
            x, y, z = self.nv_coordinates(u)
            x1, y1, z1 = self.nv_coordinates(v)

            dist = math.sqrt((x - x1) ** 2 + (y - y1) ** 2 + (z - z1) ** 2)

            self.dist_dict[u][v] = dist
            self.dist_dict[v][u] = dist

    def connected_edges(self, key):
        u = key

        for v in self.nv_neighbors(u):
            yield (u, v)

    def ok_degree(self, key, max_degree):
        if not max_degree:
            return True
        return self.nv_degree(key) <= max_degree

    def ok_degrees(self, keys, max_degree):
        keys = to_list(keys)
        for key in keys:
            if not self.ok_degree(key, max_degree):
                return False
        return True

    def ok_edges_angles(self, keys, min_angle, **kwargs):
        vkeys = to_list(keys)
        for vkey in vkeys:
            if not self.ok_edge_angles(vkey, min_angle, **kwargs):
                return False
        return True

    def ok_edge_angles(self, key, min_angle, additional_edge=None):
        if not min_angle:
            return True
        smallest_edge_angle = self.nv_smallest_edge_angle(
            key, additional_edge=additional_edge
        )
        return smallest_edge_angle > min_angle

    def nv_smallest_edge_angle(self, vkey, additional_edge=None):
        edges = list(self.connected_edges(vkey))

        if additional_edge:
            edges.append(additional_edge)

        vectors = []
        for u, v in edges:
            if u != vkey:
                u, v = v, u
            vec = rg.Vector3d(*self.edge_vector(u, v))
            vec.Unitize()
            vectors.append(vec)

        to_compare = combinations(range(len(vectors)), 2)
        vector_angles = []
        for v, v1 in list(to_compare):
            v, v1 = vectors[v], vectors[v1]
            vector_angles.append(rg.Vector3d.VectorAngle(v, v1))

        return min(vector_angles)

    def edge_length_ratios(self, fkey):
        edges = self.face_halfedges(fkey)
        edge_lengths = []

        for u, v in edges:
            edge_lengths.append([(u, v), self.edge_length(u, v)])

        edge_lengths.sort(key=lambda x: x[1])
        _, ref_length = edge_lengths[0]

        return [(ekey, length / ref_length) for ekey, length in edge_lengths]

    def ok_edge_length_ratios(self, fkey, max_ratio_diff):
        if not max_ratio_diff:
            return True

        edge_length_ratios = self.edge_length_ratios(fkey)
        _, longest_edge_ratio = edge_length_ratios[-1]

        return longest_edge_ratio <= max_ratio_diff

    def ok_edge_lengths(self, ekeys, length_domain):
        for u, v in ekeys:
            if not self.ok_edge_length(u, v, length_domain):
                return False
        return True

    def ok_edge_length(self, u, v, length_domain):
        if not length_domain:
            return True
        return is_in_domain(self.edge_length(u, v), length_domain)

    def analyse(
        self,
        max_degree=None,
        min_angle=None,
        max_edge_ratio_diff=None,
        edge_length_domain=None,
    ):
        wrong_degree = []
        wrong_angle = []
        wrong_edge_length = []
        wrong_ratio = []

        for key in self.nvs():
            if not self.ok_degree(key, max_degree):
                wrong_degree.append(key)
            if not self.ok_edge_angles(key, min_angle):
                wrong_angle.append(key)

        if hasattr(self, "faces"):
            for fkey in self.faces():
                if not self.ok_edge_length_ratios(fkey, max_edge_ratio_diff):
                    wrong_ratio.append(fkey)

        for u, v in self.edges():
            if not self.ok_edge_length(u, v, edge_length_domain):
                wrong_edge_length.append((u, v))

        return wrong_degree, wrong_angle, wrong_edge_length, wrong_ratio

    def analyse_rg(
        self,
        show_only_bad=False,
        max_degree=None,
        min_angle=None,
        edge_length_domain=None,
        max_edge_ratio_diff=None,
    ):
        """Returns Rhino.Geometry objects to visualize component properties.

        Parameters
        ----------
        show_only_bad : :obj:`bool`, optional
            Toggles between showing only components not meeting requirements
            or not. Defaults to ``False``.
        max_degree : :obj:`int`, optional
            Maximum vertex degree, defaults to ``None``.
        min_angle : :obj:`float`, optional
            Sets lower bound of allowed edge angles at any vertex. Uses radians.
            Defaults to ``None``.
        edge_length_domain : :obj:`tuple` of :obj:`float`, optional
            Lower and upper bounds of allowed edge lengths.

        Returns
        -------
        :obj:`list` of :class:`Rhino.Geometry.Point3d`
            Points representing nodes.
        :obj:`list` of :class:`Rhino.Geometry.Line`
            Lines representing edges.
        :obj:`list` of :obj:`tuple` of :obj:`Rhino.Geometry.Point3d` and :obj:`str`
            Pairs of marker placement points and marker text.
        """
        if show_only_bad:
            wrong_degree, wrong_angle, wrong_edge_length, wrong_ratio = self.analyse(
                max_degree=max_degree,
                min_angle=min_angle,
                edge_length_domain=edge_length_domain,
                max_edge_ratio_diff=max_edge_ratio_diff,
            )
            nvs = set((wrong_degree + wrong_angle))
            edges = wrong_edge_length
            faces = wrong_ratio
        else:
            nvs = self.nvs()
            edges = self.edges()
            faces = self.faces()

        nv_geos = []
        edge_geos = []
        # face_geos = []

        markers = []

        for key in nvs:
            x, y, z = self.nv_coordinates(key)
            txt_loc = [x, y, z]
            geo = rg.Point3d(x, y, z)
            txt_segments = []

            angle = math.degrees(self.nv_smallest_edge_angle(key))
            angle_txt = "A={}".format(int(angle))

            degree_txt = "D={}".format(self.nv_degree(key))

            if show_only_bad and (key in wrong_degree or key in wrong_angle):
                nv_geos.append(geo)

                if key in wrong_degree:
                    txt_segments.append(degree_txt)
                if key in wrong_angle:
                    txt_segments.append(angle_txt)
            else:
                nv_geos.append(geo)
                txt_segments = [degree_txt + angle_txt]

            if len(txt_segments) > 0:
                txt = " ".join(txt_segments)
                markers.append((txt_loc, txt))

        for u, v in edges:
            geo = self.edge_to_rgline(u, v)
            txt_loc = self.edge_midpoint(u, v)
            txt = "L={}".format(int(self.edge_length(u, v)))

            if not show_only_bad or (u, v) in wrong_edge_length:
                markers.append((txt_loc, txt))
                edge_geos.append(geo)

        for fkey in faces:
            txt_loc = self.face_center(fkey)
            edge_length_ratios = self.edge_length_ratios(fkey)
            _, smallest_longest_ratio = edge_length_ratios[-1]
            txt = "R={}".format(round(smallest_longest_ratio))

            markers.append((txt_loc, txt))

        markers = [(rg.Point3d(*loc), txt) for loc, txt in markers]

        return nv_geos, edge_geos, markers

    def correlate_lines(self, lines, data, label="correlated_data", precision="1e-3"):
        """Correlate lines to mesh edges and add attributes to correlated edges.

        Parameters
        ----------
        lines : :obj:`list` of :class:`Rhino.Geometry.Line`
            Rhino lines to correlate to mesh edges.
        data : :obj:`list`
            Data to add to edge attributes.
        label : :obj:`str`, optional
            Label for edge-data (key in `attr_dict`)
        """
        nv_dict = {}
        for key in self.nvs():
            xyz = geometric_key(self.nv_coordinates(key), precision=precision)

            nv_dict[xyz] = key

        for line, line_data in zip(lines, data):
            start_xyz = geometric_key(
                [line.FromX, line.FromY, line.FromZ], precision=precision
            )

            end_xyz = geometric_key([line.ToX, line.ToY, line.ToZ], precision=precision)

            start_nv = nv_dict.get(start_xyz)
            end_nv = nv_dict.get(end_xyz)

            try:
                self.edge_attribute((start_nv, end_nv), label, line_data)
            except KeyNotFoundException:
                exists = self.has_edge((start_nv, end_nv))
                reverse_exists = self.has_edge((end_nv, start_nv))
                print(
                    "Edge: {}-{}. Exists: {}. Reverse exists: {}".format(
                        start_nv, end_nv, exists, reverse_exists
                    )
                )

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
