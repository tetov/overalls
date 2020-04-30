import Rhino.Geometry as rg


def plane_normal_to_rgplane(plane_normal):
    pt, normal = plane_normal

    return rg.Plane(rg.Point3d(*pt), rg.Vector3d(*normal))


def rgplane_to_plane_normal(rgplane):
    pt = rgplane.Origin
    normal = rgplane.Normal
    return ((pt.X, pt.Y, pt.Z), (normal.X, normal.Y, normal.Z))


class BundleProfile(object):
    def __init__(
        self,
        identifier="standard",
        origin_plane=((0, 0, 0), (0, 0, 1)),
        subprofile_planes=(((0, 0, 0), (0, 0, 1))),
        radii=10,
    ):
        self.identifier = identifier
        self.origin_plane = origin_plane
        self.subprofile_planes = subprofile_planes
        self.radii = radii

    @property
    def n_subprofiles(self):
        return len(self.subprofile_planes)

    def rgxform_to_worldxy(self):
        rgplane_from = plane_normal_to_rgplane(self.origin_plane)
        rgplane_to = rg.Plane.WorldXY

        return rg.Transform.PlaneToPlane(rgplane_from, rgplane_to)

    def as_rgcircles(self):
        circles = []
        for plane_normal in self.subprofile_planes:
            plane = plane_normal_to_rgplane(plane_normal)
            circle = rg.Circle(plane, self.radii)
            circle.Transform(self.rgxform_to_worldxy())
            circles.append(circle)
        return circles

    @classmethod
    def from_rgplanes(cls, identifier, origin_rgplane, subprofile_rgplanes, radii):
        kwargs = {}
        if origin_rgplane:
            kwargs["origin_plane"] = rgplane_to_plane_normal(origin_rgplane)

        if len(subprofile_rgplanes) > 0:
            kwargs["subprofile_planes"] = []
            for rgplane in subprofile_rgplanes:
                plane = rgplane_to_plane_normal(rgplane)
                kwargs["subprofile_planes"].append(plane)

        if identifier:
            kwargs["identifier"] = identifier

        if radii:
            kwargs["radii"] = radii

        return cls(**kwargs)


class Member(object):
    def __init__(self, start, end, bundle_profile=None):
        self.start = start
        self.end = end
        if bundle_profile:
            self.bundle_profile = bundle_profile
        else:
            self.bundle_profile = BundleProfile()

    @property
    def start_rg(self):
        return rg.Point3d(*self.start)

    @property
    def end_rg(self):
        return rg.Point3d(*self.end)

    def get_length(self):
        return self.start_rg.DistanceTo(self.end_rg)

    def _get_subprofile_xform(self):
        profile_origin = rg.Plane(
            plane_normal_to_rgplane(self.bundle_profile.origin_plane)
        )
        member_plane = rg.Plane(self.start_rg, self.as_rgvector())

        return rg.Transform.PlaneToPlane(profile_origin, member_plane)

    def as_rgline(self):
        return rg.Line(self.start_rg, self.end_rg)

    def as_rgvector(self):
        return rg.Vector3d(self.end_rg - self.start_rg)

    def get_bundle_rglines(self):
        T = self._get_subprofile_xform()

        lines = []
        for plane_normal in self.bundle_profile.subprofile_planes:
            plane = plane_normal_to_rgplane(plane_normal)
            plane.Transform(T)
            lines.append(rg.Line(plane.Origin, plane.Normal, self.get_length()))
        return lines

    def as_breps(self, tol=1e-3, angle_tol=0.015):
        args = (
            self.bundle_profile.radii,
            True,
            rg.PipeCapMode.Flat,
            True,
            tol,
            angle_tol,
        )

        breps = []
        for line in self.get_bundle_rglines():
            crv = line.ToNurbsCurve()
            breps.append(rg.Brep.CreatePipe(crv, *args))

    def as_meshes(self, segments=6, tol=1e-3):
        args = (
            self.bundle_profile.radii,
            segments,
            tol,
            rg.MeshPipeCapStyle.Flat,
            False,
        )
        meshes = []
        for line in self.get_bundle_rglines():
            crv = line.ToNurbsCurve()
            meshes.append(rg.Mesh.CreateFromCurvePipe(crv, *args))
        return meshes

    @classmethod
    def from_line(cls, line, bundle_profile=None):
        if isinstance(line, rg.Line):
            start_pt = line.From
            end_pt = line.To

            start_coord = [start_pt.X, start_pt.Y, start_pt.Z]
            end_coord = [end_pt.X, end_pt.Y, end_pt.Z]
        else:
            raise NotImplementedError(
                "Conversion from {} not implemented.".format(type(line))
            )
        return cls(start_coord, end_coord, bundle_profile=bundle_profile)
