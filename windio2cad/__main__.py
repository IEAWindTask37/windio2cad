from typing import Dict, List, Any, Optional
import argparse
import yaml
import numpy as np
from scipy.interpolate import PchipInterpolator as spline
import geometry_tools as geom
import solid
import subprocess
from numpy.linalg import norm
from math import sin, cos


# TODO: Attempt to use ruamel


class Blade:
    """
    This class renders one blade for the rotor.
    """

    def __init__(self, yaml_filename: str):
        """
        The constructor opens the YAML file and extracts the blade
        and airfoil information into instance attributes.

        Parameters
        ----------
        yaml_filename: str
            Filename that contains the geometry for the rotor.
        """
        geometry = yaml.load(open(yaml_filename, "r"), yaml.FullLoader)
        self.outer_shape = geometry["components"]["blade"]["outer_shape_bem"]
        self.airfoils = geometry["airfoils"]

    @staticmethod
    def myinterp(xi, x, f) -> np.array:
        myspline = spline(x, f)
        return myspline(xi)

    def generate_lofted(self, n_span_min=100, n_xy=400) -> np.array:
        """
        Creates the lofted shape of a blade and returns a NumPy array
        of the polygons at each cross section.

        Parameters
        ----------
        n_span_min: int
            Number of cross sections to create across span of
            blade.

        n_xy: int
            The number of x, y points in the polygons at each slice of
            the blade.

        Returns
        -------
        np.array
            An array of the polygons at each cross section of the blade.
        """
        # Use yaml grid points and others that we add
        r_span = np.unique(
            np.r_[
                np.linspace(0.0, 1.0, n_span_min),
                self.outer_shape["chord"]["grid"],
                self.outer_shape["twist"]["grid"],
                self.outer_shape["pitch_axis"]["grid"],
                self.outer_shape["reference_axis"]["x"]["grid"],
                self.outer_shape["reference_axis"]["y"]["grid"],
                self.outer_shape["reference_axis"]["z"]["grid"],
            ]
        )
        n_span = len(r_span)

        # Read in blade spanwise geometry values and put on common grid
        chord = self.myinterp(
            r_span,
            self.outer_shape["chord"]["grid"],
            self.outer_shape["chord"]["values"],
        )
        twist = self.myinterp(
            r_span,
            self.outer_shape["twist"]["grid"],
            self.outer_shape["twist"]["values"],
        )
        pitch_axis = self.myinterp(
            r_span,
            self.outer_shape["pitch_axis"]["grid"],
            self.outer_shape["pitch_axis"]["values"],
        )
        ref_axis = np.c_[
            self.myinterp(
                r_span,
                self.outer_shape["reference_axis"]["x"]["grid"],
                self.outer_shape["reference_axis"]["x"]["values"],
            ),
            self.myinterp(
                r_span,
                self.outer_shape["reference_axis"]["y"]["grid"],
                self.outer_shape["reference_axis"]["y"]["values"],
            ),
            self.myinterp(
                r_span,
                self.outer_shape["reference_axis"]["z"]["grid"],
                self.outer_shape["reference_axis"]["z"]["values"],
            ),
        ]

        # Get airfoil names and thicknesses
        af_position = self.outer_shape["airfoil_position"]["grid"]
        af_used = self.outer_shape["airfoil_position"]["labels"]
        n_af_span = len(af_position)
        n_af = len(self.airfoils)
        name = n_af * [""]
        r_thick = np.zeros(n_af)
        for i in range(n_af):
            name[i] = self.airfoils[i]["name"]
            r_thick[i] = self.airfoils[i]["relative_thickness"]

        # Create common airfoil coordinates grid
        coord_xy = np.zeros((n_af, n_xy, 2))
        for i in range(n_af):
            points = np.c_[
                self.airfoils[i]["coordinates"]["x"],
                self.airfoils[i]["coordinates"]["y"],
            ]

            # Check that airfoil points are declared from the TE suction side to TE pressure side
            idx_le = np.argmin(points[:, 0])
            if np.mean(points[:idx_le, 1]) > 0.0:
                points = np.flip(points, axis=0)

            # Remap points using class AirfoilShape
            af = geom.AirfoilShape(points=points)
            af.redistribute(n_xy, even=False, dLE=True)
            af_points = af.points

            # Add trailing edge point if not defined
            if [1, 0] not in af_points.tolist():
                af_points[:, 0] -= af_points[np.argmin(af_points[:, 0]), 0]
            c = max(af_points[:, 0]) - min(af_points[:, 0])
            af_points[:, :] /= c

            coord_xy[i, :, :] = af_points

        # Reconstruct the blade relative thickness along span with a pchip
        r_thick_used = np.zeros(n_af_span)
        coord_xy_used = np.zeros((n_af_span, n_xy, 2))
        coord_xy_interp = np.zeros((n_span, n_xy, 2))
        coord_xy_dim = np.zeros((n_span, n_xy, 2))

        for i in range(n_af_span):
            for j in range(n_af):
                if af_used[i] == name[j]:
                    r_thick_used[i] = r_thick[j]
                    coord_xy_used[i, :, :] = coord_xy[j, :, :]

        r_thick_interp = self.myinterp(r_span, af_position, r_thick_used)

        # Spanwise interpolation of the profile coordinates with a pchip
        r_thick_unique, indices = np.unique(r_thick_used, return_index=True)
        coord_xy_interp = np.flip(
            self.myinterp(
                np.flip(r_thick_interp), r_thick_unique, coord_xy_used[indices, :, :]
            ),
            axis=0,
        )
        for i in range(n_span):
            # Correction to move the leading edge (min x point) to (0,0)
            af_le = coord_xy_interp[i, np.argmin(coord_xy_interp[i, :, 0]), :]
            coord_xy_interp[i, :, 0] -= af_le[0]
            coord_xy_interp[i, :, 1] -= af_le[1]
            c = max(coord_xy_interp[i, :, 0]) - min(coord_xy_interp[i, :, 0])
            coord_xy_interp[i, :, :] /= c
            # If the rel thickness is smaller than 0.4 apply a trailing ege smoothing step
            if r_thick_interp[i] < 0.4:
                coord_xy_interp[i, :, :] = geom.trailing_edge_smoothing(
                    coord_xy_interp[i, :, :]
                )

        # Offset by pitch axis and scale for chord
        coord_xy_dim = coord_xy_interp.copy()
        coord_xy_dim[:, :, 0] -= pitch_axis[:, np.newaxis]
        coord_xy_dim = coord_xy_dim * chord[:, np.newaxis, np.newaxis]

        # Rotate to twist angle
        coord_xy_dim_twisted = np.zeros(coord_xy_interp.shape)
        for i in range(n_span):
            x = coord_xy_dim[i, :, 0]
            y = coord_xy_dim[i, :, 1]
            coord_xy_dim_twisted[i, :, 0] = x * np.cos(twist[i]) - y * np.sin(twist[i])
            coord_xy_dim_twisted[i, :, 1] = y * np.cos(twist[i]) + x * np.sin(twist[i])

        # Assemble lofted shape along reference axis
        lofted_shape = np.zeros((n_span, n_xy, 3))
        for i in range(n_span):
            for j in range(n_xy):
                lofted_shape[i, j, :] = (
                    np.r_[
                        coord_xy_dim_twisted[i, j, 1],
                        coord_xy_dim_twisted[i, j, 0],
                        0.0,
                    ]
                    + ref_axis[i, :]
                )

        return lofted_shape

    def blade_hull(self, downsample_z: int = 1) -> solid.OpenSCADObject:
        """
        This creates an OpenSCAD hull object around cross sections of a blade,
        thereby rendering the complete geometry for a single blade.

        Parameters
        ----------
        downsample_z: int
            Skips to every nth sample across the z axis of the blade. For
            example, 10 uses only every tenth cross section.

        Returns
        -------
        solid.OpenSCADObject
            The OpenSCAD object that is ready to render to code.
        """

        # Get the lofted shape and the number of sections across its span
        lofted_shape = self.generate_lofted()
        n_span = lofted_shape.shape[0]

        # Find the distance between each cross section. Find the minimum of
        # these distances and multiply by 0.1. This will be the height of each
        # extrusion for each cross section.

        diff_z = []
        for k in range(n_span - 1):
            diff_z.append(lofted_shape[k + 1, 0, 2] - lofted_shape[k, 0, 2])
        dz = 0.1 * min(diff_z)

        # Make the range to sample the span of the blade. If downsample_z
        # is 1, that means every cross section will be plotted. If it is
        # greater than 1, samples will be skipped. This is reflected in
        # the range to sample the span.

        if downsample_z == 1:
            n_span_range = range(n_span)
        else:
            n_span_range = range(0, n_span, downsample_z)

        # Create one extrusion per cross section.
        extrusions = []
        for k in n_span_range:
            bottom = lofted_shape[k, 0, 2]
            points = tuple((row[0], row[1]) for row in lofted_shape[k, :, :])
            polygon = solid.polygon(points)
            extrusion = solid.linear_extrude(dz)(polygon)
            translated_extrusion = solid.translate((0.0, 0.0, bottom))(extrusion)
            extrusions.append(translated_extrusion)

        # Create a hull around all the cross sections and return it.
        hull_of_extrusions = solid.hull()(extrusions)
        return hull_of_extrusions


class RNA:
    """
    This class generates a rotor nacelle assembly (RNA)
    """

    def __init__(self, yaml_filename: str):
        """
        This loads the YAML file and extracts the nacelle, tower, and
        hub geometry information and puts those data in instance
        attributes.

        Parameters
        ----------
        yaml_filename: str
            The absolute path to the file containing the YAML
        """
        geometry = yaml.load(open(yaml_filename, "r"), yaml.FullLoader)
        self.nacelle_dict = geometry["components"]["nacelle"]
        self.tower_dict = geometry["components"]["tower"]
        self.hub_dict = geometry["components"]["hub"]

    def rna_union(
        self, blade_object: Optional[solid.OpenSCADObject] = None
    ) -> solid.OpenSCADObject:
        """
        This creates a union of the components of the RNA. In this method,
        the hub and nacelle are created. Optionally, a single blade can
        be passed as an object. If a blade is present in the arguments, then
        three copies are assembled into a rotor.

        Parameters
        ----------
        blade_object: Optional[solid.OpenSCADObject]
            Optional. If specified, this is a single blade, three copies
            of which are assembled into a rotor.

        Returns
        -------
        solid.OpenSCADObject
            RNA union object ready for rendering to OpenSCAD.
        """

        # Get tower height and nacelle dimensions.
        tower_height = self.tower_dict["outer_shape_bem"]["reference_axis"]["z"][
            "values"
        ][-1]
        nacelle_length = 2.0 * self.nacelle_dict["drivetrain"]["overhang"]
        nacelle_height = 2.2 * self.nacelle_dict["drivetrain"]["distance_tt_hub"]
        nacelle_width = nacelle_height
        nacelle_z_height = 0.5 * nacelle_height + tower_height

        # The "cube" is the nacelle
        cube = solid.cube(
            size=[nacelle_length, nacelle_width, nacelle_height], center=True
        )

        # The hub is a sphere
        hub_center_y = 0.0
        hub_center_x = 1.5 * self.nacelle_dict["drivetrain"]["overhang"]
        hub_radius = self.hub_dict["diameter"]
        hub = solid.translate((hub_center_x, hub_center_y, 0.0))(
            solid.sphere(hub_radius)
        )

        # If blades are present, create a rotor; otherwise just return the
        # hub attached to the nacelle

        if blade_object is None:
            union = solid.union()((hub, cube))
        else:
            translated_blades = [
                solid.translate((0.0, 0.0, hub_radius))(blade_object)
                for _ in range(3)
            ]

            rotor = solid.union()([
                solid.rotate((theta, 0.0, 0.0))(blade)
                for blade, theta in zip(translated_blades,[0.0, 120.0, -120.0])
            ])

            translate_rotor = solid.translate((hub_center_x, hub_center_y, 0.0))(rotor)
            union = solid.union()((hub, cube, translate_rotor))
        rna = solid.translate((0.0, 0.0, nacelle_z_height))(union)
        return rna


class Tower:
    """
    This class generates OpenSCAD code for a tower specified in a YAML file.
    """

    def __init__(self, yaml_filename: str):
        """
        This reads the tower YAML file and extracts the data needed
        for the geometry of the tower into instance attributes.

        The tower is assumed to extend along the z axis, and the tower
        height is taken to be the height of the last element on the
        outer_diameter list.

        Parameters
        ----------
        yaml_filename: str
            The name of the YAML file containing the geometry specification
        """
        geometry = yaml.load(open(yaml_filename, "r"), yaml.FullLoader)
        self.height = geometry["components"]["tower"]["outer_shape_bem"][
            "reference_axis"
        ]["z"]["values"][-1]
        self.grid = geometry["components"]["tower"]["outer_shape_bem"][
            "outer_diameter"
        ]["grid"]
        self.values = geometry["components"]["tower"]["outer_shape_bem"][
            "outer_diameter"
        ]["values"]

    def tower_union(self) -> solid.OpenSCADObject:
        """
        This method creates the union of conical cylinders for a tower.

        Returns
        -------
        solid.OpenSCADObject
            The union of all the tower sections.
        """
        sections = []
        for i in range(1, len(self.grid)):
            bottom = self.grid[i - 1] * self.height
            section_height = (self.grid[i] - self.grid[i - 1]) * self.height
            r1 = self.values[i - 1] / 2.0
            r2 = self.values[i] / 2.0
            cylinder = solid.cylinder(r1=r1, r2=r2, h=section_height)
            translation = solid.translate((0.0, 0.0, bottom))(cylinder)
            sections.append(translation)
        section_union = solid.union()(sections)
        return section_union


class FloatingPlatform:
    """
    This class generates OpenSCAD code for an arbitrary set of floating
    members and joints in an ontology file.
    """

    def __init__(self, yaml_filename: str):
        """
        This takes the YAML input name and sets up some instance
        variables. floating_platform_dict holds the dictionary
        for the floating platform ontology portion of the YAML
        file after the YAML file is ready. joints_dict holds
        the cartesian coordinates of all the joints after the
        YAML file is read in.

        Parameters
        ----------
        yaml_filename: str
            The absolute path to the ontology YAML file as a string.
        """
        self.yaml_filename = yaml_filename
        geometry = yaml.load(open(self.yaml_filename, "r"), yaml.FullLoader)
        self.floating_platform = geometry["components"]["floating_platform"]

    def joints(self):
        """
        Converts all normal and axial joints into cartesian coordinates.

        This method uses memoization. The first time it calculates all the
        cartesian coordinates it caches these conversions in self.joints_dict.
        Subsequent calls return this cached value.

        Returns
        -------
        Dict[str, np.array]
            The names of joints mapped to their cartesian coordinates
        """
        joints_dict = {}
        for member in self.floating_platform["joints"]:
            if "cylindrical" in member and member["cylindrical"]:
                r = member["location"][0]
                theta = member["location"][1]
                x = r * cos(theta)
                y = r * sin(theta)
                z = member["location"][2]
                joints_dict[member["name"]] = np.array([x, y, z])
            else:
                joints_dict[member["name"]] = np.array(member["location"])
        for member in self.floating_platform["members"]:
            joint1 = joints_dict[member["joint1"]]
            joint2 = joints_dict[member["joint2"]]
            direction = joint2 - joint1
            if "axial_joints" in member:
                for axial_joint in member["axial_joints"]:
                    grid = axial_joint["grid"]
                    axial_cartesian = joint1 + direction * grid
                    joints_dict[axial_joint["name"]] = axial_cartesian
        return joints_dict

    def members_union(self) -> solid.OpenSCADObject:
        """
        Creates a union of all the members on the floating platform.

        Returns
        -------
        solid.OpenSCADObject
            Returns and OpenSCAD object that is the union of all members.
        """
        members = []
        joints_dict = self.joints()
        for member in self.floating_platform["members"]:
            joint1 = joints_dict[member["joint1"]]
            joint2 = joints_dict[member["joint2"]]
            grid = member["outer_shape"]["outer_diameter"]["grid"]
            values = member["outer_shape"]["outer_diameter"]["values"]
            members.append(self.member(joint1, joint2, grid, values))
        return solid.union()(members)

    def member(
        self, joint1: np.array, joint2: np.array, grid: List[float], values: List[float]
    ) -> solid.OpenSCADObject:
        """
        Creates a member between two points.

        See more information at:
        http://forum.openscad.org/Rods-between-3D-points-td13104.html

        Parameters
        ----------
        joint1: np.array
            Cartesian coordinates of the first joint.

        joint2: np.array
            Cartesian coordinates of second joint

        grid: List[float]
            Grid of axial positions on the member, as a fraction of the
            length of the member.

        values: List[float]
            Diameters at each position on the grid.

        Returns
        -------
        solid.OpenSCADObject
            Returns the transformed member ready to put into the union.
        """
        direction = joint2 - joint1
        height = norm(direction)

        member_shapes = []
        for i in range(len(grid) - 1):
            section_height = (height * grid[i + 1]) - (height * grid[i])
            bottom = (0.0, 0.0, height * grid[i])
            r1 = values[i] / 2.0
            r2 = values[i + 1] / 2.0
            cylinder = solid.cylinder(r1=r1, r2=r2, h=section_height)
            translation = solid.translate(bottom)(cylinder)
            member_shapes.append(translation)

        member_union = solid.union()(member_shapes)

        if direction[0] == 0 and direction[1] == 0:
            return solid.translate((joint1[0], joint1[1], min(joint1[2], joint2[2])))(
                member_union
            )
        else:
            w = direction / height
            u0 = np.cross(w, [0, 0, 1])
            u = u0 / norm(u0)
            v0 = np.cross(w, u)
            v = v0 / norm(v0)

            # The mulmatrix must be a list of lists
            multmatrix = (
                (u[0], v[0], w[0], joint1[0]),
                (u[1], v[1], w[1], joint1[1]),
                (u[2], v[2], w[2], joint1[2]),
                (0, 0, 0, 1),
            )
            return solid.multmatrix(m=multmatrix)(member_union)


if __name__ == "__main__":

    # Create a command line parser
    parser = argparse.ArgumentParser(
        description="Translate a yaml definition of a semisubmersible platform into an OpenSCAD source file."
    )
    parser.add_argument("--input", help="Input .yaml file", required=True)
    parser.add_argument(
        "--output",
        help="Output .stl file. If the file exists it will be overwritten.",
        required=True,
    )
    parser.add_argument("--openscad", help="Path to OpenSCAD executable", required=True)
    parser.add_argument("--downsample", default=1, type=int, help="Defaults to 1, meaning every cross section is rendered.")
    args = parser.parse_args()

    intermediate_openscad = "intermediate.scad"

    print(f"Input yaml: {args.input}")
    print(f"Output .stl: {args.output}")
    print(f"Blade downsampling: {args.downsample}")
    print(f"Intermediate OpenSCAD: {intermediate_openscad}")
    print(f"Path to OpenSCAD: {args.openscad}")
    print("Parsing .yaml ...")

    blade = Blade(args.input)
    blade_object = blade.blade_hull(downsample_z=args.downsample)
    fp = FloatingPlatform(args.input)
    tower = Tower(args.input)
    rna = RNA(args.input)

    with open(intermediate_openscad, "w") as f:
        f.write("$fn = 25;\n")
        big_union = solid.union()(
            [fp.members_union(), tower.tower_union(), rna.rna_union(blade_object)]
        )
        f.write(solid.scad_render(big_union))

    print("Creating .stl ...")
    subprocess.run([args.openscad, "-o", args.output, intermediate_openscad])
    print("Done!")
