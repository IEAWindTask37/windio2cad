from typing import Dict, List, Any
import argparse
import yaml
import numpy as np
import solid
import subprocess
from numpy.linalg import norm
from math import sin, cos


class Tower:
    def __init__(self, yaml_filename: str):
        self.yaml_filename = yaml_filename
        geometry = yaml.load(open(self.yaml_filename, "r"), yaml.FullLoader)
        self.height = geometry["components"]["tower"]["outer_shape_bem"]["reference_axis"]["z"]["values"][-1]
        self.grid = geometry["components"]["tower"]["outer_shape_bem"]["outer_diameter"]["grid"]
        self.values = geometry["components"]["tower"]["outer_shape_bem"]["outer_diameter"]["values"]

    def tower_union(self) -> solid.OpenSCADObject:
        sections = []
        for i in range(1, len(self.grid)):
            bottom = self.grid[i - 1] * self.height
            section_height = (self.grid[i] - self.grid[i - 1]) * self.height
            r1 = self.values[i - 1]
            r2 = self.values[i]
            cylinder = solid.cylinder(r1=r1, r2=r2, h=section_height)
            translation = solid.translate((0.0, 0.0, bottom))(cylinder)
            sections.append(translation)
        section_union = solid.union()(sections)
        return section_union


class FloatingPlatform:
    """
    This class generates OpenSCAD code for an arbitrary set of members and
    joints in an ontology YAML file.
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
        self.floating_platform_dict = None
        self.joints_dict = None

    @property
    def floating_platform(self) -> Dict[str, Any]:
        """
        Reads the floating platform portion of the YAML ontology and
        returns it as a dictionary.

        This method uses memoization to cache the results of reading the
        floating platform. For the first time the floating platform is
        read, self.floating_platform_dict is None, and the file is read
        and parsed as YAML. This parsed result is stored in
        self.floating_platform_dict. Subsequent calls return this cached
        value without reading and parsing the YAML again.

        Returns
        -------
        Dict[str, Any]
            The dictionary that describes the floating platform.
        """
        if self.floating_platform_dict is None:
            geometry = yaml.load(open(self.yaml_filename, "r"), yaml.FullLoader)
            self.floating_platform_dict = geometry["components"]["floating_platform"]
        return self.floating_platform_dict

    @property
    def joints(self) -> Dict[str, np.array]:
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
        if self.joints_dict is None:
            self.joints_dict = {}
            for member in self.floating_platform["joints"]:
                if "cylindrical" in member and member["cylindrical"]:
                    r = member["location"][0]
                    theta = member["location"][1]
                    x = r * cos(theta)
                    y = r * sin(theta)
                    z = member["location"][2]
                    self.joints_dict[member["name"]] = np.array([x, y, z])
                else:
                    self.joints_dict[member["name"]] = np.array(member["location"])
            for member in self.floating_platform["members"]:
                joint1 = self.joints_dict[member["joint1"]]
                joint2 = self.joints_dict[member["joint2"]]
                direction = joint2 - joint1
                if "axial_joints" in member:
                    for axial_joint in member["axial_joints"]:
                        grid = axial_joint["grid"]
                        axial_cartesian = joint1 + direction * grid
                        self.joints_dict[axial_joint["name"]] = axial_cartesian
        return self.joints_dict

    def members_union(self) -> solid.OpenSCADObject:
        """
        Creates a union of all the members on the floating platform.

        Returns
        -------
        solid.OpenSCADObject
            Returns and OpenSCAD object that is the union of all members.
        """
        members = []
        for member in self.floating_platform["members"]:
            joint1 = self.joints[member["joint1"]]
            joint2 = self.joints[member["joint2"]]
            grid = member["outer_shape"]["outer_diameter"]["grid"]
            values = member["outer_shape"]["outer_diameter"]["values"]
            members.append(self.member(joint1, joint2, grid, values))
        return solid.union()(members)

    def member(self, joint1: np.array, joint2: np.array, grid: List[float], values: List[float]) -> solid.OpenSCADObject:
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

        member_union = solid.union()(tuple(member_shapes))

        if direction[0] == 0 and direction[1] == 0:
            return solid.translate((joint1[0], joint1[1], min(joint1[2], joint2[2])))(member_union)
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
                (0, 0, 0, 1)
            )
            return solid.multmatrix(m=multmatrix)(member_union)


if __name__ == "__main__":
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
    args = parser.parse_args()

    intermediate_openscad = "intermediate.scad"

    print(f"Input yaml: {args.input}")
    print(f"Output .stl: {args.output}")
    print(f"Intermediate OpenSCAD: {intermediate_openscad}")
    print(f"Path to OpenSCAD: {args.openscad}")
    print("Parsing .yaml ...")

    fp = FloatingPlatform(args.input)
    tower = Tower(args.input)

    with open(intermediate_openscad, "w") as f:
        f.write("$fn = 25;\n")
        f.write(solid.scad_render(fp.members_union()))
        f.write(solid.scad_render(tower.tower_union()))

    print("Creating .stl ...")
    subprocess.run([args.openscad, "-o", args.output, intermediate_openscad])
    print("Done!")
