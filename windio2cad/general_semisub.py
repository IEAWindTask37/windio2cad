import yaml
import numpy as np
import solid
from numpy.linalg import norm


class FloatingPlatform:
    def __init__(self, yaml_filename: str):
        self.yaml_filename = yaml_filename
        self.floating_platform_dict = None
        self.joints_dict = None

    @property
    def floating_platform(self):
        if self.floating_platform_dict is None:
            geometry = yaml.load(open(self.yaml_filename, "r"), yaml.FullLoader)
            self.floating_platform_dict = geometry["components"]["floating_platform"]
        return self.floating_platform_dict

    @property
    def joints(self):
        if self.joints_dict is None:
            self.joints_dict = {}
            for member in self.floating_platform["components"]["joints"]:
                self.joints_dict[member["name"]] = np.array(member["location"])
            for member in self.floating_platform["components"]["members"]:
                joint1 = self.joints_dict[member["joint1"]]
                joint2 = self.joints_dict[member["joint2"]]
                direction = joint2 - joint1
                if "axial_joints" in member:
                    for axial_joint in member["axial_joints"]:
                        grid = axial_joint["grid"]
                        axial_cartesian = joint1 + direction * grid
                        self.joints_dict[axial_joint["name"]] = axial_cartesian
        return self.joints_dict

    def members_union(self):
        members = []
        for member in self.floating_platform["components"]["members"]:
            joint1 = self.joints[member["joint1"]]
            joint2 = self.joints[member["joint2"]]
            grid = member["outer_shape"]["outer_diameter"]["grid"]
            values = member["outer_shape"]["outer_diameter"]["values"]
            members.append(self.rod(joint1, joint2, grid, values))
        return solid.union()(members)

    def rod(self, a, b, grid, values):
        direction = b - a
        height = norm(direction)

        member_shapes = []
        for i in range(len(grid) - 1):
            section_height = (height * grid[i + 1]) - (height * grid[i])
            bottom = (0.0, 0.0, height * grid[i])
            cylinder = solid.cylinder(r1=values[i], r2=values[i + 1], h=section_height)
            translation = solid.translate(bottom)(cylinder)
            member_shapes.append(translation)

        member_union = solid.union()(tuple(member_shapes))

        if direction[0] == 0 and direction[1] == 0:
            return solid.translate((a[0], a[1], min(a[2], b[2])))(member_union)
        else:
            w = direction / height
            u0 = np.cross(w, [0, 0, 1])
            u = u0 / norm(u0)
            v0 = np.cross(w, u)
            v = v0 / norm(v0)

            # The mulmatrix must be a list of lists
            multmatrix = (
                (u[0], v[0], w[0], a[0]),
                (u[1], v[1], w[1], a[1]),
                (u[2], v[2], w[2], a[2]),
                (0, 0, 0, 1)
            )
            return solid.multmatrix(m=multmatrix)(member_union)


if __name__ == "__main__":
    fp = FloatingPlatform("semisubmersible.yaml")
    print(solid.scad_render(fp.members_union()))
