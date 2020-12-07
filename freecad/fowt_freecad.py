from typing import Dict, List, Any
import yaml
import numpy as np
from numpy.linalg import norm


class FloatingPlatform:
    def __init__(self, yaml_filename: str):
        self.doc = App.newDocument("Doc")
        self.yaml_filename = yaml_filename
        self.floating_platform_dict = None
        self.joints_dict = None
        self.shape_serial = 0

    def floating_platform(self):
        if self.floating_platform_dict is None:
            geometry = yaml.load(open(self.yaml_filename, "r"), yaml.FullLoader)
            self.floating_platform_dict = geometry["components"]["floating_platform"]
        return self.floating_platform_dict

    # def joints(self):
    #     if self.joints_dict is None:
    #         self.joints_dict = {}
    #         for member in self.floating_platform["components"]["joints"]:
    #             self.joints_dict[member["name"]] = np.array(member["location"])
    #         for member in self.floating_platform["components"]["members"]:
    #             joint1 = self.joints_dict[member["joint1"]]
    #             joint2 = self.joints_dict[member["joint2"]]
    #             direction = joint2 - joint1
    #             if "axial_joints" in member:
    #                 for axial_joint in member["axial_joints"]:
    #                     grid = axial_joint["grid"]
    #                     axial_cartesian = joint1 + direction * grid
    #                     self.joints_dict[axial_joint["name"]] = axial_cartesian
    #     return self.joints_dict
    #
    # def members_union(self):
    #     members = []
    #     for member in self.floating_platform["components"]["members"]:
    #         joint1 = self.joints[member["joint1"]]
    #         joint2 = self.joints[member["joint2"]]
    #         grid = member["outer_shape"]["outer_diameter"]["grid"]
    #         values = member["outer_shape"]["outer_diameter"]["values"]
    #         members.append(self.member(joint1, joint2, grid, values))
    #     return members
    #
    # def member(self, joint1: np.array, joint2: np.array, grid: List[float], values: List[float]):
    #     direction = joint2 - joint1
    #     height = norm(direction)
    #
    #     member_shapes = []
    #     for i in range(len(grid) - 1):
    #         section_height = (height * grid[i + 1]) - (height * grid[i])
    #         bottom = (0.0, 0.0, height * grid[i])
    #         r1 = values[i] / 2.0
    #         r2 = values[i + 1] / 2.0
    #         # cylinder = solid.cylinder(r1=r1, r2=r2, h=section_height)
    #         # translation = solid.translate(bottom)(cylinder)
    #
    #         cylinder = Part.makeCone(r1, r2, section_height)
    #         cylinder.translate((0.0, 0.0, bottom))
    #         member_shapes.append(cylinder)
    #
    #     cumulative_member_shape = member_shapes[0]
    #     for shape in member_shapes[1:]:
    #         cumulative_member_shape = cumulative_member_shape.fuse(shape)
    #
    #     if direction[0] == 0 and direction[1] == 0:
    #         return cumulative_member_shape.translate((joint1[0], joint1[1], min(joint1[2], joint2[2])))
    #     else:
    #         w = direction / height
    #         u0 = np.cross(w, [0, 0, 1])
    #         u = u0 / norm(u0)
    #         v0 = np.cross(w, u)
    #         v = v0 / norm(v0)
    #         m = FreeCAD.Base.Matrix(((u[0], v[0], w[0], joint1[0]), (u[1], v[1], w[1], joint1[1]), (u[2], v[2], w[2], joint1[2]), (0, 0, 0, 1)))
    #         return cumulative_member_shape.transformShape(m)


if __name__ == "__main__":
    fp = FloatingPlatform("semisubmersible.yaml")
    fp.members_union()
