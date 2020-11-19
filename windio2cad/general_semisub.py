from typing import Dict, List, Union
import argparse
from os import path
import yaml
import numpy as np

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

    def members(self):
        for member in self.floating_platform["components"]["members"]:
            name = member["name"]
            joint1 = self.joints[member["joint1"]]
            joint2 = self.joints[member["joint2"]]
            print(f"{name}: {joint1} to {joint2}")


if __name__ == "__main__":
    fp = FloatingPlatform("semisubmersible.yaml")
    print(fp.members())
