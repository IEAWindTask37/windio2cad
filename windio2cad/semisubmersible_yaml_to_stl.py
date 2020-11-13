from typing import Dict, List, Union
import argparse
from os import path
import yaml
import subprocess


def translate_yaml_to_scad_variables(
    yaml_filename: str,
) -> Dict[str, Union[List[float], float]]:
    """
    This takes a yaml filename with a semisubmersible platform
    and covnerts it to a dictionary with key value pairs. The key
    names of this output dictionary are the variable names and
    the values are the values of those variables.

    The members are assumed to be circular.

    Parameters
    ----------
    yaml_filename: str
        The filename that contains the .yaml to be transformed into
        OpenSCAD variables.

    Returns
    -------
    Dict[str, Union[List[float], float]]
        The key/value pairs of variable names and their values
    """
    geometry = yaml.load(open(yaml_filename, "r"), yaml.FullLoader)
    floating_platform = geometry["components"]["floating_platform"]

    # This dictionary will handle all the joints, normal and axial
    joints = {}

    # Extract normal joints which have 3 dimensional vectors as locations
    for x in floating_platform["components"]["joints"]:
        joints[x["name"]] = x["location"]

    # Extract axial joints which have scalar values as locations
    for x in floating_platform["components"]["members"]:
        if "axial_joints" in x:
            for y in x["axial_joints"]:
                joints[y["name"]] = y["grid"]

    # Now parse out all members, their locations, and their outer diameters
    # For now assume circular members.
    members = {}
    for x in floating_platform["components"]["members"]:
        key = x["name"]
        value = {}
        value["joint1"] = joints[x["joint1"]]
        value["joint2"] = joints[x["joint2"]]
        # This is where I assume circular members
        value["outer_grid"] = x["outer_shape"]["outer_diameter"]["grid"]
        value["outer_values"] = x["outer_shape"]["outer_diameter"]["values"]
        members[key] = value

    # Now flatten the dictionary of dictionaries
    flattened = {}
    for outer_key, outer_value in members.items():
        for inner_key, flat_value in outer_value.items():
            flat_key = f"{outer_key}_{inner_key}"
            flattened[flat_key] = flat_value

    return flattened


def create_scad_file(
    scad_variables: Dict[str, Union[List[float], float]],
    renderer_filename: str,
    output_filename: str,
) -> None:
    """
    This takes a list variables as created by
    translate_yaml_to_scad_variables() and outputs it to a file
    that can be read by OpenSCAD. At the end it includes a
    "renderer" file, which is simply OpenSCAD code that knows
    how to render these variables.

    Parameters
    ----------
    scad_variables: Dict[str, Union[List[float], float]]:
        The variables to be rendered with OpenSCAD.

    renderer_filename: str
        The filename that contains the OpenSCAD modules that
        create each member of the platform.

    output_filename: str
        The output filename that the renderable OpenSCAD
        source code should go into.
    """
    with open(output_filename, "w") as f:
        for var_name, var_value in scad_variables.items():
            f.write(f"{var_name} = {var_value};\n")
        f.write(f"include <{renderer_filename}>;\n")
        f.write("$fn = 25;")


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
    parser.add_argument(
        "--modules", help="OpenSCAD modules that enable rendering", required=True
    )
    parser.add_argument("--openscad", help="Path to OpenSCAD executable", required=True)
    args = parser.parse_args()

    intermediate_openscad = "intermediate.scad"

    print(f"Input yaml: {args.input}")
    print(f"Output .stl: {args.output}")
    print(f"Intermediate OpenSCAD: {intermediate_openscad}")
    print(f"OpenSCAD modules: {args.modules}")
    print(f"Path to OpenSCAD: {args.openscad}")
    print("Parsing .yaml ...")
    scad_variables = translate_yaml_to_scad_variables(args.input)
    create_scad_file(scad_variables=scad_variables, output_filename="intermediate.scad", renderer_filename=args.modules)
    print("Creating .stl ...")
    subprocess.run([args.openscad, "-o", args.output, intermediate_openscad])
    print("Done!")
