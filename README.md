# windio2cad

## What is `windio2cad`?

The `windio2cad` module is a command line tool to convert `.yaml` files that follow the [windio ontology](https://windio.readthedocs.io/en/latest/) into `.stl` CAD files. These CAD files can then be used as input to other codes or for visualization.

Two example files for floating offshore wind turbine (FOWT) platforms are included in this source repository:

- `nrel5mw-spar_oc4.yaml`: A semisubmersible floating platform.

- `nrel5mw-spar_oc3.yaml`: A spar-type floating platform.

Read below to learn more about rendering these `.yaml` files into `.stl` files.

## Installation

This instllation guide assumes you are using `conda` for management of virtual environments. Clone the repository and change into the repository's directory:

```
conda create -n windio2cad python=3.8
conda activate windio2cad
pip install -e .
```

This module depends on the OpenSCAD package to be installed to render the final `.stl` file. See [http://www.openscad.org/](http://www.openscad.org/)

## Usage

Use the following command to create an `.stl` file for a turbine atop a semisubmersible foundation:

```
python -m windio2cad --input nrel5mw-semi_oc4.yaml --output turbine.stl --openscad [PATH TO OPENSCAD]
```

Replace `[PATH TO OPENSCAD]` with the path to OpenSCAD on your operating system.

### Rendering just the blade

In some applications, only a rendering of the blade is needed. In those cases, simply add the `--blade` option to the command line to render just the blade.
