#!/usr/bin/env python

from setuptools import setup, find_packages

# Top-level setup
setup(
    name             = 'windio2cad',
    version          = '0.0.1',
    description      = 'windio2cad',
    long_description =  '''windio2cad takes windio ontology YAML files and translates them into .stl files.''',
    url              = 'https://github.com/WISDEM/WISDEM',
    author           = 'NREL WISDEM Team',
    author_email     = 'systems.engineering@nrel.gov',
    install_requires = [
        'pyyaml',
        'scipy',
        'matplotlib',
        'numpy',
        'ruamel_yaml',
    ],
    python_requires  = '>=3.8',
    packages         = find_packages(exclude=['docs', 'tests', 'ext']),
    license          = 'Apache License, Version 2.0',
    zip_safe         = False
)
