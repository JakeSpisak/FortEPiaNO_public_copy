#!/usr/bin/env python
"""Script for installing dependencies for the python scripts of FortEPiaNO"""

from setuptools import setup


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="fortepiano",
    description="support scripts for FortEPiaNO",
    long_description_content_type="text/markdown",
    long_description=readme(),
    author="S. Gariazzo",
    author_email="stefano.gariazzo@gmail.com",
    url="https://bitbucket.org/ahep_cosmo/piano/",
    # license="GPL-3.0",
    # keywords="",
    install_requires=[
        "argparse",
        'matplotlib(>2,<3);python_version<"3"',
        'matplotlib(>3);python_version>"3"',
        "numpy",
        "scipy",
        "python-ternary",
    ],
)
