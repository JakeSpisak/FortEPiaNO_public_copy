#!/usr/bin/env python
"""Script for installing dependencies for the python scripts of FortEPiaNO"""

from setuptools import setup


def readme():
    with open("README.md") as f:
        return f.read()


def version():
    with open("VERSION") as f:
        return f.read().split(" ")[0]


setup(
    name="fortepiano",
    version=version(),
    description="support scripts for FortEPiaNO",
    long_description_content_type="text/markdown",
    long_description=readme(),
    author="S. Gariazzo",
    author_email="stefano.gariazzo@gmail.com",
    url="https://bitbucket.org/ahep_cosmo/fortepiano_public/",
    license="GPL-3.0",
    keywords="neutrino decoupling, neutrino oscillations, early universe, Neff, effective number of relativistic species",
    install_requires=[
        "argparse",
        'matplotlib(>2,<3);python_version<"3"',
        'matplotlib(>3);python_version>"3"',
        "numpy",
        "scipy",
        "python-ternary",
        'mock;python_version<"3"',
        'unittest2;python_version<"3"',
    ],
    data_files=[("python", ["LICENSE.txt", "VERSION", "icon.png"])],
)
