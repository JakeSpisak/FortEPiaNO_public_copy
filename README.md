# [FortEPiaNO](http://www.bitbucket.org?/...)
# FORTran-Evolved PrImordiAl Neutrino Oscillations
# by S. Gariazzo (stefano.gariazzo@gmail.com) and P.F. de Salas ()

FortEPiaNO is a Fortran code to compute the evolution of neutrino oscillations in the early universe.
The code is written to flexible, it can work with two to six neutrinos (active or sterile).

If you use this code for scientific publications, please cite the paper:
S. Gariazzo, P.F. de Salas, S. Pastor
...
arxiv:...
submitted to...
https://inspire....

## 1.Install and compile
Clone the repository ... or download its content ..., go to folder, compile with `make`, run with `bin/fortepiano.exe ini/explanatory.ini`

The `ini/explanatory.ini` contains a detailled description of the various arguments that the fortran code accepts.

### 1.1.More on compiling
Compiled by default with `ifort`, if present.
To force compilation with a `gfortran`, use `make F90=gfortran`
Compiled modules are stored in the `build/` folder by default. To change this, compile with `BUILD_DIR`,
for example `make BUILD_DIR=newbuildfolder`.
`make clean` to remove all compiled files.
`make tests && bin/tests` will run a set of numerical tests. If you modify the code, they will tell you if everything is still working as expected.

## 2.Python scripts
Some useful python scripts are delivered together with the main fortran code.

### 2.1.`prepareIni.py`
Tool that will help to generate an ini file, using default values for the non-specified parameters.
Use `python python/prepareIni.py -h` to get a list of accepted arguments and their description.
Currently accepts only up to 4 neutrinos (3 active + 1 sterile)

### 2.2.`readOutput.py`
Functions to read the output files generated by the fortran code and generate plots.
It can be used importing the main class `...`, for example:
```
from nuDensOutput import NuDensRun
run = NuDensRun("output/folder/", nnu=4, plots=True)
```
For possible plotting functions include:
* ...

### 2.3.`dogrid3p1.py`
Functions to generate the ini files, submit runs, produce plots for a specific grid.
See `python python/dogrid3p1.py -h` for a detailled description of the accepted arguments.

## 3.Source code
Approximate content of the files:
* `config.f90`: routines for the initialization
* `const.f90`: constants, variables, interfaces
* `cosmology.f90`: energy density of the various species
* `equations.f90`: evolution equations of the various quantities
* `interactions`: functions that define the neutrino interactions and the collision terms
* ``: main program

Auxiliary files:
* `bspline*`: spline interpolation utilities
* `errors.f90`: error and log management
* `HEigensystem.f90`: complex matrix diagonalization
* `IniFile.f90`: functions to read the input file
* `linear_interpolation_module.f90`: linear interpolation utilities
* `matrix_utils.f90`: utilities for creating and manipulating matrices
* `odepack*`: DLSODA and related utilities
* `stuff.f90`: old functions that were used in previous versions of the code and now enter only the tests
* `tests.f90`: numerical tests for the software
* `test_utils.f90`: assertion functions and test counts
* `utilities.f90`: utilities for integration, interpolation, checkpointing, log file manipulation, time measurements

## 4.Acknowledgments
This software is part of a project that has received funding from the European Union's Horizon 2020 research and innovation programme, under the Marie Skłodowska-Curie grant agreement No 796941.
