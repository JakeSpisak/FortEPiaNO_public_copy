# [FortEPiaNO](https://bitbucket.org/ahep_cosmo/fortepiano_public/)
FORTran-Evolved PrImordiAl Neutrino Oscillations  
by S. Gariazzo (gariazzo@ific.uv.es) and P.F. de Salas (pablo.fernandez@fysik.su.se)

FortEPiaNO is a Fortran code to compute the evolution of neutrino oscillations in the early universe.  
The code is written to flexible, it can work with two to six neutrinos (active or sterile).  
At the moment, no lepton asymmetries nor non-standard interactions are implemented.

If you use this code for scientific publications, please cite the paper:  
**Thermalisation of sterile neutrinos in the early Universe in the 3+1 scheme with full mixing matrix**  
_S. Gariazzo, P.F. de Salas, S. Pastor_  
JCAP 07 (2019) 014.  
[arxiv:1905.11290](https://arxiv.org/abs/1905.11290),  
see also on [INSPIRE](https://inspirehep.net/record/1736955).  
```
@article{Gariazzo:2019gyi,
      author         = "Gariazzo, S. and de Salas, P. F. and Pastor, S.",
      title          = "{Thermalisation of sterile neutrinos in the early Universe in the 3+1 scheme with full mixing matrix}",
      journal        = "JCAP",
      volume         = "07",
      year           = "2019",
      number         = "07",
      pages          = "014",
      doi            = "10.1088/1475-7516/2019/07/014",
      eprint         = "1905.11290",
      archivePrefix  = "arXiv",
      primaryClass   = "astro-ph.CO",
      SLACcitation   = "%%CITATION = ARXIV:1905.11290;%%"
}
```


## 1.Install and compile
In order to install the code, the suggested option is to clone the git repository.  
To do so, go to the webpage [https://bitbucket.org/ahep_cosmo/fortepiano_public/](https://bitbucket.org/ahep_cosmo/fortepiano_public/), and use the button on the top right to obtain the clone link (`git clone https://bitbucket.org/ahep_cosmo/fortepiano_public.git`).  
You can also download a compressed folder with all the content: in such case, uncompress it in a new folder.

In either cases, when you have the source code, open the folder which contains it in a terminal.  
The main code can be compiled simply using `make` (see also the section below on more options for compiling), which will generate an executable `bin/fortepiano`.  
Your first test can be to use the explanatory configuration file, as: `bin/fortepiano ini/explanatory.ini`.

The `ini/explanatory.ini` contains a detailled description of the various arguments that the fortran code accepts, and shows how to pass the neutrino mixing parameters.
You are encouraged to look inside this file to understand what parameters you can vary.

### 1.1.More on compiling
The Makefile should automatically recognize which compiler to use.
If it is present, `ifort` will be used by default, otherwise the GNU fortran compiler `gfortran` is used.
To force compilation with a `gfortran` when `ifort` exists, use `make F90=gfortran`.

Compiled modules are stored in the `build/` folder by default. You can however specify the folder where to save the builds using `BUILD_DIR`, for example `make BUILD_DIR=newbuildfolder`.
If you set the build folder to `build_debug`, a set of options useful for debug will be used for compiling.

Additional commands for the makefile include:

* `make clean` to remove all compiled files.
* `make tests` to compile a set of numerical tests (run them with `bin/tests`). If you modify the code, they will tell you if everything is still working as expected.

Moreover, some parts of the code can be enabled or disabled using precompilation flags.
For example (add the option to the `make` command):

* `NO_MUONS=1` disables the contribution of muons to the energy budget of the universe;
* `INTERP_DIV_FD=1` changes the interpolation of the neutrino density matrix to divide by the Fermi Dirac distribution before the interpolation (and multiply back after);
* `FULL_F_NU=1` allows to use the full matrix product in the F_nu phase space functions that compute nunu scattering and pair annihilation. If not used, only diagonal elements of the neutrino density matrix will be used;
* `FULL_F_AB=1` allows to use the full matrix product in the F_ab functions that appear in the collision integrals.

**WARNING**: the test suite will not work if the flag `NO_MUONS=1` is activated, or even if some modules have been compiled with that option. You will need to use `make clean` before `make tests` in order to be sure that everything works.

### 1.2.Interpolations
The code, in the default compilation setup, is designed to avoid computing several integrals at each step and to use an interpolation instead.
This mostly concerns integrals of energy densities (photons, charged leptons) and of the functions that describe electromagnetic corrections, including the coefficients in the dz/dx and dw/dx expressions.

In order to change this default behaviour and use the full expressions at each step, at the expense of slightly slowing down the calculation, one can compile the code using `make NOINTERPOLATION=1`, so that the interpolations will be disabled.
For the default configuration in the `ini/explanatory.ini` file, the actual change in the final effective number of neutrinos is smaller than `1e-4` and the execution takes approximately 15% longer.

**WARNING**: the test suite will not compile if the flag `NOINTERPOLATION=1` is activated, or if some modules have been compiled with that option. You will need to use `make clean` before `make tests` in order to be sure that everything works.


## 2.Python scripts
Some useful python scripts are delivered together with the main fortran code.

In order to be sure that all the python scripts work properly, you should install the required dependencies using `python setup.py install --user`.

### 2.1.`prepareIni.py`
Tool that will help to generate an ini file, using default values for the non-specified parameters.  
Use `python python/prepareIni.py -h` to get a list of accepted arguments and their description.  
Currently accepts only up to 4 neutrinos (3 active + 1 sterile).

### 2.2.`fortepianoOutput.py`
Functions to read the output files generated by the fortran code and generate plots.
It can be used importing the main class `FortEPiaNORun`, for example:
```
from fortepianoOutput import FortEPiaNORun
run = FortEPiaNORun("output/folder/", nnu=4, plots=True)
```

For possible plotting functions include:

* `plotFD`: plot the Fermi-Dirac distribution computed in the current momentum grid;
* `plotZ`: plot the evolution of the photon temperature `z` as a function of `x`;
* `plotW`: plot the evolution of the effective neutrino temperature `w` as a function of `x`;
* `plotZoverW`: plot the evolution of the ratio `z/w` as a function of `x`;
* `plotDeltaZ`: plot the relative difference of the evolution of the photon temperature `z` as a function of `x`, with respect to some other given run;
* `plotRhoDiagY`, `plotRhoOffDiagY`: plot the evolution of an element (diagonal or off-diagonal) of the neutrino density matrix at a fixed momentum `y`;
* `plotdRhoOffDiagY`: plot a numerical estimate of the `x` derivative of an (off-diagonal) element of the neutrino density matrix at a fixed momentum `y`;
* `plotRhoFin`: plot the momentum-dependence of the requested element of the neutrino density matrix at the final `x`;
* `plotRhoX`: plot the momentum-dependence of the requested element of the neutrino density matrix, at a given `x`;
* `plotNeff`: plot the evolution of the effective number of neutrinos, which is correctly normalized only at very early times or today;
* `plotEnergyDensity`: plot the evolution of the energy density of the different components and their sum;
* `plotEntropy`: plot the evolution of the entropy density of the different components and their sum;
* `doAllPlots`: it will create a series of pdf plots in the output folder, calling
`plotZ`, `plotW`, `plotRhoFin`, `plotRhoDiagY`, `plotRhoOffDiagY`, `plotdRhoOffDiagY`,
for flavor and mass eigenstates if available.

Additionally, the function `integrateRho_yn` provides a fast way to compute the integral of `y^n f(y)`, interpolating over the final energy density that was computed by `FortEPiaNO`.

### 2.3.`tests.py`
Testing tools to verify that everything works in the python part, from reading the output folder to plotting.
Run it if you want to test that you have all the required packages and everything is working in your current setup.


## 3.Source code
Should you need to edit the source codes, this is more or less the content of each file:

* `config.f90`: routines for the initialization;
* `const.f90`: constants, variables, interfaces;
* `cosmology.f90`: energy density of the various species;
* `equations.f90`: evolution equations of the various quantities;
* `fortepiano.f90`: main program;
* `ftqed.f90`: functions and utilities to compute finite-temperature QED corrections;
* `interactions.f90`: functions that define the neutrino interactions and the collision terms;
* `tests.f90`: numerical tests for the software.

Auxiliary files, which in principle you should not need to edit:

* `bspline*`: spline interpolation utilities;
* `errors.f90`: error and log management;
* `heigensystem.f90`: complex matrix diagonalization;
* `iniFile.f90`: functions to read the input file;
* `linear_interpolation_module.f90`: linear interpolation utilities;
* `matrix_utils.f90`: utilities for creating and manipulating matrices;
* `odepack*`: DLSODA and related utilities;
* `stuff.f90`: old functions that were used in previous versions of the code and now enter only the tests;
* `test_utils.f90`: assertion functions and test counts;
* `utilities.f90`: utilities for integration, interpolation, checkpointing, log file manipulation, time measurements.


## 4.Acknowledgments
This software is part of a project that has received funding from the European Union's Horizon 2020 research and innovation programme, under the Marie Skłodowska-Curie grant agreement No 796941.
