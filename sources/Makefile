# FLAGS
# use G_L and G_R values at zero-momentum transfer instead of the default ones
GLR_ZERO_MOMENTUM ?=
# define FULL_F_AB=1 to compute full matrix product in F_AB functions (by default assumes diagonal G matrices)
FULL_F_AB ?=
# define FULL_F_NU=1 to compute full matrix product in F_nu functions. If not set, use only diagonal elements
FULL_F_NU ?=
# define NO_INTERPOLATION=1 to use the full expressions instead of interpolations for the energy densities and electromagnetic corrections
NO_INTERPOLATION ?=
# set NO_MUONS=1 to disable muon contributions
NO_MUONS ?=
# set NO_NUE_ANNIHILATION=1 to disable contribution to collision integrals from nunu<->e+e- annihilation processes
NO_NUE_ANNIHILATION ?=
# set RHO_OFFDIAG_INTERP_DIV_FD=1 to interpolate all the entries of the neutrino density matrix after dividing by a Fermi-Dirac distribution (by default, this is done only for diagonal entries)
RHO_OFFDIAG_INTERP_DIV_FD ?=
#use this flag to set a custom value for the weak mixing angle (for example SINSQTHW=0.23)
SINSQTHW ?= 
# define TESTSPEED=1 to compute the time required for the first 1000 derivatives
TESTSPEED ?=

F90 ?= ifort

BUILD_DIR ?= ../build

EXECNAME = fortepiano

ifeq ("$(F90)","gfortran")
ifortErr=1
else
ifortErr=$(shell which ifort >/dev/null; echo $$?)
endif

ifeq "$(ifortErr)" "0"
#ifort
F90=ifort
  F90FLAGS=-O3 -fpp -L/usr/lib -I$(BUILD_DIR)/ -module $(BUILD_DIR)/ -openmp -parallel -par-report1 -no-prec-div
DEBUGFLAGS=-O0 -fpp -L/usr/lib -I$(BUILD_DIR)/ -module $(BUILD_DIR)/ -p -g -traceback -openmp -fpe0 -check all

else
#gfortran
F90=gfortran
  F90FLAGS=-O3 -cpp -L/usr/lib -J$(BUILD_DIR)/ -I$(BUILD_DIR)/ -ffast-math -ffree-line-length-none -fopenmp
DEBUGFLAGS=-O0 -cpp -L/usr/lib -J$(BUILD_DIR)/ -I$(BUILD_DIR)/ -p -g -fbacktrace -ffast-math -ffree-line-length-none -fopenmp -fcheck=all
endif

ifeq ($(GLR_ZERO_MOMENTUM), 1)
	F90FLAGS += -DGLR_ZERO_MOMENTUM=1
endif
ifeq ($(FULL_F_AB), 1)
	F90FLAGS += -DFULL_F_AB=1
endif
ifeq ($(FULL_F_NU), 1)
	F90FLAGS += -DFULL_F_NU=1
endif
ifeq ($(NO_INTERPOLATION), 1)
	F90FLAGS += -DNO_INTERPOLATION=1
endif
ifeq ($(NO_MUONS), 1)
	F90FLAGS += -DNO_MUONS=1
endif
ifeq ($(NO_NUE_ANNIHILATION), 1)
	F90FLAGS += -DNO_NUE_ANNIHILATION=1
endif
ifeq ($(RHO_OFFDIAG_INTERP_DIV_FD), 1)
	F90FLAGS += -DRHO_OFFDIAG_INTERP_DIV_FD=1
endif
ifneq ($(SINSQTHW),)
	F90FLAGS += -DSINSQTHW=$(SINSQTHW)
endif
ifeq ($(TESTSPEED), 1)
	F90FLAGS += -DTESTSPEED=1
endif

OBJ_FILES=$(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o $(BUILD_DIR)/config.o \
	$(BUILD_DIR)/iniFile.o $(BUILD_DIR)/utilities.o $(BUILD_DIR)/matrix_utils.o \
	$(BUILD_DIR)/ftqed.o \
	$(BUILD_DIR)/interactions.o $(BUILD_DIR)/cosmology.o $(BUILD_DIR)/equations.o \
	$(BUILD_DIR)/heigensystem.o \
	$(BUILD_DIR)/odepack.o $(BUILD_DIR)/odepack-sub1.o $(BUILD_DIR)/odepack-sub2.o \
	$(BUILD_DIR)/bspline_module.o $(BUILD_DIR)/bspline_oo_module.o $(BUILD_DIR)/bspline_sub_module.o \
	$(BUILD_DIR)/linear_interpolation_module.o

OBJ_TESTS=$(OBJ_FILES) $(BUILD_DIR)/stuff.o $(BUILD_DIR)/test_utils.o

ifeq ($(BUILD_DIR),../build_debug)
FFLAGS=$(DEBUGFLAGS)
else
FFLAGS=$(F90FLAGS)
endif

default: fortepiano

$(BUILD_DIR)/bspline_module.o: $(BUILD_DIR)/bspline_oo_module.o \
	$(BUILD_DIR)/bspline_sub_module.o
$(BUILD_DIR)/bspline_oo_module.o: $(BUILD_DIR)/bspline_sub_module.o
$(BUILD_DIR)/odepack.o: $(BUILD_DIR)/odepack-sub1.o $(BUILD_DIR)/odepack-sub2.o
$(BUILD_DIR)/iniFile.o: $(BUILD_DIR)/const.o
$(BUILD_DIR)/matrix_utils.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o
$(BUILD_DIR)/heigensystem.o: $(BUILD_DIR)/const.o
$(BUILD_DIR)/ftqed.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o $(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/config.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/interactions.o \
	$(BUILD_DIR)/iniFile.o $(BUILD_DIR)/matrix_utils.o $(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/equations.o
$(BUILD_DIR)/interactions.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/matrix_utils.o $(BUILD_DIR)/utilities.o $(BUILD_DIR)/bspline_module.o \
	$(BUILD_DIR)/ftqed.o \
	$(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/cosmology.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/utilities.o $(BUILD_DIR)/bspline_module.o \
	$(BUILD_DIR)/ftqed.o \
	$(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/equations.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/cosmology.o $(BUILD_DIR)/interactions.o $(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/heigensystem.o \
	$(BUILD_DIR)/ftqed.o \
	$(BUILD_DIR)/bspline_module.o $(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/fortepiano.o: $(OBJ_FILES)
$(BUILD_DIR)/stuff.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/cosmology.o $(BUILD_DIR)/interactions.o $(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/equations.o $(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/utilities.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o
$(BUILD_DIR)/test_utils.o: $(BUILD_DIR)/const.o
$(BUILD_DIR)/tests.o: $(OBJ_TESTS)

PREP_NODES_OBJ = $(BUILD_DIR)/utilities.o $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o
$(BUILD_DIR)/prepare_gl_nodes.o: $(PREP_NODES_OBJ)

all: fortepiano tests

readnodes:
	$(F90) -o ../bin/read_gl_nodes read_gl_nodes.f90

preparenodes: $(BUILD_DIR)/prepare_gl_nodes.o
	$(F90) -o ../bin/prepare_gl_nodes $(PREP_NODES_OBJ) $(BUILD_DIR)/prepare_gl_nodes.o $(FFLAGS) $(TESTFLAGS)

tests: TESTFLAGS = -DDO_TESTS=1

fortepiano: Makefile $(OBJ_FILES) $(BUILD_DIR)/fortepiano.o
	$(F90) -o ../bin/$(EXECNAME) $(OBJ_FILES) $(BUILD_DIR)/fortepiano.o $(FFLAGS) $(TESTFLAGS)

tests: Makefile $(OBJ_TESTS) $(BUILD_DIR)/tests.o
	$(F90) -o ../bin/tests $(OBJ_TESTS) $(BUILD_DIR)/tests.o $(FFLAGS) $(TESTFLAGS)

$(BUILD_DIR)/%.o: %.f90 Makefile
	$(F90) $(FFLAGS) $(TESTFLAGS) -c $*.f90 -o $(BUILD_DIR)/$*.o

$(BUILD_DIR)/%.o: %.f Makefile
	$(F90) $(FFLAGS) $(TESTFLAGS) $(stdFlag) -c $*.f -o $(BUILD_DIR)/$*.o
