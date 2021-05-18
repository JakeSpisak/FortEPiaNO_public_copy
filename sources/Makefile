# FLAGS
# use G_L and G_R values at zero-momentum transfer instead of the default ones
GLR_ZERO_MOMENTUM ?=
# define FULL_F_AB=1 to compute full matrix product in F_AB functions (by default assumes diagonal G matrices)
FULL_F_AB ?=
# define FULL_F_NU=1 to compute full matrix product in F_nu functions. If not set, use only diagonal elements
FULL_F_NU ?=
# enable low reheating
LOW_REHEATING ?=
# define NO_INTERPOLATION=1 to use the full expressions instead of interpolations for the energy densities and electromagnetic corrections
NO_INTERPOLATION ?=
# set DO_MUONS=1 to enable muon contributions
DO_MUONS ?=
# set NO_NUE_ANNIHILATION=1 to disable contribution to collision integrals from nunu<->e+e- annihilation processes
NO_NUE_ANNIHILATION ?=
# set RHO_OFFDIAG_INTERP_DIV_FD=1 to interpolate all the entries of the neutrino density matrix after dividing by a Fermi-Dirac distribution (by default, this is done only for diagonal entries)
RHO_OFFDIAG_INTERP_DIV_FD ?=
#use this flag to set a custom value for the weak mixing angle (for example SINSQTHW=0.23)
SINSQTHW ?= 
# define TESTSPEED=1 to compute the time required for the first 1000 derivatives
TESTSPEED ?=

F90 ?= gfortran

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
  F90FLAGS=-O3 -fPIC -fpp -L/usr/lib -I$(BUILD_DIR)/ -openmp -parallel -par-report1 -no-prec-div
DEBUGFLAGS=-O0 -fpp -L/usr/lib -I$(BUILD_DIR)/ -p -g -traceback -openmp -fpe0 -check all
MOD_DIR=-module $(BUILD_DIR)/
else
#gfortran
F90=gfortran
  F90FLAGS=-O3 -fPIC -cpp -L/usr/lib -I$(BUILD_DIR)/ -ffast-math -ffree-line-length-none -fopenmp
DEBUGFLAGS=-O0 -cpp -L/usr/lib- -I$(BUILD_DIR)/ -p -g -fbacktrace -ffast-math -ffree-line-length-none -fopenmp -fcheck=all
MOD_DIR=-J$(BUILD_DIR)/
endif

ifeq ($(BUILD_DIR),../build_debug)
FFLAGS=$(DEBUGFLAGS)
else
FFLAGS=$(F90FLAGS)
endif

ifeq ($(GLR_ZERO_MOMENTUM), 1)
	FFLAGS += -DGLR_ZERO_MOMENTUM=1
endif
ifeq ($(FULL_F_AB), 1)
	FFLAGS += -DFULL_F_AB=1
endif
ifeq ($(FULL_F_NU), 1)
	FFLAGS += -DFULL_F_NU=1
endif
ifeq ($(LOW_REHEATING), 1)
	FFLAGS += -DLOW_REHEATING=1 -DNO_INTERPOLATION=1
endif
ifeq ($(NO_INTERPOLATION), 1)
	FFLAGS += -DNO_INTERPOLATION=1
endif
ifeq ($(DO_MUONS), 1)
	FFLAGS += -DDO_MUONS=1
endif
ifeq ($(NO_NUE_ANNIHILATION), 1)
	FFLAGS += -DNO_NUE_ANNIHILATION=1
endif
ifeq ($(RHO_OFFDIAG_INTERP_DIV_FD), 1)
	F90FLAGS += -DRHO_OFFDIAG_INTERP_DIV_FD=1
endif
ifneq ($(SINSQTHW),)
	F90FLAGS += -DSINSQTHW=$(SINSQTHW)
endif
ifeq ($(TESTSPEED), 1)
	FFLAGS += -DTESTSPEED=1
endif

# used in prepare_interpolations:
INTERPFLAGS=
ifneq ($(XIN),)
	INTERPFLAGS += -DXIN=$(XIN)
endif
ifneq ($(XFIN),)
	INTERPFLAGS += -DXFIN=$(XFIN)
endif
ifneq ($(YMIN),)
	INTERPFLAGS += -DYMIN=$(YMIN)
endif
ifneq ($(YMAX),)
	INTERPFLAGS += -DYMAX=$(YMAX)
endif
ifneq ($(STARTX),)
	INTERPFLAGS += -DSTARTX=$(STARTX)
endif

OBJ_FILES=$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/config.o \
	$(BUILD_DIR)/iniFile.o \
	$(BUILD_DIR)/input.o \
	$(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/matrix_utils.o \
	$(BUILD_DIR)/test_utils.o \
	$(BUILD_DIR)/ftqed.o \
	$(BUILD_DIR)/interactions.o \
	$(BUILD_DIR)/cosmology.o \
	$(BUILD_DIR)/matter.o \
	$(BUILD_DIR)/output.o \
	$(BUILD_DIR)/equations.o \
	$(BUILD_DIR)/heigensystem.o \
	$(BUILD_DIR)/odepack.o \
	$(BUILD_DIR)/odepack-sub1.o \
	$(BUILD_DIR)/odepack-sub2.o \
	$(BUILD_DIR)/bspline_module.o \
	$(BUILD_DIR)/bspline_oo_module.o \
	$(BUILD_DIR)/bspline_sub_module.o \
	$(BUILD_DIR)/linear_interpolation_module.o

OBJ_TESTS=$(OBJ_FILES) $(BUILD_DIR)/stuff.o

default: fortepiano

$(BUILD_DIR)/bspline_module.o: \
	$(BUILD_DIR)/bspline_oo_module.o \
	$(BUILD_DIR)/bspline_sub_module.o
$(BUILD_DIR)/bspline_oo_module.o: \
	$(BUILD_DIR)/bspline_sub_module.o
$(BUILD_DIR)/odepack.o: \
	$(BUILD_DIR)/odepack-sub1.o \
	$(BUILD_DIR)/odepack-sub2.o
$(BUILD_DIR)/iniFile.o: \
	$(BUILD_DIR)/const.o
$(BUILD_DIR)/matrix_utils.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/errors.o
$(BUILD_DIR)/heigensystem.o: \
	$(BUILD_DIR)/const.o
$(BUILD_DIR)/ftqed.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/test_utils.o \
	$(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/config.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/interactions.o \
	$(BUILD_DIR)/iniFile.o \
	$(BUILD_DIR)/matrix_utils.o \
	$(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/input.o \
	$(BUILD_DIR)/equations.o
$(BUILD_DIR)/input.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/iniFile.o \
	$(BUILD_DIR)/errors.o
$(BUILD_DIR)/interactions.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/matrix_utils.o \
	$(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/bspline_module.o \
	$(BUILD_DIR)/ftqed.o \
	$(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/cosmology.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/bspline_module.o \
	$(BUILD_DIR)/test_utils.o \
	$(BUILD_DIR)/ftqed.o \
	$(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/equations.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/cosmology.o \
	$(BUILD_DIR)/interactions.o \
	$(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/test_utils.o \
	$(BUILD_DIR)/heigensystem.o \
	$(BUILD_DIR)/ftqed.o \
	$(BUILD_DIR)/matter.o \
	$(BUILD_DIR)/output.o \
	$(BUILD_DIR)/bspline_module.o \
	$(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/matter.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/cosmology.o \
	$(BUILD_DIR)/heigensystem.o
$(BUILD_DIR)/output.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/cosmology.o \
	$(BUILD_DIR)/matter.o \
	$(BUILD_DIR)/utilities.o
$(BUILD_DIR)/stuff.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/cosmology.o \
	$(BUILD_DIR)/interactions.o \
	$(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/equations.o \
	$(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/utilities.o: \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/errors.o
$(BUILD_DIR)/test_utils.o: \
	$(BUILD_DIR)/const.o
$(BUILD_DIR)/fortepiano.o: $(OBJ_FILES)
$(BUILD_DIR)/tests.o: $(OBJ_TESTS)

PREP_NODES_OBJ = \
	$(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/const.o \
	$(BUILD_DIR)/errors.o
$(BUILD_DIR)/prepare_gl_nodes.o: $(PREP_NODES_OBJ)

$(BUILD_DIR)/prepare_interpolations.o: $(OBJ_FILES)

all: fortepiano tests

readnodes:
	$(F90) -o ../bin/read_gl_nodes read_gl_nodes.f90

preparenodes: $(BUILD_DIR)/prepare_gl_nodes.o
	$(F90) -o ../bin/prepare_gl_nodes $(PREP_NODES_OBJ) $(BUILD_DIR)/prepare_gl_nodes.o $(FFLAGS) $(MOD_DIR) $(TESTFLAGS)

prepareinterp: $(PREP_NODES_OBJ) $(BUILD_DIR)/prepare_interpolations.o
	$(F90) -o ../bin/prepare_interpolations $(OBJ_FILES) $(BUILD_DIR)/prepare_interpolations.o $(FFLAGS) $(MOD_DIR) $(TESTFLAGS) $(INTERPFLAGS)

tests: TESTFLAGS = -DDO_TESTS=1

fortepiano: Makefile $(OBJ_FILES) $(BUILD_DIR)/fortepiano.o
	$(F90) -o ../bin/$(EXECNAME) $(OBJ_FILES) $(BUILD_DIR)/fortepiano.o $(FFLAGS) $(MOD_DIR) $(TESTFLAGS)

pythonwrapper: Makefile $(OBJ_FILES)
	python -m numpy.f2py -c --fcompiler=$(F90) --opt='$(FFLAGS)' $(OBJ_FILES) -lgomp -m fortepianoWrapper pythonWrapper.f90
	python testwrapper.py
	cp fortepianoWrapper*.so ../python

tests: Makefile $(OBJ_TESTS) $(BUILD_DIR)/tests.o
	$(F90) -o ../bin/tests $(OBJ_TESTS) $(BUILD_DIR)/tests.o $(FFLAGS) $(MOD_DIR) $(TESTFLAGS)

$(BUILD_DIR)/%.o: %.f90 Makefile
	$(F90) $(FFLAGS) $(MOD_DIR) $(TESTFLAGS) $(INTERPFLAGS) -c $*.f90 -o $(BUILD_DIR)/$*.o

$(BUILD_DIR)/%.o: %.f Makefile
	$(F90) $(FFLAGS) $(MOD_DIR) $(TESTFLAGS) $(INTERPFLAGS) $(stdFlag) -c $*.f -o $(BUILD_DIR)/$*.o
