# FLAGS
BUILD_DIR ?= build
EXECNAME = fortepiano
# define TESTSPEED=1 to compute the time required for the first 1000 derivatives
TESTSPEED ?=
# define FULL_F_AB=1 to compute full matrix product in F_AB functions (by default assumes diagonal G matrices)
FULL_F_AB ?=
# define NOINTERPOLATION=1 to use the full expressions instead of interpolations for the energy densities and electromagnetic corrections
NOINTERPOLATION ?=
F90 ?= ifort

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

ifeq ($(TESTSPEED), 1)
	F90FLAGS += -DTESTSPEED=1
endif
ifeq ($(FULL_F_AB), 1)
	F90FLAGS += -DFULLFAB=1
endif
ifeq ($(NOINTERPOLATION), 1)
	F90FLAGS += -DNOINTERPOLATION=1
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

ifeq ($(BUILD_DIR),build_debug)
FFLAGS=$(DEBUGFLAGS)
else
FFLAGS=$(F90FLAGS)
endif

default: all

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
	$(BUILD_DIR)/interactions.o $(BUILD_DIR)/utilities.o $(BUILD_DIR)/bspline_module.o \
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
$(BUILD_DIR)/test_utils.o: $(BUILD_DIR)/const.o
$(BUILD_DIR)/tests.o: $(OBJ_TESTS)

all: fortepiano

directories:
	mkdir -p bin/ log/ $(BUILD_DIR)

fortepiano: directories Makefile $(OBJ_FILES) $(BUILD_DIR)/fortepiano.o
	$(F90) -o bin/$(EXECNAME) $(OBJ_FILES) $(BUILD_DIR)/fortepiano.o $(FFLAGS)

tests: directories Makefile $(OBJ_TESTS) $(BUILD_DIR)/tests.o
	$(F90) -o bin/tests $(OBJ_TESTS) $(BUILD_DIR)/tests.o $(FFLAGS)

clean: 
	rm -rf bin/* $(BUILD_DIR)*/ build*/

$(BUILD_DIR)/%.o: sources/%.f90 Makefile
	$(F90) $(FFLAGS) -c sources/$*.f90 -o $(BUILD_DIR)/$*.o

$(BUILD_DIR)/%.o: sources/%.f Makefile
	$(F90) $(FFLAGS) $(stdFlag) -c sources/$*.f -o $(BUILD_DIR)/$*.o
