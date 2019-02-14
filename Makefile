# FLAGS
BUILD_DIR ?= build
EXECNAME ?= nuDens.exe
# use e.g. with USER_DEFINED=-DTESTSPEED=1 and EXECNAME=nuDens_speed.exe
USER_DEFINED ?= 
LOGY ?= 
#LOGY ?= -DLOGY=1

#ifort
F90=ifort
  F90FLAGS=-O3 -fpp -L/usr/lib -I$(BUILD_DIR)/ -module $(BUILD_DIR)/ -p -g -traceback -openmp -parallel -par-report1 -no-prec-div $(USER_DEFINED) $(LOGY)
DEBUGFLAGS=-O0 -fpp -L/usr/lib -I$(BUILD_DIR)/ -module $(BUILD_DIR)/ -p -g -traceback -openmp -fpe0 -check all $(USER_DEFINED) $(LOGY)
# -check all -check noarg_temp_created
# -stand f03  -check all -warn all -fstack-protector -assume protect_parens -implicitnone
# -openmp 
#-openmp-stubs 
#-openmp       
stdFlag=

#gfortran
#F90=gfortran
#F90FLAGS=-O3 -L/usr/lib -J$(BUILD_DIR)/ -I$(BUILD_DIR)/ -g -traceback -fopenmp
#stdFlat= --std=legacy

OBJ_FILES=$(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o $(BUILD_DIR)/config.o \
	$(BUILD_DIR)/IniFile.o $(BUILD_DIR)/utilities.o $(BUILD_DIR)/matrix_utils.o \
	$(BUILD_DIR)/interactions.o $(BUILD_DIR)/cosmology.o $(BUILD_DIR)/equations.o \
	$(BUILD_DIR)/odepack.o $(BUILD_DIR)/odepack-sub1.o $(BUILD_DIR)/odepack-sub2.o \
	$(BUILD_DIR)/bspline_module.o $(BUILD_DIR)/bspline_oo_module.o $(BUILD_DIR)/bspline_sub_module.o \
	$(BUILD_DIR)/linear_interpolation_module.o
#$(BUILD_DIR)/opkdmain.o $(BUILD_DIR)/opkda1.o $(BUILD_DIR)/opkda2.o

ifeq ($(BUILD_DIR),builddeb)
FFLAGS=$(DEBUGFLAGS)
else
FFLAGS=$(F90FLAGS)
endif

default: all

$(BUILD_DIR)/bspline_module.o: $(BUILD_DIR)/bspline_oo_module.o \
	$(BUILD_DIR)/bspline_sub_module.o
$(BUILD_DIR)/bspline_oo_module.o: $(BUILD_DIR)/bspline_sub_module.o
#$(BUILD_DIR)/opkdmain.o: $(BUILD_DIR)/opkda1.o $(BUILD_DIR)/opkda2.o
$(BUILD_DIR)/odepack.o: $(BUILD_DIR)/odepack-sub1.o $(BUILD_DIR)/odepack-sub2.o
$(BUILD_DIR)/IniFile.o: $(BUILD_DIR)/const.o
$(BUILD_DIR)/matrix_utils.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o
$(BUILD_DIR)/config.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/interactions.o \
	$(BUILD_DIR)/IniFile.o $(BUILD_DIR)/matrix_utils.o $(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/equations.o
$(BUILD_DIR)/interactions.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/matrix_utils.o $(BUILD_DIR)/utilities.o $(BUILD_DIR)/bspline_module.o \
	$(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/cosmology.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/interactions.o $(BUILD_DIR)/utilities.o $(BUILD_DIR)/bspline_module.o \
	$(BUILD_DIR)/linear_interpolation_module.o
$(BUILD_DIR)/equations.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/errors.o \
	$(BUILD_DIR)/cosmology.o $(BUILD_DIR)/interactions.o $(BUILD_DIR)/utilities.o \
	$(BUILD_DIR)/bspline_module.o $(BUILD_DIR)/linear_interpolation_module.o
#$(BUILD_DIR)/mc.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/config.o
#$(BUILD_DIR)/minimize.o: $(BUILD_DIR)/const.o $(BUILD_DIR)/config.o
$(BUILD_DIR)/nuDens.o: $(OBJ_FILES)
$(BUILD_DIR)/tests.o: $(OBJ_FILES)

all: nudens

directories:
	mkdir -p $(BUILD_DIR)

objects: $(OBJ_FILES) $(BUILD_DIR)/nuDens.o

nudens: directories objects Makefile
	$(F90) -o bin/$(EXECNAME) $(OBJ_FILES) $(BUILD_DIR)/nuDens.o $(FFLAGS)

nudens_debug: directories objects Makefile
	$(F90) -o bin/nuDens_debug.exe $(OBJ_FILES) $(BUILD_DIR)/nuDens.o $(FFLAGS)

tests: directories objects Makefile $(BUILD_DIR)/tests.o
	$(F90) -o bin/tests $(OBJ_FILES) $(BUILD_DIR)/tests.o $(DEBUGFLAGS)

clean: 
	rm -f bin/* $(BUILD_DIR)/*o $(BUILD_DIR)/*mod

$(BUILD_DIR)/%.o: sources/%.f90 Makefile
	$(F90) $(FFLAGS) -c sources/$*.f90 -o $(BUILD_DIR)/$*.o

$(BUILD_DIR)/%.o: sources/%.f Makefile
	$(F90) $(FFLAGS) $(stdFlag) -c sources/$*.f -o $(BUILD_DIR)/$*.o
