# FLAGS
#ifort
F90=ifort
F90FLAGS=-O3 -L/usr/lib -Ibuild/ -module build/ -g -p -traceback -no-prec-div
# -openmp
stdFlag=

#gfortran
#F90=gfortran
#F90FLAGS=-O3 -L/usr/lib -Jbuild/ -Ibuild/ -g -traceback -fopenmp
#stdFlat= --std=legacy

OBJ_FILES=build/const.o build/errors.o build/config.o build/IniFile.o \
	build/utilities.o build/matrix_utils.o \
	build/interactions.o build/cosmology.o build/equations.o build/functions.o \
	build/odepack.o build/odepack-sub1.o build/odepack-sub2.o \
	build/bspline_module.o build/bspline_oo_module.o build/bspline_sub_module.o
#build/opkdmain.o build/opkda1.o build/opkda2.o

default: all

build/%.o: sources/%.f90 Makefile
	$(F90) $(F90FLAGS) -c sources/$*.f90 -o build/$*.o

build/%.o: sources/%.f Makefile
	$(F90) $(F90FLAGS) $(stdFlag) -c sources/$*.f -o build/$*.o

build/bspline_module.o: build/bspline_oo_module.o build/bspline_sub_module.o
build/bspline_oo_module.o: build/bspline_sub_module.o
#build/opkdmain.o: build/opkda1.o build/opkda2.o
build/odepack.o: build/odepack-sub1.o build/odepack-sub2.o
build/IniFile.o: build/const.o
build/matrix_utils.o: build/const.o build/errors.o
build/config.o: build/const.o build/interactions.o build/IniFile.o build/matrix_utils.o build/errors.o
build/interactions.o: build/const.o build/errors.o build/matrix_utils.o build/utilities.o build/bspline_module.o
build/cosmology.o: build/const.o build/errors.o build/interactions.o build/utilities.o build/bspline_module.o
build/equations.o: build/const.o build/errors.o build/cosmology.o build/interactions.o build/utilities.o
build/functions.o: build/const.o build/errors.o build/cosmology.o build/interactions.o build/utilities.o build/equations.o
#build/mc.o: build/const.o build/config.o
#build/minimize.o: build/const.o build/config.o
build/nuDens.o: $(OBJ_FILES)

all: nudens

nudens: $(OBJ_FILES) build/nuDens.o Makefile
	$(F90) -o bin/nuDens.exe $(OBJ_FILES) build/nuDens.o $(F90FLAGS)

clean: 
	rm -f bin/* build/*o build/*mod
