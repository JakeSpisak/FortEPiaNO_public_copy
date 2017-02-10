# FLAGS
F90=gfortran
F90FLAGS=-O3 -L/usr/lib -Jbuild/ -Ibuild/

OBJ_FILES=build/const.o build/config.o build/IniFile.o build/matrix_utils.o

default: all

build/opk%.o: sources/%.f90 Makefile
	$(F90) $(F90FLAGS) -c sources/$*.f -o build/$*.o

build/%.o: sources/%.f90 Makefile
	$(F90) $(F90FLAGS) -c sources/$*.f90 -o build/$*.o

build/IniFile.o: build/const.o
build/matrix_utils.o: build/const.o
build/config.o: build/const.o build/IniFile.o build/matrix_utils.o
#build/data.o: build/const.o build/config.o
#build/flux.o: build/const.o build/config.o build/mc.o
#build/likelihood.o: build/const.o build/config.o
#build/mc.o: build/const.o build/config.o
#build/minimize.o: build/const.o build/config.o
build/nuDens.o: $(OBJ_FILES)

all: nudens

nudens: $(OBJ_FILES) build/nuDens.o Makefile
	$(F90) -o bin/nuDens.exe $(OBJ_FILES) build/nuDens.o $(F90FLAGS)

clean: 
	rm -f bin/* sources/*o sources/*mod
