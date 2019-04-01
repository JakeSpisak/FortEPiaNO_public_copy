#!/bin/bash
if [ "$1" == "" ]
then
	exe=nuDens
else
	exe=$1
fi
gprof bin/$exe.exe | gprof2dot | dot -Tpng -o output.png

#other way:
#valgrind --tool=callgrind bin/nuDens.exe ...
#gprof2dot -f callgrind callgrind.out.X | dot -Tpng -o output.png
