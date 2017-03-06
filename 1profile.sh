#!/bin/bash
if [ "$1" == "" ]
then
	exe=nuDens
else
	exe=$1
fi
gprof bin/$exe.exe | gprof2dot | dot -Tpng -o output.png
