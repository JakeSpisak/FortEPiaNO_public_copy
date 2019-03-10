#!/bin/bash
lab=$1
prec=$2
python python/prepareIni.py ini/3nu/20_${lab}.ini OUT/prec_3nu/20_${lab}/ 3nu damping --dlsoda_rtol=$prec --dlsoda_atol=$prec --Ny=20 --Nylog=1 --y_cen=0.01 --x_in=0.05
python python/prepareIni.py ini/3nu/40_${lab}.ini OUT/prec_3nu/40_${lab}/ 3nu damping --dlsoda_rtol=$prec --dlsoda_atol=$prec --Ny=40 --Nylog=5 --y_cen=1 --x_in=0.05
python python/prepareIni.py ini/3nu/70_${lab}.ini OUT/prec_3nu/70_${lab}/ 3nu damping --dlsoda_rtol=$prec --dlsoda_atol=$prec --Ny=70 --Nylog=10 --y_cen=1 --x_in=0.05
python python/prepareIni.py ini/3nu/100_${lab}.ini OUT/prec_3nu/100_${lab}/ 3nu damping --dlsoda_rtol=$prec --dlsoda_atol=$prec --Ny=100 --Nylog=10 --y_cen=1 --x_in=0.05
for ny in 20 40 70 100;
do
	bin/nuDens.exe ini/3nu/${ny}_${lab}.ini
done
