#!/bin/bash
lab=$1
prec=$2
python python/prepareIni.py ini/3p1/15_${lab}.ini OUT/prec_3p1/15_${lab}/ 3+1 damping --dlsoda_rtol=$prec --dlsoda_atol=$prec --Nx=2000 --Ny=15 --Nylog=1 --y_cen=0.01
python python/prepareIni.py ini/3p1/20_${lab}.ini OUT/prec_3p1/20_${lab}/ 3+1 damping --dlsoda_rtol=$prec --dlsoda_atol=$prec --Nx=2000 --Ny=20 --Nylog=1 --y_cen=0.01
python python/prepareIni.py ini/3p1/30_${lab}.ini OUT/prec_3p1/30_${lab}/ 3+1 damping --dlsoda_rtol=$prec --dlsoda_atol=$prec --Nx=2000 --Ny=30 --Nylog=3 --y_cen=1
python python/prepareIni.py ini/3p1/40_${lab}.ini OUT/prec_3p1/40_${lab}/ 3+1 damping --dlsoda_rtol=$prec --dlsoda_atol=$prec --Nx=2000 --Ny=40 --Nylog=5 --y_cen=1
python python/prepareIni.py ini/3p1/70_${lab}.ini OUT/prec_3p1/70_${lab}/ 3+1 damping --dlsoda_rtol=$prec --dlsoda_atol=$prec --Nx=2000 --Ny=70 --Nylog=10 --y_cen=1
for ny in 15 20 30 40 70;
do
	bin/nuDens.exe ini/3p1/${ny}_${lab}.ini
done
