#!/bin/bash
make BUILD_DIR=buildpar/ EXECNAME=parallel.exe

export OMP_NUM_THREADS=8
ini=ini/debug.ini
rm -r `grep outputFolder $ini |sed "s#outputFolder = ##g"`/*
sleep 1s
./bin/parallel.exe $ini
