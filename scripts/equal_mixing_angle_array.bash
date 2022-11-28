#!/bin/bash

prefix=/home/jakespisak/software/fortepiano_public

for i in {3..29}; 
    do 
        mkdir -p $prefix/output/equal_mixing_angle_array/run$i
        cd $prefix
        nohup $prefix/bin/fortepiano $prefix/ini/equal_mixing_angle_array/run$i.ini &> $prefix/output/equal_mixing_angle_array/run$i/nohup$i.out &
        sleep 15m
    done
    
