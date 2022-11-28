#!/bin/bash

prefix=/home/jakespisak/software/fortepiano_public

for i in {3..9};
    do
        for j in {0..2};
            do 
                echo "i is" $i
                echo "j is" $j
                idx=$(( 3*$i + $j))
                mkdir -p $prefix/output/equal_mixing_extended/run$idx
                cd $prefix
                nohup $prefix/bin/fortepiano $prefix/ini/equal_mixing_extended/run$idx.ini &> $prefix/output/equal_mixing_extended/run$idx/nohup$idx.out &
            done
        sleep 1h
    done
    
