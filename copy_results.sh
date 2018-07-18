#!/bin/bash

if test ! -d results; then
    mkdir results
fi

for lattice in 11.0 11.1 11.2 # 11.3 11.4 11.5 11.5014556585 11.6 11.7 11.8 11.9 12.0
do
    if test -f $lattice/relax.out;
    then
        if test ! -d results/$lattice;
        then
            mkdir results/$lattice
        fi
        cp $lattice/*.in results/$lattice/
        cp $lattice/*.out results/$lattice/
        cp $lattice/*.bands.* results/$lattice/
    fi
done
