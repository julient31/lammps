#! /usr/bin/env bash

mkdir initial_spirals
cd initial_spirals
../src/create_spirals.py 100 4 
cd ..

mkdir run_cyclo
cd run_cyclo
../src/lammps_cyclo.py 100 4
cd ..
