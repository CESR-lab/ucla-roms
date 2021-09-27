#!/bin/bash

echo "run model..."

cd ../Flux_frc
make
mpirun -n 6 ./roms sample.in
../Lmd_kpp_old/join.sh

cd ../Lmd_kpp_old
make
mpirun -n 6 ./roms sample.in
./join.sh

ncdiff sample_his.0000.nc ../Flux_frc/sample_his.0000.nc diff.nc
ncview diff.nc

echo "complete!"
