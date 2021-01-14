#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms wec_netcdf_benchmark.in > benchmarks/netcdf/wec_netcdf_test.log
rm sample_wec*.nc
echo "complete!"
