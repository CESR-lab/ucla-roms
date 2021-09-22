#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms river_ana_benchmark.in > benchmarks/river/river_ana_test.log
rm grid.*.nc
echo "complete!"
