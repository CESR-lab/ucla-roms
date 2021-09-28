#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms sample_wec_benchmark.in > benchmarks/wec/sample_wec_test.log
echo "complete!"
