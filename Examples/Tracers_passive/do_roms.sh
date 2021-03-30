#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms sample_tracers.in
echo "complete!"
