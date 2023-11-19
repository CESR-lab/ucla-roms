#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms Flux_frc.in
echo "complete!"
