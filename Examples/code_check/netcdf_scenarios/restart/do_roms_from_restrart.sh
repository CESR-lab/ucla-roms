#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms sample_from_restart.in
echo "complete!"
