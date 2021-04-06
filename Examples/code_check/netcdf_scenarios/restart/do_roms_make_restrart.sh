#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms sample_make_restart.in
echo "complete!"
