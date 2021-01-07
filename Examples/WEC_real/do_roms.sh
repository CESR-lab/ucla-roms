#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms sample_wec.in
echo "complete!"
