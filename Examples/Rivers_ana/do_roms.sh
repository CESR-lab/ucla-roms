#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms river_ana.in
echo "complete!"
