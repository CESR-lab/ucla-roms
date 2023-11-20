#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms passive.in
echo "complete!"
