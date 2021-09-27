#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms sample.in
echo "complete!"
