#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms rivers.in
echo "complete!"
