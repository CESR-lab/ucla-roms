#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms pipes.in
echo "complete!"
