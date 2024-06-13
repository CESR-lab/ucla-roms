#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms bgc.in
echo "complete!"
