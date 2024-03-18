#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms extract.in
echo "complete!"
