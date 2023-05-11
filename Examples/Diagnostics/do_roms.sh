#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms diags.in
echo "complete!"
