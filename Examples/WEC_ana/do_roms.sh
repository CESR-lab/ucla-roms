#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms USWC_wec_ana.in
echo "complete!"
