#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms USWC_flux.in
echo "complete!"
