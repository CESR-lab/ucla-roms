#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms ksink.in
echo "complete!"
