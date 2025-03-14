#!/bin/bash

echo "run model..."
mpiexec -n 6 ./roms cdr.in
echo "complete!"
