#!/bin/bash

echo "run model..."
mpiexec -n 6 ./roms rivers.in
echo "complete!"
