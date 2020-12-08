#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms sample_pipe.in
echo "complete!"
