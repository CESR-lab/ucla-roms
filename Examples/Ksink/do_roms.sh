#!/bin/bash

echo "run model..."
mpirun -n 6 ./roms ksink.in
ncjoin --delete ksink_his*.?.nc
ncjoin --delete ksink_dia*.?.nc
ncview ksink_his.*.nc &
echo "complete!"
