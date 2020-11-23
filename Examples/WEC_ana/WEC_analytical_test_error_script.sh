#!/bin/bash

echo "Running roms analytical test to ensure error values have not been altered..."
# remove wpp_his.0000.nc if it already exists else ncjoin creates funny files.
rm wpp_his.0000.nc
mpirun -n 8 ./roms
echo "Running ncjoin to join netcdf result files!"
rm grid*
ncjoin -d wpp_his.0000.*
echo "Running python error script: wave_packet_offline.py!"
python wave_packet_offline.py
echo "WEC analtyical test complete!"
read -rsp $'Press any key or wait a few seconds to continue...\n' -n 1 -t 5
