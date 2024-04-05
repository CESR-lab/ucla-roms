#!/bin/bash

echo "run from initial"

mpirun -n 6 ./roms restart1.in < /dev/null >  jobout1

ncjoin -d restart_his.20121209140755.?.nc
ncjoin -d restart_his.20121209144115.?.nc
mv restart_his.20121209140755.nc restart_his1a.nc
mv restart_his.20121209144115.nc restart_his1b.nc

echo "run from restart"

mpirun -n 6 ./roms restart2.in < /dev/null >  jobout2
ncjoin -d restart_his.20121209140755.?.nc
ncjoin -d restart_his.20121209144115.?.nc
mv restart_his.20121209140755.nc restart_his2a.nc
mv restart_his.20121209144115.nc restart_his2b.nc

ncdiff -O restart_his1a.nc restart_his2a.nc dfa.nc
ncdiff -O restart_his1b.nc restart_his2b.nc dfb.nc

ncjoin -d restart_rst.20121209140755.?.nc

rm restart*.?.nc

echo "complete!"
