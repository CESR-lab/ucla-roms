#!/bin/bash

echo "run model & netcdf delete, ncjoin, ncview..."
mpirun -n 6 ./roms USWC_sample.in
rm uswc_his.0000.nc uswc_surf_flux.0000.nc
ncjoin -d uswc_his.0000.*.nc
ncjoin -d uswc_surf_flux.*.nc
ncview uswc_his.0000.nc &
echo "netcdf delete, ncjoin and ncview complete!"
