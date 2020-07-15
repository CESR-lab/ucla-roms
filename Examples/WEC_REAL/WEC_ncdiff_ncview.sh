#!/bin/bash

echo "run model & netcdf delete, ncjoin, ncdiff, ncview..."
make
mpirun -n 6 ./roms
rm usw3_testcase.0000.nc diff_usw3_f48M_vs_f47M_ncdf.nc
ncjoin -d usw3_testcase.0000.*
ncdiff usw3_testcase.0000.nc /home/ddevin/WEC_maya/0611_USWC_WEC/47M_Makdefs_FFLAGS_oldcode/usw3_testcase.0000.nc diff_usw3_f48M_vs_f47M_ncdf.nc
ncview diff_usw3_f48M_vs_f47M_ncdf.nc
echo "netcdf delete, ncjoin and ncdiff complete!"
