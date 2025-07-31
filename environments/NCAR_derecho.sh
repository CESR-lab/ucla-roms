module restore
module load ncarenv/23.09
module load intel-oneapi/2024.2.1
module load cray-mpich/8.1.29
module load netcdf/4.9.2

export MPIHOME=${CRAY_MPICH_PREFIX}/
export NETCDFHOME=${NETCDF}/
export LIBRARY_PATH=${LIBRARY_PATH}:${NETCDFHOME}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${NETCDFHOME}/lib

export ROMS_ROOT=$(cd ../ && pwd)
export PATH=$PATH:$ROMS_ROOT/Tools-Roms
