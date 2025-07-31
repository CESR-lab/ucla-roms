# If on NERSC's perlmutter system, use `source NERSC_perlmutter.sh` to configure your environment for ROMS
module restore
module load cpu/1.0
module load cray-hdf5/1.12.2.9
module load cray-netcdf/4.9.0.9

export MPIHOME=${CRAY_MPICH_PREFIX}/
export NETCDFHOME=${CRAY_NETCDF_PREFIX}/
export PATH=${PATH}:${NETCDFHOME}/bin
export LIBRARY_PATH=${LIBRARY_PATH}:${NETCDFHOME}/lib

export ROMS_ROOT=$(cd ../ && pwd) #Set parent directory to this file as ROMS_ROOT
export PATH=$PATH:$ROMS_ROOT/Tools-Roms
