module reset
module load slurm sdsc DefaultModules shared
module load cpu/0.15.4 intel/19.1.1.217 mvapich2/2.3.4
module load netcdf-c/4.7.4

export NETCDFHOME=${NETCDF_FORTRANHOME}/
export MPIHOME=${MVAPICH2HOME}/
export MPIROOT=${MVAPICH2HOME}/

export ROMS_ROOT=$(cd ../ && pwd)


