CONDA_ENV=roms-ci

export ROMS_ROOT=.
export MPIHOME=${CONDA_PREFIX} #e.g. /Users/you/miniconda3/envs/roms_marbl
export NETCDFHOME=${CONDA_PREFIX}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NETCDFHOME/lib"
export PATH="./:$PATH"
export PATH=$PATH:$ROMS_ROOT/Tools-Roms
