#!/bin/bash

if [ -z "${ACCOUNT_KEY}" ];then
    echo "ACCOUNT_KEY environment variable empty. Set using export ACCOUNT_KEY=<your_account_key>."
    exit 1
fi


# Flags needed for mvapich2:
export MV2_USE_RDMA_CM=0
export MV2_IBA_HCA=mlx5_2
export MV2 DEFAULT PORT=1

module purge
module load slurm
module load cpu/0.15.4  intel/19.1.1.217  mvapich2/2.3.4
module load netcdf-c/4.7.4
module load netcdf-fortran/4.5.3

sbatch --job-name="Rivers_ana" \
       --output="Rivers_ana.out" \
       --partition="debug" \
       --nodes=1 \
       --ntasks-per-node=6 \
       --account=${ACCOUNT_KEY} \
       --export=ALL \
       --mail-type=ALL \
       -t 00:10:00 \
       --wrap="srun --mpi=pmi2 -n 6 roms ./river_ana.in"
