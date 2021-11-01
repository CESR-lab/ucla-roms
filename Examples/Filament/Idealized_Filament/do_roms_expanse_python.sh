#!/bin/bash
#SBATCH --job-name="roms"

# %j=job_number and %N gives nodelist output="wec_real.%j.%N.out"
#SBATCH --output="python_ttw_log.%j.%N.run.out"
#SBATCH --partition=debug

# Can only use a max of 2 nodes on 'debug' partition:
#SBATCH --nodes=1

# Expanse has 128 cores per node:
#SBATCH --ntasks-per-node=1
#SBATCH --account=cla119
#SBATCH --export=ALL

# Max run time on 'debug' is 30 minutes:
#SBATCH -t 00:02:00

# Flags needed for mvapich2:
export MV2_USE_RDMA_CM=0
export MV2_IBA_HCA=mlx5_2
export MV2 DEFAULT PORT=1

module reset
module load gcc/10.2.0
module load python
module load py-numpy
module load py-cftime
module load netcdf-c
module load py-netcdf4

python3 setup_grid_init_bry.py params_fil.py



