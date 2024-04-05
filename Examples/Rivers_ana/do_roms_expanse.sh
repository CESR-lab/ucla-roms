#!/bin/bash
#SBATCH --job-name="rivers_ana"
#SBATCH --output="rivers_ana.out"
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --account=edf100
#SBATCH --export=ALL
#SBATCH --mail-type=ALL
#SBATCH -t 00:10:00

echo "run model..."
srun --mpi=pmi2 -n 6 ./roms river_ana.in
echo "complete!"
