#!/bin/bash
#SBATCH --job-name="Rivers_real"
#SBATCH --output="Rivers_real.out"
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --account=edf100
#SBATCH --export=ALL
#SBATCH --mail-type=ALL
#SBATCH -t 00:10:00

source ~/.bashrc

echo "run model..."
srun --mpi=pmi2 -n 6 ./roms rivers.in
echo "complete!"
