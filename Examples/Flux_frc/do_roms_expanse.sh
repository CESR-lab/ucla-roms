#!/bin/bash
#SBATCH --job-name="Flux_frc"
#SBATCH --output="Flux_frc.out"
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --account=edf100
#SBATCH --export=ALL
#SBATCH --mail-type=ALL
#SBATCH -t 00:10:00

source ~/.bashrc

echo "run model..."
srun --mpi=pmi2 -n 6 ./roms Flux_frc.in
echo "complete!"
