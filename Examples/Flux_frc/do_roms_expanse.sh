#!/bin/bash

if [ -z "${ACCOUNT_KEY}" ];then
    echo "ACCOUNT_KEY environment variable empty. Set using export ACCOUNT_KEY=<your_account_key>."
    exit 1
fi

sbatch --job-name="Flux_frc" \
       --output="Flux_frc.out" \
       --partition="debug" \
       --nodes=1 \
       --ntasks-per-node=6 \
       --account=${ACCOUNT_KEY} \
       --export=ALL \
       --mail-type=ALL \
       -t 00:10:00 \
       --wrap="srun --mpi=pmi2 -n 6 ./roms Flux_frc.in"
