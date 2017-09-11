#!/bin/bash

#SBATCH --job-name=xy
#SBATCH --output=xy.out
#SBATCH --array=24,40,48,56,72

#SBATCH --time=0-7:30

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --mem=32000M

#SBATCH --mail-user=mbeach@phas.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load julia

julia -p 64 generate-xy.jl $SLURM_ARRAY_TASK_ID
