#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --job-name=south
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=ecknight@ualberta.ca
#SBATCH --mail-type=ALL
module load StdEnv/2020
module load r/4.2.1
export NODESLIST=$(echo $(srun hostname))
Rscript --vanilla 04B.south.R
