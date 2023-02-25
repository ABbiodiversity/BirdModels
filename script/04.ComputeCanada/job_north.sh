#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --mem=0
#SBATCH --time=36:00:00
#SBATCH --job-name=north
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=ecknight@ualberta.ca

module load nixpkgs/16.09
module load gcc/7.3.0
module load openmpi/3.1.2
module load r/3.5.1

export NODESLIST=$(echo $(srun hostname))
Rscript --vanilla 04A.north.R
