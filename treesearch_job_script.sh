#!/bin/bash
#SBATCH --job-name=treesearch
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:0:0
#SBATCH --output=treesearch.out
#SBATCH --error=treesearch.err
#SBATCH --mail-user=millsbro@oregonstate.edu
#SBATCH --mail-type=END


#Bash commands to run within job
python PhylogenomicPipeline.py



