#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:4
#SBATCH --time=20:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=solv-1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ca2356@nyu.edu
#SBATCH --output=slurm_%j.out

module purge
module load cuda/9.2.88
parallel ./run_solvation.sh ::: 0 1 2 3
