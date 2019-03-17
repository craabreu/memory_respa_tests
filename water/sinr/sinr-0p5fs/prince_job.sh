#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:2
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=langevin
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ca2356@nyu.edu
#SBATCH --output=slurm_%j.out
 
module purge
module load cuda/9.2.88
python simulate.py --nsteps 7200000 --device 0 --secdev 1

