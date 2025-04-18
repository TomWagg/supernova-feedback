#!/bin/bash
## Job Name
#SBATCH --job-name=alpha-vir-8.75
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=12:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/alpha-vir-8.75_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "alpha-vir-8.75"

python /mnt/home/twagg/supernova-feedback/src/variations-galaxy.py -f "alpha-vir-8.75.h5" -a 8.75