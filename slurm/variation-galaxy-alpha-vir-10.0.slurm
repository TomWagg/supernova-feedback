#!/bin/bash
## Job Name
#SBATCH --job-name=alpha-vir-10.0
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=1:30:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/alpha-vir-10.0_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "alpha-vir-10.0"

python /mnt/home/twagg/supernova-feedback/src/variations-galaxy.py -f "alpha-vir-10.0.h5" -a 10.0