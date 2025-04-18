#!/bin/bash
## Job Name
#SBATCH --job-name=ce-0.1
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=15:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/ce-0.1_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "ce-0.1"

python /mnt/home/twagg/supernova-feedback/src/variations.py -f "ce-0.1.h5" -c 0.1