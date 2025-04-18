#!/bin/bash
## Job Name
#SBATCH --job-name=ce-10.0
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=15:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/ce-10.0_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "ce-10.0"

python /mnt/home/twagg/supernova-feedback/src/variations.py -f "ce-10.0.h5" -c 10.0