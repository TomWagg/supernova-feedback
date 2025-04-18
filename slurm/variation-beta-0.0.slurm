#!/bin/bash
## Job Name
#SBATCH --job-name=beta-0.0
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=15:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/beta-0.0_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "beta-0.0"

python /mnt/home/twagg/supernova-feedback/src/variations.py -f "beta-0.0.h5" -b 0.0