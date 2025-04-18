#!/bin/bash
## Job Name
#SBATCH --job-name=no-fallback
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=15:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/no-fallback_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "no-fallback"

python /mnt/home/twagg/supernova-feedback/src/variations.py -f "no-fallback.h5" -k 3