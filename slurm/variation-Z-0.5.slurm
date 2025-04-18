#!/bin/bash
## Job Name
#SBATCH --job-name=Z-0.5
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=15:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/Z-0.5_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "Z-0.5"

python /mnt/home/twagg/supernova-feedback/src/variations.py -f "Z-0.5.h5" -Z 0.5