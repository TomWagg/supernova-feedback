#!/bin/bash
## Job Name
#SBATCH --job-name=ccsn-20
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=15:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/ccsn-20_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "ccsn-20"

python /mnt/home/twagg/supernova-feedback/src/variations.py -f "ccsn-20.h5" -C 20