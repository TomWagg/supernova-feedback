#!/bin/bash
## Job Name
#SBATCH --job-name=ecsn-265
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=15:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/ecsn-265_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "ecsn-265"

python /mnt/home/twagg/supernova-feedback/src/variations.py -f "ecsn-265.h5" -e -265