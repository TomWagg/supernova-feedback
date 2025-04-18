#!/bin/bash
## Job Name
#SBATCH --job-name=cogsworth
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=15:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/logs_%a_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "base sim"

python /mnt/home/twagg/supernova-feedback/src/feedback_simulation.py -p 128 -t 200 -f "fiducial"
