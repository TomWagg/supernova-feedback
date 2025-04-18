#!/bin/bash
## Job Name
#SBATCH --job-name=r442
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=22:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/r442_%a_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "base sim"

python /mnt/home/twagg/supernova-feedback/src/changa_feedback.py -p 128 -t 200
