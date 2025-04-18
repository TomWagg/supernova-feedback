#!/bin/bash
## Job Name
#SBATCH --job-name=qcritB-0.0
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=15:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/qcritB-0.0_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "qcritB-0.0"

python /mnt/home/twagg/supernova-feedback/src/variations.py -f "qcritB-0.0.h5" -q 0.0001