#!/bin/bash
## Job Name
#SBATCH --job-name=q-plus1
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=20:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/q-plus1_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "q-plus1"

python /mnt/home/twagg/supernova-feedback/src/feedback_simulation.py -p 128 -t 200 -f "q-plus1.h5" -q 1