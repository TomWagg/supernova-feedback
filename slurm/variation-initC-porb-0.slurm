#!/bin/bash
## Job Name
#SBATCH --job-name=porb-0
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=20:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/porb-0_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "porb-0"

python /mnt/home/twagg/supernova-feedback/src/feedback_simulation.py -p 128 -t 200 -f "porb-0.h5" -P 0