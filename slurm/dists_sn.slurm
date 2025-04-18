#!/bin/bash
## Job Name
#SBATCH --job-name=sn_dists
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=1:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/logs_%a_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "Getting supernova distances"

python /mnt/home/twagg/supernova-feedback/src/get_supernova_distances.py
