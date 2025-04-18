basic = """#!/bin/bash
## Job Name
#SBATCH --job-name=NAMEHERE
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=12:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/NAMEHERE_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "NAMEHERE"

python /mnt/home/twagg/supernova-feedback/src/variations-galaxy.py -f "NAMEHERE.h5" -FLAGHERE VALHERE"""

for name, flag, value in zip(["alpha-vir-0.0875", "alpha-vir-8.75"],
                             ['a', 'a'],
                             [0.0875, 8.75]):
    with open(f"variation-galaxy-{name}.slurm", "w") as f:
        write_this = basic.replace("NAMEHERE", name).replace("FLAGHERE", flag).replace("VALHERE", str(value))
        f.write(write_this)
