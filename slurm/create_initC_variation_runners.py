basic = """#!/bin/bash
## Job Name
#SBATCH --job-name=NAMEHERE
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=20:00:00
#SBATCH -o /mnt/home/twagg/supernova-feedback/slurm/logs/NAMEHERE_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tomwagg@uw.edu
#SBATCH --export=all

source /mnt/home/twagg/.bashrc
conda activate cogsworth

echo "NAMEHERE"

python /mnt/home/twagg/supernova-feedback/src/feedback_simulation.py -p 128 -t 200 -f "NAMEHERE.h5" -FLAGHERE VALHERE"""

for name, flag, value in zip(["porb-minus1", "porb-0", "q-minus1", "q-plus1", "singles", "imf-1.9", "imf-2.7", "porb-max-3"],
                             ['P', 'P', 'q', 'q', "b", "m", "m", "M"],
                             ['-0.999', '0', '-0.999', '1', "0", "-1.9", "-2.7", "3.0"]):
    with open(f"variation-initC-{name}.slurm", "w") as f:
        write_this = basic.replace("NAMEHERE", name).replace("FLAGHERE", flag).replace("VALHERE", str(value))
        f.write(write_this)
