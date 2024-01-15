#!/bin/bash -l
#SBATCH -p normal
#SBATCH -J Snake_sub
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --tasks 1
. /etc/profile

conda activate /mnt/home/lkirsch/anaconda3/envs/phylofisher

python3 run_pipeline.py toRUN3.txt