#!/bin/bash
#SBATCH --ntasks=4
#SBATCH --mem=16G
#
#
#

#conda envs
# load conda env
conda activate hla

THREADS="4"
# Params
#   FILE_IN_BAM
FILE_IN_BAM=$1
WINSIZE=$2
# WINSIZE="1000"

echo "[INFO] Calc depth with mosdepth"
echo "[INFO] CMD: mosdepth --by ${WINSIZE} --no-per-base --mapq 20 coverage ${FILE_IN_BAM}"
mosdepth \
    --by ${WINSIZE} \
    --threads ${THREADS} \
    --no-per-base \
    --mapq 10 \
    coverage \
    ${FILE_IN_BAM}

# usage:
# sbatch --chdir="workdir" --output="log-mosdepth.out" --error="log-mosdepth.err" myjob.sh
