#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --mem=8G
#
#
#
#

#conda envs
# load conda env
conda activate hla

THREADS="2"
# BIN

SAMPLE_ID=$1
MODBAM_IN=$2
TYPE_BAM_IN=$3
HLABIN=$4

# TODO: extract reads names and position for read extension and extraction of methylation called reads
# TODO: split read methylation analysis
${HLABIN}/utils/split_read_methyl.py \
    --bam ${TYPE_BAM_IN}
    --modbam ${MODBAM_IN}
    --sample ${SAMPLE_ID}