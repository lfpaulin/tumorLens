#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --mem=16G
#
#
#

#conda envs
# load conda env
conda activate hla

schlons_bin=$1
results_dir=$2
annot_bed=$3
reference=$4
sample_id=$5
hla_rare=$6
tumor_json=$7
control_json=$8
hla_annot_dir=$9

${schlons_bin} loss \
    --work-dir         "${results_dir}" \
    --annot            "${annot_bed}" \
    --fasta            "${reference}" \
    --output           "${sample_id}" \
    --hla-rare         "${hla_rare}" \
    --tumor            "${tumor_json}" \
    --control          "${control_json}" \
    --hla-annot        "${hla_annot_dir}" \
    --debug
