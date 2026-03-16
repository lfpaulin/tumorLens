#!/bin/bash
#SBATCH --ntasks=4
#SBATCH --mem=12G
#
#
#

#conda envs
# load conda env
conda activate clair

# tasks
THREADS="4"
# move to annot
MODEL_GUPPY6="data/guppy_models/clair3_model-guppy_6"
REF_FAS=$3

input_bam=$1
chunk=$2
output_dir="snv_$(basename ${chunk} | cut -d "." -f 1)"

run_clair3.sh \
    --ref_fn="${REF_FAS}" \
    --model_path="${MODEL_GUPPY6}" \
    --threads="${THREADS}"  --platform="ont" \
    --bam_fn="${input_bam}" \
    --output="${output_dir}" \
    --bed_fn="${chunk}"
