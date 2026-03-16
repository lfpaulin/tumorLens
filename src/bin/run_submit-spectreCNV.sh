#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=4G
#
#
#

#conda envs
# load conda env
conda activate hla

# Params
SAMPLE_ID=$1
COVERAGE_FILE=$2
FASTA_REF=$3
REF_META=$4
BLACKLIST=$5
SNV=$6
TUMOR_PURITY=$7
CANCER=$8

echo "[INFO] Starting SPECTRE"
if [[ "" == "${CANCER}" ]];
then
    echo "[INFO] single sample mode"
    spectre CNVCaller \
        --sample-id    "${SAMPLE_ID}" \
        --coverage     "${COVERAGE_FILE}" \
        --snv          "${SNV}" \
        --output-dir   . \
        --reference    "${FASTA_REF}" \
        --metadata     "${REF_META}" \
        --blacklist    "${BLACKLIST}"
else
    echo "[INFO] cancer mode"
    COVERAGE_FILES=${COVERAGE_FILE//,/ }
    SNV_FILES=${SNV//,/ }
    spectre Cancer \
        --sample-id     "${SAMPLE_ID}" \
        --coverage      "${COVERAGE_FILES}" \
        --snv           "${SNV_FILES}" \
        --output-dir    . \
        --reference     "${FASTA_REF}" \
        --metadata      "${REF_META}" \
        --blacklist     "${BLACKLIST}" \
        --tumor-content "${TUMOR_PURITY}" \
        --ploidy 2
#        --ploidy-chr chrX:1,chrY:1 \
fi
