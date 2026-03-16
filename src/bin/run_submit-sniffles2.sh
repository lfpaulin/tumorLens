#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --mem=16G
#
#
#

#conda envs
# load conda env
conda activate hla

THREADS="8"
# Params
#   FILE_IN_BAM
#   FASTA_REF
FILE_IN_BAM=$1
FASTA_REF=$2
SAMPLE_ID=$3

FILE_SV="sv.vcf.gz"
FILE_SNF="sv.snf"
FILE_SV_MOSAIC="sv_mosaic.vcf.gz"
FILE_SNF_MOSAIC="sv_mosaic.snf"
if [[ ! -f "${FILE_SV}" ]]; then
    echo "[INFO] Running sniffles2 ..."
    echo "[INFO] CMD: sniffles2 --input ${FILE_IN_BAM} --vcf ${FILE_SV} --snf ${FILE_SNF} --reference ${FASTA_REF} --threads ${THREADS} --allow-overwrite --sample-id ${SAMPLE_ID}"
    # germline
    sniffles --input "${FILE_IN_BAM}" \
        --vcf "${FILE_SV}"  \
        --snf "${FILE_SNF}"  \
        --reference "${FASTA_REF}"  \
        --threads ${THREADS} \
        --allow-overwrite --sample-id "${SAMPLE_ID}"
    # mosaic
    sniffles --input "${FILE_IN_BAM}" \
        --vcf "${FILE_SV_MOSAIC}"  \
        --snf "${FILE_SNF_MOSAIC}"  \
        --reference "${FASTA_REF}"  \
        --threads ${THREADS} \
        --mosaic \
        --allow-overwrite --sample-id "${SAMPLE_ID}"
else
    echo "[INFO] File '${FILE_SV}' already exists: skipping sniffles"
fi
