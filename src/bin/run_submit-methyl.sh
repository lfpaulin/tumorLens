#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --mem=32G
#
#
#
#OUT_ERR_CHDIR

# load conda env
conda activate hla

## arguments
REFERENCE=$1
BAM_ALN=$2
GENES_BED=$3
FINAL_DIR=$4
THREADS_USE=8
OUT="align_mod"

# extraction of modified bases with modkit
# genome-wide
modkit pileup \
    --cpg --ref "${REFERENCE}" \
    --threads ${THREADS_USE} \
    "${BAM_ALN}"  "${FINAL_DIR}/${OUT}_modkit.bed"
bgzip --force "${FINAL_DIR}/${OUT}_modkit.bed"
tabix --preset bed "${FINAL_DIR}/${OUT}_modkit.bed.gz"

# for some genes
modkit pileup \
    --cpg --ref "${REFERENCE}" \
    --include-bed "${GENES_BED}"  \
    --threads ${THREADS_USE} \
    "${BAM_ALN}"  "${FINAL_DIR}/${OUT}_modkit_genes.bed"
bgzip --force "${FINAL_DIR}/${OUT}_modkit_genes.bed"
tabix --preset bed "${FINAL_DIR}/${OUT}_modkit_genes.bed.gz"
