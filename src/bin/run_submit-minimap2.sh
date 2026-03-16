#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --mem=32G
#
#
#

# load conda env
conda activate hla

## arguments
# REFERENCE = reference version 37, 38, t2t
# READS = full path of the reads
# OUT = prefix of the results file (bam)

PRESET="ont"
MM_PRESET="map-ont"
REFERENCE=$1
READS=$2
OUT="align.bam"
THREADS=16

minimap_version=$(minimap2 --version)
samtools_version=$(samtools --version | grep samtools | head -n 1 | cut -d " " -f 2)
echo "[INFO] mapping ${PRESET} reads:"
echo "       minimap2 version ${minimap_version}"
echo "       samtools version ${samtools_version}"

echo "[INFO]: CMD: minimap2 -ax ${MM_PRESET} -t ${THREADS} ${REFERENCE} ${READS} | samtools sort -m 2G - > ${OUT}"

minimap2 \
  -ax ${MM_PRESET} \
  -t ${THREADS} \
  ${REFERENCE} \
  ${READS} | samtools sort -m 2G - > ${OUT}

echo "[INFO] indexing ${OUT} (samtools)"
samtools index ${OUT}
