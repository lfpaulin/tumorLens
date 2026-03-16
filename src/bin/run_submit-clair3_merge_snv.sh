#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --mem=16G
#
#
#

#conda envs
# load conda env
conda activate hla

SNV_OUTPUT=$1
FILES_MERGE=$2

bcftools concat \
    --allow-overlaps \
    --remove-duplicates \
    --file-list "${FILES_MERGE}" \
    --threads 8 | \
      bcftools view --apply-filters PASS --type snps --include "AF > 0.1" | \
      bcftools sort | bgzip -c > ${SNV_OUTPUT}

tabix --force ${SNV_OUTPUT}
