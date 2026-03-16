#!/bin/bash
#SBATCH --ntasks=4
#SBATCH --mem=16G
#
#
#

#conda envs
# load conda env
conda activate phase

# tasks
THREADS="4"
# parameters
REF_FAS=${1}
BED_EXTRACT=${2}

TUMOR_DIR=${3}
TUMOR_BAM="${TUMOR_DIR}/methylation/align_mod.bam"
TUMOR_SNV="${TUMOR_DIR}/snv/snv.vcf.gz"
TUMOR_SNV_PHASED="${TUMOR_DIR}/snv/snv_phased.vcf.gz"
TUMOR_BAM_GENES="${TUMOR_DIR}/methylation/align_mod_genes.bam"
TUMOR_BAM_GENES_HP="${TUMOR_DIR}/methylation/align_mod_genes_hp.bam"

CONTROL_DIR=${4}
CONTROL_BAM="${CONTROL_DIR}/methylation/align_mod.bam"
CONTROL_SNV="${CONTROL_DIR}/snv/snv.vcf.gz"
CONTROL_SNV_PHASED="${CONTROL_DIR}/snv/snv_phased.vcf.gz"
CONTROL_BAM_GENES="${CONTROL_DIR}/methylation/align_mod_genes.bam"
CONTROL_BAM_GENES_HP="${CONTROL_DIR}/methylation/align_mod_genes_hp.bam"

SOURCE_DIR=${5}

# phase control with BAM+SNV
    if [[ ! -f "${CONTROL_SNV_PHASED}" ]];
then
    echo "Running control phasing"
    whatshap phase \
        -o "${CONTROL_SNV_PHASED}" \
        --reference="${REF_FAS}" \
        --ignore-read-groups \
        "${CONTROL_SNV}"  "${CONTROL_BAM}"
else
    echo "${CONTROL_SNV_PHASED} exists, skipped"
fi

# phase tumor with BAM+phasedControl
if [[ ! -f "${TUMOR_SNV_PHASED}" ]];
then
    echo "Running tumor phasing"
    whatshap phase \
        -o "${TUMOR_SNV_PHASED}" \
        --reference="${REF_FAS}" \
        --ignore-read-groups \
        "${TUMOR_SNV}"  "${CONTROL_SNV_PHASED}"
else
    echo "${TUMOR_SNV_PHASED} exists, skipped"
fi

# indexing
if [[ ! -f "${CONTROL_SNV_PHASED}.tbi" ]];
then
    echo "Running control indexing"
    bcftools index --tbi "${CONTROL_SNV_PHASED}"
else
    echo "${CONTROL_SNV_PHASED} index exists, skipped"
fi

if [[ ! -f "${TUMOR_SNV_PHASED}.tbi" ]];
then
    echo "Running tumor indexing"
    bcftools index --tbi "${TUMOR_SNV_PHASED}"
else
    echo "${TUMOR_SNV_PHASED} index exists, skipped"
fi

# extract reads control
cd "${CONTROL_DIR}/methylation" || exit
if [[ ! -f "${CONTROL_BAM_GENES}" ]];
then
    echo "Reads extraction control"
    bash ${SOURCE_DIR}/bin/extract_bam_reads_bed.sh  "${BED_EXTRACT}"  "${CONTROL_BAM}"  "${CONTROL_BAM_GENES}"
    # cleaning
    if [[ -f "tmp.sam" ]];       then  rm tmp.sam; fi
    if [[ -f "tmp_rmdup.bam" ]]; then  rm tmp_rmdup.bam; fi
else
    echo "${CONTROL_BAM_GENES} exists, skipped"
fi

# extract reads tumor
cd "${TUMOR_DIR}/methylation" || exit
if [[ ! -f "${TUMOR_BAM_GENES}" ]];
then
    echo "Reads extraction tumor"
    bash ${SOURCE_DIR}/bin/extract_bam_reads_bed.sh  "${BED_EXTRACT}"  "${TUMOR_BAM}"    "${TUMOR_BAM_GENES}"
    # cleaning
    if [[ -f "tmp.sam" ]];       then  rm tmp.sam; fi
    if [[ -f "tmp_rmdup.bam" ]]; then  rm tmp_rmdup.bam; fi
else
    echo "${TUMOR_BAM_GENES} exists, skipped"
fi

# haplotagging
if [[ ! -f "${CONTROL_BAM_GENES_HP}" ]];
then
    cd "${CONTROL_DIR}/snv" || exit
    if [[ -f "${CONTROL_SNV_PHASED}" && -f "${CONTROL_BAM_GENES}" ]];
    then
        echo "Haplotag control"
        whatshap haplotag \
            -o "${CONTROL_BAM_GENES_HP}" \
            --output-threads ${THREADS} \
            --reference "${REF_FAS}" \
            --ignore-read-groups \
            "${CONTROL_SNV_PHASED}"  "${CONTROL_BAM_GENES}"
        samtools index "${CONTROL_BAM_GENES_HP}"
    else
        echo "Missing input ${CONTROL_SNV_PHASED} or ${CONTROL_BAM_GENES}"
    fi
else
    echo "${CONTROL_BAM_GENES_HP} exists, skipped"
fi

if [[ ! -f "${TUMOR_BAM_GENES_HP}" ]];
then
    cd "${TUMOR_DIR}/snv" || exit
    if [[ -f "${TUMOR_SNV_PHASED}" && -f "${TUMOR_BAM_GENES}" ]];
    then
        echo "Haplotag tumor"
        whatshap haplotag \
            -o "${TUMOR_BAM_GENES_HP}" \
            --reference "${REF_FAS}" \
            --output-threads ${THREADS} \
            --ignore-read-groups \
            "${TUMOR_SNV_PHASED}"  "${TUMOR_BAM_GENES}"
        samtools index "${TUMOR_BAM_GENES_HP}"
    else
        echo "Missing input ${TUMOR_SNV_PHASED} or ${TUMOR_BAM_GENES}"
    fi
else
    echo "${TUMOR_BAM_GENES_HP} exists, skipped"
fi
