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
REF_FAS=$1
BASEDIR=$2
TUMOR_DIR=$3
TUMOR_BAM="${BASEDIR}/${TUMOR_DIR}/methylation/align_mod.bam"
TUMOR_SNV="${BASEDIR}/${TUMOR_DIR}/snv/snv.vcf.gz"
TUMOR_SNV_PHASED="${BASEDIR}/${TUMOR_DIR}/snv/snv_phased.vcf.gz"

CONTROL_DIR=$4
CONTROL_BAM="${BASEDIR}/${CONTROL_DIR}/methylation/align_mod.bam"
CONTROL_SNV="${BASEDIR}/${CONTROL_DIR}/snv/snv.vcf.gz"
CONTROL_SNV_PHASED="${BASEDIR}/${CONTROL_DIR}/snv/snv_phased.vcf.gz"

TUMOR_ID=$5
CONTROL_ID=$6
BED_EXTRACT=$7
MULTI_DIR=$8

SAMPLE_BASE=$(echo "${TUMOR_DIR}" | cut -d "/" -f 1,2)
if [[ "" != "${MULTI_DIR}" ]];
then
    RESULTS_DIR="phase_methyl_${MULTI_DIR}"
else
    RESULTS_DIR="phase_methyl"
fi

echo "${BASEDIR}/${SAMPLE_BASE}/${RESULTS_DIR}"
if [[ ! -d "${BASEDIR}/${SAMPLE_BASE}/${RESULTS_DIR}" ]];
then
    mkdir "${BASEDIR}/${SAMPLE_BASE}/${RESULTS_DIR}"
fi

TUMOR_BAM_GENES="${BASEDIR}/${SAMPLE_BASE}/${RESULTS_DIR}/${TUMOR_ID}.bam"
CONTROL_BAM_GENES="${BASEDIR}/${SAMPLE_BASE}/${RESULTS_DIR}/${CONTROL_ID}.bam"
TUMOR_BAM_GENES_HP="${BASEDIR}/${SAMPLE_BASE}/${RESULTS_DIR}/${TUMOR_ID}_hp.bam"
CONTROL_BAM_GENES_HP="${BASEDIR}/${SAMPLE_BASE}/${RESULTS_DIR}/${CONTROL_ID}_hp.bam"

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

cd "${BASEDIR}/${SAMPLE_BASE}/${RESULTS_DIR}/" || exit

# bed, bamin, bamout
if [[ ! -f "${TUMOR_BAM_GENES}" ]];
then
    echo "Reads extraction tumor"
    extract_bam_reads_bed.sh  "${BED_EXTRACT}"  "${TUMOR_BAM}"    "${TUMOR_BAM_GENES}"
else
    echo "${TUMOR_BAM_GENES} exists, skipped"
fi

if [[ ! -f "${CONTROL_BAM_GENES}" ]];
then
    echo "Reads extraction control"
    extract_bam_reads_bed.sh  "${BED_EXTRACT}"  "${CONTROL_BAM}"  "${CONTROL_BAM_GENES}"
else
    echo "${CONTROL_BAM_GENES} exists, skipped"
fi

# cleaning
if [[ -f "tmp.sam" ]];       then  rm tmp.sam; fi
if [[ -f "tmp_rmdup.bam" ]]; then  rm tmp_rmdup.bam; fi

# haplotagging
if [[ ! -f "${CONTROL_BAM_GENES_HP}" ]];
then
    if [[ -f "${CONTROL_SNV_PHASED}" && -f "${CONTROL_BAM_GENES}" ]];
    then
        echo "Haplotag control"
        whatshap haplotag \
            -o "${CONTROL_BAM_GENES_HP}" \
            --output-threads ${THREADS} \
            --reference "${REF_FAS}" \
            --ignore-read-groups \
            "${CONTROL_SNV_PHASED}"  "${CONTROL_BAM_GENES}"
    else
        echo "Missing input ${CONTROL_SNV_PHASED} or ${CONTROL_BAM_GENES}"
    fi
else
    echo "${CONTROL_BAM_GENES_HP} exists, skipped"
fi

if [[ ! -f "${TUMOR_BAM_GENES_HP}" ]];
then
    if [[ -f "${TUMOR_SNV_PHASED}" && -f "${TUMOR_BAM_GENES}" ]];
    then
        echo "Haplotag tumor"
        whatshap haplotag \
            -o "${TUMOR_BAM_GENES_HP}" \
            --reference "${REF_FAS}" \
            --output-threads ${THREADS} \
            --ignore-read-groups \
            "${TUMOR_SNV_PHASED}"  "${TUMOR_BAM_GENES}"
    else
        echo "Missing input ${TUMOR_SNV_PHASED} or ${TUMOR_BAM_GENES}"
    fi
else
    echo "${TUMOR_BAM_GENES_HP} exists, skipped"
fi
