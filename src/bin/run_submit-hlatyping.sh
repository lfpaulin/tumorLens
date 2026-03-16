#!/bin/bash
#SBATCH --job-name=specHLA
#SBATCH --ntasks=8
#SBATCH --mem=32G
#
#

handle_error() {
    echo "An error occurred on line $1"
    exit 1
}

trap 'handle_error $LINENO' ERR

#conda envs
# load conda env
conda activate base

THREADS="8"
WORKDIR=$1    # sample_dir/typing
BAM_IN=$2     # aligned bam input
SPEC_HLA_BIN="bin/SpecHLA/script/long_read_typing.py"

# const
HLA_REGION="chr6:26000000-36000000"
FQ_HLA="typing_hla_region.fastq.gz"


#Required arguments:
#  -r        Long-read fastq file. PacBio or Nanopore. (default: None)
#  -n        Sample ID (default: None) <== they make a dir in the output dir, so results
#  -o        The output folder to store the typing results. (default: ./output)
#Optional arguments (used):
#  -j        Number of threads. (default: 5)
#Optional arguments (not used):
#  -p        The population of the sample [Asian, Black, Caucasian, Unknown, nonuse] for annotation. Unknown means use mean allele frequency in all populations. nonuse indicates only adopting mapping score
#              and considering zero-frequency alleles. (default: Unknown)
#  -d        Minimum score difference to assign a read to a gene. (default: 0.001)
#  -g        Whether use G group resolution annotation [0|1]. (default: 0)
#  -m        1 represents typing, 0 means only read assignment (default: 1)
#  -a        Prefix of filtered fastq file. (default: long_read)


# PART1: TYPING
if [[ "${WORKDIR}" != "" && ${BAM_IN} != "" ]]
then
    # STEP0: create work/output DIR
    HLA_WORK_DIR="${WORKDIR}/typing"
    if [[ ! -d ${HLA_WORK_DIR} ]]
    then
        mkdir "${HLA_WORK_DIR}"
    else
        echo "${HLA_WORK_DIR} exists"
    fi
    cd "${HLA_WORK_DIR}" || exit  # <== very important
    # STEP1: extract hla region from bam
    if [[ ! -f "${FQ_HLA}" ]];
    then
        conda activate hla
        samtools view -hb "${BAM_IN}" ${HLA_REGION} | \
            samtools fastq - | \
            bgzip -c > "${HLA_WORK_DIR}/${FQ_HLA}"
    else
        echo "[LOG] skipping hla extraction, ${FQ_HLA} exists"
    fi
    # STEP2: hla typing
    conda activate spechla_env
    # Run SpecHLA
    echo "CMD: python ${SPEC_HLA_BIN}  -r ${HLA_WORK_DIR}/${FQ_HLA}  -n  typing  -o  ${WORKDIR}  -j  ${THREADS}"
    python "${SPEC_HLA_BIN}" \
          -r  "${HLA_WORK_DIR}/${FQ_HLA}" \
          -n  typing \
          -o  "${WORKDIR}" \
          -j  ${THREADS}
else
  echo "[ERROR] Missing params ${WORKDIR} | ${BAM_IN}"
  exit 1
fi

# TODO: check this
#UTILS_DIR=$3   # <== path to utils
#ANNOT_DIR=$4   # <== annotation
# PART 2 personal reference
# STEP1: extract fasta files from typing results (with digits) and make personal reference
#cd ${HLA_WORK_DIR}
#conda activate hla
#PERSONAL_REF="personal_hla"
#HLA_REFERENCE_PERSONAL="${PERSONAL_REF}/my_hla.fasta"
#[[ ! -d "${PERSONAL_REF}" ]] && mkdir ${PERSONAL_REF} || echo "[LOG] directory '${PERSONAL_REF}' exists"
#if [[ ! -f "${HLA_REFERENCE_PERSONAL}" ]]; then
#    ${UTILS_DIR}/hlatyping_fasta.py hlafa \
#      --fasta "${ANNOT_DIR}/hla_gen.fasta.gz "\
#      --hla-names "${ANNOT_DIR}/hla_types_2digit_representative.csv" \
#      --hla-locus-masked $"{ANNOT_DIR}/chr6_hla_masked_grch38.fasta.gz" \
#      --type "${HLA_WORK_DIR}/results/${TODO} > ${HLA_REFERENCE_PERSONAL}"
#      # TODO: add file output
#fi
#samtools faidx ${HLA_REFERENCE_PERSONAL}
