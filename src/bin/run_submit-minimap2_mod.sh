#!/bin/bash
#SBATCH --ntasks=12
#SBATCH --mem=32G
#
#
#
#OUT_ERR_CHDIR

# load conda env
conda activate hla

# increase number of opened files
ulimit -n 32768

## arguments
PRESET="ont"
MINIMAP_PRESET="map-ont"
REFERENCE=$1
READS_DIR=$2
FINAL_DIR=$3
TMP_DIR="/space1/tmp/slurm.$SLURM_JOB_ID"
#MODBAM_IS_FILE=$4

minimap_version=$(minimap2 --version)
samtools_version=$(samtools --version | grep samtools | head -n 1 | cut -d " " -f 2)
echo "[INFO] mapping ${PRESET} reads:"
echo "       minimap2 version ${minimap_version}"
echo "       samtools version ${samtools_version}"

# threads used by program, total is 12
THREADS_SAM_CAT=1  #  2
THREADS_MINIMAP=8  #  8
THREADS_SAM_SRT=1  #  2

# work files
LIST_MODBAM="modb_bams.list"
OUT="align_mod"

# get all the files to be mapped in a list
ls ${READS_DIR}/*bam > ${TMP_DIR}/${LIST_MODBAM}

# cat | fastq | mapping | sort
if [[ ! -f "${FINAL_DIR}/${OUT}.bam" ]];
then
    cd "${TMP_DIR}" || exit
    samtools cat -b ${TMP_DIR}/${LIST_MODBAM} -@ ${THREADS_SAM_CAT} | samtools fastq -T Mm,Ml,MM,ML - | \
        minimap2 -ax ${MINIMAP_PRESET} -t ${THREADS_MINIMAP} -y ${REFERENCE} - | \
        samtools sort -m 2G -@ ${THREADS_SAM_SRT} - > ${TMP_DIR}/${OUT}.bam

    echo "[INFO]: CMD: samtools cat -b ${LIST_MODBAM} -@ ${THREADS_SAM_CAT} | samtools fastq -T Mm,Ml,MM,ML - | \ "
    echo "                 minimap2 -ax ${MINIMAP_PRESET} -t ${THREADS_MINIMAP} ${REFERENCE} - | \ "
    echo "                 samtools sort -m 2G -@ ${THREADS_SAM_SRT} - > ${TMP_DIR}/${OUT}.bam "
    # index
    echo "[INFO] indexing ${OUT} samtools"
    samtools index -@ ${THREADS_MINIMAP}  ${TMP_DIR}/${OUT}.bam
    # BAM
    mv ${TMP_DIR}/${OUT}.bam     ${FINAL_DIR}/${OUT}.bam
    mv ${TMP_DIR}/${OUT}.bam.bai ${FINAL_DIR}/${OUT}.bam.bai
fi
