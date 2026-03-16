#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --mem=16G
#
#
#
#OUT_ERR_CHDIR

# load conda env
conda activate hla

## arguments
REFERENCE=$1
BAM_MOD_PHASED_TUMOR=$2  # we need to split it tumor
BAM_MOD_PHASED_CONTROL=$3  # we need to split it tumor

OUTPUT_DIR="methyl_phase"
BAM_MOD_HPx="align_mod_genes_hp"
THREADS_USE=8

if [[ ! -d  "${OUTPUT_DIR}" ]];
then
    mkdir "${OUTPUT_DIR}"
fi
cd  "${OUTPUT_DIR}" || exit

# expected haplotypes are 1 and 2
for hp in "1" "2";
do
    # TUMOR
    bam_out="${BAM_MOD_HPx}_tumor_${hp}.bam"
    echo "${BAM_MOD_PHASED_TUMOR} => ${bam_out}";
    samtools view -H ${BAM_MOD_PHASED_TUMOR} > tmp_bam_header.txt;
    samtools view ${BAM_MOD_PHASED_TUMOR} | grep "HP:i:${hp}" > tmp_bam_reads.txt;
    cat tmp_bam_header.txt tmp_bam_reads.txt | samtools view -hb > ${bam_out}
    samtools index ${bam_out}
    #modkit
    methhp="${BAM_MOD_HPx}_tumor_${hp}.bed.gz"
    modkit pileup \
      --cpg \
      --ref ${REFERENCE} \
      --threads ${THREADS_USE} \
      --log-filepath "log_${BAM_MOD_HPx}_tumor_${hp}.txt" \
      ${bam_out} - | bgzip -c > ${methhp}
    tabix --preset bed ${methhp}

    # CONTROL
    bam_out="${BAM_MOD_HPx}_control_${hp}.bam"
    echo "${BAM_MOD_PHASED_CONTROL} => ${bam_out}";
    samtools view -H ${BAM_MOD_PHASED_CONTROL} > tmp_bam_header.txt;
    samtools view ${BAM_MOD_PHASED_CONTROL} | grep "HP:i:${hp}" > tmp_bam_reads.txt;
    cat tmp_bam_header.txt tmp_bam_reads.txt | samtools view -hb > ${bam_out}
    samtools index ${bam_out}
    #modkit
    methhp="${BAM_MOD_HPx}_control_${hp}.bed.gz"
    modkit pileup \
      --cpg \
      --ref ${REFERENCE} \
      --threads ${THREADS_USE} \
      --log-filepath "log_${BAM_MOD_HPx}_control_${hp}.txt" \
      ${bam_out} - | bgzip -c > ${methhp}
    tabix --preset bed ${methhp}
done

# Diff Methyl Analysis
GENES=$4
#SAMPLE_NAMES=$5  # comma sep, in the following order: tumor, control
TUMOR_NAME="tumor" #"$(echo   "${SAMPLE_NAMES}" | cut -d "," -f 1)"
CONTROL_NAME="control" #"$(echo "${SAMPLE_NAMES}" | cut -d "," -f 2)"

if [[ ! -f "${BAM_MOD_HPx}_control_1.bed.gz.tbi" ]];
then
    tabix --preset bed "${BAM_MOD_HPx}_control_1.bed.gz"
fi
if [[ ! -f "${BAM_MOD_HPx}_control_2.bed.gz.tbi" ]];
then
    tabix --preset bed "${BAM_MOD_HPx}_control_2.bed.gz"
fi
if [[ ! -f "${BAM_MOD_HPx}_tumor_1.bed.gz.tbi" ]];
then
    tabix --preset bed "${BAM_MOD_HPx}_tumor_1.bed.gz"
fi
if [[ ! -f "${BAM_MOD_HPx}_tumor_2.bed.gz.tbi" ]];
then
    tabix --preset bed "${BAM_MOD_HPx}_tumor_2.bed.gz"
fi
modkit dmr multi \
    --sample "${BAM_MOD_HPx}_control_1.bed.gz"  "${CONTROL_NAME}_hap1" \
    --sample "${BAM_MOD_HPx}_control_2.bed.gz"  "${CONTROL_NAME}_hap2" \
    --sample "${BAM_MOD_HPx}_tumor_1.bed.gz"    "${TUMOR_NAME}_hap1" \
    --sample "${BAM_MOD_HPx}_tumor_2.bed.gz"    "${TUMOR_NAME}_hap2" \
    --out-dir dmr_by_hap \
    --regions-bed  "${GENES}" \
    --ref  "${REFERENCE}" \
    --base C --threads ${THREADS_USE} \
    --log-filepath _log_dmr_multi.txt

if [[ -f "results_stats.txt" ]];
then
 rm results_stats.txt
fi
if [[ -f "results_genes.txt" ]];
then
 rm results_genes.txt
fi

for bed_res in dmr_results_by_hap/*bed;
do
 echo "${bed_res}" >> results_stats.txt
 echo "${bed_res}" >> results_genes.txt
 cut -f 1,2,4,5,12,13 "${bed_res}" | awk '{if ($5-$6 > 0.2 || $5-$$6<-0.2) print $_"\t"$5-$6}' | wc -l >> results_stats.txt
 cut -f 1,2,4,5,12,13 "${bed_res}" | awk '{if ($5-$6 > 0.2 || $5-$$6<-0.2) print $_"\t"$5-$6}'         >> results_genes.txt
done