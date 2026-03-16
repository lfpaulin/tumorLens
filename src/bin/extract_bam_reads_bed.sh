# UPDATE THE NAME OF THE FILE ...
BED_USE=$1
BAM_IN=$2    # Here goes the bam INPUT
BAM_OUT=$3   # Here goes the bam OUTPUT

# add header
samtools view -H "${BAM_IN}" > tmp.sam;
# add every single candidate from the bed
cat "${BED_USE}" | perl -lane 'print "$F[0]:$F[1]-$F[2]"' | \
    xargs -n1 -t -I{} samtools view "${BAM_IN}" {} | sort | uniq >> tmp.sam;

# make it bam and sort it
samtools view -bSh tmp.sam | samtools sort -o "${BAM_OUT}" -
samtools index "${BAM_OUT}"
