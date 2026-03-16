import sys

chunk_size = int(1e4)  # 10kb
counter = 1

outbed = open(f'methyl_chunk.bed', "w")
for line in sys.stdin:
    [chrom, chr_len, _, _, _] = line.rstrip("\n").split("\t")
    chr_len = int(chr_len)
    if chr_len > 1e6 and "X" not in chrom and "Y" not in chrom:
        for start in range(1,chr_len,chunk_size):
            counter_string = f'{counter:03d}'
            end = min(start+chunk_size-1, chr_len)
            #outbed = open(f'chunk_{counter_string}_{chrom}_{start}_{end}.bed', "w")
            outbed.write(f'{chrom}\t{start}\t{end}\n')
            #outbed.close()
            counter += 1

outbed.close()

