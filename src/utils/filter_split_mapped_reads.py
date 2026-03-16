#!/usr/bin/env python3
import argparse
import pysam
from logger_debug import setup_log


logger = setup_log(__name__, True)
# min query aln 500
min_query_aln = 500
min_match_proportion = 0.80
# min match = 10%
max_edit_dist = 0.10
# min mapq
min_mapq_read = 10
# exclude flags -> read unmapped (0x4) + not primary alignment (0x100) + supplementary alignment (0x800)
read_flag_exclude = 4  # 2308


def read_analysis(sam_line):
    include_read = False
    if sam_line.mapping_quality > min_mapq_read and sam_line.flag & read_flag_exclude == 0:
        ([_M, _I, _D, _N, _S, _H, _P, _eq, _X, _B, _NM], _) = sam_line.get_cigar_stats()
        # proper_match = _M - _NM
        # sam_line_query_length = sam_line.query_length if sam_line.query_length != 0 else sam_line.infer_query_length()
        # proper_match_proportion = proper_match/sam_line_query_length
        edit_dist = _NM/sam_line.query_alignment_length
        add_star = "*" if sam_line.query_length == 0 else ""
        include_read = edit_dist <= max_edit_dist and sam_line.query_alignment_length >= min_query_aln
        logger.debug(f'{sam_line.query_name}: '
                     f'edit dist: {round(edit_dist, 3)}, '
                     f'query_alignment_length: {sam_line.query_alignment_length} '
                     f'{include_read}: {edit_dist <= max_edit_dist} & {sam_line.query_alignment_length >= min_query_aln}'
                     f'{add_star}')
    return include_read


def read_name_rm_coordinates(read_name):
    # read name is <name>_start_end, remove the last two in case "_" is part of the name
    split_name = read_name.split("_")[:-2]
    return "_".join(split_name)


# TODO: for reads spanning multiple genes, check distance
def filter_reads(user_args):
    def filter_stdin():
        _total_reads, _accepted_reads = 0, 0
        for ref in bam_file.references:
            for sam_line in bam_file.fetch(ref):
                read_included = read_analysis(sam_line)
                _total_reads += 1
                if read_included:
                    _accepted_reads += 1
                    bam_list.append(sam_line)
        return [_total_reads, _accepted_reads]

    def filter_files_out_list():
        for ref in bam_file.references:
            for sam_line in bam_file.fetch(ref):
                read_included = read_analysis(sam_line)
                if read_included:
                    raw_read_name = read_name_rm_coordinates(sam_line.query_name)
                    if raw_read_name not in bam_dict:
                        logger.debug("1..")
                        bam_dict[raw_read_name] = [ref]
                    else:
                        if ref not in bam_dict[raw_read_name]:
                            bam_dict[raw_read_name].append(ref)
                        else:
                            pass
        logger.debug(f'# ## FILTER ## #')
        allele_reads = {}
        reads_discarded = []
        for read_name in bam_dict.keys():
            if len(bam_dict[read_name]) == 1:
                [ref] = bam_dict[read_name]
                if ref != "chr6_hla_locus_26mb-36mb":
                    if ref not in allele_reads:
                        allele_reads[ref] = [read_name]
                    else:
                        allele_reads[ref].append(read_name)
            elif len(bam_dict[read_name]) == 2:
                [ref1, ref2] = bam_dict[read_name]
                if ref1 == "chr6_hla_locus_26mb-36mb" or ref2 == "chr6_hla_locus_26mb-36mb":
                    if ref1 == "chr6_hla_locus_26mb-36mb":
                        ref = ref2
                    else:
                        ref = ref1
                    if ref not in allele_reads:
                        allele_reads[ref] = [read_name]
                    else:
                        allele_reads[ref].append(read_name)
                else:
                    reads_discarded.append(read_name)
            else:
                reads_discarded.append(read_name)
                # more than 2 or none are discarded
        # we expect only 2 ref
        if len(allele_reads.keys()) == 2:
            [reads_allele1, reads_allele2] = [allele_reads[k] for k in allele_reads.keys()]
            # allele 1
            file_handler = open(f'typing-mapped/{user_args.out_list}_allele1.txt', "wt")
            for each_read in reads_allele1:
                file_handler.write(f'{each_read}\n')
            file_handler.close()
            # allele 2
            file_handler = open(f'typing-mapped/{user_args.out_list}_allele2.txt', "wt")
            for each_read in reads_allele2:
                file_handler.write(f'{each_read}\n')
            file_handler.close()
            # reads discarded
            file_handler = open(f'typing-mapped/{user_args.out_list}_discarded.txt', "wt")
            for each_read in reads_discarded:
                file_handler.write(f'{each_read}\n')
            file_handler.close()
        else:
            logger.warning(f'Error during the analysis: more than 2 alleles found per gene')
            logger.debug(f'Alleles: {allele_reads}')

    def filter_files_out_bam():
        _total_reads, _accepted_reads, _uniq_chr_reads = 0, 0, 0
        for ref in bam_file.references:
            for sam_line in bam_file.fetch(ref):
                read_included = read_analysis(sam_line)
                _total_reads += 1
                if read_included:
                    logger.debug(f'# ## READ INCLUDE ## #')
                    _accepted_reads += 1
                    raw_read_name = read_name_rm_coordinates(sam_line.query_name)
                    if raw_read_name not in bam_dict:
                        logger.debug("1..")
                        bam_dict[raw_read_name] = {ref: [sam_line]}
                    else:
                        logger.debug(f'2..  {raw_read_name}')
                        if ref not in bam_dict[raw_read_name]:
                            bam_dict[raw_read_name][ref] = [sam_line]
                            logger.debug(f'2.1  {bam_dict[raw_read_name].keys()}')
                        else:
                            bam_dict[raw_read_name][ref].append(sam_line)
                            logger.debug(f'2.2  {bam_dict[raw_read_name].keys()}')
        logger.debug(f'# ## FILTER ## #')
        tmp_2_alleles = {}
        for read_name in bam_dict.keys():
            read_chrom_overlap = bam_dict[read_name].keys()
            # read only covers 1 HLA "chromosome"
            if len(read_chrom_overlap) == 1:
                logger.debug(f'YES  =>> {read_name} : {read_chrom_overlap}')
                log_reads.write(f'YES  =>> {read_name} : {read_chrom_overlap}\n')
                _uniq_chr_reads = write_filtered_reads(read_chrom_overlap, read_name)
            # read only covers 2 HLA "chromosomes"
            else:
                # 2nd chromosome is decoy, good
                if len(read_chrom_overlap) == 2 and "chr6_hla_locus_26mb-36mb" in read_chrom_overlap:
                    logger.debug(f'YES  =>> {read_name} : {read_chrom_overlap}')
                    log_reads.write(f'YES  =>> {read_name} : {read_chrom_overlap}\n')
                    _uniq_chr_reads = write_filtered_reads(read_chrom_overlap, read_name)
                elif (len(read_chrom_overlap) == 2 and "chr6_hla_locus_26mb-36mb" not in read_chrom_overlap) or \
                        (len(read_chrom_overlap) == 3 and "chr6_hla_locus_26mb-36mb" in read_chrom_overlap):
                    [chr1, chr2] = [x for x in read_chrom_overlap if x != "chr6_hla_locus_26mb-36mb"]
                    # example of HLA format: 'HLA00037__A_03__allele1'
                    # example of split indices '(0:HLA00037)_(1:)_(2:A)_(3:03)_(4:)_(5:allele1)'
                    hla_index = 2
                    hla1, hla2 = chr1.split("_")[hla_index], chr2.split("_")[hla_index]
                    # hla chromosomes are not the same (not 2 alleles), good
                    if hla1 != hla2:
                        tmp2hla = f'{chr1}|{chr2}'
                        if tmp2hla not in tmp_2_alleles:
                            tmp_2_alleles[tmp2hla] = 1
                        else:
                            tmp_2_alleles[tmp2hla] += 1
                        logger.debug(f'YES  =>> {read_name} : {read_chrom_overlap}')
                        log_reads.write(f'YES  =>> {read_name} : {read_chrom_overlap}\n')
                        _uniq_chr_reads = write_filtered_reads(read_chrom_overlap, read_name)
                    else:
                        logger.debug(f'NOPE =>> {read_name} : {len(read_chrom_overlap)}|{read_chrom_overlap}')
                        log_reads.write(f'NOPE =>> {read_name} : {len(read_chrom_overlap)}|{read_chrom_overlap}\n')
                else:
                    logger.debug(f'NOPE =>> {read_name} : {len(read_chrom_overlap)}|{read_chrom_overlap}')
                    log_reads.write(f'NOPE =>> {read_name} : {len(read_chrom_overlap)}|{read_chrom_overlap}\n')
        logger.debug(tmp_2_alleles)
        return [_total_reads, _accepted_reads, _uniq_chr_reads]

    def write_filtered_reads(_read_chrom_overlap, _read_name):
        _uniq_chr_reads = 0
        for chrom in _read_chrom_overlap:
            for this_sam_line in bam_dict[_read_name][chrom]:
                _uniq_chr_reads += 1
                bam_list.append(this_sam_line)
        return _uniq_chr_reads

    # open and close files
    out_is_bam = False
    out_is_list = False
    if user_args.out_bam == "" and user_args.out_list == "":
        logger.warning(f'Need to provide at least one output, making default bam stdout')
        user_args.out_bam = "-"
        out_is_bam = True
    elif user_args.out_bam == "" and user_args.out_list != "":
        out_is_list = True
    elif user_args.out_bam != "" and user_args.out_list != "":
        logger.warning(f'Only one output is needed, making default user selected bam')
        out_is_bam = True
    else:
        logger.warning(f'Something went wrong, making default output bam stdout')
        user_args.out_bam = "-"
        out_is_bam = True
    bam_file = pysam.AlignmentFile(user_args.in_bam, "rb")
    log_reads = open(user_args.log_reads, "wt")
    input_is_stdin = user_args.in_bam == "-"
    bam_list = []
    bam_dict = {}
    uniq_chr_reads = None
    total_reads, accepted_reads, uniq_chr_reads = 0, 0, 0
    if input_is_stdin:
        logger.info("little output is given")
        [total_reads, accepted_reads] = filter_stdin()
    else:
        if out_is_bam:
            out_bam = pysam.AlignmentFile(user_args.out_bam, "wb", template=bam_file)
            [total_reads, accepted_reads, uniq_chr_reads] = filter_files_out_bam()
            for line in bam_list:
                out_bam.write(line)
            out_bam.close()
        elif out_is_list:
            filter_files_out_list()
    if input_is_stdin or out_is_bam:
        logger.info(f'Total reads: {total_reads}')
        logger.info(f'Reads pass filter 1: {accepted_reads}')
        logger.info(f'Reads pass filter 2: {uniq_chr_reads}')
    bam_file.close()


def get_arguments():
    filter_split_reads_help = """
    Filter_split_reads:
            --in-bam/-i          File input BAM with ONT reads (if - can be stdin)
            --out-bam/-o         Output BAM file (if - can be stdout)
            --out-list/-l        Prefix for output read name in list bay allele, will create allele1 and allele2 files
            --reads-used-log/-r  Reads used in the personal genome
            --version/-v         Shows current version/build
    """
    parser = argparse.ArgumentParser(description="", usage=filter_split_reads_help)

    # ############################################################################################ #
    parser.add_argument("-v", "--version", action="version", version=f"0.1")
    parser.add_argument("-i", "--in-bam", type=str, required=True, dest='in_bam', default="", help='Input bam file')
    parser.add_argument("-o", "--out-bam", type=str, required=False, dest='out_bam', default="",
                        help='Output filtered bam file')
    parser.add_argument("-l", "--out-list", type=str, required=False, dest='out_list', default="",
                        help='Prefix for output read name lists per allele (1 and 2)')
    parser.add_argument("-r", "--reads-used-log", type=str, required=True, dest='log_reads', default="",
                        help='Output for reads used log')
    return parser.parse_args()


def main():
    args = get_arguments()
    filter_reads(args)


# main
if __name__ == '__main__':
    main()

