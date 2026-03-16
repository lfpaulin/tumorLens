#!/usr/bin/env python3
import argparse
import pysam
from logger_debug import setup_log


logger = setup_log(__name__, True)


# TODO: from here below, update
def read_analysis(sam_line, out_bam, max_edit_dist, min_match_proportion):
    # min_match_proportion
    _accepted_reads = 0
    _total_reads = 0
    ([_M, _I, _D, _N, _S, _H, _P, _eq, _X, _B, _NM], _) = sam_line.get_cigar_stats()
    proper_match = _M - _NM
    edit_dist = _NM/sam_line.query_alignment_length
    proper_match_proportion = proper_match/sam_line.query_length
    if edit_dist <= max_edit_dist and proper_match_proportion >= min_match_proportion:
        # logger.warning(f'{sam_line.query_name}: '
        #             f'edit dist: {round(edit_dist, 3)}|{1-round(edit_dist, 3)}, '
        #             f'query_alignment_length: {sam_line.query_alignment_length}')
        out_bam.write(sam_line)
        _accepted_reads += 1
    _total_reads += 1
    return [_total_reads, _accepted_reads]


def filter_reads_stdin():
    bam_file = pysam.AlignmentFile("-", "rb")
    out_bam = pysam.AlignmentFile("-",  "wb", template=bam_file)
    # min match = 80%
    min_match_proportion = 0.80
    # min match = 5%
    max_edit_dist = 0.05
    total_reads = 0
    accepted_reads = 0
    for sam_line in bam_file:
        [_total_reads, _accepted_reads] = read_analysis(sam_line, out_bam, max_edit_dist, min_match_proportion)
        total_reads += _total_reads
        accepted_reads += _accepted_reads
    logger.warning(f'Total reads: {total_reads}')
    logger.warning(f'Reads pass: {accepted_reads}')
    bam_file.close()
    out_bam.close()


# TODO: for reads spanning multiple genes, check distance
def filter_reads_files(user_args):
    bam_file = pysam.AlignmentFile(user_args.in_bam, "rb")
    out_bam = pysam.AlignmentFile(user_args.out_bam, "wb", template=bam_file)
    # min match = 80%
    min_match_proportion = 0.80
    total_reads = 0
    accepted_reads = 0
    for ref in bam_file.references:
        for sam_line in bam_file.fetch(ref):
            [_total_reads, _accepted_reads] = read_analysis(sam_line, out_bam, min_match_proportion)
            total_reads += _total_reads
            accepted_reads += _accepted_reads
    logger.info(f'Total reads: {total_reads}')
    logger.info(f'Reads pass: {accepted_reads}')
    bam_file.close()
    out_bam.close()


def get_arguments():
    filter_split_reads_help = """
    Filter_split_reads:
        from STDIN: stdin
        from FILES: files
            --in-bam/-i     File input with ONT reads
            --out-bam/-o    Directory with the modbam files (unmapped bam)
            --version/-v    Shows current version/build
    """
    parser = argparse.ArgumentParser(description="", usage=filter_split_reads_help)
    subparsers = parser.add_subparsers(help=filter_split_reads_help, dest="command")

    # ############################################################################################ #
    # stdin
    stdin_help = "Input and output are STDIN and STDOUT respectively"
    subparser_stdin = subparsers.add_parser("stdin", help=stdin_help)
    # files
    files_help = "..."
    subparser_files = subparsers.add_parser("files", help=files_help)
    subparser_files.add_argument("-v", "--version", action="version", version=f"0.1")
    subparser_files.add_argument("-i", "--in-bam", type=str, required=True, dest='in_bam', default="",
                                 help='Input bam file')
    subparser_files.add_argument("-o", "--out-bam", type=str, required=True, dest='out_bam', default="",
                                 help='Output filtered bam file')
    return parser.parse_args()


def main():
    args = get_arguments()
    if args.command == "stdin":
        filter_reads_stdin()
    else:
        filter_reads_files(args)


# main
if __name__ == '__main__':
    main()

