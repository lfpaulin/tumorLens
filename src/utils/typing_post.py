#!/usr/bin/env python3
import logging
import os
import gzip
import sys

import pysam
import argparse
import subprocess
from .logger_debug import setup_log
from .term_colors import TermColors


logger = setup_log(__name__, False)
tumor = "tumor"
control = "control"


class TypeArgs:
    sample: str = ""
    work_dir: str = ""
    min_query_aln: int = 0
    max_edit_dist: float = 0.0
    use_edit_dist: bool = False
    annot_data: str = ""
    data_dir: str = ""

    def __init__(self, sample: str, work_dir: str, annot_data: str, min_query_aln: int, max_edit_dist: float,
                 use_edit_dist: bool, data_dir: str):
        self.sample = sample
        self.work_dir = work_dir
        self.min_query_aln = min_query_aln
        self.max_edit_dist = max_edit_dist
        self.use_edit_dist = use_edit_dist
        self.annot_data = annot_data
        self.data_dir = data_dir


class HLAAlignmentPreprocess:
    name: str = ""
    fastq: str = ""
    dir_work: str = ""
    hla_fasta4 = "hla_gen.format.filter.extend.DRB.no26789.v2.fasta.gz"
    hla_fasta = "hla_gen.format.filter.extend.DRB.no26789.fasta.gz"
    dir_type_post: str = "post_analysis"
    file_type: str = "hla.result.details.txt"
    gene_names = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
    # specHLA type: 'A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'
    typed_result = {
        'A': ["", ""],
        'B': ["", ""],
        'C': ["", ""],
        'DPA1': ["", ""],
        'DPB1': ["", ""],
        'DQA1': ["", ""],
        'DQB1': ["", ""],
        'DRB1': ["", ""]
    }
    typing_fastq_lr = {
        'A': "A.long_read.fq.gz",
        'B': "B.long_read.fq.gz",
        'C': "C.long_read.fq.gz",
        'DPA1': "DPA1.long_read.fq.gz",
        'DPB1': "DPB1.long_read.fq.gz",
        'DQA1': "DQA1.long_read.fq.gz",
        'DQB1': "DQB1.long_read.fq.gz",
        'DRB1': "DRB1.long_read.fq.gz"
    }
    typing_fastq_lr_split = {
        'A': "A_long_read_split.fastq.gz",
        'B': "B_long_read_split.fastq.gz",
        'C': "C_long_read_split.fastq.gz",
        'DPA1': "DPA1_long_read_split.fastq.gz",
        'DPB1': "DPB1_long_read_split.fastq.gz",
        'DQA1': "DQA1_long_read_split.fastq.gz",
        'DQB1': "DQB1_long_read_split.fastq.gz",
        'DRB1': "DRB1_long_read_split.fastq.gz"
    }

    def __init__(self, sample_name: str, work_dir: str, digits: int = 4):
        self.name = sample_name
        self.dir_work = work_dir
        self.digits = digits
        self.digits_index = int(digits/2)
        # dirs
        self.dir_type_post = f'{self.dir_work}/{self.dir_type_post}'

    def __repr__(self):
        return (f'name: {self.name}\nfastq: {self.fastq}\nhla_fasta: {self.hla_fasta}\ndir_work: {self.dir_work}'
                f'\n{self.typed_result}')

    def get_type(self, hla_gene):
        return self.typed_result[hla_gene]

    def get_long_reads(self, hla_gene):
        return self.typing_fastq_lr[hla_gene]

    def get_long_reads_split(self, hla_gene):
        return self.typing_fastq_lr_split[hla_gene]


class FastqSplit:
    def __init__(self, new_read_len=3000):
        # expects 3 tab-sep cols: read_name, seq, qual
        # this is done by samtools view FASTQ-file | cut -f 1,10,11
        self.new_read_len = new_read_len
        self.read_name, self.seq, self.seq_quality, self.read_len = "", "", [], 0
        self.output = []

    def process_read(self, fastq_record):
        self.set_read_info(fastq_record)
        self.split_read()

    def set_read_info(self, fastq_record):
        self.read_name, self.seq, self.seq_quality = fastq_record.name, fastq_record.sequence, fastq_record.quality
        self.read_len = len(self.seq)

    def split_read(self):
        for start_pos in range(0, self.read_len, self.new_read_len):
            end_pos = start_pos+self.new_read_len if start_pos+self.new_read_len < self.read_len else self.read_len
            self.make_fastq(start_pos, end_pos)

    def make_fastq(self, read_start, read_end):
        new_read = self.seq[read_start:read_end]
        new_qual = self.seq_quality[read_start:read_end]
        if len(new_read) != 0:
            self.output.append(f'@{self.read_name}_{read_start+1}_{read_end}\n{new_read}\n+\n{new_qual}')


def get_types(sample_name: str, work_dir: str):
    hla_types = HLAAlignmentPreprocess(sample_name, work_dir)
    hla_typing_read = open(f'{work_dir}/{hla_types.file_type}', "r")
    hla_genes = []
    for line in hla_typing_read:
        line = line.rstrip("\n")
        if line.startswith("H"):
            gene, hla_type = line.split("\t")[:2]
            hla_type = use_type_digits(hla_type, hla_types.digits_index)
            [_, gene, hla_allele] = gene.split("_")
            hla_allele = int(hla_allele) - 1  # make it index 1=>0 | 2=>1
            hla_types.typed_result[gene][hla_allele] = hla_type
            if gene not in hla_genes:
                hla_genes.append(gene)
    # check if empty
    for gene in hla_genes:
        if "" == hla_types.typed_result[gene][0] and "" != hla_types.typed_result[gene][1]:
            hla_types.typed_result[gene][0] = hla_types.typed_result[gene][1]
        if "" == hla_types.typed_result[gene][1] and "" != hla_types.typed_result[gene][0]:
            hla_types.typed_result[gene][1] = hla_types.typed_result[gene][0]
    return hla_types


def use_type_digits(hla_type, digits_index):
    return ":".join(hla_type.split(":")[:digits_index])


def not_same_type(hla_types: list[str], digits: int = 4):
    digits_index = int(digits/2)
    hla1, hla2 = hla_types
    if digits_index == 4:
        return hla1 != hla2
    else:
        hla1 = ":".join(hla1.split(":")[:digits_index])
        hla2 = ":".join(hla2.split(":")[:digits_index])
        return hla1 != hla2


def extract_needed_hla(hla_fasta_handler: pysam.FastaFile, dir_work_post: str, hla_gene_name: str, hla_types: list[str],
                       hla_fasta_handler_v2: pysam.FastaFile):
    seq, pretty_type = "", ""
    """
    TODO: update the analysis:
        1. make a single fasta reference with both alleles (and include mask now)
        2. improve hla type fetch with => v1 first, v2 second, 2 digits last
    """
    hla_type_fasta_out = open(f'{dir_work_post}/{hla_gene_name}.fasta', "w")
    for hla_type in hla_types:
        try:
            logger.debug(f'{hla_type} Using HLA-db v1')
            seq = hla_fasta_handler.fetch(hla_type)
        except KeyError:
            try:
                logger.warning(f'{hla_type} not found in HLA-db v1.. Using HLA-db v2')
                seq = hla_fasta_handler_v2.fetch(hla_type)
            except KeyError:
                pass
            test_type = [ref for ref in hla_fasta_handler.references + hla_fasta_handler_v2.references if hla_type in ref]
            if len(test_type) > 0:
                # use first entry
                use_hla = test_type[0]
                if use_hla in hla_fasta_handler.references:
                    logger.warning(f'{hla_type} not exact match, using {use_hla}, from v1')
                    seq = hla_fasta_handler.fetch(use_hla)
                elif use_hla in hla_fasta_handler_v2.references:
                    logger.warning(f'{hla_type} not exact match, using {use_hla}, from v2')
                    seq = hla_fasta_handler_v2.fetch(use_hla)
                else:
                    logger.error(f'{hla_type} not exact match, using {use_hla} from ?.. exiting with error')
                    sys.exit(1)
            else:
                logger.warning(f'No sequence found for {hla_type}')

        pretty_type = "_".join("_".join(hla_type.split(":")).split("*"))
        logger.debug(f'{hla_type} => {pretty_type}')
        hla_type_fasta_out.write(f'>{hla_type}\n{seq}\n')
    hla_type_fasta_out.close()
    logger.debug(f'Extraction success')
    return f'{dir_work_post}/{hla_gene_name}.fasta'


def align_hla(hla_name: str, dir_work_post: str, params: list[str]):
    logger.debug(f'Starting alignment process for {hla_name}')
    # aln.sh REFERENCE=$1 READS=$2 OUT=$3
    [hla_gene_fasta, reads, mapped_reads] = params
    # TODO: improve parameters based on spechHLA
    if not os.path.exists(f'{dir_work_post}/{mapped_reads}'):
        logger.debug(f'Alignment...')
        cmd_aln = (f'minimap2 -ax map-ont -t 4 {hla_gene_fasta} {reads} | '
                   f'samtools view -hb | samtools sort > {dir_work_post}/{mapped_reads}; '
                   f'samtools index {dir_work_post}/{mapped_reads}')
        logger.debug(f'CMD: {cmd_aln}')
        run_cmd = subprocess.run(cmd_aln, shell=True, capture_output=True, text=True)
        if run_cmd.stderr != "":
            if "rror" in run_cmd.stderr:
                logger.error(run_cmd.stderr)
        logger.debug(f'Alignment success')
    else:
        logger.debug(f'No alignment...')
    return f"{dir_work_post}/{mapped_reads}"


def read_analysis(sam_file: str, min_query_aln: int, max_edit_dist: float, hla_alleles: dict,
                  use_edit_dist: bool) -> tuple[int, dict]:
    # min_query_aln = 600
    # min_match_proportion = 0.80
    # max_edit_dist = 0.05
    # min mapq
    min_mapq_read = 10
    # exclude flags -> read unmapped (0x4) + not primary alignment (0x100) + supplementary alignment (0x800)
    read_flag_exclude = 2308
    # output
    sam_file_handler = pysam.AlignmentFile(sam_file, "r")
    for sam_line in sam_file_handler:
        if sam_line.flag & read_flag_exclude == 0:
            if use_edit_dist:
                include_read = False
                if sam_line.mapping_quality > min_mapq_read:
                    ([_M, _I, _D, _N, _S, _H, _P, _eq, _X, _B, _NM], _) = sam_line.get_cigar_stats()
                    # proper_match = _M - _NM
                    # sam_line_query_length = sam_line.query_length if sam_line.query_length != 0 else (
                    #     sam_line.infer_query_length())
                    # proper_match_proportion = proper_match/sam_line_query_length
                    edit_dist = _NM/sam_line.query_alignment_length
                    # add_star = "*" if sam_line.query_length == 0 else ""
                    include_read = edit_dist <= max_edit_dist and sam_line.query_alignment_length >= min_query_aln
                    # logger.debug(f'{sam_line.query_name}: '
                    #              f'edit dist: {round(edit_dist, 3)}, '
                    #              f'query_alignment_length: {sam_line.query_alignment_length} '
                    #               (f'{include_read}: {edit_dist <= max_edit_dist} & '
                    #                f'{sam_line.query_alignment_length >= min_query_aln}')
                    #              f'{add_star}')
            else:
                include_read = True
            if include_read:
                hla_alleles[sam_line.reference_name].append(sam_line.query_name)
    # Remove cross mapping
    allele1, allele2 = hla_alleles.keys()
    reads_allele1 = set([read.split("_")[0] for read in hla_alleles[allele1]])
    reads_allele2 = set([read.split("_")[0] for read in hla_alleles[allele2]])
    reads_allele1_keep = [read for read in reads_allele1 if read not in reads_allele2]
    reads_allele2_keep = [read for read in reads_allele2 if read not in reads_allele1]
    logger.debug(f'{len(hla_alleles[allele1])}|{len(reads_allele1)}|{len(reads_allele1_keep)} : '
                 f'{len(hla_alleles[allele2])}|{len(reads_allele2)}|{len(reads_allele2_keep)}')
    include_reads_count = len(reads_allele1_keep) + len(reads_allele2_keep)
    return include_reads_count, {allele1: reads_allele1_keep, allele2: reads_allele2_keep}


def methylation(tumor_dir: str, tumor_hla_reads: dict, control_dir: str, control_hla_reads: dict,
                ref_fasta: str, results_dir: str, annot_bed: pysam.TabixFile, annot_type: str) -> dict:

    def extract_name_form_path(file_path: str):
        return "_".join(os.path.basename(file_path).split("_")[2:4])

    def hla_pretty_name(hla_type_format: str):
        return "_".join("".join(hla_type_format.split(":")).split("*"))

    aln_mod = pysam.TabixFile(f'{tumor_dir}/../methylation/align_mod_modkit.bed.gz')
    if len(aln_mod.contigs) == 0:
        return {}
    extra_padding = 10000
    results_allele_meth = {}
    tumor_dir_methyl = f'{results_dir}/typing/post_methyl_tumor'
    control_dir_methyl = f'{results_dir}/typing/post_methyl_control'
    hla_type_dmr_dir = f'{results_dir}/typing/dmr'
    logger.setLevel(logging.DEBUG)
    os.mkdir(f'{results_dir}/typing') if not os.path.exists(f'{results_dir}/typing') else None
    os.mkdir(tumor_dir_methyl) if not os.path.exists(tumor_dir_methyl) else None
    os.mkdir(control_dir_methyl) if not os.path.exists(control_dir_methyl) else None
    os.mkdir(hla_type_dmr_dir) if not os.path.exists(hla_type_dmr_dir) else None
    tumor_hla_bam = pysam.AlignmentFile(f'{tumor_dir}/_hla_region.bam')
    control_hla_bam = pysam.AlignmentFile(f'{control_dir}/_hla_region.bam')
    hla_bed, reg_bam, reg_bed, reg_fetch = "", "", "", ""
    logger.info(f'Annotation element: "{annot_type}"')
    for hla_gene in tumor_hla_reads.keys():
        logger.info(f'{TermColors.bold} {hla_gene} {TermColors.end}')
        if hla_gene not in results_allele_meth:
            results_allele_meth[hla_gene] = {tumor: None, control: None, "bed": ""}
        annot_bed.seek(0)
        for bed_reg in annot_bed.fetch():
            c, s, e, g, _, _ = bed_reg.split("\t")
            sp = int(s) - extra_padding
            ep = int(e) + extra_padding
            if f'HLA-{hla_gene}' == g:
                reg_bam = f'{c}:{sp}-{ep}'
                logger.info(f'HLA-{hla_gene} => {reg_bam} ')
                reg_bed = f'{c}\t{s}\t{e}'
                reg_fetch = f'{c}:{s}-{e}'
                hla_bed = open(f'{tumor_dir_methyl}/typing_meth_hla_{hla_gene}_{annot_type}.bed', "w")
                hla_bed.write(f'{reg_bed}\n')
                hla_bed.close()
                hla_bed = f'{tumor_dir_methyl}/typing_meth_hla_{hla_gene}_{annot_type}.bed'
                results_allele_meth[hla_gene]["bed"] = hla_bed
                break
        if 2 == len(tumor_hla_reads[hla_gene]):
            [hla1, hla2] = tumor_hla_reads[hla_gene].keys()
            reads_a1 = tumor_hla_reads[hla_gene][hla1]
            reads_a2 = tumor_hla_reads[hla_gene][hla2]
            hla1_pretty = hla_pretty_name(hla1)
            hla2_pretty = hla_pretty_name(hla2)
            logger.info(f'{hla1_pretty} / {hla2_pretty}')
            tumor_allele1_bam = f'{tumor_dir_methyl}/typing_meth_{hla1_pretty}_{annot_type}.bam'
            tumor_allele1_bam_out = pysam.AlignmentFile(tumor_allele1_bam, "wb", template=tumor_hla_bam)
            tumor_allele2_bam = f'{tumor_dir_methyl}/typing_meth_{hla2_pretty}_{annot_type}.bam'
            tumor_allele2_bam_out = pysam.AlignmentFile(tumor_allele2_bam, "wb", template=tumor_hla_bam)
            a1, a2 = 0, 0
            for read in tumor_hla_bam.fetch():
                if read.query_name in reads_a1:
                    tumor_allele1_bam_out.write(read)
                    a1 += 1
                elif read.query_name in reads_a2:
                    tumor_allele2_bam_out.write(read)
                    a2 += 1
                else:
                    pass
            logger.info(f'Tumor: allele1 ({len(reads_a1)}, {a1}) | allele2 ({len(reads_a2)}, {a2})')
            tumor_allele1_bam_out.close()
            tumor_allele2_bam_out.close()
            pysam.index(tumor_allele1_bam)
            pysam.index(tumor_allele2_bam)
            tumor_allele1_bed = f'{tumor_dir_methyl}/typing_meth_{hla1_pretty}_{annot_type}_modkit.bed'
            tumor_allele1_log = f'{tumor_dir_methyl}/typing_meth_{hla1_pretty}_{annot_type}_modkit.log'
            cmd_modkit = (f'modkit pileup  --cpg --ref {ref_fasta} --include-bed {hla_bed} --ignore h '
                          f'{tumor_allele1_bam}  {tumor_allele1_bed} 1> {tumor_allele1_log} 2> {tumor_allele1_log};'
                          f'bgzip {tumor_allele1_bed} --force; '
                          f'tabix --force --preset bed {tumor_allele1_bed}.gz')
            logger.debug(cmd_modkit)
            _ = subprocess.run(cmd_modkit, shell=True, capture_output=False)
            modkit_log = open(tumor_allele1_log)
            for line in modkit_log:
                if "rror" in line and "not enough data points" not in line:
                    logger.error(f'Modkit terminated or cannot be run {cmd_modkit}, exiting with error')
                    sys.exit(1)
            modkit_log.close()
            tumor_allele2_bed = f'{tumor_dir_methyl}/typing_meth_{hla2_pretty}_{annot_type}_modkit.bed'
            tumor_allele2_log = f'{tumor_dir_methyl}/typing_meth_{hla2_pretty}_{annot_type}_modkit.log'
            cmd_modkit = (f'modkit pileup  --cpg --ref {ref_fasta} --include-bed {hla_bed} --ignore h '
                          f'{tumor_allele2_bam}  {tumor_allele2_bed}  1> {tumor_allele2_log} 2> {tumor_allele2_log}; '
                          f'bgzip {tumor_allele2_bed} --force; '
                          f'tabix --force --preset bed {tumor_allele2_bed}.gz')
            _ = subprocess.run(cmd_modkit, shell=True, capture_output=False)
            modkit_log = open(tumor_allele2_log)
            for line in modkit_log:
                if "rror" in line and "not enough data points" not in line:
                    logger.error(f'Modkit terminated or cannot be run {cmd_modkit}, exiting with error')
                    sys.exit(1)
            modkit_log.close()
            results_allele_meth[hla_gene][tumor] = [f'{tumor_allele1_bed}.gz', f'{tumor_allele2_bed}.gz']
        else:
            # only one allele pass, extract the methyl results
            hla1 = tumor_hla_reads[hla_gene]
            hla1_pretty = hla_pretty_name(hla1)
            tumor_allele1_bed = f'{tumor_dir_methyl}/typing_meth_{hla1_pretty}_{annot_type}_modkit.bed'
            aln_mod = pysam.TabixFile(f'{tumor_dir}/../methylation/align_mod_modkit.bed.gz')
            file_handler = open(tumor_allele1_bed, "wt")
            for line in aln_mod.fetch(region=reg_fetch):
                file_handler.write(f'{line}\n')
            file_handler.close()
            pysam.tabix_compress(filename_in=tumor_allele1_bed, filename_out=f'{tumor_allele1_bed}.gz', force=True)
            pysam.tabix_index(tumor_allele1_bed, force=True, preset="bed")
            results_allele_meth[hla_gene][tumor] = [f'{tumor_allele1_bed}.gz']
        if 2 == len(control_hla_reads[hla_gene]):
            [hla1, hla2] = control_hla_reads[hla_gene].keys()
            reads_a1 = control_hla_reads[hla_gene][hla1]
            reads_a2 = control_hla_reads[hla_gene][hla2]
            hla1_pretty = hla_pretty_name(hla1)
            hla2_pretty = hla_pretty_name(hla2)
            logger.info(f'{hla1_pretty} / {hla2_pretty}')
            control_allele1_bam = f'{control_dir_methyl}/typing_meth_{hla1_pretty}_{annot_type}.bam'
            control_allele1_bam_out = pysam.AlignmentFile(control_allele1_bam, "wb", template=control_hla_bam)
            control_allele2_bam = f'{control_dir_methyl}/typing_meth_{hla2_pretty}_{annot_type}.bam'
            control_allele2_bam_out = pysam.AlignmentFile(control_allele2_bam, "wb", template=control_hla_bam)
            a1, a2 = 0, 0
            for read in control_hla_bam.fetch():
                if read.query_name in reads_a1:
                    control_allele1_bam_out.write(read)
                    a1 += 1
                elif read.query_name in reads_a2:
                    control_allele2_bam_out.write(read)
                    a2 += 1
                else:
                    pass
            logger.info(f'Control: allele1 ({len(reads_a1)}, {a1}) | allele2 ({len(reads_a2)}, {a2})')
            control_allele1_bam_out.close()
            control_allele2_bam_out.close()
            pysam.index(control_allele1_bam)
            pysam.index(control_allele2_bam)
            control_allele1_bed = f'{control_dir_methyl}/typing_meth_{hla1_pretty}_{annot_type}_modkit.bed'
            cmd_modkit = (f'modkit pileup  --cpg --ref {ref_fasta} --include-bed {hla_bed} --ignore h '
                          f'{control_allele1_bam}  {control_allele1_bed}; '
                          f'bgzip {control_allele1_bed} --force; '
                          f'tabix --force --preset bed {control_allele1_bed}.gz')
            run_cmd = subprocess.run(cmd_modkit, shell=True, capture_output=True, text=True)
            if run_cmd.stderr != "":
                logger.warning(run_cmd.stderr)
                if ("rror" in run_cmd.stderr and "not enough data points" not in run_cmd.stderr and
                        not "zero reads found in bam index"):
                    logger.error(f'Cannot run {cmd_modkit}, exiting with error')
                    sys.exit(1)
            control_allele2_bed = f'{control_dir_methyl}/typing_meth_{hla2_pretty}_{annot_type}_modkit.bed'
            cmd_modkit = (f'modkit pileup  --cpg --ref {ref_fasta} --include-bed {hla_bed} --ignore h '
                          f'{control_allele2_bam}  {control_allele2_bed}; '
                          f'bgzip {control_allele2_bed} --force; '
                          f'tabix --force --preset bed {control_allele2_bed}.gz')
            run_cmd = subprocess.run(cmd_modkit, shell=True, capture_output=True, text=True)
            if run_cmd.stderr != "":
                logger.warning(run_cmd.stderr)
                if ("rror" in run_cmd.stderr and "not enough data points" not in run_cmd.stderr and
                        not "zero reads found in bam index"):
                    logger.error(f'Cannot run {cmd_modkit}, exiting with error')
                    sys.exit(1)
            results_allele_meth[hla_gene][control] = [f'{control_allele1_bed}.gz', f'{control_allele2_bed}.gz']
        else:
            # only one allele pass, extract the methyl results
            hla1 = control_hla_reads[hla_gene]
            hla1_pretty = hla_pretty_name(hla1)
            control_allele1_bed = f'{control_dir_methyl}/typing_meth_{hla1_pretty}_{annot_type}_modkit.bed'
            aln_mod = pysam.TabixFile(f'{control_dir}/../methylation/align_mod_modkit.bed.gz')
            file_handler = open(control_allele1_bed, "wt")
            for line in aln_mod.fetch(region=reg_fetch):
                file_handler.write(f'{line}\n')
            file_handler.close()
            pysam.tabix_compress(filename_in=control_allele1_bed, filename_out=f'{control_allele1_bed}.gz', force=True)
            pysam.tabix_index(control_allele1_bed, force=True, preset="bed")
            results_allele_meth[hla_gene][control] = [f'{control_allele1_bed}.gz']
    logger.info(results_allele_meth)

    for hla_gene in results_allele_meth.keys():
        bed_file = results_allele_meth[hla_gene]["bed"]
        # comparison
        tumor_modkit = results_allele_meth[hla_gene][tumor]
        if 2 == len(tumor_modkit):
            [ta1, ta2] = tumor_modkit
            ta1_name = extract_name_form_path(ta1)
            ta2_name = extract_name_form_path(ta2)
            tumor_methyl = f'--sample {ta1} T_{ta1_name} -s {ta2} T_{ta2_name}'
        else:
            [ta1] = tumor_modkit
            ta1_name = extract_name_form_path(ta1)
            tumor_methyl = f'--sample {ta1} T_{ta1_name}'
        control_modkit = results_allele_meth[hla_gene][control]
        if 2 == len(control_modkit):
            [ca1, ca2] = control_modkit
            ca1_name = extract_name_form_path(ca1)
            ca2_name = extract_name_form_path(ca2)
            control_methyl = f'--sample {ca1} C_{ca1_name} -s {ca2} C_{ca2_name}'
        else:
            [ca1] = control_modkit
            ca1_name = "_".join(ca1.split("/")[-1].split("_")[2:4])
            control_methyl = f'--sample {ca1} C_{ca1_name}'
        modkit_dmr_cmd = (f"modkit dmr multi {control_methyl}  {tumor_methyl}  "
                          f"--regions-bed {bed_file} --ref {ref_fasta} --base C "
                          f"--out-dir {hla_type_dmr_dir}/_compare_methyl_by_type_hla{hla_gene}_{annot_type} "
                          f"--log-filepath {hla_type_dmr_dir}/_compare_methyl_by_type_hla{hla_gene}_{annot_type}.log "
                          f"--force")
        logger.debug(f'Running CMD: {modkit_dmr_cmd}')
        run_cmd = subprocess.run(f'{modkit_dmr_cmd}', shell=True, capture_output=True, text=True)
        if run_cmd.stderr != "":
            logger.warning(run_cmd.stderr)
            if "rror" in run_cmd.stderr:
                logger.error(f'Error while running modkit dmr multi: "{modkit_dmr_cmd}"')
                sys.exit(1)

    return results_allele_meth


def get_arguments() -> TypeArgs:
    filter_split_reads_help = """
    Filter_split_reads:
            --sample/-s          Sample name
            --work-dir/-w        Working directory path (typing)
    """
    parser = argparse.ArgumentParser(description="", usage=filter_split_reads_help)

    # ############################################################################################ #
    parser.add_argument("-s", "--sample", type=str, required=True, dest='sample', default="", help='')
    parser.add_argument("-w", "--work-dir", type=str, required=True, dest='work_dir', default="",
                        help='Working directory path')
    parser.add_argument("-d", "--data-dir", type=str, required=True, dest='data_dir', default="",
                        help='Directory path of typing results')
    parser.add_argument("-q", "--min-query-aln", type=int, required=True, dest='min_query_aln',
                        default=0, help='')
    parser.add_argument("-m", "--max-edit-dist", type=float, required=True, dest='max_edit_dist',
                        default=0.0, help='')
    parser.add_argument("-e", "--use-edit-dist", action='store_true', required=False, dest='use_edit_dist',
                        default=False, help='')
    user_args_parsed = parser.parse_args()
    return TypeArgs(user_args_parsed.sample, user_args_parsed.work_dir, user_args_parsed.min_query_aln,
                    user_args_parsed.max_edit_dist, user_args_parsed.use_edit_dist, user_args_parsed.data_dir)


def main(user_arguments: TypeArgs = None) -> tuple[dict, dict]:
    if user_arguments is not None:
        use_args = user_arguments
    else:
        use_args = get_arguments()
    hla_types = get_types(sample_name=use_args.sample, work_dir=use_args.data_dir)

    results = {}
    hla_reads_for_methl = {}

    logger.debug(f'Read fasta with HLA types + seq')
    hla_fasta = pysam.FastaFile(f'{use_args.annot_data}/{hla_types.hla_fasta}')
    hla_fasta_v2 = pysam.FastaFile(f'{use_args.annot_data}/{hla_types.hla_fasta4}')

    if not os.path.exists(user_arguments.work_dir):
        logger.info(f'making: {user_arguments.work_dir}')
        os.mkdir(user_arguments.work_dir)
    fastq_split_read = FastqSplit()

    for hla_gene in hla_types.gene_names:
        if hla_gene not in results:
            results[hla_gene] = None
            hla_reads_for_methl[hla_gene] = None
        logger.info(f'Gene: HLA-{hla_gene}')
        hla_genes_lst = hla_types.get_type(hla_gene)
        hla_genes_str = ",".join(hla_genes_lst)
        hla1, hla2 = hla_genes_lst
        if not_same_type(hla_genes_lst):
            logger.debug(f'Split LR into {fastq_split_read.new_read_len}bp')
            hla_lr = pysam.FastxFile(f'{hla_types.dir_work}/{hla_types.get_long_reads(hla_gene)}')
            hla_reads_split_fq = f'{user_arguments.work_dir}/{hla_types.get_long_reads_split(hla_gene)}'
            if not os.path.exists(hla_reads_split_fq):
                hla_lr_split = gzip.open(hla_reads_split_fq, "wt")
                for each_record in hla_lr:
                    fastq_split_read.process_read(each_record)
                result_split_reads = "\n".join(fastq_split_read.output)
                hla_lr_split.write(f'{result_split_reads}\n')
                hla_lr.close()
                hla_lr_split.close()

            logger.debug(f'Extract HLA type fasta split LR for: {hla_genes_str}')
            hla_gene_fasta = extract_needed_hla(hla_fasta, user_arguments.work_dir, hla_gene,
                                                hla_genes_lst, hla_fasta_v2)
            logger.debug(f'Map split LR {hla_genes_lst}')
            hla_bam = align_hla(f'HLA-{hla_gene}', user_arguments.work_dir,
                                [f'{hla_gene_fasta}', hla_reads_split_fq, f'{hla_gene}.bam'])
            reads_mapped, reads_by_type = read_analysis(hla_bam, use_args.min_query_aln, use_args.max_edit_dist,
                                                        {hla1: [], hla2: []}, use_args.use_edit_dist)
            allele1 = len(reads_by_type[hla1])
            allele2 = len(reads_by_type[hla2])
            if reads_mapped > 0:
                logger.debug(f'{hla1}: {(allele1/reads_mapped):.4f} ({allele1}/{reads_mapped})')
                logger.debug(f'{hla2}: {(allele2/reads_mapped):.4f} ({allele2}/{reads_mapped})')
            else:
                logger.debug(f'{hla1}, {hla2}: no reads mapped')
            results[hla_gene] = (hla1, allele1, hla2, allele2, "HET")
            hla_reads_for_methl[hla_gene] = {hla1: reads_by_type[hla1], hla2: reads_by_type[hla2]}
        else:
            logger.debug(f'{hla_genes_str} are the same type')
            results[hla_gene] = (hla1, 1, hla2, 0, "HOM")
            hla_reads_for_methl[hla_gene] = hla1
    logger.debug(results, hla_reads_for_methl)
    return results, hla_reads_for_methl


# main
if __name__ == '__main__':
    main()
