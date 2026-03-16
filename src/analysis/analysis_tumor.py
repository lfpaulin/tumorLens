#!/usr/bin/env python3
import os
import sys
import pysam
import numpy as np
from matplotlib import pyplot as plot_engine
from matplotlib import gridspec
from utils.logger_debug import setup_log
from classes.analysis_func import load_hla_annotation
from classes.analysis_func import HLAType


# Globals
logger = setup_log(__name__, debug=True)
hla_region = {"chromosome": "chr6", "start": 29e6, "end": 34e6}


class AnalysisFiles(object):
    # init
    def __init__(self):
        # files/dir
        self.work_dir = ""    # Directory where the analysis is performed, must exist
        self.fasta = ""       # Reference genome
        self.annot = ""       # Annotations file (genes) for same version of genome
        self.sv = ""          # SV file with results (vcf)
        self.snv = ""         # SNV file with results (vcf)
        self.cnv = ""         # CNV file with results (vcf)
        self.coverage = ""    # Coverage file with results (bed)
        self.methyl = ""      # Methylation calls (?)
        self.hlatype = ""     # Typing results for HLA genes
        self.out_prefix = ""  # For output
        self.output_report = ""
        self.output_file = ""
        self.names = ""
        # file handlers
        self.handler_report = None   # open W
        self.handler_output = None   # open W
        self.handler_fasta = None    # pysam.FAstaFile
        self.handler_annot = None    # open R
        self.handler_sv = None       # pysam.VariantFile
        self.handler_snv = None      # pysam.VariantFile
        self.handler_cnv = None      # pysam.VariantFile
        self.handler_cov = None      # pysam.Tabix
        self.handler_methyl = None   # pysam.Tabix
        self.handler_hlatype = None  # open R

    def args_check_tumor_only(self, user_args):
        self_cmd = "AnalysisFiles.args_check"
        # errors 0x194 -> not found
        #        0x1F7 -> could not open
        # Work dir ###########################################################
        if os.path.exists(user_args.work_dir):
            self.work_dir = user_args.work_dir
        else:
            logger.error(f'{user_args.work_dir} not found and needed. ID:0x194-1')
        # Fasta    ###########################################################
        if os.path.exists(user_args.fasta):
            self.fasta = user_args.fasta
            logger.info(f'using {self.fasta} reference')
            try:
                # .references & .lengths
                self.handler_fasta = pysam.FastaFile(self.fasta)
            except OSError as os_err:
                logger.error(f'could not open {self.fasta} reference.\n{os_err}\nID:0x1F7-2')
        else:
            logger.error(f'{user_args.fasta} not found and needed. ID:0x194-2')
        # Annot    ###########################################################
        if os.path.exists(user_args.annot):
            self.annot = user_args.annot
            logger.info(f'using {self.annot} annotation')
            try:
                # data is bed, chr,start,end,gene_name,strand,gene_id
                self.handler_annot = open(self.annot, "r")
            except OSError as os_err:
                logger.error(f'could not open {self.annot} annotation.\n{os_err}\nID:0x1F7-3')
        else:
            logger.error(f'{user_args.annot} not found and needed. ID:0x194-3')
        # SV       ###########################################################
        if os.path.exists(user_args.sv):
            self.sv = user_args.sv
            logger.info(f'using {self.sv} SV')
            try:
                self.handler_sv = pysam.VariantFile(self.sv)
            except OSError as os_err:
                logger.error(f'could not open {self.sv}.\n{os_err}\nID:0x1F7-4')
        else:
            logger.error(f'{user_args.sv} not found and needed. ID:0x194-4')
        # SNV      ###########################################################
        if os.path.exists(user_args.snv):
            self.snv = user_args.snv
            logger.info(f'using {self.snv} SNV')
            try:
                self.handler_snv = pysam.VariantFile(self.snv)
            except OSError as os_err:
                logger.error(f'could not open {self.snv}.\n{os_err}\nID:0x1F7-5')
        else:
            logger.error(f'{user_args.snv} not found and needed. ID:0x194-5')
        # CNV      ###########################################################
        if os.path.exists(user_args.cnv):
            self.cnv = user_args.cnv
            logger.info(f'using {self.cnv} CNV')
            try:
                self.handler_cnv = pysam.VariantFile(self.cnv)
            except OSError as os_err:
                logger.error(f'could not open {self.cnv}.\n{os_err}\nID:0x1F7-6')
        else:
            logger.error(f'{user_args.cnv} not found and needed. ID:0x194-6')
        # Coverage ###########################################################
        if os.path.exists(user_args.cov):
            self.coverage = user_args.cov
            logger.info(f'using {self.coverage} coverage data')
            try:
                self.handler_cov = pysam.TabixFile(self.coverage)
            except OSError as os_err:
                logger.error(f'could not open {self.coverage} coverage.\n{os_err}\nID:0x1F7-7')
        else:
            logger.error(f'{user_args.cov} not found and needed. ID:0x194-7')
        # Methyl   ###########################################################
        if os.path.exists(user_args.methyl):
            self.methyl = user_args.methyl
            logger.info(f'using {self.methyl} methylation')
            try:
                self.handler_methyl = pysam.TabixFile(self.methyl)
            except OSError as os_err:
                logger.error(f'could not open {self.methyl} methylation.\n{os_err}\nID:0x1F7-8')
        else:
            logger.error(f'{user_args.methyl} not found and needed. ID:0x194-8')
        # HLA typing   #######################################################
        if os.path.exists(user_args.hlatype):
            self.hlatype = user_args.hlatype
            logger.info(f'using {self.hlatype} HLA typing')
            try:
                self.handler_hlatype = open(self.hlatype, "r")
            except OSError as os_err:
                logger.error(f'could not open {self.hlatype}  HLA typing.\n{os_err}\nID:0x1F7-9')
        else:
            logger.error(f'{user_args.hlatype} not found and needed. ID:0x194-9')
        # Output   ###########################################################
        # directory check
        if not os.path.isdir(f'{user_args.work_dir}'):
            os.mkdir(f'{user_args.work_dir}')
        self.out_prefix = user_args.out_prefix
        # Report
        self.output_report = f'{user_args.work_dir}/{self.out_prefix}_report.html'
        self.handler_report = open(self.output_report, "wt")
        logger.info(f'using {self.output_report} as output')
        # Output
        self.output_file = f'{user_args.work_dir}/{self.out_prefix}_results.txt'
        self.handler_output = open(self.output_file, "wt")
        logger.info(f'using {self.output_file} as output')
        # End      ###########################################################

    def close_handlers_tumor_only(self):
        # Fasta    ###########################################################
        self.handler_fasta.close() if (self.handler_fasta is not None and not self.handler_fasta.closed) else None
        # Annot    ###########################################################
        self.handler_annot.close() if (self.handler_annot is not None and not self.handler_annot.closed) else None
        # SV       ###########################################################
        self.handler_sv.close() if (self.handler_sv is not None and not self.handler_sv.closed) else None
        # SNV      ###########################################################
        self.handler_snv.close() if (self.handler_snv is not None and not self.handler_snv.closed) else None
        # CNV      ###########################################################
        self.handler_cnv.close() if (self.handler_cnv is not None and not self.handler_cnv.closed) else None
        # Coverage ###########################################################
        self.handler_cov.close() if (self.handler_cov is not None and not self.handler_cov.closed) else None
        # Methyl   ###########################################################
        self.handler_methyl.close() if (self.handler_methyl is not None and not self.handler_methyl.closed) else None
        # HLA typing   #######################################################
        self.handler_hlatype.close() if (self.handler_hlatype is not None and not self.handler_hlatype.closed) else None
        # Report   ###########################################################
        self.handler_report.close() if (self.handler_report is not None and not self.handler_report.closed) else None
        # Output   ###########################################################
        self.handler_output.close() if (self.handler_output is not None and not self.handler_output.closed) else None
        # End      ###########################################################


class HLAPlot:
    def __init__(self):
        # the plot
        self.figure = plot_engine.figure(figsize=(9, 9))
        gs = gridspec.GridSpec(4, 1, height_ratios=[3, 3, 2, 1])
        self.cov_plot = plot_engine.subplot(gs[0])
        self.meth_plot = plot_engine.subplot(gs[1])
        self.variants = plot_engine.subplot(gs[2])
        self.gene_body = plot_engine.subplot(gs[3])
        self.gene_body.axes.get_yaxis().set_visible(False)
        # colors
        self.coverage_color = "#377eb8"
        self.methyl_color = "#ff7f00"
        self.variant_color = {"loh": "#de4e83", "cnv": "#009999", "sv": "#6600cc"}
        self.gene_color = "#000000"
        self.gene_color_body = {"body": "#252525", "promoter": "#969696"}
        self.gene_spacing_color = "#000000"
        # legends and axis
        self.axis_ylim = {"bottom": 0, "top": 0}  # to init
        self.file_prefix = "Analysis"
        self.output_directory = "./"

    def plot_hla_region(self, chr_region="", region_info=None, hla_genes=None):
        self_cmd = "HLAPlot.plot_hla_region"
        if region_info is None:
            logger.error("'region_info' parameter is needed")
        logger.info(f'plotting results for {chr_region}')
        # reference region
        [region_start, region_end] = region_info.position
        # coverage plot
        coverage = region_info.cov["plot"]
        self.axis_ylim["top"] = np.max(coverage["cov"]) + 1
        self.cov_plot.axes.set_ylabel(f'Coverage')
        self.cov_plot.plot(np.array(coverage["pos"]), np.array(coverage["cov"]),
                           color=self.coverage_color, linewidth='1')
        self.cov_plot.axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
        self.cov_plot.axes.set_xlim(left=region_start, right=region_end)
        # methylation plot
        methylation = region_info.meth["plot"]
        self.axis_ylim["bottom"] = -10
        self.axis_ylim["top"] = 110
        self.meth_plot.axes.set_ylabel(f'Methylation proportion')
        self.meth_plot.scatter(np.array(methylation["pos"]), np.array(methylation["meth"]),
                               color=self.methyl_color, marker="o", s=2)
        self.meth_plot.axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
        self.meth_plot.axes.set_xlim(left=region_start, right=region_end)
        # Variants plot
        variants = {"loh": region_info.loh["full"], "cnv": region_info.cnv["full"], "sv": region_info.sv["full"]}
        self.axis_ylim["bottom"] = 0
        self.axis_ylim["top"] = 4
        self.variants.axes.set_ylabel(f'Variants')
        variant_y_values = {"sv": 1, "cnv": 2, "loh": 3}
        for each_variant_type in variants.keys():
            logger.info(f'Number of {each_variant_type}: {len(variants[each_variant_type])}')
            var_counter = 0
            for each_variant in variants[each_variant_type]:
                [start, end] = each_variant.start, each_variant.stop
                var_y_value = variant_y_values[each_variant_type] + (0 if var_counter % 3 == 0 else
                                                                     (0.2 if var_counter % 3 == 1 else -0.2))
                self.variants.plot(np.array([start, end]), np.array([var_y_value, var_y_value]), linewidth='2',
                                   color=self.variant_color[each_variant_type])
                var_counter += 1
        self.variants.axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
        self.variants.axes.set_xlim(left=region_start, right=region_end)
        # Genetic element body
        self.gene_body.axes.set_ylabel(f'HLA')
        if hla_genes is not None:
            genn_counter = 0
            for each_gene in hla_genes.keys():
                logger.info(f'Result of gene: {each_gene}')
                [start, end] = hla_genes[each_gene].position
                self.gene_body.plot(np.array([start, end]), np.array([0, 0]), linewidth='2', color=self.gene_color)
                genn_counter += 1
            self.gene_body.axes.set_xlabel(f'Chromosome position')
            self.gene_body.axes.set_xlim(left=region_start, right=region_end)
        # save and close
        chr_region_pretty = "_".join("_".join(chr_region.split(":")).split("-"))
        self.figure.suptitle(f'{self.file_prefix} HLA Locus {chr_region}')
        self.figure.savefig(f'{self.output_directory}/plot-hla-{self.file_prefix}-{chr_region_pretty}.png', dpi=300)
        logger.info(f'Plot saved: plot-hla-{self.file_prefix}-{chr_region_pretty}.png')
        self.figure.clf()

    def plot_hla_result(self, gene_name="", gene_info=None, gene_padding=2e3, promoter=1e3):
        self_cmd = "HLAPlot.plot_hla_result"
        if gene_info is None:
            logger.error("'gene_info' parameters are needed")
        logger.info(f'plotting results for {gene_name}')
        # gene info
        [gene_start, gene_end] = gene_info.position
        gene_start -= gene_padding
        gene_end += gene_padding
        coverage = gene_info.cov["plot"]
        methylation = gene_info.meth["plot"]
        # coverage plot
        self.axis_ylim["top"] = np.max(coverage["cov"]) + 1
        self.axis_ylim["bottom"] = 0
        self.cov_plot.axes.set_ylabel(f'Coverage')
        self.cov_plot.plot(np.array(coverage["pos"]), np.array(coverage["cov"]),
                           color=self.coverage_color, linewidth='1')
        self.cov_plot.axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
        self.cov_plot.axes.set_xlim(left=gene_start, right=gene_end)
        # meth plot
        self.axis_ylim["bottom"] = -5
        self.axis_ylim["top"] = 105
        self.meth_plot.axes.set_ylabel(f'Methylation proportion')
        self.meth_plot.scatter(np.array(methylation["pos"]), np.array(methylation["meth"]),
                               color=self.methyl_color, marker="o", s=4)
        self.meth_plot.axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
        self.meth_plot.axes.set_xlim(left=gene_start, right=gene_end)
        # Drawing the genetic element body
        [start, end] = gene_info.position
        logger.info(f"Result of gene: {gene_name}, Gene body: {gene_info.position}")
        self.gene_body.plot(np.array([start, end]), np.array([0, 0]), linewidth='5', color=self.gene_color_body["body"])
        # Ikb upstream as promoter
        start, end = gene_info.start - promoter, gene_info.start
        self.gene_body.axes.set_xlabel(f'Chromosome position')
        self.gene_body.axes.set_xlim(left=gene_start, right=gene_end)
        self.gene_body.plot(np.array([start, end]), np.array([0, 0]), linewidth='3',
                            color=self.gene_color_body["promoter"])
        # Arrow for strand  -> or <-
        logger.info(f'{gene_info.strand}  {gene_info.strand == "+"}', f'{self_cmd} err')
        if gene_info.strand == "+":
            start = gene_info.end
            end = gene_padding/8
        else:
            start = gene_info.start - promoter
            end = -gene_padding/8
        self.gene_body.arrow(x=start, y=0, dx=end, dy=0, head_width=0.02, head_length=(gene_padding/8)*0.95,
                             length_includes_head=True, color=self.gene_spacing_color)
        # variants per gene
        variants = {"loh": gene_info.loh["full"], "cnv": gene_info.cnv["full"], "sv": gene_info.sv["full"]}
        self.axis_ylim["bottom"] = 0
        self.axis_ylim["top"] = 4
        self.variants.axes.set_ylabel(f'Variants')
        variant_y_values = {"sv": 1, "cnv": 2, "loh": 3}
        for each_variant_type in variants.keys():
            logger.info(f'Number of {each_variant_type}: {len(variants[each_variant_type])}')
            var_counter = 0
            for each_variant in variants[each_variant_type]:
                [start, end] = each_variant.start, each_variant.stop
                var_y_value = variant_y_values[each_variant_type] + (0 if var_counter % 3 == 0 else
                                                                     (0.2 if var_counter % 3 == 1 else -0.2))
                self.variants.plot(np.array([start, end]), np.array([var_y_value, var_y_value]), linewidth='2',
                                   color=self.variant_color[each_variant_type])
                var_counter += 1
        self.variants.axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
        self.variants.axes.set_xlim(left=gene_start, right=gene_end)
        # save and close
        self.figure.suptitle(f'{self.file_prefix} gene: {gene_name}')
        self.figure.savefig(f'{self.output_directory}/plot-hla-{self.file_prefix}-{gene_name}.png', dpi=300)
        logger.info(f'Plot saved: plot-hla-{self.file_prefix}-{gene_name}.png')
        self.figure.clf()


class HLALocusResults(object):
    def __init__(self, hla_region_dict):
        self.chromosome = hla_region_dict["chromosome"]
        self.start = int(hla_region_dict["start"])
        self.end = int(hla_region_dict["end"])
        self.position = [self.start, self.end]
        self.region_coordinates = f'{self.chromosome}:{self.start}-{self.end}'
        # also add results
        self.snv = None
        self.loh = None
        self.cov = None
        self.cnv = None
        self.sv = None
        self.meth = None

    def set_results(self, gene_results):
        self.snv = gene_results["snv"]
        self.loh = gene_results["loh"]
        self.cov = gene_results["coverage"]
        self.cnv = gene_results["cnv"]
        self.sv = gene_results["sv"]
        self.meth = gene_results["methylation"]

    def summary_hla_region(self, output_handler):
        pretty_cov = np.round(self.cov["summary"], 2) if self.cov["summary"] is not None else None
        pretty_mth = np.round(self.meth["summary"], 2) if self.meth["summary"] is not None else None
        sv = f'{self.sv["summary"]}/{len(self.sv["full"])}'
        loh = f'{self.loh["summary"]}/{len(self.loh["full"])}'
        cnv = f'{self.cnv["summary"]}/{len(self.cnv["full"])}'
        log_line = f'{self.region_coordinates} => Coverage: {pretty_cov}\tMethylation: {pretty_mth}\t' \
                   f'SV: {sv}\tCNV: {cnv}\tLoH: {loh}'
        out_line = f'HLA_locus\t{pretty_cov}\t{pretty_mth}\t{sv}\t{cnv}\t{loh}'
        logger.info(log_line, "HLALocusResults.out")
        output_handler.write(f'{out_line}\n')


def load_results(analysis_files=None, chromosome=None, start=None, end=None, padding=0, gene_name=""):
    analysis_results = {"methylation": None, "snv": None, "sv": None, "cnv": None, "coverage": None}
    start = start - padding
    end = end + padding
    region_fetch = f'{chromosome}:{start}-{end}'

    # Structural Variants (SV)
    # Note: analysis_files.handler_sv = None        # pysam.VariantFile
    var_list = []
    for variant in analysis_files.handler_sv.fetch(region=region_fetch):
        var_list.append(variant)
    sv_summary = True if len(var_list) > 0 else False
    analysis_results["sv"] = {"summary": sv_summary, "full": var_list}

    # Single Nucleotide Variants (SNV)
    # Note: analysis_files.handler_snv = None       # pysam.VariantFile
    # Note: This one is used for LoH inside Spectre (CNV)

    # Copy Number Variants (CNV)
    # Note: analysis_files.handler_cnv = None       # pysam.VariantFile
    var_list = []
    loh_list = []
    for variant in analysis_files.handler_cnv.fetch(region=region_fetch):
        if "LOH" in variant.alts[0]:
            loh_list.append(variant)
        else:
            var_list.append(variant)
    cnv_summary = True if len(var_list) > 0 else False
    loh_summary = True if len(loh_list) > 0 else False
    # CNV and LoH based on CNV
    analysis_results["cnv"] = {"summary": cnv_summary, "full": var_list}
    analysis_results["loh"] = {"summary": loh_summary, "full": loh_list}

    # Coverage (Used for CNV, also plot)
    # Note: analysis_files.handler_cov = None       # pysam.TabixFile
    var_list = []
    for_plot = {"pos": [], "cov": []}
    for variant in analysis_files.handler_cov.fetch(region=region_fetch):
        [conitg, start, end, cover] = variant.rstrip("\n").split("\t")
        cover_float = float(cover)
        var_list.append([conitg, int(start), int(end), cover_float])
        for_plot["pos"].append(int(start))
        for_plot["cov"].append(cover_float)
    analysis_results["coverage"] = {"full": var_list, "summary": np.nanmean(for_plot["cov"]), "plot": for_plot} \
        if len(var_list) > 0 else {"full": None, "summary": None, "plot": for_plot}

    # Methylation modification 5mC (proxy for gene expression)
    # Note: analysis_files.handler_methyl = None    # pysam.TabixFile
    # 1	chrom	name of reference sequence from BAM header	str
    # 2	start position	0-based start position	int
    # 3	end position	0-based exclusive end position	int
    # 10	Nvalid_cov	See definitions above.	int
    # 11	fraction modified	Nmod / Nvalid_cov	float
    # 12	Nmod	See definitions above.	int
    # 13	Ncanonical	See definitions above.	int
    var_list = []
    for_plot = {"pos": [], "meth": []}
    min_coverage_use_meth = 10
    for variant in analysis_files.handler_methyl.fetch(region=region_fetch):
        [conitg, start, _, _, _, _, _, _, _, mth_info] = variant.rstrip("\n").split("\t")
        [valid, mod_frac, _, _, _, _, _, _, _] = mth_info.split(" ")
        if int(valid) > min_coverage_use_meth:
            mod_frac_float = float(mod_frac)
            var_list.append([conitg, int(start), mod_frac_float])
            for_plot["pos"].append(int(start))
            for_plot["meth"].append(mod_frac_float)
    analysis_results["methylation"] = {"full": var_list, "summary": np.nanmean(for_plot["meth"]), "plot": for_plot} \
        if len(var_list) > 0 else {"full": None, "summary": None, "plot": for_plot}

    # return results
    return analysis_results


# Note: TUMOR ONLY
# TODO:
def tumor(use_args):
    analysis_files = AnalysisFiles()
    analysis_files.args_check_tumor_only(use_args)
    logger.info(f'Sample: {use_args.out_prefix}')
    analysis_files.handler_output.write(f'Sample: {use_args.out_prefix}\n')
    """
    Note: all files inputs used
    analysis_files.handler_output   -> <OUT>_report.html
    analysis_files.handler_fasta    -> <IN>.fasta.gz | <IN>.fasta.gz.fai
    analysis_files.handler_annot    -> IN>.bed
    analysis_files.handler_sv       -> <IN>-sv.vcf.gz
    analysis_files.handler_snv      -> <IN>-snv-clair3.vcf.gz
    analysis_files.handler_cnv      -> <IN>-cnv-spectre.vcf.gz
    analysis_files.handler_cov      -> <IN>-coverage-mosdepth.bed.gz
    analysis_files.handler_methyl   -> <IN>_methylation.bed.gz
    analysis_files.handler_hlatype  -> <IN>_hla_typing.bed
    """
    # Section: Summary ###########################
    # ### Section: HLA region ####################
    hla_region_results = HLALocusResults(hla_region)
    logger.info("HLA Cov\tHLA 5mC %\tHLA SV 50+/n\tHLA CNV 50kb+/n\tHLA LoH 1Mb+/n")
    analysis_files.handler_output.write("HLA Cov\tHLA 5mC %\tHLA SV 50+/n\tHLA CNV 50kb+/n\tHLA LoH 1Mb+/n\n")
    hla_region_results.set_results(load_results(analysis_files, hla_region_results.chromosome, hla_region_results.start,
                                                hla_region_results.end))
    hla_region_results.summary_hla_region(analysis_files.handler_output)
    # ### Section: Genes of interest #############
    # Genes in the bed file list
    # load bed gene info
    gene_info, hla_info = load_hla_annotation(analysis_files.handler_annot)
    # HLA typing (proxy for LoH)
    # analysis_files.handler_hlatype = None   # open
    hlatype = HLAType()
    hlatype.set_types(analysis_files.handler_hlatype)
    # padding
    padding_use = 1000
    windows_size = 1000
    # for each gene
    logger.info("HLA Cov\tHLA 5mC %\tHLA SV 50+\tHLA CNV 50kb+\tHLA LoH SNV50kb+/Type")
    analysis_files.handler_output.write("HLA Cov\tHLA 5mC %\tHLA SV 50+\tHLA CNV 50kb+\tHLA LoH SNV50kb+/Type\n")
    for use_gene in gene_info.keys():
        [gene_start, gene_end] = gene_info[use_gene].position
        gene_results = load_results(analysis_files, gene_info[use_gene].chromosome, gene_start, gene_end, padding_use, use_gene)
        gene_info[use_gene].set_results(gene_results)
        if "HLA" in use_gene:
            [_, hla_gene] = use_gene.split("-")
            if hla_gene in hlatype.gene_names:
                gene_info[use_gene].set_hla(hlatype.hla[hla_gene])
        # summary
        gene_info[use_gene].summary_hla(analysis_files.handler_output, hla_only=True)
    # Section: END Summary #######################

    # Section: Plots #############################
    plot_region = HLAPlot()
    plot_region.output_directory = use_args.work_dir
    # Region
    plot_region.plot_hla_region(chr_region=hla_region_results.region_coordinates, region_info=hla_region_results,
                                hla_genes=hla_info)
    # HLA genes
    for each_hla in hla_info.keys():
        plot_hla = HLAPlot()
        plot_hla.output_directory = use_args.work_dir
        plot_hla.plot_hla_result(gene_name=each_hla, gene_info=hla_info[each_hla])
    # Section: END Plots #########################

    # DEV
    analysis_files.close_handlers_tumor_only()
    sys.exit(0)

    # Aqui ####################
    if False:
        gene_start -= padding_use
        gene_end += padding_use
        # get coverage
        [_cov, _pos] = [[], []]
        for cov in analysis_files.handler_cov.fetch(gene_info[use_gene].chromosome, gene_start, gene_end):
            [_, _cov_start, _, _cov_val] = cov.split("\t")
            _pos.append(int(_cov_start))
            _cov.append(float(_cov_val))
        coverage_gene = {"pos": np.array(_pos), "cov": np.array(_cov)}
        # get meth
        [_meth, _pos] = [[], []]
        for me in analysis_files.handler_methyl.fetch(gene_info[use_gene].chromosome, gene_start, gene_end):
            [_, pos, _, _, _, _, _, _, _, _, p_meth, _, _, _] = me.split("\t")
            _pos.append(int(pos))
            _meth.append(float(p_meth))
        #methylation = {"pos": _pos, "ratio": _meth}
        # plot results
        plot_hla = HLAPlot()
        plot_hla.output_directory = use_args.work_dir
        plot_hla.plot_hla_result(gene_name=use_gene, gene_info=gene_info[use_gene])
    # Region analysis
    # get region coverage
    [_cov, _pos] = [[], []]
    for cov in analysis_files.handler_cov.fetch(hla_region["chromosome"], hla_region["start"], hla_region["end"]):
        [_, _cov_start, _, _cov_val] = cov.split("\t")
        _pos.append(int(_cov_start))
        _cov.append(float(_cov_val))
    coverage_region = {"pos": np.array(_pos), "cov": np.array(_cov)}
    # get region methylation
    [_meth, _pos] = [[], []]
    for me_pos in range(int(hla_region["start"]), int(hla_region["end"])+1, int(windows_size)):
        [tmp_meth, tmp_pos] = [[], []]
        for me in analysis_files.handler_methyl.fetch(hla_region["chromosome"], me_pos, int(me_pos+windows_size)):
            # [_, pos, _, _, qc, _, _, _, _, tot, p_meth, c_canon, c_meth, c_rm] = me.split("\t")
            [_, pos, _, _, _, _, _, _, _, _, p_meth, _, _, _] = me.split("\t")
            tmp_pos.append(int(pos))
            tmp_meth.append(float(p_meth))
        if len(tmp_pos) > 0:
            _pos.append(tmp_pos[0])
            _meth.append(np.nanmean(np.array(tmp_meth)))
        else:
            _pos.append(me_pos)
            _meth.append(np.nan)
    #methylation = {"pos": _pos, "ratio": _meth}
    hla_plot = HLAPlot()
    hla_plot.output_directory = use_args.work_dir
    hla_plot.plot_hla_region(chr_region="Chr6:29000000-34000000", gene_info=gene_info)
    # close
    analysis_files.close_handlers_tumor_only()


def tumor_api(pipeline_files, use_args):
    self_cmd = 'loss'
    analysis_files = AnalysisFiles()
    analysis_files.args_check_tumor_only(use_args)
    logger.info(f'Sample: {use_args.out_prefix}')
    analysis_files.handler_output.write(f'Sample: {use_args.out_prefix}\n')
    # Note: all files inputs use
    # analysis_files.handler_output   -> <OUT>_report.html
    # analysis_files.handler_fasta    -> <IN>.fasta.gz | <IN>.fasta.gz.fai
    # analysis_files.handler_annot    -> IN>.bed
    # analysis_files.handler_sv       -> <IN>-sv.tsv
    # analysis_files.handler_snv      -> <IN>-snv-clair3.tsv
    # analysis_files.handler_cnv      -> <IN>-cnv-spectre.vcf
    # analysis_files.handler_cov      -> <IN>-coverage-mosdepth.bed.gz
    # analysis_files.handler_methyl   -> <IN>_methylation.bed.gz
    # analysis_files.handler_hlatype  -> <IN>_hla_typing.bed
    # Section: Summary ###########################
    # ### Section: HLA region ####################
    hla_region_results = HLALocusResults(hla_region)
    logger.info("HLA Cov\tHLA 5mC %\tHLA SV 50+/n\tHLA CNV 50kb+/n\tHLA LoH SNV50kb+/n")
    analysis_files.handler_output.write("HLA Cov\tHLA 5mC %\tHLA SV 50+/n\tHLA CNV 50kb+/n\tHLA LoH SNV50kb+/n\n")
    hla_region_results.set_results(load_results(analysis_files, hla_region_results.chromosome, hla_region_results.start,
                                                hla_region_results.end))
    hla_region_results.summary_hla_region(analysis_files.handler_output)
    # ### Section: Genes of interest #############
    # Genes in the bed file list
    # load bed gene info
    gene_info, hla_info = load_hla_annotation(analysis_files.handler_annot)
    # HLA typing (proxy for LoH)
    # analysis_files.handler_hlatype = None   # open
    hlatype = HLAType()
    hlatype.set_types(analysis_files.handler_hlatype)
    # padding
    padding_use = 1000
    windows_size = 1000
    # for each gene
    logger.info("HLA Cov\tHLA 5mC %\tHLA SV 50+\tHLA CNV 50kb+\tHLA LoH SNV50kb+/Type")
    analysis_files.handler_output.write("HLA Cov\tHLA 5mC %\tHLA SV 50+\tHLA CNV 50kb+\tHLA LoH SNV50kb+/Type\n")
    for use_gene in gene_info.keys():
        [gene_start, gene_end] = gene_info[use_gene].position
        gene_results = load_results(analysis_files, gene_info[use_gene].chromosome, gene_start, gene_end, padding_use, use_gene)
        gene_info[use_gene].set_results(gene_results)
        if "HLA" in use_gene:
            [_, hla_gene] = use_gene.split("-")
            if hla_gene in hlatype.gene_names:
                gene_info[use_gene].set_hla(hlatype.hla[hla_gene])
        # summary
        gene_info[use_gene].summary_hla(analysis_files.handler_output, hla_only=True)
    # Section: END Summary #######################

    # Section: Plots #############################
    plot_region = HLAPlot()
    plot_region.output_directory = use_args.work_dir
    # Region
    plot_region.plot_hla_region(chr_region=hla_region_results.region_coordinates, region_info=hla_region_results,
                                hla_genes=hla_info)
    # HLA genes
    for each_hla in hla_info.keys():
        plot_hla = HLAPlot()
        plot_hla.output_directory = use_args.work_dir
        plot_hla.plot_hla_result(gene_name=each_hla, gene_info=hla_info[each_hla])
    # Section: END Plots #########################

    # DEV
    analysis_files.close_handlers_tumor_only()
    sys.exit(0)

    # Aqui ####################
    if False:
        gene_start -= padding_use
        gene_end += padding_use
        # get coverage
        [_cov, _pos] = [[], []]
        for cov in analysis_files.handler_cov.fetch(gene_info[use_gene].chromosome, gene_start, gene_end):
            [_, _cov_start, _, _cov_val] = cov.split("\t")
            _pos.append(int(_cov_start))
            _cov.append(float(_cov_val))
        coverage_gene = {"pos": np.array(_pos), "cov": np.array(_cov)}
        # get meth
        [_meth, _pos] = [[], []]
        for me in analysis_files.handler_methyl.fetch(gene_info[use_gene].chromosome, gene_start, gene_end):
            [_, pos, _, _, _, _, _, _, _, _, p_meth, _, _, _] = me.split("\t")
            _pos.append(int(pos))
            _meth.append(float(p_meth))
        #methylation = {"pos": _pos, "ratio": _meth}
        # plot results
        plot_hla = HLAPlot()
        plot_hla.output_directory = use_args.work_dir
        plot_hla.plot_hla_result(gene_name=use_gene, gene_info=gene_info[use_gene])
    # Region analysis
    # get region coverage
    [_cov, _pos] = [[], []]
    for cov in analysis_files.handler_cov.fetch(hla_region["chromosome"], hla_region["start"], hla_region["end"]):
        [_, _cov_start, _, _cov_val] = cov.split("\t")
        _pos.append(int(_cov_start))
        _cov.append(float(_cov_val))
    coverage_region = {"pos": np.array(_pos), "cov": np.array(_cov)}
    # get region methylation
    [_meth, _pos] = [[], []]
    for me_pos in range(int(hla_region["start"]), int(hla_region["end"])+1, int(windows_size)):
        [tmp_meth, tmp_pos] = [[], []]
        for me in analysis_files.handler_methyl.fetch(hla_region["chromosome"], me_pos, int(me_pos+windows_size)):
            # [_, pos, _, _, qc, _, _, _, _, tot, p_meth, c_canon, c_meth, c_rm] = me.split("\t")
            [_, pos, _, _, _, _, _, _, _, _, p_meth, _, _, _] = me.split("\t")
            tmp_pos.append(int(pos))
            tmp_meth.append(float(p_meth))
        if len(tmp_pos) > 0:
            _pos.append(tmp_pos[0])
            _meth.append(np.nanmean(np.array(tmp_meth)))
        else:
            _pos.append(me_pos)
            _meth.append(np.nan)
    #methylation = {"pos": _pos, "ratio": _meth}
    hla_plot = HLAPlot()
    hla_plot.output_directory = use_args.work_dir
    hla_plot.plot_hla_region(chr_region="Chr6:29000000-34000000", gene_info=gene_info)
    # close
    analysis_files.close_handlers_tumor_only()

