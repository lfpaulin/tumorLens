import os
import sys
import pysam
import json
import logging
import subprocess
import numpy as np
from scipy.signal import savgol_filter
from matplotlib import pyplot as plot_engine
from matplotlib import gridspec
import matplotlib.ticker as ticker
from config.parameters import LossParameters
from utils.term_colors import TermColors as Col
from utils.logger_debug import setup_log
from utils.html_result import HTMLResult
from utils import parse_vcf
from utils import typing_post
from .gene import GeneResults
from .hla import HLA_GENES
from .hla import HLALocusResults
from pipelines.data_process import BinPathsAnnot


# from .hla import load_hla_annotation
# from .hla import HLAType
# from .analysis import load_results

# Globals
logger = setup_log(__name__, debug=False)
tumor = "tumor"
normal = "normal"
padding = 1000  # 1kb


class AnalysisFiles(object):
    # init
    def __init__(self):
        # files/dir
        self.sv = ""  # SV file with results (vcf)
        self.snf = ""  # SV file with results (snf)
        self.snf_mosaic = ""  # SV file with results (snf)
        self.snv = ""  # SNV file with results (vcf)
        self.cnv = ""  # CNV file with results (vcf)
        self.spc = ""  # CNV file with results (spc)
        self.coverage = ""  # Coverage file with results (bed)
        self.methyl = ""  # Methylation calls (bed)
        self.methyl_hap_res = ""  # methylation analysis by hap summary (txt)
        self.hlatype = ""  # Typing results for HLA genes
        self.bam = ""    # Alignment file
        self.loh_af_threshold = 0.75
        self.name = ""
        # file handlers
        self.handler_sv = None       # pysam.VariantFile
        self.handler_snv = None      # pysam.VariantFile
        self.handler_cnv = None      # pysam.VariantFile
        self.handler_cov = None      # pysam.Tabix
        self.handler_methyl = None   # pysam.Tabix
        self.handler_methyl_hap = None   # open read
        self.handler_hlatype = None  # open Read

    def args_check_somatic(self, user_args, sample_type=""):
        if "tumor" == sample_type:
            try:
                sample_files = json.load(open(user_args.tumor_json))
                logger.info(f'{Col.bold}Using file {user_args.tumor_json} {Col.end}')
            except OSError:
                logger.error(f'{Col.bold}{Col.Fg.red}File {user_args.tumor_json} could not be opened'
                             f'{Col.end}')
                sys.exit(1)
        elif "normal" == sample_type or "control" == sample_type:
            try:
                sample_files = json.load(open(user_args.control_json))
                logger.info(f'{Col.bold}Using file {user_args.control_json} {Col.end}')
            except OSError:
                logger.error(f'{Col.bold}{Col.Fg.red}File {user_args.control_json} could not be opened'
                             f'{Col.end}')
                sys.exit(1)
        else:
            logger.error(f'{Col.bold}{Col.Fg.red}Option {sample_type} unknown{Col.end}')
            sys.exit(1)
        # errors 404 -> not found
        #        503 -> could not open
        self.name = sample_files["name"]
        logger.info(f'Sample {sample_type}: {Col.bold}{Col.Fg.blue}{self.name}{Col.end}')
        # SV       ###########################################################
        self.sv = sample_files["sv"]
        if os.path.exists(self.sv):
            logger.info(f'Using SV:         {self.sv}')
            try:
                self.handler_sv = pysam.VariantFile(self.sv)
            except OSError as os_err:
                logger.error(f'{Col.bold}{Col.Fg.red}Could not open file {self.sv}.\n'
                             f'{os_err}\nID:503-4{Col.end}')
                sys.exit(1)
        else:
            logger.error(f'{Col.bold}{Col.Fg.red}{self.sv} files not found and needed. '
                         f'ID:404-4{Col.end}')
            sys.exit(1)
        # SV+SNF       #######################################################
        self.snf = sample_files["snf"]
        if os.path.exists(self.snf):
            logger.info(f'Using SNF:        {self.snf}')
        else:
            logger.error(f'{Col.bold}{Col.Fg.red}{self.snf} files not found and needed. '
                         f'ID:404-4.2{Col.end}')
            sys.exit(1)
        self.snf_mosaic = sample_files["snf_mosaic"]
        if os.path.exists(self.snf_mosaic):
            logger.info(f'Using SNF mosaic: {self.snf_mosaic}')
        else:
            logger.error(f'{Col.bold}{Col.Fg.red}{self.snf_mosaic} files not found and needed. '
                         f'ID:404-4.2{Col.end}')
            sys.exit(1)
        # SNV      ###########################################################
        self.snv = sample_files["snv"]
        if os.path.exists(self.snv):
            logger.info(f'Using SNV:        {self.snv}')
            try:
                self.handler_snv = pysam.VariantFile(self.snv)
            except OSError as os_err:
                logger.error(f'{Col.bold}{Col.Fg.red}Could not open file {self.snv}.\n{os_err}\n'
                             f'ID:503-5{Col.end}')
                sys.exit(1)
        else:
            logger.error(f'{Col.bold}{Col.Fg.red}{self.snv} file not found and needed. '
                         f'ID:404-5{Col.end}')
            sys.exit(1)
        # CNV      ###########################################################
        self.cnv = sample_files["cnv"]
        if os.path.exists(self.cnv):
            logger.info(f'Using CNV:        {self.cnv}')
            try:
                self.handler_cnv = pysam.VariantFile(self.cnv)
            except OSError as os_err:
                logger.error(f'{Col.bold}{Col.Fg.red}Could not open {self.cnv}.\n{os_err}\n'
                             f'ID:503-6{Col.end}')
                sys.exit(1)
        else:
            logger.error(f'{Col.bold}{Col.Fg.red}{self.cnv} not found and needed. '
                         f'ID:404-6{Col.end}')
            sys.exit(1)
        # CNV+SPC      #######################################################
        self.spc = sample_files["spc"]
        """ not used right now"""
        # if os.path.exists(self.spc):
        #     logger.info(f'Using SPC:        {self.spc}')
        # else:
        #     logger.error(f'{Col.bold}{Col.Fg.red}{self.spc} not found and needed. '
        #                  f'ID:404-6.2{Col.end}')
        #     sys.exit(1)
        # Coverage ###########################################################
        self.coverage = sample_files["coverage"]
        if os.path.exists(self.coverage):
            logger.info(f'Using coverage:   {self.coverage}')
            try:
                self.handler_cov = pysam.TabixFile(self.coverage)
            except OSError as os_err:
                logger.error(f'{Col.bold}{Col.Fg.red}Could not open {self.coverage} coverage.\n'
                             f'{os_err}\nID:503-7{Col.end}')
                sys.exit(1)
        else:
            logger.error(f'{Col.bold}{Col.Fg.red}{self.coverage} not found and needed. '
                         f'ID:404-7{Col.end}')
            sys.exit(1)
        # Methyl   ###########################################################
        self.methyl = sample_files["methyl"]
        if os.path.exists(self.methyl):
            logger.info(f'Using methyl:     {self.methyl} methylation')
            try:
                self.handler_methyl = pysam.TabixFile(self.methyl)
            except OSError as os_err:
                logger.error(f'{Col.bold}{Col.Fg.red}could not open {self.methyl} methylation.\n'
                             f'{os_err}\nID:503-8{Col.end}')
                sys.exit(1)
        else:
            logger.error(f'{Col.bold}{Col.Fg.red}{self.methyl} not found and needed. '
                         f'ID:404-8{Col.end}')
            sys.exit(1)
        # Methyl hap  ########################################################
        self.methyl_hap_res = sample_files["methyl_hap"]
        if os.path.exists(self.methyl_hap_res):
            logger.info(f'Using methyl hap: {self.methyl_hap_res}')
            try:
                self.handler_methyl_hap = open(self.methyl_hap_res, "r") if "tumor" == sample_type else ""
            except OSError as os_err:
                logger.error(f'{Col.bold}{Col.Fg.red}could not open {self.methyl_hap_res} methylation.\n'
                             f'{os_err}\nID:503-9{Col.end}')
                sys.exit(1)
        else:
            logger.error(f'{Col.bold}{Col.Fg.red}{self.methyl_hap_res} not found and needed. '
                         f'ID:404-9{Col.end}')
            # sys.exit(1)
        # HLA typing   #######################################################
        self.hlatype = sample_files["typing"]
        if os.path.exists(self.hlatype):
            logger.info(f'Using HLA typing: {self.hlatype}')
            try:
                self.handler_hlatype = open(self.hlatype, "r")
            except OSError as os_err:
                logger.error(f'{Col.bold}{Col.Fg.red}could not open {self.hlatype}  HLA typing.\n'
                             f'{os_err}\nID:503-10{Col.end}')
                sys.exit(1)
        else:
            logger.error(f'{Col.bold}{Col.Fg.red}{self.hlatype} not found and needed. '
                         f'ID:404-10{Col.end}')
            sys.exit(1)
        # End      ###########################################################

    def close_handlers(self):
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
        # End      ###########################################################


class HLAPlot:
    def __init__(self, name_tumor: str = "", name_normal: str = "", region: bool = False, loh_af: float = 0.75):
        # names
        self.names = {tumor: name_tumor, normal: name_normal}
        # the plot
        self.ncol = 2
        self.nrow = 3 if region else 4
        self.width = 10 if region else 11
        self.height = 7
        self.figure = plot_engine.figure(figsize=(self.width, self.height))
        height_ratios_subplot = [2, 3, 2] if region else [2, 2, 2, 1]
        self.variant_line_width = 3
        self.coverage_point_size = 5
        gs = gridspec.GridSpec(self.nrow, self.ncol, height_ratios=height_ratios_subplot)
        self.cov_plot = {tumor: plot_engine.subplot(gs[0]), normal: plot_engine.subplot(gs[1])}
        self.variants = {tumor: plot_engine.subplot(gs[2]), normal: plot_engine.subplot(gs[3])}
        self.meth_plot = {tumor: plot_engine.subplot(gs[4]), normal: plot_engine.subplot(gs[5])}
        if not region:
            self.gene_body = {tumor: plot_engine.subplot(gs[6]), normal: plot_engine.subplot(gs[7])}
        # colors => see func
        # legends and axis
        self.axis_ylim = {"bottom": 0.0, "top": 0.0}  # to init
        self.file_prefix = "Analysis"
        self.output_directory = "./"
        self.loh_af_threshold = loh_af

    @staticmethod
    def colors(element):
        non_var_colors = {
            "coverage": "#377eb8",       # blue
            "methyl": "#ec7014",         # orange
            "methyl_smooth": "#fe9929",  # orange light
            "gene_body": "#252525",      # dark gray
            "promoter": "#969696",       # light gray
            "gene_spacing": "#111111",   # black
            "gt_vaf": "#a1a1a1",         # gray
            "gt_vaf_alpha": 0.5,
            "used_vaf": "#cb181d"
        }
        variant_color = {
            "loh_loh": "#c985bd",   # pink
            "cnv_del": "#66cc00",   # green
            "cnv_dup": "#ff2d2d",   # red
            "sv_ins": "#6600cc",    # purple
            "sv_del": "#00CC66",    # green
            "sv_inv": "#0066cc",    # blue
            "sv_dup": "#ff3355",    # red
            "sv_bnd": "#404040"     # black
        }
        if element in non_var_colors:
            return non_var_colors[element]
        if element in variant_color:
            return variant_color[element]
        logger.warning(f'{Col.Fg.orange}element {element} nor found, default color is black (#000000)'
                       f'{Col.end}')
        return '#000000'

    def plot_hla_region(self, plot_results=None):
        if plot_results is None:
            logger.error(f"{Col.bold}{Col.Fg.red}'plot_results' parameter is needed{Col.end}")
            sys.exit(1)
        for sample_type in [tumor, normal]:
            chr_region = plot_results[sample_type].region_coordinates
            logger.info(f"plotting: {sample_type} => {chr_region}")
            # reference region
            [region_start, region_end] = plot_results[sample_type].position
            # coverage plot
            self.axis_ylim["top"] = np.max([v for k in [tumor, normal] for v in plot_results[k].cov["cov"]])
            coverage = plot_results[sample_type].cov
            self.cov_plot[sample_type].set_title(self.names[sample_type])
            self.cov_plot[sample_type].axes.set_ylabel(f'Coverage')
            self.cov_plot[sample_type].scatter(np.array(coverage["pos"]), np.array(coverage["cov"]),
                                               color=self.colors("coverage"), s=self.variant_line_width)
            self.cov_plot[sample_type].axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
            self.cov_plot[sample_type].axes.set_xlim(left=region_start, right=region_end)
            # Variants plot
            variants = {"loh": plot_results[sample_type].loh,
                        "cnv": plot_results[sample_type].cnv,
                        "sv": plot_results[sample_type].sv}
            self.axis_ylim["bottom"] = 0.5
            self.axis_ylim["top"] = 3.5
            self.variants[sample_type].axes.set_ylabel(f'Variants')
            variant_y_values = {"sv": 1, "cnv": 2, "loh": 3}
            for each_variant_type in variants.keys():
                var_counter = 0
                for each_variant in variants[each_variant_type]:
                    [start, end] = each_variant.start, each_variant.stop
                    variant_color = self.colors(f'{each_variant_type}_{each_variant.info["SVTYPE"].lower()}')
                    offset = (0 if var_counter % 3 == 0 else (0.1 if var_counter % 3 == 1 else -0.1))
                    var_y_value = variant_y_values[each_variant_type] + (offset if each_variant_type == "sv" else 0)
                    self.variants[sample_type].plot(np.array([start, end]),
                                                    np.array([var_y_value, var_y_value]),
                                                    color=variant_color,
                                                    linewidth=f'{self.variant_line_width}')
                    var_counter += 1
            self.variants[sample_type].axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
            self.variants[sample_type].axes.set_xlim(left=region_start, right=region_end)
            self.variants[sample_type].tick_params(left=False)
            self.variants[sample_type].set_yticks((1, 2, 3))
            self.variants[sample_type].set_yticklabels(("SV", "CNV", "LOH"))
            # meth plot used for LOH
            self.meth_plot[sample_type].axes.set_ylabel(f'% ALT VAF')
            # snv_gt = plot_results[sample_type].snv["gt"]
            # divided by af >= lof_ah and af < lof_ah
            snv_pos, snv_vaf, snv_posu, snv_vafu = [], [], [], []
            for i in range(0, len(plot_results[sample_type].snv["pos"])):
                vaf = plot_results[sample_type].snv["vaf"][i]
                pos = plot_results[sample_type].snv["pos"][i]
                if vaf < self.loh_af_threshold:
                    snv_pos.append(pos)
                    snv_vaf.append(vaf)
                else:
                    snv_posu.append(pos)
                    snv_vafu.append(vaf)
            for each_variant in variants["loh"]:
                [start, end] = each_variant.start, each_variant.stop
                variant_color = self.colors("loh_loh")
                self.meth_plot[sample_type].plot(np.array([start, end]), np.array([0.5, 0.5]), color=variant_color,
                                                 linewidth=f'{self.variant_line_width}')
            self.meth_plot[sample_type].scatter(np.array(snv_pos), np.array(snv_vaf), color=self.colors("gt_vaf"),
                                                s=self.coverage_point_size, alpha=self.colors("gt_vaf_alpha"))
            self.meth_plot[sample_type].scatter(np.array(snv_posu), np.array(snv_vafu), color=self.colors("used_vaf"),
                                                s=self.coverage_point_size, alpha=self.colors("gt_vaf_alpha"))
            self.meth_plot[sample_type].axes.set_ylim(bottom=-0.1, top=1.1)
            self.meth_plot[sample_type].axes.set_xlim(left=region_start, right=region_end)
        # save and close
        region_pretty = "_".join("_".join(chr_region.split(":")).split("-"))
        self.figure.suptitle(f'{self.file_prefix} HLA Locus {chr_region}')
        self.figure.savefig(f'{self.output_directory}/plot-HLA-{self.file_prefix}-{region_pretty}.png', dpi=300)
        self.figure.clf()

    def plot_gene_result(self, res_data: dict, gene_name: str, gene_info: GeneResults, gene_padding=2e3, promoter=1e3):
        if gene_info is None:
            logger.error(f"{Col.bold}{Col.Fg.red}'gene_info' parameters are needed{Col.end}")
            sys.exit(1)
        logger.info(f'plotting results for {gene_name}')
        # gene info
        for sample_type in [tumor, normal]:
            [gene_start, gene_end] = gene_info.position
            gene_start -= gene_padding
            gene_end += gene_padding
            coverage = res_data[sample_type].cov
            methylation = res_data[sample_type].meth
            # coverage plot
            self.axis_ylim["top"] = np.max(coverage["cov"])*1.05
            self.axis_ylim["top"] = np.max([v for k in [tumor, normal] for v in res_data[k].cov["cov"]])
            self.axis_ylim["bottom"] = 0
            self.cov_plot[sample_type].set_title(self.names[sample_type])
            self.cov_plot[sample_type].axes.set_ylabel(f'Coverage')
            self.cov_plot[sample_type].scatter(np.array(coverage["pos"]), np.array(coverage["cov"]),
                                               color=self.colors("coverage"), marker="o", s=4)
            self.cov_plot[sample_type].axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
            self.cov_plot[sample_type].axes.set_xlim(left=gene_start, right=gene_end)
            tick_sep = np.round((gene_end - gene_start) / 3) - (np.round((gene_end - gene_start) / 3) % 100)
            self.cov_plot[sample_type].xaxis.set_major_locator(ticker.MultipleLocator(base=tick_sep))
            # meth plot
            self.axis_ylim["bottom"] = -5
            self.axis_ylim["top"] = 105
            self.meth_plot[sample_type].axes.set_ylabel(f'')
            # methylation_status = np.array([100 if m > 0.65 else 0 if m < 0.35 else m for m in methylation["mod"]])
            self.meth_plot[sample_type].scatter(np.array(methylation["pos"]), methylation["mod"],  # methylation_status,
                                                color=self.colors("methyl"), marker="o", s=4)
            try:
                meth_mod_yhat = savgol_filter(methylation["mod"], 30, 2)
                self.meth_plot[sample_type].plot(methylation["pos"], meth_mod_yhat, linewidth='2',
                                                 color=self.colors("methyl_smooth"), linestyle="solid")
            except ValueError:
                logger.warning(f"Unable to perform 'savgol_filter'")
            self.meth_plot[sample_type].axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
            self.meth_plot[sample_type].axes.set_xlim(left=gene_start, right=gene_end)
            self.meth_plot[sample_type].xaxis.set_major_locator(ticker.MultipleLocator(base=tick_sep))
            self.meth_plot[sample_type].tick_params(left=False)
            self.meth_plot[sample_type].set_yticks((0, 25, 50, 75, 100))
            # self.meth_plot[sample_type].set_yticklabels(("Non-methylated", "", "Unresolved", "", "Methylated"))
            # self.meth_plot[sample_type].tick_params(axis='y', labelrotation=45)
            # Drawing the genetic element body
            # 1kb upstream as promoter and arrow direction
            if "+" == gene_info.strand:
                promoter_start = gene_info.start - promoter
                promoter_end = gene_info.start
                arrow_start = gene_info.end
                arrow_end = gene_padding / 8
            else:
                promoter_start = gene_info.end
                promoter_end = gene_info.end + promoter
                arrow_start = gene_info.start
                arrow_end = -gene_padding / 8
            self.gene_body[sample_type].axes.set_xlabel(f'Chromosome position')
            self.gene_body[sample_type].axes.set_xlim(left=gene_start, right=gene_end)
            self.gene_body[sample_type].plot(np.array([promoter_start, promoter_end]), np.array([0, 0]),
                                             linewidth='3', color=self.colors("promoter"))
            # Arrow for strand  -> or <-
            self.gene_body[sample_type].arrow(x=arrow_start, y=0, dx=arrow_end, dy=0, head_width=0.5,
                                              head_length=(gene_padding / 8) * 0.95,
                                              length_includes_head=True, color=self.colors("gene_body"))
            # gene body
            [start, end] = gene_info.position
            self.gene_body[sample_type].plot(np.array([start, end]), np.array([0, 0]), linewidth='5',
                                             color=self.colors("gene_body"))
            self.gene_body[sample_type].xaxis.set_major_locator(ticker.MultipleLocator(base=tick_sep))
            self.gene_body[sample_type].axes.set_ylabel(f'HLA genes')
            self.gene_body[sample_type].tick_params(left=False)
            self.gene_body[sample_type].set_yticks((-1, 0, 1))
            self.gene_body[sample_type].set_yticklabels(("", "", ""))
            # variants per gene
            variants = {"loh": res_data[sample_type].loh,
                        "cnv": res_data[sample_type].cnv,
                        "sv": res_data[sample_type].sv}
            self.axis_ylim["bottom"] = 0
            self.axis_ylim["top"] = 4
            self.variants[sample_type].axes.set_ylabel(f'Variants')
            variant_y_values = {"sv": 1, "cnv": 2, "loh": 3}
            vars_unpack = {"sv": "", "cnv": "", "loh": ""}
            for each_variant_type in variants.keys():
                n_variants = len(variants[each_variant_type])
                if n_variants > 0:
                    vars_unpack_svt = [var.info["SVTYPE"] for var in variants[each_variant_type]]
                    vars_unpack_svl = [var.info["SVLEN"] for var in variants[each_variant_type]]
                    vars_details = [f'{svt}:{svl}' if abs(svl) >= 50 else ""
                                    for svt, svl in zip(vars_unpack_svt, vars_unpack_svl)]
                    vars_unpack[each_variant_type] = f'{n_variants}|{vars_details}'
                var_counter = 0
                for each_variant in variants[each_variant_type]:
                    [start, end] = each_variant.start, each_variant.stop
                    variant_color = self.colors(f'{each_variant_type}_{each_variant.info["SVTYPE"].lower()}')
                    var_y_value = variant_y_values[each_variant_type] + (0 if var_counter % 3 == 0 else
                                                                         (0.1 if var_counter % 3 == 1 else -0.1))
                    self.variants[sample_type].plot(np.array([start, end]), np.array([var_y_value, var_y_value]),
                                                    linewidth='2', color=variant_color)
                    var_counter += 1
            self.variants[sample_type].axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
            self.variants[sample_type].axes.set_xlim(left=gene_start, right=gene_end)
            self.variants[sample_type].xaxis.set_major_locator(ticker.MultipleLocator(base=tick_sep))
            self.variants[sample_type].tick_params(left=False)
            self.variants[sample_type].set_yticks((0, 1, 2, 3, 4))
            self.variants[sample_type].set_yticklabels(("", "SV", "CNV", "LOH", ""))
            logger.info(f"Result of gene: {gene_name}|{sample_type}, Gene body: {gene_info.position} | "
                        f"{vars_unpack}")
        # save and close
        self.figure.suptitle(f'{self.file_prefix} gene: {gene_name}')
        gene_is = "HLA" if "HLA" in gene_name else "GENE"
        self.figure.savefig(f'{self.output_directory}/plot-{gene_is}-{self.file_prefix}-{gene_name}.png', dpi=300)
        self.figure.clf()


def somatic_sv(tumor_snf, normal_snf, annot_bed, reference, tumor_snf_mosaic):
    logger.info("SV merge/compare somatic")
    annot_bed_handler = pysam.TabixFile(annot_bed)
    merge_tsv = open("merge_sv.tsv", "w")
    merge_tsv.write(f'{tumor_snf}\ttumor\n{normal_snf}\tnormal')
    merge_tsv.close()
    merge_tsv = open("merge_sv_mosaic.tsv", "w")
    merge_tsv.write(f'{tumor_snf_mosaic}\ttumor\n{normal_snf}\tnormal')
    merge_tsv.close()

    # SV merge germline
    cmd_sniffles = (f'sniffles --input merge_sv.tsv --vcf merge_sv.vcf.gz --allow-overwrite --threads 4 '
                    f'--reference {reference}')
    run_cmd = subprocess.run(cmd_sniffles, shell=True, capture_output=True, text=True)
    if run_cmd.stderr != "":
        logger.warning(f'{Col.Fg.orange}{run_cmd.stderr}{Col.end}')
        logger.error(f'{Col.bold}{Col.Fg.red}Cannot run {cmd_sniffles}, '
                     f'exiting with error{Col.end}')
        sys.exit(1)
    # SV merge mosaic
    cmd_sniffles = (f'sniffles --input merge_sv_mosaic.tsv --vcf merge_sv_mosaic.vcf.gz --allow-overwrite --threads 4 '
                    f'--reference {reference}')
    run_cmd = subprocess.run(cmd_sniffles, shell=True, capture_output=True, text=True)
    if run_cmd.stderr != "":
        logger.warning(f'{Col.Fg.orange}{run_cmd.stderr}{Col.end}')
        logger.error(f'{Col.bold}{Col.Fg.red}Cannot run {cmd_sniffles}, '
                     f'exiting with error{Col.end}')
        sys.exit(1)
    # populate genes that are affected by somatic SVs
    sv_affecting_genes_results = {}
    for gene_entry in annot_bed_handler.fetch():
        contig, pos, end, gene, _, _ = gene_entry.split("\t")
        reg_gene = f'{contig}:{pos}-{end}'
        sv_affecting_genes_results[gene] = {"region": reg_gene, "sv": [], "init": False}
    sv_affecting_genes_results, sv_somatic_list = parse_vcf.sv_somatic(pysam.VariantFile("merge_sv.vcf.gz"),
                                                                       sv_affecting_genes_results)
    return sv_affecting_genes_results


def somatic_cnv(tumor_vcf, annot_bed):
    logger.info("CNV somatic/compare")
    somatic_cnv_tumor = pysam.VariantFile(tumor_vcf)
    annot_bed_handler = pysam.TabixFile(annot_bed)
    cnv_loh_affecting_genes_results = {}
    for gene_entry in annot_bed_handler.fetch():
        contig, pos, end, gene, _, _ = gene_entry.split("\t")
        reg_gene = f'{contig}:{pos}-{end}'
        cnv_loh_affecting_genes_results[gene] = {"region": reg_gene, "cnv": [], "init": False}
    cnv_loh_affecting_genes_results, cnv_somatic_list = parse_vcf.cnv_somatic(somatic_cnv_tumor,
                                                                              cnv_loh_affecting_genes_results)
    return cnv_loh_affecting_genes_results


def merge_type(tumor_type, normal_type, hla_rareness_file):
    logger.info("HLA typing comparison")
    tumor_results = open(tumor_type)
    normal_results = open(normal_type)
    hla_rareness = json.load(open(hla_rareness_file))
    rare_hla = ["WD", "X"]
    # compare by pairs, use counter, if even save hla gene name if odd compare
    allele1 = {tumor: "", normal: ""}
    allele2 = {tumor: "", normal: ""}
    counter = 0
    gene = ""
    hla_type_compare = {}
    for tumor_line, normal_line in zip(tumor_results, normal_results):
        if tumor_line.startswith("#") or tumor_line.startswith("G"):
            pass
        else:
            tumor_res_list = tumor_line.split("\t")
            normal_res_list = normal_line.split("\t")
            if counter % 2 == 0:
                allele1[tumor] = ":".join(tumor_res_list[1].split(":")[:2])
                allele1[normal] = ":".join(normal_res_list[1].split(":")[:2])
                gene = "-".join(tumor_res_list[0].split("_")[:2])
            else:
                allele2[tumor] = ":".join(tumor_res_list[1].split(":")[:2])
                allele2[normal] = ":".join(normal_res_list[1].split(":")[:2])
                # check empty #####
                if allele1[tumor] == "":
                    allele1[tumor] = allele2[tumor]
                elif allele2[tumor] == "":
                    allele2[tumor] = allele1[tumor]
                else:
                    pass
                if allele1[normal] == "":
                    allele1[normal] = allele2[normal]
                elif allele2[normal] == "":
                    allele2[normal] = allele1[normal]
                else:
                    pass
                # compare #####
                # check rarity
                tumor1 = hla_rareness[allele1[tumor]] if allele1[tumor] in hla_rareness else None
                tumor2 = hla_rareness[allele2[tumor]] if allele2[tumor] in hla_rareness else None
                normal1 = hla_rareness[allele1[normal]] if allele1[normal] in hla_rareness else None
                normal2 = hla_rareness[allele2[normal]] if allele2[normal] in hla_rareness else None
                if tumor1 in rare_hla or tumor2 in rare_hla:
                    if tumor1 in rare_hla and tumor2 in rare_hla:
                        allele1[tumor] = f'{allele1[tumor]}(rare: {tumor1})'
                        allele2[tumor] = f'{allele2[tumor]}(rare: {tumor2})'
                    elif tumor1 in rare_hla and tumor2 not in rare_hla:
                        allele1[tumor] = f'{allele1[tumor]}(rare: {tumor1})'
                    elif tumor1 not in rare_hla and tumor2 in rare_hla:
                        allele2[tumor] = f'{allele2[tumor]}(rare: {tumor2})'
                    else:
                        pass
                if normal1 in rare_hla or normal2 in rare_hla:
                    if normal1 in rare_hla and normal2 in rare_hla:
                        allele1[normal] = f'{allele1[normal]}(rare: {normal1})'
                        allele2[normal] = f'{allele2[normal]}(rare: {normal2})'
                    elif normal1 in rare_hla and normal2 not in rare_hla:
                        allele1[normal] = f'{allele1[normal]}(rare: {normal1})'
                    elif normal1 not in rare_hla and normal2 in rare_hla:
                        allele2[normal] = f'{allele2[normal]}(rare: {normal2})'
                    else:
                        pass
                # output result
                if allele1[tumor] == allele2[tumor] and allele1[normal] == allele2[normal]:
                    hla_type_compare[gene] = {"class": "LOH_T_N", "tumor": allele1[tumor], "normal": allele1[normal]}
                elif allele1[tumor] == allele2[tumor] and allele1[normal] != allele2[normal]:
                    hla_type_compare[gene] = {"class": "LOH_T", "tumor": allele1[tumor],
                                              "normal": f'{allele1[normal]}/{allele2[normal]}'}
                elif allele1[tumor] != allele2[tumor] and allele1[normal] == allele2[normal]:
                    hla_type_compare[gene] = {"class": "LOH_N", "tumor": f'{allele1[tumor]}/{allele2[tumor]}',
                                              "normal": allele1[normal]}
                else:
                    control_hla = [allele1[normal], allele2[normal]]
                    if allele1[tumor] not in control_hla or allele2[tumor] not in control_hla:
                        # try 2 digits before bang
                        hla1t_2d, hla2t_2d = allele1[tumor].split(":")[0],  allele2[tumor].split(":")[0]
                        hlac_2d = [allele1[normal].split(":")[0], allele2[normal].split(":")[0]]
                        if hla1t_2d not in hlac_2d or hla2t_2d not in hlac_2d:
                            type_bang = "[!]"
                        else:
                            type_bang = "[~]"
                    else:
                        type_bang = ""
                    hla_type_compare[gene] = {"class": "HET", "tumor": f'{type_bang}{allele1[tumor]}/{allele2[tumor]}',
                                              "normal": f'{allele1[normal]}/{allele2[normal]}'}
            counter += 1
    return hla_type_compare


def coverage_finer(tumor_type: str, control_type: str, our_dir: str) -> tuple[pysam.TabixFile, pysam.TabixFile]:
    logger.info("Running coverage_finer")
    tumor_dir = "/".join(tumor_type.split("/")[:-1])
    control_dir = "/".join(control_type.split("/")[:-1])
    tumor_hla_region_bam = f'{tumor_dir}/_hla_region.bam'
    control_hla_region_bam = f'{control_dir}/_hla_region.bam'
    cmd_new_cov_tumor = (f'cd {our_dir}; mosdepth --by 10 --threads 2 --no-per-base --mapq 10 --chrom chr6 '
                         f'tumor_cov  {tumor_hla_region_bam}; '
                         f'tabix --preset bed {our_dir}/tumor_cov.regions.bed.gz')
    logger.info(f'Running tumor')
    logger.debug(f'Running: {cmd_new_cov_tumor}')
    run_cmd_tumor = subprocess.run(cmd_new_cov_tumor, shell=True, capture_output=True, text=True)
    if run_cmd_tumor.stderr != "":
        logger.warning(f'{Col.Fg.orange}{run_cmd_tumor.stderr}{Col.end}')
        if "rror" in run_cmd_tumor.stderr:
            logger.error(f'{Col.bold}{Col.Fg.red}Cannot run {cmd_new_cov_tumor}, '
                         f'exiting with error{Col.end}')
            sys.exit(1)
    cmd_new_cov_control = (f'cd {our_dir}; mosdepth --by 10 --threads 2 --no-per-base --mapq 10 --chrom chr6 '
                           f'control_cov  {control_hla_region_bam}; '
                           f'tabix --preset bed {our_dir}/control_cov.regions.bed.gz')
    logger.info(f'Running control')
    logger.debug(f'Running: {cmd_new_cov_control}')
    run_cmd_control = subprocess.run(cmd_new_cov_control, shell=True, capture_output=True, text=True)
    if run_cmd_control.stderr != "":
        logger.warning(f'{Col.Fg.orange}{run_cmd_control.stderr}{Col.end}')
        if "rror" in run_cmd_control.stderr:
            logger.error(f'{Col.bold}{Col.Fg.red}Cannot run {cmd_new_cov_control}, '
                         f'exiting with error{Col.end}')
            sys.exit(1)
    # TODO: get the new file names of new bed files
    return (pysam.TabixFile(f"{our_dir}/tumor_cov.regions.bed.gz"),
            pysam.TabixFile(f"{our_dir}/control_cov.regions.bed.gz"))


def post_type_log(tumor_hla: dict, control_hla: dict):
    logger.info(f'Gene\tTumor\tControl')
    for hla_gene in tumor_hla.keys():
        tumor1, t_allele1, tumor2, t_allele2, _ = tumor_hla[hla_gene]
        control1, c_allele1, control2, c_allele2, _ = control_hla[hla_gene]
        tumor_total = t_allele1 + t_allele2
        control_total = c_allele1 + c_allele2
        # cases: T1/T2 v C1/C2
        # * T1 == C1, T2 == C2  <= only here we have HOM v HOM
        # *         , T2 != C2
        # * T1 == C2, T2 == C1
        # *         , T2 != C1
        # * T1 != C1, T2 == C2
        # T1 != C2, T2 == C1
        # T1 != C1, T2 != C2
        if tumor1 == control1:
            logger.info(f'{tumor1}\t{(t_allele1/tumor_total):0.3f}\t{(c_allele1/control_total):0.3f}')
            if tumor2 == control2:
                if tumor1 == tumor2:
                    pass
                else:
                    logger.info(f'{tumor2}\t{(t_allele2/tumor_total):0.3f}\t{(c_allele2/control_total):0.3f}')
            else:
                if control1 == control2 and tumor1 != tumor2:
                    logger.info(f'{tumor2}\t{(t_allele2/tumor_total):0.3f}\t0')
                elif tumor1 == tumor2 and control1 != control2:
                    logger.info(f'{control2}\t0\t{(c_allele2/control_total):0.3f}')
                else:
                    logger.info(f'{control2}\t0\t{(c_allele2/control_total):0.3f}')
                    logger.info(f'{tumor2}\t{(t_allele2/tumor_total):0.3f}\t0')
        elif tumor1 == control2:
            logger.info(f'{tumor1}\t{(t_allele1/tumor_total):0.3f}\t{(c_allele2/control_total):0.3f}')
            if tumor2 == control1:
                logger.info(f'{tumor2}\t{(t_allele2/tumor_total):0.3f}\t{(c_allele1/control_total):0.3f}')
            else:
                if tumor1 == tumor2:
                    logger.info(f'{control1}\t0\t{(c_allele1/control_total):0.3f}')
                else:
                    logger.info(f'{control1}\t0\t{(c_allele1/control_total):0.3f}')
                    logger.info(f'{tumor2}\t{(t_allele2/tumor_total):0.3f}\t0')
        elif tumor1 != control1 and tumor2 == control2:
            logger.info(f'{tumor2}\t{(t_allele2/tumor_total):0.3f}\t{(c_allele2/control_total):0.3f}')
            logger.info(f'{control1}\t0\t{(c_allele1/control_total):0.3f}')
            if tumor1 != tumor2:
                logger.info(f'{tumor1}\t{(t_allele1/tumor_total):0.3f}\t0')
        elif tumor1 != control2 and tumor2 == control1:
            logger.info(f'{tumor2}\t{(t_allele2/tumor_total):0.3f}\t{(c_allele1/control_total):0.3f}')
            logger.info(f'{control2}\t0\t{(c_allele2/control_total):0.3f}')
            if tumor1 != tumor2:
                logger.info(f'{tumor1}\t{(t_allele1/tumor_total):0.3f}\t0')
        else:
            logger.info(f'{control1}\t0\t{(c_allele1/control_total):0.3f}')
            logger.info(f'{control2}\t0\t{(c_allele2/control_total):0.3f}')
            logger.info(f'{tumor1}\t{(t_allele1/tumor_total):0.3f}\t0')
            logger.info(f'{tumor2}\t{(t_allele2/tumor_total):0.3f}\t0')


def post_type_log2(tumor_hla: dict, control_hla: dict) -> list:
    save_results = []
    header = f'Gene\tAllele1\tAllele2\tAllele3\tAllele4'
    logger.info(header)
    save_results.append(f'Type\t{header}')
    logger.debug(tumor_hla)
    logger.debug(control_hla)
    for hla_gene in tumor_hla.keys():
        tumor1, t_allele1, tumor2, t_allele2, _ = tumor_hla[hla_gene]
        control1, c_allele1, control2, c_allele2, _ = control_hla[hla_gene]
        tumor_total = t_allele1 + t_allele2
        control_total = c_allele1 + c_allele2
        tumor_line = []
        control_line = []
        # cases: T1/T2 v C1/C2
        # * T1 == C1, T2 == C2  <= only here we have HOM v HOM
        # *         , T2 != C2
        # * T1 == C2, T2 == C1
        # *         , T2 != C1
        # * T1 != C1, T2 == C2
        # T1 != C2, T2 == C1
        # T1 != C1, T2 != C2
        if tumor1 == control1:
            # Allele 1/2
            tumor_line = [tumor1, f'{(t_allele1/tumor_total):0.3f}'] if tumor_total > 0 else [tumor1, "0.0"]
            control_line = [control1, f'{(c_allele1/control_total):0.3f}'] if control_total > 0 else [control1, "0.0"]
            if tumor2 == control2:
                if tumor1 == tumor2:
                    t1, t2 = tumor_line
                    logger.info(f'T|{t1},-\t{t2}\t0\t0\t0')
                    save_results.append(f'Tumor\t{t1},-\t{t2}\t0\t0\t0')
                    c1, c2 = control_line
                    logger.info(f'C|{c1},-\t{c2}\t0\t0\t0')
                    save_results.append(f'Control\t{c1},-\t{c2}\t0\t0\t0')
                else:
                    tumor_line += [tumor2, f'{(t_allele2/tumor_total):0.3f}'] if tumor_total > 0 else [tumor2, "0.0"]
                    control_line += [control2, f'{(c_allele2/control_total):0.3f}'] if control_total > 0 else [control2, "0.0"]
                    t1, t2, t3, t4 = tumor_line
                    logger.info(f'T|{t1},{t3}\t{t2}\t{t4}\t0\t0')
                    save_results.append(f'Tumor\t{t1},{t3}\t{t2}\t{t4}\t0\t0')
                    c1, c2, c3, c4 = control_line
                    logger.info(f'C|{c1},{c3}\t{c2}\t{c4}\t0\t0')
                    save_results.append(f'Control\t{c1},{c3}\t{c2}\t{c4}\t0\t0')
            else:
                if tumor1 == tumor2 and control1 != control2:
                    t1, t2 = tumor_line
                    logger.info(f'T|{t1},-\t{t2}\t0\t0\t0')
                    save_results.append(f'Tumor\t{t1},-\t{t2}\t0\t0\t0')
                    control_line += [control2, f'{(c_allele2/control_total):0.3f}'] if control_total > 0 else [control2, "0.0"]
                    c1, c2, c3, c4 = control_line
                    logger.info(f'C|{c1},{c3}\t{c2}\t{c4}\t0\t0')
                    save_results.append(f'Control\t{c1},{c3}\t{c2}\t{c4}\t0\t0')
                # Allele 1/2/3
                elif control1 == control2 and tumor1 != tumor2:
                    tumor_line += [tumor2, f'{(t_allele2/tumor_total):0.3f}'] if tumor_total > 0 else [tumor2, "0.0"]
                    t1, t2, t3, t4 = tumor_line
                    logger.info(f'T|{t1},{t3}\t{t2}\t0\t{t4}\t0')
                    save_results.append(f'Tumor\t{t1},{t3}\t{t2}\t0\t{t4}\t0')
                    c1, c2 = control_line
                    logger.info(f'C|{c1},-\t{c2}\t0\t0\t0')
                    save_results.append(f'Control\t{c1},-\t{c2}\t0\t0\t0')
                else:
                    tumor_line += [tumor2, f'{(t_allele2/tumor_total):0.3f}'] if tumor_total > 0 else [tumor2, "0.0"]
                    control_line += [control2, f'{(c_allele2/control_total):0.3f}'] if control_total > 0 else [control2, "0.0"]
                    t1, t2, t3, t4 = tumor_line
                    logger.info(f'T|{t1},{t3}\t{t2}\t0\t{t4}\t0')
                    save_results.append(f'Tumor\t{t1},{t3}\t{t2}\t0\t{t4}\t0')
                    c1, c2, c3, c4 = control_line
                    logger.info(f'C|{c1},{c3}\t{c2}\t{c4}\t0\t0')
                    save_results.append(f'Control\t{c1},{c3}\t{c2}\t{c4}\t0\t0')
        elif tumor1 == control2:
            tumor_line = [tumor1, f'{(t_allele1/tumor_total):0.3f}'] if tumor_total > 0 else [tumor1, "0.0"]
            control_line = [control2, f'{(c_allele2/control_total):0.3f}'] if control_total > 0 else [control2, "0.0"]
            if tumor2 == control1:
                tumor_line += [tumor2, f'{(t_allele2/tumor_total):0.3f}'] if tumor_total > 0 else [tumor2, "0.0"]
                control_line += [control1, f'{(c_allele1/control_total):0.3f}'] if control_total > 0 else [control1, "0.0"]
                t1, t2, t3, t4 = tumor_line
                logger.info(f'T|{t1},{t3}\t{t2}\t{t4}\t0\t0')
                save_results.append(f'Tumor\t{t1},{t3}\t{t2}\t{t4}\t0\t0')
                c1, c2, c3, c4 = control_line
                logger.info(f'C|{c1},{c3}\t{c2}\t{c4}\t0\t0')
                save_results.append(f'Control\t{c1},{c3}\t{c2}\t{c4}\t0\t0')
            else:
                if tumor1 == tumor2:
                    t1, t2 = tumor_line
                    logger.info(f'T|{t1},-\t{t2}\t0\t0\t0')
                    save_results.append(f'Tumor\t{t1},-\t{t2}\t0\t0\t0')
                    control_line += [control1, f'{(c_allele1/control_total):0.3f}'] if control_total > 0 else [control1, "0.0"]
                    c1, c2, c3, c4 = control_line
                    logger.info(f'C|{c1},{c3}\t{c2}\t{c4}\t0\t0')
                    save_results.append(f'Control\t{c1},{c3}\t{c2}\t{c4}\t0\t0')
                else:
                    tumor_line += [tumor2, f'{(t_allele2/tumor_total):0.3f}']
                    control_line += [control1, f'{(c_allele1/control_total):0.3f}'] if control_total > 0 else [control1, "0.0"]
                    t1, t2, t3, t4 = tumor_line
                    logger.info(f'T|{t1},{t3}\t{t2}\t0\t{t4}\t0')
                    save_results.append(f'Tumor\t{t1},{t3}\t{t2}\t0\t{t4}\t0')
                    c1, c2, c3, c4 = control_line
                    logger.info(f'C|{c1},{c3}\t{c2}\t{c4}\t0\t0')
                    save_results.append(f'Control\t{c1},{c3}\t{c2}\t{c4}\t0\t0')
        elif tumor1 != control1 and tumor2 == control2:
            c1, c2, c3, c4 = (control1, f'{(c_allele1/control_total):0.3f}', control2,
                              f'{(c_allele2/control_total):0.3f}') if control_total > 0 else \
                             (control1, f'0.0', control2, f'0.0')
            t1, t2, t3, t4 = (tumor1, f'{(t_allele1/tumor_total):0.3f}', tumor2,
                              f'{(t_allele2 / tumor_total):0.3f}') if tumor_total > 0 else \
                             (tumor1, f'0.0', tumor2, f'0.0')
            if tumor1 != tumor2:
                logger.info(f'T|{t1},{t3}\t{t2}\t0\t{t4}\t0')
                save_results.append(f'Tumor\t{t1},{t3}\t{t2}\t0\t{t4}\t0')
            else:
                logger.info(f'T|{t1},-\t{t2}\t0\t0\t0')
                save_results.append(f'Tumor\t{t1},-\t{t2}\t0\t0\t0')
            logger.info(f'C|{c1},{c3}\t{c2}\t{c4}\t0\t0')
            save_results.append(f'Control\t{c1},{c3}\t{c2}\t{c4}\t0\t0')
        elif tumor1 != control2 and tumor2 == control1:
            c1, c2, c3, c4 = (control1, f'{(c_allele1/control_total):0.3f}', control2,
                              f'{(c_allele2/control_total):0.3f}') if control_total > 0 else \
                             (control1, f'0.0', control2, f'0.0')
            t1, t2, t3, t4 = (tumor2, f'{(t_allele2/tumor_total):0.3f}', tumor1,
                              f'{(t_allele1/tumor_total):0.3f}') if tumor_total > 0 else \
                             (tumor2, f'0.0', tumor1, f'0.0')
            if tumor1 != tumor2:
                logger.info(f'T|{t1},{t3}\t{t2}\t0\t{t4}\t0')
                save_results.append(f'Tumor\t{t1},{t3}\t{t2}\t0\t{t4}\t0')
            else:
                logger.info(f'T|{t1},-\t0\t{t2}\t0\t0')
                save_results.append(f'Tumor\t{t1},-\t0\t{t2}\t0\t0')
            logger.info(f'C|{c1},{c3}\t{c2}\t{c4}\t0\t0')
            save_results.append(f'Control\t{c1},{c3}\t{c2}\t{c4}\t0\t0')
        else:
            c1, c2, c3, c4 = (control1, f'{(c_allele1 / control_total):0.3f}', control2,
                              f'{(c_allele2 / control_total):0.3f}') if control_total > 0 else \
                             (control1, f'0.0', control2, f'0.0')
            t1, t2, t3, t4 = (tumor2, f'{(t_allele2 / tumor_total):0.3f}', tumor1,
                              f'{(t_allele1 / tumor_total):0.3f}') if tumor_total > 0 else \
                             (tumor2, f'0.0', tumor1, f'0.0')
            if tumor1 == tumor2:
                logger.info(f'T|{t1},-\t0\t0\t{t2}\t0')
                save_results.append(f'Tumor\t{t1},-\t0\t0\t{t2}\t0')
            else:
                logger.info(f'T|{t1},{t3}\t0\t0\t{t2}\t{t4}')
                save_results.append(f'Tumor\t{t1},{t3}\t0\t0\t{t2}\t{t4}')
            logger.info(f'C|{c1},{c3}\t{c2}\t{c4}\t0\t0')
            save_results.append(f'Control\t{c1},{c3}\t{c2}\t{c4}\t0\t0')
    return save_results


def post_type(tumor_type: str, control_type: str, annot_dir: str, ref_fasta: str,
              out_dir: str) -> tuple[dict, dict, dict]:
    min_aln_len = 1000
    max_edit_dist = 0.05  # 5%
    allele_ratio_min = 3
    logger.debug(f'{tumor_type} | {control_type}')
    tumor_dir = f'{out_dir}/typing/tumor_post'
    control_dir = f'{out_dir}/typing/control_post'
    if not os.path.exists(f'{out_dir}/typing/'):
        os.mkdir(f'{out_dir}/typing/')
    tumor_type_dir = "/".join(tumor_type.split("/")[:-1])
    control_type_dir = "/".join(control_type.split("/")[:-1])
    tumor_post_args = typing_post.TypeArgs("", tumor_dir, annot_dir, min_aln_len, max_edit_dist,
                                           False, tumor_type_dir)
    control_post_args = typing_post.TypeArgs("", control_dir, annot_dir, min_aln_len, max_edit_dist,
                                             False, control_type_dir)
    control_hla_postprocess, control_reads_alleles = typing_post.main(control_post_args)
    tumor_hla_postprocess, tumor_reads_alleles = typing_post.main(tumor_post_args)
    logger.info('# #########################')
    hla_count_allele = post_type_log2(tumor_hla_postprocess, control_hla_postprocess)
    logger.info(f'HLA allele count written to {out_dir}/hla_count_allele.tsv')
    hla_count_allele_output = open(f'{out_dir}/hla_count_allele.tsv', "w")
    hla_count_allele_output.write("\n".join(hla_count_allele))
    hla_count_allele_output.close()
    hla_postprocess = {}
    for hla_gene in tumor_hla_postprocess.keys():
        hla_out = f"HLA-{hla_gene}"
        hla1_t, allele1_tc, hla2_t, allele2_tc, class_t = tumor_hla_postprocess[hla_gene]
        hla1_c, allele1_cc, hla2_c, allele2_cc, class_c = control_hla_postprocess[hla_gene]
        logger.debug(f'{tumor_hla_postprocess[hla_gene]}|{control_hla_postprocess[hla_gene]}')
        # het / alt
        if class_t == class_c:
            if class_t == "HOM":
                hla_postprocess[hla_out] = ["LOH:2", f'{hla1_c}', f'{hla1_t}']  # both loh
                pass  # not much to say
            if class_t == "HET":
                if allele1_tc > allele2_tc:
                    tum_ratio = allele1_tc/allele2_tc if allele2_tc > 0 else 1
                    tum_res = f'{hla1_t}/{hla2_t}'
                else:
                    tum_ratio = allele2_tc/allele1_tc if allele1_tc > 0 else 1
                    tum_res = f'{hla2_t}/{hla1_t}'
                if allele1_cc > allele2_cc:
                    ctr_ratio = allele1_cc/allele2_cc if allele2_cc > 0 else 1
                    ctr_res = f'{hla1_c}/{hla2_c}'
                else:
                    ctr_ratio = allele2_cc/allele1_cc if allele1_cc > 0 else 1
                    ctr_res = f'{hla2_c}/{hla1_c}'
                if tum_ratio > allele_ratio_min and ctr_ratio > allele_ratio_min:
                    hla_postprocess[hla_out] = ["loh:2", f'{tum_res}({tum_ratio:.3f})', f'{ctr_res} ({ctr_ratio:.3f})']
                    pass  # both loh
                elif tum_ratio > allele_ratio_min > ctr_ratio:
                    hla_postprocess[hla_out] = ["loh:1", f'{tum_res} ({tum_ratio:.3f})', f'{ctr_res} ({ctr_ratio:.3f})']
                    pass  # important
                elif tum_ratio < allele_ratio_min < ctr_ratio:
                    hla_postprocess[hla_out] = ["loh:0", f'{tum_res} ({tum_ratio:.3f})', f'{ctr_res} ({ctr_ratio:.3f})']
                    pass  # weird
                else:
                    hla_postprocess[hla_out] = ["HET:2", f'{tum_res} ({tum_ratio:.3f})', f'{ctr_res} ({ctr_ratio:.3f})']
        else:
            if class_t == "HOM":
                if allele1_cc > allele2_cc:
                    ctr_ratio = allele1_cc / allele2_cc if allele2_cc > 0 else 1
                    ctr_res = f'{hla1_c}/{hla2_c}'
                else:
                    ctr_ratio = allele2_cc / allele1_cc if allele1_cc > 0 else 1
                    ctr_res = f'{hla2_c}/{hla1_c}'
                hla_postprocess[hla_out] = ["LOH:1", f'{hla1_t}', f'{ctr_res} ({ctr_ratio:.3f})']
                pass  # important
            if class_t == "HET":
                if allele1_tc > allele2_tc:
                    tum_ratio = allele1_tc / allele2_tc if allele2_tc > 0 else 1
                    tum_res = f'{hla1_t}/{hla2_t}'
                else:
                    tum_ratio = allele2_tc / allele1_tc if allele1_tc > 0 else 1
                    tum_res = f'{hla2_t}/{hla1_t}'
                hla_postprocess[hla_out] = ["LOH:0", f'{tum_res} ({tum_ratio:.3f})', f'{hla1_c}']
    # methylation analysis

    bed_annot = pysam.TabixFile(f'{annot_dir}/gene_annotate_padding.bed.gz')
    methylation_type_post_gene = typing_post.methylation(tumor_type_dir, tumor_reads_alleles,
                                                         control_type_dir, control_reads_alleles,
                                                         ref_fasta, out_dir, bed_annot,
                                                         "gene")
    bed_annot = pysam.TabixFile(f'{annot_dir}/gene_annotate_promoter2k.bed.gz')
    methylation_type_post_promoter = typing_post.methylation(tumor_type_dir, tumor_reads_alleles,
                                                             control_type_dir, control_reads_alleles,
                                                             ref_fasta, out_dir, bed_annot,
                                                             "promoter")
    return hla_postprocess, methylation_type_post_gene, methylation_type_post_promoter
    # for hla in hla_postprocess.keys():
    #     loh_class, tum_res, ctr_res = hla_postprocess[hla]
    #     hla_postprocess[hla] = f'{loh_class}|{tum_res}|{ctr_res}'
    # return hla_postprocess


def compare_methyl(tumor_methyl, control_methyl, annot_bed, fasta_ref, work_dir) -> tuple:
    # TODO: deprecated
    logger.info("Methylation compare")
    annot_bed_handler = open(annot_bed)
    tumor_results = pysam.TabixFile(tumor_methyl)
    normal_results = pysam.TabixFile(control_methyl)
    # if files are empty break
    if len(tumor_results.contigs) == 0 or len(normal_results.contigs) == 0:
        from pathlib import Path
        Path(f"{work_dir}/mod_parsed.txt").touch()
        return {}, f"{work_dir}/mod_parsed.txt"
    methyl_normal, methyl_tumor, methyl_compare = {}, {}, {}
    gene_divided_regions = ["upstream", "body", "downstream"]
    # modkit output https://github.com/nanoporetech/modkit?tab=readme-ov-file#description-of-bedmethyl-output
    # modkit_results in column 10
    modkit_res_col = 9
    methyl_frac_id = 0
    methyl_site_count_id = 1
    for gene_entry in annot_bed_handler:
        contig, start, end, gene, = gene_entry.rstrip("\n").split("\t")
        methyl_normal[gene] = {"upstream": 0, "body": 0, "downstream": 0, "output": ""}
        methyl_tumor[gene] = {"upstream": 0, "body": 0, "downstream": 0, "output": ""}
        methyl_compare[gene] = {"upstream": "0|0/0", "body": "0|0/0", "downstream": "0|0/0", "output": ""}
        gene_regions = {
            "upstream": f'{contig}:{int(start)-padding}-{int(start)-1}',
            "body": f'{contig}:{start}-{end}',
            "downstream": f'{contig}:{int(end)+1}-{int(end)+padding}'
        }
        # fetch upstream, gene body and downstream the genes
        for region in gene_divided_regions:
            # Tumor
            counter, methyl_frac, mod_bases, canon_bases = 0, 0, 0, 0
            for line in tumor_results.fetch(region=gene_regions[region]):
                line_tab = line.split("\t")
                if len(line_tab) == 10:
                    # older version of modkit
                    _, mfrac, nmod, ncanon, _, _, _, _, _ = line_tab[modkit_res_col].split(" ")
                elif len(line_tab) == 18:
                    # newer version of modkit (v0.5.0)
                    _, mfrac, nmod, ncanon, _, _, _, _, _ = line_tab[modkit_res_col:]
                else:
                    logger.error(f'{len(line_tab)}, {line_tab}')
                    sys.exit(1)
                methyl_frac += float(mfrac)
                mod_bases += int(nmod)
                canon_bases += int(ncanon)
                counter += 1
            methyl_tumor[gene][region] = [methyl_frac/counter, counter] if counter > 0 else [0, 0]
            # normal
            counter, methyl_frac, mod_bases, canon_bases = 0, 0, 0, 0
            for line in normal_results.fetch(region=gene_regions[region]):
                line_tab = line.split("\t")
                if len(line_tab) == 10:
                    # older version of modkit
                    _, mfrac, nmod, ncanon, _, _, _, _, _ = line_tab[modkit_res_col].split(" ")
                elif len(line_tab) == 18:
                    # newer version of modkit (v0.5.0)
                    _, mfrac, nmod, ncanon, _, _, _, _, _ = line_tab[modkit_res_col:]
                else:
                    logger.error("")
                    sys.exit(1)
                methyl_frac += float(mfrac)
                mod_bases += int(nmod)
                canon_bases += int(ncanon)
                counter += 1
            methyl_normal[gene][region] = [methyl_frac/counter, counter] if counter > 0 else [0, 0]
            # T/N: if x>1 higher methylation in tumor, else higher methylation in normal
            # https://pubmed.ncbi.nlm.nih.gov/17334365/
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6791695/
            # strong CpG island promoters are mostly unmethylated, even when inactive (^m ~ no activity)
            # Silencing gene promoters, usually as a result of extensive DNA methylation near the transcription start
            # site (TSS) in regions rich in CpG
            # CpG-poor promoters hypermethylated in somatic cells, which does not preclude their activity. (^m ~ unaff)
            # lack methylation is related to "activity"
            methyl_bang, ratio_out = "", ""
            methyl_tumor_frac = methyl_tumor[gene][region][methyl_frac_id]
            counter_tumor = methyl_tumor[gene][region][methyl_site_count_id]
            methyl_normal_frac = methyl_normal[gene][region][methyl_frac_id]
            if methyl_frac > 0:
                methyl_ratio = methyl_tumor_frac/methyl_normal_frac
                if methyl_tumor_frac+methyl_normal_frac > 80.0:
                    if 1/2 > methyl_ratio and (counter_tumor+counter)/2.0 > 7.0:
                        methyl_bang = "[!]"
                    elif 2 < methyl_ratio and (counter_tumor+counter)/2.0 > 7.0:
                        methyl_bang = "[^]"
                ratio_out = f'{methyl_ratio:.2f}|{methyl_tumor_frac:.2f}/{methyl_normal_frac:.2f}'
            else:
                if methyl_tumor_frac > 50.0 and counter_tumor > 7:
                    methyl_bang = "^"
                ratio_out = f'{methyl_tumor_frac:.2f}/{methyl_normal_frac:.2f}'
            methyl_compare[gene][region] = [f'{methyl_bang}{ratio_out}', f'({counter_tumor}, {counter})']
        methyl_tumor[gene]["output"] = "/".join([f'{methyl_tumor[gene][region]}' for region in gene_divided_regions])
        methyl_normal[gene]["output"] = "/".join([f'{methyl_normal[gene][region]}' for region in gene_divided_regions])
        methyl_compare[gene]["output"] = "\t".join([" ".join(methyl_compare[gene][region]) for region
                                                    in gene_divided_regions])
    annot_bed_handler.close()
    tumor_results.close()
    normal_results.close()
    # run modkit dmr pair
    logger.debug(f'Running modkit dmr pair: {control_methyl} v. {tumor_methyl}')
    modkit_dmr_cmd = (f"modkit dmr pair -a {control_methyl} -b {tumor_methyl} -r {annot_bed} --ref {fasta_ref} "
                      f"--base C -o {work_dir}/mod.txt --log-filepath {work_dir}/dmr.log --force")
    run_cmd = subprocess.run(f'{modkit_dmr_cmd}', shell=True, capture_output=True, text=True)
    if run_cmd.stderr != "":
        logger.warning(f'{Col.Fg.orange}{run_cmd.stderr}{Col.end}')
    modkit_dmr_file = open(f"{work_dir}/mod.txt", "r")
    modkit_dmr_parse = open(f"{work_dir}/mod_parsed.txt", "w")
    min_mod_fracc_diff = 20  # 20%
    max_site_diff = 30  # 30%
    min_sites = 50
    modkit_dmr_parse.write(f'#gene\tregion\tlength\tscore\tmod_frac_control\tn_control\tmod_frac_tumor\tn_tumor\t'
                           f'mod_frac_diff(C-T)|av_sites_diff\n')
    for line in modkit_dmr_file:
        (chrom, start, end, gene, score, _, n_allele1, _, n_allele2, _, p_allele1, p_allele2,
         mod_frac_allele1, mod_frac_allele2, _, _, _, _) = line.rstrip("\n").split("\t")
        score = float(score)
        start, end = int(start), int(end)
        if "," in n_allele1:
            mod1, mod2 = n_allele1.split(",")
            mod_m = mod1 if "m" in mod1 else mod2
            mod1, mod2 = p_allele1.split(",")
            modf_m = mod1 if "m" in mod1 else mod2
        else:
            mod_m = n_allele1 if "m" in n_allele1 else "m:0"
            modf_m = p_allele1 if "m" in p_allele1 else "m:0.0"
        n_allele1 = int(mod_m.split(":")[1])
        p_allele1 = float(modf_m.split(":")[1]) * 100
        if "," in n_allele2:
            mod1, mod2 = n_allele2.split(",")
            mod_m = mod1 if "m" in mod1 else mod2
            mod1, mod2 = p_allele2.split(",")
            modf_m = mod1 if "m" in mod1 else mod2
        else:
            mod_m = n_allele2 if "m" in n_allele2 else "m:0"
            modf_m = p_allele2 if "m" in p_allele2 else "m:0.0"
        n_allele2 = int(mod_m.split(":")[1])
        p_allele2 = float(modf_m.split(":")[1]) * 100
        av_sites_n = (n_allele1 + n_allele2) / 2
        av_sites_diff = abs(n_allele1 - n_allele2) / av_sites_n * 100
        if min(n_allele1, n_allele2) <= 100:
            check_site_diff = av_sites_diff < max_site_diff
        else:
            check_site_diff = True
        if abs(p_allele1 - p_allele2) >= min_mod_fracc_diff and check_site_diff and av_sites_n > min_sites:
            modkit_dmr_parse.write(f'{gene}\t{chrom}:{start}-{end}\t{end - start}\t{score:0.3f}\t'
                                   f'{p_allele1:0.2f}\t{n_allele1}\t{p_allele2:0.2f}\t{n_allele2}\t'
                                   f'{(p_allele1 - p_allele2):0.2f}|{av_sites_diff:0.3f}|'
                                   f'{av_sites_n:0.0f}\n')
    modkit_dmr_file.close()
    modkit_dmr_parse.close()
    return methyl_compare, f"{work_dir}/mod_parsed.txt"


def load_modkit_parsed(modkit_drm_parsed: str) -> dict:
    # TODO: UPDATE for results file
    file_read = open(modkit_drm_parsed, "r")
    modkit_drm_parsed_dict = {}
    for line in file_read:
        line_split = line.rstrip("\n").split("\t")
        gene_name = line_split[0]
        if gene_name not in modkit_drm_parsed_dict:
            modkit_drm_parsed_dict[gene_name] = line_split
    # gene, region, len, score, mod1, cpg1, mod2, cpg2, mod_diff|cpg_diff
    return modkit_drm_parsed_dict


def loss(user_args: LossParameters, script_path):
    if user_args.as_dev:
        logger.setLevel(logging.DEBUG)
    logger.debug(f"Running debug mode")
    annot_paths = BinPathsAnnot(script_path, user_args.as_dev)
    annot_paths.fasta_reference = user_args.reference
    logger.info(f'{Col.bold}Starting {Col.Fg.cyan}SCHLoNS loss{Col.end} '
                f'{Col.bold}with parameters: {user_args.command_string}{Col.end}')
    # Work dir ###########################################################
    if not os.path.exists(user_args.work_dir):
        try:
            os.mkdir(user_args.work_dir)
            os.chdir(user_args.work_dir)
        except OSError as os_err:
            logger.error(f'{Col.bold}{Col.Fg.red}{user_args.work_dir} not found and could not be created'
                         f' {os_err}. ID:404-1{Col.end}')
            sys.exit(1)
    else:
        os.chdir(user_args.work_dir)
    handler_annot: pysam.TabixFile
    if os.path.exists(annot_paths.annotation.hla_bed_pad_gz):
        logger.info(f'Using annotation: {annot_paths.annotation.hla_bed_pad_gz}')
        try:
            # data is bed, chr,start,end,gene_name,strand,gene_id
            handler_annot = pysam.TabixFile(annot_paths.annotation.hla_bed_pad_gz)
        except OSError as os_err:
            logger.error(f'{Col.bold}{Col.Fg.red}Could not open {annot_paths.annotation.hla_bed_pad_gz}'
                         f'annotation.\n{os_err}\nID:503-3{Col.end}')
    else:
        logger.error(f'{Col.bold}{Col.Fg.red}{annot_paths.annotation.hla_bed_pad_gz} not found and needed. '
                     f'ID:404-3{Col.end}')
        sys.exit(1)
    # Tumor + Normal    ######################################################
    tumor_sample = AnalysisFiles()
    normal_sample = AnalysisFiles()
    tumor_sample.args_check_somatic(user_args, "tumor")
    normal_sample.args_check_somatic(user_args, "normal")
    # Report + Output    #####################################################
    output_report = f'{user_args.work_dir}/{user_args.output_prefix}_report.html'
    handler_report = open(output_report, "wt")
    html_report = HTMLResult(user_args.output_prefix)
    logger.info(f'using {output_report} as output')
    output_file = f'{user_args.work_dir}/{user_args.output_prefix}_results.tsv'
    handler_output = open(output_file, "wt")
    handler_output.write(f'Sample: {user_args.output_prefix}\n')
    """
    Note: all files inputs used
    handler_output   -> <OUT>_report.html
    handler_fasta    -> <IN>.fasta.gz | <IN>.fasta.gz.fai
    handler_annot    -> IN>.bed
    tumor|normal.handler_sv          -> sv.vcf.gz
    tumor|normal.handler_snf         -> sv.snf
    tumor|normal.handler_snv         -> snv.vcf.gz
    tumor|normal.handler_cnv         -> cnv.vcf.gz
    tumor|normal.handler_spc         -> cnv.spc.gz
    tumor|normal.handler_cov         -> coverage.regions.bed.gz
    tumor|normal.handler_methyl      -> align_mod_modkit.bed.gz
    tumor|normal.handler_methyl_hap  -> results_genes.txt, single file
    tumor|normal.handler_hlatype     -> typing.txt
    """
    # Section: Analysis pairs ###########################
    comparison_sv = somatic_sv(tumor_sample.snf, normal_sample.snf, annot_paths.annotation.hla_bed_pad_gz,
                               annot_paths.fasta_reference, tumor_sample.snf_mosaic)
    logger.debug(comparison_sv)
    comparison_cnv = somatic_cnv(tumor_sample.cnv, annot_paths.annotation.hla_bed_pad_gz)
    logger.debug(comparison_cnv)
    comparison_type = merge_type(tumor_sample.hlatype, normal_sample.hlatype, annot_paths.annotation.hla_ciwd)
    logger.debug(comparison_type)
    align_ratio_type, methylation_by_type_gene, methylation_by_type_promoter = post_type(tumor_sample.hlatype,
                                                                                         normal_sample.hlatype,
                                                                                         annot_paths.annotation.hla,
                                                                                         annot_paths.fasta_reference,
                                                                                         user_args.work_dir)

    logger.debug(align_ratio_type)
    logger.debug(methylation_by_type_gene)
    logger.debug(methylation_by_type_promoter)
    # TODO: fix this
    # methyl_output = open(f"methyl_results.txt", "w")
    # for line in tumor_sample.handler_methyl_hap:
    #     methyl_output.write(line)
    # methyl_output.close()
    # TODO: update analysis for phased/allele specific methylation
    comparison_methylation, modkit_dmr = compare_methyl(tumor_sample.methyl, normal_sample.methyl,
                                                        annot_paths.annotation.hla_bed_simple,
                                                        annot_paths.fasta_reference, user_args.work_dir)
    modkit_dmr_res = load_modkit_parsed(modkit_dmr)
    logger.info(f'All files loaded {Col.bold}successfully{Col.end}')
    # Section: Summary ###########################
    # ### Section: HLA region ####################
    # TODO: Update table output
    # header_out = ("Gene\tGenomic region\tSV 50+\tCNV 100kb+\tLoH 1Mb+\tTyping tumor\tTyping-mapped tumor"
    #               "\tTyping control\tTyping-mapped control\tPromoter 5mC%\tGene 5mC%\tDownstream 5mC%\n")
    # handler_output.write(header_out)
    genes_output_plot = []
    genes_output_write = []
    hla_output_write = []
    results_save = {}
    logger.debug(comparison_sv)
    for gene_entry in handler_annot.fetch():
        gene_results = GeneResults(gene_entry.rstrip("\n").split("\t"))
        if (comparison_sv[gene_results.name]["init"] or comparison_cnv[gene_results.name]["init"] or
                gene_results.name in modkit_dmr_res or gene_results.name in HLA_GENES):
            cnv = comparison_cnv[gene_results.name]["cnv"]
            sv = comparison_sv[gene_results.name]["sv"]
            me = comparison_methylation[gene_results.name]["output"] \
                if gene_results.name in comparison_methylation else "."
            gene_results.set_results({
                "snv": "NA",
                "coverage": np.nan,
                "cnv": cnv,
                "loh": "",
                "sv": sv,
                "methylation": me,
                "mdr": modkit_dmr_res[gene_results.name] if gene_results.name in modkit_dmr_res else "."
            })
            if gene_results.name in HLA_GENES:
                type_map_class, tumor_map, control_map = align_ratio_type[gene_results.name]
                gene_results.set_hla(comparison_type[gene_results.name]["class"],
                                     comparison_type[gene_results.name]["tumor"],
                                     comparison_type[gene_results.name]["normal"],
                                     type_map_class, tumor_map, control_map)
            if gene_results.name in HLA_GENES:
                hla_output_write.append(f'*{gene_results.print_output()}')
            else:
                genes_output_write.append(gene_results.print_output())
            genes_output_plot.append(gene_results.name)
            results_save[gene_results.name] = gene_results
    for gene_output_text in hla_output_write + genes_output_write:
        handler_output.write(gene_output_text)
    handler_output.close()
    # we will re-use:
    #    comparison_sv, comparison_cnv, comparison_type, comparison_methylation
    # Section: END Summary #######################

    # Section: Plots #############################
    # Data for HLA regions plot
    hla_region_plot_data = {
        tumor: HLALocusResults(),
        normal: HLALocusResults()
    }
    hla_region_plot_data[tumor].set_coverage(tumor_sample.handler_cov)
    hla_region_plot_data[normal].set_coverage(normal_sample.handler_cov)
    hla_region_plot_data[tumor].set_cnv_loh(tumor_sample.handler_cnv)
    hla_region_plot_data[normal].set_cnv_loh(normal_sample.handler_cnv)
    hla_region_plot_data[tumor].set_sv(tumor_sample.handler_sv)
    hla_region_plot_data[normal].set_sv(normal_sample.handler_sv)
    hla_region_plot_data[tumor].set_snv(tumor_sample.handler_snv)
    hla_region_plot_data[normal].set_snv(normal_sample.handler_snv)
    hla_region_plot_data[tumor].set_methylation(tumor_sample.handler_methyl)
    hla_region_plot_data[normal].set_methylation(normal_sample.handler_methyl)

    # HLA region
    plot_region = HLAPlot(tumor_sample.name, normal_sample.name, True, tumor_sample.loh_af_threshold)
    plot_region.output_directory = user_args.work_dir
    plot_region.plot_hla_region(plot_results=hla_region_plot_data)
    # HLA genes
    tumor_coverage_new, control_coverage_new = coverage_finer(tumor_sample.hlatype,
                                                              normal_sample.hlatype,
                                                              user_args.work_dir)
    hla_region_plot_data[tumor].set_coverage(tumor_coverage_new)     # tumor_sample.handler_cov)
    hla_region_plot_data[normal].set_coverage(control_coverage_new)  # normal_sample.handler_cov
    for each_hla in HLA_GENES:
        plot_hla = HLAPlot(tumor_sample.name, normal_sample.name, False, tumor_sample.loh_af_threshold)
        plot_hla.output_directory = user_args.work_dir
        plot_hla.plot_gene_result(res_data=hla_region_plot_data, gene_name=each_hla,
                                  gene_info=results_save[each_hla])

    # All genes
    handler_annot.seek(0)
    for each_gene in results_save.keys():
        if each_gene not in HLA_GENES:
            hla_region_plot_data[tumor].set_coverage(tumor_sample.handler_cov, results_save[each_gene].region)
            hla_region_plot_data[normal].set_coverage(normal_sample.handler_cov, results_save[each_gene].region)
            hla_region_plot_data[tumor].set_methylation(tumor_sample.handler_methyl, results_save[each_gene].region)
            hla_region_plot_data[normal].set_methylation(normal_sample.handler_methyl, results_save[each_gene].region)
            plot_hla = HLAPlot(tumor_sample.name, normal_sample.name)
            plot_hla.output_directory = user_args.work_dir
            plot_hla.plot_gene_result(res_data=hla_region_plot_data, gene_name=each_gene,
                                      gene_info=results_save[each_gene])
    # Section: END Plots #########################
    # Finish
    logger.info("Closing/cleaning")
    tumor_sample.close_handlers()
    normal_sample.close_handlers()
    handler_report.close()
    handler_output.close()
    logger.info("Loss analysis finished")
