import os
from pathlib import Path
from dataclasses import dataclass
from config.parameters import SomaticParameters
from utils.logger_debug import setup_log


class BinPathsAnnot(object):
    # init
    def __init__(self, self_path: Path, as_debug: bool):
        self.logger = setup_log(__name__, debug=as_debug)
        self.as_dev = as_debug
        self.schlons_py = str(self_path)
        self_path_len = len(self.schlons_py.split("/"))
        self.schlons_dir = "/".join(self.schlons_py.split("/")[:self_path_len - 1])
        self.logger.debug(self.schlons_dir)
        # annotation files
        self.annotation = Annotation()
        self.annotation.directory = f'{self.schlons_dir}/data/annotation'
        self.annotation.hla = f'{self.annotation.directory}/hla_genes'
        self.annotation.snv_bed = f'{self.annotation.directory}/chunks'
        self.annotation.cnv_blacklist = f'{self.annotation.directory}/spectre/grch38_blacklist_spectre.bed'
        self.annotation.cnv_meta_ref = f'{self.annotation.directory}/spectre/grch38_spectre.mdr'
        self.annotation.hla_bed_pad_gz = f'{self.annotation.directory}/hla_genes/gene_annotate_padding.bed.gz'
        self.annotation.hla_bed_pad = f'{self.annotation.directory}/hla_genes/gene_annotate_padding.bed'
        self.annotation.hla_bed = f'{self.annotation.directory}/hla_genes/gene_annotate.bed'
        self.annotation.hla_bed_simple = f'{self.annotation.directory}/hla_genes/gene_annotate_simple.bed'
        self.annotation.hla_ciwd = f'{self.annotation.directory}/hla_genes/hla_CIWD3.json'
        # bin/scripts
        self.bin = BinRun(f'{self.schlons_dir}/bin')
        # HLA locus
        self.HLA_REGION = "chr6:26000000-36000000"
        # reference genome and CNV blacklist
        self.fasta_reference = ""

    def print_starting(self, as_debug):
        self.logger.info(f'Reference genome is {self.fasta_reference} ...')
        self.logger.info(f'Working in debug mode ...') if as_debug else None
        self.logger.info(f'Bin dir is {self.bin.directory}')
        self.logger.info(f'HLA genes are in {self.annotation.hla_bed}')

    def snv_bed_regions(self):
        self.logger.debug(f'snv_bed for clair3: {self.annotation.snv_bed}')
        bed_chunks_list = os.listdir(self.annotation.snv_bed)
        bed_chunks_list.sort()
        bed_chunks_for_snv = [f'{self.annotation.snv_bed}/{bed_file}'
                              for bed_file in bed_chunks_list if "bed" in bed_file]
        return bed_chunks_for_snv

    def __repr__(self):
        debug_status = f'Working in debug mode...\n' if self.as_dev else ""
        return (f'Reference genome is {self.fasta_reference}\n{debug_status}'
                f'HLA genes are in {self.annotation.hla_bed}')


class PipelineFiles(object):
    # init
    def __init__(self, user_args: SomaticParameters):
        self.logger = setup_log(__name__, user_args.as_dev)
        self.sample_name = user_args.sample_name
        self.analysis_type = user_args.subcommand
        # working dir
        if not os.path.exists(user_args.work_dir):
            self.logger.info(f'Input dir {user_args.work_dir} do not exists, creating...')
            os.mkdir(user_args.work_dir)
        self.work_dir = user_args.work_dir
        # Samples
        # tumor by default
        self.tumor = ResultsPreprocess("tumor", f'{self.work_dir}/tumor')
        self.tumor.set_input_type(user_args, True)
        if not os.path.exists(self.tumor.dir_work):
            os.mkdir(self.tumor.dir_work)
        # check for control being use
        if self.analysis_type == "somatic":
            self.control = ResultsPreprocess("control", f'{self.work_dir}/control')
            self.control.set_input_type(user_args, False)
            if not os.path.exists(self.control.dir_work):
                os.mkdir(self.control.dir_work)
        # guppy version
        self.guppy_ver = user_args.guppy_ver

    def print_starting(self):
        self.logger.info(f'Working directory is {self.work_dir}/')
        if self.analysis_type == "somatic":
            self.logger.info(f'Somatic analysis of {self.sample_name}')
        elif self.analysis_type == "tumor":
            self.logger.info(f'Tumor-only analysis of {self.sample_name}')
        self.logger.info(f'Using guppy model version {self.guppy_ver}')

    def mod_as_bam(self, sample: str = ""):
        if sample == "tumor":
            self.tumor.file_bam = self.tumor.file_bam_mod
            self.tumor.dir_aln = self.tumor.dir_methyl
        elif sample == "control":
            self.control.file_bam = self.control.file_bam_mod
            self.control.dir_aln = self.control.dir_methyl
        else:
            self.logger.warning(f'{sample} not recognized')


@dataclass
class Annotation:
    directory: str
    snv_bed: str
    cnv_blacklist: str
    cnv_meta_ref: str
    hla_bed: str
    hla_ciwd: str
    hla_bed_pad: str
    hla_bed_pad_gz: str
    hla_bed_simple: str

    def __init__(self):
        pass


class BinRun:
    directory: str
    aln = f'run_submit-minimap2.sh'
    cov = f'run_submit-mosdepth.sh'
    snv = f'run_submit-clair3.sh'
    snv_merge = f'run_submit-clair3_merge_snv.sh'
    snv_phase = f'run_submit-phasing.sh'
    sv = f'run_submit-sniffles2.sh'
    cnv = f'run_submit-spectreCNV.sh'
    mthaln = f'run_submit-minimap2_mod.sh'
    methl = f'run_submit-methyl.sh'
    typing = f'run_submit-hlatyping.sh'
    me_hap = f'run_submit-methyl_hap.sh'
    loss = f'run_submit-loss.sh'

    def __init__(self, bin_dir):
        self.directory = bin_dir
        self.aln = f'{self.directory}/{self.aln}'
        self.cov = f'{self.directory}/{self.cov}'
        self.snv = f'{self.directory}/{self.snv}'
        self.snv_merge = f'{self.directory}/{self.snv_merge}'
        self.snv_phase = f'{self.directory}/{self.snv_phase}'
        self.sv = f'{self.directory}/{self.sv}'
        self.cnv = f'{self.directory}/{self.cnv}'
        self.mthaln = f'{self.directory}/{self.mthaln}'
        self.methl = f'{self.directory}/{self.methl}'
        self.typing = f'{self.directory}/{self.typing}'
        self.me_hap = f'{self.directory}/{self.me_hap}'
        self.loss = f'{self.directory}/{self.loss}'


class ResultsPreprocess:
    name: str = ""
    # inputs
    fastq: str = ""
    modbam: str = ""
    modbam_file: bool = False
    skip_fastq: bool = False
    skip_modbam: bool = False
    # directories
    dir_work: str = ""
    dir_aln: str = "align"
    dir_cov: str = "coverage"
    dir_cnv: str = "cnv"
    dir_snv: str = "snv"
    dir_sv: str = "sv"
    dir_methyl: str = "methylation"
    dir_me_phase: str = "methylation/methyl_phase"
    dir_type: str = "typing"
    dir_analysis: str = "analysis"
    # files to be used
    file_bam: str = "align.bam"
    file_cov: str = "coverage.regions.bed.gz"
    file_cnv: str = "cnv.vcf.gz"
    file_spc: str = "cnv.spc.gz"
    file_sv: str = "sv.vcf.gz"
    file_sv_mosaic: str = "sv_mosaic.vcf.gz"
    file_snf: str = "sv.snf"
    file_snf_mosaic: str = "sv_mosaic.snf"
    file_snv: str = "snv.vcf.gz"
    file_snv_phased: str = "snv_phased.vcf.gz"
    file_type: str = "hla.result.details.txt"
    file_methyl_genes: str = "align_mod_modkit_genes.bed.gz"
    file_methyl: str = "align_mod_modkit.bed.gz"
    file_bam_mod: str = "align_mod.bam"
    file_bam_mod_hap: str = "align_mod_genes_hp.bam"
    file_me_phase_dmr_res: str = "results_genes.txt"
    file_me_phase_dmr_stats: str = "results_stats.txt"

    def __init__(self, sample_name, work_dir):
        self.name = sample_name
        self.dir_work = work_dir
        # dirs
        self.dir_aln = f'{self.dir_work}/{self.dir_aln}'
        self.dir_cov = f'{self.dir_work}/{self.dir_cov}'
        self.dir_cnv = f'{self.dir_work}/{self.dir_cnv}'
        self.dir_snv = f'{self.dir_work}/{self.dir_snv}'
        self.dir_sv = f'{self.dir_work}/{self.dir_sv}'
        self.dir_methyl = f'{self.dir_work}/{self.dir_methyl}'
        self.dir_type = f'{self.dir_work}/{self.dir_type}'
        self.dir_analysis = f'{self.dir_work}/{self.dir_analysis}'

    def set_input_type(self, user_args: SomaticParameters, is_tumor: bool = True):
        file_input = user_args.tumor if is_tumor else user_args.control
        file_name, file_ext = os.path.splitext(file_input) if is_tumor else os.path.splitext(file_input)
        if "bam" == file_ext:
            self.modbam = file_input
            self.modbam_file = True
            self.skip_fastq = True
        elif "gz" == file_ext:
            _, file_ext2 = os.path.splitext(file_name)
            if "fq" == file_ext2 or "fastq" == file_ext2:
                self.skip_modbam = True
                self.fastq = file_input
        elif "fq" == file_ext or "fastq" == file_ext:
            self.skip_modbam = True
            self.fastq = file_input
