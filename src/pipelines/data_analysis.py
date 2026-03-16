#!/usr/bin/env python3
import os
import datetime
import subprocess
from pathlib import Path
from config.parameters import SomaticParameters
from config.parameters import LossParameters
from utils.term_colors import TermColors as Col
from utils.rng_string import rng_string_gen
from utils.logger_debug import setup_log
from .data_process import PipelineFiles
from .data_process import BinPathsAnnot
from .data_process import ResultsPreprocess
from utils.submit_jobs import SubmitJobsSlurm as SubmitJobs


class TumorLensProcessing(object):
    def __init__(self, user_args: SomaticParameters, script_path: Path):
        # #########################################################################
        # logger
        self.logger = setup_log(__name__, user_args.as_dev)
        self.logger.info(f"Starting {Col.Fg.lightblue}process steps {Col.end} ...")
        # Starting ....
        self.user_args = user_args
        self.sample_jobs_id = rng_string_gen(12)
        self.time_start = datetime.datetime.now()
        # #########################################################################
        self.bin_paths = BinPathsAnnot(script_path, user_args.as_dev)
        self.bin_paths.fasta_reference = user_args.reference
        self.bin_paths.print_starting(user_args.as_dev)
        # file_use.tumor & file_use.control
        self.files_use = PipelineFiles(user_args)
        self.files_use.print_starting()
        self.logger.info(f'Sample ID: {self.sample_jobs_id}')
        # for cov and cnv
        self.window_size = "1000"
        # dev
        self.dry_run = user_args.dry_run

    @staticmethod
    def get_jobid_from_stdout(job_stdout: str = "", as_list: bool = False) -> str:
        job_id_list = []
        if job_stdout != "":
            file_jobs_list = job_stdout if as_list else [job_stdout]
            for each_job in file_jobs_list:
                try:
                    [_, _, _, job_id] = each_job.rstrip("\n").split(" ")
                except ValueError:
                    job_id = ""
                job_id_list.append(job_id)
        else:
            job_id_list.append("")
        if len(job_id_list) == 1:
            return job_id_list[0]
        else:
            return ",".join(job_id_list)

    def align_bam(self, sample: ResultsPreprocess) -> str:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}alignment{Col.end} process '
                         f'for {Col.bold}{Col.Fg.cyan}{sample.name}{Col.end}')
        os.chdir(sample.dir_work)
        if not os.path.exists(sample.dir_aln):
            self.logger.info(f'Making directory "{sample.dir_aln}"')
            os.mkdir(sample.dir_aln)
        os.chdir(sample.dir_aln)
        if os.path.exists(sample.file_bam):
            self.logger.info(f'{Col.bold}{Col.Fg.red}Skipping mapping, {sample.file_bam} '
                             f'exists{Col.end}')
            return ""
        # aln.sh REFERENCE=$1 READS=$2 OUT=$3
        cmd_aln = f'{self.bin_paths.bin.aln} {self.bin_paths.fasta_reference} {sample.fastq} '
        if self.dry_run:
            Path(f'{sample.dir_aln}/align.bam').touch()
            return "0001"
        job_aln = SubmitJobs(self.logger)
        job_aln.set_details(name=f'aln{self.sample_jobs_id}',
                            output=f'log_minimap2_{self.sample_jobs_id}.out',
                            error=f'log_minimap2_{self.sample_jobs_id}.err',
                            chdir=sample.dir_aln)
        run_cmd = subprocess.run(job_aln.make_submit_job(cmd_aln), shell=True, capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_aln.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        return job_aln.jobid

    def methylation_aln(self, sample: ResultsPreprocess) -> str:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}methylation{Col.end} pre-process '
                         f'for {Col.bold}{Col.Fg.cyan}{sample.name}{Col.end}')
        os.chdir(sample.dir_work)
        if not os.path.exists(sample.dir_methyl):
            self.logger.info(f'Making directory "{sample.dir_methyl}"')
            os.mkdir(sample.dir_methyl)
        os.chdir(sample.dir_methyl)
        self.logger.debug(f'Checking if {sample.file_bam_mod} exists')
        if os.path.exists(sample.file_bam_mod):
            self.logger.info(f'{Col.bold}{Col.Fg.red}Skipping modbam mapping, {sample.file_bam_mod} '
                             f'exists{Col.end}')
            return ""
        # 1. cat all modbam + map
        # modbam_map.sh REFERENCE=$1  READS_DIR=$2 GENES_BED=$3 FINAL_DIR=$4
        cmd_modbam_map = f'{self.bin_paths.bin.mthaln}  {self.bin_paths.fasta_reference} ' + \
                         f'{sample.modbam}  {sample.dir_methyl}'
        if self.dry_run:
            Path(f'{sample.dir_methyl}/align_mod.bam').touch()
            return "0002"
        job_modbam_map = SubmitJobs(self.logger)
        job_modbam_map.set_details(name=f'mth{self.sample_jobs_id}',
                                   output=f'log_modbam_minimap2_{self.sample_jobs_id}.out',
                                   error=f'log_modbam_minimap2_{self.sample_jobs_id}.err',
                                   chdir=sample.dir_methyl)
        run_cmd = subprocess.run(job_modbam_map.make_submit_job(cmd_modbam_map), shell=True,
                                 capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_modbam_map.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        self.logger.debug(f'Job ID: {job_modbam_map.jobid}')
        return job_modbam_map.jobid

    def methylation_stats(self, sample: ResultsPreprocess, dependency_job: str = "") -> str:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}methylation{Col.end} stats '
                         f'for {Col.bold}{Col.Fg.cyan}{sample.name}{Col.end}')
        os.chdir(sample.dir_work)
        if not os.path.exists(sample.dir_methyl):
            self.logger.info(f'Making directory "{sample.dir_methyl}"')
            os.mkdir(sample.dir_methyl)
        os.chdir(sample.dir_methyl)
        self.logger.debug(f'Checking if {sample.file_methyl_genes} exists')
        if os.path.exists(sample.file_methyl_genes):
            self.logger.info(f'{Col.bold}{Col.Fg.red}Skipping methyl stats, {sample.file_methyl_genes} '
                             f'exists{Col.end}')
            return ""
        # 1. cat all modbam + map
        # modbam_map.sh REFERENCE=$1  BAM_ALN=$2 GENES_BED=$3 FINAL_DIR=$4
        cmd_methyl_stats = (f'{self.bin_paths.bin.methl}  {self.bin_paths.fasta_reference} {sample.file_bam_mod}  '
                            f'{self.bin_paths.annotation.hla_bed}  {sample.dir_methyl}')
        if self.dry_run:
            Path(f'{sample.dir_methyl}/align_mod_modkit.bed.gz').touch()
            return "0003"
        job_methyl_stats = SubmitJobs(self.logger)
        job_methyl_stats.set_details(name=f'mth{self.sample_jobs_id}',
                                     output=f'log_methyl_stats{self.sample_jobs_id}.out',
                                     error=f'log_methyl_stats{self.sample_jobs_id}.err',
                                     chdir=sample.dir_methyl)
        if not dependency_job:
            self.logger.info(f'Expecting indexed bam ({sample.dir_aln}/{sample.file_bam_mod}) to be present')
        else:
            job_methyl_stats.set_dependencies(f'afterok:{dependency_job}')
        run_cmd = subprocess.run(job_methyl_stats.make_submit_job(cmd_methyl_stats), shell=True,
                                 capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_methyl_stats.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        self.logger.debug(f'Job ID: {job_methyl_stats.jobid}')
        return job_methyl_stats.jobid

    def coverage(self, sample: ResultsPreprocess, dependency_job: str = "") -> str:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}coverage{Col.end} process '
                         f'for {Col.bold}{Col.Fg.cyan}{sample.name}{Col.end}')
        os.chdir(sample.dir_work)
        if not os.path.exists(sample.dir_cov):
            self.logger.info(f'Making directory "{sample.dir_cov}"')
            os.mkdir(sample.dir_cov)
        os.chdir(sample.dir_cov)
        if os.path.exists(sample.file_cov):
            self.logger.info(f'{Col.bold}{Col.Fg.red}Skipping coverage, {sample.file_cov} '
                             f'exists{Col.end}')
            return ""
        # cov.sh FILE_IN_BAM=$1 WINSIZE=$2
        cmd_cov = f'{self.bin_paths.bin.cov}  {sample.dir_aln}/{sample.file_bam}  {self.window_size}'
        if self.dry_run:
            Path(f'{sample.dir_cov}/coverage.regions.bed.gz').touch()
            return "0004"
        job_cov = SubmitJobs(self.logger)
        job_cov.set_details(name=f'cov{self.sample_jobs_id}',
                            output=f'log_mosdepth_{self.sample_jobs_id}.out',
                            error=f'log_mosdepth_{self.sample_jobs_id}.err',
                            chdir=sample.dir_cov)
        if not dependency_job:
            self.logger.info(f'Expecting indexed bam ({sample.dir_aln}/{sample.file_bam}) to be present')
        else:
            job_cov.set_dependencies(f'afterok:{dependency_job}')
        run_cmd = subprocess.run(job_cov.make_submit_job(cmd_cov), shell=True, capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_cov.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        return job_cov.jobid

    def sv(self, sample: ResultsPreprocess, dependency_job: str = "") -> str:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}SV{Col.end} process '
                         f'for {Col.bold}{Col.Fg.cyan}{sample.name}{Col.end}')
        os.chdir(sample.dir_work)
        if not os.path.exists(sample.dir_sv):
            self.logger.info(f'Making directory "{sample.dir_sv}"')
            os.mkdir(sample.dir_sv)
        os.chdir(sample.dir_sv)
        if os.path.exists(sample.file_sv) and os.path.exists(sample.file_snf):
            self.logger.info(f'{Col.bold}{Col.Fg.red}Skipping sv-calling, {sample.file_sv} '
                             f'exists{Col.end}')
            return ""
        # sv.sh  FILE_IN_BAM=$1  FASTAREF=$2
        cmd_sv_call = (f'{self.bin_paths.bin.sv}  {sample.dir_aln}/{sample.file_bam} '
                       f'{self.bin_paths.fasta_reference}  {sample.name}')
        if self.dry_run:
            Path(f'{sample.dir_sv}/sv.vcf.gz').touch()
            Path(f'{sample.dir_sv}/sv.snf').touch()
            return "0005"
        job_sv_call = SubmitJobs(self.logger)
        job_sv_call.set_details(name=f'sv{self.sample_jobs_id}',
                                output=f'log_sniffles2_{self.sample_jobs_id}.out',
                                error=f'log_sniffles2_{self.sample_jobs_id}.err',
                                chdir=sample.dir_sv)
        if not dependency_job:
            self.logger.info(f'Expecting indexed bam ({sample.dir_aln}/{sample.file_bam}) to be present')
        else:
            job_sv_call.set_dependencies(f'afterok:{dependency_job}')
        run_cmd = subprocess.run(job_sv_call.make_submit_job(cmd_sv_call), shell=True, capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_sv_call.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        return job_sv_call.jobid

    def snv(self, sample: ResultsPreprocess, dependency_job: str = "") -> tuple:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}SNV{Col.end} process '
                         f'for {Col.bold}{Col.Fg.cyan}{sample.name}{Col.end}')
        os.chdir(sample.dir_work)
        if not os.path.exists(sample.dir_snv):
            self.logger.info(f'Making directory "{sample.dir_snv}"')
            os.mkdir(sample.dir_snv)
        os.chdir(sample.dir_snv)
        if os.path.exists(sample.file_snv):
            self.logger.info(f'{Col.bold}{Col.Fg.red}Skipping snv-calling, {sample.file_snv} '
                             f'exists{Col.end}')
            return "", ""
        # snv.sh PREFIX_FILE=$1  FILE_IN_BAM=$2  FASTAREF=$3  REGIONS_FILE=$4  OUT_SNV_CLAIR3=$5  GUPPY_USE=$6
        os.mkdir(f'{sample.dir_snv}/chunks') if not os.path.exists(f'{sample.dir_snv}/chunks') else None
        job_snv_call_ids = []
        snv_call_outdirs = []
        genome_snv_chunks = self.bin_paths.snv_bed_regions()
        self.logger.info(f'HLA bed {genome_snv_chunks} ...')
        if self.dry_run:
            Path(f'{sample.dir_snv}/chunk.vcf.gz').touch()
            return "0006", [""]
        for use_chunk in genome_snv_chunks:
            chunk_use_out = os.path.basename(use_chunk).split(".")[0]
            cmd_snv_call = (f'{self.bin_paths.bin.snv}  {sample.dir_aln}/{sample.file_bam}  {use_chunk}  '
                            f'{self.bin_paths.fasta_reference}')
            job_snv_call = SubmitJobs(self.logger)
            job_snv_call.set_details(name=f'snv{self.sample_jobs_id}',
                                     output=f'log_clair3_{self.sample_jobs_id}_{chunk_use_out}.out',
                                     error=f'log_clair3_{self.sample_jobs_id}_{chunk_use_out}.err',
                                     chdir=f'{sample.dir_snv}/chunks')
            if not dependency_job:
                self.logger.info(f'Expecting indexed bam ({sample.dir_aln}/{sample.file_bam}) to be present')
            else:
                job_snv_call.set_dependencies(f'afterok:{dependency_job}')
            run_cmd = subprocess.run(job_snv_call.make_submit_job(cmd_snv_call), shell=True,
                                     capture_output=True, text=True)
            if run_cmd.stderr != "":
                self.logger.error(run_cmd.stderr)
            job_snv_call_ids.append(self.get_jobid_from_stdout(run_cmd.stdout))
            snv_call_outdirs.append(f'snv_{chunk_use_out}')
        job_snv_jobid = ",".join(job_snv_call_ids)
        return job_snv_jobid, snv_call_outdirs

    def snv_merge(self, sample: ResultsPreprocess, dependency_job: str = "", clair_chunk_outdirs=None) -> str:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}SNV merging{Col.end} process '
                         f'for {Col.bold}{Col.Fg.cyan}{sample.name}{Col.end}')
        if clair_chunk_outdirs is None:
            clair_chunk_outdirs = []
        os.chdir(sample.dir_snv)
        if os.path.exists(sample.file_snv):
            self.logger.info(f'{Col.bold}{Col.Fg.red}Skipping snv-calling, {sample.file_snv} '
                             f'exists{Col.end}')
            return ""
        if self.dry_run:
            Path(f'{sample.dir_snv}/snv.vcf.gz').touch()
            return "0007"
        clair_chunk_outdirs.sort()
        chunk_list = open("filelist.txt", "w")
        for each_chunk in clair_chunk_outdirs:
            chunk_list.write(f'{sample.dir_snv}/chunks/{each_chunk}/merge_output.vcf.gz\n')
        chunk_list.close()
        # snv_merge.sh  SNV_OUTPUT=$1  files_merge=$2
        cmd_snv_merge = f'{self.bin_paths.bin.snv_merge}  {sample.file_snv}  filelist.txt'
        job_snv_merge = SubmitJobs(self.logger)
        job_snv_merge.set_details(name=f'mrg{self.sample_jobs_id}',
                                  output=f'log_snv_merge_{self.sample_jobs_id}.out',
                                  error=f'log_snv_merge_{self.sample_jobs_id}.err',
                                  chdir=sample.dir_snv)
        if not dependency_job:
            self.logger.info("Expecting Clair3 calls in 'chunks' to be present")
        else:
            job_snv_merge.set_dependencies(f'afterok:{dependency_job}')
        run_cmd = subprocess.run(job_snv_merge.make_submit_job(cmd_snv_merge), shell=True,
                                 capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_snv_merge.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        return job_snv_merge.jobid

    def snv_phase(self, tumor: ResultsPreprocess, control: ResultsPreprocess, dependency_job: str = "") -> str:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}SNV phasing{Col.end} '
                         f'for {Col.bold}{Col.Fg.red}{tumor.name}{Col.end}/'
                         f'{Col.bold}{Col.Fg.green}{control.name}{Col.end}')
        """
        self.bin_paths.fasta_reference
        self.bin_paths.annotation.hla_bed
        control.dir_work
        tumor.dir_work
        """
        # snv_phasing.sh  ref=$1 bed_annot=$2 tumor=$3 control=$4  script_path=$5
        cmd_phasing = (f'{self.bin_paths.bin.snv_phase}  {self.bin_paths.fasta_reference}  '
                       f'{self.bin_paths.annotation.hla_bed_pad} {tumor.dir_work} {control.dir_work} '
                       f'{self.bin_paths.schlons_dir}')
        if self.dry_run:
            Path(f'{tumor.dir_snv}/snv_phase.vcf.gz').touch()
            return "0008"
        job_snv_phase = SubmitJobs(self.logger)
        job_snv_phase.set_details(name=f'phase{self.sample_jobs_id}',
                                  output=f'log_snv_phase_{self.sample_jobs_id}.out',
                                  error=f'log_snv_phase_{self.sample_jobs_id}.err',
                                  chdir=f'{tumor.dir_work}/snv')
        if not dependency_job:
            self.logger.info("alignment and SNV data need to be available")
        else:
            job_snv_phase.set_dependencies(f'afterok:{dependency_job}')
        run_cmd = subprocess.run(job_snv_phase.make_submit_job(cmd_phasing), shell=True,
                                 capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_snv_phase.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        return job_snv_phase.jobid

    def typing(self, sample: ResultsPreprocess, dependency_job: str = "") -> str:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}HLA typing{Col.end} process '
                         f'for {Col.bold}{Col.Fg.cyan}{sample.name}{Col.end}')
        os.chdir(sample.dir_work)
        if not os.path.exists(sample.dir_type):
            self.logger.info(f'Making directory "{sample.dir_type}"')
            os.mkdir(sample.dir_type)
        os.chdir(sample.dir_type)
        if os.path.exists(sample.file_type):
            self.logger.info(f'{Col.bold}{Col.Fg.red}Skipping typing, {sample.file_type} '
                             f'exists{Col.end}')
            return ""
        # typing.sh  WORKDIR=$1  BAM_IN=$2 | => optional HLABIN=$3  HLATYPE=$4
        cmd_hla_typing = f'{self.bin_paths.bin.typing}  {sample.dir_work}  {sample.dir_aln}/{sample.file_bam} '
        if self.dry_run:
            Path(f'{sample.dir_type}/hla.result.details.txt').touch()
            return "0009"
        job_hla_typing = SubmitJobs(self.logger)
        job_hla_typing.set_details(name=f'typ{self.sample_jobs_id}',
                                   output=f'log_hla_typing_{self.sample_jobs_id}.out',
                                   error=f'log_hla_typing_{self.sample_jobs_id}.err',
                                   chdir=sample.dir_type)
        if not dependency_job:
            self.logger.info("Expecting alignment (bam files) to be present")
        else:
            job_hla_typing.set_dependencies(f'afterok:{dependency_job}')
        run_cmd = subprocess.run(job_hla_typing.make_submit_job(cmd_hla_typing), shell=True,
                                 capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_hla_typing.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        return job_hla_typing.jobid

    def cnv(self, sample: ResultsPreprocess, dependency_job: str = "") -> str:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}CNV{Col.end} process '
                         f'for {Col.bold}{Col.Fg.cyan}{sample.name}{Col.end}')
        os.chdir(sample.dir_work)
        if not os.path.exists(sample.dir_cnv):
            self.logger.info("Making directory 'cnv'")
            os.mkdir(sample.dir_cnv)
        os.chdir(sample.dir_cnv)
        if os.path.exists(sample.file_cnv):
            self.logger.info(f'{Col.bold}{Col.Fg.red}Skipping snv-calling, {sample.file_cnv} '
                             f'exists{Col.end}')
            return ""
        # cnv.sh SAMPLE_ID=$1  COVERAGE_FILE=$2  FASTAREF=$3  REF_META=$4  BLACKLIST=$5  SNV=$6
        if self.dry_run:
            Path(f'{sample.dir_cnv}/cnv.vcf.gz').touch()
            return "0010"
        cmd_cnv_call = f'{self.bin_paths.bin.cnv}  {sample.name}  {sample.dir_cov}/{sample.file_cov}  ' \
                       f'{self.bin_paths.fasta_reference}  {self.bin_paths.annotation.cnv_meta_ref}' \
                       f'{self.bin_paths.annotation.cnv_blacklist}  {sample.dir_snv}/{sample.file_snv} '
        job_cnv_call = SubmitJobs(self.logger)
        job_cnv_call.set_details(name=f'cnv{self.sample_jobs_id}',
                                 output=f'log_spectre_cnv_{self.sample_jobs_id}.out',
                                 error=f'log_spectre_cnv_{self.sample_jobs_id}.err',
                                 chdir=sample.dir_cnv)
        if not dependency_job:
            self.logger.info("Expecting coverage directory (bed files) to be present")
        else:
            job_cnv_call.set_dependencies(f'afterok:{dependency_job}')
        run_cmd = subprocess.run(job_cnv_call.make_submit_job(cmd_cnv_call), shell=True,
                                 capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_cnv_call.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        return job_cnv_call.jobid

    def cnv_cancer(self, tumor: ResultsPreprocess, control: ResultsPreprocess, tumor_purity: float,
                   dependency_job: str = "") -> str:
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}CNV cancer{Col.end} process '
                         f'for {Col.bold}{Col.Fg.red}{tumor.name}{Col.end}/'
                         f'{Col.bold}{Col.Fg.green}{control.name}{Col.end} pair')
        os.chdir(tumor.dir_work)
        if not os.path.exists(tumor.dir_cnv):
            self.logger.info("Making directory 'cnv'")
            os.mkdir(tumor.dir_cnv)
        os.chdir(tumor.dir_cnv)
        if os.path.exists(tumor.file_cnv):
            self.logger.info(f'{Col.bold}{Col.Fg.red}Skipping snv-calling, {tumor.file_cnv} '
                             f'exists{Col.end}')
            return ""
        # cnv.sh SAMPLE_ID=$1  COVERAGE_FILE=$2  FASTAREF=$3  REF_META=$4  BLACKLIST=$5
        #        SNV=$6  TUMOR_PURITY=$7  CANCER=$8
        if self.dry_run:
            Path(f'{tumor.dir_cnv}/cnv.vcf.gz').touch()
            return "0011"
        cmd_cnv_call = (f'{self.bin_paths.bin.cnv}  {tumor.name}  '
                        f'{tumor.dir_cov}/{tumor.file_cov},{control.dir_cov}/{control.file_cov}    '
                        f'{self.bin_paths.fasta_reference}  {self.bin_paths.annotation.cnv_meta_ref}  '
                        f'{self.bin_paths.annotation.cnv_blacklist}  '
                        f'{tumor.dir_snv}/{tumor.file_snv},{control.dir_snv}/{control.file_snv}' 
                        f'{tumor_purity}  1')
        job_cnv_call = SubmitJobs(self.logger)
        job_cnv_call.set_details(name=f'cnvc{self.sample_jobs_id}',
                                 output=f'log_spectre_cnv_{self.sample_jobs_id}.out',
                                 error=f'log_spectre_cnv_{self.sample_jobs_id}.err',
                                 chdir=tumor.dir_cnv)
        if not dependency_job:
            self.logger.info("Expecting files: coverage file for tumor, SNV for tumor and SNV for control")
        else:
            job_cnv_call.set_dependencies(f'afterok:{dependency_job}')
        run_cmd = subprocess.run(job_cnv_call.make_submit_job(cmd_cnv_call), shell=True,
                                 capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_cnv_call.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        return job_cnv_call.jobid

    def methylation_hap(self, tumor: ResultsPreprocess, control: ResultsPreprocess, dependency_job: str = ""):
        self.logger.info(f'Starting {Col.bold}{Col.Fg.blue}Methylation by hap{Col.end} '
                         f'for {Col.bold}{Col.Fg.red}{tumor.name}{Col.end}/'
                         f'{Col.bold}{Col.Fg.green}{control.name}{Col.end}')
        """
        self.bin_paths.fasta_reference
        self.bin_paths.annotation.hla_bed_pad_gz
        control.dir_work
        tumor.dir_work
        """
        # methyl_hap.sh  ref=$1 bed_annot=$2 tumor_hp_bam=$3 control_hopp_bam=$4
        cmd_meth_phasing = (f'{self.bin_paths.bin.me_hap}  {self.bin_paths.fasta_reference}  '
                            f'{tumor.dir_methyl}/{tumor.file_bam_mod_hap}  '
                            f'{control.dir_methyl}/{control.file_bam_mod_hap}  '
                            f'{self.bin_paths.annotation.hla_bed_simple}  {tumor.name},{control.name}')
        if self.dry_run:
            Path(f'{tumor.dir_methyl}/methyl_hps.bed.gz').touch()
            return "0012"
        job_meth_phase = SubmitJobs(self.logger)
        job_meth_phase.set_details(name=f'methap{self.sample_jobs_id}',
                                   output=f'log_meth_phase_{self.sample_jobs_id}.out',
                                   error=f'log_meth_phase_{self.sample_jobs_id}.err',
                                   chdir=tumor.dir_methyl)
        if not dependency_job:
            self.logger.info("alignment and SNV data need to be available")
        else:
            job_meth_phase.set_dependencies(f'afterok:{dependency_job}')
        run_cmd = subprocess.run(job_meth_phase.make_submit_job(cmd_meth_phasing), shell=True,
                                 capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_meth_phase.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        return job_meth_phase.jobid


class TumorLensLoss(object):
    def __init__(self, user_args_loss: LossParameters, job_id: str = "None"):
        # jobID
        self.job_id = job_id
        # args
        self.loss_param: LossParameters = user_args_loss
        # logger
        self.logger = setup_log(__name__, user_args_loss.as_dev)
        self.logger.info("Starting loss analysis")
        self.dry_run = user_args_loss.dry_run

    @staticmethod
    def get_jobid_from_stdout(job_stdout: str = "", as_list: bool = False) -> str:
        job_id_list = []
        if job_stdout != "":
            file_jobs_list = job_stdout if as_list else [job_stdout]
            for each_job in file_jobs_list:
                try:
                    [_, _, _, job_id] = each_job.rstrip("\n").split(" ")
                except ValueError:
                    job_id = ""
                job_id_list.append(job_id)
        else:
            job_id_list.append("")
        if len(job_id_list) == 1:
            return job_id_list[0]
        else:
            return ",".join(job_id_list)

    def submit_loss(self, tumor_json: str, control_json: str, dependencies: list[str],
                    analysis_dir: str):
        """
        ${schlonsbin} loss \
            --work-dir         ${results_dir} \
            --fasta            /mnt/ssd_ubuntu/reference_genomes/grch38/genome/grch38.fasta \
            --output           ${sampleid} \
            --tumor            ${results_dir}/tumor.json \
            --control          ${results_dir}/control.json \
            --debug
            schlons_bin=$1      results_dir=$2
            reference=$3        sample_id=$4
            tumor_json=$5       control_json=$6
            """
        cmd_sv_call = (f'{self.loss_param.command_string}  {self.loss_param.subcommand} {analysis_dir}  '
                       f'{self.loss_param.reference}  {self.loss_param.output_prefix}  '
                       f'{tumor_json}  {control_json}')
        if self.loss_param.dry_run:
            Path(f'{analysis_dir}/analysis.txt').touch()
            return "0013"
        job_sv_call = SubmitJobs(self.logger)
        job_sv_call.set_details(name=f'loss',
                                output=f'log_loss_{self.job_id}.out',
                                error=f'log_loss_{self.job_id}.err',
                                chdir=analysis_dir)
        dependency_job = ",".join(dependencies)
        job_sv_call.set_dependencies(f'afterok:{dependency_job}')
        run_cmd = subprocess.run(job_sv_call.make_submit_job(cmd_sv_call), shell=True, capture_output=True, text=True)
        if run_cmd.stderr != "":
            self.logger.error(run_cmd.stderr)
        job_sv_call.set_jobid(self.get_jobid_from_stdout(run_cmd.stdout))
        return job_sv_call.jobid

