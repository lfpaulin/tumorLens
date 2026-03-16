#!/usr/bin/env python3
import os.path
import sys
import time
from pathlib import Path
from config.config import GetArguments
from config.parameters import SomaticParameters
from config.parameters import LossParameters
from utils.term_colors import TermColors as Col
from pipelines.data_analysis import TumorLensProcessing
from pipelines.data_analysis import TumorLensLoss
from analysis.loss_somatic import loss
from utils.logger_debug import setup_log


script_path = Path(__file__).absolute()
version = "1.0.250801"


def tumor_lens(user_args: SomaticParameters):
    tl_logger = setup_log(__name__, user_args.as_dev)
    tl_logger.info(f'{Col.bold}Starting {Col.Fg.cyan}TumorLens {Col.end}{Col.bold}with parameters:\n'
                   f'>   {user_args.command_string}{Col.end}')
    pipeline = TumorLensProcessing(user_args, script_path)
    time.sleep(2)
    if "somatic" == user_args.subcommand:
        tl_logger.info(f'Running in {Col.bold}somatic analysis{Col.end} for {Col.bold}{user_args.sample_name}{Col.end}')
        tl_logger.info(f'Tumor sample processing')
        tumor_jobs = sample_preprocess(user_args, pipeline, True, {})
        tl_logger.info(f'Control sample processing')
        control_jobs = sample_preprocess(user_args, pipeline, False, tumor_jobs)
        loss_args, tumor_json_path, control_json_path, analysis_dir = make_loss_script(user_args, pipeline)
        loss_submit(loss_args, tumor_json_path, control_json_path, tumor_jobs, control_jobs,
                    analysis_dir, pipeline.sample_jobs_id)
    if "tumor" == user_args.subcommand:
        tl_logger.info(f'Running in {Col.bold}tumor-only mode{Col.end} for {Col.bold}{user_args.sample_name}{Col.end}')
        sample_preprocess(user_args, pipeline, True, {})


def sample_preprocess(user_args: SomaticParameters, pipeline: TumorLensProcessing, is_tumor: bool,
                      tracked_jobs: dict) -> dict:
    preprocess_logger = setup_log(__name__, user_args.as_dev)
    use_sample = pipeline.files_use.tumor if is_tumor else pipeline.files_use.control
    sample_type = "tumor" if is_tumor else "control"
    tumor_data = pipeline.files_use.tumor if not is_tumor else None
    # step 1 preprocess
    # job tracker
    job_tracker = {"sample_id": pipeline.sample_jobs_id}
    preprocess_logger.info(f'Running HLA loss pre-processing with job ID {Col.bold}{pipeline.sample_jobs_id}{Col.end}')
    time.sleep(2)
    # STEP 1: Alignment
    map_dependency = ""
    job_tracker["minimap"] = ""
    if not use_sample. skip_fastq:
        job_tracker["minimap"] = pipeline.align_bam(use_sample)
        map_dependency = job_tracker["minimap"]
        preprocess_logger.debug(f'using minimap dependency: {map_dependency}')
    else:
        preprocess_logger.debug(f'using modbam as aln')
        pipeline.files_use.mod_as_bam(sample_type)
    time.sleep(5)
    # methylation aln merge
    job_tracker["meth_aln"] = ""
    if not use_sample.skip_modbam:
        job_tracker["meth_aln"] = pipeline.methylation_aln(use_sample)
        if map_dependency == "":
            map_dependency = job_tracker["meth_aln"]
    # alignment dependencies and file assign
    time.sleep(10)
    # continue analysis ...
    # Coverage
    job_tracker["mosdepth"] = pipeline.coverage(sample=use_sample, dependency_job=map_dependency)
    time.sleep(5)
    # SV
    job_tracker["sniffles2"] = pipeline.sv(sample=use_sample, dependency_job=map_dependency)
    time.sleep(5)
    # SNV
    [job_tracker["clair3"], clair_chunk_out_dirs] = pipeline.snv(sample=use_sample, dependency_job=map_dependency)
    time.sleep(5)
    # SNV.merge
    snv_dependency = job_tracker["clair3"]
    job_tracker["snv_merge"] = pipeline.snv_merge(sample=use_sample, dependency_job=snv_dependency,
                                                  clair_chunk_outdirs=clair_chunk_out_dirs)
    time.sleep(5)
    # Methylation => TODO: update for phasing
    job_tracker["methyl"] = pipeline.methylation_stats(sample=use_sample, dependency_job=map_dependency)
    time.sleep(5)
    # typing
    job_tracker["hla_typing"] = pipeline.typing(sample=use_sample, dependency_job=map_dependency)
    time.sleep(5)
    # CNV
    if "control" == sample_type:
        # Tumor sample CNV calling <== coverage tumor, SNV tumor + SNV control
        cnv_dependency = ",".join([job_id for job_id in [tracked_jobs["mosdepth"], tracked_jobs["snv_merge"],
                                                         job_tracker["snv_merge"],] if job_id != ""])
        job_tracker["spectre_t"] = pipeline.cnv_cancer(tumor=pipeline.files_use.tumor,
                                                       tumor_purity=user_args.tumor_purity,
                                                       control=pipeline.files_use.control,
                                                       dependency_job=cnv_dependency)
        # Control sample CNV calling is not needed, cancer mode of spectre already gives somatic calling
        # SNV.phase
        phase_dependency = ",".join([job_id for job_id in [tracked_jobs["snv_merge"], job_tracker["snv_merge"],
                                                           tracked_jobs["meth_aln"],  job_tracker["meth_aln"]]
                                     if job_id != ""])
        job_tracker["snv_phase"] = pipeline.snv_phase(tumor=tumor_data, control=use_sample,
                                                      dependency_job=phase_dependency)
        time.sleep(5)
        job_tracker["methyl_phase"] = pipeline.methylation_hap(tumor=tumor_data, control=use_sample,
                                                               dependency_job=job_tracker["snv_phase"])
        time.sleep(5)

    elif "tumor" == sample_type and "tumor" == user_args.subcommand:
        # Tumor only CNV calling
        cnv_dependency = ",".join([job_id for job_id in [job_tracker["mosdepth"], job_tracker["snv_merge"]]
                                   if job_id != ""])
        job_tracker["spectre"] = pipeline.cnv(sample=pipeline.files_use.tumor, dependency_job=cnv_dependency)
    else:
        pass
    # wrapping up
    time.sleep(2)
    preprocess_logger.info(f'All jobs sent: {job_tracker}')
    return job_tracker


def loss_submit(loss_param: LossParameters, tumor_json: str, control_json: str, tumor_jobs: dict, control_jobs: dict,
                analysis_dir: str, job_id: str):
    analysis_logger = setup_log(__name__, loss_param.as_dev)
    analysis_logger.info(f'{Col.bold}{Col.Fg.cyan}TumorLens loss analysis: {Col.end}'
                         f'{Col.bold}{Col.Fg.red}TODO{Col.end}')
    loss_analysis_job = TumorLensLoss(loss_param, job_id)
    # TODO: submit via sbatch to be able to use dependencies flow control
    dependencies = [tumor_jobs[tj] for tj in tumor_jobs] + [control_jobs[cj] for cj in control_jobs]
    loss_analysis_job.submit_loss(tumor_json, control_json, dependencies, analysis_dir)


def make_loss_script(use_args: SomaticParameters, pipeline_data: TumorLensProcessing) -> tuple[LossParameters, str, str, str]:
    # pipeline_data.files_use.work_dir
    loss_param = LossParameters()
    # from user arguments
    loss_param.command_string = pipeline_data.bin_paths.schlons_py
    loss_param.subcommand = "loss"
    loss_param.work_dir = use_args.work_dir
    loss_param.fasta = use_args.reference
    loss_param.as_dev = use_args.as_dev
    loss_param.dry_run = use_args.dry_run
    loss_param.output_prefix = f'{use_args.sample_name}_pair'
    # from piepline data
    # loss_param.annot = pipeline_data.bin_paths.annotation.hla_bed_pad_gz
    # loss_param.hla_rare = pipeline_data.bin_paths.annotation.hla_ciwd
    # loss_param.hla_annot = pipeline_data.bin_paths.annotation.directory
    # analysis
    analysis_dir = f'{loss_param.work_dir}/analysis'
    if not os.path.exists(f'{analysis_dir}'):
        os.mkdir(analysis_dir)
    # tumor
    tumor_json = open(f'{analysis_dir}/tumor.json', 'w')
    tumor_json.write(f'''{{
        "name": "tumor",
        "sv": "{pipeline_data.files_use.tumor.dir_sv}/{pipeline_data.files_use.tumor.file_sv}",
        "snf": "{pipeline_data.files_use.tumor.dir_sv}/{pipeline_data.files_use.tumor.file_snf}",
        "sv_mosaic": "{pipeline_data.files_use.tumor.dir_sv}/{pipeline_data.files_use.tumor.file_sv_mosaic}",
        "snf_mosaic": "{pipeline_data.files_use.tumor.dir_sv}/{pipeline_data.files_use.tumor.file_snf_mosaic}",
        "snv": "{pipeline_data.files_use.tumor.dir_snv}/{pipeline_data.files_use.tumor.file_snv}",
        "cnv": "{pipeline_data.files_use.tumor.dir_cnv}/{pipeline_data.files_use.tumor.file_cnv}",
        "spc": "{pipeline_data.files_use.tumor.dir_cnv}/{pipeline_data.files_use.tumor.file_spc}",
        "coverage": "{pipeline_data.files_use.tumor.dir_cov}/{pipeline_data.files_use.tumor.file_cov}",
        "methyl": "{pipeline_data.files_use.tumor.dir_methyl}/{pipeline_data.files_use.tumor.file_methyl}",
        "typing": "{pipeline_data.files_use.tumor.dir_type}/{pipeline_data.files_use.tumor.file_type}",
        "methyl_hap": "{pipeline_data.files_use.tumor.dir_type}/{pipeline_data.files_use.tumor.file_me_phase_dmr_res}"
    }}\n''')
    tumor_json.close()
    # control
    control_json = open(f'{analysis_dir}/control.json', 'w')
    control_json.write(f'''{{
        "name": "control",
        "sv": "{pipeline_data.files_use.control.dir_sv}/{pipeline_data.files_use.control.file_sv}",
        "snf": "{pipeline_data.files_use.control.dir_sv}/{pipeline_data.files_use.control.file_snf}",
        "sv_mosaic": "{pipeline_data.files_use.control.dir_sv}/{pipeline_data.files_use.control.file_sv_mosaic}",
        "snf_mosaic": "{pipeline_data.files_use.control.dir_sv}/{pipeline_data.files_use.control.file_snf_mosaic}",
        "snv": "{pipeline_data.files_use.control.dir_snv}/{pipeline_data.files_use.control.file_snv}",
        "cnv": "{pipeline_data.files_use.control.dir_cnv}/{pipeline_data.files_use.control.file_cnv}",
        "spc": "{pipeline_data.files_use.control.dir_cnv}/{pipeline_data.files_use.control.file_spc}",
        "coverage": "{pipeline_data.files_use.control.dir_cov}/{pipeline_data.files_use.control.file_cov}",
        "methyl": "{pipeline_data.files_use.tumor.dir_me_phase}/{pipeline_data.files_use.control.file_methyl}",
        "typing": "{pipeline_data.files_use.control.dir_type}/{pipeline_data.files_use.control.file_type}",
        "methyl_hap": "Single file"
    }}\n''')
    control_json.close()

    return loss_param, f'{analysis_dir}/tumor.json', f'{analysis_dir}/control.json', analysis_dir


def main():
    get_user_args = GetArguments()
    args, main_help = get_user_args.tumor_lens_args()
    if not args.subcommand:
        print(main_help)
    else:
        if "somatic" == args.subcommand:
            tumor_lens(args)
        elif "tumor" == args.subcommand:
            tumor_lens(args)
        elif "loss" == args.subcommand:
            loss(args, script_path)
        else:
            print(main_help)
            sys.exit(1)


# main
if __name__ == '__main__':
    main()
