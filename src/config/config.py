#!/usr/bin/env python3
import sys
import argparse
from utils.term_colors import TermColors as Col
from .parameters import ConfigParameters

VERSION = "1.0"
PROCESSING = "1.0"
LOSS = "1.0"


def preprocess_version():
    print(f'{Col.bold}{Col.Fg.green}Pre-process{Col.end} version: {PROCESSING}')


def loss_version():
    print(f'{Col.bold}{Col.Fg.cyan}Analysis loss{Col.end} version: {LOSS}')


def tumorlens_version():
    print(f'{Col.bold}{Col.Fg.blue}TumorLens{Col.end} version: {VERSION}')


def version_print(ver = ""):
    tumorlens_version()
    preprocess_version()
    loss_version()


class GetArguments(object):
    def __init__(self):
        self.tumor_lens_help = (f'{Col.bold}TumorLens{Col.end}, a systematic somatic-variant analysis of cancer sample.'
                                f'\nThe pipline can be run in two modes: {Col.bold}{Col.Fg.cyan}somatic{Col.end} for '
                                f'tumor/normal comparison and {Col.bold}{Col.Fg.green}tumor{Col.end} for tumor-only')

    def tumor_lens_args(self):

        tl_help = self.tumor_lens_help
        
        parser = argparse.ArgumentParser(
                 description=f"{Col.bold}TumorLens{Col.end}: Systematic somatic-variant analysis of cancer sample.",
                 usage=tl_help
        )
        subparsers = parser.add_subparsers(help=tl_help, dest="command")

        # ############################################################################################################ #
        somatic_help = "Somatic analysis"
        somatic = subparsers.add_parser("somatic", help=somatic_help)
        # Tumor
        somatic.add_argument('-1', '--tumor', type=str, required=True, dest='tumor_in',
                             default="", help='Tumor sample input file (fastq or bam), long reads')
        somatic.add_argument('-s', '--sample-name', type=str, required=True, dest='sample_name',
                             default="", help='Sample name')
        somatic.add_argument('-2', '--control', type=str, required=True, dest='control_in',
                             default="", help='Control sample input file (fastq or bam), long reads')
        somatic.add_argument('-d', '--work-dir', type=str, required=True, dest='work_dir', default="",
                             help='Directory where the analysis is performed, can be created')
        somatic.add_argument('-r', '--reference', type=str, required=True, dest='reference', default="",
                             help=f'Human reference genome in fasta. GRCh38 is needed as the used annotation are in '
                                  f'GRCh38 coordinates')
        # Other
        somatic.add_argument('-p', '--tumor-purity', type=float, required=False,
                             dest='tumor_purity', default=1.0, help='Tumor purity')
        somatic.add_argument('-g', '--guppy-ver', type=int, required=False, dest='guppy_ver', default=6,
                             help='Guppy version for Clair3, default=6')
        somatic.add_argument('-!', '--debug', action='store_true', required=False, dest='as_dev',
                             default=False, help='')
        somatic.add_argument('-0', '--dry-run', action='store_true', required=False, dest='dry_run',
                             default=False, help='')

        # ############################################################################################################ #
        tumor_help = "Tumor only analysis"
        tumor = subparsers.add_parser("tumor", help=tumor_help)
        # from fastq/modbam dir to all files for analysis
        tumor.add_argument('-1', '--tumor', type=str, required=True, dest='tumor_in', default="",
                           help='Tumor sample input file (fastq or bam), long reads')
        tumor.add_argument('-s', '--sample-name', type=str, required=True, dest='sample_name',
                           default="", help='Prefix for all results files')
        tumor.add_argument('-p', '--tumor-purity', type=float, required=False,  dest='tumor_purity', 
                           default=1.0, help='Tumor purity')
        tumor.add_argument('-d', '--work-dir', type=str, required=True, dest='work_dir', default="",
                           help='Directory where the analysis is performed, can be created')
        tumor.add_argument('-r', '--reference', type=str, required=True, dest='reference', default="",
                           help=f'Human reference genome in fasta. GRCh38 is needed as the used annotation are in '
                                f'GRCh38 coordinates')
        tumor.add_argument('-g', '--guppy-ver', type=int, required=False, dest='guppy_ver', default=6,
                           help='Guppy version for Clair3, default=6')
        tumor.add_argument('-!', '--debug', action='store_true', required=False, dest='as_dev',
                           default=False, help='')
        tumor.add_argument('-0', '--dry-run', action='store_true', required=False, dest='dry_run',
                           default=False, help='')

        # ############################################################################################################ #
        # LoH pipeline for paired (tumor-control) samples
        loss_help = "Perform somatic analysis and HLA loss"
        loss = subparsers.add_parser("loss", help=loss_help)
        # basic arguments
        loss.add_argument('-d', '--work-dir', type=str, required=True, dest='work_dir',
                          default="", help='Directory where the analysis is performed, must exist')
        loss.add_argument('-f', '--reference', type=str, required=False, dest='reference', default="",
                          help='fasta reference')
        loss.add_argument('-t', '--tumor-json', type=str, required=True, dest='tumor_json', default="",
                          help='JSON file with required files/parameters')
        loss.add_argument('-n', '--control-json', type=str, required=True, dest='control_json', default="",
                          help='JSON file with required files/parameters')
        loss.add_argument('-o', '--output', type=str, required=True, dest='output_prefix', default="paired",
                          help='Output prefix for result files')
        # dev
        loss.add_argument('-!', '--debug', action='store_true', required=False, dest='as_dev',
                          default=False, help='as dev')
        loss.add_argument('-0', '--dry-run', action='store_true', required=False, dest='dry_run',
                          default=False, help='')

        # ############################################################################################################ #
        # Version
        show_version = subparsers.add_parser("version")

        # ############################################################################################################ #
        args = parser.parse_args()
        full_command = " ".join(sys.argv)
        if "version" == args.command:
            version_print()
            sys.exit(0)
        elif "somatic" == args.command:
            return ConfigParameters.somatic(args, full_command), tl_help
        elif "tumor" == args.command:
            return ConfigParameters.tumor(args, full_command), tl_help
        elif "loss" == args.command:
            return ConfigParameters.loss(args, full_command), tl_help
        else:
            print(tl_help)
            sys.exit(1)
