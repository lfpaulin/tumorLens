#!/usr/bin/env python


class SomaticParameters:
    # info
    command_string: str = ""
    subcommand: str = ""
    # params
    tumor: str = ""
    sample_name: str = ""
    tumor_purity: float = 1.0
    work_dir: str = ""
    guppy_ver: int = 6
    reference: str = ""
    as_dev: bool = False
    dry_run: bool = False
    control: str = ""
    # extras
    tumor_file: str
    control_file: str


class LossParameters:
    # info
    command_string: str = ""
    subcommand: str = ""
    # params
    work_dir: str = ""
    reference: str = ""
    tumor_json: str = ""
    control_json: str = ""
    output_prefix: str = ""
    as_dev: bool = False
    dry_run: bool = False


class NoneParameters:
    subcommand: str = ""


class ConfigParameters:
    @staticmethod
    def somatic(user_args, full_command_str) -> SomaticParameters:
        args = SomaticParameters()
        # info
        args.command_string = full_command_str
        args.subcommand = user_args.command
        # params
        args.tumor = user_args.tumor_in
        args.sample_name = user_args.sample_name
        args.tumor_purity = user_args.tumor_purity
        args.work_dir = user_args.work_dir
        args.reference = user_args.reference
        args.guppy_ver = user_args.guppy_ver
        args.as_dev = user_args.as_dev
        args.dry_run = user_args.dry_run
        args.control = user_args.control_in
        return args

    @staticmethod
    def tumor(user_args, full_command_str) -> SomaticParameters:
        args = SomaticParameters()
        # info
        args.command_string = full_command_str
        args.subcommand = user_args.command
        # params
        args.tumor = user_args.tumor_in
        args.sample_name = user_args.sample_name
        args.tumor_purity = user_args.tumor_purity
        args.work_dir = user_args.work_dir
        args.reference = user_args.reference
        args.guppy_ver = user_args.guppy_ver
        args.as_dev = user_args.as_dev
        args.dry_run = user_args.dry_run
        return args

    @staticmethod
    def loss(user_args, full_command_str) -> LossParameters:
        args = LossParameters()
        # info
        args.command_string = full_command_str
        args.subcommand = user_args.command
        # params
        args.work_dir = user_args.work_dir
        args.reference = user_args.reference
        args.tumor_json = user_args.tumor_json
        args.control_json = user_args.control_json
        args.output_prefix = user_args.output_prefix
        args.as_dev = user_args.as_dev
        args.dry_run = user_args.dry_run
        return args

    @staticmethod
    def none(user_args_str) -> NoneParameters:
        params_none = NoneParameters()
        params_none.subcommand = user_args_str
        return params_none
