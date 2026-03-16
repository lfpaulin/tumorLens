#!/usr/bin/env python3
import argparse
import pysam
import csv
from logger_debug import setup_log


logger = setup_log(__name__, True)


# source https://github.com/ANHIG/IMGTHLA/
# uses git lfs
# Version
def the_version():
    print(f'hlatyping_fasta.py version 0.1\nIMGTHLA version 3.51 26.01.2023')


class HLATYPEGeneric(object):
    # init
    def __init__(self, input_line, digits=8):
        # the results file has 11 columns, only number 3 and 4 are needed
        #   Locus
        #   Chromosome
        #   Allele *needed
        #   Q1     *needed
        #   Q2 | AverageCoverage | CoverageFirstDecile | MinimumCoverage | proportionkMersCovered
        #   LocusAvgColumnError | NColumns_UnaccountedAllele_fGT0.2 | (perfectG)
        tab_sep_fields = input_line.split("\t")
        self.entry = input_line
        self.MIN_QUAL = 0.5  # coin toss
        self.allele = ""
        # extract needed fields
        needed_fields_index = 4
        self.hla, self.chromosome, allele, qual_score = tab_sep_fields[:needed_fields_index]
        self.qual_score = float(qual_score)
        self.digits = digits
        self.hla_use_n_digits(allele)

    def hla_use_n_digits(self, allele):
        # examples
        #   2digit  A*01          -> [:1] 2/2
        #   4digit  A*01:01       -> [:2] 4/2
        #   6digit  A*01:01:01    -> [:3] 6/2
        #   8digit  A*01:01:01:01 -> [:4] 8/2
        digit_to_index = int(self.digits/2)
        self.allele = ":".join(allele.split(":")[:digit_to_index])


class HLATYPEBest(object):
    # init
    def __init__(self, input_line, digits=8):
        # the results file has 11 columns, only number 3 and 4 are needed
        #   Locus
        #   Chromosome
        #   Allele *needed
        #   Q1     *needed
        #   Q2
        #   AverageCoverage
        #   CoverageFirstDecile
        #   MinimumCoverage
        #   proportionkMersCovered
        #   LocusAvgColumnError
        #   NColumns_UnaccountedAllele_fGT0.2
        # split and analysis
        self.N_FIELDS = 11
        tab_sep_fields = input_line.split("\t")
        self.ERROR = len(tab_sep_fields) != self.N_FIELDS
        self.MIN_QUAL = 0.5  # coin toss
        if not self.ERROR:
            # extract needed fields
            field_chromosome_index = 1  #2
            field_allele_index = 2  #3
            field_q1_index = 3      #4
            self.chromosome = tab_sep_fields[field_chromosome_index]
            self.qual_score = float(tab_sep_fields[field_q1_index])
            self.allele = hla_use_n_digits(tab_sep_fields[field_allele_index], digits) if self.qual_score > self.MIN_QUAL else ""


class HLATYPEGBest(object):
    # init
    def __init__(self, input_line):
        # the results file has 12 columns, only number 3 and 4 are needed
        #   Locus
        #   Chromosome
        #   Allele *needed
        #   Q1     *needed
        #   Q2
        #   AverageCoverage
        #   CoverageFirstDecile
        #   MinimumCoverage
        #   proportionkMersCovered
        #   LocusAvgColumnError
        #   NColumns_UnaccountedAllele_fGT0.2
        #   perfectG
        # split and analysis
        self.N_FIELDS = 12
        tab_sep_fields = input_line.split("\t")
        self.ERROR = len(tab_sep_fields) != self.N_FIELDS
        self.MIN_QUAL = 0.5  # coin toss
        if not self.ERROR:
            # extract needed fields
            field_chromosome_index = 1  #2
            field_allele_index = 2      #3
            field_q1_index = 3          #4
            self.chromosome = tab_sep_fields[field_chromosome_index]
            self.qual_score = float(tab_sep_fields[field_q1_index])
            self.allele = tab_sep_fields[field_allele_index]  # Not used: if self.qual_score > self.MIN_QUAL else ""


class HLACIWD(object):
    # init
    def __init__(self):
        # need the HLA G group and classification, ID is optional, but if it has it will be added
        # 3 columns
        #   G group
        #   ID
        #   classification:
        #       C  = common,          (≥1 in 10 000)
        #       I  = intermediate     (≥1 in 100 000)
        #       WD = Well-Documented  (≥5 occurrences)
        #       not-CIWD otherwise
        # split and analysis
        self.ciwd_by_hla_g_group = {}
        self.id_by_hla_g_group = {}
        self.hla_genes_in_list = ["A", "B", "C", "DRB1", "DQB1", "DPB1"]

    def add_ciwd_allele(self, line):
        line_to_list = line.split(",")
        if len(line_to_list) != 3:
            logger.warning("Could not be added")
        [g_group, id, cwid] = line_to_list
        if id != "":
            self.id_by_hla_g_group[g_group] = id
        self.ciwd_by_hla_g_group[g_group] = cwid

    def get_ciwd_allele(self, g_group):
        if g_group in self.id_by_hla_g_group:
            return self.ciwd_by_hla_g_group[g_group], self.id_by_hla_g_group[g_group]
        else:
            return self.ciwd_by_hla_g_group[g_group], None

    def is_hla_rare(self, hla_g_group):
        # will only give answer for genes in the list, if not assume not rare
        hla_gene = hla_g_group.split("*")[0]
        if hla_gene in self.hla_genes_in_list:
            return hla_g_group not in self.ciwd_by_hla_g_group


# from hla-names (csv) to dict
def get_hla_names(hla_csv, digit=2):
    hla_names = open(hla_csv, newline='')
    hla_names_csv_reader = csv.reader(hla_names, delimiter=",")
    hla_dict = {}
    for row in hla_names_csv_reader:
        if "#" not in row:
            [hla_type, hla_name] = row
            hla_type_by_digit = hla_use_n_digits(hla_type, digit)
            if hla_type_by_digit not in hla_dict:
                hla_dict[hla_type_by_digit] = [hla_name]
            else:
                hla_dict[hla_type_by_digit].append(hla_name)  # should only occur when using 8 digits(?)
    return hla_dict


def hla_use_n_digits(hla_type, digit=8):
    # examples
    #   2digit  A*01          -> [:1] 2/2
    #   4digit  A*01:01       -> [:2] 4/2
    #   6digit  A*01:01:01    -> [:3] 6/2
    #   8digit  A*01:01:01:01 -> [:4] 8/2
    digit_to_index = int(digit/2)
    return ":".join(hla_type.split(":")[:digit_to_index])


# only if G groups is not available
def get_hla_types(hla_type_file, digits=2):
    file_handler = open(hla_type_file, "r")
    hla_la = []        # final result
    hla_la_tmp = []    # only the alleles is counted
    for line in file_handler:
        if "Allele" not in line:
            hlatype_entry = HLATYPEGeneric(line.rstrip("\n"), digits)
            logger.debug(hlatype_entry.allele)
            if hlatype_entry.allele not in hla_la_tmp:
                hla_la.append(f'{hlatype_entry.allele}|{hlatype_entry.chromosome}')
                hla_la_tmp.append(hlatype_entry.allele)
    return hla_la


# if available, use
def make_hla_ciwd_obj(ciwd_file):
    file_handler = open(ciwd_file, "r")
    hla_ciwd_class = HLACIWD()
    for line in file_handler:
        if not line.startswith("#"):
            hla_ciwd_class.add_ciwd_allele(line.rstrip("\n"))
    return hla_ciwd_class


def hlatype_to_fasta(args):
    # args:
    #   'fasta'    big fasta file with sequences for each type (8 digits)
    #   'hlatype'    the results from HLA-LA can be R1_bestguess_G.txt or R1_bestguess.txt it uses the first result
    #              and only n digits, default digits=2
    #   'hla_csv'  HLA type to ID in the DB/fasta file, make dict
    # hla type to ID. Mandatory, used to fetch the fasta seq
    hla_names = get_hla_names(args.hla_csv, args.hla_digit)
    # Optional: get the hla rarity classification (Optional)
    hla_ciwd = None if args.hla_ciwd == "" else make_hla_ciwd_obj(args.hla_ciwd)
    # get the hla-la results
    hla_types = get_hla_types(args.hlatype, args.hla_digit)
    logger.debug(hla_types)
    # read the fasta file wilt hla genes by type, mandatory for fasta seq retrieval
    fas_genomic, fas_cds, fasta_cds = "", "", None
    if "," in args.fasta:
        fasta_list = args.fasta.split(",")
        if len(fasta_list) == 2:
            [fas_genomic, fas_cds] = fasta_list
        else:
            [fas_genomic, fas_cds] = fasta_list[:2]
        fasta_cds = pysam.FastaFile(fas_cds)
    else:
        fas_genomic = args.fasta
    fasta = pysam.FastaFile(fas_genomic)
    # use this genes
    # hla_genes_use = {"A*", "B*", "C*", "DQA1*", "DQB1*", "DRB1*", "DPA1*", "DPB1*"}
    hla_genes_use = ["A", "B", "C", "DQA1", "DQB1", "DRB1", "DPA1", "DPB1"]
    # get the sequences that are in the hla-la results
    fasta_seq = ""  # work var
    for each_type in hla_types:
        [hla_type, this_chr] = each_type.split("|")
        hla_type_no_start = hla_type.split("*")[0]
        this_hla = None
        if hla_type in hla_names and hla_type_no_start in hla_genes_use:
            [this_hla] = hla_names[hla_type]
            if this_hla is not None:
                try:
                    my_reference = f'HLA:{this_hla}'
                    fasta_seq = fasta.fetch(reference=my_reference)
                    logger.debug(f'[*] HAS-SEQ  {hla_type} | {this_chr} -> {this_hla}')
                    print(f'>{this_hla}__{pretty_hlatype(hla_type)}__allele{this_chr}', end="\n")
                    print(f'{fasta_seq}', end="\n")
                except KeyError:
                    logger.warning(f'[!] NO-FAS   {hla_type} | {this_chr} -> {this_hla}')
                    logger.warning(f'[!] NO-FAS   using CDS instead ...')
            else:
                logger.debug(f'[!] NO-NAME  {hla_type} | {this_chr} -> {this_hla}')
    # add the masked hla locus 26Mb - 36Mb
    fasta_hla_masked = pysam.FastaFile(args.hla_locus_masked)
    [hla_mask] = fasta_hla_masked.references
    fasta_hla_masked_seq = fasta_hla_masked.fetch(reference=hla_mask)
    print(f'>{hla_mask}', end="\n")
    print(f'{fasta_hla_masked_seq}', end="\n")


def pretty_hlatype(hla_type_digits):
    if ":" in hla_type_digits:
        hla_type_digits = "_".join(hla_type_digits.split(":"))
    if "*" in hla_type_digits:
        hla_type_digits = "_".join(hla_type_digits.split("*"))
    return hla_type_digits


# NOTE: not used
# from hla-groups (csv) to dict
def get_hla_groups(hla_csv):
    hla_groups = open(hla_csv, newline='')
    hla_groups_reader = csv.reader(hla_groups, delimiter=",")
    hla_dict = dict()
    for row in hla_groups_reader:
        # [key, val]
        [hla_group, hla_types] = row
        if "/" in hla_types:
            """
            From Chris Hammer:
            Yes, I guess that's an issue that comes up due to our HLA typing tool choice.
            We used HLA-HD internally for ASHLI development, which provides 6-digit resolution, not G group.
            That said, it shouldn't be a big issue. G groups are identical in their genomic sequence in the
            exons of HLA genes that define the peptide binding domain, and usually don't have a lot of additional
            variation. Thinking pragmatically, I'd say just use the first representative of a g group, which very
            often will be the one very common allele in that group, among many very rare ones.
            # will only use the first one
            """
            only_first = 0  # index of list
            use_first_hlatype = hla_types.split("/")[only_first]
            hla_gene = hla_group.split("*")[only_first]
            hla_dict[hla_group] = f'{hla_gene}*{use_first_hlatype}'  # only use the first one in the list
        else:
            hla_dict[hla_group] = [hla_types]  # they are the same
    return hla_dict


# preferably use this
def get_hla_types_g_groups(hla_type_g_group_file):
    file_handler = open(hla_type_g_group_file, "r")
    hla_la = []
    hla_la_tmp = []
    for line in file_handler:
        if "Allele" in line:
            pass  # skip header
        else:
            hlatype_entry = HLATYPEGBest(line.rstrip("\n"))
            # if hlatype_entry.qual_score > hlatype_entry.MIN_QUAL and hlatype_entry.allele not in hla_la:
            if hlatype_entry.allele not in hla_la_tmp:
                hla_la.append(f'{hlatype_entry.allele}|{hlatype_entry.chromosome}')
                hla_la_tmp.append(hlatype_entry.allele)
    return hla_la


def get_arguments():
    hlaextract_help = """
    hlaextract <command> [<args>]
        # Extract from two digit:
            hlafa    Extract fasta seq from HLA-typing result (HLA-LA)
                     Mandatory:  --fasta | --type | --hla-names | --hla-locus-masked
                     Optional:   --ciwd-groups

        # Version
            version    Shows current version/build
    """
    parser = argparse.ArgumentParser(
             description="HLA typing to fasta",
             usage=hlaextract_help
    )
    subparsers = parser.add_subparsers(help=hlaextract_help, dest="command")

    # ############################################################################################ #
    # Version
    version_help = "Gives the version number"
    subparser_version = subparsers.add_parser("version", help=version_help)

    # ############################################################################################ #
    # Help
    the_help = "Help"
    subparser_help = subparsers.add_parser("help", help=the_help)

    # ############################################################################################ #
    #  Extract from two digit:
    # Genotype SNV
    hla_fasta_help = "Extracts the fasta sequence from an HLA typing result (HLA-LA)"
    subparser_hla_fasta = subparsers.add_parser("hlafa", help=hla_fasta_help)

    subparser_hla_fasta.add_argument('-f', '--fasta', type=str, required=True, dest='fasta', default="",
                                     help='Fasta file from which the sequences will be extracted')
    subparser_hla_fasta.add_argument('-t', '--type', type=str, required=True, dest='hlatype', default="",
                                     help='HLA-LA typing result, only first result column will be used')
    subparser_hla_fasta.add_argument('-n', '--hla-names', type=str, required=True, dest='hla_csv', default="",
                                     help='HLA alleles file with the relationship between the HLA type and its "name"')
    subparser_hla_fasta.add_argument('-d', '--digit', type=int, required=False, dest='hla_digit', default=2,
                                     help='Number of digits to use from the typing results, Default = 2')
    subparser_hla_fasta.add_argument('-m', '--hla-locus-masked', type=str, required=True, dest='hla_locus_masked',
                                     default="hla_locus_26Mb_36Mb-genes_masked.fasta.gz", help='...')
    # optional
    subparser_hla_fasta.add_argument('-c', '--ciwd-groups', type=str, required=False, dest='hla_ciwd', default="",
                                     help='Common, Intermediate and Well-Documented HLA Alleles by G group')

    # ############################################################################################ #

    args = parser.parse_args()
    return args, hlaextract_help


def main():
    args, main_help = get_arguments()
    command = args.command
    logger.info(f'{command}')
    if command is None:
        print(main_help)
    # Version
    elif command == "version":
        the_version()
    elif command == "help":
        print(main_help)
    elif command == "hlafa":
        hlatype_to_fasta(args)
    else:
        pass


# main
if __name__ == '__main__':
    main()
