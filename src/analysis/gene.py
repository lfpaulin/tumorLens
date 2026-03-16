import numpy as np
from utils.logger_debug import setup_log
from utils import parse_vcf


logger = setup_log(__name__, debug=True)
GTs = {(0, 0): "HOM", (0, 1): "HET", (1, 0): "HET", (1, 1): "ALT", (None, None): "NA"}


class GeneResults(object):
    def __init__(self, gene_info):
        contig, start, end, gname, gid, strand = gene_info
        self.chromosome = contig
        self.start = int(start)
        self.end = int(end)
        self.region = f'{self.chromosome}:{self.start}-{self.end}'
        self.position = [self.start, self.end]
        self.name = gname
        self.strand = strand
        self.id = gid
        # also add results
        self.snv = None
        self.loh = None
        self.cov = None
        self.cnv = None
        self.sv = None
        self.meth = None
        self.type = ".\t.\t.\t."
        self.type_t = "NA"
        self.type_c = "NA"
        self.type_map_t = "NA"
        self.type_map_c = "NA"

    def __repr__(self):
        return (f"Gene: {self.name}\nSNV: {self.snv}\nLOH: {self.loh}\nCOV: {self.cov}\n"
                f"CNV: {self.cnv}\nSV: {self.sv}\nMeth: {self.meth}\nType: {self.type}\n"
                f"Tumor: {self.type_t}\nControl: {self.type_c}")

    def set_results(self, gene_results):
        self.snv = gene_results["snv"]
        self.loh = gene_results["loh"]
        self.cov = gene_results["coverage"]
        self.cnv = gene_results["cnv"]
        self.sv = gene_results["sv"]
        self.meth = gene_results["methylation"]

    def print_output(self):
        # header_out = "\tContig\tMethylation 5mC%\tSV 50+\tCNV 100kb+\tLoH 1Mb+\tTypingx4\n"
        sv: parse_vcf.SV
        sv_out = "|".join([f'{sv.svtype},{sv.svlen},{GTs[sv.tumor]}' for sv in self.sv])
        cnv: parse_vcf.SV
        cnv_out = "|".join([f'{cnv.svtype},{cnv.svlen},{cnv.tumor}' for cnv in self.cnv])
        output_line = "\t".join([self.name, self.region, sv_out, cnv_out, self.type, self.meth])
        return f'{output_line}\n'

    def set_hla(self, type_class, hla_t, hla_c, map_type_class, hla_mt, hla_mc):
        self.type_t = hla_t
        self.type_c = hla_c
        self.type_map_t = hla_mt
        self.type_map_c = hla_mc
        # typing classification
        if type_class == "LOH_T_N":
            hla_class_t, hla_class_c = "HOM", "HOM"
        elif type_class == "LOH_T":
            hla_class_t, hla_class_c = "LOH", "HET"
        elif type_class == "LOH_N":
            hla_class_t, hla_class_c = "HET", "HOM"
        else:
            hla_class_t, hla_class_c = "HET", "HET"
        # typing-mapped classification
        if map_type_class == "LOH:2":
            hla_class_mt, hla_class_mc = "HOM", "HOM"
        elif map_type_class == "LOH:1":
            hla_class_mt, hla_class_mc = "LOH", "HET"
        elif map_type_class == "LOH:0":
            hla_class_mt, hla_class_mc = "HET", "LOH"
        elif map_type_class == "loh:2":
            hla_class_mt, hla_class_mc = "HET/loh", "HET/loh"
        elif map_type_class == "loh:1":
            hla_class_mt, hla_class_mc = "HET/loh", "HET"
        elif map_type_class == "loh:0":
            hla_class_mt, hla_class_mc = "HET", "HET/loh"
        else:
            hla_class_mt, hla_class_mc = "HET", "HET"
        self.type = (f'{hla_class_t}:{self.type_t}\t{hla_class_mt}:{self.type_map_t}\t'
                     f'{hla_class_c}:{self.type_c}\t{hla_class_mc}:{self.type_map_c}')
