import sys
import numpy as np
import pysam

from .gene import GeneResults
from utils.logger_debug import setup_log

logger = setup_log(__name__, debug=True)

HLA_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]
HLA_REGION = {"chromosome": "chr6", "start": 26000000, "end": 36000000, "coordinates": "chr6:29000000-34000000"}
GT_ALT = (1, 1)
GT_HET = (0, 1)
ACCEPTED_GT = [(1, 1), (0, 1)]


class HLALocusResults(object):
    def __init__(self):
        self.chromosome = HLA_REGION["chromosome"]
        self.start = int(HLA_REGION["start"])
        self.end = int(HLA_REGION["end"])
        self.position = [self.start, self.end]
        self.region_coordinates = f'{self.chromosome}:{self.start}-{self.end}'
        self.window_size = 1000
        # also add results
        self.loh = None
        self.cov = None
        self.cnv = None
        self.sv = None
        self.snv = None
        # unused
        self.meth = None

    def __repr__(self):
        return f'loh:{self.loh}, cov:{self.cov}, cnv:{self.cnv}, sv:{self.sv}, snv:{self.snv}'

    def set_coverage(self, file_handler, chr_region: str = ""):
        # file expected to have contig<TAB>start<TAB>end<TAB>coverage
        # go to start
        file_handler.seek(0)
        coverage, positions = [], []
        use_region = self.region_coordinates if "" == chr_region else chr_region
        for line in file_handler.fetch(region=use_region):
            _, pos, _, cover = line.rstrip("\n").split("\t")
            coverage.append(float(cover))
            positions.append(int(pos))
        self.cov = {"pos": positions, "cov": coverage}

    def set_cnv_loh(self, file_handler):
        # file expected to be VCF
        # go to start
        file_handler.seek(0)
        self.loh = []
        self.cnv = []
        min_loh_size = 1000000
        min_cnv_size = 100000
        for var_record in file_handler.fetch(region=self.region_coordinates):
            if var_record.info["SVLEN"] >= min_loh_size and var_record.info["SVTYPE"] == "LOH":
                self.loh.append(var_record)
            elif var_record.info["SVLEN"] >= min_cnv_size and var_record.info["SVTYPE"] != "LOH":
                self.cnv.append(var_record)
            else:
                pass

    def set_sv(self, file_handler):
        # file expected to be VCF
        # go to start
        file_handler.seek(0)
        # for visualization purposes only
        self.sv = []
        min_sv_size = 100000
        for var_record in file_handler.fetch(region=self.region_coordinates):
            if var_record.info["SVTYPE"] != "BND":
                if abs(var_record.info["SVLEN"]) >= min_sv_size:
                    self.sv.append(var_record)

    def set_snv(self, file_handler):
        # file expected to be VCF
        # go to start
        file_handler.seek(0)
        self.snv = {"pos": [], "vaf": [], "gt": []}
        # summarized by window size
        for start_region in range(self.start, self.end, self.window_size):
            region_chunk = f'{self.chromosome}:{start_region}-{start_region+self.window_size-1}'
            tmp_vaf, tmp_gt = [], []
            for var_record in file_handler.fetch(region=region_chunk):
                gt = var_record.samples.values()[0]["GT"]
                if gt in ACCEPTED_GT:
                    tmp_vaf.append(var_record.samples.values()[0]["AF"])
                    tmp_gt.append(gt)
            self.snv["pos"].append(start_region)
            self.snv["vaf"].append(np.nanmean(tmp_vaf) if len(tmp_vaf) else np.nan)
            self.snv["gt"].append(sum([1 if gt == GT_ALT else 0 for gt in tmp_gt])/len(tmp_gt)
                                  if len(tmp_gt) > 0 else np.nan)

    def set_methylation(self, file_handler: pysam.TabixFile, chr_region: str = ""):
        # file expected to have contig<TAB>start<TAB>end<TAB>coverage
        # go to start
        file_handler.seek(0)
        mod_frac, positions = [], []
        use_region = self.region_coordinates if "" == chr_region else chr_region
        pos_idx = 1
        mod_idx = 3
        frac_idx = 9
        if len(file_handler.contigs) == 0:
            self.meth = {"pos": [0], "mod": [0.0]}
        else:
            logger.debug(file_handler.contigs)
            for line in file_handler.fetch(region=use_region):
                line_split = line.rstrip("\n").split("\t")
                mod = line_split[mod_idx]
                if "m" == mod:
                    pos = line_split[pos_idx]
                    positions.append(int(pos))
                    if len(line_split) == 10:
                        # older version of modkit
                        _, mfrac, nmod, ncanon, _, _, _, _, _ = line_split[frac_idx].split(" ")
                    elif len(line_split) == 18:
                        # newer version of modkit (v0.5.0)
                        _, mfrac, nmod, ncanon, _, _, _, _, _ = line_split[frac_idx:]
                    else:
                        logger.error(f'{len(line_split)}, {line_split}')
                        sys.exit(1)
                    mod_frac.append(float(mfrac))
            self.meth = {"pos": positions, "mod": mod_frac}


def load_hla_annotation(bed_file, as_dict=False):
    # chr start end gene_name strand gene_id
    hla_annot = dict()
    gene_annot = dict()
    bed_file.seek(0)
    for line in bed_file:
        my_gene = GeneResults(line.rstrip("\n").split("\t"))
        if my_gene.name not in gene_annot:
            gene_annot[my_gene.name] = my_gene
            if "HLA" in my_gene.name:
                hla_annot[my_gene.name] = my_gene
        else:
            logger.warning(f'repeated gene {my_gene.name}', f'load_annot')
    if as_dict:
        return {"gene": gene_annot, "hla": hla_annot}
    return gene_annot, hla_annot
