#!/usr/bin/env python3

import pysam

GTs = {(0, 0): "HOM", (0, 1): "HET", (1, 0): "HET", (1, 1): "ALT", (None, None): "NA"}


class SV(object):
    def __init__(self):
        self.contig: str
        self.pos: int
        self.end: int
        self.svtype: str = ""
        self.svlen: int = 0
        self.suppvec: str
        self.svid: str = ""
        self.tumor: tuple[int, int] = (0, 0)
        self.region: str = ""
        self.type: str = ""

    def __repr__(self):
        if self.svlen == 0:
            return "No SV"
        tumor = GTs[self.tumor]
        return f'{self.region}\t{self.svid}\t{self.svtype}\t{self.svlen}\t{tumor}\t{self.type}'


def sv_somatic(vcf: pysam.VariantFile, genes: dict[str, dict[str, str | list[SV] | bool]]) -> (
        tuple)[dict[str, dict[str, str | list[SV] | bool]], list[SV]]:
    sv_somatic_list = []
    if len(genes) > 0:
        for gene in genes:
            vcf.seek(0)
            for sv_entry in vcf.fetch(region=genes[gene]["region"]):
                sv_info, is_somatic = variant_info_somatic(sv_entry, "sv")
                if is_somatic:
                    sv_somatic_list.append(sv_info)
                    genes[gene]["sv"].append(sv_info)
                    genes[gene]["init"] = True
    else:
        for sv_entry in vcf.fetch():
            sv_info, is_somatic = variant_info_somatic(sv_entry, "sv")
            if is_somatic:
                sv_somatic_list.append(sv_info)
    return genes, sv_somatic_list


class CNV(object):
    def __init__(self):
        self.contig: str
        self.pos: int
        self.end: int
        self.svtype: str = ""
        self.svlen: int = 0
        self.svid: str = ""
        self.tumor: tuple[int, int] = (0, 0)
        self.region: str = ""
        self.type: str = ""

    def __repr__(self):
        if self.svlen == 0:
            return "No CNV"
        tumor = GTs[self.tumor]
        return f'{self.region}\t{self.svid}\t{self.svtype}\t{self.svlen}\t{tumor}\t{self.type}'


def cnv_somatic(vcf: pysam.VariantFile, genes: dict[str, dict[str, str | list[CNV] | bool]]) -> (
        tuple)[dict[str, dict[str, str | list[CNV] | bool]], list[CNV]]:
    sv_somatic_list = []
    if len(genes) > 0:
        for gene in genes:
            vcf.seek(0)
            for sv_entry in vcf.fetch(region=genes[gene]["region"]):
                sv_info, is_somatic = variant_info_somatic(sv_entry, "cnv")
                if is_somatic:
                    sv_somatic_list.append(sv_info)
                    genes[gene]["cnv"].append(sv_info)
                    genes[gene]["init"] = True
    else:
        for sv_entry in vcf.fetch():
            sv_info, is_somatic = variant_info_somatic(sv_entry, "cnv")
            if is_somatic:
                sv_somatic_list.append(sv_info)
    return genes, sv_somatic_list


def variant_info_somatic(variant: pysam.VariantRecord, var_type: str) -> tuple[SV | CNV, bool]:
    if "sv" == var_type:
        sv_info = SV()
        sv_samples = variant.samples.keys()
        tum: str
        nor: str
        if len(sv_samples) == 2:
            sv_info.type = "somatic"
            suppvec = variant.info.get("SUPP_VEC")
            if "10" != suppvec:
                sv_info.contig = variant.contig
                sv_info.pos = variant.pos
                sv_info.svtype = variant.info.get("SVTYPE")
                sv_info.svlen = variant.info.get("SVLEN")
                sv_end = variant.info.get("END")
                sv_info.end = sv_end if sv_end is not None else (sv_info.pos + abs(sv_info.svlen)
                                                                 if sv_info.svtype not in ["INS", "BND"]
                                                                 else sv_info.pos + 1)

                sv_info.svid = variant.id
                sv_info.suppvec = suppvec
                sv_info.region = f'{sv_info.contig}:{sv_info.pos}-{sv_info.end}'
                [tum, _] = variant.samples.keys()
                sv_info.tumor = variant.samples.get(tum).get("GT")
                return sv_info, True
            else:
                return sv_info, False
        else:
            sv_info.type = "tumor/single"
            return sv_info, False
    elif "cnv" == var_type:
        sv_info = CNV()
        [tum] = variant.samples.keys()
        sv_info.type = variant.info.get("SOMATIC")
        if sv_info.type is not None:
            sv_info.type = "somatic"
            sv_info.contig = variant.contig
            sv_info.pos = variant.pos
            sv_info.svtype = variant.info.get("SVTYPE")
            sv_info.svlen = variant.info.get("SVLEN")
            sv_end = variant.info.get("END")
            sv_info.end = sv_end if sv_end is not None else (sv_info.pos + abs(sv_info.svlen)
                                                             if sv_info.svtype not in ["INS", "BND"]
                                                             else sv_info.pos + 1)

            sv_info.svid = variant.id
            sv_info.region = f'{sv_info.contig}:{sv_info.pos}-{sv_info.end}'
            sv_info.tumor = variant.samples.get(tum).get("GT")
            return sv_info, True
        else:
            return sv_info, False
    else:
        return SV(), False
