import numpy as np


def load_results(analysis_files=None, chromosome=None, start=None, end=None, padding=0):
    analysis_results = {"snv": {}, "sv": {}, "cnv": {}, "loh": {}, "coverage": {}, "methylation": {}}
    start = start - padding
    end = end + padding
    region_fetch = f'{chromosome}:{start}-{end}'

    # Structural Variants (SV)
    # Note: analysis_files.handler_sv => pysam.VariantFile
    var_list = []
    genomic_iterator = analysis_files.handler_sv.fetch(region=region_fetch) \
        if type(analysis_files.handler_sv) is dict else analysis_files.handler_sv.fetch(region=region_fetch)
    for variant in genomic_iterator:
        var_list.append(variant)
    sv_summary = True if len(var_list) > 0 else False
    analysis_results["sv"] = {"summary": sv_summary, "full": var_list}

    # Single Nucleotide Variants (SNV)
    # Note: analysis_files.handler_snv => pysam.VariantFile
    # Note: This one is used for LoH inside Spectre (CNV)

    # Copy Number Variants (CNV)
    # Note: analysis_files.handler_cnv => pysam.VariantFile
    var_list, loh_list = [], []
    genomic_iterator = analysis_files.handler_cnv.fetch(region=region_fetch) \
        if type(analysis_files.handler_cnv) is dict else analysis_files.handler_cnv.fetch(region=region_fetch)
    for variant in genomic_iterator:
        if "LOH" in variant.alts[0]:
            loh_list.append(variant)
        else:
            var_list.append(variant)
    cnv_summary = True if len(var_list) > 0 else False
    loh_summary = True if len(loh_list) > 0 else False
    # CNV and LoH based on CNV
    analysis_results["cnv"] = {"summary": cnv_summary, "full": var_list}
    analysis_results["loh"] = {"summary": loh_summary, "full": loh_list}

    # Coverage (Used for CNV, also plot)
    # Note: analysis_files.handler_cov => pysam.TabixFile
    for_plot = {"pos": [], "cov": []}
    var_list = []
    genomic_iterator = analysis_files.handler_cov.fetch(region=region_fetch) \
        if type(analysis_files.handler_cov) is dict else analysis_files.handler_cov.fetch(region=region_fetch)
    for variant in genomic_iterator:
        [contig, start, end, cover] = variant.rstrip("\n").split("\t")
        cover_float = float(cover)
        var_list.append([contig, int(start), int(end), cover_float])
        for_plot["pos"].append(int(start))
        for_plot["cov"].append(cover_float)
    analysis_results["coverage"] = {"full": var_list, "summary": np.nanmean(for_plot["cov"])} \
        if len(var_list) > 0 else {"full": None, "summary": None}
    analysis_results["coverage"]["plot"] = for_plot

    # Methylation modification 5mC (proxy for gene expression)
    # Note: analysis_files.handler_methyl = None    # pysam.TabixFile
    # 1	contig	name of reference sequence from BAM header	str
    # 2	start position	0-based start position	int
    # 3	end position	0-based exclusive end position	int
    # 10	Nvalid_cov	See definitions above.	int
    # 11	fraction modified	Nmod / Nvalid_cov	float
    # 12	Nmod	See definitions above.	int
    # 13	Ncanonical	See definitions above.	int
    var_list = []
    for_plot = {"pos": [], "meth": []}
    min_coverage_use_meth = 10
    genomic_iterator = analysis_files.handler_methyl.fetch(region=region_fetch) \
        if type(analysis_files.handler_methyl) is dict else analysis_files.handler_methyl.fetch(region=region_fetch)
    for variant in genomic_iterator:
        [contig, start, _, _, _, _, _, _, _, mth_info] = variant.rstrip("\n").split("\t")
        [valid, mod_frac, _, _, _, _, _, _, _] = mth_info.split(" ")
        if int(valid) > min_coverage_use_meth:
            mod_frac_float = float(mod_frac)
            var_list.append([contig, int(start), mod_frac_float])
            for_plot["pos"].append(int(start))
            for_plot["meth"].append(mod_frac_float)
    analysis_results["methylation"] = {"full": var_list, "summary": np.nanmean(for_plot["meth"])} \
        if len(var_list) > 0 else {"full": None, "summary": None}
    analysis_results["methylation"]["plot"] = for_plot

    # return results
    return analysis_results


def zscore(some_np_array=None):
    if isinstance(some_np_array, np.ndarray):
        return (some_np_array - np.nanmean(some_np_array))/np.nanstd(some_np_array)
    else:
        return None
