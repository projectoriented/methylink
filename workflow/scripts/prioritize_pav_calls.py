#!/usr/bin/env python3

import pybedtools

PAV_HEADER = ["#CHROM", "POS", "END", "ID", "SVTYPE", "SVLEN", "HAP", "GT", "CLUSTER_MATCH", "CALL_SOURCE",
              "TIG_REGION", "QUERY_STRAND", "ALIGN_INDEX", "LEFT_SHIFT", "HOM_REF", "HOM_TIG", "HAP_VARIANTS", "CI",
              "HAP_RO", "HAP_OFFSET", "HAP_SZRO", "HAP_OFFSZ"]


# filter out low-complexity regions
def remove_low_complexity(insertion_bed, lc_bed):
    a = pybedtools.BedTool(insertion_bed)
    b = pybedtools.BedTool(lc_bed)
    return a.intersect(b, v=True)


# filter for coding regions
def get_cds(insertion_bed, gencode_bed):
    a = pybedtools.BedTool(insertion_bed)
    b = pybedtools.BedTool(gencode_bed)
    return a.intersect(b, u=True)


def get_sizes(filtered_bed_obj, size=10000):
    df = filtered_bed_obj.to_dataframe(header=None, keep_default_na=False, names=PAV_HEADER)
    return df.loc[(df['SVLEN'] >= size)]
