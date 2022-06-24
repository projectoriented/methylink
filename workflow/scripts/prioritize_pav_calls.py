#!/usr/bin/env python3

import pybedtools
import pandas as pd


def filter_sizes(bed_file, size=10000):
    df = pd.read_csv(bed_file, header=True, keep_default_na=False)
    return df.loc[(df['SVLEN'] >= size)]


def filter_bed(query_bed, lc_bed, gencode_bed, out_file):
    query = pybedtools.BedTool(query_bed)
    lc = pybedtools.BedTool(lc_bed)
    gencode = pybedtools.BedTool(gencode_bed)

    # remove low complexity regions
    lc_removed = query.intersect(lc, v=True)
    # filter for coding regions - this is not strand aware
    return lc_removed.intersect(gencode, u=True, output=out_file)


filter_bed(filter_sizes(snakemake.input.pav_ins_bed), snakemake.input.lc_bed, snakemake.input.gencode_bed, snakemake.output.out_file)

