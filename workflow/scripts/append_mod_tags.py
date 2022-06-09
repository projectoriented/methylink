#!/usr/bin/env python3

import pysam
import sys

# LOGGING
sys.stdout = open(snakemake.log[0], "w")


def fetch_modified_bases(modified_obj) -> dict:
    # modified_obj: assumes a bam pysam obj with just modified bases in optional tags e.g. Mm & Ml
    tags_dict = {}
    for read in modified_obj.fetch(until_eof=True):
        if read.has_tag("Mm"):
            tags = read.get_tags()
            qname = read.query_name
            tags_dict[qname] = tags
    modified_obj.close()
    print(f"Base modification tags fetched for {modified_obj.filename.decode()}")
    return tags_dict


def write_linked_tags(bam, tags_dict, out_file) -> None:
    # bam: equivalent aligned bam
    # dict_tags: {query_name: [Mm tags and possibly Ml]}
    appended_tags = pysam.AlignmentFile(out_file, "wb", template=bam)
    for read in bam.fetch():
        if read.query_name in tags_dict.keys():
            read.set_tags(read.get_tags() + tags_dict[read.query_name])
        appended_tags.write(read)
    print(f"File written to: {out_file}")
    appended_tags.close()

    # write index
    pysam.index(out_file)
    print(f"Index written for {out_file}.bai")


def collect_tags(methyl_sn_input: list) -> dict:
    # methyl_sn_input: snakemake input
    tags = {}
    if not len(methyl_sn_input) == 1:
        for bam in methyl_sn_input:
            methyl_bam = pysam.AlignmentFile(bam, "rb", check_sq=False)
            dict_of_tags_per_bam = fetch_modified_bases(methyl_bam)
            tags.update(dict_of_tags_per_bam)
    else:
        methyl_bam = pysam.AlignmentFile(methyl_sn_input, "rb", check_sq=False)
        dict_of_tags_per_bam = fetch_modified_bases(methyl_bam)
        tags.update(dict_of_tags_per_bam)
    return tags


aln_bam = pysam.AlignmentFile(*snakemake.input.aln_bam, "rb")
tags_dict = collect_tags(snakemake.input.methyl_bam)
output_file = snakemake.output.linked_bam

write_linked_tags(aln_bam, tags_dict, output_file)
