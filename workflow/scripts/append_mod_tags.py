#!/usr/bin/env python3

import pysam
import sys

# LOGGING
sys.stdout = open(snakemake.log[0], "w")

aln_bam = pysam.AlignmentFile(snakemake.input.aln_bam, "rb")
methyl_bam = pysam.AlignmentFile(snakemake.input.methyl_bam, "rb", check_sq=False)
output_file = snakemake.output.linked_bam


def fetch_modified_bases(modified_obj):
    # modified_obj: assumes a bam with just modified bases in optional tags e.g. Mm & Ml
    tags_dict = {}
    for read in modified_obj.fetch(until_eof=True):
        if read.has_tag("Mm"):
            tags = read.get_tags()
            qname = read.query_name
            tags_dict[qname] = tags
    modified_obj.close()
    print(f"Base modification tags fetched for {modified_obj.filename.decode()}")
    return tags_dict


def write_linked_tags(bam, tags_dict, out_file):
    # bam: equivalent aligned bam
    # dict_tags: {query_name: [Mm tags and possibly Ml]}
    appended_tags = pysam.AlignmentFile(out_file, "wb", template=bam)
    for read in bam.fetch(until_eof=True):
        if read.query_name in tags_dict.keys() and not read.is_unmapped:
            read.set_tags(read.get_tags() + tags_dict[read.query_name])
        appended_tags.write(read)
    print(f"File written to: {out_file}")
    appended_tags.close()

    # write index
    pysam.index(output_file)
    print(f"Index written for {out_file}.bai")


tags_dict = fetch_modified_bases(methyl_bam)

write_linked_tags(aln_bam, tags_dict, output_file)
