#!/usr/bin/env python3

import pysam

from utils import ROOT_DIR

aln_file = ROOT_DIR + "/data/CHM1_aln_test-subsampled.bam"
m5_file = ROOT_DIR + "/data/CHM1_m5_test-subsampled.bam"

aln_bam = pysam.AlignmentFile(aln_file, "rb", check_sq=False)
m5_bam = pysam.AlignmentFile(m5_file, "rb", check_sq=False)


def fetch_modified_bases(modified_obj):
    # this function assumes a bam with just modified bases in optional tags
    tags_dict = {}
    for read in modified_obj.fetch(until_eof=True):
        if read.has_tag("Mm"):
            tags = read.get_tags()
            qname = read.query_name
            tags_dict[qname] = tags
    m5_bam.close()
    return tags_dict


tags_dict = fetch_modified_bases(m5_bam)
out_file = ROOT_DIR + "/data/linked.bam"


def write_linked_tags(bam, dict_tags, out_file):
    # bam: equivalent aligned bam
    # dict_tags: {query_name: Mm tags and possibly Ml}
    appended_tags = pysam.AlignmentFile(out_file, "wb", template=bam)
    counter=0
    for read in bam.fetch(until_eof=True):
        if read.query_name in dict_tags.keys() and read.is_mapped:
            read.set_tags(read.get_tags() + dict_tags[read.query_name])
            counter += 1
        appended_tags.write(read)
    print(counter)
    print(f"File written to: {out_file}")
    appended_tags.close()


write_linked_tags(aln_bam, tags_dict, out_file)
