#!/usr/bin/env python3

import pysam
import numpy as np
import os
import sys
from multiprocessing import Pool


# LOGGING
sys.stdout = open(snakemake.log[0], "w")


def fetch_modified_bases(modified_obj) -> dict:
    """
    Fetch base modification tags Mm & Ml
    :param modified_obj: An unsorted bam pysam object with just methylation calls
    :return: A dictionary of tags where keys = query name and value = list of optional tags
    """
    tags_dict = {}
    for read in modified_obj.fetch(until_eof=True):
        if read.has_tag("Mm") or read.has_tag("MM"):
            tags = read.get_tags()
            qname = read.query_name
            tags_dict[qname] = tags
    modified_obj.close()
    print(f"Base modification tags fetched for {modified_obj.filename.decode()}")
    return tags_dict


def write_linked_tags(bam, tags_dict, out_file) -> None:
    """
    Write out merged bam with Mm tags and possibly Ml, and its index.
    :param bam: equivalent aligned bam
    :param tags_dict: a dict of {query_name: [Mm tags and possibly Ml]}
    :param out_file: merged bam file path
    :return: None
    """
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
    """
    Collect optional tags from ONT bam with methyl calls
    :param methyl_sn_input: a list of file paths pointing to methyl bam
    :return: a dict of {query_name: [Mm tags and possibly Ml]}
    """
    tags = {}
    if not len(methyl_sn_input) == 1:
        for bam in methyl_sn_input:
            methyl_bam = pysam.AlignmentFile(bam, "rb", check_sq=False)
            dict_of_tags_per_bam = fetch_modified_bases(methyl_bam)
            tags.update(dict_of_tags_per_bam)
    else:
        methyl_bam = pysam.AlignmentFile(methyl_sn_input[0], "rb", check_sq=False)
        dict_of_tags_per_bam = fetch_modified_bases(methyl_bam)
        tags.update(dict_of_tags_per_bam)
    return tags


def get_n_records(bam: pysam.AlignmentFile):
    records = list(map(lambda x: x, bam))
    return records


def make_subset_bams(bam: pysam.AlignmentFile, n_splits):
    partition = np.array_split(get_n_records(bam=bam), n_splits)

    for idx, collection_of_records in enumerate(partition):
        subset_bam = pysam.AlignmentFile(f"tmp.{idx}.bam", "wb", template=bam)
        for record in collection_of_records:
            subset_bam.write(record)

        subset_bam.close()
        pysam.index(f"tmp.{idx}.bam")

    bam.close()


def combine_the_chunked(bams: list[pysam.AlignmentFile], merge_output: str):

    aln_bams = [pysam.AlignmentFile(x, check_sq=False) for x in bams]

    out_bam = pysam.AlignmentFile(merge_output, "wb", template=aln_bams[0])
    for bam in aln_bams:
        for records in bam:
            out_bam.write(records)
        bam.close()

    out_bam.close()
    pysam.index(merge_output)


def execute_the_commands(bam_file:str, methyl_file: list, output_file: str):
    aln_bam = pysam.AlignmentFile(bam_file, "rb")
    tags_dict = collect_tags(methyl_file)
    write_linked_tags(aln_bam, tags_dict, output_file)


def clean_up_temps(files: list):
    for f in files:
        index = f + '.bai'
        try:
            os.remove(f)
            os.remove(index)
            print("removed: ", f)
            print("removed: ", index)
        except FileNotFoundError:
            pass


def main():
    # Grabbing from snakemake
    threads = snakemake.threads
    methyl_collection = snakemake.input.methyl_bam
    final_output = snakemake.output.linked_bam
    bam = pysam.AlignmentFile(snakemake.input.aln_bam, check_sq=False)

    # Make the chunks
    make_subset_bams(bam=bam, n_splits=10)

    chunked_bams = [f"tmp.{idx}.bam" for idx in range(0, 10)]
    link_bam_output_names = [f"tmp.{idx}-linked.bam" for idx in range(0, 10)]

    with Pool(threads) as p:
        p.starmap(execute_the_commands, zip(chunked_bams, link_bam_output_names))
        p.close()

    combine_the_chunked(link_bam_output_names, final_output)

    # Clean up the chunked bams and their index
    clean_up_temps(chunked_bams)

    # Clean up the linked bams and their index
    clean_up_temps(link_bam_output_names)
    

if __name__ == '__main__':
    sys.exit(main())
