#!/usr/bin/env python3

import os
import sys
import time
import math
import itertools
import gzip
import pickle
from multiprocessing import Pool

import pysam


# LOGGING
sys.stdout = open(snakemake.log[0], "w")

def get_time():
    time_rn = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    return time_rn

def fetch_modified_bases(modified_obj) -> dict:
    """
    Fetch base modification tags Mm & Ml
    :param modified_obj: An unsorted bam pysam object with just methylation calls
    :return: A dictionary of tags where keys = query name and value = list of optional tags
    """
    tags_dict = {}
    print(f"Opening {modified_obj.filename.decode()} to fetch tags. {get_time()}")

    for read in modified_obj.fetch(until_eof=True):
        if read.has_tag("Mm") or read.has_tag("MM"):
            tags = read.get_tags()
            qname = read.query_name
            tags_dict[qname] = tags
    modified_obj.close()

    print(f"Base modification tags fetched for {modified_obj.filename.decode()}. {get_time()}")
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
    for read in bam.fetch(until_eof=True):
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


def make_subset_bams(bam, n_splits, prefix):

    # Create output BAM files equivalent to n_splits
    outbams = [pysam.AlignmentFile(f"{prefix}_tmp.{i}.bam", "wb", template=bam) for i in range(n_splits)]

    print(f'Chunking up the bam into {n_splits} parts')
    # Iterate over the reads in the input BAM file
    for read in bam:
        # Determine which output BAM file to write the read to
        output_idx = hash(read.qname) % n_splits
        outbam = outbams[output_idx]

        # Write the read to the output BAM file
        outbam.write(read)

    # Close all the BAM files and write index for each
    bam.close()
    for outbam in outbams:
        outbam.close()
        pysam.index(outbam.filename.decode())


def combine_the_chunked(bams: list[str], merge_output: str):
    aln_bams = [pysam.AlignmentFile(x, check_sq=False) for x in bams]

    out_bam = pysam.AlignmentFile(f'{merge_output}.unsorted', "wb", template=aln_bams[0])
    for bam in aln_bams:
        for records in bam:
            out_bam.write(records)
        bam.close()

    out_bam.close()
    pysam.sort('-o', merge_output, f'{merge_output}.unsorted')
    pysam.index(merge_output)
    os.remove(f'{merge_output}.unsorted')


def run_pool(bam_file: str, tags_dict_pkl: str, output_file) -> None:
    aln_bam = pysam.AlignmentFile(bam_file, "rb")

    with gzip.open(tags_dict_pkl, 'rb') as f:
        tags_dict = pickle.load(f)

    write_linked_tags(aln_bam, tags_dict, output_file)

    wait_time = 1  # wait 1 second between each check

    # fait for the file to become available
    while not os.path.exists(bam_file):
        time.sleep(wait_time)

    # file is available, remove it
    os.remove(bam_file)


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
    aln_bam = snakemake.input.aln_bam
    bam = pysam.AlignmentFile(aln_bam, check_sq=False)
    prefix = os.path.join(snakemake.resources.tmpdir, snakemake.wildcards.sample)
    final_output = snakemake.output.linked_bam

    # Make sure that each chunk is roughly 100 MB
    file_size = os.path.getsize(aln_bam)
    chunk_size = 100 * 1024 * 1024  # 100 MB in bytes
    def get_n_splits():
        num_chunks = math.ceil(file_size / chunk_size)
        return num_chunks
    n = 10 if file_size < chunk_size else get_n_splits()

    # Get the meth dictionary
    tags_dict = collect_tags(methyl_collection)
    pkl_out = f'{prefix}_tags-dict.pickle.gz'
    with gzip.open(pkl_out, 'wb') as f:
        pickle.dump(tags_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    del tags_dict

    # Make the chunks
    make_subset_bams(bam=bam, n_splits=n, prefix=prefix)

    # Gather the arguments for run_pool
    chunked_bams = [f"{prefix}_tmp.{idx}.bam" for idx in range(0, n)]
    link_bam_output_names = [f"{prefix}_tmp.{idx}-linked.bam" for idx in range(0, n)]

    with Pool(threads) as p:
        p.starmap(run_pool, zip(chunked_bams, itertools.repeat(pkl_out), link_bam_output_names))
        p.close()
        p.join()

    combine_the_chunked(link_bam_output_names, final_output)

    # CLEANING UP!
    clean_up_temps(chunked_bams)
    clean_up_temps(link_bam_output_names)
    clean_up_temps([pkl_out])


if __name__ == '__main__':
    sys.exit(main())
