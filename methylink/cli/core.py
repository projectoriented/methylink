"""Module for base CLI"""
import logging
import sys
import os
import tempfile

import click

from methylink.append_mod_tags import AppendModTags, ScatterGather

import concurrent.futures as cf

LOG_LEVELS = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]


@click.command("base")
@click.option(
    "--threads",
    show_default=True,
    default=1,
    help="Number of threads to use.",
)
@click.option(
    "--tmp",
    required=False,
    type=click.Path(),
    default=None,
    help="Temp directory to use.",
)
@click.option(
    "--methyl_bams",
    required=True,
    multiple=True,
    help="Unmapped bam files with methylation tags.",
)
@click.option(
    "--aln_bam",
    required=True,
    help="Aligned bam to map the meth tags to."
)
@click.option(
    "--sample",
    required=True,
    help="Sample name.",
)
@click.option(
    "--output",
    required=True,
    help="Output file.",
)
@click.option(
    "--log_level",
    required=False,
    default="INFO",
    type=click.Choice(LOG_LEVELS),
    help="Set the level of log output.",
    show_default=True,
)
def base(sample, threads, methyl_bams, aln_bam, output, log_level, tmp=None):
    # Logging
    logging.basicConfig(stream=sys.stdout, level=log_level)

    prefix = tempfile.mkdtemp(suffix="_methylink", dir=tmp)

    methylink_obj = AppendModTags(threads=threads, methyl_collection=methyl_bams, aln_bam_path=aln_bam, output=output,
                                  prefix=prefix)

    # Create database
    methylink_obj.create_database()

    # Populate the database with methylation tags
    methylink_obj.collect_tags(methyl_collection=methyl_bams)

    # Make the chunks
    scattergather_obj = ScatterGather(
        aln_bam_path=aln_bam,
        prefix=os.path.join(prefix, sample),
        output=output
    )
    chunked_bam_names = scattergather_obj.make_subset_bams()
    link_bam_output_names = [x.replace("_tmp.", "_tmp-linked.") for x in chunked_bam_names]

    with cf.ThreadPoolExecutor(threads) as executor:
        executor.map(methylink_obj.run_pool, chunked_bam_names, link_bam_output_names)

    scattergather_obj.combine_the_chunked(linked_bam_output_fp=link_bam_output_names)

    # # CLEANING UP!
    scattergather_obj.clean_up_tempdir()
