"""Module for base CLI"""
import click

import tempfile

from methylink.append_mod_tags import AppendModTags, ScatterGather

import multiprocessing as mp

@click.group("base")
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
    "--output",
    required=True,
    help="Output file.",
)
def base(threads, methyl_bams, aln_bam, output, tmp=None):
    
    prefix = tempfile.mkdtemp(suffix="_methylink", dir=tmp)
    
    methylink_obj = AppendModTags(threads=threads, methyl_collection=methyl_bams, aln_bam=aln_bam, output=output, prefix=prefix)
    
    # Create database
    methylink_obj.create_database()

    # Populate the database with methylation tags
    methylink_obj.collect_tags(methyl_collection=methyl_bams)

    # Make the chunks
    scattergather_obj = ScatterGather(aln_bam=methylink_obj.aln_bam, prefix=methylink_obj.prefix, output=methylink_obj.output)
    scattergather_obj.make_subset_bams()

    with mp.Pool(threads) as p:
        p.starmap(
            methylink_obj.run_pool,
            zip(
                scattergather_obj.chunked_bams_names, 
                scattergather_obj.link_bam_output_names
                ),
        )
        p.close()
        p.join()

    scattergather_obj.combine_the_chunked()
    
    # CLEANING UP!
    scattergather_obj.clean_up_tempdir()