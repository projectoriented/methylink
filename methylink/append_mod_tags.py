"""
Link an aligned bam to its unmapped counterpart(s) with methylation tags.
Author: Mei Wu, github.com/projectoriented
"""

import os
import time
import pickle
import sqlite3
import logging
import shutil

import pysam

from methylink.database import crud

LOG = logging.getLogger(__name__)


class AppendModTags:
    def __init__(self, prefix, threads, methyl_collection, aln_path, output) -> None:
        self.prefix = prefix
        self.threads = threads
        self.methyl_collection = methyl_collection
        self.aln_obj = pysam.AlignmentFile(aln_path, check_sq=False)
        self.output = output
        self.db_name = os.path.join(self.prefix, crud.DB_NAME)

    def create_database(self):
        crud.create_database(db_name=self.db_name)
        return self.db

    @property
    def db(self):
        db_obj = crud.DataBaseAdapter(db_name=self.db_name)
        return db_obj

    def fetch_modified_bases(self, modified_obj) -> None:
        """
        Fetch base modification tags Mm & Ml
        """
        LOG.info(f"Opening {modified_obj.filename.decode()} to fetch tags")

        for read in modified_obj.fetch(until_eof=True):
            if read.has_tag("Mm") or read.has_tag("MM"):
                tags = read.get_tags()
                qname = read.query_name

                # serialize the tags list
                serialized_list = pickle.dumps(tags)

                crud.insert_one(
                    qname=str(qname), tag=sqlite3.Binary(serialized_list), db=self.db
                )

        modified_obj.close()
        LOG.info(f"Base modification tags fetched for {modified_obj.filename.decode()}")

    def collect_tags(self, methyl_collection: list) -> None:
        """
        Collect optional tags from ONT bam with methyl calls
        """
        for bam in methyl_collection:
            methyl_bam = pysam.AlignmentFile(bam, "rb", check_sq=False)
            self.fetch_modified_bases(modified_obj=methyl_bam)

    def write_linked_tags(self, aln_obj, out_file) -> None:
        """
        Write out merged bam with Mm tags and possibly Ml, and its index.
        """
        appended_tags = pysam.AlignmentFile(out_file, "wb", template=aln_obj)
        for read in aln_obj.fetch(until_eof=True):
            result = crud.select_one(qname=str(read.qname), db=self.db)
            if result:
                deserialized_tag = pickle.loads(result[0])
                read.set_tags(read.get_tags() + deserialized_tag)
            else:
                pass
            appended_tags.write(read)

        LOG.info(f"File written to: {out_file}")
        appended_tags.close()

        # write index
        pysam.index(out_file)
        LOG.info(f"Index written for {out_file}.bai")

    @staticmethod
    def clean_aln(aln_path: str):
        suffix = ".bai"
        files = [aln_path, aln_path + suffix]
        for f in files:
            try:
                os.remove(f)
            except FileNotFoundError:
                LOG.warning(f"{f} not found.")

    def run_pool(self, chunked_aln_fp: str, outfile: str) -> None:
        chunked_aln_obj = pysam.AlignmentFile(chunked_aln_fp, "rb")

        self.write_linked_tags(aln_obj=chunked_aln_obj, out_file=outfile)

        # wait for the file to become available
        while not os.path.exists(outfile):
            time.sleep(1)

        # file is available, remove it
        self.clean_aln(aln_path=chunked_aln_fp)


class ScatterGather:
    def __init__(self, aln_path, prefix, output) -> None:
        self.aln_obj = pysam.AlignmentFile(aln_path, check_sq=False)
        self.prefix = prefix
        self.output = output

    def make_subset_bams(self) -> list[str]:
        subset_size = 100 * 1024 * 1024  # 100MB in bytes

        if os.path.getsize(self.aln_obj.filename.decode()) < subset_size:
            subset_size = int(subset_size / 10)

        subset_idx = 0
        subset_size_bytes = 0
        current_subset = None

        bam_file_list = []

        for read in self.aln_obj:
            # If the current subset is None or its size has exceeded the subset size, create a new subset
            if current_subset is None or subset_size_bytes >= subset_size:
                # If this is not the first subset, close the previous subset file
                if current_subset is not None:
                    current_subset.close()
                    pysam.index(current_subset.filename.decode())

                # Create a new subset file with a name based on the subset index
                subset_idx += 1
                current_subset = pysam.AlignmentFile(
                    f"{self.prefix}_tmp.{subset_idx}.bam", "wb", template=self.aln_obj
                )
                bam_file_list.append(f"{self.prefix}_tmp.{subset_idx}.bam")

            # Write the current read to the current subset file
            current_subset.write(read)
            subset_size_bytes = os.path.getsize(current_subset.filename.decode())

        # Close the last subset file
        if current_subset is not None:
            current_subset.close()
            pysam.index(current_subset.filename.decode())

        self.aln_obj.close()

        return bam_file_list

    def combine_the_chunked(self, linked_bam_output_fp):
        # Read in the chunked bams
        aln_bams = [
            pysam.AlignmentFile(x, check_sq=False) for x in linked_bam_output_fp
        ]

        out_bam = pysam.AlignmentFile(self.output, "wb", template=aln_bams[0])
        for bam in aln_bams:
            for records in bam:
                out_bam.write(records)
            bam.close()

        out_bam.close()
        pysam.index(self.output)

    def clean_up_tempdir(self):
        parent_dir = os.path.dirname(self.prefix)
        try:
            shutil.rmtree(parent_dir)
            LOG.debug(f"Cleaned up {parent_dir}.")
        except OSError:
            LOG.warning(f"Cannot find {parent_dir}")
