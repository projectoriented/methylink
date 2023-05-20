# methylink
[![GitHub actions status](https://github.com/projectoriented/methylink/workflows/Tests/badge.svg?branch=main)](https://github.com/projectoriented/methylink/actions?query=branch%3Amain+workflow%3ATests)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/projectoriented/methylink/blob/main/LICENSE)

A command line tool to link methylated sites between two BAM/SAM files of the same origin.

## Installation
In a **clean** environment
```shell
pip install methylink
```

## Usage:
```shell 
Usage: methylink [OPTIONS] COMMAND [ARGS]...

  A command line tool to link methylated sites between two BAM files of the
  same origin.

Options:
  -h, --help                      Show this message and exit.
  --version                       Show the version and exit.
  --log_level [DEBUG|INFO|WARNING|ERROR|CRITICAL]
                                  Set the level of log output.  [default:
                                  INFO]
  --threads INTEGER               Number of threads to use.  [default: 1]
  --tmp PATH                      Temp directory to use.
  --methyl_bams TEXT              Unmapped bam files with methylation tags.
                                  [required]
  --aln TEXT                      Aligned bam to map the meth tags to.
                                  [required]
  --sample TEXT                   Sample name.  [required]
  --output TEXT                   Output file.  [required]

```

## Quick start
```shell
methylink --threads 2 --aln tests/data/CHM1_aln_test-subsampled.bam --sample CHM1 --methyl_bams $(ls tests/data/CHM1_methylated_test-*.bam | tr '\n' ',') --output CHM1-linked.bam
```

## Development
I'm happy with any contributions to make this code better :muscle:. You should be able to go forth with the following:
```shell
git clone [this repo]
cd methylink
python -m venv vmeth
source vmeth/bin/activate

pip install --editable .
pytest
```