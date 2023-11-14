import io
import os
import re

from setuptools import setup, find_packages


AUTHOR = "Mei Wu"
DESCRIPTION = "methylink is a tool to link methylation tags between SAM/BAM files."
HERE = os.path.abspath(os.path.dirname(__file__))
NAME = "methylink"
REQUIREMENTS = ["click >= 8.1.3", "pysam >= 0.21.0", "pytest"]

try:
    with io.open(os.path.join(HERE, "README.md"), encoding="utf-8") as f:
        LONG_DESCRIPTION = "\n" + f.read()
except FileNotFoundError:
    LONG_DESCRIPTION = DESCRIPTION

with io.open(os.path.join(HERE, "methylink", "__init__.py"), encoding="utf-8") as f:
    pattern = re.compile(r".*(?P<version>[0-9]\.[0-9]\.[0-9])")
    version_file = f.read()
    METHYLINK_VERSION = pattern.match(version_file).groupdict()["version"]

setup(
    name=NAME,
    version=METHYLINK_VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    keywords="methylation nanopore bioinformatics",
    license="MIT",
    author=AUTHOR,
    packages=find_packages(),
    include_package_data=True,
    url="https://github.com/projectoriented/methylink",
    install_requires=[REQUIREMENTS],
    python_requires=">=3.9.13",
    entry_points={
        "console_scripts": [
            "methylink = methylink.main:base",
        ],
    },
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.9",
        "Environment :: Console",
    ],
)
