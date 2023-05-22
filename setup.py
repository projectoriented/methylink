import io
import os

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

setup(
    name=NAME,
    version="0.1.0",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
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
