from setuptools import setup, find_packages

import methylink

requirements = [
    "click",
    "pysam",
    "pytest"
]

setup(
    name='methylink',
    version=methylink.__version__,
    author="Mei Wu",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        requirements
    ],
    python_requires=">=3.8.3",
    entry_points={
        'console_scripts': [
            'methylink = methylink.main:base',
        ],
    },
)