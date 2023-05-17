from setuptools import setup, find_packages

import os

os.environ["HTSLIB_CONFIGURE_OPTIONS"] = str("--enable-plugins")
print(os.environ['HTSLIB_CONFIGURE_OPTIONS'])

setup(
    name='methylink',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click >= 8.1.3', 'pysam'
    ],
    entry_points={
        'console_scripts': [
            'methylink = methylink.main:base',
        ],
    },
)