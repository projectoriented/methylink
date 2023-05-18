from setuptools import setup, find_packages

requirements = [
    "click >= 8.1.3",
    "pysam",

]
setup(
    name='methylink',
    version='0.1.0',
    author="Mei Wu",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        requirements
    ],
    entry_points={
        'console_scripts': [
            'methylink = methylink.main:base',
        ],
    },
)