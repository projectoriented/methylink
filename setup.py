from setuptools import setup, find_packages


REQUIREMENTS = [
    "click >= 8.1.3",
    "pysam >= 0.21.0",
    "pytest"
]

setup(
    name='methylink',
    version="0.1.0",
    author="Mei Wu",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        REQUIREMENTS
    ],
    python_requires=">=3.9.13",
    entry_points={
        'console_scripts': [
            'methylink = methylink.main:base',
        ],
    },
)