#!/usr/bin/env python
# coding=utf-8
import sys
from copy import copy

from setuptools import setup, find_packages

setup(
    name = "dimspy",
    version="0.1",
    author="Keiron O'Shea",
    author_email="keo7@aber.ac.uk",
    url="https://www.github.com/KeironO/DIMSpy",
    download_url="https://www.github.com/KeironO/DIMSpy/zipball/master",
    description="",
    long_description="",
    packages=["dimspy"],
    include_package_data=True,
    package_data={
        '': ['*.txt', '*.rst', '*.md'],
    },
    exclude_package_data={'': ['README.txt']},

    install_requires=[
        'numpy>=1.11.3',
        'scipy>=0.18.1',
        'pymzml>=0.7.7',
        'pandas>=0.18.1',
        "scikit-learn>=0.18.1"
    ],

    keywords='bioinformatics metabolomics research analysis science',

)