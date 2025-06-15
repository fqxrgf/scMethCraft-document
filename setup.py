#!/usr/bin/env python
#-*- coding:utf-8 -*-

from setuptools import setup, find_packages

setup(name='scMethCraft',
    version='1.0.1',
    keywords=("pip", "scMethCraft", "single-cell"),
    url="https://hub.njuu.cf/BioX-NKU/",
    author="BioX-NKU",
    packages=find_packages(),

    python_requires='>=3.8.13',


    license='MIT Licence',
    install_requires=[
        'scanpy',
        'pysam',
        'torch',
        'seaborn',
     ],
    classifiers=['Intended Audience :: Science/Research',
      'License :: OSI Approved :: MIT License',
      'Programming Language :: Python :: 3.8',
      'Operating System :: MacOS :: MacOS X',
      'Operating System :: Microsoft :: Windows',
      'Operating System :: POSIX :: Linux',
      'Topic :: Scientific/Engineering :: Bio-Informatics']
     )




