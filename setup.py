#!/usr/bin/env python

import sys

try:
   from setuptools import setup
except ImportError:
   sys.stderr.write( "Could not import 'setuptools', falling back to 'distutils'.\n" )
   from distutils.core import setup

if sys.version_info[0] < 2 or sys.version_info < 5:
   sys.stderr.write( "Error in setup script for Virana:\n" )
   sys.stderr.write( "You need at least version 2.7 of Python to use Virana.\n" )
   sys.exit( 1 )

if sys.version_info[0] >= 3:
   sys.stderr.write( "Error in setup script for Virana:\n" )
   sys.stderr.write( "Sorry, this package does not yet work with Python 3.\n" )
   sys.stderr.write( "Please use Python 2.x, x>=7.\n" )
   sys.exit( 1 )

setup(
    name='Virana',
    version='1.1.0',
    author='Sven-Eric Schelhorn',
    author_email='sven@mpi-inf.mpg.de',
    packages=['virana'],
    scripts=['bin/vmap','bin/vref','bin/vhom', 'bin/vsim'],
    license='LICENSE.txt',
    description='Python application for performing metagenomic next-generation sequencing analyses on microbial targets.',
    long_description=open('README.txt').read(),
    install_requires=[
        "plumbum",
	      "numpy",
	      "matplotlib",
        "biopython==1.61",
        "pysam",
        "ftputil>=2.4",
        "HTSeq"
    ],
)
