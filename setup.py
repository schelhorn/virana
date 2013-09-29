#!/usr/bin/env python

import sys
import os.path

try:
   from setuptools import setup, Extension
except ImportError:   
   sys.stderr.write( "Could not import 'setuptools', falling back to 'distutils'.\n" )
   from distutils.core import setup, Extension

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
    version='1.0.0',
    author='Sven-Eric Schelhorn',
    author_email='sven@mpi-inf.mpg.de',
    packages=['virana'],
    scripts=['bin/vmap','bin/vref','bin/vhom'],
    license='LICENSE.txt',
    description='Python applicaton for performing next-generation sequence analyses on viral targets.',
    long_description=open('README.txt').read(),
    install_requires=[
        "plumbum",
	"numpy",
	"matplotlib",
        "biopython==1.61",
        "pysam",
        "ftputil>=2.4",
        "HTSeq",
	"numpy"
    ],
)
