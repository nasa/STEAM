#!/usr/bin/env python
""" A directory named "mpi4py" comes with the mpi4py module and must be
linked into this directory.  This script will provide the path to that
directory so that the Makefile can capture the output of this script and 
link it here.
"""

from mpi4py import __file__ as f
from os.path import dirname

top_dir = dirname( f )
bottom_dir = 'include/mpi4py'
print( '{:s}/{:s}'.format( top_dir, bottom_dir ) )
