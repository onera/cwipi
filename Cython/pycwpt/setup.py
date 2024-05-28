from distutils.core import setup
from distutils.sysconfig import get_python_lib
import glob
import os
import sys

setup(
    name = "pycwpt",
    packages     = ['pycwpt'],
    data_files = [('', ["pycwpt.so"])],
    author = 'E. Quemerais, B. Andrieu, K. Hoogveld, N. Dellinger',
    description = 'Test utilities',
    license = 'LGPL'
    )

