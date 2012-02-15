cimport cython

import numpy as np
cimport numpy as np

cimport mpi4py.MPI as MPI
from mpi4py.mpi_c cimport *

interp_f={}
current_cpl = "" 

cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)

cdef extern from "fileobject.h":
    ctypedef class __builtin__.file [object PyFileObject]:
        pass

