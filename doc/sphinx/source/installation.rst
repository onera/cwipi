.. _installation:

Installation
############

Basic Installation
==================

mkdir build
cd build
cmake ..
make
make install

If the installation fails, use the following CMake options.



CMake general options
=====================

cmake -D<option1_name>=<option1_value> ... -D<optionn_name>=<optionn_value>

Prefix :
    CMAKE_INSTALL_PREFIX=<prefix>

Enable fortran interface :
    CWP_ENABLE_Fortran=<ON | OFF> (default : OFF)

Enable python interface :
    CWP_ENABLE_PYTHON_BINDINGS=<ON | OFF> (default : OFF)
      If a simple autodetection fails, you can use these options to find Python :
        PYTHON_LIBRARY=<path>
        PYTHON_INCLUDE_DIR=<path>
      Refere to FindPython in the CMake documentation for more informations.

Enable shared libraries :
    CWP_ENABLE_SHARED=<ON | OFF> (default : ON)

Enable static libraries :
    CWP_ENABLE_STATIC=<ON | OFF> (default : ON)

Enable the use of the library BLAS :
    CWP_ENABLE_BLAS=<ON | OFF> (default : ON)

      If a simple autodetection fails, you can use these options to find BLAS :
        BLAS_DIR=<path>      Where to find the base directory of blas
        BLAS_INCDIR=<path>   Where to find the header files
        BLAS_LIBDIR=<path>   Where to find the library files

      To force the use of a list of libraries, use :
        DBLAS_LIBRARIES="<lib_1> ... <lib_n>"

Compiler choice
===============

CC=<C compiler> CXX=<CXX compiler> FC=<Fortran compiler> cmake ...

or use the following CMake options
    CMAKE_C_COMPILER=<C compiler>
    CMAKE_CXX_COMPILER=<CXX compiler>
    CMAKE_Fortran_COMPILER=<Fortran compiler>


CMake MPI options
=================

    MPI_C_COMPILER=<C mpi wrapper>
    MPI_CXX_COMPILER=<CXX mpi wrapper>
    MPI_Fortran_COMPILER=<Fortran mpi wrapper>

If a simple autodetection fails, you can use these options to find MPI :
    MPI_<lang>_LIBRARIES
    MPI_<lang>_INCLUDE_PATH

Refere to FindMPI in the CMake documentation for more informations.
