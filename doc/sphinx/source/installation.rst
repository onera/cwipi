.. _installation:

Installation
############

Basic Installation
==================

.. code-block:: sh

  mkdir build
  cd build
  cmake ..
  make
  make install

If the installation fails, use the following CMake options.



CMake general options
=====================

.. code-block:: sh

    cmake -D<option1_name>=<option1_value> ... -D<option_name>=<option_value>

Prefix : ``CMAKE_INSTALL_PREFIX=<prefix>``

Enable Fortran interface : ``CWP_ENABLE_Fortran=<ON | OFF>`` (default : ``OFF``)

Enable Python interface : ``CWP_ENABLE_PYTHON_BINDINGS=<ON | OFF>`` (default : ``OFF``)

If a simple autodetection fails, you can use these options to find Python :

.. code-block:: sh

    PYTHON_LIBRARY=<path>
    PYTHON_INCLUDE_DIR=<path>

Refer to FindPython in the CMake documentation for more information.

Enable shared libraries : ``CWP_ENABLE_SHARED=<ON | OFF>`` (default : ``ON``)

Enable static libraries : ``CWP_ENABLE_STATIC=<ON | OFF>`` (default : ``ON``)

Enable the use of the library BLAS : ``CWP_ENABLE_BLAS=<ON | OFF>`` (default : ``ON``)

If a simple autodetection fails, you can use these options to find BLAS :

.. code-block:: sh

    BLAS_DIR=<path>     # Where to find the base directory of BLAS
    BLAS_INCDIR=<path>  # Where to find the header files
    BLAS_LIBDIR=<path>  # Where to find the library files

To force the use of a list of libraries : ``DBLAS_LIBRARIES="<lib_1> ... <lib_n>"``

Compiler choice
===============

.. code-block:: sh

    CC=<C compiler> CXX=<CXX compiler> FC=<Fortran compiler> cmake ...

or use the following CMake options

.. code-block:: sh

    CMAKE_C_COMPILER=<C compiler>
    CMAKE_CXX_COMPILER=<CXX compiler>
    CMAKE_Fortran_COMPILER=<Fortran compiler>


CMake MPI options
=================

.. code-block:: sh

    MPI_C_COMPILER=<C MPI wrapper>
    MPI_CXX_COMPILER=<CXX MPI wrapper>
    MPI_Fortran_COMPILER=<Fortran MPI wrapper>

If a simple autodetection fails, you can use these options to find MPI :

.. code-block:: sh

    MPI_<lang>_LIBRARIES
    MPI_<lang>_INCLUDE_PATH

Refer to FindMPI in the CMake documentation for more information.
