**CWIPI** (Coupling With Interpolation Parallel Interface) is a C++/C parallel coupling library under LGPL.

## Documentation  ##
 
User documentation is deployed on the Gitlab pages server: https://numerics.gitlab-pages.onera.net/coupling/cwipi/cwipi-1.0.0/index.html

## Build and install ##

### Basic Installation

Follow these steps to build CWIPI from the sources:

1. `git clone git@gitlab.onera.net:numerics/coupling/cwipi.git` (for Onera users)
2. `cd cwipi`
3. `mkdir build`
4. `cd build`
5. `cmake ..`
6. `make`
7. `make install`
8. `./cwp_run` (if you want to run the test cases)

### CMake general options

cmake -D<option1_name>=<option1_value> ... -D<option_name>=<option_value>

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

Enable client-server mode :
    CWP_ENABLE_CLIENT_SERVER=<ON | OFF> (default : OFF)

Enable documentation mode :
     CWP_ENABLE_DOCUMENTATION=<ON | OFF> (default : OFF)
     Once build, the documentation can be found in build/doc/sphinx/html and launch `index.html` file

### Compiler choice

CC=<C compiler> CXX=<CXX compiler> FC=<Fortran compiler> cmake ...

or use the following CMake options
    CMAKE_C_COMPILER=<C compiler>
    CMAKE_CXX_COMPILER=<CXX compiler>
    CMAKE_Fortran_COMPILER=<Fortran compiler>

### CMake MPI options

    MPI_C_COMPILER=<C mpi wrapper>
    MPI_CXX_COMPILER=<CXX mpi wrapper>
    MPI_Fortran_COMPILER=<Fortran mpi wrapper>

If a simple autodetection fails, you can use these options to find MPI :
    MPI_<lang>_LIBRARIES
    MPI_<lang>_INCLUDE_PATH

Refere to FindMPI in the CMake documentation for more informations.

## Issues ##

Issues can be reported directly on [the Issues section](https://gitlab.onera.net/numerics/coupling/cwipi/-/issues).

## License ##

`CWIPI` is available under the LGPLv3 license (https://www.gnu.org/licenses/lgpl-3.0.fr.html).
