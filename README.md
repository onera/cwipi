<p align="center">
    <img src="doc/images/cwipiNew.svg" alt="Logo" width="25%"/>
</p>

**CWIPI** (Coupling With Interpolation Parallel Interface) is a library for coupling parallel scientific codes via MPI communications to perform multi-physics simulations in massively parallel, distributed-memory environments, with interfaces in C, Python and Fortran.

## 📖 Documentation

The user documentation is available either on [the GitHub page](https://onera.github.io/cwipi/index.html) or on [ONERA's internal Gitlab](https://numerics.gitlab-pages.onera.net/coupling/cwipi/index.html).

## 🛠️ Build and install

### Dependencies

General dependencies for building **CWIPI** are:
- a C++ compiler (tested with `gcc 10, 12, 14, 15`, `intel 2022`, `intel-llvm 2022`)
- [CMake](https://cmake.org/) (version 3.16 or higher)
- an MPI distribution (tested with `openmpi` and `intel oneapi`)

To run **CWIPI** tests with `cwp_run`, additional dependencies are:
- `bc`
- `awk`
- [gnu-time](https://formulae.brew.sh/formula/gnu-time) on macOS

To build **CWIPI** Fortran API, additional dependencies are:
- a Fortran 90 compiler (tested with `gfortran` based on `gcc 10, 12, 14, 15`, `intel 2022`, `intel-llvm 2022`)

To build **CWIPI** Python API, additional dependencies are:
- `numpy (1.x)`
- `mpi4py`
- `Cython (0.29.x, 3.x)`

### Basic Installation

Follow these steps to build **CWIPI** from the sources:

First, clone the repository
  - either from [GitHub](https://github.com/onera/cwipi): `git clone git@github.com:onera/cwipi.git`
  - or from [GitLab](https://gitlab.onera.net/numerics/coupling/cwipi) (for ONERA users only): `git clone git@gitlab.onera.net:numerics/coupling/cwipi.git`

Then, use the following commands:
1. `cd cwipi`
2. `git submodule update --init` (needed for dependencies such as [**ParaDiGM**](https://github.com/onera/paradigm))
3. `mkdir build`
4. `cd build`
5. `cmake ..`
6. `make`
7. `make install`
8. `./cwp_run` (if you want to run the test cases, at least 10 CPU cores are required)

<details>

<summary>Configuration with CMake</summary>

### CMake general options
    cmake -D<option1_name>=<option1_value> ... -D<option2_name>=<option2_value>

#### Installation prefix
    CMAKE_INSTALL_PREFIX=<prefix>

#### Enable Fortran interface
    CWP_ENABLE_Fortran=<ON | OFF> (default : OFF)
    CWP_ENABLE_Fortran_MPI_MODULE=<ON | OFF> (default : OFF)

#### Enable Python interface
    CWP_ENABLE_PYTHON_BINDINGS=<ON | OFF> (default : OFF)

If a simple autodetection fails, you can use these options to find Python :

    PYTHON_LIBRARY=<path>
    PYTHON_INCLUDE_DIR=<path>

Refer to [FindPython](https://cmake.org/cmake/help/latest/module/FindPython.html) in the CMake documentation for more information.

#### Build shared library
    CWP_ENABLE_SHARED=<ON | OFF> (default : ON)

#### Build static library
    CWP_ENABLE_STATIC=<ON | OFF> (default : OFF)

#### Hide symbols of internal ParaDiGM library
    CWP_ENABLE_HIDE_PDM_SYMBOLS=<ON | OFF> (default : ON)

If ```CWP_ENABLE_STATIC=ON``` then ```CWP_ENABLE_HIDE_PDM_SYMBOLS=OFF``` is forced

#### Enable the use of external ParaDiGM library
    CWP_ENABLE_EXTERNAL_PDM=<ON | OFF> (default : OFF)

If ```CWP_ENABLE_EXTERNAL_PDM=ON```, you must define this variable to find ParaDiGM:

    PDM_SOURCE_DIR=<path> Where to find the base directory of ParaDiGM

If ```CWP_ENABLE_HIDE_PDM_SYMBOLS=ON``` then ```CWP_ENABLE_EXTERNAL_PDM=OFF``` is forced

</details>

## 🚀 Quick start

You can find a basic example of two-way coupling using **CWIPI** in the [quick start section of the documentation](#documentation).
Additional examples for more advanced scenarios can be found in the ``tests`` directory.

## 🤕 Issues

Issues can be reported directly in the [Issues](https://gitlab.onera.net/numerics/coupling/cwipi/-/issues) section.

## ⚖️ License

**CWIPI** is available under the [LGPL3 license](https://www.gnu.org/licenses/lgpl-3.0.en.html).

## © Copyright

Copyright 2026, ONERA The French Aerospace Lab
