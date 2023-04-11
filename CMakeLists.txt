

#==============================================================================
#
#                Top level CMakeLists.txt file for CWIPI
#
#==============================================================================

cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

enable_testing()

# ----------------------------------------------------------------------
# RPATH protection
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

#------------------------------------------------------------------------------
# Internal CMake module (shared with paradigm)
#------------------------------------------------------------------------------

#set(CWP_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/paradigm/cmake")
set(CWP_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${CWP_CMAKE_DIR}/modules")

project(CWP VERSION 0.0.0 LANGUAGES C CXX Fortran)

#------------------------------------------------------------------------------
# Enable languages
#------------------------------------------------------------------------------

set(PDM_IN_CWP 1)
set(PDM_VERSION "${PDM_VERSION_MAJOR}.${PDM_VERSION_MINOR}.${PDM_VERSION_PATCH}")

#------------------------------------------------------------------------------
# User option Definition
#------------------------------------------------------------------------------

# OpenMP
# ------

option(CWP_ENABLE_OPENMP "Enable OpenMP" ON)

# Pt-scotch
# ---------

option(CWP_ENABLE_PTSCOTCH "Compile with support for Pt-SCOTCH." ON)

# ParMetis
# --------

option(CWP_ENABLE_PARMETIS "Compile with support for ParMetis." ON)

# Shared libraries
# ----------------

option(CWP_ENABLE_SHARED "Build Shared Libraries" ON)

# Static libraries
# ----------------

option(CWP_ENABLE_STATIC "Build Static Libraries" ON)

# Check MPI Wrapper
# -----------------

option(CWP_ENABLE_MPI_CHECK "Check MPI Wrapper" ON)

# Python Bindings
# ---------------

option(CWP_ENABLE_PYTHON_BINDINGS "Build Python module" OFF)

# Fortran interface
# -----------------

option(CWP_ENABLE_Fortran  "Build Shared Fortran Interface" OFF)

option(CWP_ENABLE_Fortran_MPI_MODULE "Use Fortran MPI module" OFF)

if (CWP_ENABLE_Fortran_MPI_MODULE)
  set (CWP_HAVE_FORTRAN_MPI_MODULE 1)
endif()

# New to old API
# -----------------

option(CWP_NEW_TO_OLD_API  "Build Fake Old API refering to New one" OFF)

option(CWP_ENABLE_NEW_API  "Build new API (under construction)" ON)

# MKL, BLAS, LAPACK
# -----------------

option(CWP_ENABLE_BLAS  "Enable the use of the library blas" ON)

if (CWP_ENABLE_BLAS)
  if (CMAKE_C_COMPILER_ID STREQUAL "Intel")
    set (BLA_VENDOR "Intel10_64lp_seq")
    find_package(BLAS)
    if (BLAS_FOUND)
      set(CWP_HAVE_MKL 1)
    endif()
    unset (BLA_VENDOR)
  endif()

  if (NOT BLAS_FOUND)
    find_package(BLAS)
  endif()

  find_package(LAPACK)


endif()

if (BLAS_FOUND)
  set(CWP_HAVE_BLAS 1)
endif()


# High order (Nodal bases from SPACE)
# -----------------

#if (BLAS_FOUND AND LAPACK_FOUND)
#  set (ENABLE_SPACE_NODAL_BASIS ON)
#else()
#  set (ENABLE_SPACE_NODAL_BASIS OFF)
#endif()
#option(CWP_ENABLE_SPACE_NODAL_BASIS  "Enable nodal basis from SPACE" ${ENABLE_SPACE_NODAL_BASIS})

#if (CWP_ENABLE_SPACE_NODAL_BASIS)
# set(SNB_IN_CWP 1)
# set(HAVE_SPACE_BASIS 1)
#endif ()

# Allow Fortran call from C
# --------------------------

# Client-Server
# -------------

option(CWP_ENABLE_CLIENT_SERVER "Enable Client-Server" OFF)

# Documentation
# -------------

option(CWP_ENABLE_DOCUMENTATION "Enable sphinx documentation" OFF)

#------------------------------------------------------------------------------
# Default build : Release
#------------------------------------------------------------------------------

if (CWP_ENABLE_PYTHON_BINDINGS AND (NOT CWP_ENABLE_SHARED))
   set (CWP_ENABLE_SHARED "ON" CACHE STRING "Build Shared Libraries" FORCE )
endif()

if (CWP_ENABLE_SHARED AND CWP_ENABLE_Fortran)
   set (CWP_ENABLE_SHARED_Fortran_LIB ON)
endif()

if (CWP_ENABLE_STATIC AND CWP_ENABLE_Fortran)
   set (CWP_ENABLE_STATIC_Fortran_LIB ON)
endif()

if ((NOT CWP_ENABLE_STATIC) AND (NOT CWP_ENABLE_SHARED))
   message (FATAL_ERROR "No enabled library. Please, enable shared or static library")
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING  "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()

#------------------------------------------------------------------------------
# Default compiler flags
#------------------------------------------------------------------------------

include(${CMAKE_MODULE_PATH}/default_flags.cmake)

#------------------------------------------------------------------------------
# ADD "-DDEBUG_CLASSE" in Debug
#------------------------------------------------------------------------------

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions (-DDEBUG_CLASSE)
endif()

#------------------------------------------------------------------------------
# Check include file
#------------------------------------------------------------------------------

include(CheckIncludeFile)
check_include_file("dlfcn.h"        HAVE_DLFCN_H)
check_include_file("fcntl.h"        HAVE_FCNTL_H)
check_include_file("malloc.h"       HAVE_MALLOC_H)
check_include_file("stddef.h"       HAVE_STDDEF_H)

check_include_file("inttypes.h"     HAVE_INTTYPES_H)
check_include_file("memory.h"       HAVE_MEMORY_H)
check_include_file("stdbool.h"      HAVE_STDBOOL_H)
check_include_file("stdint.h"       HAVE_STDINT_H)
check_include_file("stdlib.h"       HAVE_STDLIB_H)
check_include_file("strings.h"      HAVE_STRINGS_H)
check_include_file("string.h"       HAVE_STRING_H)
check_include_file("unistd.h"       HAVE_UNISTD_H)
if (NOT HAVE_UNISTD_H)
  check_include_file("sys/unistd.h"   HAVE_UNISTD_H)
endif()

check_include_file("sys/stat.h"     HAVE_SYS_STAT_H)
check_include_file("sys/types.h"    HAVE_SYS_TYPES_H)
check_include_file("sys/procfs.h"   HAVE_SYS_PROCFS_H)
check_include_file("sys/signal.h"   HAVE_SYS_SIGNAL_H)
check_include_file("sys/sysinfo.h"  HAVE_SYS_SYSINFO_H)
check_include_file("sys/utsname.h"  HAVE_SYS_UTSNAME_H)
check_include_file("sys/resource.h" HAVE_SYS_RESOURCE_H)
check_include_file("sys/time.h"     HAVE_SYS_TIME_H)
check_include_file("time.h"         HAVE_TIME_H)

#------------------------------------------------------------------------------
# Check type size
#------------------------------------------------------------------------------

include(CheckTypeSize)
check_type_size("size_t"       SIZEOF_SIZE_T)
if (HAVE_SIZEOF_SIZE_T)
  set (HAVE_SIZE_T 1)
endif()
check_type_size("_Bool"        SIZEOF__BOOL)
if (HAVE_SIZEOF__BOOL)
  set (HAVE__BOOL 1)
endif()
check_type_size("ssize_t"      SIZEOF_SSIZE_T)
if (HAVE_SIZEOF_SSIZE_T)
  set (HAVE_SSIZE_T 1)
endif()
check_type_size("long long"    SIZEOF_LONG_LONG)
if (HAVE_SIZEOF_LONG_LONG)
  set (HAVE_LONG_LONG 1)
endif()
check_type_size("unsigned long long"  SIZEOF_UNSIGNED_LONG_LONG)
if (HAVE_SIZEOF_UNSIGNED_LONG_LONG)
  set (HAVE_UNSIGNED_LONG_LONG 1)
endif()
check_type_size("long"         SIZEOF_LONG)
if (HAVE_SIZEOF_LONG)
  set (HAVE_LONG 1)
endif()
check_type_size("short"        SIZEOF_SHORT)
if (HAVE_SIZEOF_SHORT)
  set (HAVE_SHORT 1)
endif()
check_type_size("int"          SIZEOF_INT)
if (HAVE_SIZEOF_INT)
  set (HAVE_INT 1)
endif()
check_type_size("int32_t"      SIZEOF_INT32_T)
if (HAVE_SIZEOF_INT32_T)
  set (HAVE_INT32_T 1)
endif()
check_type_size("int64_t"      SIZEOF_INT64_T)
if (HAVE_SIZEOF_INT64_T)
  set (HAVE_INT64_T 1)
endif()
check_type_size("uint32_t"     SIZEOF_UINT32_T)
if (HAVE_SIZEOF_UINT32_T)
  set (HAVE_UINT32_T 1)
endif()
check_type_size("uint64_t"     SIZEOF_UINT64_T)
if (HAVE_SIZEOF_UINT64_T)
  set (HAVE_UINT64_T 1)
endif()
check_type_size("unistd_t"     SIZEOF_UNISTD_T)
if (HAVE_SIZEOF_UNISTD_T)
  set (HAVE_UNISTD_T 1)
endif()
check_type_size("float"        SIZEOF_FLOAT)
if (HAVE_SIZEOF_FLOAT_T)
  set (HAVE_FLOAT_T 1)
endif()
check_type_size("double"       SIZEOF_DOUBLE)
if (HAVE_SIZEOF_DOUBLE_T)
  set (HAVE_DOUBLE_T 1)
endif()
check_type_size("long double"  SIZEOF_LONG_DOUBLE)
if (HAVE_SIZEOF_LONG_DOUBLE_T)
  set (HAVE_LONG_DOUBLE_T 1)
endif()
check_type_size("void *"       SIZEOF_VOID_P)
if (HAVE_SIZEOF_VOID_P)
  set (HAVE_VOID_P 1)
endif()
check_type_size("off_t"        SIZEOF_OFF_T)
if (HAVE_SIZEOF_OFF_T)
  set (HAVE_OFF_T 1)
endif()
check_type_size("ptrdiff_t"    SIZEOF_PTRDIFF_T)
if (HAVE_SIZEOF_PTRDIFF_T)
  set (HAVE_PTRDIFF_T 1)
endif()

#------------------------------------------------------------------------------
# Check C internal functions
#------------------------------------------------------------------------------

include(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(sysinfo HAVE_SYSINFO)
CHECK_FUNCTION_EXISTS(sbrk HAVE_SBRK)
CHECK_FUNCTION_EXISTS(mkdir HAVE_MKDIR)
CHECK_FUNCTION_EXISTS(posix_memalign HAVE_POSIX_MEMALIGN)
CHECK_FUNCTION_EXISTS(backtrace HAVE_BACKTRACE)
CHECK_FUNCTION_EXISTS(fseeko HAVE_FSEEKO)
CHECK_FUNCTION_EXISTS(getrusage HAVE_GETRUSAGE)
CHECK_FUNCTION_EXISTS(gettimeofday HAVE_GETTIMEOFDAY)
CHECK_FUNCTION_EXISTS(vprintf HAVE_VPRINTF)
CHECK_FUNCTION_EXISTS(snprintf HAVE_SNPRINTF)

#------------------------------------------------------------------------------
# Libraries to add to link
#------------------------------------------------------------------------------

set(COMMON_LIBRARIES)

if (CMAKE_BUILD_TYPE STREQUAL "Coverage")
   add_link_options("--coverage")
endif()

#------------------------------------------------------------------------------
# Check dependecies
#------------------------------------------------------------------------------

# Check MPI (TODO : check MPI wrapper function check )
# ----------------------------------------------------

if (CWP_ENABLE_Fortran)
  find_package(MPI REQUIRED COMPONENTS C Fortran)
else()
  find_package(MPI REQUIRED COMPONENTS C)
endif()

mark_as_advanced (MPI_EXTRA_LIBRARY MPI_LIBRARY)

# Check MPI Wrapper
# -----------------

if (CWP_ENABLE_MPI_CHECK)
  include(${CMAKE_MODULE_PATH}/CheckMPI.cmake)
endif()

# Checck SCOTCH an ParMETIS
# -------------------------

if (CWP_ENABLE_PARMETIS)
  find_package(ParMETIS 4.0.3)
  if (PARMETIS_FOUND)
    set(CWP_HAVE_PARMETIS 1)
  endif()
endif()

if (CWP_ENABLE_PTSCOTCH)
  find_package(PTSCOTCH 6.0.0)
  if (PTSCOTCH_FOUND)
    set(CWP_HAVE_PTSCOTCH 1)
  endif()
endif()

# Check OpenMP (TODO : check MPI wrapper function check )
# -------------------------------------------------------

if (CWP_ENABLE_OPENMP)
  find_package(OpenMP)
endif()

# Python bindings
# ---------------

if (CWP_ENABLE_PYTHON_BINDINGS)
  find_package(MPI REQUIRED COMPONENTS CXX)

  find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
  find_package(Mpi4Py REQUIRED)
  include(UseCython)
  if (NOT CMAKE_CWP_INSTALL_PYTHON_DIR)

     set (CMAKE_CWP_INSTALL_PYTHON_DIR  "${CMAKE_INSTALL_PREFIX}")
     set (CMAKE_CWP_INSTALL_PYTHON_DIR
         ${CMAKE_CWP_INSTALL_PYTHON_DIR} CACHE STRING "Install path directory for Python packages")
  endif ()

endif()

#------------------------------------------------------------------------------
# System for bft
#------------------------------------------------------------------------------

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(bftc_OS_Linux 1)
elseif (CMAKE_SYSTEM_NAME STREQUAL "AIX")
  set(bftc_OS_AIX 1)
elseif (CMAKE_SYSTEM_NAME STREQUAL "HP_UX")
  set(bftc_OS_HP_UX 1)
elseif (CMAKE_SYSTEM_NAME STREQUAL "IRIX64")
  set(bftc_OS_IRIX64 1)
elseif (CMAKE_SYSTEM_NAME STREQUAL "OSF1")
  set(bftc_OS_OSF1 1)
elseif (CMAKE_SYSTEM_NAME STREQUAL "SunOS")
  set(bftc_OS_Solaris 1)
elseif (CMAKE_SYSTEM_NAME STREQUAL "System_V")
  set(bftc_OS_System_V 1)
endif()

#------------------------------------------------------------------------------
# Call Fortran in C
#------------------------------------------------------------------------------

if (NOT (CWP_ENABLE_SHARED_Fortran_LIB OR CWP_ENABLE_STATIC_Fortran_LIB))
   set (CWP_HAVE_NOT_FORTRAN_IN_C 1)
endif()

#------------------------------------------------------------------------------
# Store variable in configure file
#------------------------------------------------------------------------------

configure_file(cwp_config_priv.h.in "${CMAKE_CURRENT_BINARY_DIR}/config_priv.h")
configure_file(bftc_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/bftc_config.h")
configure_file(fvmc_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/fvmc_config.h")
configure_file(cwp_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/cwipi_config.h")

if (CWP_ENABLE_Fortran)
  if ("${CMAKE_CURRENT_BINARY_DIR}/cwipi_config.h" IS_NEWER_THAN  "${CMAKE_CURRENT_BINARY_DIR}/cwipi_configf.h")
    include (FortranizeFile)
    fortranize(${CMAKE_CURRENT_BINARY_DIR}/cwipi_config.h)
  endif ()
endif ()

#------------------------------------------------------------------------------
# Installation de PDM
#------------------------------------------------------------------------------

# Set PDM install sub-directories
set(CWP_BIN_DIR "bin")
set(CWP_LIB_DIR "lib")
set(CWP_INCLUDE_DIR "include")
set(CWP_PKGCONFIG_DIR "lib/pkgconfig")
set(CWP_SHARE_DIR "share/cedre")
set(CWP_MAN_DIR "share/man")
set(CWP_DOC_DIR "${CWP_SHARE_DIR}/doc")
set(CWP_ETC_DIR "etc")

#------------------------------------------------------------------------------
# Create test list
#------------------------------------------------------------------------------

set(CWP_LIST_TEST "")

#------------------------------------------------------------------------------
# Ajout des sources de tous les elements du projet
# En fonction des top_targets a installer
#------------------------------------------------------------------------------

if (CWP_ENABLE_NEW_API)
  if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/external/paradigm")
  add_subdirectory(external/paradigm)
#  add_subdirectory(paradigm)
#  add_subdirectory(paradigm/doc)
  endif()
endif()

add_subdirectory(src)

if (CWP_ENABLE_PYTHON_BINDINGS)
  add_subdirectory(Cython)
endif()

if (CWP_ENABLE_SPACE_NODAL_BASIS)
 add_subdirectory(spaceBasis)
endif ()

if (CWP_ENABLE_NEW_API)
  if (CWP_ENABLE_CLIENT_SERVER)
    add_subdirectory(bin)
  endif()
endif()

add_subdirectory(tests)

add_subdirectory(doc)

#------------------------------------------------------------------------------
# Print summary of found and not found optional packages
#------------------------------------------------------------------------------

include(FeatureSummary)
FEATURE_SUMMARY(WHAT ALL)

#------------------------------------------------------------------------------
# Affichage du message de fin d'installation
#------------------------------------------------------------------------------

add_subdirectory(${CWP_CMAKE_DIR}/post-install)

#------------------------------------------------------------------------------
# Print coverage report
#------------------------------------------------------------------------------

if (CMAKE_BUILD_TYPE STREQUAL "Coverage")

  if (NOT CMAKE_C_COMPILER_ID STREQUAL "GNU") # "Clang"
    message(FATAL_ERROR "Code coverage unavailable for compilers other than GCC and Clang")
  endif()

  # remove .gcno
  # TO DO : remove in make clean
  # execute_process(COMMAND find . -type f -name '*.o' -exec rm {} +)
  # execute_process(COMMAND find . -type f -name '*.o.d' -exec rm {} +)
  # execute_process(COMMAND find . -type f -name '*.gcov' -exec rm {} +)
  # execute_process(COMMAND find . -type f -name '*.gcno' -exec rm {} +)
  # execute_process(COMMAND find . -type f -name '*.gcda' -exec rm {} +)

  add_custom_target (distclean @echo cleaning for source distribution)
  add_custom_command(TARGET distclean POST_BUILD
                     COMMAND find ARGS ${CMAKE_BINARY_DIR} -type f -name '*.gcno' -exec rm {} +
                     DEPENDS clean
                     VERBATIM )

  # list of tests
  get_property(TEST_NAMES DIRECTORY ${CMAKE_SOURCE_DIR}/tests PROPERTY TESTS)
  get_property(TUTORIAL_NAMES DIRECTORY ${CMAKE_SOURCE_DIR}/tests/tutorial PROPERTY TESTS)
  list (APPEND TEST_NAMES ${TUTORIAL_NAMES})

  # run test
  add_custom_target(cwp_test ALL
                    COMMAND ctest -R c_new_api_dcube || echo "" # || echo "" to ensure zero error code, non-zero if tests fail
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
                    DEPENDS ${TEST_NAMES}
                    COMMENT "Running tests for coverage")

  # command to generate xml coverage summary
  # TO DO: check if gcov and gcovr installed
  add_custom_target(cwp_coverage ALL
                    COMMAND python3 -m gcovr
                            --cobertura-pretty
                            --root ${CMAKE_SOURCE_DIR}/src/new
                            --object-directory ${CMAKE_BINARY_DIR}
                            # --exclude-directories='(.+/)?paradigm'
                            --output cwp_coverage.xml
                    BYPRODUCTS cwp_coverage.xml
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    DEPENDS cwp_test
                    COMMENT "Running gcovr to produce Cobertura code coverage report")

  # execute command
  if(EXISTS ${CMAKE_BINARY_DIR}/cwp_coverage.xml)
    message(FATAL_ERROR "cwp_coverage.xml already exists, make clean before make again for a new coverage report")
  else()
    add_custom_command(TARGET cwp_coverage POST_BUILD
                       COMMAND ;
                       COMMENT "Cobertura code coverage report saved in cwp_coverage.xml")
  endif()

endif()
