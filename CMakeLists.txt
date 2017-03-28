
#==============================================================================
#
#                Top level CMakeLists.txt file for CWIPI
#
#==============================================================================

cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

enable_testing()

#------------------------------------------------------------------------------
# Internal CMake module (shared with paradigm)
#------------------------------------------------------------------------------

set(CWP_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/paradigm/cmake")
set(CMAKE_MODULE_PATH "${CWP_CMAKE_DIR}/modules")

include(${CMAKE_MODULE_PATH}/version.cmake)
project(CWP VERSION 1.0.0-rc1 LANGUAGES C CXX Fortran)

set(PDM_IN_CWP 1)
set(PDM_VERSION "${PDM_VERSION_MAJOR}.${PDM_VERSION_MINOR}.${PDM_VERSION_PATCH}")

#------------------------------------------------------------------------------
# Enable languages
#------------------------------------------------------------------------------

enable_language(Fortran)
enable_language(C)
enable_language(CXX)

#------------------------------------------------------------------------------
# User option Definition
#------------------------------------------------------------------------------

# OpenMP
# ------

option(CWP_ENABLE_OPENMP "Enable OpenMP" OFF)

# Enable long for absolute number
# -------------------------------

#option(CWP_ENABLE_LONG_G_NUM "Enable long for global number" OFF)

# Pt-scotch
# ---------

option(CWP_ENABLE_PTSCOTCH "Compile with support for Pt-SCOTCH." OFF)

# ParMetis
# --------

option(CWP_ENABLE_PARMETIS "Compile with support for ParMetis." OFF)

# Shared libraries
# ----------------

option(CWP_ENABLE_SHARED "Build Shared Libraries" OFF)

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

option(CWP_ENABLE_SHARED_Fortran_INTERFACE "Build Shared Fortran Interface" OFF)

option(CWP_ENABLE_STATIC_Fortran_INTERFACE "Build Static Fortran Interface" ON)

#------------------------------------------------------------------------------
# Default build : Release
#------------------------------------------------------------------------------

if (CWP_ENABLE_PYTHON_BINDINGS AND (NOT CWP_ENABLE_SHARED))
   set (CWP_ENABLE_SHARED "ON" CACHE STRING "Build Shared Libraries" FORCE )
endif()

if ((NOT CWP_ENABLE_STATIC) AND (NOT CWP_ENABLE_SHARED))
   message (FATAL_ERROR "No enabled library. Please, enable shared or static library")
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING  "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()

#------------------------------------------------------------------------------
# ADD "-DDEBUG_CLASSE" in Debug
#------------------------------------------------------------------------------

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions (-DDEBUG_CLASSE)
endif()

#------------------------------------------------------------------------------
# Check include file
#------------------------------------------------------------------------------

check_include_file_concat("inttypes.h"        HAVE_INTTYPES_H)
check_include_file_concat("memory.h"        HAVE_MEMORY_H)
check_include_file_concat("stdbool.h"        HAVE_STDBOOL_H)
check_include_file_concat("stdint.h"        HAVE_STDINT_H)
check_include_file_concat("stdlib.h"        HAVE_STDLIB_H)
check_include_file_concat("strings.h"        HAVE_STRINGS_H)
check_include_file_concat("string.h"        HAVE_STRING_H)
check_include_file_concat("sys/stat.h"        HAVE_SYS_STAT_H)
check_include_file_concat("sys/types.h"        HAVE_SYS_TYPES_H)
check_include_file_concat("sys/unistd.h"        HAVE_UNISTD_H)

#------------------------------------------------------------------------------
# Check type size
#------------------------------------------------------------------------------

check_type_size("size_t"  SIZEOF_SIZE_T)
check_type_size("ssize_t"  SIZEOF_SSIZE_T)
check_type_size("long long"  SIZEOF_LONG_LONG)
check_type_size("long"  SIZEOF_LONG)
check_type_size("short"  SIZEOF_SHORT)
check_type_size("int"  SIZEOF_INT)
check_type_size("int32_t"  SIZEOF_INT32_T)
check_type_size("int64_t"  SIZEOF_INT64_T)
check_type_size("float"  SIZEOF_FLOAT)
check_type_size("double"  SIZEOF_DOUBLE)
check_type_size("long double"  SIZEOF_LONG_DOUBLE)
check_type_size("void *"  SIZEOF_VOID_P)
check_type_size("off_t"  SIZEOF_OFF_T)

#------------------------------------------------------------------------------
# Check C internal functions
#------------------------------------------------------------------------------

include(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(backtrace CWP_HAVE_BACKTRACE)
CHECK_FUNCTION_EXISTS(fseeko CWP_HAVE_FSEEKO)
CHECK_FUNCTION_EXISTS(getrusage CWP_HAVE_GETRUSAGE)
CHECK_FUNCTION_EXISTS(gettimeofday CWP_HAVE_GETTIMEOFDAY)
CHECK_FUNCTION_EXISTS(vprintf CWP_HAVE_VPRINTF)
CHECK_FUNCTION_EXISTS(vnprintf CWP_HAVE_VNPRINTF)

#------------------------------------------------------------------------------
# Default compiler flags 
#------------------------------------------------------------------------------

include(${CMAKE_MODULE_PATH}/default_flags.cmake)

#------------------------------------------------------------------------------
# Libraries to add to link 
#------------------------------------------------------------------------------

set(COMMON_LIBRARIES)

#------------------------------------------------------------------------------
# Check dependecies                                                            
#------------------------------------------------------------------------------

# Check MPI (TODO : check MPI wrapper function check )
# ----------------------------------------------------

find_package(MPI REQUIRED COMPONENTS MPI_C_FOUND MPI_CXX_FOUND MPI_Fortran_FOUND)

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
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Mpi4Py REQUIRED)  
  include(UseCython)
  if (NOT CMAKE_CWP_INSTALL_PYTHON_DIR)

     set (CMAKE_CWP_INSTALL_PYTHON_DIR  "${CMAKE_INSTALL_PREFIX}")
     set (CMAKE_CWP_INSTALL_PYTHON_DIR
         ${CMAKE_CWP_INSTALL_PYTHON_DIR} CACHE STRING "Install path directory for Python packages")

  endif ()

endif()

#------------------------------------------------------------------------------
# Store variable in configure file
#------------------------------------------------------------------------------

configure_file(cwp_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/cwp_config.h")
configure_file(bftc_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/bftc_config.h")
configure_file(fvmc_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/fvmc_config.h")

if (CWP_ENABLE_STATIC_Fortran_INTERFACE OR CWP_ENABLE_SHARED_Fortran_INTERFACE)

function(fortranize file)
    unset(result)
    unset(basename)
    file(STRINGS "${file}" lines)
    get_filename_component(basename ${file} NAME_WE)
    foreach(i IN LISTS lines)
        string(REGEX REPLACE "^/\\* (#undef .*) \\*/$" "! \\1" i "${i}")
        set(result "${result}${i}\n")
    endforeach()
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${basename}f.h" "${result}")
endfunction()

fortranize(${CMAKE_CURRENT_BINARY_DIR}/cwp_config.h)

endif ()

#------------------------------------------------------------------------------
# Doc
#------------------------------------------------------------------------------

find_package(Doxygen)       

#------------------------------------------------------------------------------
# Print summary of found and not found optional packages
#------------------------------------------------------------------------------

include(FeatureSummary)
FEATURE_SUMMARY(WHAT ALL)

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
# Ajout des sources de tous les elements du projet
# En fonction des top_targets a installer
#------------------------------------------------------------------------------

add_subdirectory(paradigm)

add_subdirectory(src)

#if (PDM_ENABLE_PYTHON_BINDINGS)
#  add_subdirectory(Cython)
#endif()

#add_subdirectory(tests)

#------------------------------------------------------------------------------
# Affichage du message de fin d'installation
#------------------------------------------------------------------------------

add_subdirectory(${CWP_CMAKE_DIR}/post-install)





