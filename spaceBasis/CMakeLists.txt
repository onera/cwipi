#==============================================================================
#
#                Top level CMakeLists.txt file for BASESNODALES
#
#==============================================================================

# Require CMake 3.0
cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

enable_testing()

#------------------------------------------------------------------------------
# Internal CMake module (shared with paradigm)
#------------------------------------------------------------------------------

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
project(BASESNODALES VERSION 2017.1.1 LANGUAGES C Fortran)

message("Source du projet:      ${BASESNODALES_SOURCE_DIR}")
message("Compilation du projet: ${BASESNODALES_BINARY_DIR}")

#------------------------------------------------------------------------------
# Definition de la version BASESNODALES
#------------------------------------------------------------------------------

set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

set(BASESNODALES_VERSION "${BASESNODALES_VERSION_MAJOR}.${BASESNODALES_VERSION_MINOR}.${BASESNODALES_VERSION_MICRO}")

#------------------------------------------------------------------------------
# Localisation des fonctions cmake du projet
#------------------------------------------------------------------------------

set(BASESNODALES_CMAKE_DIR "${BASESNODALES_SOURCE_DIR}/cmake" CACHE INTERNAL "")
set(CMAKE_MODULE_PATH "${BASESNODALES_CMAKE_DIR}/modules")

#------------------------------------------------------------------------------
# Activation des langages Fortran C
#------------------------------------------------------------------------------

enable_language(Fortran)
enable_language(C)


# Shared libraries
# ----------------

if(PDM_IN_CWP)

 option(SNB_ENABLE_SHARED "Build Shared Libraries" ${CWP_ENABLE_SHARED})

else()

 option(SNB_ENABLE_SHARED "Build Shared Libraries" ON)

endif()

# Static libraries
# ----------------

if(PDM_IN_CWP)

 option(SNB_ENABLE_STATIC "Build Static Libraries" ${CWP_ENABLE_STATIC})

else()

 option(SNB_ENABLE_STATIC "Build Static Libraries" OFF)

endif()

#------------------------------------------------------------------------------
# Construction par defaut : mode release
#------------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

#------------------------------------------------------------------------------
# Flags de compilation Fortran par defaut suivant le compilateur
#------------------------------------------------------------------------------


include(${CMAKE_MODULE_PATH}/default_flags.cmake)

add_definitions (-DF77API)

if (CMAKE_C_COMPILER_ID STREQUAL "Intel")
   set (BLA_VENDOR "Intel10_64lp")
endif()

find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

#------------------------------------------------------------------------------
# Ajout des sources de tous les elements du projet
# En fonction des top_targets a installer
#------------------------------------------------------------------------------

add_subdirectory(src)
add_subdirectory(tests)


#------------------------------------------------------------------------------
# Print summary of found and not found optional packages
#------------------------------------------------------------------------------

include(FeatureSummary)
FEATURE_SUMMARY(WHAT ALL)

