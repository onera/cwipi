# Top level CMakeLists.txt file for PTK

# Require CMake 3.0
cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
project(BASESNODALES VERSION 2015.0.0 LANGUAGES C Fortran)

message("Source du projet:      ${BASESNODALES_SOURCE_DIR}")
message("Compilation du projet: ${BASESNODALES_BINARY_DIR}")

#------------------------------------------------------------------------------
# Definition de la version BASESNODALES
#------------------------------------------------------------------------------
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

set(BASESNODALES_VERSION "${BASESNODALES_VERSION_MAJOR}.${BASESNODALES_VERSION_MINOR}.${BASESNODALES_VERSION_MICRO}")

#------------------------------------------------------------------------------
# Utilisation des flags de compilation par dÈfaut du projet
# (Flags par dÈfaut sont dÈfinis plus bas dans le fichier) 
#------------------------------------------------------------------------------

option(BASESNODALES_ENABLE_DEFAULT_FLAGS "Enable default flags" ON)


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


#------------------------------------------------------------------------------
# Flags de compilation C par defaut suivant le compilateur
#------------------------------------------------------------------------------

add_definitions (-DF77API)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

#------------------------------------------------------------------------------
# Ajout des sources de tous les elements du projet
# En fonction des top_targets a installer
#------------------------------------------------------------------------------

add_subdirectory(src)

