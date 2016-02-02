# Top level CMakeLists.txt file for PTK

# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

project(BASESNODALES C Fortran)

message("Source du projet:      ${BASESNODALES_SOURCE_DIR}")
message("Compilation du projet: ${BASESNODALES_BINARY_DIR}")

#------------------------------------------------------------------------------
# Definition de la version BASESNODALES
#------------------------------------------------------------------------------
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

# Defintion de la version pour le futur
# CMake pilotera l'ensemble via la generation de fichier de configuration
# Ces fichiers seront ensuite inclus dans les sources de BASESNODALES
set(BASESNODALES_VERSION_MAJOR "2015")
set(BASESNODALES_VERSION_MINOR "0")
set(BASESNODALES_VERSION_MICRO "0")
set(BASESNODALES_VERSION "${BASESNODALES_VERSION_MAJOR}.${BASESNODALES_VERSION_MINOR}.${BASESNODALES_VERSION_MICRO}")

message("Version du projet: ${BASESNODALES_VERSION}")

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

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
get_filename_component (C_COMPILER_NAME       ${CMAKE_C_COMPILER}       NAME)


if (Fortran_COMPILER_NAME MATCHES "ifort.*")

  if (BASESNODALES_ENABLE_DEFAULT_FLAGS)
    
    # Ajout des flags commun 
    set (CMAKE_Fortran_FLAGS                          "-i4 -r8"           CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -xHost"             CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fp-model fast=2"  CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -no-prec-div"       CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fpp"               CACHE STRING "Flags used by the compiler during all build types." FORCE)
    
#    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -mkl"               CACHE STRING "Flags used by the compiler during all build types." FORCE)
        
    
    #set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} --opt-report=3"  CACHE STRING "Flags used by the compiler during all build types." FORCE)
    
    # Defintion des flags release 
    set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -nogen-interfaces -zero" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    
    # Verification de la possibilite de vectorisation
    include(CheckCCompilerFlag)
    CHECK_C_COMPILER_FLAG("-axAVX,SSE4.2" HAVE_VEC)
    if (HAVE_VEC)
      # ax option generates many versions of possible optmimizations depending on target
      set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -axAVX,SSE4.2" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    endif ()
    
    
    # Defintion des flags debug   
    #    set (CMAKE_Fortran_FLAGS_DEBUG           "-O0 -g -check all -ftrapuv -gen-interfaces" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_Fortran_FLAGS_DEBUG           "-O0 -g -check all -ftrapuv -nogen-interfaces" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    
    # Defintion des flags profiling
    set (CMAKE_Fortran_FLAGS_PROFILING       "${CMAKE_Fortran_FLAGS_RELEASE} -pg" CACHE STRING "Flags used by the compiler during profiling builds." FORCE)
    
    # Defintion des flags release avec info debug
    set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "${CMAKE_Fortran_FLAGS_RELEASE} -g" CACHE STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)

    # Defintion des flags minimum release
    set (CMAKE_Fortran_FLAGS_MINSIZEREL      "-O2 -g" CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)

  endif () 

endif ()


#------------------------------------------------------------------------------
# Flags de compilation C par defaut suivant le compilateur
#------------------------------------------------------------------------------

if (C_COMPILER_NAME MATCHES "icc.*")

  if (BASESNODALES_ENABLE_DEFAULT_FLAGS)
    set (CMAKE_C_FLAGS "-std=c99 -no-multibyte-chars" CACHE STRING "Flags used by the compiler during all build types." FORCE)
    
    #libmesh7 (provisoire)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DF77API "     CACHE STRING "Flags used by the compiler during all build types." FORCE) 
    
    set (CMAKE_C_FLAGS_RELEASE "-O3" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    
    if (HAVE_VEC)
      # ax option generates many versions of possible optmimizations depending on target
      set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -axAVX,SSE4.2" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    endif ()
    #set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -vec-report1" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set (CMAKE_C_FLAGS_DEBUG           "-O0 -g -Wcheck" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_C_FLAGS_PROFILING       "-O3 -pg" CACHE STRING "Flags used by the compiler during profiling builds." FORCE)
    set (CMAKE_C_FLAGS_RELWITHDEBINFO  "-O3 -g" CACHE STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
    set (CMAKE_C_FLAGS_MINSIZEREL      "-O2 -g" CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
  endif()

endif ()



find_package(BLAS)
find_package(LAPACK)


#------------------------------------------------------------------------------
# Ajout des sources de tous les elements du projet
# En fonction des top_targets a installer
#------------------------------------------------------------------------------

add_subdirectory(src)

