# Top level CMakeLists.txt file for PTK

# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

project(BNOD C CXX Fortran)

message("Source du projet: ${BNOD_SOURCE_DIR}")
message("Compilation du projet: ${BNOD_BINARY_DIR}")

#------------------------------------------------------------------------------
# Definition de la version BNOD
#------------------------------------------------------------------------------
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

# Defintion de la version pour le futur
# CMake pilotera l'ensemble via la generation de fichier de configuration
# Ces fichiers seront ensuite inclus dans les sources de BNOD
set(BNOD_VERSION_MAJOR "2015")
set(BNOD_VERSION_MINOR "0")
set(BNOD_VERSION_MICRO "0")
set(BNOD_VERSION "${BNOD_VERSION_MAJOR}.${BNOD_VERSION_MINOR}.${BNOD_VERSION_MICRO}")

message("Version du projet: ${BNOD_VERSION}")

#------------------------------------------------------------------------------
# Utilisation des flags de compilation par dÈfaut du projet
# (Flags par dÈfaut sont dÈfinis plus bas dans le fichier) 
#------------------------------------------------------------------------------

option(BNOD_ENABLE_DEFAULT_FLAGS "Enable default flags" ON)


#------------------------------------------------------------------------------
# Localisation des fonctions cmake du projet
#------------------------------------------------------------------------------

set(BNOD_CMAKE_DIR "${BNOD_SOURCE_DIR}/cmake" CACHE INTERNAL "")
set(CMAKE_MODULE_PATH "${BNOD_CMAKE_DIR}/modules")


#------------------------------------------------------------------------------
# Activation des langages Fortran C et C++
#------------------------------------------------------------------------------

enable_language(Fortran)
enable_language(C)
enable_language(CXX)

#------------------------------------------------------------------------------
# Construction par defaut : mode release
#------------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()


# Ajout de "-DDEBUG_CLASSE" em mode DEBUG
#----------------------------------------

#if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#  add_definitions (-DDEBUG_CLASSE)
#endif()


#------------------------------------------------------------------------------
# Flags de compilation Fortran par defaut suivant le compilateur
#------------------------------------------------------------------------------

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
get_filename_component (C_COMPILER_NAME       ${CMAKE_C_COMPILER}       NAME)
get_filename_component (CXX_COMPILER_NAME     ${CMAKE_CXX_COMPILER}     NAME)


if (Fortran_COMPILER_NAME MATCHES "gfortran.*")

  # gfortran
  # --------

  if (BNOD_ENABLE_DEFAULT_FLAGS)

    set (CMAKE_Fortran_FLAGS "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none" CACHE STRING "Flags used by the compiler during all build types." FORCE)
    
    # Activation du traceback  
    
    find_package(Libunwind)
    if (NOT DVD AND LIBUNWIND_FOUND)
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace" CACHE STRING "Flags used by the compiler during all build types." FORCE)
    endif (NOT DVD AND LIBUNWIND_FOUND)
    
    set (CMAKE_Fortran_FLAGS_RELEASE         "-O3"     CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set (CMAKE_Fortran_FLAGS_DEBUG           "-O0 -g"  CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_Fortran_FLAGS_PROFILING       "-O3 -pg" CACHE STRING "Flags used by the compiler during profiling builds." FORCE)
    set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "-O3 -g"  CACHE STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
    set (CMAKE_Fortran_FLAGS_MINSIZEREL      "-O2 -g"  CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)

  endif()

  set (FORTRAN_LIBRARIES                   )
  set (FORTRAN_LIBRARIES_FLAG              )
  set (COMPIL_SPEC_SCRIPT  "compil_spec"   )

elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")

  # ifort
  # -----
  message("je passe ici")

  if (BNOD_ENABLE_DEFAULT_FLAGS)

    message("puis passe ici")
    
    # Ajout des flags commun 

    set (CMAKE_Fortran_FLAGS                          "-fpp"             CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -arch x86_64"      CACHE STRING "Flags used by the compiler during all build types." FORCE)
    
    # Activation du traceback  
    find_package(Libunwind)
    if (NOT DVD AND LIBUNWIND_FOUND)
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -traceback"      CACHE STRING "Flags used by the compiler during all build types." FORCE)
    endif (NOT DVD AND LIBUNWIND_FOUND)
    
    # Defintion des flags release 
    set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -nogen-interfaces -zero" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    include(CheckCCompilerFlag)
    CHECK_C_COMPILER_FLAG("-axAVX,SSE4.2" HAVE_VEC)
    
    # Verification de la possibilite de vectorisation
    if (HAVE_VEC)
      # ax option generates many versions of possible optmimizations depending on target
      set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -axAVX,SSE4.2" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    endif ()
    
    if (DVD)
      set (CMAKE_Fortran_FLAGS_RELEASE   "${CMAKE_Fortran_FLAGS_RELEASE} -static-intel -i-static" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    endif()

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

  find_library(FORTRAN_LIBRARIES ifcore)
  set (FORTRAN_LIBRARIES_FLAG          )
  set (COMPIL_SPEC_SCRIPT  "compil_spec"   )

elseif (Fortran_COMPILER_NAME MATCHES "xlf*")

  # xlf
  # ---

  if (BNOD_ENABLE_DEFAULT_FLAGS)

    # Ajout des flags communs  


    
    set (CMAKE_Fortran_FLAGS " -qintsize=4 -qrealsize=8" CACHE STRING "Flags used by the compiler during all build types." FORCE)
    #   set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsuffix=f=f90 -qsuffix=cpp=F90" CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsuffix=f=f90 -qsuffix=cpp=F90 -q64 -qundef -qextchk -qhalt=e -qfloat=nans" CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qextname -qflag=i:i -qsuppress=1501-245" CACHE STRING "Flags used by the compiler during all build types." FORCE) 
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qinitauto=FF -qspillsize=5000 " CACHE STRING "Flags used by the compiler during all build types." FORCE)
    # 
    # Defintion des flags release 

    set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -qstrict" CACHE STRING "Flags used by the compiler during release builds." FORCE)

    # Defintion des flags debug   

    set (CMAKE_Fortran_FLAGS_DEBUG "-g -qfullpath -O0 -C -qsigtrap -qnooptimize" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -qflttrap=enable:invalid:overflow:zerodivide" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -qwarn64" CACHE STRING "Flags used by the compiler during debug builds." FORCE)

    # Defintion des flags profiling

    set (CMAKE_Fortran_FLAGS_PROFILING       "${CMAKE_Fortran_FLAGS_RELEASE} -p" CACHE STRING "Flags used by the compiler during profiling builds." FORCE)
    set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "-O3 -g" CACHE STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
    set (CMAKE_Fortran_FLAGS_MINSIZEREL      "-O2" CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)

  endif ()
  
  set(IBM_LIB stdc++ ibmc++)
  set(XL_LIB xl xlf90 xlsmp xlopt)
  link_directories(/opt/ibmcmp/xlsmp/3.1/lib64 /opt/ibmcmp/vacpp/12.1/lib64 /opt/ibmcmp/xlf/14.1/lib64)

  set (COMPIL_SPEC_SCRIPT  "compil_spec_aix"   )

elseif (Fortran_COMPILER_NAME MATCHES "pgf90*")

  # pgi
  # ---

  if (BNOD_ENABLE_DEFAULT_FLAGS)

    # Ajout des flags communs  

    set (CMAKE_Fortran_FLAGS "-i4 -r8" CACHE STRING "Flags used by the compiler during all build types." FORCE)

    set (CMAKE_Fortran_FLAGS_RELEASE         "-fastsse" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set (CMAKE_Fortran_FLAGS_DEBUG           "-g -Minform=inform -Ktrap=inv,divz,ovf" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_Fortran_FLAGS_PROFILING       "${CMAKE_Fortran_FLAGS_RELEASE} -pg" CACHE STRING "Flags used by the compiler during profiling builds." FORCE)
    set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "${CMAKE_Fortran_FLAGS_RELEASE} -g" CACHE STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
    set (CMAKE_Fortran_FLAGS_MINSIZEREL      "-O2" CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)

  endif ()

  set (FORTRAN_LIBRARIES                   )
  set (FORTRAN_LIBRARIES_FLAG    -pgf90libs)
  set (COMPIL_SPEC_SCRIPT  "compil_spec"   )

else (Fortran_COMPILER_NAME MATCHES "GNU*")

  # si on ne trouve pas de compilateur, on essaye ca ...
  # ----------------------------------------------------

  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})

  if (BNOD_ENABLE_DEFAULT_FLAGS)
    message ("No optimized Fortran compiler flags are known, we just try -O2...")

    set (CMAKE_Fortran_FLAGS "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none" CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
  endif()

  set (FORTRAN_LIBRARIES                   )
  set (FORTRAN_LIBRARIES_FLAG              )
  set (COMPIL_SPEC_SCRIPT  "compil_spec"   )

endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


#------------------------------------------------------------------------------
# Flags de compilation C par defaut suivant le compilateur
#------------------------------------------------------------------------------

if (C_COMPILER_NAME MATCHES "gcc.*")

  # gcc
  # ---

  if (BNOD_ENABLE_DEFAULT_FLAGS)
    set (CMAKE_C_FLAGS "-std=c99 -arch x86_64" CACHE STRING "Flags used by the compiler during all build types." FORCE)

    set (CMAKE_C_FLAGS_RELEASE         "-O3" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set (CMAKE_C_FLAGS_DEBUG           "-O0 -g" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_C_FLAGS_PROFILING       "-O3 -pg" CACHE STRING "Flags used by the compiler during profiling builds." FORCE)
    set (CMAKE_C_FLAGS_RELWITHDEBINFO  "-O3 -g" CACHE STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
    set (CMAKE_C_FLAGS_MINSIZEREL      "-O2 -g" CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
  endif()

elseif (C_COMPILER_NAME MATCHES "icc.*")

  # icc
  # ---

  if (BNOD_ENABLE_DEFAULT_FLAGS)
    set (CMAKE_C_FLAGS "-std=c99 -no-multibyte-chars" CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_C_FLAGS_RELEASE "-O3" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    if (HAVE_VEC)
      # ax option generates many versions of possible optmimizations depending on target
      set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -axAVX,SSE4.2" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    endif ()
    if (DVD)
      set (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE} -static-intel" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    endif()
    set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -vec-report1" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set (CMAKE_C_FLAGS_DEBUG           "-O0 -g -Wcheck" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_C_FLAGS_PROFILING       "-O3 -pg" CACHE STRING "Flags used by the compiler during profiling builds." FORCE)
    set (CMAKE_C_FLAGS_RELWITHDEBINFO  "-O3 -g" CACHE STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
    set (CMAKE_C_FLAGS_MINSIZEREL      "-O2 -g" CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
  endif()

elseif (C_COMPILER_NAME MATCHES "xlc*")

  # xlc
  # ---

  if (BNOD_ENABLE_DEFAULT_FLAGS)
    set (CMAKE_C_FLAGS "-q64 -qflag=i:i -qsuppress=1501-245 -qlanglvl=stdc99" CACHE STRING "Flags used by the compiler during all build types." FORCE)
  #  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qextname=cedre_io_get_timer_distrib,cedre_io_ecriture_globale -qflag=i:i -qsuppress=1501-245" CACHE STRING "Flags used by the compiler during all build types." FORCE) 


    set (CMAKE_C_FLAGS_RELEASE         "-O3 -qstrict" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set (CMAKE_C_FLAGS_DEBUG           "-g" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_C_FLAGS_PROFILING       "${CMAKE_C_FLAGS_RELEASE} -pg" CACHE STRING "Flags used by the compiler during profiling builds." FORCE)
    set (CMAKE_C_FLAGS_RELWITHDEBINFO  "${CMAKE_C_FLAGS_RELEASE} -g" CACHE STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
    set (CMAKE_C_FLAGS_MINSIZEREL      "-O2 -g" CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
  endif()

elseif (C_COMPILER_NAME MATCHES "pgcc.*")

  # pgcc
  # ----


  if (BNOD_ENABLE_DEFAULT_FLAGS)
    set (CMAKE_C_FLAGS "-c99" CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_C_FLAGS_RELEASE         "-fastsse" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set (CMAKE_C_FLAGS_DEBUG           "-g" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_C_FLAGS_PROFILING       "${CMAKE_C_FLAGS_RELEASE} -pg" CACHE STRING "Flags used by the compiler during profiling builds." FORCE)
    set (CMAKE_C_FLAGS_RELWITHDEBINFO  "${CMAKE_C_FLAGS_RELEASE} -g" CACHE STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
    set (CMAKE_C_FLAGS_MINSIZEREL      "-O2 -g" CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
  endif()

else (C_COMPILER_NAME MATCHES "gcc.*")

  # si on ne trouve pas de compilateur, on essaye ca ...
  # ----------------------------------------------------

  message ("CMAKE_C_COMPILER full path: " ${CMAKE_C_COMPILER})
  message ("C compiler: " ${C_COMPILER_NAME})
  if (BNOD_ENABLE_DEFAULT_FLAGS)
    message ("No optimized C compiler flags are known, we just try -O2...")

    set (CMAKE_C_FLAGS "-std=c99" CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set (CMAKE_C_FLAGS_RELEASE         "-O3" CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set (CMAKE_C_FLAGS_DEBUG           "-O0 -g" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set (CMAKE_C_FLAGS_PROFILING       "-O3 -pg" CACHE STRING "Flags used by the compiler during profiling builds." FORCE)
    set (CMAKE_C_FLAGS_RELWITHDEBINFO  "-O3 -g" CACHE STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
    set (CMAKE_C_FLAGS_MINSIZEREL      "-O2 -g" CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
  endif()

endif (C_COMPILER_NAME MATCHES "gcc.*")

mark_as_advanced(CXX_LIBRARIES CXX_LIBRARIES_FLAG FORTRAN_LIBRARIES FORTRAN_LIBRARIES_FLAG)

if (BNOD_ENABLE_DEFAULT_FLAGS)
  mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING CMAKE_C_FLAGS_PROFILING CMAKE_Fortran_FLAGS_PROFILING)
endif()


#------------------------------------------------------------------------------
# Definition des bibliothèques communes à ajouter au link
#------------------------------------------------------------------------------


#find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)
#find_package(MKL REQUIRED)

find_package(MKL REQUIRED)

#include_directories(${MKL_INCLUDE_DIRS})
#link_directories(${MKL_LIBRARIES})
#target_link_libraries(<module>
#mkl_sequential
#)


#set(COMMON_LIBRARIES)
#set(COMMON_LIBRARIES -L/opt/intel/composer_xe_2013.0.088/mkl/lib)

# Pour libunwind : Ecriture du flag sous la forme -lunwind 
#(sinon erreur avec ifort si on met directement ${LIBUNWIND_LIBRARIES_PATH})
# --------------------------------------------------------------------------

if (NOT DVD AND LIBUNWIND_FOUND)
  set(COMMON_LIBRARIES "${COMMON_LIBRARIES} -L${LIBUNWIND_LIBRARY_PATH} -lunwind")
endif (NOT DVD AND LIBUNWIND_FOUND)



#------------------------------------------------------------------------------
# Ajout des sources de tous les elements du projet
# En fonction des top_targets a installer
#------------------------------------------------------------------------------

add_subdirectory(src)

