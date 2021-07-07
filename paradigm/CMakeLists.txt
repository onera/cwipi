#==============================================================================
#
#             Top level CMakeLists.txt file for ParaDiGM
#
#==============================================================================

cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

#------------------------------------------------------------------------------
# Internal CMake module (shared with ParaDiGMA or CWIPI)
#------------------------------------------------------------------------------

if (NOT (PDM_IN_PDMA OR PDM_IN_CWP))
   set(PDM_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
   set(CMAKE_MODULE_PATH "${PDM_CMAKE_DIR}/modules")
endif()

include(${CMAKE_MODULE_PATH}/version.cmake)
project(PDM VERSION ${PDM_DEF_VERSION_MAJOR} LANGUAGES C CXX)

if (NOT (PDM_IN_PDMA))
   set(PDM_VERSION "${PDM_VERSION_MAJOR}.${PDM_VERSION_MINOR}.${PDM_VERSION_PATCH}")
   set(PDM_IN_PDMA_BOOL 0)
endif()

# ----------------------------------------------------------------------
find_program( MEMORYCHECK_COMMAND valgrind )
set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" )
# ----------------------------------------------------------------------

include(CTest)
enable_testing ()

#------------------------------------------------------------------------------
# User option Definition
#------------------------------------------------------------------------------

if (NOT (PDM_IN_PDMA OR PDM_IN_CWP))

   # OpenMP
   # ------

   option(PDM_ENABLE_OPENMP "Enable OpenMP" OFF)

   # Enable long for absolute number
   # -------------------------------

   option(PDM_ENABLE_LONG_G_NUM "Enable long for global number" OFF)

   # Pt-scotch
   # ---------

   option(PDM_ENABLE_PTSCOTCH "Compile with support for SCOTCH." ON)

   # ParMetis
   # --------

   option(PDM_ENABLE_PARMETIS "Compile with support for SCOTCH." ON)

   # Shared libraries
   # ----------------

   option(PDM_ENABLE_SHARED "Build Shared Libraries" OFF)

   # Static libraries
   # ----------------

   option(PDM_ENABLE_STATIC "Build Static Libraries" ON)

   # Check MPI Wrapper
   # -----------------

   option(PDM_ENABLE_MPI_CHECK "Check MPI Wrapper" ON)

   # Python Bindings
   # ---------------

   option(PDM_ENABLE_PYTHON_BINDINGS "Build Python module" OFF)

   # Fortran interface
   # -----------------

   option(PDM_ENABLE_Fortran "Build Fortran Interface" ON)

   if (PDM_ENABLE_Fortran)
     enable_language(Fortran)
   endif()

   if (PDM_ENABLE_STATIC AND PDM_ENABLE_Fortran)
      set(PDM_ENABLE_STATIC_Fortran_INTERFACE ON)
   endif ()

   if (PDM_ENABLE_SHARED AND PDM_ENABLE_Fortran)
      set(PDM_ENABLE_SHARED_Fortran_INTERFACE ON)
   endif ()

   option(PDM_ENABLE_Fortran_MPI_MODULE "Use Fortran MPI module" ON)

   # Unit test mpi with doctect
   # ------
   option(PDM_ENABLE_UNIT_TEST "Enable unit test " OFF)

elseif (PDM_IN_CWP)

  # OpenMP
  # ------

  set(PDM_ENABLE_OPENMP ${CWP_ENABLE_OPENMP})

  # Enable long for absolute number
  # -------------------------------

  set(PDM_ENABLE_LONG_G_NUM ON)

  # Pt-scotch
  # ---------

  set(PDM_ENABLE_PTSCOTCH ${CWP_ENABLE_PTSCOTCH})

  # ParMetis
  # --------

  set(PDM_ENABLE_PARMETIS ${CWP_ENABLE_PARMETIS})

  # Shared libraries
  # ----------------

  set(PDM_ENABLE_SHARED ${CWP_ENABLE_SHARED})

  # Static libraries
  # ----------------

  set(PDM_ENABLE_STATIC ${CWP_ENABLE_STATIC})

  # Check MPI Wrapper
  # -----------------

  #set(PDM_ENABLE_MPI_CHECK ${CWP_ENABLE_MPI_CHECK})

  # Python Bindings
  # ---------------

  #set(PDM_ENABLE_PYTHON_BINDINGS ${CWP_ENABLE_PYTHON_BINDINGS})

  # Fortran interface
  # -----------------

  set(PDM_ENABLE_SHARED_Fortran ${CWP_ENABLE_Fortran})

  set(PDM_ENABLE_SHARED_Fortran_INTERFACE ${CWP_ENABLE_SHARED_Fortran_LIB})

  set(PDM_ENABLE_STATIC_Fortran_INTERFACE ${CWP_ENABLE_STATIC_Fortran_LIB})

  option(PDM_ENABLE_Fortran_MPI_MODULE "Use Fortran MPI module" ON)

endif ()

#------------------------------------------------------------------------------
# Default build : Release
#------------------------------------------------------------------------------

if (NOT PDM_IN_PDMA)

    if ((PDM_ENABLE_PYTHON_BINDINGS OR PDM_ENABLE_SHARED_Fortran_INTERFACE) AND (NOT PDM_ENABLE_SHARED))
      set (PDM_ENABLE_SHARED "ON" CACHE STRING "Build Shared Libraries" FORCE )
    endif()

    if ((NOT PDM_ENABLE_STATIC) AND (NOT PDM_ENABLE_SHARED))
       message (FATAL_ERROR "No enabled library. Please, enable shared or static library")
    endif()

    if (NOT CMAKE_BUILD_TYPE)
       set(CMAKE_BUILD_TYPE "Release" CACHE STRING  "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel" FORCE)
    endif()

    # ADD "-DDEBUG_CLASSE" in Debug
    #------------------------------

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
       add_definitions (-DDEBUG_CLASSE)
    endif()

    # Add "-DHAVE_GETRUSAGE" if function getrusage and gettimeofday exist
    #--------------------------------------------------------------------

    include(CheckFunctionExists)
    CHECK_FUNCTION_EXISTS(getrusage PDM_HAVE_GETRUSAGE)
    CHECK_FUNCTION_EXISTS(gettimeofday PDM_HAVE_GETTIMEOFDAY)

endif ()

#------------------------------------------------------------------------------
# Default compiler flags
#------------------------------------------------------------------------------

if (NOT (PDM_IN_PDMA OR PDM_IN_CWP))

   include(${CMAKE_MODULE_PATH}/default_flags.cmake)

endif ()

#------------------------------------------------------------------------------
# Libraries to add to link
#------------------------------------------------------------------------------

if (NOT (PDM_IN_PDMA OR PDM_IN_CWP))

   set(COMMON_LIBRARIES)

endif ()

#------------------------------------------------------------------------------
# Check dependecies
#------------------------------------------------------------------------------

if (NOT (PDM_IN_PDMA OR PDM_IN_CWP))

  # Check MPI (TODO : check MPI wrapper function check )
  # -------------

  if (PDM_ENABLE_Fortran)
    find_package(MPI REQUIRED COMPONENTS C Fortran)
  else()
    find_package(MPI REQUIRED COMPONENTS C)
  endif()

  if (PDM_ENABLE_Fortran_MPI_MODULE)
    set (PDM_HAVE_FORTRAN_MPI_MODULE 1)
  endif()

  mark_as_advanced (MPI_EXTRA_LIBRARY MPI_LIBRARY)

  # ---------------------
  # Check MPI Wrapper
  # ---------------------

  if (PDM_ENABLE_MPI_CHECK)
     include(${CMAKE_MODULE_PATH}/CheckMPI.cmake)
  endif()

  # Checck SCOTCH an ParMETIS
  # -------------------------

  if (PDM_ENABLE_PARMETIS)
     find_package(ParMETIS 4.0.3)
     if (PARMETIS_FOUND)
        set(PDM_HAVE_PARMETIS 1)
     endif()
  endif()

  if (PDM_ENABLE_PTSCOTCH)
     find_package(PTSCOTCH 6.0.0)
     if (PTSCOTCH_FOUND)
        set(PDM_HAVE_PTSCOTCH 1)
     endif()
  endif()

  # Check OpenMP (TODO : check MPI wrapper function check )
  # -------------------------------------------------------

  if (PDM_ENABLE_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
      set(PDM_HAVE_OPENMP 1)
    endif()
  endif()

  # Python bindings
  # ---------------

  if (PDM_ENABLE_PYTHON_BINDINGS)
    # NB: cython module must use mpicxx to avoid "mpi.h" not found.
    find_package(MPI REQUIRED COMPONENTS CXX)

    find_package(PythonInterp REQUIRED)
    find_package(PythonLibs REQUIRED)
    find_package(NumPy REQUIRED)
    find_package(Mpi4Py REQUIRED)
    include(UseCython)

    if (NOT CMAKE_PDM_INSTALL_PYTHON_DIR)

       set (CMAKE_PDM_INSTALL_PYTHON_DIR  "${CMAKE_INSTALL_PREFIX}")
       set (CMAKE_PDM_INSTALL_PYTHON_DIR
            ${CMAKE_PDM_INSTALL_PYTHON_DIR} CACHE STRING "Install path directory for Python packages")

    endif ()

  endif()

elseif (PDM_IN_CWP)

  # Checck SCOTCH an ParMETIS
  # -------------------------

  if (PDM_ENABLE_PARMETIS)
     if (PARMETIS_FOUND)
        set(PDM_HAVE_PARMETIS 1)
     endif()
  endif()

  if (PDM_ENABLE_PTSCOTCH)
     if (PTSCOTCH_FOUND)
        set(PDM_HAVE_PTSCOTCH 1)
     endif()
  endif()

  if (PDM_ENABLE_Fortran_MPI_MODULE)
    set (PDM_HAVE_FORTRAN_MPI_MODULE 1)
  endif()


endif ()

#------------------------------------------------------------------------------
# Store variable in configure file
#------------------------------------------------------------------------------

if (NOT (PDM_IN_PDMA))

  if (PDM_ENABLE_LONG_G_NUM)
     set(PDM_LONG_G_NUM 1)
  endif ()

  configure_file(pdm_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h")

  if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)
    if ("${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h" IS_NEWER_THAN  "${CMAKE_CURRENT_BINARY_DIR}/pdm_configf.h")
     include (FortranizeFile)
     fortranize ("${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h")
    endif ()
  endif ()

endif ()

#------------------------------------------------------------------------------
# Doc
#------------------------------------------------------------------------------

if (NOT (PDM_IN_PDMA OR PDM_IN_CWP))

   find_package(Doxygen)

endif ()

#------------------------------------------------------------------------------
# Print summary of found and not found optional packages
#------------------------------------------------------------------------------

if (NOT (PDM_IN_PDMA OR PDM_IN_CWP))

   include(FeatureSummary)
   FEATURE_SUMMARY(WHAT ALL)

endif()

#------------------------------------------------------------------------------
# Installation de PDM
#------------------------------------------------------------------------------

# Append the library version information to the library target properties
if (PDM_WITH_LIBRARY_VERSION)
  string(REPLACE "+" "" PDM_LIBRARY_VERSION ${PDM_VERSION})
  # This setting of SOVERSION assumes that any API change
  # will increment either the minor or major version number.
  set(PDM_LIBRARY_PROPERTIES ${PDM_LIBRARY_PROPERTIES}
    VERSION ${PDM_LIBRARY_VERSION}
    SOVERSION ${PDM_VERSION_MAJOR}.${PDM_VERSION_MINOR}
     )
endif()

# Set PDM install sub-directories
set(PDM_BIN_DIR "bin")
set(PDM_LIB_DIR "lib")
set(PDM_INCLUDE_DIR "include")
set(PDM_PKGCONFIG_DIR "lib/pkgconfig")
set(PDM_SHARE_DIR "share/cedre")
set(PDM_MAN_DIR "share/man")
set(PDM_DOC_DIR "${PDM_SHARE_DIR}/doc")
set(PDM_ETC_DIR "etc")

#------------------------------------------------------------------------------
# Ajout des sources de tous les elements du projet
# En fonction des top_targets a installer
#------------------------------------------------------------------------------

add_subdirectory(src)
if (PDM_ENABLE_PYTHON_BINDINGS)
  add_subdirectory(Cython)
endif()
add_subdirectory(test)

#------------------------------------------------------------------------------
if(PDM_ENABLE_UNIT_TEST)
  include(FetchContent)
  # FetchContent_Declare(
  #   doctest
  #   GIT_REPOSITORY https://github.com/onqtam/doctest.git
  #   GIT_TAG        2.4.1
  # )
  FetchContent_Declare(
    doctest
     GIT_REPOSITORY https://github.com/BerengerBerthoul/doctest.git
     GIT_TAG        f13a00607463ce211004be23b5aa0bdf24d3ab61
  )
  FetchContent_GetProperties(doctest)
  if(NOT doctest_POPULATED)
    FetchContent_Populate(doctest)
    add_subdirectory(${doctest_SOURCE_DIR} ${doctest_BINARY_DIR})
  endif()
endif()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Création d'un package pdmConfig.cmake pour importer facilement paradigm
# pour les clients utilisant CMake
#------------------------------------------------------------------------------
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/templates/pdmConfig.cmake.in"
  pdmConfig.cmake
  @ONLY
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pdmConfig.cmake"
  DESTINATION lib/cmake
)

#------------------------------------------------------------------------------
# Affichage du message de fin d'installation
#------------------------------------------------------------------------------

if (NOT (PDM_IN_PDMA OR PDM_IN_CWP))

add_subdirectory(${PDM_CMAKE_DIR}/post-install)

endif()
