# - Try to find SCOTCH
# Once done this will define
#
#  SCOTCH_FOUND        - system has found SCOTCH
#  SCOTCH_INCLUDE_DIRS - include directories for SCOTCH
#  SCOTCH_LIBARIES     - libraries for SCOTCH
#  SCOTCH_VERSION      - version for SCOTCH
#

# mettre un if si on ne peut pas compiler et mettre un message

set(SCOTCH_FOUND FALSE)

# Check for header file
find_path(SCOTCH_INCLUDE_DIRS ptscotch.h
  HINTS ${SCOTCH_DIR}/include $ENV{SCOTCH_DIR}/include
  DOC "Directory where the Pt-SCOTCH header is located"
  )

# Check for scotch
find_library(SCOTCH_LIBRARY
  NAMES scotch
  HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The SCOTCH library"
  )
find_library(SCOTCH_LIBRARY
  NAMES scotch
  DOC "The SCOTCH library"
  )

# Check for scotcherr
find_library(SCOTCHERR_LIBRARY
  NAMES scotcherr
  HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The SCOTCH-ERROR library"
  )
find_library(SCOTCHERR_LIBRARY
  NAMES scotcherr
  DOC "The SCOTCH-ERROR library"
  )

# Check for ptscotch
find_library(PTSCOTCH_LIBRARY
  NAMES ptscotch
  HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The Pt-SCOTCH library"
  )
find_library(PTSCOTCH_LIBRARY
  NAMES ptscotch
  DOC "The Pt-SCOTCH library"
  )

# Check for ptscotcherr
find_library(PTSCOTCHERR_LIBRARY
  NAMES ptscotcherr
  HINTS ${SCOTCH_DIR}/lib $ENV{SCOTCH_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The Pt-SCOTCH-ERROR library"
  )
find_library(PTSCOTCHERR_LIBRARY
  NAMES ptscotcherr
  DOC "The Pt-SCOTCH-ERROR library"
  )

set(SCOTCH_LIBRARIES ${PTSCOTCH_LIBRARY} ${PTSCOTCHERR_LIBRARY} ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY} CACHE STRING "SOCTCH parallel libraries")
set(SCOTCH_SEQ_LIBRARIES ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY} CACHE STRING "SOCTCH sequential libraries")

if (SCOTCH_INCLUDE_DIRS AND SCOTCH_LIBRARIES)

  mark_as_advanced(SCOTCH_INCLUDE_DIRS SCOTCH_LIBRARIES SCOTCH_SEQ_LIBRARIES)

  if (PTSCOTCH_LIBRARY)
    string(REGEX REPLACE "(^.*)/lib/libptscotch.*$" "\\1" SCOTCH_LIBRARY_PATH ${PTSCOTCH_LIBRARY} )
  endif (PTSCOTCH_LIBRARY)

  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES ${SCOTCH_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${SCOTCH_LIBRARIES})

  # Add MPI variables if MPI has been found
  if (MPI_C_LIBRARIES)
    set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_C_INCLUDE_PATH})
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${MPI_C_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${MPI_C_COMPILE_FLAGS}")
  endif()
    
  set(CMAKE_REQUIRED_FLAGS  "${CMAKE_REQUIRED_FLAGS} ${CMAKE_C_FLAGS}")

  set(SCOTCH_CONFIG_TEST_VERSION_C
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/scotch_config_test_version.c")
 
  file(WRITE ${SCOTCH_CONFIG_TEST_VERSION_C} "
#include <stdint.h>
#include <stdio.h>
#include <mpi.h>
#include <ptscotch.h>

int main() {
  printf(\"%i.%i.%i\\n\", SCOTCH_VERSION,
	                  SCOTCH_RELEASE,
	                  SCOTCH_PATCHLEVEL);
  return 0;
}
")

try_run(
    SCOTCH_CONFIG_TEST_VERSION_EXITCODE
    SCOTCH_CONFIG_TEST_VERSION_COMPILED
    ${CMAKE_CURRENT_BINARY_DIR}
    ${SCOTCH_CONFIG_TEST_VERSION_C}
    CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
      "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}"
    COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE OUTPUT
)


if (   (NOT SCOTCH_CONFIG_TEST_VERSION_COMPILED)
    OR (NOT (SCOTCH_CONFIG_TEST_VERSION_EXITCODE EQUAL 0))) 
    message(WARNING "Unable to determine SCOTCH version")
    set(SCOTCH_VERSION_OK TRUE)
    set(SCOTCH_VERSION "??.??.??" CACHE TYPE STRING)
    set(SCOTCH_TEST_COMPILE TRUE)
else ()
    set(SCOTCH_VERSION ${OUTPUT} CACHE TYPE STRING)
    mark_as_advanced(SCOTCH_VERSION)

   if (SCOTCH_FIND_VERSION)
  # Check if version found is >= required version
      if (NOT "${SCOTCH_VERSION}" VERSION_LESS "${SCOTCH_FIND_VERSION}")
  	set(SCOTCH_VERSION_OK TRUE)
      endif()
   else()
  # No specific version requested
       set(SCOTCH_VERSION_OK TRUE)
   endif()
   mark_as_advanced(SCOTCH_VERSION_OK)
endif ()

unset(PTSCOTCH_LIBRARY CACHE)
unset(PTSCOTCHERR_LIBRARY CACHE)
unset(SCOTCH_LIBRARY CACHE)
unset(SCOTCHERR_LIBRARY CACHE)
unset(SCOTCH_LIBRARY_PATH CACHE)

endif ()
#
# Standard package handling
#
find_package_handle_standard_args(SCOTCH
                                  "SCOTCH could not be found. Be sure to set SCOTCH_DIR."
                                  SCOTCH_LIBRARIES
                                  SCOTCH_SEQ_LIBRARIES
                                  SCOTCH_INCLUDE_DIRS
                                  SCOTCH_VERSION
                                  SCOTCH_VERSION_OK)
