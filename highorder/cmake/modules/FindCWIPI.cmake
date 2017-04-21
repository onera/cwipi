# - Try to find CWIPI    
# Once done this will define
#
#  CWIPI_FOUND        - system has found CWIPI
#  CWIPI_INCLUDE_DIR - include directories for CWIPI
#  CWIPI_LIBRARIES    - libraries for CWIPI
#  CWIPI_VERSION      - version for CWIPI
#

set(CWIPI_FOUND FALSE)

# Check for header file
find_path(CWIPI_INCLUDE_DIR cwipi_config.h  cwipi_f.f90  cwipi.h  cwipi.mod
  HINTS ${CWIPI_DIR}/include $ENV{CWIPI_DIR}/include
  DOC "Directory where the CWIPI header is located"
)

# Check for cwipi_f
find_library(CWIPIF_LIBRARY
  NAMES cwipi_f 
  HINTS ${CWIPI_DIR}/lib $ENV{CWIPI_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The CWIPI library"
  )
find_library(CWIPIF_LIBRARY
  NAMES cwipi_f 
  DOC "The CWIPI library"
  )

# Check for cwipi
find_library(CWIPI_LIBRARY
  NAMES cwipi 
  HINTS ${CWIPI_DIR}/lib $ENV{CWIPI_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The CWIPI library"
  )
find_library(CWIPI_LIBRARY
  NAMES cwipi 
  DOC "The CWIPI library"
  )

# Check for fvmc
find_library(FVMC_LIBRARY
  NAMES fvmc 
  HINTS ${CWIPI_DIR}/lib $ENV{CWIPI_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The CWIPI library"
  )
find_library(FVMC_LIBRARY
  NAMES fvmc 
  DOC "The CWIPI library"
  )

# Check for bftc
find_library(BFTC_LIBRARY
  NAMES bftc 
  HINTS ${CWIPI_DIR}/lib $ENV{CWIPI_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The CWIPI library"
  )
find_library(BFTC_LIBRARY
  NAMES bftc 
  DOC "The CWIPI library"
  )

# Check for cwipi-config
find_program(CWIPI_CONFIG_EXE
  NAMES cwipi-config    
  HINTS ${CWIPI_DIR}/bin $ENV{CWIPI_DIR}/bin
  NO_DEFAULT_PATH
  DOC "The CWIPI config executable"
  )
find_program(CWIPI_CONFIG_EXE
  NAMES cwipi-config    
  DOC "The CWIPI config executable"
  )

set(CWIPI_LIBRARIES ${CWIPIF_LIBRARY} CACHE STRING "CWIPI libraries")
set(CWIPI_LIBRARIES ${CWIPI_LIBRARIES} ${CWIPI_LIBRARY} CACHE STRING "CWIPI libraries")
set(CWIPI_LIBRARIES ${CWIPI_LIBRARIES} ${FVMC_LIBRARY} CACHE STRING "CWIPI libraries")
set(CWIPI_LIBRARIES ${CWIPI_LIBRARIES} ${BFTC_LIBRARY} CACHE STRING "CWIPI libraries")

if (CWIPI_INCLUDE_DIR AND CWIPI_LIBRARIES AND CWIPI_CONFIG_EXE)

    execute_process(COMMAND ${CWIPI_CONFIG_EXE} --version
                    OUTPUT_VARIABLE OUTPUT
                    ERROR_VARIABLE ERROR_VAR)

    if (OUTPUT)
      set(CWIPI_VERSION ${OUTPUT} CACHE TYPE STRING)
      mark_as_advanced(CWIPI_VERSION)
    endif()

    # Check if version found is >= required version

    if (CWIPI_FIND_VERSION)
      if (NOT "${CWIPI_VERSION}" VERSION_LESS "${CWIPI_FIND_VERSION}")
        set(CWIPI_VERSION_OK TRUE)
       endif()
    else()
   # No specific version requested
      set(CWIPI_VERSION_OK TRUE)
     endif()
endif()

mark_as_advanced(CWIPI_VERSION_OK)

# Standard package handling
find_package_handle_standard_args(CWIPI
                                  "CWIPI could not be found. Be sure to set CWIPI_DIR."
                                  CWIPI_LIBRARIES
                                  CWIPI_INCLUDE_DIR
                                  CWIPI_CONFIG_EXE
                                  CWIPI_VERSION
                                  CWIPI_VERSION_OK)

mark_as_advanced(CWIPI_LIBRARIES
                 CWIPI_INCLUDE_DIR
                 CWIPI_CONFIG_EXE
                 CWIPI_VERSION
                 CWIPI_VERSION_OK)
unset(CWIPIF_LIBRARY CACHE)
unset(CWIPI_LIBRARY CACHE)
unset(BFTC_LIBRARY CACHE)
unset(FVMC_LIBRARY CACHE)
