# - Try to find FVM    
# Once done this will define
#
#  FVM_FOUND        - system has found FVM
#  FVM_INCLUDE_DIRS - include directories for FVM
#  FVM_LIBARIES     - libraries for FVM
#  FVM_VERSION      - version for FVM
#

set(FVM_FOUND FALSE)

# Check for header file
find_path(FVM_INCLUDE_DIRS fvm_config.h fvm_point_location.h
  HINTS ${FVM_DIR}/include $ENV{FVM_DIR}/include
  DOC "Directory where the FVM header is located"
)

# Check for fvm
if(FVM_USE_STATIC_LIBRARIES)
  find_library(FVM_LIBRARY
    NAMES libfvm_scal.a
    HINTS ${FVM_DIR}/lib $ENV{FVM_DIR}/lib
    NO_DEFAULT_PATH
    DOC "The FVM library"
    )
  find_library(FVM_LIBRARY
    NAMES fvm_scal 
    DOC "The FVM library"
    )
else()
  find_library(FVM_LIBRARY
    NAMES fvm_scal 
    HINTS ${FVM_DIR}/lib $ENV{FVM_DIR}/lib
    NO_DEFAULT_PATH
    DOC "The FVM library"
    )
  find_library(FVM_LIBRARY
    NAMES fvm_scal 
    DOC "The FVM library"
    )
endif()

# Check for fvm-config
find_program(FVM_CONFIG_EXE
  NAMES fvm-config    
  HINTS ${FVM_DIR}/bin $ENV{FVM_DIR}/bin
  NO_DEFAULT_PATH
  DOC "The FVM config executable"
  )
find_program(FVM_CONFIG_EXE
  NAMES fvm-config    
  DOC "The FVM config executable"
  )

set(FVM_LIBRARIES ${FVM_LIBRARY} CACHE STRING "FVM libraries")

if (FVM_INCLUDE_DIRS AND FVM_LIBRARIES AND FVM_CONFIG_EXE)

    execute_process(COMMAND ${FVM_CONFIG_EXE} --version
                    OUTPUT_VARIABLE OUTPUT)

    if (OUTPUT)
      set(FVM_VERSION ${OUTPUT} CACHE TYPE STRING)
      mark_as_advanced(FVM_VERSION)
    endif()

    # Check if version found is >= required version

    if (FVM_FIND_VERSION)
      if (NOT "${FVM_VERSION}" VERSION_LESS "${FVM_FIND_VERSION}")
        set(FVM_VERSION_OK TRUE)
       endif()
    else()
   # No specific version requested
      set(FVM_VERSION_OK TRUE)
     endif()
endif()

mark_as_advanced(FVM_VERSION_OK)

# Standard package handling
find_package_handle_standard_args(FVM
                                  "FVM could not be found. Be sure to set FVM_DIR."
                                  FVM_LIBRARIES
                                  FVM_INCLUDE_DIRS
                                  FVM_CONFIG_EXE
                                  FVM_VERSION
                                  FVM_VERSION_OK)

mark_as_advanced(FVM_LIBRARIES 
                 FVM_INCLUDE_DIRS
                 FVM_CONFIG_EXE
                 FVM_VERSION
                 FVM_VERSION_OK)
unset(FVM_LIBRARY CACHE)