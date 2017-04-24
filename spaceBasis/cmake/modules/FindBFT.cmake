# - Try to find BFT    
# Once done this will define
#
#  BFT_FOUND        - system has found BFT
#  BFT_INCLUDE_DIRS - include directories for BFT
#  BFT_LIBARIES     - libraries for BFT
#  BFT_VERSION      - version for BFT
#

set(BFT_FOUND FALSE)

# Check for header file
find_path(BFT_INCLUDE_DIRS bft_backtrace.h  bft_config.h  bft_config_priv.h  bft_error.h  bft_file.h  bft_fp_trap.h  bft_mem.h  bft_mem_usage.h  bft_printf.h  bft_sys_info.h  bft_timer.h  bft_version.h
  HINTS ${BFT_DIR}/include $ENV{BFT_DIR}/include
  DOC "Directory where the BFT header is located"
)

# Check for bft
if(BFT_USE_STATIC_LIBRARIES)
  find_library(BFT_LIBRARY
    NAMES libbft.a    
    HINTS ${BFT_DIR}/lib $ENV{BFT_DIR}/lib
    NO_DEFAULT_PATH
    DOC "The BFT library"
    )
  find_library(BFT_LIBRARY
    NAMES bft    
    DOC "The BFT library"
    )
else()
  find_library(BFT_LIBRARY
    NAMES bft    
    HINTS ${BFT_DIR}/lib $ENV{BFT_DIR}/lib
    NO_DEFAULT_PATH
    DOC "The BFT library"
    )
  find_library(BFT_LIBRARY
    NAMES bft    
    DOC "The BFT library"
    )
endif()

# Check for bft_config
find_program(BFT_CONFIG_EXE
  NAMES bft-config    
  HINTS ${BFT_DIR}/bin $ENV{BFT_DIR}/bin
  NO_DEFAULT_PATH
  DOC "The BFT config executable"
  )
find_program(BFT_CONFIG_EXE
  NAMES bft-config    
  DOC "The BFT config executable"
  )

set(BFT_LIBRARIES ${BFT_LIBRARY} CACHE STRING "BFT libraries")

if (BFT_INCLUDE_DIRS AND BFT_LIBRARIES AND BFT_CONFIG_EXE)

    execute_process(COMMAND ${BFT_CONFIG_EXE} --version
                    OUTPUT_VARIABLE OUTPUT)

    if (OUTPUT)
      set(BFT_VERSION ${OUTPUT} CACHE TYPE STRING)
      mark_as_advanced(BFT_VERSION)
    endif()

    # Check if version found is >= required version

    if (BFT_FIND_VERSION)
      if (NOT "${BFT_VERSION}" VERSION_LESS "${BFT_FIND_VERSION}")
        set(BFT_VERSION_OK TRUE)
       endif()
    else()
   # No specific version requested
      set(BFT_VERSION_OK TRUE)
     endif()
endif()

mark_as_advanced(BFT_VERSION_OK)

# Standard package handling
find_package_handle_standard_args(BFT
                                  "BFT could not be found. Be sure to set BFT_DIR."
                                  BFT_LIBRARIES
                                  BFT_INCLUDE_DIRS
                                  BFT_CONFIG_EXE
                                  BFT_VERSION
                                  BFT_VERSION_OK)

 mark_as_advanced(BFT_LIBRARIES
                  BFT_INCLUDE_DIRS
                  BFT_CONFIG_EXE
                  BFT_VERSION
                  BFT_VERSION_OK)
unset(BFT_LIBRARY CACHE)