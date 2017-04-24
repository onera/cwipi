# - Try to find PALM    
# Once done this will define
#
#  PALM_FOUND        - system has found PALM
#  PALM_INCLUDE_DIR  - include directories for PALM
#  PALM_LIBRARIES     - libraries for PALM
#

set(PALM_FOUND FALSE)

# Check for header file
find_path(PALM_INCLUDE_DIR config.h palmlib.h
  HINTS ${PALM_DIR}/include $ENV{PALM_DIR}/include
  DOC "Directory where the PALM header is located"
)

# Check for palm
find_library(PALM_LIBRARY
  NAMES palm 
  HINTS ${PALM_DIR}/lib $ENV{PALM_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The PALM library"
  )
find_library(PALM_LIBRARY
  NAMES palm 
  DOC "The PALM library"
  )

# Check for server palm
find_library(PALM_SERVER_LIBRARY
  NAMES  serv_opalm
  HINTS ${PALM_SERVER_LIB_DIR} $ENV{PALM_SERVER_LIB_DIR}
  NO_DEFAULT_PATH
  DOC "The PALM server library"
  )
find_library(PALM_SERVER_LIBRARY
  NAMES serv_opalm
  DOC "The PALM server library"
  )

set(PALM_LIBRARIES ${PALM_SERVER_LIBRARY} ${PALM_LIBRARY} CACHE STRING "Palm libraries")

mark_as_advanced(PALM_LIBRARIES PALM_SERVER_LIBRARY PALM_LIBRARY PALM_INCLUDE_DIR)

# Standard package handling
find_package_handle_standard_args(PALM
                                  "PALM could not be found. Be sure to set PALM_DIR PALM_SERVER_DIR."
                                  PALM_LIBRARIES
                                  PALM_LIBRARY
                                  PALM_SERVER_LIBRARY
                                  PALM_INCLUDE_DIR)


