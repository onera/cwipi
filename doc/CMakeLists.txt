# Writing this file using:
#   - Maia's CMake
#   - https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/

find_package(Doxygen 1.8.9 REQUIRED)
find_package(Sphinx  3     REQUIRED)
find_package(Breathe       REQUIRED)

#------------------------------------------------------------------------------
# Build doxygen documentation
#------------------------------------------------------------------------------

if(NOT DOXYGEN_FOUND)
  message(FATAL_ERROR "Doxygen is needed to build the documentation.")
endif()

if(NOT BREATHE_FOUND)
  message(FATAL_ERROR "Breathe is needed to build the documentation.")
endif()

set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
set(SPHINX_IN ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/conf.py.in)
set(SPHINX_OUT ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/source/conf.py)

set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)

configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)
configure_file(${SPHINX_IN} ${SPHINX_OUT} @ONLY)

if (CWP_ENABLE_SHARED)
  set(DOXYGEN_DEP cwp_shared)
else ()
  set(DOXYGEN_DEP cwp_static)
endif()

if (CWP_ENABLE_Fortran)
  if (CWP_ENABLE_SHARED)
    set(DOXYGEN_DEP ${DOXYGEN_DEP} cwpf_shared)
  else ()
    set(DOXYGEN_DEP ${DOXYGEN_DEP} cwpf_static)
  endif()
endif()

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating API documentation with Doxygen"
                   DEPENDS ${DOXYGEN_DEP} # NB : trying to do documentation after module generation
                   VERBATIM)

add_custom_target(doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

#------------------------------------------------------------------------------
# Build sphinx documentation
#------------------------------------------------------------------------------

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/source)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx/html)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

if(DEFINED ENV{PYTHONPATH})
  set (DOC_PYTHONPATH "PYTHONPATH=${PROJECT_BINARY_DIR}/Cython:$ENV{PYTHONPATH}")
else()
  set (DOC_PYTHONPATH "PYTHONPATH=${PROJECT_BINARY_DIR}/Cython")
endif()

# Only regenerate Sphinx when:
# - Doxygen has rerun
# - Our doc files have been updated
# - The Sphinx config has been updated

file(GLOB_RECURSE SOURCES_SPHINX
  ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/source/*.rst
)

add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                   COMMAND ${CMAKE_COMMAND} -E env ${DOC_PYTHONPATH} env ${DOC_PRELOAD} env ${RUN_LSAN}
                   ${SPHINX_EXECUTABLE} -b html
                   -Q
                   -w ${CMAKE_CURRENT_BINARY_DIR}/sphinx_build_warnings.log
                   ${SPHINX_SOURCE} ${SPHINX_BUILD} > ${CMAKE_CURRENT_BINARY_DIR}/sphinx.log 2>&1
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   DEPENDS
                   # Other docs files you want to track should go here (or in some variable)
                   ${DOXYGEN_INDEX_FILE}
                   ${SOURCES_SPHINX}
                   doxygen
                   MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                   COMMENT "Generating documentation with Sphinx"
                   VERBATIM)

add_custom_target(sphinx ALL DEPENDS ${SPHINX_INDEX_FILE} doxygen)

add_custom_target(doc ALL DEPENDS doxygen sphinx)

if(CWP_ENABLE_PYTHON_BINDINGS)
  add_dependencies(sphinx pycwp pycwpt)
  if(CWP_ENABLE_CLIENT_SERVER)
    add_dependencies(sphinx pycwpclt)
  endif()
endif()

install(DIRECTORY ${SPHINX_BUILD}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/cwp)
