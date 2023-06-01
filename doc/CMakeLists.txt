# Writing this file using:
#   - Maia's CMake
#   - https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/

if(CWP_ENABLE_DOCUMENTATION)
    find_package(Doxygen 1.8.9 REQUIRED) 
    find_package(Sphinx 3 REQUIRED)

    #------------------------------------------------------------------------------
    # Build doxygen documentation
    #------------------------------------------------------------------------------

    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)

    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                       COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                       MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                       COMMENT "Generating API documentation with Doxygen"
#                       DEPENDS pycwp # NB : trying to do documentation after module generation
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

    if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
      execute_process(COMMAND gcc -print-file-name=libasan.so OUTPUT_VARIABLE PRELOAD_ASAN OUTPUT_STRIP_TRAILING_WHITESPACE)
      set(DOC_PRELOAD "LD_PRELOAD=${PRELOAD_ASAN}")
    endif()

    # Only regenerate Sphinx when:
    # - Doxygen has rerun
    # - Our doc files have been updated
    # - The Sphinx config has been updated
    add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                       COMMAND ${CMAKE_COMMAND} -E env ${DOC_PYTHONPATH} env ${DOC_PRELOAD}
                         ${SPHINX_EXECUTABLE} -Q -b html # -Q to suppress warning output
                         #-w sphinx_build_warnings.log
                         # Tell Breathe where to find the Doxygen output
                         -Dbreathe_projects.cwipi=${CMAKE_CURRENT_BINARY_DIR}/xml
                         ${SPHINX_SOURCE} ${SPHINX_BUILD}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       DEPENDS
                       # Other docs files you want to track should go here (or in some variable)
                       ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/source/index.rst
                       ${DOXYGEN_INDEX_FILE}
                       MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                       COMMENT "Generating documentation with Sphinx")

    add_custom_target(doc ALL DEPENDS ${SPHINX_INDEX_FILE})

    add_dependencies(doc doxygen pycwp)

    add_custom_command(TARGET doc POST_BUILD
                       COMMAND ;
                       COMMENT "Sphinx documentation can be found at ${CMAKE_BINARY_DIR}/doc/sphinx/html. Open index.html")

    install(DIRECTORY ${SPHINX_BUILD}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/cwp)
endif()
