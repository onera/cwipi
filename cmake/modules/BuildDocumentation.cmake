# Adapted from Maia
macro(build_sphinx_documentation)
# 1. Sphinx
  find_package(Sphinx 3 REQUIRED)
  set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/source)
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx/html)
  set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

  file(GLOB_RECURSE doc_files ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/source/*)
  add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                     COMMAND ${SPHINX_EXECUTABLE} -b html ${SPHINX_SOURCE} ${SPHINX_BUILD}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     DEPENDS ${doc_files}
                     COMMENT "Generating Sphinx documentation")

  add_custom_target(cwp_sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})

# 2. Install
  install(DIRECTORY ${SPHINX_BUILD}
          DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/cwp)
endmacro()
