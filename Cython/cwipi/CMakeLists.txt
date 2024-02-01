file(GLOB SOURCE_PXI
     CONFIGURE_DEPENDS
     *.pxi)

file(GLOB SOURCE_PY
     CONFIGURE_DEPENDS
     *.py)

foreach(F IN LISTS SOURCE_PXI SOURCE_PY)
  get_filename_component(F_BASENAME "${F}" NAME )
  configure_file(${F_BASENAME} ${CMAKE_CURRENT_BINARY_DIR}/${F_BASENAME} )
endforeach()

configure_file(cwipi.pyx.in cwipi.pyx)

cython_add_module (cwipi
                   PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/cwipi.pyx"
                   INCLUDE_DIRECTORIES ${include_dirs}
                   DEPENDS ${SOURCE_PXI} cwp_shared)

set_target_properties(cwipi PROPERTIES PREFIX "") # do not prefix by lib

get_target_property(CWP_INC cwipi INCLUDE_DIRECTORIES)
target_include_directories(cwipi PUBLIC ${CWP_INC})

target_link_libraries(cwipi Python::NumPy MPI::MPI_CXX)
target_link_libraries(cwipi ${LINK_LIBRARIES})
target_link_libraries(cwipi cwp_shared)

if (CWP_ENABLE_EXTERNAL_PDM)
  target_link_libraries(cwipi pdm::pdm_shared)
endif ()

if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang") OR
    (CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (cwipi PROPERTIES COMPILE_FLAGS "-w")
endif()

set (SITE_PACKAGES_CWP_OLD "${CMAKE_CWP_INSTALL_PYTHON_DIR}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/cwipi")

install(DIRECTORY DESTINATION ${SITE_PACKAGES_CWP_OLD})
install(FILES  "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
        DESTINATION            "${SITE_PACKAGES_CWP_OLD}")
install(PROGRAMS  "${CMAKE_CURRENT_BINARY_DIR}/cwipi.so"
        DESTINATION            "${SITE_PACKAGES_CWP_OLD}")
