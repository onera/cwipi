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

configure_file(pycwpclt.pyx.in pycwpclt.pyx)

cython_add_module (pycwpclt
                   PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/pycwpclt.pyx"
                   INCLUDE_DIRECTORIES ${include_dirs}
                   DEPENDS ${SOURCE_PXI})

set_target_properties(pycwpclt PROPERTIES PREFIX "") # do not prefix by lib

get_target_property(CWP_INC pycwpclt INCLUDE_DIRECTORIES)
target_include_directories(pycwpclt PUBLIC ${CWP_INC})

target_link_libraries(pycwpclt Python::NumPy MPI::MPI_CXX)
target_link_libraries(pycwpclt ${LINK_LIBRARIES})
target_link_libraries(pycwpclt cwp_shared)

if (CWP_ENABLE_EXTERNAL_PDM)
  target_link_libraries(pycwpclt pdm::pdm_shared)
endif ()

if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang") OR
    (CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (pycwpclt PROPERTIES COMPILE_FLAGS "-w")
endif()

set (SITE_PACKAGES_CWP "${CMAKE_CWP_INSTALL_PYTHON_DIR}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/pycwpclt")

install(DIRECTORY DESTINATION ${SITE_PACKAGES_CWP})
install(FILES  "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
        DESTINATION            "${SITE_PACKAGES_CWP}")
install(PROGRAMS  "${CMAKE_CURRENT_BINARY_DIR}/pycwpclt.so"
        DESTINATION            "${SITE_PACKAGES_CWP}")
