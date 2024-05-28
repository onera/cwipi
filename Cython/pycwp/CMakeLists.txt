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

configure_file(pycwp.pyx.in pycwp.pyx)

cython_add_module (pycwp
                   PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/pycwp.pyx"
                   INCLUDE_DIRECTORIES ${include_dirs}
                   DEPENDS ${SOURCE_PXI})

set_target_properties(pycwp PROPERTIES PREFIX "") # do not prefix by lib

get_target_property(CWP_INC pycwp INCLUDE_DIRECTORIES)
target_include_directories(pycwp PUBLIC ${CWP_INC})

target_link_libraries(pycwp Python::NumPy MPI::MPI_CXX)
target_link_libraries(pycwp ${LINK_LIBRARIES})

if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang") OR
    (CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (pycwp PROPERTIES COMPILE_FLAGS "-w")
endif()

set (SITE_PACKAGES_CWP "${CMAKE_CWP_INSTALL_PYTHON_DIR}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/pycwp")

install(DIRECTORY DESTINATION ${SITE_PACKAGES_CWP})
install(FILES  "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
        DESTINATION            "${SITE_PACKAGES_CWP}")
install(PROGRAMS  "${CMAKE_CURRENT_BINARY_DIR}/pycwp.so"
        DESTINATION            "${SITE_PACKAGES_CWP}")
