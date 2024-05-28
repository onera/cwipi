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

configure_file(Pypdm.pyx.in Pypdm.pyx)
configure_file(pdm_part_coarse.pxi.in pdm_part_coarse.pxi)

cython_add_module (Pypdm
                   PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/Pypdm.pyx"
                   INCLUDE_DIRECTORIES ${include_dirs}
                   DEPENDS ${SOURCE_PXI})

set_target_properties(Pypdm PROPERTIES PREFIX "") # do not prefix by lib

get_target_property(PDM_INC Pypdm INCLUDE_DIRECTORIES)
target_include_directories(Pypdm PUBLIC ${PDM_INC})

target_link_libraries(Pypdm Python::NumPy MPI::MPI_CXX)
target_link_libraries(Pypdm ${LINK_LIBRARIES})

if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang") OR
    (CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (Pypdm PROPERTIES COMPILE_FLAGS "-w")
endif()

set (SITE_PACKAGES_PDM "${CMAKE_CWP_INSTALL_PYTHON_DIR}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/Pypdm")

install(DIRECTORY DESTINATION ${SITE_PACKAGES_PDM})
install(FILES  "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
        DESTINATION            "${SITE_PACKAGES_PDM}")
install(PROGRAMS  "${CMAKE_CURRENT_BINARY_DIR}/Pypdm.so"
        DESTINATION            "${SITE_PACKAGES_PDM}")
