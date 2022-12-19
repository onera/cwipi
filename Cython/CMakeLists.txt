
if (MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif() 

if (MPI_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif(MPI_Fortran_COMPILER) 

if ((Mpi4Py_VERSION VERSION_GREATER 2.0.0) OR 
    (Mpi4Py_VERSION VERSION_EQUAL 2.0.0))
   set (mpi4pylibmpi libmpi)
else ()
   set (mpi4pylibmpi mpi_c)
endif()

# configure_file(cwipi.pyx.in cwipi.pyx)
configure_file(cwp.pyx.in cwp.pyx)

file (COPY . DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" FILES_MATCHING PATTERN "*.pxi")
file (COPY . DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" FILES_MATCHING PATTERN "*.py")

set (include_dirs "")
list (APPEND include_dirs ${Mpi4Py_INCLUDE_DIR})
list (APPEND include_dirs ${NumPy_INCLUDE_DIR})
list (APPEND include_dirs ${CMAKE_CURRENT_BINARY_DIR})
list (APPEND include_dirs ${CMAKE_CURRENT_SOURCE_DIR})
list (APPEND include_dirs ${CMAKE_BINARY_DIR})

include(UseCython)
cython_add_module (cwipi 
                   # PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/cwipi.pyx"
                   PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/cwp.pyx"
                   OTHER_SOURCES $<TARGET_OBJECTS:cwp>
                   INCLUDE_DIRECTORIES ${include_dirs})

get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)                          
target_include_directories(cwipi PUBLIC ${CWP_INC})

unset(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
if (${PARMETIS_FOUND})
  list(APPEND LINK_LIBRARIES ${PARMETIS_LIBRARIES})
  list(APPEND LINK_LIBRARIES ${METIS_LIBRARIES})
endif ()
if (${PTSCOTCH_FOUND})
  list(APPEND LINK_LIBRARIES ${PTSCOTCH_LIBRARIES})
  list(APPEND LINK_LIBRARIES ${SCOTCH_LIBRARIES})
endif ()

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/new") 
 target_link_libraries(cwipi ${LINK_LIBRARIES} pdm_shared pdm_mpi_shared)
endif()

target_link_libraries(cwipi cwp_shared)
if (CWP_ENABLE_Fortran)
  target_link_libraries(cwipi cwpf_shared)
endif()

if (CWP_ENABLE_SPACE_NODAL_BASIS)
  target_link_libraries(cwipi spaceBasis_shared)
endif()


if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
   set_target_properties(cwipi PROPERTIES
                              COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/cwipi/__init__.py")
set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/timestamp")

add_custom_command(OUTPUT ${OUTPUT}
                       COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet build
                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                       DEPENDS ${DEPS})

set (SITE_PACKAGES_CWP "${CMAKE_CWP_INSTALL_PYTHON_DIR}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/cwipi")

add_custom_target(target ALL DEPENDS ${OUTPUT})


if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (cwipi PROPERTIES COMPILE_FLAGS "-w")
endif()

install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet install --prefix=${CMAKE_CWP_INSTALL_PYTHON_DIR} --install-data=${SITE_PACKAGES_CWP} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
