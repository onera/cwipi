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

file ( COPY . DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" FILES_MATCHING PATTERN "*.pxi")
file ( COPY . DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" FILES_MATCHING PATTERN "*.py")

# old API
configure_file(cwipi.pyx.in cwipi.pyx)
# new API
configure_file(pycwp.pyx.in pycwp.pyx)
# client-server
if (CWP_ENABLE_CLIENT_SERVER)
  configure_file(client_server/pycwpclt.pyx.in client_server/pycwpclt.pyx)
endif()

set (include_dirs "")
list (APPEND include_dirs ${Mpi4Py_INCLUDE_DIR})
list (APPEND include_dirs ${NumPy_INCLUDE_DIR})
list (APPEND include_dirs ${CMAKE_CURRENT_BINARY_DIR})
list (APPEND include_dirs ${CMAKE_CURRENT_SOURCE_DIR})
list (APPEND include_dirs ${CMAKE_BINARY_DIR})


include(UseCython)

# old API
cython_add_module (cwipi
                   PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/cwipi.pyx"
                   INCLUDE_DIRECTORIES ${include_dirs})
# new API
cython_add_module (pycwp
                   PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/pycwp.pyx"
                   INCLUDE_DIRECTORIES ${include_dirs})
# client-server
if (CWP_ENABLE_CLIENT_SERVER)
  cython_add_module (pycwpclt
                     PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/client_server/pycwpclt.pyx"
                     INCLUDE_DIRECTORIES ${include_dirs})
endif()

# old API
get_target_property(CWP_INC cwipi INCLUDE_DIRECTORIES)
target_include_directories(cwipi PUBLIC ${CWP_INC})

# new API
get_target_property(CWP_INC pycwp INCLUDE_DIRECTORIES)
target_include_directories(pycwp PUBLIC ${CWP_INC})
# client-server
if (CWP_ENABLE_CLIENT_SERVER)
  get_target_property(CWP_INC pycwpclt INCLUDE_DIRECTORIES)
  target_include_directories(pycwpclt PUBLIC ${CWP_INC})
endif()

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
  if (CWP_ENABLE_CLIENT_SERVER)
    target_link_libraries(pycwpclt ${LINK_LIBRARIES} pdm_shared pdm_mpi_shared)
  endif()
  target_link_libraries(pycwp ${LINK_LIBRARIES} pdm_shared pdm_mpi_shared)
  target_link_libraries(cwipi ${LINK_LIBRARIES} pdm_shared pdm_mpi_shared)
endif()

# old API
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

# new API
target_link_libraries(pycwp cwp_shared)
if (CWP_ENABLE_Fortran)
  target_link_libraries(pycwp cwpf_shared)
endif()

if (CWP_ENABLE_SPACE_NODAL_BASIS)
  target_link_libraries(pycwp spaceBasis_shared)
endif()

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
   set_target_properties(pycwp PROPERTIES
                         COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

# client-server
if (CWP_ENABLE_CLIENT_SERVER)
  target_link_libraries(pycwpclt cwp_shared)
  if (CWP_ENABLE_Fortran)
    target_link_libraries(pycwpclt cwpf_shared)
  endif()

  if (CWP_ENABLE_SPACE_NODAL_BASIS)
    target_link_libraries(pycwpclt spaceBasis_shared)
  endif()

  if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
     set_target_properties(pycwpclt PROPERTIES
                           COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
  endif()
endif()

set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/timestamp")

# old API
set(DEPS_OLD    "${CMAKE_CURRENT_SOURCE_DIR}/cwipi/__init__.py")
add_custom_command(OUTPUT ${OUTPUT}
                   COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet build
                   COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                   DEPENDS ${DEPS_OLD})

# new API
set(DEPS_NEW    "${CMAKE_CURRENT_SOURCE_DIR}/pycwp/__init__.py")
add_custom_command(OUTPUT ${OUTPUT}
                   COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet build
                   COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                   DEPENDS ${DEPS_NEW})

# client-server
if (CWP_ENABLE_CLIENT_SERVER)
  set(DEPS_NEW    "${CMAKE_CURRENT_SOURCE_DIR}/client_server/pycwpclt/__init__.py")
  add_custom_command(OUTPUT ${OUTPUT}
                     COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet build
                     COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                     DEPENDS ${DEPS_NEW})
endif()

# old API
set (SITE_PACKAGES_CWP_OLD "${CMAKE_CWP_INSTALL_PYTHON_DIR}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/cwipi")

if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang") OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (cwipi PROPERTIES COMPILE_FLAGS "-w")
endif()

# new API
set (SITE_PACKAGES_CWP_NEW "${CMAKE_CWP_INSTALL_PYTHON_DIR}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/pycwp")

if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang") OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (pycwp PROPERTIES COMPILE_FLAGS "-w")
endif()

# client-server
if (CWP_ENABLE_CLIENT_SERVER)
  set (SITE_PACKAGES_CWP_CS "${CMAKE_CWP_INSTALL_PYTHON_DIR}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/pycwpclt")

  if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
      (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
      (CMAKE_C_COMPILER_ID STREQUAL "AppleClang") OR
      (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
    set_target_properties (pycwpclt PROPERTIES COMPILE_FLAGS "-w")
  endif()
endif()

# old API
install(DIRECTORY DESTINATION ${SITE_PACKAGES_CWP_OLD})
install(FILES  "${CMAKE_CURRENT_SOURCE_DIR}/cwipi/__init__.py"
        DESTINATION            "${SITE_PACKAGES_CWP_OLD}")
install(PROGRAMS  "${CMAKE_CURRENT_BINARY_DIR}/cwipi.so"
        DESTINATION            "${SITE_PACKAGES_CWP_OLD}")

# new API
install(DIRECTORY DESTINATION ${SITE_PACKAGES_CWP_NEW})
install(FILES  "${CMAKE_CURRENT_SOURCE_DIR}/pycwp/__init__.py"
        DESTINATION            "${SITE_PACKAGES_CWP_NEW}")
install(PROGRAMS  "${CMAKE_CURRENT_BINARY_DIR}/pycwp.so"
        DESTINATION            "${SITE_PACKAGES_CWP_NEW}")

# client-server
if (CWP_ENABLE_CLIENT_SERVER)
  install(DIRECTORY DESTINATION ${SITE_PACKAGES_CWP_CS})
  install(FILES  "${CMAKE_CURRENT_SOURCE_DIR}/client_server/pycwpclt/__init__.py"
          DESTINATION            "${SITE_PACKAGES_CWP_CS}")
  install(PROGRAMS  "${CMAKE_CURRENT_BINARY_DIR}/pycwpclt.so"
          DESTINATION            "${SITE_PACKAGES_CWP_CS}")
endif()
