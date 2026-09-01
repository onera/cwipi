# Sources list
# ------------

file(GLOB SOURCES
     CONFIGURE_DEPENDS
     bft/*.[ch]
     fvm/*.[ch]
     new/*.[ch]
     new/*.[ch]xx
     *.[ch]
     *.[ch]xx)

if (CWP_ENABLE_CLIENT_SERVER)
  list(APPEND SOURCES new/client_server/client.cxx)
  list(APPEND SOURCES new/client_server/client.h)
  list(APPEND SOURCES new/client_server/message.h)
  list(APPEND SOURCES new/client_server/server.cxx)
  list(APPEND SOURCES new/client_server/server.h)
  list(APPEND SOURCES new/client_server/struct.hxx)
  list(APPEND SOURCES new/client_server/transfer.c)
  list(APPEND SOURCES new/client_server/transfer.h)
endif()

if (CWP_ENABLE_Fortran)
  file(GLOB SOURCESF
       CONFIGURE_DEPENDS
       # --> old fortran
       fortran/cwipi_cf.cxx
       fortran/cwipi_f.f90
       fortran/cwipi_printfort_cf.cxx
       fortran/cwipi_printfort.f90
       fortran/cwipi_set_output_logical_unit.f90
       fortran/printfort.f90
       fortran/mod_fvmc_ho_basis.f90
       # --> new fortran
       fortran/new/cwp_cf.cxx
       fortran/new/cwp_f.f90
       fortran/new/cwp_printfort_cf.cxx
       fortran/new/cwp_printfort.f90)
  list(APPEND SOURCES fortran/callfortinterpfct.f90)
  list(APPEND SOURCES fortran/callforthointerpfct.f90)
  list(APPEND SOURCES fortran/callfortlocinterpfct.f90)
endif()

if (NOT CWP_ENABLE_EXTERNAL_PDM)
  file(GLOB_RECURSE SOURCESPDM
      ${PDM_SOURCE_DIR}/src/*.[ch])
  add_library(pdm OBJECT ${SOURCESPDM})
  list(APPEND PDM_INC ${PDM_SOURCE_DIR}/src
                      ${PDM_SOURCE_DIR}/src/algorithm
                      ${PDM_SOURCE_DIR}/src/exchange_protocol
                      ${PDM_SOURCE_DIR}/src/ext_wrapper
                      ${PDM_SOURCE_DIR}/src/geometry
                      ${PDM_SOURCE_DIR}/src/gnum_utils
                      ${PDM_SOURCE_DIR}/src/gpu
                      ${PDM_SOURCE_DIR}/src/io
                      ${PDM_SOURCE_DIR}/src/mesh
                      ${PDM_SOURCE_DIR}/src/mesh_struct
                      ${PDM_SOURCE_DIR}/src/mesh_topology
                      ${PDM_SOURCE_DIR}/src/meshgen
                      ${PDM_SOURCE_DIR}/src/mpi_wrapper
                      ${PDM_SOURCE_DIR}/src/ppart
                      ${PDM_SOURCE_DIR}/src/struct
                      ${PDM_SOURCE_DIR}/src/test_utils
                      ${PDM_SOURCE_DIR}/src/util)
  target_include_directories(pdm
                             PRIVATE   $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
                             PUBLIC    $<INSTALL_INTERFACE:include/>
                             PRIVATE  "$<BUILD_INTERFACE:${PDM_INC}>")
else()
  set(PDM_INC "")
endif()

add_library(cwp OBJECT ${SOURCES})

target_include_directories(cwp
                           INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           PRIVATE   $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           PRIVATE   $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
                           PUBLIC    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                           PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/new>
                           PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fvm>
                           PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bft>
                           PUBLIC    $<INSTALL_INTERFACE:include/>
                           PRIVATE  "$<BUILD_INTERFACE:${PDM_INC}>")

if (CWP_ENABLE_Fortran)

  add_library(cwpf OBJECT ${SOURCESF})

  target_include_directories(cwpf
                             INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                             PRIVATE   $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                             PRIVATE   $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
                             PUBLIC    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                             PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                             PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/new>
                             PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fvm>
                             PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bft>
                             PUBLIC    $<INSTALL_INTERFACE:include/>
                             PRIVATE  "$<BUILD_INTERFACE:${PDM_INC}>")

  add_dependencies(cwpf cwp)

endif()

if (BLAS_FOUND)
  string(REPLACE " " ";" BLAS_COMPILER_FLAGS "${BLAS_COMPILER_FLAGS}")
  target_compile_options (cwp PRIVATE ${BLAS_COMPILER_FLAGS})
  if (BLAS_INCLUDE_DIRS)
    target_include_directories(cwp PRIVATE ${BLAS_INCLUDE_DIRS})
  endif()
endif()

# Liste des header
# ----------------

set (LIST_HEADER "")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/cwipi.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/new/cwp.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/new/cwp_logging.h")
list(APPEND LIST_HEADER "${CMAKE_BINARY_DIR}/cwipi_config.h")

if (CWP_ENABLE_Fortran)
  set (LIST_HEADERF "")
  list(APPEND LIST_HEADERF "${CMAKE_BINARY_DIR}/cwipi_configf.h")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/cwipi.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/cwp.mod")
endif ()

# Add Fortran dependencies if Fortran is called in C
# --------------------------------------------------

unset(LINK_LIBRARIES)

if ((NOT CWP_HAVE_NOT_FORTRAN_IN_C) AND FORTRAN_LIBRARIES)
  list(APPEND LINK_LIBRARIES ${FORTRAN_LIBRARIES})
endif()

if ((NOT CWP_HAVE_NOT_FORTRAN_IN_C) AND FORTRAN_LIBRARIES_FLAG)
  set_target_properties(cwp PROPERTIES COMPILE_FLAGS ${FORTRAN_LIBRARIES_FLAG})
endif ()

# Add BLAS libraries
# --------------------------------------------------

if (BLAS_FOUND)
  list(APPEND LINK_LIBRARIES ${BLAS_LIBRARIES})
endif()

# Add cwp libraries
# --------------------------------------------------

if (CWP_ENABLE_SHARED)
  set_property(TARGET ${cwp} PROPERTY POSITION_INDEPENDENT_CODE 1)
  if (CWP_ENABLE_EXTERNAL_PDM)
    add_library(cwp_shared SHARED $<TARGET_OBJECTS:cwp>)
    target_link_libraries (cwp        pdm::pdm_shared)
    target_link_libraries (cwp_shared pdm::pdm_shared)
  else ()
    add_library(cwp_shared SHARED $<TARGET_OBJECTS:cwp> $<TARGET_OBJECTS:pdm>)
    add_dependencies(cwp_shared pdm)
    if (CWP_ENABLE_HIDE_PDM_SYMBOLS)
      if (APPLE)
        set(EXPORTED_SYMBOLS_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../cwp_apple.map)
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -exported_symbols_list ${EXPORTED_SYMBOLS_FILE}")
        set_target_properties(cwp_shared PROPERTIES LINK_DEPENDS ${EXPORTED_SYMBOLS_FILE})
      else()
        set_target_properties(cwp_shared PROPERTIES LINK_FLAGS "-Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/../cwp.map")
      endif()
    endif()
  endif ()
  add_dependencies(cwp_shared cwp)
  set_target_properties(cwp_shared PROPERTIES OUTPUT_NAME "cwp")
  set_target_properties(cwp_shared PROPERTIES VERSION ${CWP_VERSION})
  target_link_libraries(cwp_shared ${LINK_LIBRARIES})
  set_target_properties(cwp_shared PROPERTIES
                        PUBLIC_HEADER "${LIST_HEADER}")
  get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)
  target_include_directories(cwp_shared PUBLIC "${CWP_INC}")

  if (CWP_ENABLE_Fortran)
    set_property(TARGET ${cwpf} PROPERTY POSITION_INDEPENDENT_CODE 1)
    if (CWP_ENABLE_EXTERNAL_PDM)
      add_library(cwpf_shared SHARED $<TARGET_OBJECTS:cwpf>)
      target_link_libraries (cwpf        pdm::pdmf_shared)
      target_link_libraries (cwpf_shared pdm::pdmf_shared)
    else ()
      if (CWP_ENABLE_HIDE_PDM_SYMBOLS)
        add_library(cwpf_shared SHARED $<TARGET_OBJECTS:cwpf> $<TARGET_OBJECTS:pdm>)
        add_dependencies(cwpf_shared pdm)
        if (APPLE)
          set(EXPORTED_SYMBOLS_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../cwp_apple.map)
          set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -exported_symbols_list ${EXPORTED_SYMBOLS_FILE}")
          set_target_properties(cwpf_shared PROPERTIES LINK_DEPENDS ${EXPORTED_SYMBOLS_FILE})
        else()
          set_target_properties(cwpf_shared PROPERTIES LINK_FLAGS "-Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/../cwp.map")
        endif()
      else ()
        add_library(cwpf_shared SHARED $<TARGET_OBJECTS:cwpf>)
      endif()
    endif ()
    add_dependencies(cwpf_shared cwpf cwp_shared)
    set_target_properties(cwpf_shared PROPERTIES OUTPUT_NAME "cwpf")
    set_target_properties(cwpf_shared PROPERTIES VERSION ${CWP_VERSION})
    target_link_libraries(cwpf_shared cwp_shared ${CXX_LIBRARIES})
    set_target_properties(cwpf_shared PROPERTIES
                          PUBLIC_HEADER "${LIST_HEADERF}")
    get_target_property(CWP_INC cwpf INCLUDE_DIRECTORIES)
    target_include_directories(cwpf_shared PUBLIC "${CWP_INC}")
  endif()

endif()

if (CWP_ENABLE_STATIC)
  set_property(TARGET ${cwp} PROPERTY POSITION_INDEPENDENT_CODE 1)
  if (CWP_ENABLE_EXTERNAL_PDM)
    add_library(cwp_static STATIC $<TARGET_OBJECTS:cwp>)
    target_link_libraries (cwp        pdm::pdm_static)
    target_link_libraries (cwp_static pdm::pdm_static)
  else ()
    add_library(cwp_static STATIC $<TARGET_OBJECTS:cwp> $<TARGET_OBJECTS:pdm>)
    add_dependencies(cwp_static pdm)
  endif ()
  add_dependencies(cwp_static cwp)
  set_target_properties(cwp_static PROPERTIES OUTPUT_NAME "cwp")
  set_target_properties(cwp_static PROPERTIES VERSION ${CWP_VERSION})
  target_link_libraries(cwp_static ${LINK_LIBRARIES})
  set_target_properties(cwp_static PROPERTIES
                        PUBLIC_HEADER "${LIST_HEADER}")
  get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)
  target_include_directories(cwp_static PUBLIC "${CWP_INC}")

  if (CWP_ENABLE_Fortran)
    set_property(TARGET ${cwpf} PROPERTY POSITION_INDEPENDENT_CODE 1)
    if (CWP_ENABLE_EXTERNAL_PDM)
      add_library(cwpf_static STATIC $<TARGET_OBJECTS:cwpf>)
      target_link_libraries (cwpf        pdm::pdmf_static)
      target_link_libraries (cwpf_static pdm::pdmf_static)
    else ()
      add_library(cwpf_static STATIC $<TARGET_OBJECTS:cwpf>)
    endif ()
    add_dependencies(cwpf_static cwpf cwp_static)
    set_target_properties(cwpf_static PROPERTIES OUTPUT_NAME "cwpf")
    set_target_properties(cwpf_static PROPERTIES VERSION ${CWP_VERSION})
    target_link_libraries(cwpf_static cwp_static ${MPI_CXX_LIBRARIES})
    set_target_properties(cwpf_static PROPERTIES
                          PUBLIC_HEADER "${LIST_HEADERF}")
    get_target_property(CWP_INC cwpf INCLUDE_DIRECTORIES)
    target_include_directories(cwpf_static PUBLIC "${CWP_INC}")
  endif()

endif()

# Spécification de l'installation
# -------------------------------

set (_targets "")

list (APPEND _targets cwp)
if (CWP_ENABLE_Fortran)
  list (APPEND _targets cwpf)
endif()

if (CWP_ENABLE_SHARED)
  list (APPEND _targets cwp_shared)
  if (CWP_ENABLE_Fortran)
    list (APPEND _targets cwpf_shared)
  endif()
endif()

if (CWP_ENABLE_STATIC)
  list (APPEND _targets cwp_static)
  if (CWP_ENABLE_Fortran)
    list (APPEND _targets cwpf_static)
  endif()
endif()

foreach(_target IN ITEMS ${_targets})
  install(TARGETS ${_target} EXPORT ${_target}Targets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include
    ARCHIVE DESTINATION lib
  )
  install(EXPORT ${_target}Targets
    FILE ${_target}Targets.cmake
    NAMESPACE cwp::
    DESTINATION lib/cmake/${_target}
  )
  add_library(cwp::${_target} ALIAS ${_target})
endforeach()
