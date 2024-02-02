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

set(SOURCESPDM  "")
set(SOURCESPDMF "")
if (NOT CWP_ENABLE_EXTERNAL_PDM)
  file(GLOB_RECURSE SOURCESPDM
      ${PDM_SOURCE_DIR}/src/*.[ch])
  list(FILTER SOURCESPDM EXCLUDE REGEX "pdm_no_mpi*")
endif()

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
  if (NOT CWP_ENABLE_EXTERNAL_PDM)
    list(APPEND SOURCESPDMF ${PDM_SOURCE_DIR}/src/util/pdm_fortran.F90)
    list(APPEND SOURCESPDMF ${PDM_SOURCE_DIR}/src/pdm.F90)
    list(APPEND SOURCESPDMF ${PDM_SOURCE_DIR}/src/struct/pdm_pointer_array.F90)
    if (CWP_ENABLE_TESTS)
      list(APPEND SOURCESPDMF ${PDM_SOURCE_DIR}/src/meshgen/pdm_generate_mesh.F90)
      list(APPEND SOURCESPDMF ${PDM_SOURCE_DIR}/src/meshgen/pdm_sphere_surf_gen.F90)
    endif()
  endif()
endif()

add_library(cwp OBJECT ${SOURCES} ${SOURCESPDM} ${SOURCESF} ${SOURCESPDMF})

target_include_directories(cwp
                           INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           PRIVATE   $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           PRIVATE   $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
                           PUBLIC    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                           PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/new>
                           PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fvm>
                           PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bft>
                           PUBLIC    $<INSTALL_INTERFACE:include/>)

if (NOT CWP_ENABLE_EXTERNAL_PDM)
  target_include_directories(cwp
                             PRIVATE $<BUILD_INTERFACE:${PDM_SOURCE_DIR}/src>
                             PRIVATE $<BUILD_INTERFACE:${PDM_SOURCE_DIR}/src/meshgen>
                             PRIVATE $<BUILD_INTERFACE:${PDM_SOURCE_DIR}/src/mesh>
                             PRIVATE $<BUILD_INTERFACE:${PDM_SOURCE_DIR}/src/struct>
                             PRIVATE $<BUILD_INTERFACE:${PDM_SOURCE_DIR}/src/io>
                             PRIVATE $<BUILD_INTERFACE:${PDM_SOURCE_DIR}/src/pario>
                             PRIVATE $<BUILD_INTERFACE:${PDM_SOURCE_DIR}/src/ppart>
                             PRIVATE $<BUILD_INTERFACE:${PDM_SOURCE_DIR}/src/util>
                             PRIVATE $<BUILD_INTERFACE:${PDM_SOURCE_DIR}/src/ext_wrapper>
                             PRIVATE $<BUILD_INTERFACE:${PDM_SOURCE_DIR}/src/mpi_wrapper>)
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
list(APPEND LIST_HEADER "${CMAKE_BINARY_DIR}/cwipi_config.h")

if (CWP_ENABLE_Fortran)
  list(APPEND LIST_HEADER "${CMAKE_BINARY_DIR}/cwipi_configf.h")
  list(APPEND LIST_HEADER "${CMAKE_CURRENT_BINARY_DIR}/cwipi.mod")
  list(APPEND LIST_HEADER "${CMAKE_CURRENT_BINARY_DIR}/cwp.mod")
endif ()

# Options du link
# ---------------

unset(LINK_LIBRARIES)

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  list(APPEND LINK_LIBRARIES "m")
endif()

# Add Fortran dependencies if Fortran is called in C
# --------------------------------------------------

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
  add_library(cwp_shared SHARED $<TARGET_OBJECTS:cwp>)
  set_target_properties(cwp_shared PROPERTIES LINK_FLAGS "-Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/../cwp.map")
  add_dependencies(cwp_shared cwp)
  set_target_properties(cwp_shared PROPERTIES OUTPUT_NAME "cwp")
  set_target_properties(cwp_shared PROPERTIES VERSION ${CWP_VERSION})
  target_link_libraries(cwp_shared ${LINK_LIBRARIES})
  set_target_properties(cwp_shared PROPERTIES
                        PUBLIC_HEADER "${LIST_HEADER}")
  if (CWP_ENABLE_EXTERNAL_PDM)
    target_link_libraries (cwp        pdm::pdm_shared)
    target_link_libraries (cwp_shared pdm::pdm_shared)
  endif ()
  get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)
  target_include_directories(cwp_shared PUBLIC ${CWP_INC})
  if (CWP_ENABLE_TESTS AND NOT CWP_ENABLE_EXTERNAL_PDM) #add library with non-private paradigm for tests
    add_library(cwp_test_shared SHARED $<TARGET_OBJECTS:cwp>)
    add_dependencies(cwp_test_shared cwp)
    target_link_libraries(cwp_test_shared ${LINK_LIBRARIES})
    if (CWP_ENABLE_EXTERNAL_PDM)
      target_link_libraries (cwp_test_shared pdm::pdm_shared)
    endif ()
    get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)
    target_include_directories(cwp_test_shared PUBLIC ${CWP_INC})
  endif ()
endif()

if (CWP_ENABLE_STATIC)
  set_property(TARGET ${cwp} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(cwp_static STATIC $<TARGET_OBJECTS:cwp>)
  add_dependencies(cwp_static cwp)
  set_target_properties(cwp_static PROPERTIES OUTPUT_NAME "cwp")
  set_target_properties(cwp_static PROPERTIES VERSION ${CWP_VERSION})
  target_link_libraries(cwp_static ${LINK_LIBRARIES})
  set_target_properties(cwp_static PROPERTIES
                        PUBLIC_HEADER "${LIST_HEADER}")
  if (CWP_ENABLE_EXTERNAL_PDM)
    target_link_libraries (cwp        pdm::pdm_static)
    target_link_libraries (cwp_static pdm::pdm_static)
  endif ()
  get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)
  target_include_directories(cwp_static PUBLIC ${CWP_INC})
endif()

# Spécification de l'installation
# -------------------------------

set (_targets "")

if (CWP_ENABLE_SHARED)
  list (APPEND _targets cwp_shared)
endif()

if (CWP_ENABLE_STATIC)
  list (APPEND _targets cwp_static)
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
