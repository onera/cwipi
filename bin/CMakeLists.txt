# ------------
# Build server
# ------------

unset(LINK_LIBRARIES)

if (CWP_ENABLE_SHARED)
  list(APPEND LINK_LIBRARIES cwp_shared)
else ()
  list(APPEND LINK_LIBRARIES cwp_static)
endif ()

set(SOURCESPDM "")
if (CWP_ENABLE_EXTERNAL_PDM)
  if (CWP_ENABLE_SHARED)
    list(APPEND LINK_LIBRARIES
          pdm::pdm_mpi_shared
          pdm::pdm_shared)
  else ()
    list(APPEND LINK_LIBRARIES
          pdm::pdm_mpi_static
          pdm::pdm_static)
  endif ()
else ()
  file(GLOB_RECURSE SOURCESPDM
    ${PDM_SOURCE_DIR}/src/*/pdm_io.c
    ${PDM_SOURCE_DIR}/src/*/pdm_timer.c
    ${PDM_SOURCE_DIR}/src/*/pdm_mpi.c
    ${PDM_SOURCE_DIR}/src/*/pdm_array.c
    ${PDM_SOURCE_DIR}/src/*/pdm_file_seq.c
    ${PDM_SOURCE_DIR}/src/*/pdm_file_par.c
    ${PDM_SOURCE_DIR}/src/*/pdm_mpi_node_first_rank.c
    ${PDM_SOURCE_DIR}/src/*/pdm_printf.c
    ${PDM_SOURCE_DIR}/src/*/pdm_error.c)
endif ()

add_executable (cwp_server server_main.c ${SOURCESPDM})
target_include_directories(cwp_server PRIVATE ${CMAKE_SOURCE_DIR}
                           PRIVATE ${CMAKE_BINARY_DIR}
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(cwp_server ${LINK_LIBRARIES})
install(TARGETS cwp_server DESTINATION bin)
