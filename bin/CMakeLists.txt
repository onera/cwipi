# ------------
# Build server
# ------------

unset(LINK_LIBRARIES)

if (NOT CWP_ENABLE_SHARED)
  list(APPEND LINK_LIBRARIES cwp_static)
endif ()

if (CWP_ENABLE_EXTERNAL_PDM)
  if (NOT CWP_ENABLE_SHARED)
    list(APPEND LINK_LIBRARIES
          pdm::pdm_mpi_static
          pdm::pdm_static)
  else ()
    list(APPEND LINK_LIBRARIES cwp_shared)
    list(APPEND LINK_LIBRARIES
          pdm::pdm_mpi_shared
          pdm::pdm_shared)
  endif ()
else ()
  if (CWP_ENABLE_SHARED)
    list(APPEND LINK_LIBRARIES cwp_test_shared)
  endif ()
endif ()

add_executable (cwp_server server_main.c)
target_include_directories(cwp_server PRIVATE ${CMAKE_SOURCE_DIR}
                           PRIVATE ${CMAKE_BINARY_DIR}
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(cwp_server ${LINK_LIBRARIES})
install(TARGETS cwp_server DESTINATION bin)
