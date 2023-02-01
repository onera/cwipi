# ------------
# Build server
# ------------

if (MPI_C_COMPILER)
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

if (MPI_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

unset(LINK_LIBRARIES)

if (NOT MPI_C_COMPILER)
  if (MPI_C_LIBRARIES)
    list(APPEND LINK_LIBRARIES ${MPI_C_LIBRARIES})
  endif()
endif ()

if (NOT MPI_CXX_COMPILER)
  if (MPI_CXX_LIBRARIES)
    list(APPEND LINK_LIBRARIES ${MPI_C_LIBRARIES})
  else ()
    if (MPI_C_LIBRARIES)
     list(APPEND LINK_LIBRARIES ${MPI_C_LIBRARIES})
    endif()
  endif()
endif ()

if (CWP_ENABLE_STATIC)
  if (CWP_ENABLE_NEW_API)
     if (EXISTS ${CWP_SOURCE_DIR}/paradigm)
       list(APPEND LINK_LIBRARIES
            pdm_mpi_static
            pdm_static)
     endif()
   endif()
   list(APPEND LINK_LIBRARIES cwp_static)
   list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
else ()
  if (CWP_ENABLE_NEW_API)
     if (EXISTS ${CWP_SOURCE_DIR}/paradigm)
       list(APPEND LINK_LIBRARIES
            pdm_mpi_shared
            pdm_shared)
     endif()
   endif()
   list(APPEND LINK_LIBRARIES cwp_shared)
   list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
endif()

# add executable

if (CWP_ENABLE_NEW_API)
  if (CWP_ENABLE_CLIENT_SERVER)
    add_executable (server_main server_main.c)
    target_include_directories(server_main PRIVATE ${CMAKE_SOURCE_DIR}
                               PRIVATE ${CMAKE_BINARY_DIR}
                               PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    target_link_libraries(server_main ${LINK_LIBRARIES})
    install(TARGETS server_main DESTINATION bin)
  endif()
endif()
