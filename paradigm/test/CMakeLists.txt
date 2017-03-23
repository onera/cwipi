# --------------------------
# Construction des Tests
# --------------------------

if (MPI_C_COMPILER)
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif() 

if (MPI_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif() 

unset(LINK_LIBRARIES)

if (PDM_ENABLE_STATIC)
 list(APPEND LINK_LIBRARIES
      pdm_static
      pdm_mpi_static
      ${COMMON_LIBRARIES})
else ()
 list(APPEND LINK_LIBRARIES
      pdm_shared
      pdm_mpi_shared ${COMMON_LIBRARIES})
endif()

get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)                          


# Test pdm_t_part_cube
# --------------

file(GLOB_RECURSE SOURCES 
     pdm_t_part_dcube.c)

add_executable(pdm_t_part_dcube ${SOURCES})

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(pdm_t_part_dcube
    PROPERTIES
    COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

if ((NOT MPI_Fortran_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(pdm_t_part_dcube
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
endif()

target_include_directories(pdm_t_part_dcube PRIVATE ${CMAKE_SOURCE_DIR}
                                     PRIVATE ${CMAKE_BINARY_DIR}
                                     PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(pdm_t_part_dcube PRIVATE ${PDM_INC})

target_link_libraries(pdm_t_part_dcube ${LINK_LIBRARIES})

install(TARGETS pdm_t_part_dcube
        RUNTIME DESTINATION bin)

add_test (pdm_t_part_cube ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
   ${MPIEXEC_PREFLAGS}
   ${CMAKE_CURRENT_BINARY_DIR}/pdm_t_part_dcube
   ${MPIEXEC_POSTFLAGS})

message(STATUS "Add test pdm_t_part_cube - Succes")

# Test 3D (plusieurs sous-domaines par proc)
# ------------------------------------------

unset(SOURCES)
file(GLOB_RECURSE SOURCES 
     pdm_t_writer_ensight_3d.c)

add_executable(pdm_t_writer_ensight_3d ${SOURCES})

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(pdm_t_writer_ensight_3d
    PROPERTIES
    COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

if ((NOT MPI_Fortran_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(pdm_t_writer_ensight_3d
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
endif()

target_include_directories(pdm_t_writer_ensight_3d PRIVATE ${CMAKE_SOURCE_DIR}
                                      PRIVATE ${CMAKE_BINARY_DIR}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR/../src/msg}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(pdm_t_writer_ensight_3d PRIVATE ${PDM_INC})
target_link_libraries(pdm_t_writer_ensight_3d ${LINK_LIBRARIES})

install(TARGETS pdm_t_writer_ensight_3d
        RUNTIME DESTINATION bin)

# Command to launch by CTest.

add_test (pdm_t_writer_ensight_3d ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
   ${MPIEXEC_PREFLAGS}
   ${CMAKE_CURRENT_BINARY_DIR}/pdm_t_writer_ensight_3d
   ${MPIEXEC_POSTFLAGS})

# Test Coarse mesh surf
# ------------------------------------------

unset(SOURCES)
file(GLOB_RECURSE SOURCES 
     pdm_t_plane_agglo.c)

add_executable(pdm_t_plane_agglo ${SOURCES})

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(pdm_t_plane_agglo
    PROPERTIES
    COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

if ((NOT MPI_Fortran_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(pdm_t_plane_agglo
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
endif()

target_include_directories(pdm_t_plane_agglo PRIVATE ${CMAKE_SOURCE_DIR}
                                      PRIVATE ${CMAKE_BINARY_DIR}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR/../src/msg}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(pdm_t_plane_agglo PRIVATE ${PDM_INC})
target_link_libraries(pdm_t_plane_agglo ${LINK_LIBRARIES})

install(TARGETS pdm_t_plane_agglo
        RUNTIME DESTINATION bin)

# Command to launch by CTest.

add_test (pdm_t_plane_agglo ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
   ${MPIEXEC_PREFLAGS}
   ${CMAKE_CURRENT_BINARY_DIR}/pdm_t_plane_agglo
   ${MPIEXEC_POSTFLAGS})

# Test Gnum
# ------------------------------------------

unset(SOURCES)
file(GLOB_RECURSE SOURCES 
     pdm_t_gnum.c)

add_executable(pdm_t_gnum ${SOURCES})

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(pdm_t_gnum
    PROPERTIES
    COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

if ((NOT MPI_Fortran_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(pdm_t_gnum
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
endif()

target_include_directories(pdm_t_gnum PRIVATE ${CMAKE_SOURCE_DIR}
                                      PRIVATE ${CMAKE_BINARY_DIR}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR/../src/msg}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(pdm_t_gnum PRIVATE ${PDM_INC})
target_link_libraries(pdm_t_gnum ${LINK_LIBRARIES})

install(TARGETS pdm_t_gnum
        RUNTIME DESTINATION bin)

# Command to launch by CTest.

add_test (pdm_t_gnum ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
   ${MPIEXEC_PREFLAGS}
   ${CMAKE_CURRENT_BINARY_DIR}/pdm_t_gnum
   ${MPIEXEC_POSTFLAGS})
