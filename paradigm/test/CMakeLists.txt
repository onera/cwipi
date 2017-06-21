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

if (NOT MPI_C_COMPILER)
  if (MPI_C_LIBRARIES)
    list(APPEND LINK_LIBRARIES ${MPI_C_LIBRARIES})   
  endif()
endif ()

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

include(TestCreate)

test_c_create (pdm_t_part_dcube ${MPIEXEC_MAX_NUMPROCS})
test_c_create (pdm_t_writer_ensight_3d ${MPIEXEC_MAX_NUMPROCS})
test_c_create (pdm_t_plane_agglo ${MPIEXEC_MAX_NUMPROCS})
test_c_create (pdm_t_gnum ${MPIEXEC_MAX_NUMPROCS})
