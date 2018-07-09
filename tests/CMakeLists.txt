# ------------
# Build tests
# ------------

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

if (CWP_ENABLE_SHARED)

   if (CWP_ENABLE_Fortran)
     if (EXISTS ${CWP_SOURCE_DIR}/paradigm)
       list(APPEND LINK_LIBRARIES pdmf_shared)
     endif()  
     list(APPEND LINK_LIBRARIES cwpf_shared)
   endif()
   if (EXISTS ${CWP_SOURCE_DIR}/paradigm)
     list(APPEND LINK_LIBRARIES
          pdm_mpi_shared
          pdm_shared)
   endif()          
   list(APPEND LINK_LIBRARIES cwp_shared)
   if (CWP_ENABLE_Fortran)
    list(APPEND LINK_LIBRARIES cwpf_shared)
   endif()
   list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})

# list(APPEND LINK_LIBRARIES
#       cwp_shared)
#  if (EXISTS ${CWP_SOURCE_DIR}/paradigm)
#      list (APPEND LINK_LIBRARIES 
#            pdm_shared
#            pdm_mpi_shared)
#  endif()
#  list(APPEND LINK_LIBRARIES
#       ${COMMON_LIBRARIES})
#
# if (CWP_ENABLE_Fortran)
#   list(APPEND LINK_LIBRARIES cwpf_shared)
# endif()

else ()
   if (CWP_ENABLE_Fortran)
     if (EXISTS ${CWP_SOURCE_DIR}/paradigm)
       list(APPEND LINK_LIBRARIES pdmf_static)
     endif()  
     list(APPEND LINK_LIBRARIES cwpf_static)
   endif()
   if (EXISTS ${CWP_SOURCE_DIR}/paradigm)
     list(APPEND LINK_LIBRARIES
          pdm_mpi_static
          pdm_static)
   endif()          
   list(APPEND LINK_LIBRARIES cwp_static)
   if (CWP_ENABLE_Fortran)
    list(APPEND LINK_LIBRARIES cwpf_static)
   endif()
   list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})

endif()


get_target_property(TEST_INC cwp INCLUDE_DIRECTORIES)                          

file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/meshes DESTINATION ${CMAKE_CURRENT_BINARY_DIR}) 

# Test list
# ---------

set (LIST_TESTS "")
list(APPEND LIST_TESTS "c_linear_coupling")

# Add test
# --------

if (CWP_ENABLE_PYTHON_BINDINGS)
  SET (CTEST_ENVIRONMENT
       "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
       "PYTHONPATH=${CWP_BINARY_DIR}/Cython:${PYTHONPATH}")
else()
  SET (CTEST_ENVIRONMENT
       "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}")
endif()

include(TestCreate)

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/c_new_api.c")
  test_c_create (c_new_api 10)
endif ()  
test_c_create (c_linear_coupling ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_location_triaP2 2)
test_c_create (c_surf_location_curvetriaP2 2)
test_c_create (c_surf_coupling_P1P1_micron_polygon ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_P1P1_polygon ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_P1P1_polygon2 ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_P1P1_polygon_with_external_points ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_P1P1 ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_P1P0_P0P1 ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_P1P0_P0P1_savedisk ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_P1P0_P0P1_savememory ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_P1P1_asynchronous ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_user_target_points ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_user_interpolation ${MPIEXEC_MAX_NUMPROCS})
#c_surf_coupling_user_interpolation2_SOURCES = c_surf_coupling_user_interpolation2.c
test_c_create (c_surf_coupling_vector_exchange ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_output_error ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_coupling_simple_location ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_without_partitioning_coupling ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_without_partitioning_coupling2 ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_surf_sequential_coupling ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_vol_coupling_P1P1 ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_vol_coupling_P1P1_with_external_points ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_vol_coupling_micro_polyhedra_P1P1 ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_vol_coupling_polyhedra_P1P1 ${MPIEXEC_MAX_NUMPROCS})
test_c_create (c_vol_coupling_polyhedra_P1P1_with_external_points ${MPIEXEC_MAX_NUMPROCS})

if (CWP_ENABLE_Fortran)
  test_fortran_create (fortran_surf_P1P1            ${MPIEXEC_MAX_NUMPROCS})
  test_fortran_create (fortran_surf_P1P1_savememory ${MPIEXEC_MAX_NUMPROCS})
  
   list(APPEND LINK_LIBRARIES spaceBasis)
  test_fortran_create  (fortran_surf_TriaP1_PiPj   ${MPIEXEC_MAX_NUMPROCS})
  target_link_libraries(fortran_surf_TriaP1_PiPj   spaceBasis)
   add_dependencies(fortran_surf_TriaP1_PiPj   spaceBasis)

  test_fortran_create  (fortran_surf_TriaP2_PiPj   ${MPIEXEC_MAX_NUMPROCS})
  target_link_libraries(fortran_surf_TriaP2_PiPj   spaceBasis)
endif()

if (CWP_ENABLE_PYTHON_BINDINGS)
  test_python_create (python_api ${MPIEXEC_MAX_NUMPROCS})
endif()
