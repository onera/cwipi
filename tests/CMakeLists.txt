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

if (CWP_ENABLE_SHARED)
 list(APPEND LINK_LIBRARIES
      cwp_shared
      pdm_shared
      pdm_mpi_shared ${COMMON_LIBRARIES})
 if (CWP_ENABLE_Fortran)
    list(APPEND LINK_LIBRARIES cwpf_shared)
 endif()

else ()
   if (CWP_ENABLE_Fortran)
      list(APPEND LINK_LIBRARIES pdmf_static cwpf_static)
   endif()
   list(APPEND LINK_LIBRARIES
      pdm_mpi_static
      pdm_static
      cwp_static)
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
test_c_create (c_linear_coupling)
test_c_create (c_surf_coupling_P1P1_micron_polygon)
test_c_create (c_surf_coupling_P1P1_polygon)
test_c_create (c_surf_coupling_P1P1_polygon2)
test_c_create (c_surf_coupling_P1P1_polygon_with_external_points)
test_c_create (c_surf_coupling_P1P1)
test_c_create (c_surf_coupling_P1P0_P0P1)
test_c_create (c_surf_coupling_P1P0_P0P1_savedisk)
test_c_create (c_surf_coupling_P1P0_P0P1_savememory)
test_c_create (c_surf_coupling_P1P1_asynchronous)
test_c_create (c_surf_coupling_user_target_points)
test_c_create (c_surf_coupling_user_interpolation)
#c_surf_coupling_user_interpolation2_SOURCES = c_surf_coupling_user_interpolation2.c
test_c_create (c_surf_coupling_vector_exchange)
test_c_create (c_surf_coupling_output_error)
test_c_create (c_surf_coupling_simple_location)
test_c_create (c_surf_without_partitioning_coupling)
test_c_create (c_surf_without_partitioning_coupling2)
test_c_create (c_surf_sequential_coupling)
test_c_create (c_vol_coupling_P1P1)
test_c_create (c_vol_coupling_P1P1_with_external_points)
test_c_create (c_vol_coupling_micro_polyhedra_P1P1)
test_c_create (c_vol_coupling_polyhedra_P1P1)
test_c_create (c_vol_coupling_polyhedra_P1P1_with_external_points)

if (CWP_ENABLE_Fortran)
  test_fortran_create (fortran_surf_P1P1)
  test_fortran_create (fortran_surf_P1P1_savememory)
endif()

if (CWP_ENABLE_PYTHON_BINDINGS)
  test_python_create (python_api)
endif()