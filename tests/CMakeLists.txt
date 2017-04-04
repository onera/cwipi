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

if (CWP_ENABLE_STATIC)
 list(APPEND LINK_LIBRARIES
      cwp_static
      pdm_static
      pdm_mpi_static
      ${COMMON_LIBRARIES})

 if (CWP_ENABLE_Fortran)
    list(APPEND LINK_LIBRARIES cwpf_static)
  endif()
else ()
 list(APPEND LINK_LIBRARIES
      cwp_shared
      pdm_shared
      pdm_mpi_shared ${COMMON_LIBRARIES})
 if (CWP_ENABLE_Fortran)
    list(APPEND LINK_LIBRARIES cwpf_shared)
 endif()
endif()

get_target_property(TEST_INC cwp INCLUDE_DIRECTORIES)                          


# Test list
# ---------

set (LIST_TESTS "")
list(APPEND LIST_TESTS "c_linear_coupling")

# Add test
# --------

include(TestCreate)
test_c_create (c_linear_coupling)
#test_c_create (c_surf_coupling_P1P1_micron_polygon) (remplacer "cree maillage polygone" par la fonction dans pdm




