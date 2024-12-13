file(READ ${CMAKE_SOURCE_DIR}/cwp_Build.defs.in cwp_build_defs_file)

# TODO modifier fichier en fonction des options cmake
if (CWP_ENABLE_EXTERNAL_PDM)
  string(REPLACE "#PDM_DIR" "PDM_DIR = ${PDM_DIR}" cwp_build_defs_file ${cwp_build_defs_file})
else()
  string(REPLACE "#PDM_DIR" "" cwp_build_defs_file ${cwp_build_defs_file})
endif()


if (MPI_C_COMPILER OR MPI_CXX_COMPILER OR MPI_Fortran_COMPILER)
  set(wrapper_mpi "# MPI Wrapper used\n#-----------------------------------------------------------------------")
  if (MPI_C_COMPILER)
    set(wrapper_mpi "${wrapper_mpi}\nMPI_C       = ${MPI_C_COMPILER}")
  endif()
  if (MPI_CXX_COMPILER)
    set(wrapper_mpi "${wrapper_mpi}\nMPI_CXX     = ${MPI_CXX_COMPILER}")
  endif()
  if (MPI_Fortran_COMPILER)
    set(wrapper_mpi "${wrapper_mpi}\nMPI_Fortran = ${MPI_Fortran_COMPILER}")
  endif()
  string(REPLACE "#MPI_Wrapper" ${wrapper_mpi} cwp_build_defs_file ${cwp_build_defs_file})
else()
  string(REPLACE "#MPI_Wrapper" "" cwp_build_defs_file ${cwp_build_defs_file})
endif ()


file(WRITE ${CMAKE_BINARY_DIR}/cwp_Build.defs.in "${cwp_build_defs_file}")

configure_file(${CMAKE_BINARY_DIR}/cwp_Build.defs.in "${CMAKE_CURRENT_BINARY_DIR}/cwp_Build.defs")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cwp_Build.defs
        DESTINATION include)
