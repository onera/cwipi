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

get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)                          

# Test create function
# --------------------

function(test_c_create name)
   add_executable(${name} "${name}.c")
   if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
     set_target_properties(${name}
                           PROPERTIES
                           COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
   endif()
   target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}
                                      PRIVATE ${CMAKE_BINARY_DIR}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
   target_include_directories(${name} PRIVATE ${CWP_INC})
   target_link_libraries(${name} ${LINK_LIBRARIES})
   install(TARGETS ${name} RUNTIME DESTINATION bin)
   add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
             ${MPIEXEC_PREFLAGS}
             ${CMAKE_CURRENT_BINARY_DIR}/${name}
             ${MPIEXEC_POSTFLAGS})
endfunction()

function(test_fortran_create name)
   add_executable(${name} "${name}.f90")
   if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
     set_target_properties(${name}
                           PROPERTIES
                           COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
   endif()
   if ((NOT MPI_Fortran_COMPILER) AND MPI_C_COMPILE_FLAGS)
     set_target_properties(pdm_t_writer_ensight_3d
                           PROPERTIES
                           COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
   endif()
   set_target_properties(${name} PROPERTIES COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
   target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}
                                      PRIVATE ${CMAKE_BINARY_DIR}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
   target_include_directories(${name} PRIVATE ${CWP_INC})
   target_link_libraries(${name} ${LINK_LIBRARIES})
   install(TARGETS ${name} RUNTIME DESTINATION bin)
   add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
             ${MPIEXEC_PREFLAGS}
             ${CMAKE_CURRENT_BINARY_DIR}/${name}
             ${MPIEXEC_POSTFLAGS})
endfunction()


# Test list
# ---------

set (LIST_TESTS "")
list(APPEND LIST_TESTS "c_linear_coupling")

# Add test
# --------

test_c_create (c_linear_coupling)
#foreach(ctest IN ${LIST_TESTS})
#  test_c_create (ctest)
#endforeach()



