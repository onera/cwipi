################################################################################
#  
# test_c_create and test_fortran_create add a Fortran test or a C test
#
# They uses LINK_LIBRARIES and TEST_INC variables 
#
################################################################################

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
   target_include_directories(${name} PRIVATE ${TEST_INC})
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
   target_include_directories(${name} PRIVATE ${TEST_INC})
   target_link_libraries(${name} ${LINK_LIBRARIES})
   install(TARGETS ${name} RUNTIME DESTINATION bin)
   add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
             ${MPIEXEC_PREFLAGS}
             ${CMAKE_CURRENT_BINARY_DIR}/${name}
             ${MPIEXEC_POSTFLAGS})
endfunction()

