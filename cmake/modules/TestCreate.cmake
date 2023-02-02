################################################################################
#  
# test_c_create and test_fortran_create add a Fortran test or a C test
#
# They uses LINK_LIBRARIES and TEST_INC variables 
#
################################################################################

function(test_c_create name n_proc)
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
   add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
             ${MPIEXEC_PREFLAGS}
             ${CMAKE_CURRENT_BINARY_DIR}/${name}
             ${MPIEXEC_POSTFLAGS})

  if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    set_property(TEST ${name} PROPERTY ENVIRONMENT "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/scripts/asan/asan.supp")
  endif()

endfunction()

function(test_fortran_create name n_proc)
   if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${name}.f90")
     add_executable(${name} "${name}.f90")
   else ()
     add_executable(${name} "${name}.F90")
   endif()
   if ((NOT MPI_Fortran_COMPILER) AND MPI_C_COMPILE_FLAGS)
     set_target_properties(${name}
                           PROPERTIES
                           COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
     target_include_directories(${name} PRIVATE ${MPI_Fortran_INCLUDE_PATH})                      
   endif()
   target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}
                                      PRIVATE ${CMAKE_BINARY_DIR}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
   target_include_directories(${name} PRIVATE ${TEST_INC})
   target_link_libraries(${name} ${LINK_LIBRARIES})
   set_target_properties(${name} PROPERTIES LINKER_LANGUAGE "Fortran")
   install(TARGETS ${name} RUNTIME DESTINATION bin)
   add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
             ${MPIEXEC_PREFLAGS}
             ${CMAKE_CURRENT_BINARY_DIR}/${name}
             ${MPIEXEC_POSTFLAGS})

  if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    set_property(TEST ${name} PROPERTY ENVIRONMENT "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/scripts/asan/asan.supp")
  endif()
endfunction()

function(test_python_create name n_proc)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${name}.py ${CMAKE_CURRENT_BINARY_DIR}/${name}.py)  
  file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/${name}.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
            ${MPIEXEC_PREFLAGS}
            ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/${name}.py
            ${MPIEXEC_POSTFLAGS})

  set (LIST_PYTHON_TEST_ENV "")

  if (PYTHON_TEST_ENV)
    list(APPEND LIST_PYTHON_TEST_ENV "${PYTHON_TEST_ENV}")
  endif()

  if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    list(APPEND LIST_PYTHON_TEST_ENV "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/script/asan/asan.supp")
  endif()

  if(DEFINED ENV{PYTHONPATH})
    list(APPEND LIST_PYTHON_TEST_ENV "PYTHONPATH=${CMAKE_BINARY_DIR}/Cython:$ENV{PYTHONPATH}")
  else()
    list(APPEND LIST_PYTHON_TEST_ENV "PYTHONPATH=${CMAKE_BINARY_DIR}/Cython")
  endif()

  if (LIST_PYTHON_TEST_ENV)
    set_property(TEST ${name} PROPERTY ENVIRONMENT "${LIST_PYTHON_TEST_ENV}")
  endif()

endfunction()
