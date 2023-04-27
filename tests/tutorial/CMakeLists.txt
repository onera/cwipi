
# C

include(TestCreate)

if (CWP_ENABLE_Fortran)
  set (exe_tuto_ex "")
  set (exe_tuto_sol "")
  set (nproc_tuto "")

  list (APPEND exe_tuto_ex c_new_api_deformable_ex.c fortran_new_api_deformable_ex.F90)

  list (APPEND exe_tuto_sol c_new_api_deformable_sol.c fortran_new_api_deformable_sol.F90)

  list (APPEND nproc_tuto 2 2) 


  message("names1 : ${exe_tuto_ex}")
  message("n_procs1 : ${nproc_tuto}")

  test_create ("${exe_tuto_ex}" "${nproc_tuto}" CWP_LIST_TEST CWP_LIST_NRANK)
  test_create ("${exe_tuto_sol}" "${nproc_tuto}" CWP_LIST_TEST CWP_LIST_NRANK)
endif()

test_c_create (c_new_api_polygon_sol     2 CWP_LIST_TEST CWP_LIST_NRANK)
test_c_create (c_new_api_polygon_ex      2 CWP_LIST_TEST CWP_LIST_NRANK)

# fortran
if (CWP_ENABLE_Fortran)
  test_fortran_create (fortran_new_api_polygon_sol    2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_fortran_create (fortran_new_api_polygon_ex     2 CWP_LIST_TEST CWP_LIST_NRANK)
endif()

# python
if (CWP_ENABLE_PYTHON_BINDINGS)
  test_python_create (python_new_api_polygon_sol 2 CWP_LIST_TEST NRANK)
  test_python_create (python_new_api_polygon_ex  2 CWP_LIST_TEST NRANK)

  if (CWP_ENABLE_Fortran)
    test_c_create       (c_new_api_callback_sendrecv_sol       2 CWP_LIST_TEST NRANK)
    test_fortran_create (fortran_new_api_callback_sendrecv_sol 2 CWP_LIST_TEST NRANK)
    test_python_create  (python_new_api_callback_sendrecv_sol  2 CWP_LIST_TEST NRANK)
  endif()
endif()

set (CWP_LIST_TEST  ${CWP_LIST_TEST} PARENT_SCOPE )
set (CWP_LIST_NRANK ${CWP_LIST_NRANK} PARENT_SCOPE )
