
# C
test_c_create (c_new_api_deformable_sol  2 CWP_LIST_TEST CWP_LIST_NRANK)
test_c_create (c_new_api_polygon_sol     2 CWP_LIST_TEST CWP_LIST_NRANK)

test_c_create (c_new_api_deformable_ex   2 CWP_LIST_TEST CWP_LIST_NRANK)
test_c_create (c_new_api_polygon_ex      2 CWP_LIST_TEST CWP_LIST_NRANK)

# fortran
if (CWP_ENABLE_Fortran)
  test_fortran_create (fortran_new_api_deformable_sol 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_fortran_create (fortran_new_api_polygon_sol    2 CWP_LIST_TEST CWP_LIST_NRANK)

  test_fortran_create (fortran_new_api_deformable_ex  2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_fortran_create (fortran_new_api_polygon_ex     2 CWP_LIST_TEST CWP_LIST_NRANK)
endif()

# python
if (CWP_ENABLE_PYTHON_BINDINGS)
  test_python_create (python_new_api_polygon_sol 2 CWP_LIST_TEST CWP_LIST_NRANK)

  test_python_create (python_new_api_polygon_ex  2 CWP_LIST_TEST CWP_LIST_NRANK)
endif()

set (CWP_LIST_TEST  ${CWP_LIST_TEST} PARENT_SCOPE )
set (CWP_LIST_NRANK ${CWP_LIST_NRANK} PARENT_SCOPE )
