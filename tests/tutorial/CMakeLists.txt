if (CWP_ENABLE_NEW_API)

  # C
  # test_c_create (c_new_api_deformable_sol  2)
  # test_c_create (c_new_api_polygon_sol     2)

  # test_c_create (c_new_api_polygon_ex      2)

  # fortran
  if (CWP_ENABLE_Fortran)
    # test_fortran_create (fortran_new_api_deformable_sol 2)
    # test_fortran_create (fortran_new_api_polygon_sol    2)

    # test_fortran_create (fortran_new_api_polygon_ex     2)
  endif()

  # python
  if (CWP_ENABLE_PYTHON_BINDINGS)
    # test_python_create (python_new_api_polygon_sol 2)

    # test_python_create (python_new_api_polygon_ex  2)
  endif()

endif()
