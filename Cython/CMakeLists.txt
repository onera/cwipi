if ((Mpi4Py_VERSION VERSION_GREATER 2.0.0) OR
    (Mpi4Py_VERSION VERSION_EQUAL 2.0.0))
  set (mpi4pylibmpi libmpi)
else ()
  set (mpi4pylibmpi mpi_c)
endif()

set (include_dirs "")
list (APPEND include_dirs ${Mpi4Py_INCLUDE_DIR})
list (APPEND include_dirs ${NumPy_INCLUDE_DIR})
list (APPEND include_dirs ${CMAKE_CURRENT_BINARY_DIR})
list (APPEND include_dirs ${CMAKE_CURRENT_SOURCE_DIR})
list (APPEND include_dirs ${CMAKE_BINARY_DIR})

unset(LINK_LIBRARIES)

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  list(APPEND LINK_LIBRARIES "m")
endif()

set(PYX_NPY_G_NUM np.int64_t)
set(PYX_NPY_G_NUM_DTYPE np.int64)
set(CWP_G_NUM_NPY_INT np.NPY_INT64)
set(PYX_G_NUM long)
set(PYX_L_NUM int)

add_subdirectory(cwipi)

add_subdirectory(pycwp)

if (CWP_ENABLE_CLIENT_SERVER)
  add_subdirectory(pycwpclt)
endif()
