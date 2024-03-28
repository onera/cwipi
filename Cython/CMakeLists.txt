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

list(APPEND LINK_LIBRARIES cwp_shared)

set(PYX_NPY_G_NUM np.int64_t)
set(PYX_NPY_G_NUM_DTYPE np.int64)
set(CWP_G_NUM_NPY_INT np.NPY_INT64)
set(PYX_G_NUM long)
set(PYX_L_NUM int)

add_subdirectory(cwipi)

add_subdirectory(pycwp)

add_subdirectory(pycwpt)

if (CWP_ENABLE_CLIENT_SERVER)
  add_subdirectory(pycwpclt)
endif()

if (NOT CWP_ENABLE_EXTERNAL_PDM AND (CWP_ENABLE_TESTS OR CWP_ENABLE_TRAINING))
  list(APPEND LINK_LIBRARIES pdm_shared)
  set(PYX_NPY_G_NUM NPY.int64_t)
  set(PYX_NPY_G_NUM_DTYPE NPY.int64)
  set(PDM_G_NUM_NPY_INT NPY.NPY_INT64)
  set(PYX_PDM_HAS_PARMETIS False)
  set(PYX_PDM_HAS_PTSCOTCH False)
  set(PYX_PDM_HAS_PDMA     False)
  set(COMMENT_NOT_PDMA      "\#")
  file(COPY        ${PDM_SOURCE_DIR}/Cython/Pypdm
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY        Pypdm/CMakeLists.txt
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Pypdm)
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/Pypdm
                   ${CMAKE_CURRENT_BINARY_DIR}/Pypdm)
endif()
