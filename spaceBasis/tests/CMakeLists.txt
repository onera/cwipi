# ------------
# Build tests
# ------------

#if (MPI_C_COMPILER)
#    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
#endif() 

#if (MPI_Fortran_COMPILER)
#    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
#endif() 

unset(LINK_LIBRARIES)

if (SNB_ENABLE_SHARED)
  list(APPEND LINK_LIBRARIES spaceBasis_shared)
else()
  list(APPEND LINK_LIBRARIES spaceBasis_static)
endif()
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})


get_target_property(TEST_INC spaceBasis INCLUDE_DIRECTORIES)                          

file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/Mesh2D DESTINATION ${CMAKE_CURRENT_BINARY_DIR}) 
file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/Mesh3D DESTINATION ${CMAKE_CURRENT_BINARY_DIR}) 


# Test list
# ---------

set (LIST_TESTS "")
list(APPEND LIST_TESTS "spaceBasis")

# Add test
# --------

SET (CTEST_ENVIRONMENT
       "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}")


include(TestCreate)
#test_fortran_create (basesNodales)
test_fortran_create (test_QuadraturesGauss)
add_dependencies(test_QuadraturesGauss spaceBasis)

test_fortran_create (test_ConnectivitiesOfSides)
add_dependencies(test_ConnectivitiesOfSides spaceBasis)

test_fortran_create (test_Jacobi)
add_dependencies(test_Jacobi spaceBasis)

test_fortran_create (test_Edges00)
add_dependencies(test_Edges00 spaceBasis)

test_fortran_create (test_Edges01)
add_dependencies(test_Edges01 spaceBasis)

test_fortran_create (test_Edges02)
add_dependencies(test_Edges02 spaceBasis)

test_fortran_create (test_Triangles00)
add_dependencies(test_Triangles00 spaceBasis)

test_fortran_create (test_Triangles01)
add_dependencies(test_Triangles01 spaceBasis)

test_fortran_create (test_Triangles02)
add_dependencies(test_Triangles02 spaceBasis)

test_fortran_create (test_Quad00)
add_dependencies(test_Quad00 spaceBasis)

test_fortran_create (test_Tetra00)
add_dependencies(test_Tetra00 spaceBasis)

test_fortran_create (test_Tetra01)
add_dependencies(test_Tetra01 spaceBasis)

test_fortran_create (test_Pyramid00)
add_dependencies(test_Pyramid00 spaceBasis)

test_fortran_create (test_Pyramid01)
add_dependencies(test_Pyramid01 spaceBasis)

#test_fortran_create (test_Pyramid02)
test_fortran_create (test_Pyramid03)
add_dependencies(test_Pyramid03 spaceBasis)



