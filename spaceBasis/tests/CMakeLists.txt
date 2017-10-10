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

list(APPEND LINK_LIBRARIES spaceBasis)
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
test_fortran_create (test_ConnectivitiesOfSides)
test_fortran_create (test_Jacobi)

test_fortran_create (test_Edges00)
test_fortran_create (test_Edges01)

test_fortran_create (test_Triangles00)
test_fortran_create (test_Triangles01)
test_fortran_create (test_Triangles02)

test_fortran_create (test_Quad00)

test_fortran_create (test_Tetra00)
test_fortran_create (test_Tetra01)

test_fortran_create (test_Pyramid00)
test_fortran_create (test_Pyramid01)
#test_fortran_create (test_Pyramid02)
test_fortran_create (test_Pyramid03)


#add_test(basesNodales   basesNodales)
#add_test(test_Triangles test_Triangles)
