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

list(APPEND LINK_LIBRARIES spaceBasis)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})


get_target_property(TEST_INC spaceBasis INCLUDE_DIRECTORIES)                          

file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/Mesh2D DESTINATION ${CMAKE_CURRENT_BINARY_DIR}) 

# Test list
# ---------

set (LIST_TESTS "")
list(APPEND LIST_TESTS "spaceBasis")

# Add test
# --------

SET (CTEST_ENVIRONMENT
       "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}")


include(TestCreate)
test_fortran_create (basesNodales)

#add_test( basesNodales basesNodales 2)
