# Liste de tous les sources
# -------------------------

include_directories(${PROJECT_SOURCE_DIR}/src)
file(GLOB_RECURSE SOURCES
	${PROJECT_SOURCE_DIR}/src/*.F90
	${PROJECT_SOURCE_DIR}/src/*.f90
	${PROJECT_SOURCE_DIR}/src/*.c
)



add_executable(BasesNodales ${SOURCES})

set_target_properties(BasesNodales PROPERTIES LINKER_LANGUAGE Fortran)

target_include_directories(BasesNodales
	INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
	PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
	PRIVATE   ${CMAKE_BINARY_DIR}
	PUBLIC    ${CMAKE_CURRENT_SOURCE_DIR}
)

if (BLAS_FOUND)
	target_include_directories(BasesNodales PRIVATE ${BLAS_INCLUDE_DIRS})
	target_link_libraries(BasesNodales ${BLAS_LIBRARIES})
endif (BLAS_FOUND)

if (LAPACK_FOUND)
	target_include_directories(BasesNodales PRIVATE ${LAPACK_INCLUDE_DIRS})
	target_link_libraries(BasesNodales ${LAPACK_LIBRARIES})
endif (LAPACK_FOUND)


install(TARGETS BasesNodales RUNTIME DESTINATION .)



