# Liste de tous les sources
# -------------------------

include_directories(${PROJECT_SOURCE_DIR}/src)

#Bring the headers, such as Student.h into the project
#include_directories(include)


file(GLOB_RECURSE SOURCES
	${PROJECT_SOURCE_DIR}/src/*.F90
	${PROJECT_SOURCE_DIR}/src/*.f90
	${PROJECT_SOURCE_DIR}/src/*.c
)


#add_executable(BasesNodales ${SOURCES})

#set_target_properties(BasesNodales PROPERTIES LINKER_LANGUAGE Fortran)

#target_include_directories(BasesNodales
#	INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
#	PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
#	PRIVATE   ${CMAKE_BINARY_DIR}
#	PUBLIC    ${CMAKE_CURRENT_SOURCE_DIR}
#)




#add_library(spaceBasis STATIC ${SOURCES})

#add_library(spaceBasis SHARED ${SOURCES})

add_library(spaceBasis OBJECT ${SOURCES})



set (LIST_HEADERF "")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/basewedge.mod")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/basesimplex1d.mod")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/basesimplex2d.mod")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/basesimplex3d.mod")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pyramidrule.mod")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/basewedge.mod")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/basesimplextools.mod")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/moddeterminant.mod")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/space_ldlt.mod")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/table_tet_mesh.mod")
list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/spacebasis.h")

#message("LIST_HEADERF:       ${LIST_HEADERF}")


set_target_properties(spaceBasis PROPERTIES PUBLIC_HEADER "${LIST_HEADERF}")

get_target_property(spaceBasis_INC spaceBasis INCLUDE_DIRECTORIES)
target_include_directories(spaceBasis PUBLIC ${CWP_INC})



target_include_directories(spaceBasis
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE  ${CMAKE_BINARY_DIR}
                           PUBLIC  ${CMAKE_SOURCE_DIR}
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories (spaceBasis
                           INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
#                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})





#target_include_directories(BasesNodales PRIVATE ${BLAS_INCLUDE_DIRS})
#target_link_libraries(BasesNodales ${BLAS_LIBRARIES})

#target_include_directories(BasesNodales PRIVATE ${LAPACK_INCLUDE_DIRS})
#target_link_libraries(BasesNodales ${LAPACK_LIBRARIES})

target_include_directories(spaceBasis PRIVATE ${BLAS_INCLUDE_DIRS})

target_include_directories(spaceBasis PRIVATE ${LAPACK_INCLUDE_DIRS})

set (_targets "")

if (SNB_ENABLE_SHARED)
  set_property(TARGET ${spaceBasis} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(spaceBasis_shared SHARED  $<TARGET_OBJECTS:spaceBasis>)
  add_dependencies(spaceBasis_shared spaceBasis)
  set_target_properties(spaceBasis_shared PROPERTIES OUTPUT_NAME "spaceBasis")
  set_target_properties(spaceBasis_shared PROPERTIES VERSION ${BASISNODALES_VERSION})
  if (MKL_FOUND)
    target_link_libraries(spaceBasis_shared ${MKL_LIBRARIES})
  else ()
  if (BLAS_FOUND)
    target_link_libraries(spaceBasis_shared ${BLAS_LIBRARIES})
  endif()
  if (LAPACK_FOUND)
    target_link_libraries(spaceBasis_shared ${LAPACK_LIBRARIES})
   endif()
  endif()
  list (APPEND _targets spaceBasis_shared)
endif ()

if (SNB_ENABLE_STATIC)
  set_property(TARGET ${spaceBasis} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(spaceBasis_static STATIC  $<TARGET_OBJECTS:spaceBasis>)
  add_dependencies(spaceBasis_static spaceBasis)
  set_target_properties(spaceBasis_static PROPERTIES OUTPUT_NAME "spaceBasis")
  set_target_properties(spaceBasis_static PROPERTIES VERSION ${BASISNODALES_VERSION})
  if (MKL_FOUND)
    target_link_libraries(spaceBasis_static ${MKL_LIBRARIES})
  else ()
  if (BLAS_FOUND)
    target_link_libraries(spaceBasis_static ${BLAS_LIBRARIES})
  endif()
  if (LAPACK_FOUND)
    target_link_libraries(spaceBasis_static ${LAPACK_LIBRARIES})
   endif()
  endif()
  list (APPEND _targets spaceBasis_static)
endif ()

#list (APPEND _targets BasesNodales)

install(TARGETS ${_targets}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        ARCHIVE DESTINATION lib)
