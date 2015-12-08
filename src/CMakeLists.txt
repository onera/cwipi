# Liste de tous les sources
# -------------------------


include_directories(${PROJECT_SOURCE_DIR}/src)
file(GLOB_RECURSE SOURCES
     ${PROJECT_SOURCE_DIR}/src/*.f90
     ${PROJECT_SOURCE_DIR}/src/*.c
)



add_executable(BasesNodales ${SOURCES})


target_link_libraries(BasesNodales "-static-intel")


#target_link_libraries(space "-arch x86_64"})
#target_link_libraries(space ${MPI_Fortran_LIBRARIES})


set_target_properties(BasesNodales
                      PROPERTIES 
                      LINKER_LANGUAGE Fortran)





# Recherhe des header
# -------------------

target_include_directories(BasesNodales
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_BINARY_DIR}
                           PUBLIC    ${CMAKE_CURRENT_SOURCE_DIR}
                           )



install(TARGETS BasesNodales RUNTIME DESTINATION bin)



