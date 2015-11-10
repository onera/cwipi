# Liste de tous les sources
# -------------------------

file(GLOB_RECURSE SOURCES 
     *.[chfF]
     *.[fF]90
)


add_executable(BasesNodales ${SOURCES})


# Recherhe des header
# -------------------

target_include_directories(BasesNodales
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_BINARY_DIR}
                           PUBLIC    ${CMAKE_CURRENT_SOURCE_DIR}
                           )


if (MKL_FOUND)
   target_include_directories(BasesNodales PRIVATE ${MKL_INCLUDE_DIRS})
   target_link_libraries(BasesNodales ${MKL_LIBRARIES})
endif (MKL_FOUND)


install(TARGETS BasesNodales RUNTIME DESTINATION bin)



