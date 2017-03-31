
# Sources list
# ------------

file(GLOB_RECURSE SOURCES 
     bft/*.[ch]
     fvm/*.[ch]
     *.[ch]
     *.[ch]xx)

if (NOT CWP_HAVE_NOT_FORTRAN_IN_C)
   list(APPEND SOURCES 
        ${CMAKE_SOURCE_DIR}/src/printfort.f90
        ${CMAKE_SOURCE_DIR}/src/callfortinterpfct.f90
        ${CMAKE_SOURCE_DIR}/src/cwipi_set_output_logical_unit.f90)
endif()

message ("${SOURCES}")

file(GLOB_RECURSE SOURCESF 
     cwipi_f.f90)

add_library(cwp OBJECT ${SOURCES})

add_library(cwpf OBJECT ${SOURCESF})

target_include_directories(cwp
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE  ${CMAKE_BINARY_DIR}
                           PUBLIC  ${CMAKE_SOURCE_DIR}
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/fvm
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/bft)

target_include_directories(cwpf
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_BINARY_DIR}
                           PUBLIC    ${CMAKE_SOURCE_DIR}
                           PUBLIC    ${CMAKE_CURRENT_SOURCE_DIR})
# Liste des header
# ----------------

set (LIST_HEADER "")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/cwipi.h")
list(APPEND LIST_HEADER "${CMAKE_BINARY_DIR}/cwipi_config.h")

if (CWP_ENABLE_STATIC_Fortran_INTERFACE OR CWP_ENABLE_SHARED_Fortran_INTERFACE)

  set (LIST_HEADERF "")
  list(APPEND LIST_HEADERF "${CMAKE_BINARY_DIR}/cwipi_configf.h")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/cwipi.mod")

endif ()

# Options du link
# ---------------

unset(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})

message(STATUS "Add library cwp - Succes")

if (CWP_ENABLE_SHARED)
  set_property(TARGET ${cwp} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(cwp_shared SHARED  $<TARGET_OBJECTS:cwp>)
  set_target_properties(cwp_shared PROPERTIES OUTPUT_NAME "cwp")
  set_target_properties(cwp_shared PROPERTIES VERSION ${CWP_VERSION})
  target_link_libraries(cwp_shared ${LINK_LIBRARIES} pdm_shared)
  target_link_libraries(cwp_shared ${LINK_LIBRARIES} pdm_mpi_shared)
  set_target_properties(cwp_shared PROPERTIES
                        PUBLIC_HEADER "${LIST_HEADER}")

  get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)                          
  target_include_directories(cwp_shared PUBLIC ${CWP_INC})

endif()

if (CWP_ENABLE_STATIC)
  set_property(TARGET ${cwp} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(cwp_static STATIC  $<TARGET_OBJECTS:cwp>)
  set_target_properties(cwp_static PROPERTIES OUTPUT_NAME "cwp")
  set_target_properties(cwp_static PROPERTIES VERSION ${CWP_VERSION})
  target_link_libraries(cwp_static ${LINK_LIBRARIES} pdm_static)
  target_link_libraries(cwp_static ${LINK_LIBRARIES} pdm_mpi_static)
  set_target_properties(cwp_static PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADER}")
  get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)
  target_include_directories(cwp_static PUBLIC ${CWP_INC})

endif()

if (CWP_ENABLE_SHARED_Fortran_INTERFACE)
  set_property(TARGET ${cwpf} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(cwpf_shared SHARED  $<TARGET_OBJECTS:cwpf>)
  set_target_properties(cwpf_shared PROPERTIES OUTPUT_NAME "cwpf")
  set_target_properties(cwpf_shared PROPERTIES VERSION ${CWP_VERSION})
  target_link_libraries(cwpf_shared cwp_shared)
  set_target_properties(cwpf_shared PROPERTIES
                            PUBLIC_HEADER "${LIST_HEADERF}")

  get_target_property(CWP_INC cwpf INCLUDE_DIRECTORIES)
  target_include_directories(cwpf_shared PUBLIC ${CWP_INC})

endif()

if (CWP_ENABLE_STATIC_Fortran_INTERFACE)
  set_property(TARGET ${cwpf} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(cwpf_static STATIC  $<TARGET_OBJECTS:cwpf>)
  set_target_properties(cwpf_static PROPERTIES OUTPUT_NAME "cwpf")
  set_target_properties(cwpf_static PROPERTIES VERSION ${CWP_VERSION})

  target_link_libraries(cwpf_static cwp_static)
  set_target_properties(cwpf_static PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADERF}")

  get_target_property(CWP_INC cwpf INCLUDE_DIRECTORIES)
  target_include_directories(cwpf_static PUBLIC ${CWP_INC})

endif()

# Spécification de l'installation
# -------------------------------

set (_targets "")

if (CWP_ENABLE_SHARED)
  list (APPEND _targets cwp_shared)
endif()

if (CWP_ENABLE_STATIC)
  list (APPEND _targets cwp_static)
endif()

if (CWP_ENABLE_SHARED_Fortran_INTERFACE)
  list (APPEND _targets cwpf_shared)
endif()

if (CWP_ENABLE_STATIC_Fortran_INTERFACE)
  list (APPEND _targets cwpf_static)
endif()

install(TARGETS ${_targets}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        ARCHIVE DESTINATION lib)
