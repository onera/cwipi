if (MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif() 

if (MPI_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif() 

if (MPI_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif(MPI_Fortran_COMPILER) 

# Flags pour MPI
# --------------

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
   set (CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
endif()

if ((NOT MPI_CXX_COMPILER) AND MPI_CXX_COMPILE_FLAGS)
   set (CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
endif()

if ((NOT MPI_Fortran_COMPILER) AND MPI_Fortran_COMPILE_FLAGS)
   set (CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
endif()


# Sources list
# ------------

if (CWP_ENABLE_NEW_API)
 if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/new") 
  file(GLOB SOURCES 
       bft/*.[ch]
       fvm/*.[ch]
       new/*.[ch]
       new/*.[ch]xx
       *.[ch]
       *.[ch]xx)

  if (CWP_NEW_TO_OLD_API)
    list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/cwipi.cxx")
    list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/cwipi.h")
  else()
    list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/new/cwipi_newtoold.cxx")
    list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/new/cwipi_newtoold.h")
  endif()

  get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)
 endif()
else()
  file(GLOB SOURCES 
       bft/*.[ch]
       fvm/*.[ch]
       *.[ch]
       *.[ch]xx)
endif()


if (CWP_ENABLE_Fortran)
   list(APPEND SOURCES fortran/callfortinterpfct.f90)
   list(APPEND SOURCES fortran/callforthointerpfct.f90)
   list(APPEND SOURCES fortran/callfortlocinterpfct.f90)
  endif()

file(GLOB SOURCESF
    fortran/cwipi_cf.cxx
    fortran/cwipi_f.f90
    fortran/cwipi_printfort_cf.cxx
    fortran/cwipi_printfort.f90
    fortran/cwipi_set_output_logical_unit.f90
    fortran/cwp_output_fortran_unit_set.f90
    fortran/cwp_printfort_cf.cxx
    fortran/printfort.f90
    fortran/mod_fvmc_ho_basis.f90
    new/cwp_f.f90)

add_library(cwp OBJECT ${SOURCES})

target_include_directories(cwp
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE  ${CMAKE_BINARY_DIR}
                           PUBLIC  ${CMAKE_SOURCE_DIR}
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/new
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/fvm
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/bft)


if (CWP_ENABLE_NEW_API)
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/new") 
    add_dependencies(cwp pdm)
    target_include_directories(cwp
                               PRIVATE ${PDM_INC}) 
  endif ()
endif ()

#if (CWP_ENABLE_SPACE_NODAL_BASIS) 
#  add_dependencies(cwp spaceBasis)
#  target_include_directories(cwp PRIVATE ${CMAKE_SOURCE_DIR}/spaceBasis/src)
#  target_include_directories(cwp PRIVATE ${CMAKE_SOURCE_DIR}/spaceBasis)
#endif()

if (BLAS_FOUND)
  string(REPLACE " " ";" BLAS_COMPILER_FLAGS "${BLAS_COMPILER_FLAGS}")
  target_compile_options (cwp PRIVATE ${BLAS_COMPILER_FLAGS})
  if (BLAS_INCLUDE_DIRS) 
    target_include_directories(cwp PRIVATE ${BLAS_INCLUDE_DIRS})
  endif()  
endif()

if (CWP_ENABLE_Fortran)
 add_library(cwpf OBJECT ${SOURCESF})

 target_include_directories(cwpf
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_BINARY_DIR}
                           PUBLIC    ${CMAKE_SOURCE_DIR}
                           PUBLIC    ${MPI_Fortran_ADDITIONAL_INCLUDE_DIRS}
                           PUBLIC    ${CMAKE_CURRENT_SOURCE_DIR}
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/new
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/fvm
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/bft) 

  add_dependencies(cwpf cwp)        

if (CWP_ENABLE_SPACE_NODAL_BASIS)
  add_dependencies(cwpf spaceBasis)        
endif()

if (CWP_ENABLE_NEW_API)
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/new")
    add_dependencies(cwpf pdmf)
    target_include_directories(cwpf PRIVATE ${PDM_INC}) 
  endif()
endif()
endif()

# Liste des header
# ----------------

set (LIST_HEADER "")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/cwipi.h")
list(APPEND LIST_HEADER "${CMAKE_BINARY_DIR}/cwipi_config.h")

if (CWP_ENABLE_Fortran)

  set (LIST_HEADERF "")
  list(APPEND LIST_HEADERF "${CMAKE_BINARY_DIR}/cwipi_configf.h")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/cwipi.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/cwp.mod")

endif ()

# Options du link
# ---------------

unset(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})

# Add Fortran dependencies if Fortran is called in C 
# --------------------------------------------------

if ((NOT CWP_HAVE_NOT_FORTRAN_IN_C) AND FORTRAN_LIBRARIES)
    list(APPEND LINK_LIBRARIES ${FORTRAN_LIBRARIES})
endif()

if ((NOT CWP_HAVE_NOT_FORTRAN_IN_C) AND FORTRAN_LIBRARIES_FLAG)
   set_target_properties(cwp PROPERTIES COMPILE_FLAGS ${FORTRAN_LIBRARIES_FLAG})
endif ()

if (MKL_FOUND)
  list(APPEND LINK_LIBRARIES ${MKL_LIBRARIES})
else ()
if (BLAS_FOUND)
  list(APPEND LINK_LIBRARIES ${BLAS_LIBRARIES})
endif()

if (LAPACK_FOUND)
  list(APPEND LINK_LIBRARIES ${LAPACK_LIBRARIES})
endif()
endif()


if (CWP_ENABLE_SHARED)
  set_property(TARGET ${cwp} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(cwp_shared SHARED  $<TARGET_OBJECTS:cwp>)
  add_dependencies(cwp_shared cwp)
  set_target_properties(cwp_shared PROPERTIES OUTPUT_NAME "cwp")
  set_target_properties(cwp_shared PROPERTIES VERSION ${CWP_VERSION})
  target_link_libraries(cwp_shared ${LINK_LIBRARIES}) 
  if (CWP_ENABLE_NEW_API)
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/new")
      target_link_libraries(cwp_shared pdm_shared)
      target_link_libraries(cwp_shared pdm_mpi_shared)
    endif()
  endif()
  set_target_properties(cwp_shared PROPERTIES
                        PUBLIC_HEADER "${LIST_HEADER}")

  get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)                          
  target_include_directories(cwp_shared PUBLIC ${CWP_INC})
  if (CWP_ENABLE_SPACE_NODAL_BASIS)
     message("add dependencies spacebasis_shared")
     target_link_libraries(cwp_shared spaceBasis_shared)
  endif()

  if (CWP_ENABLE_Fortran)
    set_property(TARGET ${cwpf} PROPERTY POSITION_INDEPENDENT_CODE 1)
    add_library(cwpf_shared SHARED  $<TARGET_OBJECTS:cwpf>)
    add_dependencies(cwpf_shared cwpf)
    set_target_properties(cwpf_shared PROPERTIES OUTPUT_NAME "cwpf")
    set_target_properties(cwpf_shared PROPERTIES VERSION ${CWP_VERSION})
    target_link_libraries(cwpf_shared cwp_shared ${CXX_LIBRARIES})
    set_target_properties(cwpf_shared PROPERTIES
                            PUBLIC_HEADER "${LIST_HEADERF}")

    get_target_property(CWP_INC cwpf INCLUDE_DIRECTORIES)
    target_include_directories(cwpf_shared PUBLIC ${CWP_INC})
    if (CWP_ENABLE_SPACE_NODAL_BASIS)
     target_link_libraries(cwpf_shared spaceBasis_shared)        
    endif()

  endif()
endif()

if (CWP_ENABLE_STATIC)
  set_property(TARGET ${cwp} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(cwp_static STATIC  $<TARGET_OBJECTS:cwp>)
  add_dependencies(cwp_static cwp)
  set_target_properties(cwp_static PROPERTIES OUTPUT_NAME "cwp")
  set_target_properties(cwp_static PROPERTIES VERSION ${CWP_VERSION})
  target_link_libraries(cwp_static ${LINK_LIBRARIES}) 
  if (CWP_ENABLE_NEW_API)
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/new")
      target_link_libraries(cwp_static pdm_static)
      target_link_libraries(cwp_static pdm_mpi_static)
    endif()
  endif()  
  set_target_properties(cwp_static PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADER}")
  get_target_property(CWP_INC cwp INCLUDE_DIRECTORIES)
  target_include_directories(cwp_static PUBLIC ${CWP_INC})
  if (CWP_ENABLE_SPACE_NODAL_BASIS)
     target_link_libraries(cwp_static spaceBasis_static)        
  endif()

  if (CWP_ENABLE_Fortran)
    set_property(TARGET ${cwpf} PROPERTY POSITION_INDEPENDENT_CODE 1)
    add_library(cwpf_static STATIC  $<TARGET_OBJECTS:cwpf>)
    add_dependencies(cwpf_static cwpf)
    set_target_properties(cwpf_static PROPERTIES OUTPUT_NAME "cwpf")
    set_target_properties(cwpf_static PROPERTIES VERSION ${CWP_VERSION})

    target_link_libraries(cwpf_static cwp_static)
    set_target_properties(cwpf_static PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADERF}")

    get_target_property(CWP_INC cwpf INCLUDE_DIRECTORIES)
    target_include_directories(cwpf_static PUBLIC ${CWP_INC})
    if (CWP_ENABLE_SPACE_NODAL_BASIS)
      target_link_libraries(cwp_static spaceBasis_static)        
    endif()

  endif()

endif()

# Spécification de l'installation
# -------------------------------

set (_targets "")

if (CWP_ENABLE_SHARED)
  list (APPEND _targets cwp_shared)
  if (CWP_ENABLE_Fortran)
    list (APPEND _targets cwpf_shared)
  endif()
endif()

if (CWP_ENABLE_STATIC)
  list (APPEND _targets cwp_static)
  if (CWP_ENABLE_Fortran)
    list (APPEND _targets cwpf_static)
  endif()
endif()

install(TARGETS ${_targets}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        ARCHIVE DESTINATION lib)
