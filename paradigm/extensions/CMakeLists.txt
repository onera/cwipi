if (PDM_ENABLE_EXTENSION_PDMA)
  if (EXISTS ${CMAKE_SOURCE_DIR}/extensions/paradigma/CMakeLists.txt)
    add_subdirectory(paradigma)
  else ()
     message (FATAL_ERROR "Sources of ParaDiGMA extension are not into the extensions/paradigma directory. Please, donwload them from the git command :
     git submodule update --init")
  endif ()
endif ()
