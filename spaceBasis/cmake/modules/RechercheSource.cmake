#------------------------------------------------------------------------------
#  Création de macro pour la détection des sources
#------------------------------------------------------------------------------
function(recherche_source ARG)
  unset(A_COMPILER)

  # Parsing du fichier des chemins
  set(file Zcmake/src_path${ARG}.txt)
  file(READ "${file}" CONTENTS)
  string(REPLACE "\n" ";" CONTENTS  ${CONTENTS})

  foreach(ITEM IN ITEMS ${CONTENTS})
    # Détermination du mode (exhaustif ou recursif)
    if( ${ITEM} MATCHES "MODE_RECHERCHE" )
      if( ${ITEM} MATCHES "RECURSIF" )
	set(MODE "RECURSIF")
      else( ${ITEM} MATCHES "EXHAUSTIF" ) 
	set(MODE "EXHAUSTIF")
      endif()
    endif()
    
    # Suppression des commentaires et du mode.
    # Il ne doit rester en principe que les chemins et des blancs...
    if(     (NOT (${ITEM} MATCHES "\#")) 
	AND (NOT (${ITEM} MATCHES "MODE_RECHERCHE"))
	)
      list(APPEND DIR_SOURCES ${ITEM})
    endif()
  endforeach()

  # Recherche des sources dans les répertoires spécifiés. 
  unset(ITEM)
  if (${MODE} MATCHES "RECURSIF" )
    foreach(ITEM IN ITEMS ${DIR_SOURCES})
      file(GLOB_RECURSE FIC_SOURCE 
	${CMAKE_CURRENT_SOURCE_DIR}/${ITEM}/*.[chfF] 
	${CMAKE_CURRENT_SOURCE_DIR}/${ITEM}/*.[fF]90
	${CMAKE_CURRENT_SOURCE_DIR}/${ITEM}/*.cpp
	)
      list(APPEND A_COMPILER ${FIC_SOURCE})
    endforeach()
  else()
    foreach(ITEM IN ITEMS ${DIR_SOURCES})
      file(GLOB FIC_SOURCE 
	${CMAKE_CURRENT_SOURCE_DIR}/${ITEM}/*.[chfF] 
	${CMAKE_CURRENT_SOURCE_DIR}/${ITEM}/*.[fF]90
	${CMAKE_CURRENT_SOURCE_DIR}/${ITEM}/*.cpp
	)
      list(APPEND A_COMPILER ${FIC_SOURCE})
    endforeach()
  endif()

  list(REMOVE_DUPLICATES A_COMPILER)

  # Parsing du fichier des exceptions
  set(file Zcmake/src_ignore${ARG}.txt)
  file(READ "${file}" CONTENTS)
  string(REPLACE "\n" ";" CONTENTS  ${CONTENTS})

  foreach(ITEM IN ITEMS ${CONTENTS})
    # Suppression des commentaires 
    # Il ne doit rester en principe que des fichiers
    if((NOT (${ITEM} MATCHES "\#")))
      list(REMOVE_ITEM A_COMPILER  ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
    endif()
  endforeach()
  
  set(A_COMPILER ${A_COMPILER} PARENT_SCOPE)
endfunction(recherche_source)


