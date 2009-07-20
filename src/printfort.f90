subroutine printfort (chaine, taille )

!***********************************************************************
!  fonction  : 
!  ---------
!      impression d'une chaine de caracteres (issue d'une fonction c)
!
!-----------------------------------------------------------------------
  use cwipi
  implicit none
  
  character, dimension(*) :: chaine
  integer   ::    taille
  character(len = 16384) :: chloc
  character(len = 64) ::  nom
  integer       ii


  taille = min(taille, 16384 - 1)
  !
  do ii = 1, taille
     chloc(ii:ii) = chaine(ii)
  enddo
  !
  write(ifile, 1000) chloc(1:taille)
  !
  return
  
1000 format(a, $)
  !
end subroutine printfort

