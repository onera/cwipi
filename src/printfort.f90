!-----------------------------------------------------------------------------
! This file is part of the CWIPI library. 
!
! Copyright (C) 2011  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

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
  write(ifile, 1000, advance='no') chloc(1:taille)
  !
  return
  
1000 format(a)
  !
end subroutine printfort

