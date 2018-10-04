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
 

!  mpirun -n 1 ./fortran_surf_TriaPi_PiPj : -n 1 ./fortran_surf_TriaPi_PiPj
!  mpirun -n 1 tests/fortran_surf_TriaPi_PiPj : -n 1 tests/fortran_surf_TriaPi_PiPj



module variablesCommunes
 !logical :: visu=.true.
  logical :: visu=.false.
  integer :: commWorld,rankWorld,sizeWorld
  integer :: commLocal,rankLocal,sizeLocal
  integer :: compOrder
  integer :: meshOrder
  
  real(8), pointer :: vand(:,:)
end module variablesCommunes


!> numerotation des types de cellules
!> compatible vtk
module spaceCellTypes
  integer, parameter :: none         =  0
  integer, parameter :: line         =  3  !> bar2      (VTK_LINE                )
  integer, parameter :: triangle     =  5  !> tria3     (VTK_TRIANGLE            )
  integer, parameter :: quad         =  9  !> quad4     (VTK_QUAD                )
  integer, parameter :: tetra        = 10  !> tetra4    (VTK_TETRA               )
  integer, parameter :: hexahedron   = 12  !> hexa8     (VTK_HEXAHEDRON          )
  integer, parameter :: wedge        = 13  !> penta6    (VTK_WEDGE               ) 06 nodes pentahedron
  integer, parameter :: pyramid      = 14  !> pyramid5  (VTK_PYRAMID             )  
  integer, parameter :: line2        = 21  !> bar3      (VTK_QUADRATIC_EDGE      )
  integer, parameter :: triangle2    = 22  !> tria6     (VTK_QUADRATIC_TRIANGLE  )
  integer, parameter :: quad2        = 23  !> quad8     (VTK_QUADRATIC_QUAD      )
  integer, parameter :: tetra2       = 24  !> tetra10   (VTK_QUADRATIC_TETRA     )
  integer, parameter :: hexahedron2  = 25  !> hexa20    (VTK_QUADRATIC_HEXAHEDRON) 
  integer, parameter :: wedge2       = 26  !> penta15   (VTK_QUADRATIC_WEDGE     ) 15 nodes pentahedron
  integer, parameter :: pyramid2     = 27  !> pyramid13 (VTK_QUADRATIC_PYRAMID   ) 13 nodes pyramides

  integer, parameter :: hexahedron3  =101
  integer, parameter :: wedge3       =102
  integer, parameter :: pyramid3     =103
  integer, parameter :: tetra3       =104
  integer, parameter :: quad3        =105
  integer, parameter :: triangle3    =106
  integer, parameter :: line3        =107 
  
  integer, parameter :: hexahedron4  =201
  integer, parameter :: wedge4       =202
  integer, parameter :: pyramid4     =203
  integer, parameter :: tetra4       =204
  integer, parameter :: quad4        =205
  integer, parameter :: triangle4    =206
  integer, parameter :: line4        =207 
  
  integer, parameter :: hexahedron5  =301
  integer, parameter :: wedge5       =302
  integer, parameter :: pyramid5     =303
  integer, parameter :: tetra5       =304
  integer, parameter :: quad5        =305
  integer, parameter :: triangle5    =306
  integer, parameter :: line5        =307 
  
  contains
  
  function cellTypeChar(iCell) result(char)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>
    implicit none
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer       :: iCell
    character(32) :: char
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>
    select case(iCell)
    case( 12)    ; char="hexahedron"
    case( 25)    ; char="hexahedron2"
    case(101)    ; char="hexahedron3"
    case(201)    ; char="hexahedron4"
    case(301)    ; char="hexahedron5"
    case( 13)    ; char="wedge"
    case( 26)    ; char="wedge2"
    case(102)    ; char="wedge3"
    case(202)    ; char="wedge4"
    case(302)    ; char="wedge5"
    case( 14)    ; char="pyramid"
    case( 27)    ; char="pyramid2"
    case(103)    ; char="pyramid3"
    case(203)    ; char="pyramid4"
    case(303)    ; char="pyramid5"
    case( 10)    ; char="tetra"
    case( 24)    ; char="tetra2"
    case(104)    ; char="tetra3"
    case(204)    ; char="tetra4"
    case(304)    ; char="tetra5"
    case(  9)    ; char="quad"
    case( 23)    ; char="quad2"
    case(105)    ; char="quad3"
    case(205)    ; char="quad4"
    case(305)    ; char="quad5"
    case(  5)    ; char="triangle"
    case( 22)    ; char="triangle2"
    case(106)    ; char="triangle3"
    case(206)    ; char="triangle4"
    case(306)    ; char="triangle5"
    case(  3)    ; char="line"
    case( 21)    ; char="line2"
    case(107)    ; char="line3"
    case(207)    ; char="line4"
    case(307)    ; char="line5"
    case default ; char="unknown"
    end select
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end function cellTypeChar
    
end module spaceCellTypes

module spaceMessages

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use iso_fortran_env
  use mpi
  use variablesCommunes
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  

  implicit none
  
contains

  subroutine msg0(msg)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*)                :: msg
    !>
    integer                     :: iRank,iErr
    real(8), allocatable        :: dTab0(:)
    integer, allocatable        :: iTab0(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(iTab0(0:sizeWorld-1))
    
    call mpi_gather(               &
    &    rankWorld, 1, mpi_integer,&
    &    iTab0(0) , 1, mpi_integer,&
    &    0                        ,&
    &    commWorld                ,&
    &    iErr                      )
    
    if( rankWorld==0 )then
     if( sizeWorld>1 )write(*,'()')
      do iRank=0,sizeWorld-1
       !write(*,'("Trace: ",a,t130,"@rkw",i3)')trim(msg),iTab0(iRank)
        write(*,'(a,t130,"@rkw",i3)')trim(msg),iTab0(iRank)
      enddo
    endif
    
    deallocate(iTab0)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine msg0
  
  subroutine msg1(buffer)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*)                   :: buffer
    !>
    integer                        :: length
    integer                        :: iRank,iErr
    character(len=:), allocatable  :: cTab0(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    length=len(buffer)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    allocate( character(len=length) :: cTab0(1:sizeWorld) )
    
    call mpi_gather(                      &
    &    buffer   , length, mpi_character,&
    &    cTab0(1) , length, mpi_character,&
    &    0                               ,&
    &    commWorld                       ,&
    &    iErr                             )
    
    if( rankWorld==0 )then
      !if( sizeWorld>1 )write(*,'()')
      do iRank=1,sizeWorld
        if( .not.len(trim(cTab0(iRank)))==0 )then   !> on retire le cas d'une chaine vide (utilitse par plotInit => micro)
          write(*,'(a)')trim(cTab0(iRank))
        endif
      enddo
    endif
    
    deallocate(cTab0)
    
    call mpi_barrier(commWorld,iErr)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef msg1
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine msg1
  
  subroutine msg2(buffer)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*)                   :: buffer
    !>
    integer                        :: length
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !length=len(buffer)
   !print '(">>> msg2 length=",i6,t130,"@rkw",i3)',length,rankWorld
    
    if( rankWorld==0 )then
     !if( sizeWorld>1 )write(*,'()')
      write(*,'(a)')trim(buffer)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !print '("<<< msg2",t130,"@rkw",i3)',rankWorld
    
    return
  end subroutine msg2
  
  subroutine stopAlert(msg)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*) :: msg
    !>
    integer      :: iErr,iErrCode
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(output_unit ,'(/"ALERT: ",a,"  rkw",i3," called mpi_abort")')trim(msg),rankWorld
    call mpi_abort(commWorld,iErr,iErrCode)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine stopAlert

end module spaceMessages


module additionnal_Functions

contains
    
  subroutine setT3MeshBasis_P1(uv,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Numerotation des sommets
    !>   3
    !>   1 2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! delcaration des variables passees en argument
    real(8), intent(in)    :: uv(:,:)
    real(8), intent(inout) :: ai(:,:)
    !>
    integer                :: iNod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do iNod=1,size(uv,2)
      ai(1,iNod)=1d0-uv(1,iNod)-uv(2,iNod)
      ai(2,iNod)=    uv(1,iNod)
      ai(3,iNod)=               uv(2,iNod)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine setT3MeshBasis_P1
  
  subroutine setT3MeshBasis_P2(uv,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Numerotation des sommets
    !>   3
    !>   6 5
    !>   1 4 2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! delcaration des variables passees en argument
    real(8), intent(in)    :: uv(:,:)
    real(8), intent(inout) :: ai(:,:)
    !>
    integer                :: iNod
    real(8)                :: u,v,w
    real(8)                :: u2,v2,w2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do iNod=1,size(uv,2)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      u=uv(1,iNod)
      v=uv(2,iNod)
      w=1d0-u-v
      u2=2d0*u ; v2=2d0*v ; w2=2d0*w
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ai(1,iNod)=w*(-1d0+w2)    !> (i,j,k)=(0,0,2)
      ai(2,iNod)=u*(-1d0+u2)    !> (i,j,k)=(2,0,0)
      ai(3,iNod)=v*(-1d0+v2)    !> (i,j,k)=(0,2,0)
      ai(4,iNod)=u2*w2          !> (i,j,k)=(1,0,1)
      ai(5,iNod)=u2*v2          !> (i,j,k)=(1,1,0)
      ai(6,iNod)=v2*w2          !> (i,j,k)=(0,1,1)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine setT3MeshBasis_P2
    
end module additionnal_Functions

subroutine  userInterpolation                        ( &
  &           entitiesDim                             ,&
  &           order                                   ,&
  &           nLocalVertex                            ,&
  &           nLocalElement                           ,&
  &           nLocalPolhyedra                         ,&
  &           nDistantPoint                           ,&
  &                                                    &
  &           localCoordinates                        ,&
  &                                                    &
  &           localConnectivityIndex                  ,&
  &           localConnectivity                       ,&
  &           localPolyFaceIndex                      ,&
  &           localPolyCellToFaceConnec               ,&
  &           localPolyFaceConnecIdx                  ,&
  &           localPolyFaceConnec                     ,&
  &                                                    &
  &           disPtsCoordinates                       ,&
  &           disPtsLocation                          ,&
  &           disPtsDistance                          ,&
  &           disPtsBaryCoordIdx                      ,&
  &           distantPointsBarycentricCoordinates     ,&
  &           dist_uvw                                ,& !> new cwipi
  &                                                    &
  &           stride                                  ,& ! =ker(calc)
  &           solverType                              ,&
  &           localField                              ,& !   mySolu
  &           distantField                             ) ! linkSolu
  !---
  use iso_c_binding, only: c_loc,c_f_pointer
  use cwipi
  use baseSimplex2D, only: setT3BasisEqui,setT3MeshIJK
  
  use variablesCommunes
  use additionnal_Functions
  use spaceMessages
  !---
  implicit none
  !---
  integer :: entitiesDim
  integer :: order
  integer :: nLocalVertex
  integer :: nLocalElement
  integer :: nLocalPolhyedra
  integer :: nDistantPoint
  real(8) :: localCoordinates                        (*)
  integer :: localConnectivityIndex                  (*)
  integer :: localConnectivity                       (*)
  integer :: localPolyFaceIndex                      (*)
  integer :: localPolyCellToFaceConnec               (*)
  integer :: localPolyFaceConnecIdx                  (*)
  integer :: localPolyFaceConnec                     (*)
  real(8) :: disPtsCoordinates                       (*)
  integer :: disPtsLocation                          (*)
  real(4) :: disPtsDistance                          (*)
  integer :: disPtsBaryCoordIdx                      (*)
  real(8) :: distantPointsBarycentricCoordinates     (*)
  real(8) :: dist_uvw                                (*)
  integer :: stride
  integer :: solverType
  real(8) ::   localField                            (*)
  real(8) :: distantField                            (*)
  !>
  integer          :: i,j,k,iRank,iErr
  integer          :: iNod,nNod,iMod,nMod
  integer          :: iCell
  real(8), pointer :: uv(:,:)!,a(:),b(:)!c(:)
  integer          :: iDistantPoint
  integer          :: iVert
  integer, pointer :: ij(:,:)
  real(8), pointer :: lagrangeMesh(:,:)
  integer          :: nod(1:16)
  real(8)          :: delta(1:3)
  character(2048)  :: buffer
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(buffer,'(a,">>> userInterpolation (callback)")')char(10) ; call msg2(buffer)  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> dist_uvw(:) -> uv(1:2,:) (sans dupliquer le bloc mémoire)
  !> la commmande c_f_pointer a besoin du use iso_c_binding, only: c_loc,c_f_pointer
  
  call c_f_pointer(cptr=c_loc(dist_uvw), fptr=uv, shape=[2,nDistantPoint])  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Calcul de distantField
  
  nMod=(meshOrder+1)*(meshOrder+2)/2
  allocate(lagrangeMesh(1:nMod,1:nDistantPoint))
  
  select case(meshOrder)
  case(1) ; call setT3MeshBasis_P1(uv=uv,ai=lagrangeMesh) !> base Triangle Geometrique P2
 !case(2) ; call setT3MeshBasis_P2(uv=uv,ai=lagrangeMesh) !> base Triangle Geometrique P2
  case default
    
    allocate(ij(1:2,1:nMod))
    call setT3MeshIJK(meshOrder=meshOrder,ij=ij)
    call setT3BasisEqui(ord=meshOrder,ijk=ij,uvw=uv,ai=lagrangeMesh)
    deallocate(ij)
  end select
  
  j=0
  do iDistantPoint=1,nDistantPoint
    iCell=disPtsLocation(iDistantPoint)
    nMod=localConnectivityIndex(iCell+1)-localConnectivityIndex(iCell)
    nod(1:nMod)=localConnectivity(localConnectivityIndex(iCell)+1:localConnectivityIndex(iCell+1))
    
    distantField(j+1:j+stride)=0d0
    do iMod=1,nMod
      i=stride*(nod(iMod)-1)
      distantField(j+1:j+stride)=distantField(j+1:j+stride)+lagrangeMesh(iMod,iDistantPoint)*localField(i+1:i+stride)
      
      !> Controling localField
      k=3     *(nod(iMod)-1)
      delta(1:3)=localField(i+1:i+3)-localCoordinates(k+1:k+3)
      if( .not.(delta(1)==0d0.and.delta(2)==0d0.and.delta(3)==0d0) )then
        print '("iCell=",i6," Delta=",3(e22.15,1x),t130,"@rkw",i3)',iCell,delta(1:3),rankWorld
      endif
    enddo
    
    j=j+stride
    k=k+3
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Visu
  if( nDistantPoint== 1)then
    iDistantPoint=1
    j=stride*(iDistantPoint-1)
    iCell=disPtsLocation(iDistantPoint)
    nMod=localConnectivityIndex(iCell+1)-localConnectivityIndex(iCell)
    nod(1:nMod)=localConnectivity(localConnectivityIndex(iCell)+1:localConnectivityIndex(iCell+1))    
    
    i=stride*(iDistantPoint-1)
    j=     3*(iDistantPoint-1)
    delta(1:3)=disPtsCoordinates(j+1:j+3)-distantField(i+1:i+3)
    
    if( meshOrder==1 )then
      
      write(buffer,'(                                                   &
      &                                                              a, &
      &              3x,"iCell=",i6," nod: ",3(i3,1x),t130,"@rkw",i3,a, &
      &              6x,"xyz1=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz2=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz3=",3(e22.15,1x)                        ,a, &
      &                                                              a, &
      &              6x,"uv  =",2(e22.15,1x)                        ,a, &
      &              6x,"disPtsCoordinates=",3(e22.15,1x)           ,a, &
      &              6x,"distantField=     ",3(e22.15,1x)           ,a, &
      &              6x,"delta=            ",3(e22.15,1x)               &
      &                                                              )')&
      &                                                  char(10),&
      &  iCell,nod(1:nMod),rankWorld                    ,char(10),&
      &  localCoordinates(3*(nod(1)-1)+1:3*(nod(1)-1)+3),char(10),&
      &  localCoordinates(3*(nod(2)-1)+1:3*(nod(2)-1)+3),char(10),&
      &  localCoordinates(3*(nod(3)-1)+1:3*(nod(3)-1)+3),char(10),&
      &                                                  char(10),&
      &  uv(1:2,iDistantPoint)                          ,char(10),&
      & disPtsCoordinates(j+1:j+3)                      ,char(10),&
      & distantField     (i+1:i+3)                      ,char(10),&
      & disPtsCoordinates(j+1:j+3)-distantField(i+1:i+3)
      
    elseif( meshOrder==2 )then
      
      write(buffer,'(                                                   &
      &                                                              a, &
      &              3x,"iCell=",i6," nod: ",6(i6,1x),t130,"@rkw",i3,a, &
      &              6x,"xyz1=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz2=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz3=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz4=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz5=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz6=",3(e22.15,1x)                        ,a, &
      &                                                              a, &
      &              6x,"uv  =",2(e22.15,1x)                        ,a, &
      &              6x,"disPtsCoordinates=",3(e22.15,1x)           ,a, &
      &              6x,"distantField=     ",3(e22.15,1x)           ,a, &
      &              6x,"delta=            ",3(e22.15,1x)               &
      &                                                              )')&
      &                                                  char(10),&
      &  iCell,nod(1:nMod),rankWorld                    ,char(10),&
      &  localCoordinates(3*(nod(1)-1)+1:3*(nod(1)-1)+3),char(10),&
      &  localCoordinates(3*(nod(2)-1)+1:3*(nod(2)-1)+3),char(10),&
      &  localCoordinates(3*(nod(3)-1)+1:3*(nod(3)-1)+3),char(10),&
      &  localCoordinates(3*(nod(4)-1)+1:3*(nod(4)-1)+3),char(10),&
      &  localCoordinates(3*(nod(5)-1)+1:3*(nod(5)-1)+3),char(10),&
      &  localCoordinates(3*(nod(6)-1)+1:3*(nod(6)-1)+3),char(10),&
      &                                                  char(10),&
      &  uv(1:2,iDistantPoint)                          ,char(10),&
      & disPtsCoordinates(j+1:j+3)                      ,char(10),&
      & distantField     (i+1:i+3)                      ,char(10),&
      & disPtsCoordinates(j+1:j+3)-distantField(i+1:i+3)
      
    elseif( meshOrder==3 )then
      
      write(buffer,'(                                                    &
      &                                                               a, &
      &              3x,"iCell=",i6," nod: ",10(i6,1x),t130,"@rkw",i3,a, &
      &              6x,"xyz01=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz02=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz03=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz04=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz05=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz06=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz07=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz08=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz09=",3(e22.15,1x)                        ,a, &
      &              6x,"xyz10=",3(e22.15,1x)                        ,a, &
      &                                                               a, &
      &              6x,"uv  =",2(e22.15,1x)                         ,a, &
      &              6x,"disPtsCoordinates=",3(e22.15,1x)            ,a, &
      &              6x,"distantField=     ",3(e22.15,1x)            ,a, &
      &              6x,"delta=            ",3(e22.15,1x)                &
      &                                                               )')&
      &                                                    char(10),&
      &  iCell,nod(1:nMod),rankWorld                      ,char(10),&
      &  localCoordinates(3*(nod(01)-1)+1:3*(nod(01)-1)+3),char(10),&
      &  localCoordinates(3*(nod(02)-1)+1:3*(nod(02)-1)+3),char(10),&
      &  localCoordinates(3*(nod(03)-1)+1:3*(nod(03)-1)+3),char(10),&
      &  localCoordinates(3*(nod(04)-1)+1:3*(nod(04)-1)+3),char(10),&
      &  localCoordinates(3*(nod(05)-1)+1:3*(nod(05)-1)+3),char(10),&
      &  localCoordinates(3*(nod(06)-1)+1:3*(nod(06)-1)+3),char(10),&
      &  localCoordinates(3*(nod(07)-1)+1:3*(nod(07)-1)+3),char(10),&
      &  localCoordinates(3*(nod(08)-1)+1:3*(nod(08)-1)+3),char(10),&
      &  localCoordinates(3*(nod(09)-1)+1:3*(nod(09)-1)+3),char(10),&
      &  localCoordinates(3*(nod(10)-1)+1:3*(nod(10)-1)+3),char(10),&
      &                                                   char(10),&
      &  uv(1:2,iDistantPoint)                           ,char(10),&
      & disPtsCoordinates(j+1:j+3)                       ,char(10),&
      & distantField     (i+1:i+3)                       ,char(10),&
      & disPtsCoordinates(j+1:j+3)-distantField(i+1:i+3)
      
    endif
    call msg1(buffer)    
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !call mpi_barrier(commWorld,iErr)
  !if( rankWorld==0 )print'(3x,"<<< userInterpolation")'
  !call mpi_barrier(commWorld,iErr)
  
  write(buffer,'(a,"<<< userInterpolation (callback)")')char(10) ; call msg2(buffer)  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine userInterpolation


program fortran_surf_TriaPi_PiPj
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use iso_fortran_env
  use iso_c_binding, only: c_loc,c_f_pointer
  
  use mpi
  use cwipi
  
  use variablesCommunes
  use additionnal_Functions
  use spaceCellTypes
  use spaceMessages
  use baseSimplex1D, only: nodes1D,setL2BasisEqui,setQ4MeshIJK
  use baseSimplex2D, only: nodes2D,setT3BasisEqui,setT3MeshIJK
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  interface
    subroutine  userInterpolation                      ( &
    &           entitiesDim                             ,&
    &           order                                   ,&
    &           nLocalVertex                            ,&
    &           nLocalElement                           ,&
    &           nLocalPolhyedra                         ,&
    &           nDistantPoint                           ,&
    &                                                    &
    &           localCoordinates                        ,&
    &                                                    &
    &           localConnectivityIndex                  ,&
    &           localConnectivity                       ,&
    &           localPolyFaceIndex                      ,&
    &           localPolyCellToFaceConnec               ,&
    &           localPolyFaceConnecIdx                  ,&
    &           localPolyFaceConnec                     ,&
    &                                                    &
    &           disPtsCoordinates                       ,&
    &           disPtsLocation                          ,&
    &           disPtsDistance                          ,&
    &           disPtsBaryCoordIdx                      ,&
    &           distantPointsBarycentricCoordinates     ,&
    &           dist_uvw                                ,&
    &                                                    &
    &           stride                                  ,&  ! =ker(calc)
    &           solverType                              ,&
    &           localField                              ,&  !   mySolu
    &           distantField                             )  ! linkSolu
    !---
    integer :: entitiesDim
    integer :: order
    integer :: nLocalVertex
    integer :: nLocalElement
    integer :: nLocalPolhyedra
    integer :: nDistantPoint
    real(8) :: localCoordinates                        (*)
    integer :: localConnectivityIndex                  (*)
    integer :: localConnectivity                       (*)
    integer :: localPolyFaceIndex                      (*)
    integer :: localPolyCellToFaceConnec               (*)
    integer :: localPolyFaceConnecIdx                  (*)
    integer :: localPolyFaceConnec                     (*)
    real(8) :: disPtsCoordinates                       (*)
    integer :: disPtsLocation                          (*)
    real(4) :: disPtsDistance                          (*)
    integer :: disPtsBaryCoordIdx                      (*)
    real(8) :: distantPointsBarycentricCoordinates     (*)
    real(8) :: dist_uvw(*)
    integer :: stride
    integer :: solverType
    real(8) ::   localField                            (*)
    real(8) :: distantField                            (*)
    end subroutine  userInterpolation
  end interface
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  include 'libmeshb7.ins'
  
  character(5)       :: codeName,codeCoupledName  
  character(128)     :: meshName
  
  integer            :: iVert,nVert
  real(8), pointer   :: vertx(:)
  integer, pointer   :: vertM(:)
  integer, pointer   :: cells(:),cellsIdx(:),mark(:),types(:)
  integer, pointer   :: nod(:)  !> ensemble des noeuds pour iCell : nod=>cells(cellsIdx(iCell)+1:cellsIdx(iCell+1))

  
  character(256)     :: key
  integer, parameter :: iFile=100
  logical            :: test
  
  integer(8)         :: InpMsh
  integer(8)         :: ad0
  integer            :: ad1,ad2
  integer            :: dim,res,ver
  integer            :: iCell0
  
  
  integer            :: i,j,k,l,iu,iv
  integer            :: iMod,nMod
  integer            :: iNod,nNod
  integer            :: iCell,nCell
  integer            :: nQ4,nT3
  
  integer, pointer :: ij(:,:),ijCwipi(:)
  real(8), pointer :: lagrangeMeshQ4(:,:)
  real(8), pointer :: lagrangeMeshT3(:,:)
  real(8), pointer :: lagrangeL2    (:,:)
  
  real(8)          :: tol
  real(8)          :: xyz(1:3)
  integer          :: linkVertSize
  real(8), pointer :: linkVert(:)
  integer          :: notLocatedPoints
  
  integer          :: stride
  real(8), pointer ::   myValues(:)
  real(8), pointer :: linkValues(:)
  
  real(8), pointer :: u (  :)  !> pour les quad (segments tensorisés)
  real(8), pointer :: uv(:,:)  !> pour les triangles
  integer          :: numberOfUnlocatedPoints,numberOfUnlocatedPointsGlob
  
  integer          :: iRank,iErr
  
  integer          :: iVertMax
  real(8)          :: delta,deltaMin,deltaMax,sumDelta
  character(1024)  :: buffer
  real(8)          :: t0,t1
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_init(iErr)
  commWorld=mpi_comm_world
  
  call mpi_comm_rank(commWorld, rankWorld, iErr)
  call mpi_comm_size(commWorld, sizeWorld, iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !if( rankWorld==0) print '(/"START: fortran_surf_TriaPi_PiPj")'
  write(buffer,'("START: fortran_surf_TriaPi_PiPj")') ; call msg2(trim(buffer))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do meshOrder=1,1
    
    call cpu_time(t0)
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Initialisation de l'interface de couplage
    select case(rankWorld)
    case(0)
       codeName        = "code1"
       codeCoupledName = "code2"
    case(1)
       codeName        = "code2"
       codeCoupledName = "code1"
    end select
    
    !> permet de recuperer un commLocal
    call cwipi_init_f(           &
    &    globalComm=commWorld   ,&
    &    appliName=codeName     ,&
    &    appliComm=commLocal     )
    
    call mpi_comm_rank(commLocal,rankLocal,iErr)
    call mpi_comm_size(commLocal,sizeLocal,iErr)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    select case(meshOrder)
    case(1) ; tol=5d-1
    case(2) ; tol=1d-1
    case(3) ; tol=1d-1
    end select
    
    write(buffer,'("")')                                                         ; call msg2(trim(buffer))
    write(buffer,'("Create coupling tol=",e22.15,t130,"@rkw",i3)')tol,rankWorld  ; call msg1(trim(buffer))
      
    call cwipi_create_coupling_f(                  &
    &    couplingName="testPiPj"                  ,&
    &    couplingType=cwipi_cpl_parallel_with_part,&
    &    cplAppli=codeCoupledName                 ,&
    &    entitiesDim=2                            ,& !> Nature du couplage surfacique
    &    tolerance=tol                            ,& !> Tolerance geometrique 1d-1 par defaut
    &    meshT=cwipi_static_mesh                  ,&
    &    solvert=cwipi_solver_cell_vertex         ,&
    &    outputfreq=-1                            ,& !> Frequence du post-traitement
    &    outputfmt="Ensight Gold"                 ,&
    &    outputfmtopt="binary"                     )  
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    select case(rankWorld)
    case(0) ; compOrder=1 !10 !07 !07
    case(1) ; compOrder=1 !10 !07 !10
    end select
    
    write(buffer,'("")')                                                   ; call msg2(trim(buffer))
    write(buffer,'(3x,"compOrder=",i3,t130,"@rkw",i3)')compOrder,rankWorld ; call msg1(trim(buffer))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Création des maillages avec gmsh
    
    write(buffer,'("")')                                                   ; call msg2(trim(buffer))
    write(buffer,'(3x,"meshOrder=",i3,t130,"@rkw",i3)')meshOrder,rankWorld ; call msg1(trim(buffer))
        
    select case(rankWorld)
    case(0)
      write(key,'(6x,"gmsh spaceBasis/tests/Mesh2D/sphere01.geo -2 -format msh -order ",i1," > sphere01.log")')meshOrder      
      call msg1(trim(key))
      call execute_command_line (key, exitstat=iErr)
      write(key,'(6x,"~/Maillages/mshTomesh spaceBasis/tests/Mesh2D/sphere01.msh >> sphere01.log")')
      call msg1(trim(key))
      call execute_command_line (key, exitstat=iErr)
    case(1)
      write(key,'(6x,"gmsh spaceBasis/tests/Mesh2D/sphere02.geo -2 -format msh -order ",i1," > sphere02.log")')meshOrder
      call msg1(trim(key))
      call execute_command_line (key, exitstat=iErr)
      write(key,'(6x,"~/Maillages/mshTomesh spaceBasis/tests/Mesh2D/sphere02.msh >> sphere02.log")')
      call msg1(trim(key))
      call execute_command_line (key, exitstat=iErr)
    end select
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Reading Geometric Mesh with INRIA libMesh7
    !> attention avec gmfGetBlock certains entiers sont integer(8)
    
    write(buffer,'("")')                                               ; call msg2(trim(buffer))
    write(buffer,'("Reading Geometric Mesh",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))
    select case(rankWorld)
    case(0) ; meshName="./spaceBasis/tests/Mesh2D/sphere01.mesh"
    case(1) ; meshName="./spaceBasis/tests/Mesh2D/sphere02.mesh"
    end select
    
    !>>>>>>>
    !> Opening File
    InpMsh = gmfOpenMesh(trim(meshName),GmfRead,ver,dim)
    !<<<<<<<
    
    !>>>>>>>
    !> Mesh Sizes
    nVert = gmfstatkwd(InpMsh, GmfVertices)
    
    nCell=0 ; dim=0
    select case(meshOrder)
    case(1) ; nQ4=gmfstatkwd(InpMsh, GmfQuadrilaterals  )
    case(2) ; nQ4=gmfstatkwd(InpMsh, GmfQuadrilateralsQ2)
    case(3) ; nQ4=gmfstatkwd(InpMsh, GmfQuadrilateralsQ3)
    case(4) ; nQ4=gmfstatkwd(InpMsh, GmfQuadrilateralsQ4)
    case default ; call stopAlert("meshOrder>4")
    end select
    nCell=nCell+nQ4
    dim=dim    +nQ4*(meshOrder+1)*(meshOrder+1)
    
    select case(meshOrder)
    case(1) ; nT3=gmfstatkwd(InpMsh, GmfTriangles  )
    case(2) ; nT3=gmfstatkwd(InpMsh, GmfTrianglesP2)
    case(3) ; nT3=gmfstatkwd(InpMsh, GmfTrianglesP3)
    case(4) ; nT3=gmfstatkwd(InpMsh, GmfTrianglesP4)
    case default ; call stopAlert("meshOrder>4")
    end select
    nCell=nCell+nT3
    dim  =dim  +nT3*(meshOrder+1)*(meshOrder+2)/2
    
    write(buffer,'(3x,"meshOrder=",i2,"  nVert=",i6," nQ4=",i6," nT3=",i6,t130,"@rkw",i3)')meshOrder,nVert,nQ4,nT3,rankWorld ; call msg1(trim(buffer))  
    !<<<<<<<
    
    !>>>>>>>
    !> Reading Vertices
    
    !> allocation vertx,vertM
    allocate( vertx(1:3*nVert),vertM(1:nVert) )
    
    !> read block
    ad0=int(1,kind=8)
    ad1=0
    ad2=3*(nVert-1)
    !write(buffer,'(3x,"Vert: ad1:ad2=",i6,":",i6,t130,"@rkw",i3)')ad1,ad2,rankWorld ; call msg1(trim(buffer))  
    res = gmfGetBlock(                         &
    &     InpMsh                              ,&
    &     GmfVertices                         ,&
    &     ad0                                 ,&
    &     int(nVert,kind=8)                   ,&
    &     0, %val(0), %val(0)                 ,&
    &     GmfDouble,vertx(ad1+1),vertx(ad2+1) ,&
    &     GmfDouble,vertx(ad1+2),vertx(ad2+2) ,&
    &     GmfDouble,vertx(ad1+3),vertx(ad2+3) ,&
    &                                          &
    &     GmfInt   ,vertM(1)    ,vertM(nVert)  )  
    
  !  if( rankWorld==0 )then
  !    do iVert=1,10
  !      xyz(1:3)=vertx(3*(iVert-1)+1:3*iVert)
  !      print '("iVert=",i6," xyz=",3(e22.15,1x)," mark=",i12)',iVert,xyz(1:3),vertM(iVert)
  !    enddo
  !  endif
    !<<<<<<<
    
    !>>>>>>>
    !> Setting cellsIdx and Reading cells,cellsRef
    
    !> allocation cellsIdx,cells,mark
    allocate(cellsIdx(1:nCell+1),cells(1:dim),mark(1:nCell),types(1:nCell))
    
    !> Initialization (iCell0,cellsIdx(1))
    iCell0=0 ; cellsIdx(1)=0
    
    !> Reading Quadrilaterals
    ReadingQuadrilaterals: if( .not.nQ4==0 )then
      nCell=nQ4                                         ! <=
      nNod=(meshOrder+1)*(meshOrder+1)                  ! <=
      !> nodesIdx
      do i=1,nCell
        cellsIdx(iCell0+1+i)=cellsIdx(iCell0+1)+nNod*i
      enddo
      
      !> types
      select case(meshOrder)
      case(1) ; types(iCell0+1:iCell0+nCell)=quad       ! <=
      case(2) ; types(iCell0+1:iCell0+nCell)=quad2      ! <=
      case(3) ; types(iCell0+1:iCell0+nCell)=quad3      ! <=
      case(4) ; types(iCell0+1:iCell0+nCell)=quad4      ! <=
      case default ; call stopAlert("meshOrder>4")
      end select
      
      !> Adr  Block
      ad0=1_8                    !> debut
      ad1=cellsIdx(iCell0    +1) !> iCell0+1
      ad2=cellsIdx(iCell0+nCell) !> iCell0+nCell
      !print '("Quad: ad1:ad2=",i6,":",i6)',ad1,ad2
      
      !> Read Block
      select case(meshOrder)
      case(1)
        res=gmfGetBlock(                               &
        &   InpMsh                                    ,&
        &   GmfQuadrilaterals                         ,&  ! <=
        &   ad0                                       ,&
        &   int(nCell,kind=8)                         ,&
        &   0, %val(0), %val(0)                       ,&
        &   GmfInt, cells(ad1+ 1) , cells(ad2+ 1)     ,&
        &   GmfInt, cells(ad1+ 2) , cells(ad2+ 2)     ,&
        &   GmfInt, cells(ad1+ 3) , cells(ad2+ 3)     ,&
        &   GmfInt, cells(ad1+ 4) , cells(ad2+ 4)     ,&
        &                                              &
        &   GmfInt, mark(iCell0+1), mark(iCell0+nCell) )
      case(2)
        res=GmfGetBlock(                               &
        &   InpMsh                                    ,&
        &   GmfQuadrilateralsQ2                       ,&  ! <=
        &   ad0                                       ,&
        &   int(nCell,kind=8)                         ,&
        &   0, %val(0), %val(0)                       ,&
        &   GmfInt, cells(ad1+ 1) , cells(ad2+ 1)     ,&
        &   GmfInt, cells(ad1+ 2) , cells(ad2+ 2)     ,&
        &   GmfInt, cells(ad1+ 3) , cells(ad2+ 3)     ,&
        &   GmfInt, cells(ad1+ 4) , cells(ad2+ 4)     ,&
        &   GmfInt, cells(ad1+ 5) , cells(ad2+ 5)     ,&
        &   GmfInt, cells(ad1+ 6) , cells(ad2+ 6)     ,&
        &   GmfInt, cells(ad1+ 7) , cells(ad2+ 7)     ,&
        &   GmfInt, cells(ad1+ 8) , cells(ad2+ 8)     ,&
        &   GmfInt, cells(ad1+ 9) , cells(ad2+ 9)     ,&
        &                                              &
        &   GmfInt, mark(iCell0+1), mark(iCell0+nCell) )
      
      case default ; call stopAlert("meshOrder>4")
      end select
      
      iCell0=iCell0+nCell
    endif ReadingQuadrilaterals
    
    !> Reading Triangles
    ReadingTriangle: if( .not.nT3==0 )then
      nCell=nT3                                         ! <=
      nNod=(meshOrder+1)*(meshOrder+2)/2                ! <=
      !> nodesIdx
      do i=1,nCell
        cellsIdx(iCell0+1+i)=cellsIdx(iCell0+1)+nNod*i
      enddo
      
      !> types
      select case(meshOrder)                            
      case(1) ; types(iCell0+1:iCell0+nCell)=triangle   ! <=
      case(2) ; types(iCell0+1:iCell0+nCell)=triangle2  ! <=
      case(3) ; types(iCell0+1:iCell0+nCell)=triangle3  ! <=
      case(4) ; types(iCell0+1:iCell0+nCell)=triangle4  ! <=
      case default ; call stopAlert("meshOrder>4")
      end select
      
      !> Adr  Block
      ad0=1_8                    !> debut
      ad1=cellsIdx(iCell0    +1) !> iCell0+1
      ad2=cellsIdx(iCell0+nCell) !> iCell0+nCell
      !write(buffer,'(3x,"Tria: ad1:ad2=",i6,":",i6,t130,"@rkw",i3)')ad1,ad2,rankWorld ; call msg1(trim(buffer))
      
      !> Read Block
      select case(meshOrder)
      case(1)    
        res=GmfGetBlock(                               &
        &   InpMsh                                    ,&
        &   GmfTriangles                              ,&  ! <=
        &   ad0                                       ,&
        &   int(nCell,kind=8)                         ,&
        &   0, %val(0), %val(0)                       ,&
        &   GmfInt, cells(ad1+ 1) , cells(ad2+ 1)     ,&
        &   GmfInt, cells(ad1+ 2) , cells(ad2+ 2)     ,&
        &   GmfInt, cells(ad1+ 3) , cells(ad2+ 3)     ,&
        &                                              &
        &   GmfInt, mark(iCell0+1), mark(iCell0+nCell) )
      case(2)
        res=GmfGetBlock(                               &
        &   InpMsh                                    ,&
        &   GmfTrianglesP2                            ,&  ! <=
        &   ad0                                       ,&
        &   int(nCell,kind=8)                         ,&
        &   0, %val(0), %val(0)                       ,&
        &   GmfInt, cells(ad1+ 1) , cells(ad2+ 1)     ,&
        &   GmfInt, cells(ad1+ 2) , cells(ad2+ 2)     ,&
        &   GmfInt, cells(ad1+ 3) , cells(ad2+ 3)     ,&
        &   GmfInt, cells(ad1+ 4) , cells(ad2+ 4)     ,&
        &   GmfInt, cells(ad1+ 5) , cells(ad2+ 5)     ,&
        &   GmfInt, cells(ad1+ 6) , cells(ad2+ 6)     ,&
        &                                              &
        &   GmfInt, mark(iCell0+1), mark(iCell0+nCell) )
      case(3)
        res=GmfGetBlock(                               &
        &   InpMsh                                    ,&
        &   GmfTrianglesP3                            ,&  ! <=
        &   ad0                                       ,&
        &   int(nCell,kind=8)                         ,&
        &   0, %val(0), %val(0)                       ,&
        &   GmfInt, cells(ad1+ 1) , cells(ad2+ 1)     ,&
        &   GmfInt, cells(ad1+ 2) , cells(ad2+ 2)     ,&
        &   GmfInt, cells(ad1+ 3) , cells(ad2+ 3)     ,&
        &   GmfInt, cells(ad1+ 4) , cells(ad2+ 4)     ,&
        &   GmfInt, cells(ad1+ 5) , cells(ad2+ 5)     ,&
        &   GmfInt, cells(ad1+ 6) , cells(ad2+ 6)     ,&
        &   GmfInt, cells(ad1+ 7) , cells(ad2+ 7)     ,&
        &   GmfInt, cells(ad1+ 8) , cells(ad2+ 8)     ,&
        &   GmfInt, cells(ad1+ 9) , cells(ad2+ 9)     ,&
        &   GmfInt, cells(ad1+10) , cells(ad2+10)     ,&
        &                                              &
        &   GmfInt, mark(iCell0+1), mark(iCell0+nCell) )
        
      case default ; call stopAlert("reading Triangles meshOrder>3")
      end select
      
      iCell0=iCell0+nCell
    endif ReadingTriangle
    
    !if( rankWorld==0 )then
    !  do iCell=1,10
    !    nod=>cells(cellsIdx(iCell)+1:cellsIdx(iCell+1))
    !    print '("iCell=",i6," nod=",*(i6,1x))',iCell,nod(:)
    !  enddo
    !endif
    !<<<<<<<
    
    !>>>>>>>
    !> Closing File
    res = gmfclosemesh(InpMsh)
    !<<<<<<<  
    
    write(buffer,'(                                &
    &                                           a, &
    &              3x,"mesh=",a,t130,"@rkw",i3 ,a, &
    &              6x,"meshOrder=",i1          ,a, &
    &              6x,"nVert=",i6              ,a, &
    &              6x,"nQ4  =",i6              ,a, &
    &              6x,"nT3  =",i6                  &
    &                                           )')&
    &                          char(10),&
    & trim(meshName),rankWorld,char(10),&
    & meshOrder               ,char(10),&
    & nVert                   ,char(10),&
    & nQ4                     ,char(10),&
    & nT3
    
    call msg1(trim(buffer))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> On se couple sur les sphères maillées
    
    write(buffer,'("")')                                              ; call msg2(trim(buffer))
    write(buffer,'("Sending Mesh to Cwipi",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))
    
    nCell=nQ4+nT3
    
    !> Transmission des maillages à cwipi
    call cwipi_ho_define_mesh_f(  & !> NEW Cwipi
    &   couplingName="testPiPj"  ,&
    &   nVertex     =nVert       ,&
    &   nElts       =nQ4+nT3     ,&
    &   order       =meshOrder   ,&  
    &   coords      =vertx       ,&
    &   connecIndex =cellsIdx    ,&
    &   connec      =cells        )  
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Definition du Quad Qi (ijk)
    if( .not.nQ4==0 )then
      nMod=(meshOrder+1)*(meshOrder+1)
      allocate(ij(1:2,nMod))
      
      call setQ4MeshIJK(meshOrder=meshOrder,ij=ij)
      call c_f_pointer(cptr=c_loc(ij), fptr=ijCwipi, shape=[2*nMod])  
      
      call cwipi_ho_ordering_from_IJK_set_f( & !> NEW Cwipi
      &   couplingName ="testPiPj"          ,&
      &   tElt         = CWIPI_FACE_QUADHO  ,&
      &   nNodes       = nMod               ,&
      &   IJK          = ijCwipi             )
      
      deallocate(ij)
    endif
    
    !> Definition du Triangle géométrique Pi (ijk)
    if( .not.nT3==0 )then
      nMod=(meshOrder+1)*(meshOrder+2)/2
      allocate(ij(1:2,nMod))
      
      call setT3MeshIJK(meshOrder=meshOrder,ij=ij)
      call c_f_pointer(cptr=c_loc(ij), fptr=ijCwipi, shape=[2*nMod])  
      
      call cwipi_ho_ordering_from_IJK_set_f( & !> NEW Cwipi
      &   couplingName ="testPiPj"          ,&
      &   tElt         = CWIPI_FACE_TRIAHO  ,&
      &   nNodes       = nMod               ,&
      &   IJK          = ijCwipi             )
      
      deallocate(ij)
    endif
    
    write(buffer,'("")')                                                             ; call msg2(trim(buffer))
    write(buffer,'("Geometric HO Cell ij P",i1,t130,"@rkw",i3)')meshOrder,rankWorld  ; call msg1(trim(buffer))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Initialisation of myValues(:)
    
    write(buffer,'("")')                                                  ; call msg2(trim(buffer))
    write(buffer,'("Init myValues (x,y,z,rkw)",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))  
    
    stride=4 !> x,y,z,real(rankWorld,kind=8)
    allocate(myValues(1:nVert*stride))
    i=0 ; j=0
    do iVert=1,nVert
      myValues(j+1:j+stride)=[vertx(i+1),vertx(i+2),vertx(i+3),real(rankWorld,kind=8)]
      i=i+3
      j=j+4
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Points de couplage linkVertSize,linkVert(:)
    
    write(buffer,'("")')                                        ; call msg2(trim(buffer))
    write(buffer,'("Set linkVert(:)",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))
    
#if 0==0
    
    if( .not.nQ4==0 )then
      !> calcul lagrangeMesh
      
      nMod=(meshOrder+1)                                 !> Edge meshOrder (entree)
      nNod=(compOrder+1)                                 !> Edge compOrder (sortie)
      call nodes1D(ord=compOrder,uvw=u,display=.false.)  !> ordre du calcul
      
      allocate(ij(1:1,1:nMod)) 
      do iMod=1,meshOrder+1              !> 01 02  ... meshOrder meshOrder+1 (rangement regulier)
        ij(1,iMod)=(iMod-1)
      enddo
      
      allocate(lagrangeL2(1:nMod,1:nNod))
      call setL2BasisEqui(ord=meshOrder,ijk=ij,uvw=u,ai=lagrangeL2)
      
      block
      real(8) :: toto(11)
      print '(/"u=",*(e12.5,1x))',u(:)
      print '(/"lagrangeL2")'
      do iMod=1,nMod
        toto(1:nMod)=lagrangeL2(iMod,1:nNod)
        print '(3x,"ai(",i2,",1:nNod)=",*(e12.5,1x) )',iMod,toto(1:nMod)
      enddo
      print '("")'
      end block
      
      deallocate(ij)
      
      allocate(ij(1:2,1:nMod*nNod))      
      call setQ4MeshIJK(meshOrder=meshOrder,ij=ij)
      
      allocate(lagrangeMeshQ4(1:nMod*nMod,1:nNod*nNod))
      
      iMod=0
      do j=1,nMod ; do i=1,nMod
          iMod=iMod+1 ! print '(3x,"iMod=",i3,3x,"i,j=",2(i3,1x))',iMod,ij(1:2,iMod)
          iu=ij(1,iMod)+1
          iv=ij(2,iMod)+1
          !>
          print '(3x,"iMod=",i3,3x,"(i,j)=",2(i3,1x))',iMod,iu,iv          
          iNod=0
          do l=1,nNod ; do k=1,nNod
            iNod=iNod+1
            print '(6x,"iNod=",i3,3x,"(u,v)=",2(e12.5,1x),3x,"aiL2=",2(e12.5,1x))',iNod,u(k),u(l),lagrangeL2(iu,k),lagrangeL2(iv,l)
            
            lagrangeMeshQ4(iMod,iNod)= lagrangeL2(iu,k) &
            &                         *lagrangeL2(iv,l)
          enddo ; enddo
      enddo ; enddo
      
      deallocate(u)
      deallocate(lagrangeL2)
      
      
      block
      real(8) :: toto(11)
      print '(/"lagrangeMeshQ4")'
      do iMod=1,nMod*nMod
        toto(1:nMod*nMod)=lagrangeMeshQ4(iMod,1:nMod*nMod)
        print '(3x,"ai(",i2,",1:nNod)=",*(e12.5,1x) )',iMod,toto(1:nMod*nMod)
      enddo
      print '("")'
      end block
      
      
               
    endif
    
    if( .not.nT3==0 )then
      !> calcul lagrangeMesh
      nMod=(meshOrder+1)*(meshOrder+2)/2   !> Triangle meshOrder (entree)
      call nodes2D(ord=compOrder,uvw=uv,display=.false.) !> ordre du calcul  
      nNod=size(uv,2)
      allocate(lagrangeMeshT3(1:nMod,1:nNod))
      select case(meshOrder)
      case(01) ; call setT3MeshBasis_P1(uv=uv,ai=lagrangeMeshT3)
     !case(02) ; call setT3MeshBasis_P2(uv=uv,ai=lagrangeMeshT3)
      case default
        allocate(ij(1:2,1:nMod))
        call setT3MeshIJK(meshOrder=meshOrder,ij=ij)
        call setT3BasisEqui(ord=meshOrder,ijk=ij,uvw=uv,ai=lagrangeMeshT3)
        deallocate(ij)
      end select
      deallocate(uv)
    endif
    
    
    !> calcul linkVert
    linkVertSize=0
    if( .not.nQ4==0 )then
      nNod=(compOrder+1)*(compOrder+1)                   !> Quad meshOrder (entree)
      linkVertSize=linkVertSize+nQ4*nNod                 !> nombre total de point de couplages
    endif
    if( .not.nT3==0 )then
      nNod=(compOrder+1)*(compOrder+2)/2                 !> Triangle meshOrder (entree)
      linkVertSize=linkVertSize+nT3*nNod                 !> nombre total de point de couplages
    endif    
    allocate(linkVert(1:3*linkVertSize))                 !> 3 coordonnées par point de couplage
    
    
    !> Initialization
    iCell0=0 ; j=0
    
    !> Quadrilaterals
    if( .not.nQ4==0 )then
      nMod=(meshOrder+1)*(meshOrder+1)                   !> Quad meshOrder (entree)
      nNod=(compOrder+1)*(compOrder+1)                   !> Quad compOrder (sortie)
      do i=1,nQ4
        iCell=iCell0+i
        nod=>cells(cellsIdx(iCell)+1:cellsIdx(iCell+1))  !> nMod=size(node)
        
        if( iCell==1 )then
          do iMod=1,nMod
            k=3*(nod(iMod)-1)
            xyz(1:3)=vertx(k+1:k+3)
            print '(3x"nod(",i8,")=",3(e22.15,1x))',nod(iMod),xyz(1:3)
          enddo
        endif
        
        do iNod=1,nNod
          !> linkVert
          linkVert(j+1:j+3)=0d0
          do iMod=1,nMod
            k=3*(nod(iMod)-1)
            linkVert(j+1:j+3)=linkVert(j+1:j+3)+lagrangeMeshQ4(iMod,iNod)*vertx(k+1:k+3)            
          enddo
          
          if( iCell==1 )then
            xyz(1:3)=linkVert(j+1:j+3)
            print '(3x"linkVert(",i3,")=",3(e22.15,1x))',iNod,xyz(1:3)
          endif
          
          j=j+3      
        enddo
      enddo
      !>
      deallocate(lagrangeMeshQ4)
      iCell0=iCell0+nQ4
    endif
    
    !> Triangles
    if( .not.nT3==0 )then
      nMod=(meshOrder+1)*(meshOrder+2)/2                 !> Triangle meshOrder (entree)
      nNod=(compOrder+1)*(compOrder+2)/2                 !> Triangle compOrder (sortie)
      do i=1,nT3
        iCell=iCell0+i
        nod=>cells(cellsIdx(iCell)+1:cellsIdx(iCell+1))        
        do iNod=1,nNod
          !> linkVert
          linkVert(j+1:j+3)=0d0
          do iMod=1,nMod
            k=3*(nod(iMod)-1)
            linkVert(j+1:j+3)=linkVert(j+1:j+3)+lagrangeMeshT3(iMod,iNod)*vertx(k+1:k+3)
          enddo
          j=j+3      
        enddo
      enddo
      !>
      deallocate(lagrangeMeshT3)
      iCell0=iCell0+nT3
    endif
    
    
    !if( rankWorld==0 )then
    !  print '(/"linkVert")'
    !  j=0
    !  do iVert=1,11 !linkVertSize
    !    xyz(1:3)=linkVert(j+1:j+3)
    !    print '(3x"linkVert(",i2,")=",3(e22.15,1x))',iVert,xyz(1:3)
    !    j=j+3
    !  enddo
    !endif
    
    
    
#else
    
    linkVertSize=1                       !> nombre total de point de couplages
    allocate(linkVert(1:3*linkVertSize)) !> 3 coordonnées par point de couplage
    
    select case(rankWorld)
    case(0) ; linkVert(1:3)=[-0.919633675189875E+00,-0.250163898379974E-01, 0.392131352367653E+00]  
    case(1) ; linkVert(1:3)=[-0.848807623678876E+00,-0.213274154008489E-01, 0.424121411593377E+00]
    end select
    
    write(buffer,'(6x,"linkVert(1:3)=    ",3(e22.15,1x),t130,"@rkw",i3)')linkVert(1:3),rankWorld     ; call msg1(trim(buffer))
    
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(buffer,'("")')                                                                              ; call msg2(trim(buffer))
    write(buffer,'("Allocation of linkValues size=",i6,t130,"@rkw",i3)')stride*linkVertSize,rankWorld ; call msg1(trim(buffer))
    
    allocate(linkValues(1:stride*linkVertSize))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transmission à cwipi des coordonnees de couplage
    
    write(buffer,'("")')                                                         ; call msg2(trim(buffer))
    write(buffer,'("Transmission de linkVert a cwipi",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))  
      
    call cwipi_set_points_to_locate_f( &
    &    couplingName="testPiPj"      ,&
    &    nPts  =linkVertSize          ,&
    &    coords=linkVert               )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Localisation par cwipi des coordonnees de couplage
    
    write(buffer,'("")')                                              ; call msg2(trim(buffer))
    write(buffer,'("Localisation by Cwipi",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))  
    
    call cwipi_locate_f(couplingName="testPiPj")
    
    call cwipi_get_n_not_located_pts_f(couplingName="testPiPj", nNotLocatedPoints=numberOfUnlocatedPoints)
    
    call mpi_allreduce(numberOfUnlocatedPoints,numberOfUnlocatedPointsGlob,1,mpi_integer,mpi_sum,commWorld,iErr)
    
    if( .not.numberOfUnlocatedPointsGlob==0 )then
      write(buffer,'(3x,"nNotLocatedPoints=",i10,t130,"@rkw",i3)')numberOfUnlocatedPoints,rankWorld ; call msg1(trim(buffer))  
      write(buffer,'("fortran_surf_TriaPi_PiPj stop line: ",i6,t130,"@rkw",i3)')__LINE__+1,rankWorld
      call stopAlert(trim(buffer))  
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(buffer,'("")')                                                 ; call msg2(trim(buffer))
    write(buffer,'("Echange cwipi_exchange_f",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))  
    
    call cwipi_exchange_f(                       &
    &    couplingName=          "testPiPj"      ,&
    &    exchangeName="exch1_"//"testPiPj"      ,&
    &    meshOrder=meshOrder                    ,&  !> NEW cwipi
    &    exchangeDim=stride                     ,&  !> scalar
    &    ptHoInterpolationFct=userInterpolation ,&  !> utilisation de la procedure plug
    &                                            &
    &    sendingFieldName="mySolu"              ,&  !> solution calculee localement
    &    sendingField=myValues                  ,&
    &                                            &
    &    receivingFieldName="linkSolu"          ,&
    &    receivingField=linkValues              ,&  !> solution de raccord
    &                                            &
    &    nStep=1                                ,&  !> pas utilisee juste pour visu cwipi
    &    timeValue=0d0                          ,&  !> pas utilisee juste pour visu cwipi
    &    nNotLocatedPoints=notLocatedPoints     ,&
    &    status=iErr                             )  
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(buffer,'("")')                                        ; call msg2(trim(buffer))
    write(buffer,'("Delete coupling",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))  
    
    call cwipi_delete_coupling_f("testPiPj")
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(buffer,'("")')                                           ; call msg2(trim(buffer))
    write(buffer,'("Controling Results",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))
    
    iVertMax=1
    sumDelta= 0d0
    deltaMax=-1d50
    deltaMin= 1d50
    j=0
    k=0
    do iVert=1,linkVertSize
      delta=norm2(linkVert(j+1:j+3)-linkValues(k+1:k+3))
      sumDelta=sumDelta+delta
      if( deltaMax<delta )then
        iVertMax=iVert
        deltaMax=delta
      endif
      if( delta<deltaMin )deltaMin=delta
      j=j+3
      k=k+4
    enddo
    
    sumDelta=sumDelta/real(linkVertSize,kind=8)  
    
    j=(iVertMax-1)*3
    k=(iVertMax-1)*stride
    write(buffer,'(                                  &
    &                                             a, &
    &              3x,"Control",t130,"@rkw",i3   ,a, &
    &              6x,"meshOrder=",i1            ,a, &
    &              6x,"linkVertSize=",i6         ,a, &
    &                                             a, &
    &              6x,"deltaMin  =",e22.15       ,a, &
    &              6x,"sumDelta  =",e22.15       ,a, &
    &              6x,"deltaMax  =",e22.15       ,a, &
    &                                             a, &
    &              6x,"linkVert  =",3(e22.15,1x) ,a, &
    &              6x,"linkValues=",3(e22.15,1x) ,a, &
    &              6x,"Delta     =",3(e22.15,1x)     &
    &                                             )')&
    &                          char(10),&
    & rankWorld               ,char(10),&
    & meshOrder               ,char(10),&
    & linkVertSize            ,char(10),&
    &                          char(10),&
    & deltaMin                ,char(10),&
    & sumDelta                ,char(10),&
    & deltaMax                ,char(10),&
    &                          char(10),&
    & linkVert  (j+1:j+3)     ,char(10),&
    & linkValues(k+1:k+3)     ,char(10),&
    & linkVert  (j+1:j+3)-linkValues(k+1:k+3)
    
    call msg1(trim(buffer))      
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(vertx,vertM)
    deallocate(cellsIdx,cells,mark,types)
    deallocate(myValues,linkVert,linkValues)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call cwipi_finalize_f()
    
    call cpu_time(t1)
    write(buffer,'("")')                                                       ; call msg2(trim(buffer))
    write(buffer,'(3x,"cpu_time=",f12.5," s",t130,"@rkw",i3)')t1-t0,rankWorld  ; call msg1(trim(buffer))
    
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_finalize(iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
end program fortran_surf_TriaPi_PiPj
