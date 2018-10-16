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
    
  subroutine mshToMesh(mshName)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> http://geuz.org/gmsh/doc/texinfo/#MSH-ASCII-file-format
    !>
    !> 01:   2-node 1st order line.
    !> 02:   3-node 1st order triangle.
    !> 03:   4-node 1st order quadrangle.
    !> 04:   4-node 1st order tetrahedron.
    !> 05:   8-node 1st order hexahedron.
    !> 06:   6-node 1st order prism.
    !> 07:   5-node 1st order pyramid.
    !> 08:   3-node 2nd order line         (2 nodes associated with the vertices and 1 with the edge).
    !> 09:   6-node 2nd order triangle     (3 nodes associated with the vertices and 3 with the edges).
    !> 10:   9-node 2nd order quadrangle   (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
    !> 11:  10-node 2nd order tetrahedron  (4 nodes associated with the vertices and 6 with the edges).
    !> 12:  27-node 2nd order hexahedron   (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
    !> 13:  18-node 2nd order prism        (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
    !> 14:  14-node 2nd order pyramid      (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
    !> 15:   1-node point.
    !> 16:   8-node 2nd order quadrangle   (4 nodes associated with the vertices and 4 with the edges).
    !> 17:  20-node 2nd order hexahedron   (8 nodes associated with the vertices and 12 with the edges).
    !> 18:  15-node 2nd order prism        (6 nodes associated with the vertices and 9 with the edges).
    !> 19:  13-node 2nd order pyramid      (5 nodes associated with the vertices and 8 with the edges).
    !> 20:   9-node 3rd order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
    !> 21:  10-node 3rd order triangle     (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
    !> 22:  12-node 4th order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
    !> 23:  15-node 4th order triangle     (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
    !> 24:  15-node 5th order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
    !> 25:  21-node 5th order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
    !> 26:   4-node 3rd order edge (2 nodes associated with the vertices, 2 internal to the edge)
    !> 27:   5-node 4th order edge (2 nodes associated with the vertices, 3 internal to the edge)
    !> 28:   6-node 5th order edge (2 nodes associated with the vertices, 4 internal to the edge)
    !> 29:  20-node 3rd order tetrahedron  (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
    !> 30:  35-node 4th order tetrahedron  (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
    !> 31:  56-node 5th  order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
    !> 92:  64-node 3rd  order hexahedron  (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)
    !> 93: 125-node 4th order hexahedron   (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume)
    !>  
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*)    , intent(in)  :: mshName
    character(len=:), allocatable :: meshName
    integer                       :: iArg,nArg,length,l
    
    integer              :: i,j,cellType,nbParam,entite,numPhysicalNames,dimEntity,mark
    integer              :: iVert,nVert,iCell,nCell,iNod
    integer              :: nL2  ,nT3  ,nQ4  ,nT4  ,nP5  ,nW5  ,nH6
    integer              :: nL2P2,nT3P2,nQ4P2,nT4P2,nP5P2,nW5P2,nH6P2
    integer              :: nL2P3,nT3P3,nQ4P3,nT4P3,nP5P3,nW5P3,nH6P3
    integer              :: nL2P4,nT3P4,nQ4P4,nT4P4,nP5P4,nW5P4,nH6P4
    integer              :: iVert1,nVert1
    real(8), allocatable :: vert(:,:),vert1(:,:)
    logical, allocatable :: connected(:)
    integer, allocatable :: idx(:)
    character(256)       :: ligne
    character( 32)       :: version
    integer, allocatable :: hexas(:,:),hexas2(:,:)
    integer, allocatable :: wedge(:,:),wedge2(:,:)
    integer, allocatable :: pyras(:,:),pyras2(:,:)
    integer, allocatable :: tetra(:,:),tetra2(:,:)
    integer, allocatable :: quadr(:,:),quadr2(:,:),quadr3(:,:),quadr4(:,:)
    integer, allocatable :: trian(:,:),trian2(:,:),trian3(:,:),trian4(:,:)
    integer, allocatable :: edges(:,:),edges2(:,:),edges3(:,:),edges4(:,:)
    integer              :: tab(1:27)
    integer, parameter   :: hexaQ2Gmsh2Inria (1:27)=[01,02,03,04,05,06,07,08,09,12,14,10,11,13,15,16,17,19,20,18,21,23,22,26,24,25,27] ! <=
    integer, parameter   :: hexaQ2Inria2Medit(1:27)=[01,02,03,04,05,06,07,08,09,10,11,12,17,18,19,20,13,14,15,16,21,24,23,25,26,22,27]
   !integer, parameter   :: hexaQ2test (1:27)=[01,02,03,04,05,06,07,08,09,12,14,10,17,19,20,18,11,13,15,16,21,26,22,24,25,23,27]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    length=index(mshName,'.',.true.)-1 ! print '("length=",i3)',length
    
    allocate(character(len=length+5) :: meshName)
    write(meshName,'(a,".mesh")')mshName(1:length)
    
    !print '("input:  ",a)',trim( mshName)
    !print '("output: ",a)',trim(meshName)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nH6=0 ; nH6P2=0
    nW5=0 ; nW5P2=0
    nP5=0 ; nP5P2=0
    nT4=0 ; nT4P2=0
    nQ4=0 ; nQ4P2=0 ; nQ4P3=0 ; nQ4P4=0
    nT3=0 ; nT3P2=0 ; nT3P3=0 ; nT3P4=0
    nL2=0 ; nL2P2=0 ; nL2P3=0 ; nL2P4=0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !print '(">>> Allocation")'
    open(unit=10,file=mshName,action='read')
    lecture0: do
      read(10,*)ligne
      
      if( trim(ligne)=="$MeshFormat" )then
        
        read(10,'(a8)')version ; !print '("msh version: ",a)',version
        !> End Item
        read(10,*)ligne
        if( .not.trim(ligne)=="$EndMeshFormat" )then
          print '("Problème durant lecture: ",a)',trim(ligne)
        endif
        
      elseif( trim(ligne)=="$PhysicalNames" )then
        
        read(10,*)numPhysicalNames ; !print '(3x,"numPhysicalNames=",i3)',numPhysicalNames
        do i=1,numPhysicalNames    
          read(10,*)dimEntity,mark,ligne ; !print '(6x,"mark=",i3," -> ",a)',mark,trim(ligne)
        enddo
        
        read(10,*)ligne
        if( .not.trim(ligne)=="$EndPhysicalNames" )then
          print '("Problème durant lecture: ",a)',trim(ligne)
        endif
        
      elseif( trim(ligne)=="$Nodes" )then
        read(10,*)nVert 
        do iVert=1,nVert
          read(10,*)!i,x,y,z
        enddo
        !> End Item
        read(10,*)ligne
        if( .not.trim(ligne)=="$EndNodes" )then
          print '("Problème durant lecture: ",a)',trim(ligne)
        endif
        !print '("nVert:",i10)',nVert
        
      elseif( trim(ligne)=="$Elements" )then
        
        !> $Elements
        !> number−of−elements
        !> elm−number elm−type number−of−tags <tags> node−number−
        !> list ...
        !> $EndElements
         
        read(10,*)nCell ; !print '("nCell:",i10)',nCell
        do iCell=1,nCell
          read(10,'(a256)')ligne ! print '(a)',ligne
          
          read(ligne,*)i,cellType,NbParam
          select case(cellType)
          case( 1) ; nL2  =nL2  +1 !> Edge
          case( 2) ; nT3  =nT3  +1 !> Triangle
          case( 3) ; nQ4  =nQ4  +1 !> Quad
          case( 4) ; nT4  =nT4  +1 !> Tetra
          case( 5) ; nH6  =nH6  +1 !> Hexa
          case( 6) ; nW5  =nW5  +1 !> Wedge
          case( 7) ; nP5  =nP5  +1 !> Pyramid
          !> P2
          case( 8) ; nL2P2=nL2P2+1 !> Edge2
          case( 9) ; nT3P2=nT3P2+1 !> Triangle2
          case(10) ; nQ4P2=nQ4P2+1 !> Quad2
          case(11) ; nT4P2=nT4P2+1 !> Tetra2
          case(12) ; nH6P2=nH6P2+1 !> Hexa2
          case(13) ; nW5P2=nW5P2+1 !> Wedge2
          case(14) ; nP5P2=nP5P2+1 !> Pyramid2
          !>
          case(15) !> node point on s'en fiche
          !> P3
          case(21) ; nT3P3=nT3P3+1 !> Triangle3
          case(26) ; nL2P3=nL2P3+1 !> Edge3
          case(36) ; nQ4P3=nQ4P3+1 !> Quad3
          !> P4
          case(23) ; nT3P4=nT3P4+1 !> Triangle4 -> Triangle
          case(27) ; nL2P4=nL2P4+1 !> Edge4     -> Edge
          case(37) ; nQ4P4=nQ4P4+1 !> Quad4     -> Quad
         !case(23) ; nT3P2=nT3P2+1 !> Triangle4 -> Triangle2
         !case(27) ; nL2P2=nL2P2+1 !> Edge4     -> Edge2
         !case(37) ; nQ4P2=nQ4P2+1 !> Quad4     -> Quad2
          case default ; print '("iCell=",i6," kind of cell not implemented cellType=",i10)',iCell,cellType ; stop
          end select
        enddo
        !> End Item
        read(10,*)ligne
        if( .not.trim(ligne)=="$EndElements" )then
          print '("Problème durant lecture: ",a)',trim(ligne)
        else
          exit lecture0
        endif
        
      endif
    enddo lecture0
    
    close(10)
    
    
!     if( .not.nVert==0 )print '(3x,"nVert=",i10)',nVert
!     
!     if( .not.nH6  ==0 )print '(3x,"nH6  =",i10)',nH6
!     if( .not.nH6P2==0 )print '(3x,"nH6P2=",i10)',nH6P2
!     
!     if( .not.nW5  ==0 )print '(3x,"nW5  =",i10)',nW5
!     if( .not.nW5P2==0 )print '(3x,"nW5P2=",i10)',nW5P2
!     
!     if( .not.nP5  ==0 )print '(3x,"nP5  =",i10)',nP5
!     if( .not.nP5P2==0 )print '(3x,"nP5P2=",i10)',nP5P2
!     
!     if( .not.nT4  ==0 )print '(3x,"nT4  =",i10)',nT4
!     if( .not.nT4P2==0 )print '(3x,"nT4P2=",i10)',nT4P2
!     
!     if( .not.nQ4  ==0 )print '(3x,"nQ4  =",i10)',nQ4
!     if( .not.nQ4P2==0 )print '(3x,"nQ4P2=",i10)',nQ4P2
!     if( .not.nQ4P3==0 )print '(3x,"nQ4P3=",i10)',nQ4P3
!     if( .not.nQ4P4==0 )print '(3x,"nQ4P4=",i10)',nQ4P4
!     
!     if( .not.nT3  ==0 )print '(3x,"nT3  =",i10)',nT3
!     if( .not.nT3P2==0 )print '(3x,"nT3P2=",i10)',nT3P2
!     if( .not.nT3P3==0 )print '(3x,"nT3P3=",i10)',nT3P3
!     if( .not.nT3P4==0 )print '(3x,"nT3P4=",i10)',nT3P4
!     
!     if( .not.nL2  ==0 )print '(3x,"nL2  =",i10)',nL2
!     if( .not.nL2P2==0 )print '(3x,"nL2P2=",i10)',nL2P2
!     if( .not.nL2P3==0 )print '(3x,"nL2P3=",i10)',nL2P3
!     if( .not.nL2P4==0 )print '(3x,"nL2P4=",i10)',nL2P4
    
    allocate(vert(1:3,1:nVert))
    
    if( .not.nH6  ==0 )allocate( hexas (1:09,1:nH6  ) ) !>  8 + 1
    if( .not.nH6P2==0 )allocate( hexas2(1:28,1:nH6P2) ) !> 27 + 1
    
    if( .not.nW5  ==0 )allocate( wedge (1:07,1:nW5  ) ) !>  6 + 1
    if( .not.nW5P2==0 )allocate( wedge2(1:19,1:nW5P2) ) !> 18 + 1
    
    if( .not.nP5  ==0 )allocate( pyras (1:06,1:nP5  ) ) !>  5 + 1
    if( .not.nP5P2==0 )allocate( pyras2(1:15,1:nP5P2) ) !> 14 + 1
    
    if( .not.nT4  ==0 )allocate( tetra (1:05,1:nT4  ) ) !>  4 + 1
    if( .not.nT4P2==0 )allocate( tetra2(1:11,1:nT4P2) ) !> 10 + 1
    
    if( .not.nQ4  ==0 )allocate( quadr (1:05,1:nQ4  ) ) !>  4 + 1
    if( .not.nQ4P2==0 )allocate( quadr2(1:10,1:nQ4P2) ) !>  9 + 1
    if( .not.nQ4P3==0 )allocate( quadr3(1:17,1:nQ4P3) ) !> 16 + 1
    if( .not.nQ4P4==0 )allocate( quadr4(1:26,1:nQ4P4) ) !> 25 + 1
    
    if( .not.nT3  ==0 )allocate( trian (1:04,1:nT3  ) ) !>  3 + 1
    if( .not.nT3P2==0 )allocate( trian2(1:07,1:nT3P2) ) !>  6 + 1
    if( .not.nT3P3==0 )allocate( trian3(1:11,1:nT3P3) ) !> 10 + 1
    if( .not.nT3P4==0 )allocate( trian4(1:16,1:nT3P4) ) !> 15 + 1
    
    if( .not.nL2  ==0 )allocate( edges (1:03,1:nL2  ) ) !>  2 + 1
    if( .not.nL2P2==0 )allocate( edges2(1:04,1:nL2P2) ) !>  3 + 1
    if( .not.nL2P3==0 )allocate( edges3(1:05,1:nL2P3) ) !>  4 + 1
    if( .not.nL2P4==0 )allocate( edges4(1:06,1:nL2P4) ) !>  5 + 1
    
    !print '("<<< Allocation")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !print '(">>> Reading")'
    
    nH6=0 ; nH6P2=0
    nW5=0 ; nW5P2=0
    nP5=0 ; nP5P2=0
    nT4=0 ; nT4P2=0
    nQ4=0 ; nQ4P2=0 ; nQ4P3=0 ; nQ4P4=0
    nT3=0 ; nT3P2=0 ; nT3P3=0 ; nT3P4=0
    nL2=0 ; nL2P2=0 ; nL2P3=0 ; nL2P4=0
    
    open(unit=10,file=mshName ,action='read' )
    lecture: do
      read(10,*)ligne
      
      if( trim(ligne)=="$MeshFormat" )then
        
        read(10,'(a8)')version
        read(10,*)ligne
        if( .not.trim(ligne)=="$EndMeshFormat" )then
          print '("Problème durant lecture: ",a)',trim(ligne)
        endif
        
      elseif( trim(ligne)=="$PhysicalNames" )then
        
        read(10,*)numPhysicalNames ! print '(3x,"numPhysicalNames=",i3)',numPhysicalNames
        do i=1,numPhysicalNames    
          read(10,*)dimEntity,mark,ligne ! print '(6x,"mark=",i3," -> ",a)',mark,trim(ligne)
        enddo
        
        read(10,*)ligne
        if( .not.trim(ligne)=="$EndPhysicalNames" )then
          print '("Problème durant lecture: ",a)',trim(ligne)
        endif
        
      elseif( trim(ligne)=="$Nodes" )then
        
        read(10,*)nVert
        do iVert=1,nVert
          read(10,*)i,vert(1:3,i)
        enddo
        
        read(10,*)ligne
        if( .not.trim(ligne)=="$EndNodes" )then
          print '("Problème durant lecture: ",a)',trim(ligne)
        endif
        
      elseif( trim(ligne)=="$Elements" )then
        
        read(10,*)nCell ! print '("nCell=",i10)',nCell
        do iCell=1,nCell
          read(10,'(a256)')ligne ! print '(a)',ligne
          read(ligne,*)i,cellType,NbParam
          select case(cellType)
          case( 1) ; nL2  =nL2  +1 ; read(ligne,*)i,cellType,NbParam,edges (03,nL2  ),entite,edges (1:02,nL2  )
          case( 2) ; nT3  =nT3  +1 ; read(ligne,*)i,cellType,NbParam,trian (04,nT3  ),entite,trian (1:03,nT3  )
          case( 3) ; nQ4  =nQ4  +1 ; read(ligne,*)i,cellType,NbParam,quadr (05,nQ4  ),entite,quadr (1:04,nQ4  ) ! print '(5(i3,1x))',quadr (1:5,nQ4  )
          case( 4) ; nT4  =nT4  +1 ; read(ligne,*)i,cellType,NbParam,tetra (05,nT4  ),entite,tetra (1:04,nT4  )
          case( 5) ; nH6  =nH6  +1 ; read(ligne,*)i,cellType,NbParam,hexas (09,nH6  ),entite,hexas (1:08,nH6  )
          case( 6) ; nW5  =nW5  +1 ; read(ligne,*)i,cellType,NbParam,wedge (07,nW5  ),entite,wedge (1:06,nW5  )
          case( 7) ; nP5  =nP5  +1 ; read(ligne,*)i,cellType,NbParam,pyras (06,nP5  ),entite,pyras (1:05,nP5  )
          !>
          case( 8) ; nL2P2=nL2P2+1 ; read(ligne,*)i,cellType,NbParam,edges2(04,nL2P2),entite,edges2(1:03,nL2P2)
          case( 9) ; nT3P2=nT3P2+1 ; read(ligne,*)i,cellType,NbParam,trian2(07,nT3P2),entite,trian2(1:06,nT3P2)
          case(10) ; nQ4P2=nQ4P2+1 ; read(ligne,*)i,cellType,NbParam,quadr2(10,nQ4P2),entite,quadr2(1:09,nQ4P2)
          case(11) ; nT4P2=nT4P2+1 ; read(ligne,*)i,cellType,NbParam,tetra2(11,nT4P2),entite,tetra2(1:08,nT4P2),tetra2(10,nT4P2),tetra2(9,nT4P2) !> Retournement 9-10
          case(12) ; nH6P2=nH6P2+1 ; read(ligne,*)i,cellType,NbParam,hexas2(28,nH6P2),entite,hexas2(1:27,nH6P2)
          !>
          case(21) ; nT3P3=nT3P3+1 ; read(ligne,*)i,cellType,NbParam,trian3(11,nT3P3),entite,trian3(1:10,nT3P3)
          case(26) ; nL2P3=nL2P3+1 ; read(ligne,*)i,cellType,NbParam,edges3(05,nL2P3),entite,edges3(1:04,nL2P3)
          case(36) ; nQ4P3=nQ4P3+1 ; read(ligne,*)i,cellType,NbParam,quadr3(17,nQ4P3),entite,quadr3(1:16,nQ4P3)
          !>
          case(23) ; nT3P4=nT3P4+1 ; read(ligne,*)i,cellType,NbParam,trian4(16,nT3P4),entite,trian4(1:15,nT3P4)
          case(27) ; nL2P4=nL2P4+1 ; read(ligne,*)i,cellType,NbParam,edges4(06,nL2P4),entite,edges4(1:05,nL2P4)
          case(37) ; nQ4P4=nQ4P4+1 ; read(ligne,*)i,cellType,NbParam,quadr4(26,nQ4P4),entite,quadr4(1:25,nQ4P4)
          !>
          case(15) !> node point on s'en fiche
          case default ; print '("iCell=",i6," kind of cell not implemented cellType=",i10)',iCell,cellType ; stop
          end select
        enddo
        
        read(10,*)ligne
        if( .not.trim(ligne)=="$EndElements" )then
          print '("Problème durant lecture: ",a)',trim(ligne)
        else
          exit lecture
        endif
        
      endif
    enddo lecture
    
    close(10)
    !print '("<<< Reading")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Supression des points non connectés
    
    !print '(">>> Supressing unconnected vertices")'
    allocate(connected(1:nVert)) ; connected(1:nVert)=.false.
    
    do iCell=1,nH6
      do iNod=1,8
        connected(hexas (iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nH6P2
      do iNod=1,27
        connected(hexas2(iNod,iCell))=.true.
      enddo
    enddo
    
    do iCell=1,nW5
      do iNod=1,6
        connected(wedge (iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nW5P2
      do iNod=1,18
        connected(wedge2(iNod,iCell))=.true.
      enddo
    enddo
    
    do iCell=1,nP5
      do iNod=1,5
        connected(pyras (iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nP5P2
      do iNod=1,14
        connected(pyras2(iNod,iCell))=.true.
      enddo
    enddo
    
    do iCell=1,nT4
      do iNod=1,4
        connected(tetra (iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nT4P2
      do iNod=1,10
        connected(tetra2(iNod,iCell))=.true.
      enddo
    enddo
    
    do iCell=1,nQ4
      do iNod=1,4
        connected(quadr (iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nQ4P2
      do iNod=1,9
        connected(quadr2(iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nQ4P3
      do iNod=1,16
        connected(quadr3(iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nQ4P4
      do iNod=1,25
        connected(quadr4(iNod,iCell))=.true.
      enddo
    enddo
    
    do iCell=1,nT3
      do iNod=1,3
        connected(trian (iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nT3P2
      do iNod=1,6
        connected(trian2(iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nT3P3
      do iNod=1,10
        connected(trian3(iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nT3P4
      do iNod=1,15
        connected(trian4(iNod,iCell))=.true.
      enddo
    enddo
    
    do iCell=1,nL2
      do iNod=1,2
        connected(edges (iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nL2P2
      do iNod=1,3
        connected(edges2(iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nL2P3
      do iNod=1,4
        connected(edges3(iNod,iCell))=.true.
      enddo
    enddo
    do iCell=1,nL2P4
      do iNod=1,5
        connected(edges4(iNod,iCell))=.true.
      enddo
    enddo
    
    nVert1=count(connected)
    if(.not.nVert1==nVert )then
      !print '(3x,"nVert=",i10," -> ",i10)',nVert,nVert1
      
      allocate(idx(1:nVert))
      allocate(vert1(1:3,1:nVert1))
      iVert1=0
      do iVert=1,nVert
        if( connected(iVert) )then
          iVert1=iVert1+1
          vert1(1:3,iVert1)=vert(1:3,iVert)
          idx(iVert)=iVert1
         !print '("iVert=",i10,2x,"iVert1=",i10,2x,"xyz=",3(f12.5,1x))',iVert,iVert1,vert1(1:3,iVert1)
        endif
      enddo
      
      nVert=nVert1 ; call move_alloc(from=vert1, to=vert)
      
      do iCell=1,nH6
        do iNod=1,8
          hexas (iNod,iCell)=idx(hexas (iNod,iCell))
        enddo
      enddo
      do iCell=1,nH6P2
        do iNod=1,27
          hexas2(iNod,iCell)=idx(hexas2(iNod,iCell))
        enddo
      enddo
      
      do iCell=1,nW5
        do iNod=1,6
          wedge (iNod,iCell)=idx(wedge (iNod,iCell))
        enddo
      enddo
      do iCell=1,nW5P2
        do iNod=1,18
          wedge2(iNod,iCell)=idx(wedge2(iNod,iCell))
        enddo
      enddo
      
      do iCell=1,nP5
        do iNod=1,5
          pyras (iNod,iCell)=idx(pyras (iNod,iCell))
        enddo
      enddo
      do iCell=1,nP5P2
        do iNod=1,14
          pyras2(iNod,iCell)=idx(pyras2(iNod,iCell))
        enddo
      enddo
      
      do iCell=1,nT4
        do iNod=1,4
          tetra (iNod,iCell)=idx(tetra (iNod,iCell))
        enddo
      enddo
      do iCell=1,nT4P2
        do iNod=1,10
          tetra2(iNod,iCell)=idx(tetra2(iNod,iCell))
        enddo
      enddo
      
      do iCell=1,nQ4
        do iNod=1,4
          quadr (iNod,iCell)=idx(quadr (iNod,iCell))
        enddo
      enddo
      do iCell=1,nQ4P2
        do iNod=1,9
          quadr2(iNod,iCell)=idx(quadr2(iNod,iCell))
        enddo
      enddo
      do iCell=1,nQ4P3
        do iNod=1,16
          quadr3(iNod,iCell)=idx(quadr3(iNod,iCell))
        enddo
      enddo
      do iCell=1,nQ4P4
        do iNod=1,25
          quadr4(iNod,iCell)=idx(quadr4(iNod,iCell))
        enddo
      enddo
      
      do iCell=1,nT3
        do iNod=1,3
          trian (iNod,iCell)=idx(trian (iNod,iCell))
        enddo
      enddo
      do iCell=1,nT3P2
        do iNod=1,6
          trian2(iNod,iCell)=idx(trian2(iNod,iCell))
        enddo
      enddo
      do iCell=1,nT3P3
        do iNod=1,10
          trian3(iNod,iCell)=idx(trian3(iNod,iCell))
        enddo
      enddo
      do iCell=1,nT3P4
        do iNod=1,15
          trian4(iNod,iCell)=idx(trian4(iNod,iCell))
        enddo
      enddo
      
      do iCell=1,nL2
        do iNod=1,2
          edges (iNod,iCell)=idx(edges (iNod,iCell))
        enddo
      enddo
      do iCell=1,nL2P2
        do iNod=1,3
          edges2(iNod,iCell)=idx(edges2(iNod,iCell))
        enddo
      enddo
      do iCell=1,nL2P3
        do iNod=1,4
          edges3(iNod,iCell)=idx(edges3(iNod,iCell))
        enddo
      enddo
      do iCell=1,nL2P4
        do iNod=1,5
          edges4(iNod,iCell)=idx(edges4(iNod,iCell))
        enddo
      enddo
      
      deallocate(idx)
      
    endif
    deallocate(connected)
    
    !print '(3x,"x \in [",f12.5,",",f12.5,"]")',minval(vert(1,:)),maxval(vert(1,:))
    !print '(3x,"y \in [",f12.5,",",f12.5,"]")',minval(vert(2,:)),maxval(vert(2,:))
    !print '(3x,"z \in [",f12.5,",",f12.5,"]")',minval(vert(3,:)),maxval(vert(3,:))
    
    !print '("<<< Supressing unconnected vertices")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !print '(">>> Writing")'
    
    open(unit=20,file=meshName,action='write')
    write(20,'("MeshVersionFormatted 2")')
    write(20,'(/"Dimension"/"3"/)')
    
    write(20,'(/"Vertices")')
    write(20,*)nVert
    do iVert=1,nVert
       write(20,'(3(e22.15,1x),1x,i1)')vert(1:3,iVert),0
      !write(20,'(3(f3.0,1x),1x,i2)')vert(1:3,iVert),iVert
    enddo
    
    if( .not.nH6==0 )then
      write(20,'(/"Hexahedra")')
      write(20,'(i10)')nH6
      do iCell=1,nH6
        write(20,'(8(i10,1x),2x,i3)')hexas(:,iCell)
      enddo
    endif
    if( .not.nH6P2==0 )then
      write(20,'(/"HexahedraQ2")')
      write(20,'(i6)')nH6P2
      do iCell=1,nH6P2
        
        !if( iCell==1 )then
        !print '(/"Cell:",i6)',iCell
        !do i=1,27
        !  if( .not. hexaQ2test(i)==hexaQ2Gmsh2Inria(hexaQ2Inria2Medit(i)) )then
        !    stop 'bad conversion'
        !  endif
        !enddo
        !endif
        
        write(20,'(27(i10,1x),2x,i2)')hexas2(hexaQ2Gmsh2Inria(hexaQ2Inria2Medit(1:27)),iCell),hexas2(28,iCell)
      enddo
    endif
    
    
    if( .not.nW5==0 )then
      write(20,'(/"Prisms")')
      write(20,'(i10)')nW5
      do iCell=1,nW5
        write(20,'(6(i10,1x),2x,i3)')wedge(:,iCell)
      enddo
    endif
    
    if( .not.nP5==0 )then
      write(20,'(/"Pyramids")')
      write(20,'(i10)')nP5
      do iCell=1,nP5
        write(20,'(5(i10,1x),2x,i3)')pyras(:,iCell)
      enddo
    endif
    
    if( .not.nT4==0 )then
      write(20,'(/"Tetrahedra")')
      write(20,'(i10)')nT4
      do iCell=1,nT4
        write(20,'(4(i10,1x),2x,i3)')tetra(:,iCell)
      enddo
    endif
    if( .not.nT4P2==0 )then
      write(20,'(/"TetrahedraP2")')
      write(20,'(i10)')nT4P2
      do iCell=1,nT4P2
        write(20,'(10(i10,1x),2x,i3)')tetra2(:,iCell)
      enddo
    endif
    
    if( .not.nQ4==0 )then
      write(20,'(/"Quadrilaterals")')
      write(20,'(i10)')nQ4
      do iCell=1,nQ4
        write(20,'(4(i10,1x),2x,i3)')quadr(:,iCell)
      enddo
    endif
    if( .not.nQ4P2==0 )then
      write(20,'(/"QuadrilateralsQ2")')
      write(20,'(i10)')nQ4P2
      do iCell=1,nQ4P2
        write(20,'(9(i10,1x),2x,i3)')quadr2(:,iCell)
      enddo
    endif
    if( .not.nQ4P3==0 )then
      write(20,'(/"QuadrilateralsQ3")')
      write(20,'(i10)')nQ4P3
      do iCell=1,nQ4P3
        write(20,'(17(i10,1x),2x,i3)')quadr3(:,iCell)
      enddo
    endif
    if( .not.nQ4P4==0 )then
      write(20,'(/"QuadrilateralsQ4")')
      write(20,'(i10)')nQ4P4
      do iCell=1,nQ4P4
        write(20,'(26(i10,1x),2x,i3)')quadr4(:,iCell)
      enddo
    endif
    
    if( .not.nT3==0 )then
      write(20,'(/"Triangles")')
      write(20,'(i10)')nT3
      do iCell=1,nT3
        write(20,'(3(i10,1x),2x,i3)')trian(:,iCell)
      enddo
    endif
    if( .not.nT3P2==0 )then
      write(20,'(/"TrianglesP2")')
      write(20,'(i10)')nT3P2
      do iCell=1,nT3P2
        write(20,'(6(i10,1x),2x,i3)')trian2(:,iCell)
      enddo
    endif
    if( .not.nT3P3==0 )then
      write(20,'(/"TrianglesP3")')
      write(20,'(i10)')nT3P3
      do iCell=1,nT3P3
        write(20,'(11(i10,1x),2x,i3)')trian3(:,iCell)
      enddo
    endif
    if( .not.nT3P4==0 )then
      write(20,'(/"TrianglesP4")')
      write(20,'(i10)')nT3P4
      do iCell=1,nT3P4
        write(20,'(16(i10,1x),2x,i3)')trian4(:,iCell)
      enddo
    endif
    
    if( .not.nL2==0 )then
      write(20,'(/"Edges")')
      write(20,'(i10)')nL2
      do iCell=1,nL2
        write(20,'(2(i10,1x),2x,i3)')edges(:,iCell)
      enddo
    endif
    if( .not.nL2P2==0 )then
      write(20,'(/"EdgesP2")')
      write(20,'(i10)')nL2P2
      do iCell=1,nL2P2
        write(20,'(3(i10,1x),2x,i3)')edges2(:,iCell)
      enddo
    endif
    if( .not.nL2P3==0 )then
      write(20,'(/"EdgesP3")')
      write(20,'(i10)')nL2P3
      do iCell=1,nL2P3
        write(20,'(4(i10,1x),2x,i3)')edges3(:,iCell)
      enddo
    endif
    if( .not.nL2P4==0 )then
      write(20,'(/"EdgesP4")')
      write(20,'(i10)')nL2P4
      do iCell=1,nL2P4
        write(20,'(5(i10,1x),2x,i3)')edges4(:,iCell)
      enddo
    endif
    
    write(20,'(/"End")')
    
    close(20)
    !print '("<<< Writing")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(vert)
    if( allocated(hexas ) )deallocate(hexas )
    if( allocated(hexas2) )deallocate(hexas2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  end subroutine mshToMesh  
  
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
  use baseSimplex1D, only: setQ4BasisEqui_uv,setQ4MeshIJK
  use baseSimplex2D, only: setT3BasisEqui_uv,setT3MeshIJK,setT3MeshBasis_P1,setT3MeshBasis_P2,setT3MeshBasis_P3
  
  use variablesCommunes
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
  integer          :: nQ4,nT3
  integer          :: iCell
  real(8), pointer :: uv0(:,:),uQ4(:),vQ4(:),uvT3(:,:)
  integer          :: iDistantPoint
  integer          :: iVert
  integer, pointer :: ij(:,:)
  real(8), pointer :: lagrangeMeshQ4(:,:)
  real(8), pointer :: lagrangeMeshT3(:,:)
  integer          :: nod(1:16)
  real(8)          :: delta(1:3)
  character(2048)  :: buffer
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !write(buffer,'(a,">>> userInterpolation (callback)")')char(10) ; call msg2(buffer)  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> dist_uvw(:) -> uv0(1:2,:) (sans dupliquer le bloc mémoire)
  !> la commmande c_f_pointer a besoin du use iso_c_binding, only: c_loc,c_f_pointer
  call c_f_pointer(cptr=c_loc(dist_uvw), fptr=uv0, shape=[2,nDistantPoint])  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Setting nQ4 nT3
  nQ4=0 ; nT3=0
  do iDistantPoint=1,nDistantPoint
    iCell=disPtsLocation(iDistantPoint)
    nMod=localConnectivityIndex(iCell+1)-localConnectivityIndex(iCell)
    if    ( nMod==(order+1)*(order+1)   )then ; nQ4=nQ4+1
    elseif( nMod==(order+1)*(order+2)/2 )then ; nT3=nT3+1
    else ; call stopAlert("userInterpolation")
    endif
  enddo  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Allocating uvQ4,uvT3
  if( .not.nQ4==0 )allocate(uQ4(1:nQ4),vQ4(1:nQ4))
  if( .not.nT3==0 )allocate(uvT3(1:2,1:nT3))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Setting uvQ4 uvT3
  nQ4=0 ; nT3=0
  do iDistantPoint=1,nDistantPoint
    iCell=disPtsLocation(iDistantPoint)
    nMod=localConnectivityIndex(iCell+1)-localConnectivityIndex(iCell)
    if    ( nMod==(order+1)*(order+1)   )then ; nQ4=nQ4+1 ; uQ4(nQ4)=2d0*uv0(1,iDistantPoint)-1d0 ; vQ4(nQ4)=2d0*uv0(2,iDistantPoint)-1d0
    elseif( nMod==(order+1)*(order+2)/2 )then ; nT3=nT3+1 ; uvT3(1:2,nT3)=uv0(1:2,iDistantPoint)
    endif
  enddo    
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Computing lagrangeMeshQ4
  if( .not.nQ4==0 )then
    nMod=(order+1)*(order+1)
    allocate(ij(1:2,1:nMod))
    call setQ4MeshIJK(meshOrder=order,ij=ij)
    allocate(lagrangeMeshQ4(1:nMod,1:nQ4))
    call setQ4BasisEqui_uv(ord=order,ijk=ij,u=uQ4,v=vQ4,ai=lagrangeMeshQ4)
    deallocate(ij)
    deallocate(uQ4,vQ4)    
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Computing lagrangeMeshT3
  if( .not.nT3==0 )then
    nMod=(order+1)*(order+2)/2
    allocate(lagrangeMeshT3(1:nMod,1:nT3))
    
    select case(order)
    case(1) ; call setT3MeshBasis_P1(uv=uvT3,ai=lagrangeMeshT3) !> base Triangle Geometrique P1
    case(2) ; call setT3MeshBasis_P2(uv=uvT3,ai=lagrangeMeshT3) !> base Triangle Geometrique P2
    case(3) ; call setT3MeshBasis_P3(uv=uvT3,ai=lagrangeMeshT3) !> base Triangle Geometrique P3
    case default
      allocate(ij(1:2,1:nMod))
      call setT3MeshIJK(meshOrder=order,ij=ij)
      call setT3BasisEqui_uv(ord=order,ijk=ij,uvw=uvT3,ai=lagrangeMeshT3)
      deallocate(ij)
    end select
    deallocate(uvT3)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Computing distantField
  
  nQ4=0 ; nT3=0 ; j=0
  do iDistantPoint=1,nDistantPoint
    distantField(j+1:j+stride)=0d0
    
    iCell=disPtsLocation(iDistantPoint)
    nMod=localConnectivityIndex(iCell+1)-localConnectivityIndex(iCell)
    nod(1:nMod)=localConnectivity(localConnectivityIndex(iCell)+1:localConnectivityIndex(iCell+1))
    
    if    ( nMod==(order+1)*(order+1)   )then ; nQ4=nQ4+1
      
      do iMod=1,nMod
        i=stride*(nod(iMod)-1)
        distantField(j+1:j+stride)=distantField(j+1:j+stride)+lagrangeMeshQ4(iMod,nQ4)*localField(i+1:i+stride)
        
        !> Controling localField
        k=3*(nod(iMod)-1)
        delta(1:3)=localField(i+1:i+3)-localCoordinates(k+1:k+3)
        if( .not.(delta(1)==0d0.and.delta(2)==0d0.and.delta(3)==0d0) )then
          print '("iCell=",i6," Delta=",3(e22.15,1x),t130,"@rkw",i3)',iCell,delta(1:3),rankWorld
        endif
      enddo
      
    elseif( nMod==(order+1)*(order+2)/2 )then ; nT3=nT3+1
      
      do iMod=1,nMod
        i=stride*(nod(iMod)-1)
        distantField(j+1:j+stride)=distantField(j+1:j+stride)+lagrangeMeshT3(iMod,nT3)*localField(i+1:i+stride)
        
        !> Controling localField
        k=3*(nod(iMod)-1)
        delta(1:3)=localField(i+1:i+3)-localCoordinates(k+1:k+3)
        if( .not.(delta(1)==0d0.and.delta(2)==0d0.and.delta(3)==0d0) )then
          print '("iCell=",i6," Delta=",3(e22.15,1x),t130,"@rkw",i3)',iCell,delta(1:3),rankWorld
        endif
      enddo
      
    endif
    j=j+stride
    k=k+3
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Deallocating lagrangeMeshQ4,lagrangeMeshT3
  if( .not.nQ4==0 )deallocate(lagrangeMeshQ4)
  if( .not.nT3==0 )deallocate(lagrangeMeshT3)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
#if 0==1
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Visu focusing on one specific point
  if( nDistantPoint== 1)then
    iDistantPoint=1
    j=stride*(iDistantPoint-1)
    iCell=disPtsLocation(iDistantPoint)
    nMod=localConnectivityIndex(iCell+1)-localConnectivityIndex(iCell)
    nod(1:nMod)=localConnectivity(localConnectivityIndex(iCell)+1:localConnectivityIndex(iCell+1))    
    
    i=stride*(iDistantPoint-1)
    j=     3*(iDistantPoint-1)
    delta(1:3)=disPtsCoordinates(j+1:j+3)-distantField(i+1:i+3)
    
    if( order==1 )then
      
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
      
    elseif( order==2 )then
      
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
      
    elseif( order==3 )then
      
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
#endif
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !call mpi_barrier(commWorld,iErr)
  !if( rankWorld==0 )print'(3x,"<<< userInterpolation")'
  !call mpi_barrier(commWorld,iErr)
  
  !write(buffer,'(a,"<<< userInterpolation (callback)")')char(10) ; call msg2(buffer)  
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
  use additionnal_Functions, only: mshToMesh
  use spaceCellTypes
  use spaceMessages
  use baseSimplex1D, only: nodes1D,setQ4BasisEqui_u ,setQ4MeshIJK
  use baseSimplex2D, only: nodes2D,setT3BasisEqui_uv,setT3MeshIJK, setT3MeshBasis_P1,setT3MeshBasis_P2,setT3MeshBasis_P3
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
  
  character(16)      :: couplingName
  character(5)       :: codeName,codeCoupledName  
  character(128)     :: meshName
  
  integer             :: compOrder
  integer             :: meshOrder
  
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
  write(couplingName,'(a)')"testPiPj"
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !if( rankWorld==0) print '(/"START: fortran_surf_TriaPi_PiPj")'
  write(buffer,'("START: fortran_surf_TriaPi_PiPj")') ; call msg2(trim(buffer))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do meshOrder=4,4
    
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
!    select case(meshOrder)
!    case(1) ; tol=5d-1
!    case(2) ; tol=1d-1
!    case(3) ; tol=1d-1
!    end select
    
    select case(meshOrder)
    case(1) ; tol=2d-1 !7.8d-2
   !case(2) ; tol=1d-2
   ! case(2) ; tol=1.75d-3
    case(2) ; tol=1.d-2
    case(3) ; tol=1.d-2
    case(4) ; tol=1.d-1
    end select
    
    write(buffer,'("")')                                                                                        ; call msg2(trim(buffer))
    write(buffer,'("Code: ",a," creates coupling: ",a," with code: ",a," and tol=",e22.15,t130,"@rkw",i3)') &
    &  trim(codeName),trim(couplingName),trim(codeCoupledName),tol,rankWorld  ; call msg1(trim(buffer))
      
    call cwipi_create_coupling_f(                  &
    &    couplingName=trim(couplingName)          ,&
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
    case(0) ; compOrder=10 !07 !07
    case(1) ; compOrder=10 !07 !10
    end select
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Création des maillages avec gmsh
    
   !write(buffer,'("")')                                                   ; call msg2(trim(buffer))
   !write(buffer,'(3x,"meshOrder=",i3,t130,"@rkw",i3)')meshOrder,rankWorld ; call msg1(trim(buffer))
        
    select case(rankWorld)
    case(0)
      write(key,'(6x,"gmsh ../spaceBasis/tests/Mesh2D/sphere01.geo -2 -format msh -order ",i1," > sphere01.log")')meshOrder      
     !call msg1(trim(key))
      call execute_command_line (key, exitstat=iErr)      
      call mshToMesh(" ../spaceBasis/tests/Mesh2D/sphere01.msh")
      
      !write(key,'(6x,"~/Maillages/mshTomesh spaceBasis/tests/Mesh2D/sphere01.msh >> sphere01.log")')
      !call msg1(trim(key))
      !call execute_command_line (key, exitstat=iErr)     
    case(1)
      write(key,'(6x,"gmsh ../spaceBasis/tests/Mesh2D/sphere02.geo -2 -format msh -order ",i1," > sphere02.log")')meshOrder
     !call msg1(trim(key))
      call execute_command_line (key, exitstat=iErr)       
      call mshToMesh(" ../spaceBasis/tests/Mesh2D/sphere02.msh")
     
      !write(key,'(6x,"~/Maillages/mshTomesh spaceBasis/tests/Mesh2D/sphere02.msh >> sphere01.log")')
      !call msg1(trim(key))
      !call execute_command_line (key, exitstat=iErr)
    end select
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Reading Geometric Mesh with INRIA libMesh7
    !> attention avec gmfGetBlock certains entiers sont integer(8)
    
    write(buffer,'("")')                                               ; call msg2(trim(buffer))
    write(buffer,'("Reading Geometric Mesh",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))
    select case(rankWorld)
    case(0) ; meshName="../spaceBasis/tests/Mesh2D/sphere01.mesh"
    case(1) ; meshName="../spaceBasis/tests/Mesh2D/sphere02.mesh"
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
    
   !write(buffer,'(3x,"meshOrder=",i2,"  nVert=",i6," nQ4=",i6," nT3=",i6,t130,"@rkw",i3)')meshOrder,nVert,nQ4,nT3,rankWorld ; call msg1(trim(buffer))  
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
      case(3)
        res=GmfGetBlock(                               &
        &   InpMsh                                    ,&
        &   GmfQuadrilateralsQ3                       ,&  ! <=
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
        &   GmfInt, cells(ad1+11) , cells(ad2+11)     ,&
        &   GmfInt, cells(ad1+12) , cells(ad2+12)     ,&
        &   GmfInt, cells(ad1+13) , cells(ad2+13)     ,&
        &   GmfInt, cells(ad1+14) , cells(ad2+14)     ,&
        &   GmfInt, cells(ad1+15) , cells(ad2+15)     ,&
        &   GmfInt, cells(ad1+16) , cells(ad2+16)     ,&
        &                                              &
        &   GmfInt, mark(iCell0+1), mark(iCell0+nCell) )
      case(4)
        res=GmfGetBlock(                               &
        &   InpMsh                                    ,&
        &   GmfQuadrilateralsQ4                       ,&  ! <=
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
        &   GmfInt, cells(ad1+11) , cells(ad2+11)     ,&
        &   GmfInt, cells(ad1+12) , cells(ad2+12)     ,&
        &   GmfInt, cells(ad1+13) , cells(ad2+13)     ,&
        &   GmfInt, cells(ad1+14) , cells(ad2+14)     ,&
        &   GmfInt, cells(ad1+15) , cells(ad2+15)     ,&
        &   GmfInt, cells(ad1+16) , cells(ad2+16)     ,&
        &   GmfInt, cells(ad1+17) , cells(ad2+17)     ,&
        &   GmfInt, cells(ad1+18) , cells(ad2+18)     ,&
        &   GmfInt, cells(ad1+19) , cells(ad2+19)     ,&
        &   GmfInt, cells(ad1+20) , cells(ad2+20)     ,&
        &   GmfInt, cells(ad1+21) , cells(ad2+21)     ,&
        &   GmfInt, cells(ad1+22) , cells(ad2+22)     ,&
        &   GmfInt, cells(ad1+23) , cells(ad2+23)     ,&
        &   GmfInt, cells(ad1+24) , cells(ad2+24)     ,&
        &   GmfInt, cells(ad1+25) , cells(ad2+25)     ,&
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
      case(4)
        res=GmfGetBlock(                               &
        &   InpMsh                                    ,&
        &   GmfTrianglesP4                            ,&  ! <=
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
        &   GmfInt, cells(ad1+11) , cells(ad2+11)     ,&
        &   GmfInt, cells(ad1+12) , cells(ad2+12)     ,&
        &   GmfInt, cells(ad1+13) , cells(ad2+13)     ,&
        &   GmfInt, cells(ad1+14) , cells(ad2+14)     ,&
        &   GmfInt, cells(ad1+15) , cells(ad2+15)     ,&
        &                                              &
        &   GmfInt, mark(iCell0+1), mark(iCell0+nCell) )
        
      case default ; call stopAlert("reading Triangles meshOrder>4")
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
    
   !write(buffer,'("")')                                              ; call msg2(trim(buffer))
   !write(buffer,'("Sending Mesh to Cwipi",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))
    
    nCell=nQ4+nT3
    
    !> Transmission des maillages à cwipi
    call cwipi_ho_define_mesh_f(         & !> NEW Cwipi
    &   couplingName=trim(couplingName) ,&
    &   nVertex     =nVert              ,&
    &   nElts       =nQ4+nT3            ,&
    &   order       =meshOrder          ,&  
    &   coords      =vertx              ,&
    &   connecIndex =cellsIdx           ,&
    &   connec      =cells               )  
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Definition du Quad Qi (ijk)
    if( .not.nQ4==0 )then
      nMod=(meshOrder+1)*(meshOrder+1)
      allocate(ij(1:2,nMod))
      
      call setQ4MeshIJK(meshOrder=meshOrder,ij=ij)
      call c_f_pointer(cptr=c_loc(ij), fptr=ijCwipi, shape=[2*nMod])  

      call cwipi_ho_ordering_from_IJK_set_f( & !> NEW Cwipi
      &   couplingName =trim(couplingName)  ,&
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
      &   couplingName =trim(couplingName)  ,&
      &   tElt         = CWIPI_FACE_TRIAHO  ,&
      &   nNodes       = nMod               ,&
      &   IJK          = ijCwipi             )
      
      deallocate(ij)
    endif
    
    !write(buffer,'("")')                                                             ; call msg2(trim(buffer))
    !write(buffer,'("Geometric HO Cell ij P",i1,t130,"@rkw",i3)')meshOrder,rankWorld  ; call msg1(trim(buffer))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Initialisation of myValues(:)
    
    !write(buffer,'("")')                                                  ; call msg2(trim(buffer))
    !write(buffer,'("Init myValues (x,y,z,rkw)",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))  
    
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
    
    
#if 0==0
    
    !> calcul lagrangeMeshQ4
    if( .not.nQ4==0 )then      
      nMod=(meshOrder+1)*(meshOrder+1)                  !> Quad meshOrder (entree)
      call nodes1D(ord=compOrder,uvw=u,display=.false.) !> ordre du calcul  
      nNod=size(u) ; nNod=nNod*nNod                     !> tensorise
      allocate(lagrangeMeshQ4(1:nMod,1:nNod))
      allocate(ij(1:2,1:nMod))
      call setQ4MeshIJK(meshOrder=meshOrder,ij=ij)
      call setQ4BasisEqui_u(ord=meshOrder,ijk=ij,u=u,ai=lagrangeMeshQ4)
      deallocate(ij,u)
      
!      block
!      real(8) :: toto(11)
!      print '(/"lagrangeMeshQ4")'
!      do iMod=1,nMod
!        toto(1:nNod)=lagrangeMeshQ4(iMod,1:nNod)
!        print '(3x,"ai(",i2,",1:nNod)=",*(e12.5,1x) )',iMod,toto(1:nNod)
!      enddo
!      print '("")'
!      end block
      
    endif
    
    !> calcul lagrangeMeshT3
    if( .not.nT3==0 )then
      nMod=(meshOrder+1)*(meshOrder+2)/2                 !> Triangle meshOrder (entree)
      call nodes2D(ord=compOrder,uvw=uv,display=.false.) !> ordre du calcul  
      nNod=size(uv,2)
      allocate(lagrangeMeshT3(1:nMod,1:nNod))
      select case(meshOrder)
      case(01) ; call setT3MeshBasis_P1(uv=uv,ai=lagrangeMeshT3)
      case(02) ; call setT3MeshBasis_P2(uv=uv,ai=lagrangeMeshT3)
      case(03) ; call setT3MeshBasis_P3(uv=uv,ai=lagrangeMeshT3)
      case default
        allocate(ij(1:2,1:nMod))
        call setT3MeshIJK(meshOrder=meshOrder,ij=ij)
        call setT3BasisEqui_uv(ord=meshOrder,ijk=ij,uvw=uv,ai=lagrangeMeshT3)
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
        
        !if( iCell==1 )then
        !  do iMod=1,nMod
        !    k=3*(nod(iMod)-1)
        !    xyz(1:3)=vertx(k+1:k+3)
        !    print '(3x"nod(",i8,")=",3(e22.15,1x))',nod(iMod),xyz(1:3)
        !  enddo
        !endif
        
        do iNod=1,nNod
          !> linkVert
          linkVert(j+1:j+3)=0d0
          do iMod=1,nMod
            k=3*(nod(iMod)-1)
            linkVert(j+1:j+3)=linkVert(j+1:j+3)+lagrangeMeshQ4(iMod,iNod)*vertx(k+1:k+3)            
          enddo
          
          !if( iCell==1 )then
          !  xyz(1:3)=linkVert(j+1:j+3)
          !  print '(3x"linkVert(",i3,")=",3(e22.15,1x))',iNod,xyz(1:3)
          !endif
          
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
    
    write(buffer,'("")') ; call msg2(trim(buffer))
    write(buffer,'("linkVert compOrder=",i3," -> linkVertSize=",i6,t130,"@rkw",i3)')compOrder,linkVertSize,rankWorld ; call msg1(trim(buffer))
    
#else
    
    linkVertSize=1                       !> nombre total de point de couplages
    allocate(linkVert(1:3*linkVertSize)) !> 3 coordonnées par point de couplage
    
    select case(rankWorld)
!    case(0) ; linkVert(1:3)=[-0.919633675189875E+00,-0.250163898379974E-01, 0.392131352367653E+00]  
 !   case(1) ; linkVert(1:3)=[-0.848807623678876E+00,-0.213274154008489E-01, 0.424121411593377E+00]
    case(0) ; linkVert(1:3)=[-0.962125889866290E+00, -0.353822407491257E-01, -0.270719711354757E+00]  
    case(1) ; linkVert(1:3)=[0.690377557847734E+00,  0.515357183278473E+00, -0.519234683615825E+00]
    end select
    
    write(buffer,'(6x,"linkVert(1:3)=    ",3(e22.15,1x),t130,"@rkw",i3)')linkVert(1:3),rankWorld     ; call msg1(trim(buffer))
    
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !write(buffer,'("")')                                                                              ; call msg2(trim(buffer))
    !write(buffer,'("Allocation of linkValues size=",i6,t130,"@rkw",i3)')stride*linkVertSize,rankWorld ; call msg1(trim(buffer))
    
    allocate(linkValues(1:stride*linkVertSize))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transmission à cwipi des coordonnees de couplage
    
    !write(buffer,'("")')                                                         ; call msg2(trim(buffer))
    !write(buffer,'("Transmission de linkVert a cwipi",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))  
      
    call cwipi_set_points_to_locate_f(    &
    &    couplingName=trim(couplingName) ,&
    &    nPts  =  linkVertSize           ,&
    &    coords=linkVert                  )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Localisation par cwipi des coordonnees de couplage
    
    write(buffer,'("")')                                              ; call msg2(trim(buffer))
    write(buffer,'("Localisation by Cwipi",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))  
    
    call cwipi_locate_f(couplingName=trim(couplingName))
    
    call cwipi_get_n_not_located_pts_f(couplingName=trim(couplingName), nNotLocatedPoints=numberOfUnlocatedPoints)
    
    call mpi_allreduce(numberOfUnlocatedPoints,numberOfUnlocatedPointsGlob,1,mpi_integer,mpi_sum,commWorld,iErr)
    
    if( .not.numberOfUnlocatedPointsGlob==0 )then
      write(buffer,'(3x,"nNotLocatedPoints=",i10,t130,"@rkw",i3)')numberOfUnlocatedPoints,rankWorld ; call msg1(trim(buffer))  
      write(buffer,'("fortran_surf_TriaPi_PiPj stop line: ",i6,t130,"@rkw",i3)')__LINE__+1,rankWorld
      call stopAlert(trim(buffer))  
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !write(buffer,'("")')                                                 ; call msg2(trim(buffer))
    !write(buffer,'("Echange cwipi_exchange_f",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))  
    
    call cwipi_exchange_f(                          &
    &    couplingName=          trim(couplingName) ,&
    &    exchangeName="exch1_"//trim(couplingName) ,&
    &    meshOrder=meshOrder                       ,&  !> NEW cwipi
    &    exchangeDim=stride                        ,&  !> scalar
    &    ptHoInterpolationFct=userInterpolation    ,&  !> utilisation de la procedure plug
    &                                               &
    &    sendingFieldName="mySolu"                 ,&  !> solution calculee localement
    &    sendingField=myValues                     ,&
    &                                               &
    &    receivingFieldName="linkSolu"             ,&
    &    receivingField=linkValues                 ,&  !> solution de raccord
    &                                               &
    &    nStep=1                                   ,&  !> pas utilisee juste pour visu cwipi
    &    timeValue=0d0                             ,&  !> pas utilisee juste pour visu cwipi
    &    nNotLocatedPoints=notLocatedPoints        ,&
    &    status=iErr                                )  
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !write(buffer,'("")')                                        ; call msg2(trim(buffer))
    !write(buffer,'("Delete coupling",t130,"@rkw",i3)')rankWorld ; call msg1(trim(buffer))  
    
    call cwipi_delete_coupling_f(trim(couplingName))
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
    write(buffer,'(                                     a, &
    &              3x,"Control",t130,"@rkw",i3         ,a, &
    &              6x,"mesh: ",a," meshOrder=",i1      ,a, &
    &              6x,"nVert=",i6," nQ4=",i6," nT3=",i6,a, &
    &                                                   a, &
    &              6x,"linkVertSize=",i6               ,a, &
    &                                                   a, &
    &              6x,"deltaMin  =",e22.15             ,a, &
    &              6x,"sumDelta  =",e22.15             ,a, &
    &              6x,"deltaMax  =",e22.15             ,a, &
    &                                                   a, &
    &              6x,"linkVert  =",3(e22.15,1x)       ,a, &
    &              6x,"linkValues=",3(e22.15,1x)       ,a, &
    &              6x,"Delta     =",3(e22.15,1x)           &
    &                                                   )')&
    &                          char(10),&
    & rankWorld               ,char(10),&
    & trim(meshName)                   ,&    
    & meshOrder               ,char(10),&
    & nVert,nQ4,nT3           ,char(10),&
    &                          char(10),&
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
