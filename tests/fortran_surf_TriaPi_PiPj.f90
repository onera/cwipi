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


module spaceMessages

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use iso_fortran_env
  use mpi
  use variablesCommunes
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  

  implicit none
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
 !include 'mpif.h'
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
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
    write(output_unit ,'(/"SPACE ALERT: ",a,"  rkw",i3," called mpi_abort")')trim(msg),rankWorld
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
  &           uvw_size                                ,& !> new cwipi
  &           dist_uvw                                ,& !> new cwipi
  &                                                    &
  &           stride                                  ,&  ! =ker(calc)
  &           solverType                              ,&
  &           localField                              ,&  !   mySolu
  &           distantField                             )  ! linkSolu
  !---
  use iso_c_binding, only: c_loc,c_f_pointer
  use cwipi
  use baseSimplex2D
  
  use variablesCommunes
  use additionnal_Functions
  use spaceMessages
  !---
  implicit none
  !---
  integer :: entitiesDim
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
  integer :: uvw_size
  real(8) :: dist_uvw                                (*)
  integer :: stride
  integer :: solverType
  real(8) ::   localField                            (*)
  real(8) :: distantField                            (*)
  !>
  integer          :: i,j,k,iRank,iErr
  integer          :: iNod,nNod,iMod,nMod
  integer          :: iCell
  real(8), pointer :: uv(:,:),a(:),b(:)!c(:)
  integer          :: iDistantPoint
  integer          :: iVert
  real(8), pointer :: lagrangeMesh(:,:)
  integer          :: nod(1:16)
  real(8)          :: delta(1:3)
  character(1024)  :: buffer
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
  case(2) ; call setT3MeshBasis_P2(uv=uv,ai=lagrangeMesh) !> base Triangle Geometrique P2
  case default
    call nodes2Duv2ab(uv=uv,a=a,b=b ,display=.false.) !> rs(1:2,:)=2d0*uv(1:2,:)-1d0 && a=2 (1+r)/(1-s)-1 && b=s
    call lagrange2Dv(      &
    &    ord=meshOrder    ,&
    &    vand=vand        ,&
    &    a=a,b=b          ,&
    &    lx=lagrangeMesh  ,&
    &    transpose=.true.  ) !> lagrangeMesh(1:nMod,1:nNod) nNod=size(u)
    deallocate(a,b)    
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
    delta(1:3)=distantField(i+1:i+3)-disPtsCoordinates(k+1:k+3)
    
!    write(buffer,'(3x,"disPtsCoordinates=",3(e22.15,1x),t130,"@rkw",i3)')disPtsCoordinates(j+1:j+3),rankWorld ; call msg1(buffer)
!    write(buffer,'(3x,"distantField=     ",3(e22.15,1x),t130,"@rkw",i3)')distantField     (i+1:i+3),rankWorld ; call msg1(buffer)
!    write(buffer,'(3x,"delta=            ",3(e22.15,1x),t130,"@rkw",i3)')disPtsCoordinates(j+1:j+3)-distantField(i+1:i+3),rankWorld ; call msg1(buffer)
    
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
  
  use mpi
  use cwipi  
  use baseSimplex2D
  
  use variablesCommunes
  use additionnal_Functions
  use spaceMessages
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  interface
    subroutine  userInterpolation                      ( &
    &           entitiesDim                             ,&
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
    &           uvw_size                                ,&
    &           dist_uvw                                ,&
    &                                                    &
    &           stride                                  ,&  ! =ker(calc)
    &           solverType                              ,&
    &           localField                              ,&  !   mySolu
    &           distantField                             )  ! linkSolu
    !---
    integer :: entitiesDim
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
    integer :: uvw_size
    real(8) :: dist_uvw(*)
    integer :: stride
    integer :: solverType
    real(8) ::   localField                            (*)
    real(8) :: distantField                            (*)
    end subroutine  userInterpolation
  end interface
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  character(5)     :: codeName,codeCoupledName  
  character(128)   :: meshName
  
  integer          :: iVert,nVert
  real(8), pointer :: vertx(:)
  integer, pointer :: cells(:),cellsIdx(:)
  integer, pointer :: nod(:)  !> ensemble des noeuds pour iCell : cel=>cells(cellsIdx(iCell)+1:cellsIdx(iCell+1))

  
  character(80)      :: key
  integer, parameter :: iFile=100
  logical            :: test
  
  
  integer          :: i,j,k
  integer          :: iMod,nMod
  integer          :: iNod,nNod
  integer          :: iCell,nCell
  
  integer, pointer :: ijk(:)
  real(8), pointer :: lagrangeMesh(:,:)
  real(8), pointer :: xi(:,:)
  real(8), pointer :: a(:),b(:)
  
  real(8)          :: xyz(1:3)
  integer          :: linkVertSize
  real(8), pointer :: linkVert(:)
  integer          :: notLocatedPoints
  
  integer          :: stride
  real(8), pointer ::   myValues(:)
  real(8), pointer :: linkValues(:)
  
  real(8), pointer :: uv(:,:)
  
  integer          :: iRank,iErr
  
  integer          :: iVertMax
  real(8)          :: delta,deltaMin,deltaMax,sumDelta
  character(1024)  :: buffer
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_init(iErr)
  commWorld=mpi_comm_world
  
  call mpi_comm_rank(commWorld, rankWorld, iErr)
  call mpi_comm_size(commWorld, sizeWorld, iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !if( rankWorld==0) print '(/"START: fortran_surf_TriaPi_PiPj")'
  write(buffer,'("START: fortran_surf_TriaPi_PiPj")') ; call msg2(buffer)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
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
  !> Create coupling
!  call mpi_barrier(commWorld,iErr)
!  if( rankWorld==0 )print '(/"Creation du couplage testPiPj")'  
!  call mpi_barrier(commWorld,iErr)
  
  write(buffer,'("")')                                         ; call msg2(buffer)
  write(buffer,'("Create coupling",t130,"@rkw",i3)')rankWorld  ; call msg1(buffer)
    
  call cwipi_create_coupling_f(                  &
  &    couplingName="testPiPj"                  ,&
  &    couplingType=cwipi_cpl_parallel_with_part,&
  &    cplAppli=codeCoupledName                 ,&
  &    entitiesDim=2                            ,& !> Nature du couplage surfacique
  &    tolerance=1d-3                           ,& !> Tolerance geometrique 1d-1 par defaut
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
  
  write(buffer,'("")')                                                   ; call msg2(buffer)
  write(buffer,'(3x,"compOrder=",i3,t130,"@rkw",i3)')compOrder,rankWorld ; call msg1(buffer)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Création des maillages avec gmsh
  select case(rankWorld)
  case(0) ; meshOrder=2 !07 !07
  case(1) ; meshOrder=2 !07 !10
  end select
  
  write(buffer,'("")')                                                   ; call msg2(buffer)
  write(buffer,'(3x,"meshOrder=",i3,t130,"@rkw",i3)')meshOrder,rankWorld ; call msg1(buffer)
  
  select case(rankWorld)
  case(0)
    write(key,'("gmsh spaceBasis/tests/Mesh2D/sphere01.geo -2 -format msh -order ",i1," > sphere01.log")')meshOrder
    call execute_command_line (key, exitstat=iErr)
    write(key,'("~/Maillages/mshTomesh spaceBasis/tests/Mesh2D/sphere01.msh >> sphere01.log")')
    call execute_command_line (key, exitstat=iErr)
  case(1)
    write(key,'("gmsh spaceBasis/tests/Mesh2D/sphere02.geo -2 -format msh -order ",i1," > sphere02.log")')meshOrder
    call execute_command_line (key, exitstat=iErr)
    write(key,'("~/Maillages/mshTomesh spaceBasis/tests/Mesh2D/sphere02.msh >> sphere02.log")')
    call execute_command_line (key, exitstat=iErr)
  end select
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Reading Geometric Mesh
  
  write(buffer,'("")')                                                                               ; call msg2(buffer)
  write(buffer,'("Reading Geometric Mesh",t130,"@rkw",i3)')rankWorld ; call msg1(buffer)
  select case(rankWorld)
  case(0) ; meshName="./spaceBasis/tests/Mesh2D/sphere01.mesh"
  case(1) ; meshName="./spaceBasis/tests/Mesh2D/sphere02.mesh"
  end select
    
  open(unit=iFile,file=trim(meshName),status='old',action='read')  
  lecture: do
    read(unit=iFile,fmt=*)key ! print '("key: ",a,t130,"@rkw",i3)',trim(key),rankWorld
    
    select case(trim(key))
    case("Vertices","vertices")
      
      read(iFile,*)nVert
      allocate( vertx(1:3*nVert) )  !> sommets
      do iVert=1,nVert
        read(iFile,*)xyz(1:3)!x,y,z
        vertx(3*(iVert-1)+1:3*iVert)=xyz(1:3)
      enddo
      
    case("Quadrilaterals","QuadrilateralsQ2","QuadrilateralsQ3","Triangles","TrianglesP2","TrianglesP3")
      
      select case(trim(key))
      case("Quadrilaterals"  ) ; meshOrder=1 ; nMod=(meshOrder+1)*(meshOrder+1)
      case("QuadrilateralsQ2") ; meshOrder=2 ; nMod=(meshOrder+1)*(meshOrder+1)
      case("QuadrilateralsQ3") ; meshOrder=3 ; nMod=(meshOrder+1)*(meshOrder+1)
      case("Triangles"       ) ; meshOrder=1 ; nMod=(meshOrder+1)*(meshOrder+2)/2
      case("TrianglesP2"     ) ; meshOrder=2 ; nMod=(meshOrder+1)*(meshOrder+2)/2
      case("TrianglesP3"     ) ; meshOrder=3 ; nMod=(meshOrder+1)*(meshOrder+2)/2
      end select
      
      read(iFile,*)nCell
      allocate( cells   (1:nMod*nCell) )
      allocate( cellsIdx(1:nCell+1   ) )
      cellsIdx(1)=0
      do iCell=1,nCell
        cellsIdx(iCell+1)=cellsIdx(iCell)+nMod
        read(iFile,*)cells(cellsIdx(iCell)+1:cellsIdx(iCell+1))
      enddo
      
    case("End","end") ; exit lecture
    end select
  enddo lecture
  close(unit=iFile)
  
  write(buffer,'(                                &
  &                                           a, &
  &              3x,"mesh=",a,t130,"@rkw",i3 ,a, &
  &              6x,"meshOrder=",i1          ,a, &
  &              6x,"nVert=",i6              ,a, &
  &              6x,"nCell=",i6                  &
  &                                           )')&
  &                          char(10),&
  & trim(meshName),rankWorld,char(10),&
  & meshOrder               ,char(10),&
  & nVert                   ,char(10),&
  & nCell
  
  call msg1(buffer)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> On se couple sur les sphères maillées
  
  write(buffer,'("")')                                                         ; call msg2(buffer)
  write(buffer,'("Sending Mesh to Cwipi",t130,"@rkw",i3)')rankWorld ; call msg1(buffer)
  
  !> Transmission des maillages à cwipi
  call cwipi_ho_define_mesh_f(  & !> NEW Cwipi
  &   couplingName="testPiPj"  ,&
  &   nVertex     =nVert       ,&
  &   nElts       =nCell       ,&
  &   order       =meshOrder   ,&  
  &   coords      =vertx       ,&
  &   connecIndex =cellsIdx    ,&
  &   connec      =cells        )  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Definition du Triangle géométrique Pi
  
  write(buffer,'("")')                                                             ; call msg2(buffer)
  write(buffer,'("Geometric HO Triangle P",i1,t130,"@rkw",i3)')meshOrder,rankWorld ; call msg1(buffer)
    
  nMod=(meshOrder+1)*(meshOrder+2)/2
  allocate(ijk(1:2*nMod))
  
  if    ( meshOrder==1 )then !> TriangleP1
    
    ! 03
    ! 01 02
    
    ijk( 1: 2)=[0,0] !> 1
    ijk( 3: 4)=[1,0] !> 2
    ijk( 5: 6)=[0,1] !> 3
    
  elseif( meshOrder==2 )then !> TriangleP2
    
    ! 03
    ! 06 05
    ! 01 04 02
    
    ijk( 1: 2)=[0,0] !> 1
    ijk( 3: 4)=[2,0] !> 2
    ijk( 5: 6)=[0,2] !> 3
    ijk( 7: 8)=[1,0] !> 4
    ijk( 9:10)=[1,1] !> 5
    ijk(11:12)=[0,1] !> 6
    
  elseif(meshOrder==3 )then !> TriangleP3
    
    ! 03
    ! 08 07
    ! 09 10 06
    ! 01 04 05 02
    
    ijk( 1: 2)=[0,0] !> 01
    ijk( 3: 4)=[3,0] !> 02
    ijk( 5: 6)=[0,3] !> 03
    ijk( 7: 8)=[1,0] !> 04
    ijk( 9:10)=[2,0] !> 05
    ijk(11:12)=[2,1] !> 06
    ijk(13:14)=[1,2] !> 07
    ijk(15:16)=[0,2] !> 08
    ijk(17:18)=[0,1] !> 09
    ijk(19:20)=[1,1] !> 10
    
  elseif(meshOrder==4 )then !> TriangleP4
  
    !> 03
    !> 10 09
    !> 11 15 08
    !> 12 13 14 07
    !> 01 04 05 06 02
  
    ijk( 1: 2)=[0,0] !> 01
    ijk( 3: 4)=[4,0] !> 02
    ijk( 5: 6)=[0,4] !> 03
    ijk( 7: 8)=[1,0] !> 04
    ijk( 9:10)=[2,0] !> 05
    ijk(11:12)=[3,0] !> 06
    ijk(13:14)=[3,1] !> 07
    ijk(15:16)=[2,2] !> 08
    ijk(17:18)=[1,3] !> 09
    ijk(19:20)=[0,3] !> 10
    ijk(21:22)=[0,2] !> 11
    ijk(23:24)=[0,1] !> 12
    ijk(25:26)=[1,1] !> 13
    ijk(27:28)=[2,1] !> 14
    ijk(29:30)=[1,2] !> 15
    
  else ; stop "meshOrder>4 not implemented"
  endif
  
  call cwipi_ho_ordering_from_IJK_set_f( & !> NEW Cwipi
  &   couplingName ="testPiPj"          ,&
  &   tElt         = CWIPI_FACE_TRIAHO  ,&
  &   nNodes       = nMod               ,&
  &   IJK          = ijk                 )
  
  if( meshOrder>2 )then
    allocate(xi(1:2,1:nMod))
    j=0
    do iNod=1,nMod
      j=j+1 ; xi(1,iNod)=real(ijk(j),kind=8)/real(meshOrder,kind=8)
      j=j+1 ; xi(2,iNod)=real(ijk(j),kind=8)/real(meshOrder,kind=8)
    enddo
    call nodes2Duv2ab (uv=xi,a=a,b=b,display=.false.) ! rs(1:2,:)=2d0*uv(1:2,:)-1d0
    call vandermonde2D(ord=meshOrder,a=a,b=b,vand=vand)
    deallocate(a,b,xi)
  endif
  
  deallocate(ijk)  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Initialisation of myValues(:)
  
  write(buffer,'("")')                                                  ; call msg2(buffer)
  write(buffer,'("Init myValues (x,y,z,rkw)",t130,"@rkw",i3)')rankWorld ; call msg1(buffer)  
  
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
  
  write(buffer,'("")')                                        ; call msg2(buffer)
  write(buffer,'("Set linkVert(:)",t130,"@rkw",i3)')rankWorld ; call msg1(buffer)
  
#if 0==1
  
  !> calcul lagrangeMesh
  nMod=(meshOrder+1)*(meshOrder+2)/2   !> Triangle meshOrder (entree)
  call nodes2D(ord=compOrder,uvw=uv,display=.false.) !> ordre du calcul  
  nNod=size(uv,2)
  allocate(lagrangeMesh(1:nMod,1:nNod))
  select case(meshOrder)
  case(01) ; call setT3MeshBasis_P1(uv=uv,ai=lagrangeMesh)
  case(02) ; call setT3MeshBasis_P2(uv=uv,ai=lagrangeMesh)
  case default      
    call nodes2Duv2ab(uv=uv,a=a,b=b ,display=.false.) !> rs(1:2,:)=2d0*uv(1:2,:)-1d0 && a=2 (1+r)/(1-s)-1 && b=s
    call lagrange2Dv(      &
    &    ord=meshOrder    ,&
    &    vand=vand        ,&
    &    a=a,b=b          ,&
    &    lx=lagrangeMesh  ,&
    &    transpose=.true.  ) !> lagrangeMesh(1:nMod,1:nNod) nNod=size(u)
    deallocate(a,b)
  end select
  deallocate(uv)
  
  !> calcul linkVert
  linkVertSize=nCell*nNod              !> nombre total de point de couplages
  allocate(linkVert(1:3*linkVertSize)) !> 3 coordonnées par point de couplage
  
  j=0
  do iCell=1,nCell
    nod=>cells(cellsIdx(iCell)+1:cellsIdx(iCell+1))
    do iNod=1,nNod
      !> linkVert
      linkVert(j+1:j+3)=0d0
      do iMod=1,nMod
        i=3*(nod(iMod)-1)
        linkVert(j+1:j+3)=linkVert(j+1:j+3)+lagrangeMesh(iMod,iNod)*vertx(i+1:i+3)
      enddo
      j=j+3      
    enddo
  enddo
  
  deallocate(lagrangeMesh)
  
#else
  
  linkVertSize=1              !> nombre total de point de couplages
  allocate(linkVert(1:3*linkVertSize)) !> 3 coordonnées par point de couplage
  
  select case(rankWorld)
  case(0) ; linkVert(1:3)=[-0.151188007905265E+00, 0.353840776353232E-01, 0.987777264278830E+00]  
  case(1) ; linkVert(1:3)=[-0.151188007905265E+00, 0.353840776353232E-01, 0.987777264278830E+00]
  end select
  
  write(buffer,'(6x,"linkVert(1:3)=    ",3(e22.15,1x),t130,"@rkw",i3)')linkVert(1:3),rankWorld ; call msg1(buffer)
  
#endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(buffer,'("")')                                                                              ; call msg2(buffer)
  write(buffer,'("Allocation of linkValues size=",i6,t130,"@rkw",i3)')stride*linkVertSize,rankWorld ; call msg1(buffer)
  
  allocate(linkValues(1:stride*linkVertSize))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Transmission à cwipi des coordonnees de couplage
  
  write(buffer,'("")')                                                               ; call msg2(buffer)
  write(buffer,'("Transmission de linkVert a cwipi",t130,"@rkw",i3)')rankWorld ; call msg1(buffer)  
    
  call cwipi_set_points_to_locate_f( &
  &    couplingName="testPiPj"      ,&
  &    nPts  =linkVertSize          ,&
  &    coords=linkVert               )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Localisation par cwipi des coordonnees de couplage
  
  write(buffer,'("")')                                                               ; call msg2(buffer)
  write(buffer,'("Localisation by Cwipi",t130,"@rkw",i3)')rankWorld ; call msg1(buffer)  
  
  call cwipi_locate_f(couplingName="testPiPj")
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  write(buffer,'("")')                                                 ; call msg2(buffer)
  write(buffer,'("Echange cwipi_exchange_f",t130,"@rkw",i3)')rankWorld ; call msg1(buffer)  
  
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
  write(buffer,'("")')                                        ; call msg2(buffer)
  write(buffer,'("Delete coupling",t130,"@rkw",i3)')rankWorld ; call msg1(buffer)  
  
  call cwipi_delete_coupling_f("testPiPj")
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(buffer,'("")')                                           ; call msg2(buffer)
  write(buffer,'("Controling Results",t130,"@rkw",i3)')rankWorld ; call msg1(buffer)
  
  iVertMax=1
  sumDelta= 0d0
  deltaMax=-1d50
  deltaMin= 1d50
  j=0
  k=0
!  do iRank=0,sizeWorld-1
!    if( iRank==rankWorld )then
!      print '(/3x,"meshOrder =",i2," controling linkValues - linkVertSize=",i6,t130,"@rkw",i3)',meshOrder,linkVertSize,rankWorld
      do iVert=1,linkVertSize
        delta=norm2(linkVert(j+1:j+3)-linkValues(k+1:k+3)) !+( real(rankWorld,kind=8)-linkValues(k+4) )**2
        sumDelta=sumDelta+delta
        
       !if( deltaMax<delta )deltaMax=delta
        
        if( deltaMax<delta )then
          iVertMax=iVert
          deltaMax=delta
          !if( iRank==0 )then
          !  print '("deltaMax=",e22.15," iVert=",i10,t130,"@rkw",i3)',deltaMax,iVert,rankWorld
          !  print '( 3x,"linkVert  (j+1:j+3)=",3(e22.15,1x),t130,"@rkw",i3)',linkVert  (j+1:j+3),rankWorld
          !  print '( 3x,"linkValues(k+1:k+3)=",3(e22.15,1x),t130,"@rkw",i3)',linkValues(k+1:k+3),rankWorld
          !endif
        endif
        if( delta<deltaMin )deltaMin=delta
        j=j+3
        k=k+4
      enddo      
!    endif
!    call mpi_barrier(commWorld,iErr)
!  enddo
  
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
  
  call msg1(buffer)
  
    
  !call mpi_barrier(commWorld,iErr)  
  !do iRank=0,sizeWorld-1
  !  if( iRank==rankWorld )then
  !    print '(/3x,"deltaMin=min( |linkVert-linkValues|^2 )             =",e22.15,t130,"@rkw",i3)',deltaMin,rankWorld
  !    print '( 3x,"sumDelta=sum( |linkVert-linkValues|^2 )/linkVertSize=",e22.15,t130,"@rkw",i3)',sumDelta,rankWorld
  !    print '( 3x,"deltaMax=max( |linkVert-linkValues|^2 )             =",e22.15,t130,"@rkw",i3)',deltaMax,rankWorld
  !  endif
  !  call mpi_barrier(commWorld,iErr)
  !enddo
  
  
  
 !call mpi_allreduce(sumDelta,deltaMax,1,mpi_real8,mpi_max,commWorld,iErr)
 !if( sumDelta<1d-08 )then
 !  if( rankWorld==0 )print '(/"SUCCESS: fortran_surf_TriaPi_PiPj"/)'
 !else
 !  if( rankWorld==0 )print '(/"FAILED: fortran_surf_TriaPi_PiPj"/)'
 !  stop
 !endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(vertx,cellsIdx,cells,myValues)
  deallocate(linkVert,linkValues)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call cwipi_finalize_f()
  call mpi_finalize(iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
end program fortran_surf_TriaPi_PiPj
