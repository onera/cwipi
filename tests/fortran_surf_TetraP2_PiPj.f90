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


!  mpirun -n 1 ./fortran_surf_PiPj : -n 1 ./fortran_surf_PiPj

module additionnal_Functions

contains

  subroutine setT3MeshBasis_P1(u,v,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Numerotation des sommets
    !>   3
    !>   1 2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! delcaration des variables passees en argument
    real(8), intent(in)    :: u,v
    real(8), intent(inout) :: ai(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ai(1)=1d0-u-v
    ai(2)=    u
    ai(3)=      v
   !write(*,'("u,v=",2(f12.5,1x),"li=",3(f12.5,1x))')u,v,ai(1:3)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine setT3MeshBasis_P1
  
  subroutine setT3MeshBasis_P2(u,v,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Numerotation des sommets
    !>   3
    !>   6 5
    !>   1 4 2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! delcaration des variables passees en argument
    real(8), intent(in)    :: u,v
    real(8), intent(inout) :: ai(:)
    !>
    real(8)                :: w
    real(8)                :: u2,v2,w2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    w=1d0-u-v
    u2=2d0*u ; v2=2d0*v ; w2=2d0*w
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ai(1)=w*(-1d0+w2)    !> (i,j,k)=(0,0,2)
    ai(2)=u*(-1d0+u2)    !> (i,j,k)=(2,0,0)
    ai(3)=v*(-1d0+v2)    !> (i,j,k)=(0,2,0)
    ai(4)=u2*w2          !> (i,j,k)=(1,0,1)
    ai(5)=u2*v2          !> (i,j,k)=(1,1,0)
    ai(6)=v2*w2          !> (i,j,k)=(0,1,1)
   !write(*,'("u,v=",2(f12.5,1x),"li=",6(f12.5,1x))')u,v,ai(1:6)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine setT3MeshBasis_P2

  subroutine setT4MeshBasisP1(u,v,w,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Numerotation des sommets
    !> 01 03  04
    !> 02
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! delcaration des variables passees en argument
    real(8), intent(in)    :: u,v,w
    real(8), intent(inout) :: ai(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ai(1)=1d0-u-v-w
    ai(2)=    u
    ai(3)=      v
    ai(4)=        w
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine setT4MeshBasisP1
  
  subroutine setT4MeshBasisP2(u,v,w,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Numerotation des sommets
    !> 01 07 03  08 10  04
    !> 05 06     09
    !> 02
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! delcaration des variables passees en argument
    real(8), intent(in)    :: u,v,w
    real(8), intent(inout) :: ai(:)
    !>
    real(8)                :: x
    real(8)                :: u2,v2,w2,x2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    x=1d0-u-v-w
    u2=2d0*u ; v2=2d0*v ; w2=2d0*w ; x2=2d0*x
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
    ai(01)=(-1d0+x2)*x  !> (i,j,k)=(000)
    ai(02)=(-1d0+u2)*u  !> (i,j,k)=(200)
    ai(03)=(-1d0+v2)*v  !> (i,j,k)=(020) 
    ai(04)=(-1d0+w2)*w  !> (i,j,k)=(002)
    ai(05)=u2      *x2  !> (i,j,k)=(100)
    ai(06)=u2*v2        !> (i,j,k)=(110)
    ai(07)=   v2   *x2  !> (i,j,k)=(010)
    ai(08)=      w2*x2  !> (i,j,k)=(001)
    ai(09)=u2   *w2     !> (i,j,k)=(101)    
    ai(10)=   v2*w2     !> (i,j,k)=(011) 
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine setT4MeshBasisP2  

end module additionnal_Functions


module variablesCommunes
  logical :: visu=.true.
  integer :: commWorld,rankWorld,sizeWorld
  integer :: commLocal,rankLocal,sizeLocal
  integer :: order
  integer :: meshOrder=2
end module variablesCommunes


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
  &                                                    &
  &           stride                                  ,&  ! =ker(calc)
  &           solverType                              ,&
  &           localField                              ,&  !   mySolu
  &           distantField                             )  ! linkSolu
  !---
  use cwipi
  use modDeterminant
  use baseSimplex3D
  
  use additionnal_Functions
  
  use variablesCommunes
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
  integer :: stride
  integer :: solverType
  real(8) ::   localField                            (*)
  real(8) :: distantField                            (*)
  !>
  integer          :: i,j,k,iErr
  integer          :: iNod,nNod,iMod,nMod
  real(8), pointer :: uvw (:,:),rst(:,:),a(:),b(:),c(:),vand(:,:)
  integer          :: iDistantPoint
  integer          :: iBary,iVert
  real(8), pointer :: uvwOut(:,:),lagrange(:,:)
  real(8)          :: lagrangeMesh(1:10)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0 )print'(/">>> userInterpolation rankWorld=",i2)',rankWorld
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( visu .and. rankWorld==0 )then
    print '(/"Mon maillage surfacique de couplage")'
    iVert=0
    do i=1,nLocalVertex
      print '("localCoordinates (",i3,")=",3(f12.5,1x))',i,localCoordinates(iVert+1:iVert+3)
      iVert=iVert+3
    enddo
    do i=1,nLocalElement
      print '("localConnectivity(",i3,")=",*(i3,1x))',i,localConnectivity(localConnectivityIndex(i)+1:localConnectivityIndex(i+1))
    enddo
  endif
  
  if( visu .and. rankWorld==0 )then
    nMod=(order+1)*(order+2)*(order+3)/6
    print '(/"Mon Champ Volumique")'
    j=0
    do iMod=1,nMod
      print '("iMod=",i3," localField       =",4(f12.5,1x),t100,"@rkw",i3)',iMod,localField(j+1:j+stride),rankWorld
      j=j+stride
    enddo
  endif
  
  if( visu .and. rankWorld==0 )then
    print '()'
    iVert=0
    do iDistantPoint=1,nDistantPoint
      print '("iDis=",i3," disPtsCoordinates=",3(f12.5,1x)," inside Cell: ",i3,t100,"@rkw",i3)',&
      & iDistantPoint,disPtsCoordinates(iVert+1:iVert+3),disPtsLocation(iDistantPoint),rankWorld
      iVert=iVert+3
    enddo
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Calcul des coordonnées barycentriques dans le triangle P2
  print '()'
  block
    real(8) :: xyz0(1:3)
    
    allocate(uvwOut(1:2,1:nDistantPoint))
    iBary=0
    do iDistantPoint=1,nDistantPoint
      uvwOut(1:2,iDistantPoint)=distantPointsBarycentricCoordinates(iBary+2:iBary+3) ! <= Attention 2:3
      iBary=iBary+3
    enddo
    
    if( visu .and. rankWorld==0 )then
      print '(/"Coordonnées barycentriques")'
      do iDistantPoint=1,nDistantPoint
        print '("uvwOut=",*(f12.5,1x))',uvwOut(1:2,iDistantPoint)
      enddo
      
      print '(/"Coordonnées calculées")'
      do iDistantPoint=1,nDistantPoint
        
        !> Avec maillage dégradé ordre 1
        nMod=4                               !> TriangleP1
        !> Avec maillage ordre 2
        !nMod=(meshOrder+1)*(meshOrder+2)/2  !> TriangleP2
        
        call setT3MeshBasis_P2(u=uvwOut(1,iDistantPoint),v=uvwOut(2,iDistantPoint),ai=lagrangeMesh)
        
        xyz0(1:3)=0d0
        iVert=0
        do iMod=1,nMod
          xyz0(1:3)=xyz0(1:3)+lagrangeMesh(iMod)*localCoordinates(iVert+1:iVert+3)
          iVert=iVert+3
        enddo
        
        print '("iDis=",i3," xyz0=",3(f12.5,1x),t100,"@rkw",i3)',&
        & iDistantPoint,xyz0(1:3),rankWorld
        
      enddo
    endif
    
    deallocate(uvwOut)

  end block
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
call mpi_barrier(commWorld,iErr)
call cwipi_finalize_f()
call mpi_finalize(iErr)
stop 'A POURSUIVRE'
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Affectation des coordonées barycentriques uvwOut
  allocate(uvwOut(1:3,1:nDistantPoint))
  
  iBary=0
  do iDistantPoint=1,nDistantPoint
    uvwOut(1:3,iDistantPoint)=distantPointsBarycentricCoordinates(iBary+2:iBary+4) ! <= Attention 2:4
    iBary=iBary+4
  enddo
  
  if( visu .and. rankWorld==0 )then
    print '()'
    do iDistantPoint=1,nDistantPoint
      print '("iDis=",i3," uvwOut           =",3(f12.5,1x),t100,"@rkw",i3)',&
      & iDistantPoint,uvwOut(1:3,iDistantPoint),rankWorld
    enddo
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Base fonctionelle d'ordre order
  
  nMod=(order+1)*(order+2)*(order+3)/6
  nNod=size(uvwOut,2)
  
  !> transpose = .true. => lagrange(1:nMod,1:nNod)
  call lagrange3Dv(ord=order,uvwOut=uvwOut,lagrange=lagrange,transpose=.true.)
  
!  !> Points d'interpolation
!  call nodes3D   (ord=order,uvw=uvw,display=.false.)
!  call nodes3Dopt(ord=order,uvw=uvw,display=.false.)
!  !> Calcul de Vand(:,:)
!  call nodes3Duvw2abc(uvw=uvw,a=a,b=b,c=c,display=.false.)
!  call vandermonde3D(ord=order,a=a,b=b,c=c,vand=vand)
!  deallocate(uvw,a,b,c)
!  !> Calcul des polonômes de Lagrange d'ordre order en uvwOut
!  allocate(lagrange(1:nMod,1:nNod))
!  call nodes3Duvw2abc(uvw=uvwOut,a=a,b=b,c=c,display=.false.)
!  call lagrange3Dv(ord=order,vand=vand,a=a,b=b,c=c,lx=lagrange,transpose=.true.)  !> lagrange= Inverse[Transpose[Vand]].Psi[xyzOut] lxOut(nPt,np)
  
  if( visu .and. rankWorld==0 )then
    print '()'
    do iNod=1,nNod
      print '("iNod=",i3," lagrange         =",*(f12.5,1x))',iNod,lagrange(1:nMod,iNod)
    enddo
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Calcul de distantField
  j=0
  do iDistantPoint=1,nDistantPoint
    distantField(j+1:j+stride)=0d0
    k=0
    do iMod=1,nMod
      distantField(j+1:j+stride)= distantField(j+1:j+stride)                          &
      &                          +lagrange(iMod,iDistantPoint)*localField(k+1:k+stride)
      k=k+stride
    enddo
    j=j+stride
  enddo
  
  !> Visu de distantField
  if( visu .and. rankWorld==0 )then
    print '()'
    j=0
    do iDistantPoint=1,nDistantPoint
      print '("iDis=",i3," distantField     =",4(f12.5,1x),t100,"@rkw",i3)',&
      & iDistantPoint,distantField(j+1:j+stride),rankWorld
      j=j+stride
    enddo
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(lagrange)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0 )print'("<<< userInterpolation rankWorld=",i2)',rankWorld
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine userInterpolation


program testf
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use iso_fortran_env
  
  use mpi
  use cwipi
  
  use modDeterminant
  use baseSimplex2D
  use baseSimplex3D
  use table_tet_mesh
  
  use additionnal_Functions
  
  use variablesCommunes
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
    integer :: stride
    integer :: solverType
    real(8) ::   localField                            (*)
    real(8) :: distantField                            (*)
    end subroutine  userInterpolation
  end interface
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  character(5)     :: codeName,codeCoupledName
  
  character(64)    :: meshName
  character(80)    :: fileName
  
  integer          :: iVert,nVert
  real(8), pointer :: vertx(:,:)
  integer          :: iTetra,nTetra
  integer, pointer :: tetra(:,:)
  integer          :: iTrian,nTrian
  integer, pointer :: trian(:,:)
  
  integer          :: j,k
  integer          :: iMod,nMod
  integer          :: iNod,nNod
  integer          :: iCell,nCell
  real(8), pointer :: vertices   (:)
  integer, pointer :: connec     (:)
  integer, pointer :: connecIndex(:)
  integer, pointer :: tetraNodes(:,:)
  
  real(8), pointer :: lagrange(:,:)
  real(8)          :: lagrangeMesh(1:10)
  
  real(8)          :: xyz(1:3)
  integer          :: linkVertSize
  real(8), pointer :: linkVert(:)
  integer          :: notLocatedPoints
  
  integer          :: stride
  real(8), pointer ::   myValues(:)
  real(8), pointer :: linkValues(:)
  
  real(8), pointer :: uvw  (:,:),a(:),b(:),c(:)
  real(8), pointer :: uv   (:,:),rs (:,:)
  real(8), pointer :: vand (:,:)
  
  real(8)          :: node_xy(1:2,1:3) !> Triangle
  
  integer          :: iRank,iErr
  
  real(8)          :: delta,deltaMax
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_init(iErr)
  commWorld=mpi_comm_world
  
  call mpi_comm_rank(commWorld, rankWorld, iErr)
  call mpi_comm_size(commWorld, sizeWorld, iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0) print '(/"START: fortran_surf_TetraP2_PiPj")'
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Initialisation de l'interface de couplage
  select case(rankWorld)
  case(0)
     codeName        = "code1"
     codeCoupledName = "code2"
  case(1)
     codeName        = "code2"
     codeCoupledName = "code1"
  end select
  
  call cwipi_init_f(           &
  &    globalComm=commWorld   ,&
  &    appliName=codeName     ,&
  &    appliComm=commLocal     )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_comm_rank(commLocal,rankLocal,iErr)
  call mpi_comm_size(commLocal,sizeLocal,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  select case(rankWorld)
  case(0) ; order=1 !07
  case(1) ; order=1 !10
  end select
  print '("fortran_surf_PiPj : meshOrder=",i2," Order=",i2,t100,"@rkw",i3)',meshOrder,order,rankWorld
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Create coupling
  
  if( rankWorld==0 )print '(/"Create coupling")'
  
  call cwipi_create_coupling_f(                  &
  &    couplingName="testPiPj"                  ,&
  &    couplingType=cwipi_cpl_parallel_with_part,&
  &    cplAppli=codeCoupledName                 ,&
  &    entitiesDim=2                            ,& !> Nature du couplage
  &    tolerance=1d-1                           ,& !> Tolerance geometrique 1d-1 par defaut
  &    meshT=cwipi_static_mesh                  ,&
  &    solvert=cwipi_solver_cell_vertex         ,&
  &    outputfreq=1                             ,& !> Frequence du post-traitement
  &    outputfmt="Ensight Gold"                 ,&
  &    outputfmtopt="binary"                     )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Create Geometric Mesh
  
  if( rankWorld==0 )print '(/"Create Geometric Mesh (inria mesh format")'
  
  !  04        niveau2
  !
  !  10
  !  08 09     niveau1
  !
  !  03
  !  07 06
  !  01 05 02  niveau0
  
  select case(rankWorld)
  case(0)
    !> Vertices
    nVert=10
    allocate(vertx(1:3,1:nVert))    
    vertx(1:3,01)=[0.00, 0.00, 0.00]
    vertx(1:3,02)=[0.00,-1.00, 0.00]
    vertx(1:3,03)=[1.00, 0.00, 0.00]
    vertx(1:3,04)=[0.00, 0.00, 1.00]
    !>
    vertx(1:3,05)=[0.00,-0.50, 0.00]
    vertx(1:3,06)=[0.50,-0.50, 0.00]
    vertx(1:3,07)=[0.50, 0.00, 0.00]
    vertx(1:3,08)=[0.00, 0.00, 0.50]
    vertx(1:3,09)=[0.00,-0.50, 0.50]
    vertx(1:3,10)=[0.50, 0.00, 0.50]
    !> Tetrahedra
    nTetra=1
    allocate(tetra(1:11,1:nTetra)) !> 10 sommets + 1 marqueur
    tetra(1:11,1)=[01,02,03,04,05,06,07,08,09,10, 1]
    !>  Triangles
    nTrian=4
    allocate(trian(1:7,1:nTrian)) !> 6 sommets + 1 marqueur
    trian(1:7,1)=[02,03,04,06,10,09,1]
    trian(1:7,2)=[01,03,02,07,06,05,1]
    trian(1:7,3)=[01,04,03,08,10,07,3] !> Couplage
    trian(1:7,4)=[01,02,04,05,09,08,1]
  case(1)
    !> Vertices
    nVert=10
    allocate(vertx(1:3,1:nVert))
    vertx(1:3,01)=[0.00, 0.00, 0.00]
    vertx(1:3,02)=[1.00, 0.00, 0.00]
    vertx(1:3,03)=[0.00, 1.00, 0.00]
    !>
    vertx(1:3,04)=[0.00, 0.00, 1.00]
    vertx(1:3,05)=[0.50, 0.00, 0.00]
    vertx(1:3,06)=[0.50, 0.50, 0.00]
    vertx(1:3,07)=[0.00, 0.50, 0.00]
    vertx(1:3,08)=[0.00, 0.00, 0.50]
    vertx(1:3,09)=[0.50, 0.00, 0.50]
    vertx(1:3,10)=[0.00, 0.50, 0.50]
    !> Tetrahedra
    nTetra=1
    allocate(tetra(1:11,1:nTetra)) !> 10 sommets + 1 marquer
    tetra(1:11,1)=[1,2,3,4,5,6,7,8,9,10, 1]
    !>  Triangles
    nTrian=4
    allocate(trian(1:7,1:nTrian)) !> 6 sommets + 1 marquer
    trian(1:7,1)=[02,03,04,06,10,09, 1]
    trian(1:7,2)=[01,03,02,07,06,05, 1]
    trian(1:7,3)=[01,04,03,08,10,07, 1]
    trian(1:7,4)=[01,02,04,05,09,08, 3] !> Couplage
  end select
  
  if( visu )then
    !> Ecriture des maillages au format mesh de l'inria
    if( rankWorld==0)then
      do iRank=0,sizeWorld-1
        print '(/"Writing mesh file: Tetra",i1,".mesh")',iRank
      enddo
    endif
    write(meshName,'("Tetra",i1,".mesh")')rankWorld
    open(unit=100,file=trim(meshName),action='write',status='unknown')
    write(100,'("MeshVersionFormatted 1"/)')
    write(100,'("Dimension 3"/)')
    write(100,'("Vertices")')
    write(100,'(i2)')nVert
    do iVert=1,nVert
      write(100,'(3(e22.15,1x),i2)')vertx(1:3,iVert),0
    enddo
    write(100,'(/"TetrahedraP2")')
    write(100,'(i1)')nTetra
    do iTetra=1,nTetra
      write(100,'(*(i6,1x))')tetra(:,iTetra)
    enddo
    write(100,'(/"TrianglesP2")')
    write(100,'(i1)')nTrian
    do iTrian=1,nTrian
      write(100,'(*(i6,1x))')trian(:,iTrian)
    enddo
    write(100,'(/"End")')
    close(100)
  endif
  
  !> On se couple sur le triangle commun aux deux tetraP2 (y=0) qui va servir de maillage pour le couplage
  !
  !  03
  !  06 05
  !  01 04 02
  !
  
  !> rkw0:Triangle3 <=> rkw1:Triangle4
  select case(rankWorld)
  case(0) ; iTrian=3  !> on se couple sur le triangle 3
  case(1) ; iTrian=4  !> on se couple sur le triangle 4
  end select
  
  if( 1==0 )then
    
    !> On degrade le maillage à l'ordre 1
    nVert=03
    nCell=01
    allocate( vertices   (1:3*nVert)    )  !> sommets
    allocate( connec     (1:3*nCell)    )  !> triangle P1
    allocate( connecIndex(1:nCell+1)    )  !> triangle
    connec(1:3)=[1,2,3]
    connecIndex(1:2)=[0,3]
    
    j=0
    do iNod=1,nVert
      vertices(j+1:j+3)=vertx(1:3,trian(iNod,iTrian))
      j=j+3
    enddo
  
  !> Transmission des maillages à cwipi
  call cwipi_define_mesh_f(     &
  &   couplingName="testPiPj"  ,&
  &   nVertex     =nVert       ,&
  &   nElts       =nCell       ,&
  &   coords      =vertices    ,&
  &   connecIndex =connecIndex ,&
  &   connec      =connec       )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    if( visu )then
      write(meshName,'("Triangle",i1,".mesh")')rankWorld
      open(unit=100,file=trim(meshName),action='write',status='unknown')
      write(100,'("MeshVersionFormatted 1"/)')
      write(100,'("Dimension 3"/)')
      write(100,'("Vertices")')
      write(100,'(i2)')nVert
      j=0
      do iVert=1,nVert
        write(100,'(3(e22.15,1x),i2)')vertices(j+1:j+3),0
        j=j+3
      enddo
      write(100,'(/"Triangles")')
      write(100,'(i1)')nCell
      do iCell=1,nCell
        write(100,'(*(i6,1x))')connec( connecIndex(iCell)+1:connecIndex(iCell+1) ),0
      enddo
      write(100,'(/"End")')
      close(100)
    endif
    
  else

   
  
    
    !> On conserve un maillage d'ordre 2
    nVert=06
    nCell=01
    allocate( vertices   (1:3*nVert)    )  !> sommets
    allocate( connec     (1:6*nCell)    )  !> triangle P2
    allocate( connecIndex(1:nCell+1)    )  !> triangle
    connec(1:6)=[1,2,3,4,5,6]
    connecIndex(1:2)=[0,6]
    
    j=0
    do iNod=1,nVert
      vertices(j+1:j+3)=vertx(1:3,trian(iNod,iTrian))
      j=j+3
    enddo
    
  !> Transmission des maillages à cwipi
  call cwipi_define_ho_mesh_f(  &
  &   couplingName="testPiPj"  ,&
  &   nVertex     =nVert       ,&
  &   nElts       =nCell       ,&
  &   order       =2           ,&  
  &   coords      =vertices    ,&
  &   connecIndex =connecIndex ,&
  &   connec      =connec       )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  if( visu )then
      write(meshName,'("Triangle",i1,".mesh")')rankWorld
      open(unit=100,file=trim(meshName),action='write',status='unknown')
      write(100,'("MeshVersionFormatted 1"/)')
      write(100,'("Dimension 3"/)')
      write(100,'("Vertices")')
      write(100,'(i2)')nVert
      j=0
      do iVert=1,nVert
        write(100,'(3(e22.15,1x),i2)')vertices(j+1:j+3),0
        j=j+3
      enddo
      write(100,'(/"TrianglesP2")')
      write(100,'(i1)')nCell
      do iCell=1,nCell
        write(100,'(*(i6,1x))')connec( connecIndex(iCell)+1:connecIndex(iCell+1) ),0
      enddo
      write(100,'(/"End")')
      close(100)
    endif
    
  endif
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Points de couplage situés Triangle avec le marqueur de peau 1 => face commune aux deux tetras <=
  
  if( rankWorld==0 ) print'(/"Calcul des coordonnees des points de couplage")'
  
  !> Calcul des coordonnees barycentriques sur face trianglulaire
  call nodes3Dopt_2D(ord=order,uvw=uvw,display=.false.) !> ordre du calcul
  
  !> Calculs des coordonnées des points de couplage
  linkVertSize=size(uvw,2)
  allocate(linkVert(1:6*linkVertSize))
  
  nMod=(meshOrder+1)*(meshOrder+2)/2  !> TriangleP2
  nNod=size(uvw,2)
  
  !> transpose = .true. => lagrange(1:nMod,1:nNod)
  !call lagrange2Dv(ord=meshOrder,uvwOut=uvw,lagrange=lagrange,transpose=.true.)
  
  j=0
  do iVert=1,linkVertSize
    !> Fonction
    call setT3MeshBasis_P2(u=uvw(1,iVert),v=uvw(2,iVert),ai=lagrangeMesh)
    linkVert(j+1:j+3)=0d0
    do iMod=1,nMod
      linkVert(j+1:j+3)=linkVert(j+1:j+3)+lagrangeMesh(iMod)*vertx(1:3,trian(iMod,iTrian))
    enddo
    j=j+3
  enddo
  deallocate(uvw)
  
  
  !> Visu des coordonnees barycentriques dans le triangle unité
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then
        print '()'
        j=0
        do iVert=1,linkVertSize
          print '("linkVert(",i2,")=",3(f12.5,1x),t100,"@rkw",i3)',iVert,linkVert(j+1:j+3),rankWorld
          j=j+3
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
  endif
  
  call cwipi_set_points_to_locate_f( &
  &    couplingName="testPiPj"      ,&
  &    nPts  =linkVertSize          ,&
  &    coords=linkVert               )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Localisation
  
  if( rankWorld==0 )print '(/"Localisation")'
  
  call cwipi_locate_f(couplingName="testPiPj")
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Initialisation of myValues
  
  if( rankWorld==0 )print '(/"Initialisation of myValues")'
  
  call nodes3D   (ord=order,uvw=uvw,display=.false.)
  if( visu )then
    call driverTetMesh(ord=order,node_xyz=uvw,tetra_node=tetraNodes)
  endif
  
  call nodes3Dopt(ord=order,uvw=uvw,display=.false.)
  if( visu )then
    call saveTetMesh  (ord=order,node_xyz=uvw,tetra_node=tetraNodes)
    deallocate(tetraNodes)
  endif
  
  stride=4 ; nNod=size(uvw,2)
  allocate(myValues(1:stride*nNod))
  
  nMod=(meshOrder+1)*(meshOrder+2)*(meshOrder+3)/6 !> Tetra P2
  j=0
  do iNod=1,nNod
    call setT4MeshBasisP2(u=uvw(1,iNod),v=uvw(2,iNod),w=uvw(3,iNod),ai=lagrangeMesh)    
    xyz(1:3)=0d0
    do iMod=1,nMod
      xyz(1:3)=xyz(1:3)+lagrangeMesh(iMod)*vertx(1:3,tetra(iMod,1))
    enddo
    myValues(j+1:j+stride)=[xyz(1),xyz(2),xyz(3),real(rankWorld,kind=8)]
    j=j+stride
  enddo
  deallocate(uvw)
  
  if( visu .and. rankWorld==0 )then
    print '("nMod=",i3,2x,"nNod=",i3)',nMod,nNod
    j=0
    do iNod=1,nNod
      print '("iMod=",i3," myValues         =",4(f12.5,1x),t100,"@rkw",i3)',iNod,myValues(j+1:j+stride),rankWorld
      j=j+stride
    enddo
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  allocate(linkValues(1:stride*linkVertSize))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0 )print '(/"cwipi_exchange_f")'
  
  call cwipi_exchange_f(                       &
  &    couplingName=          "testPiPj"      ,&
  &    exchangeName="exch1_"//"testPiPj"      ,&
  &    exchangeDim=stride                     ,&  ! scalar
  &    ptInterpolationFct=userInterpolation   ,&  ! utilisation de la procedure plug
  &                                            &
  &    sendingFieldName="mySolu"              ,&  ! solution calculee localement
  &    sendingField=myValues                  ,&
  &                                            &
  &    receivingFieldName="linkSolu"          ,&
  &    receivingField=linkValues              ,&  ! solution de raccord
  &                                            &
  &    nStep=1                                ,&  ! pas utilisee juste pour visu cwipi
  &    timeValue=0d0                          ,&  ! pas utilisee juste pour visu cwipi
  &    nNotLocatedPoints=notLocatedPoints     ,&
  &    status=iErr                             )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0 )print '(/"Delete coupling")'
  call cwipi_delete_coupling_f("testPiPj")
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0 )print '(/"Controling Coupling Results")'
  
  deltaMax=-1d0
  j=0
  k=0
  do iVert=1,linkVertSize
    delta=norm2(linkVert(j+1:j+3)-linkValues(k+1:k+3)) !+( real(rankWorld,kind=8)-linkValues(k+4) )**2
    if( deltaMax<delta )deltaMax=delta
   !print'("Delta=",e22.15)',delta
    j=j+3
    k=k+4
  enddo
  delta=deltaMax
  
  call mpi_allreduce(delta,deltaMax,1,mpi_real8,mpi_max,commWorld,iErr)
  
  if( rankWorld==0 )then
    print '(/"deltaMax=",e22.15/)',deltaMax
  endif
  
  if( deltaMax<1d-12 )then
    if( rankWorld==0 )print '(/"SUCCESS: fortran_surf_PiPj"/)'
  else
    if( rankWorld==0 )print '(/"FAILED: fortran_surf_PiPj"/)'
    stop
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(vertx,tetra,trian)
  deallocate(linkVert,linkValues)
  deallocate(vertices,connecIndex,connec,myValues)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call cwipi_finalize_f()
  call mpi_finalize(iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
end program testf