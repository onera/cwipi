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

module additionnal_Functions

contains
  
  subroutine triangleP2UV_to_TetraP2UVW(iTrian,uv,uvw)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Cette procedure retourne les coordonnées barycentriques dans le tetraP2 des sommets de la face iTrian du tetraP2
    !> entree uv (1,2,;) triangleP2
    !> sortie uvw(1,3,;) tetraP2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer         , intent(in)  :: iTrian
    real(8), pointer, intent(in)  :: uv (:,:)
    real(8), pointer, intent(out) :: uvw(:,:)
    !>
    integer                       :: iVert
    integer                       :: iNod,jNod
    integer                       :: nodes(1:6)
    real(8)                       :: TetraP2(1:3,1:10)
    real(8)                       :: TrianP2(1:3,1:06)
    real(8)                       :: ai(1:6),u,v
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> uvw du TetraP2
    TetraP2(1:3,01)=[0d+0,0d+0,0d+0]
    TetraP2(1:3,02)=[1d+0,0d+0,0d+0]
    TetraP2(1:3,03)=[0d+0,1d+0,0d+0]
    TetraP2(1:3,04)=[0d+0,0d+0,1d+0]
    !>
    TetraP2(1:3,05)=[5d-1,0d+0,0d+0]
    TetraP2(1:3,06)=[5d-1,5d-1,0d+0]
    TetraP2(1:3,07)=[0d+0,5d-1,0d+0]
    TetraP2(1:3,08)=[0d+0,0d+0,5d-1]
    TetraP2(1:3,09)=[5d-1,0d+0,5d-1]
    TetraP2(1:3,10)=[0d+0,5d-1,5d-1]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> slection des noeuds
    select case(iTrian)
    case(1) ; nodes(1:6)=[02,03,04, 06,10,09]
    case(2) ; nodes(1:6)=[01,03,02, 07,06,05]
    case(3) ; nodes(1:6)=[01,04,03, 08,10,07]
    case(4) ; nodes(1:6)=[01,02,04, 05,09,08]
    case default ; stop "stop @ triangleP2UV_to_TetraP2UVW"
    end select
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> uvw du TriangleP2 face iTrian
    do iNod=1,6
      jNod=nodes(iNod)
      TrianP2(1:3,iNod)=TetraP2(1:3,jNod)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Calcul de uvw dans le TetraP2 correspond a uv dans le TriangleP2
    allocate(uvw(1:3,size(uv,2)))
    do iVert=1,size(uv,2)
      u=uv(1,iVert)
      v=uv(2,iVert)
      call setT3MeshBasis_P2(u=u,v=v,ai=ai)
      uvw(1:3,iVert)=0d0
      do iNod=1,6
        uvw(1:3,iVert)=uvw(1:3,iVert)+ai(iNod)*TrianP2(1:3,iNod)
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine triangleP2UV_to_TetraP2UVW

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
  integer :: compOrder
  integer :: meshOrder
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
  integer :: uvw_size
  real(8) :: dist_uvw                                (*)
  integer :: stride
  integer :: solverType
  real(8) ::   localField                            (*)
  real(8) :: distantField                            (*)
  !>
  integer          :: i,j,k,iRank,iErr
  integer          :: iNod,nNod,iMod,nMod
  real(8), pointer :: uvw (:,:),rst(:,:),a(:),b(:),c(:),vand(:,:)
  integer          :: iDistantPoint
  integer          :: iTrian,iVert,iBary
  real(8), pointer :: uv(:,:),uvwOut(:,:),lagrange(:,:)
  real(8)          :: lagrangeMesh(1:10)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,">>> userInterpolation")'  
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Control localCoordinates localConnectivity
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Control: localCoordinates, localConnectivity")'  
  call mpi_barrier(commWorld,iErr)
  if( 0==1 )then
 !if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        iVert=0
        do i=1,nLocalVertex
          print '(6x,"localCoordinates (",i2,")=",3(e22.15,1x))',i,localCoordinates(iVert+1:iVert+3)
          iVert=iVert+3
        enddo
        do i=1,nLocalElement
          print '(6x,"localConnectivity(",i2,")=",*(i3,1x))',i,localConnectivity(localConnectivityIndex(i)+1:localConnectivityIndex(i+1))
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
  endif
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Control localField
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Control: localField")'
  call mpi_barrier(commWorld,iErr)
 !if( visu )then
  if( 0==1 )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        nMod=(compOrder+1)*(compOrder+2)*(compOrder+3)/6
        j=0
        do iMod=1,nMod
          print '(6x,"localField(",i3,")=",*(e22.15,1x))',iMod,localField(j+1:j+stride)
          j=j+stride
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Control disPtsCoordinates,disPtsLocation  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Control: disPtsCoordinates,disPtsLocation")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then    
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        iVert=0
        do iDistantPoint=1,nDistantPoint
          if(  iRank==1 .and. iDistantPoint==29 )print '(6x,"disPtsCoordinates(",i3,")=",3(e22.15,1x)," inside Cell: ",i3)',&
          & iDistantPoint,disPtsCoordinates(iVert+1:iVert+3),disPtsLocation(iDistantPoint)
          iVert=iVert+3
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> dist_uvw(:) -> uv(1:2,:) (sans dupliquer le bloc mémoire)
  !> la commmande c_f_pointer a besoin du use iso_c_binding, only: c_loc,c_f_pointer
  
  call c_f_pointer(cptr=c_loc(dist_uvw), fptr=uv, shape=[2,nDistantPoint])
  
  
  !if( rankWorld==1 )then
  !  print '(6x,"uv (1:2,11)=",2(e22.15,1x), " -> "$)',uv(1:2,11)  
  !  uv(1:2,11)=[0.966503472726928d-01,0.248894314819442d+00]
  !  print '(2(e22.15,1x))',uv(1:2,11)
  !  print '(6x,"uv (1:2,29)=",2(e22.15,1x), " -> "$)',uv(1:2,29)  
  !  uv(1:2,29)=[0.654455337907865d+00,0.248894314819442d+00]
  !  print '(2(e22.15,1x))',uv(1:2,29)
  !endif

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Control dist_uvw,uv
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Control: dist_uvw,disPtsLocation")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then    
        print '(/3x,"meshOrder=",i1," compOrder=",i2," stride_uvw=",i2,t120,"@rkw",i3)',meshOrder,compOrder,uvw_size,rankWorld
        call mpi_barrier(commWorld,iErr)
        iVert=0
        do iDistantPoint=1,nDistantPoint
         !print '(6x,"dist_uvw(",i3,")=",*(e22.15,1x))',iDistantPoint,dist_uvw(iVert+1:iVert+uvw_size)
          if( iRank==1 .and. iDistantPoint==29 )print '(6x,"uv (1:2,",i3,")=",2(e22.15,1x))',iDistantPoint,uv(1:2,iDistantPoint)
          iVert=iVert+uvw_size
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> uv TriangleP2 -> uvw tetraP2
  !> On passe de la face TriangleP2 au TetraP2 afin de pouvoir calculer des gradients du champs
  !> cette étape n'est pas nécessaire si aucun gradient est calculé
  
  allocate(uvw(1:3,1:nDistantPoint))
  
  !> rkw0:Triangle3 <=> rkw1:Triangle4
  select case(rankWorld)
  case(0) ; iTrian=3  !> on se couple sur le triangle 3
  case(1) ; iTrian=4  !> on se couple sur le triangle 4
  end select
  
  call triangleP2UV_to_TetraP2UVW(iTrian=iTrian,uv=uv,uvw=uvw)
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Calcul: uv TriangleP2 -> uvw tetraP2")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then    
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        call mpi_barrier(commWorld,iErr)
        do iDistantPoint=1,nDistantPoint
          if( iRank==1 .and. iDistantPoint==29 )print '(6x,"uvw(",i3,")=",3(e22.15,1x))',iDistantPoint,uvw(1:3,iDistantPoint)
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)        
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Base fonctionelle d'ordre compOrder
  
  nMod=(compOrder+1)*(compOrder+2)*(compOrder+3)/6  !> Tetra
  
  !> transpose = .true. => lagrange(1:nMod,1:nDistantPoint)
  call lagrange3Dv(ord=compOrder,uvwOut=uvw,lagrange=lagrange,transpose=.true.)
  
!  !> Points d'interpolation
!  call nodes3D   (ord=compOrder,uvw=uvw,display=.false.)
!  call nodes3Dopt(ord=compOrder,uvw=uvw,display=.false.)
!  !> Calcul de vand(:,:)
!  call nodes3Duvw2abc(uvw=uvw,a=a,b=b,c=c,display=.false.)
!  call vandermonde3D(ord=compOrder,a=a,b=b,c=c,vand=vand)
!  deallocate(uvw,a,b,c)
!  !> Calcul des polonômes de Lagrange d'ordre compOrder en uvwOut
!  allocate(lagrange(1:nMod,1:nNod))
!  call nodes3Duvw2abc(uvw=uvwOut,a=a,b=b,c=c,display=.false.)
!  call lagrange3Dv(ord=compOrder,vand=vand,a=a,b=b,c=c,lx=lagrange,transpose=.true.)  !> lagrange= Inverse[Transpose[Vand]].Psi[xyzOut] lagrange(nPt,np)  
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Calcul: Bases de lagrange")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then    
        !print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        call mpi_barrier(commWorld,iErr)
        do iDistantPoint=1,nDistantPoint
          !print '(6x,"lagrange(",i3,")=",*(e22.15,1x))',iDistantPoint,lagrange(1:nMod,iDistantPoint)
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)  
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
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Calcul: distantField")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then    
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        call mpi_barrier(commWorld,iErr)
        j=0
        do iDistantPoint=1,nDistantPoint
          if( iRank==1 .and. iDistantPoint==29 )print '(6x,"distantField(",i3,")=",4(e22.15,1x),t120,"@rkw",i3)',&
          & iDistantPoint,distantField(j+1:j+stride),rankWorld
          j=j+stride
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)    
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(lagrange)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_barrier(commWorld,iErr)
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'("<<< userInterpolation rankWorld=",i2)',rankWorld
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine userInterpolation


program fortran_surf_TriaP2_PiPj
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
  
  logical          :: droit=.false.
 !logical          :: droit=.true.
  integer          :: iVert,nVert
  real(8), pointer :: vertx(:)
  integer, pointer :: cells(:),cellsIdx(:)
  integer, pointer :: cel(:)  !> ensemble des noeuds pour iCell : cel=>cells(cellsIdx(iCell)+1:cellsIdx(iCell+1))

  
  character(80)      :: key
  integer, parameter :: iFile=100
  logical            :: test
  
  
  integer          :: iTetra,nTetra
  integer, pointer :: tetra(:,:)
  integer          :: iTrian,nTrian
  integer, pointer :: trian(:,:)
  
  integer          :: i,j,k
  integer          :: iMod,nMod
  integer          :: iNod,nNod
  integer          :: iCell,nCell
  integer, pointer :: ijk     (:)
  integer, pointer :: tetraNodes(:,:)
  
  real(8), pointer :: lagrange(:,:)
  real(8), pointer :: lagrangeMesh(:)
  
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
  if( rankWorld==0) print '(/"START: fortran_surf_TriaPi_PiPj")'
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
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_comm_rank(commLocal,rankLocal,iErr)
  call mpi_comm_size(commLocal,sizeLocal,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  select case(rankWorld)
  case(0) ; compOrder=07 !07 !07
  case(1) ; compOrder=10 !07 !10
  end select
  
  call mpi_barrier(commWorld,iErr)
  print '("fortran_surf_TriaPi_PiPj: meshOrder=",i2," compOrder=",i2,t130,"@rkw",i3)',meshOrder,compOrder,rankWorld
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Reading Geometric Mesh
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Lecture des maillages géométriques (inria mesh format)")'
  call mpi_barrier(commWorld,iErr)
  
  select case(rankWorld)
  case(0) ; meshName="~/Desktop/sphere01.mesh"
  case(1) ; meshName="~/Desktop/sphere02.mesh"
  end select
  
  print '("meshName: ",a,t130,"@rkw",i3)',trim(meshName),rankWorld
  
  inquire(file=trim(meshName),exist=test) ; print '("exist0:  ",l)',test
  open(unit=iFile,file=trim(meshName))
  inquire(unit=iFile,opened=test)         ; print '("opened0: ",l)',test
  do i=1,25
!    key="coucoucoucoucoucou"
    write(unit=iFile,fmt='(i10)')2*i 
  enddo
  close(unit=iFile)
  
  
  print '("meshName: ",a,t130,"@rkw",i3)',trim(meshName),rankWorld
  inquire(file=trim(meshName),exist=test) ; print '("exist1:  ",l)',test
  inquire(unit=iFile,opened=test)         ; print '("opened1: ",l)',test
  
 !open(unit=iFile,file=trim(meshName),status='old',action='read')  
  open(unit=iFile,file=trim(meshName))
  
  inquire(unit=iFile,opened=test)         ; print '("opened2: ",l)',test
  
  do i=1,10
   !write(unit=iFile,fmt='(i10)')2*i+1 
    read(unit=iFile,fmt='(i10)')j 
    !read(iFile,fmt=*)key
    !print '(i2," key: ",a)',i,key
  enddo
    
  
!  endOfFile=0
!  do while(endOfFile==0)
!    print '("endOfFile=",i10)',endOfFile
!    read(unit=iFile,fmt=*,iostat=endOfFile)key
!    print '("key=",a)',trim(key)
!  enddo
  close(unit=iFile)
  inquire(unit=iFile,opened=test)
  print '("opened2: ",l)',test
  
  
  call cwipi_finalize_f()
  call mpi_finalize(iErr)
  stop

  
  open(unit=iFile,file=trim(meshName),status='old',action='read', iostat=iErr)  
  print '("iostat=",i10)',iErr
  
  do i=1,10
    read(unit=iFile,fmt=*)key
    print '(i2," key=",a)',i,trim(key)
  enddo
  
  close(unit=iFile)
  
  call cwipi_finalize_f()
  call mpi_finalize(iErr)
  stop
  
  
  
  lecture: do
    read(unit=iFile,fmt=*)key ! print '("key: ",a,t130,"@rkw",i3)',trim(key),rankWorld
    
   !print '("key: ",a,t130,"@rkw",i3)',trim(key),rankWorld
    print '("key=",a)',trim(key)
    
    select case(trim(key))
    case("Vertices","vertices")
      
      read(iFile,*)nVert
      allocate( vertx(1:3*nVert) )  !> sommets
      do iVert=1,nVert
        read(iFile,*)vertx(3*(iVert-1)+1:3*iVert)
      enddo
      
    case("Quadrilaterals","QuadrilateralsQ2","QuadrilateralsQ3","Triangles","TrianglesP2","TrianglesP3")
      
      select case(trim(key))
      case("Quadrilaterals"  ) ; meshOrder=1 ; nMod=(meshOrder+1)*(meshOrder+1)
      case("QuadrilateralsQ2") ; meshOrder=2 ; nMod=(meshOrder+1)*(meshOrder+1)
      case("QuadrilateralsQ3") ; meshOrder=2 ; nMod=(meshOrder+1)*(meshOrder+1)
      case("Triangles"       ) ; meshOrder=1 ; nMod=(meshOrder+1)*(meshOrder+2)/2
      case("TrianglesP2"     ) ; meshOrder=2 ; nMod=(meshOrder+1)*(meshOrder+2)/2
      case("TrianglesP3"     ) ; meshOrder=3 ; nMod=(meshOrder+1)*(meshOrder+2)/2
      end select
      
      read(iFile,*)nCell
      allocate( cells   (1:nMod*nCell) )
      allocate( cellsIdx(1:nCell+1   ) )
      cellsIdx(1)=0
      do iCell=1,nCell
        cellsIdx(iCell)=cellsIdx(iCell-1)+nMod
        read(iFile,*)cells(cellsIdx(iCell)+1:cellsIdx(iCell+1))
      enddo
      
    case("End","end") ; exit lecture
    end select
  enddo lecture
  
  close(unit=iFile)
  
  call cwipi_finalize_f()
  call mpi_finalize(iErr)
  stop
  
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Fin Lecture des maillages géométriques (inria mesh format)")'
  call mpi_barrier(commWorld,iErr)

stop
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Initialisation of myValues
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 ) print'(/"Initialisation de myValues (x,y,z,rkw)")'
  call mpi_barrier(commWorld,iErr)
  
  stride=4    !> x,y,z,real(rankWorld,kind=8)
  allocate(myValues(1:nCell*stride*nMod))
  i=0 ; j=0
  do iCell=1,nCell
    cel=>cells(cellsIdx(iCell)+1:cellsIdx(iCell+1))
    do iMod=1,nMod
      i=3*cel(iMod)
      myValues(j+1:j+stride)=[vertx(i+1),vertx(i+2),vertx(i+3),real(rankWorld,kind=8)]
    enddo
    i=i+3
    j=j+stride
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Points de couplage (on utilise compOrder>meshOrder pour ajouter des points)
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 ) print'(/"Calcul des coordonnees de linkVert")'
  call mpi_barrier(commWorld,iErr)
  
  !> Calcul des coordonnees barycentriques optimisées sur face triangulaire compOrder
  call nodes3Dopt_2D(ord=compOrder,uvw=uvw,display=.false.) !> ordre du calcul
    
  !> Calculs des coordonnees geometriques correspondantes aux coordonnees barycentriques
  linkVertSize=nCell*size(uvw,2)       !> nombre total de point de couplages
  allocate(linkVert(1:3*linkVertSize)) !> 3 coordonnées par point de couplage
  
  nNod=size(uvw,2)                     !> Triangle compOrder
  
  allocate(lagrangeMesh(1:nMod))
  
  j=0
  do iCell=1,nCell
    do iNod=1,nNod
      !> Fonctions
      call setT3MeshBasis_P2(u=uvw(1,iNod),v=uvw(2,iNod),ai=lagrangeMesh) !> base Triangle Geometrique P2
      !> linkVert
      linkVert(j+1:j+3)=0d0
      do iMod=1,nMod
        i=3*trian(iMod,iCell)
        linkVert(j+1:j+3)=linkVert(j+1:j+3)+lagrangeMesh(iMod)*vertx(i+1:i+3)
      enddo
      j=j+3      
    enddo
  enddo
  
  deallocate(lagrangeMesh)
  deallocate(uvw)
  
  !> Visu des coordonnees de couplage
  if( visu )then
    do iRank=0,sizeWorld-1
      call mpi_barrier(commWorld,iErr)
      if( iRank==rankWorld )then
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        j=0
        do iCell=1,nCell
          if( iCell==1 )then
            do iNod=1,nNod
              print '(6x,"linkVert(",i2,")=",3(e22.15,1x))',iNod,linkVert(j+1:j+3)
              j=j+3      
            enddo
          endif
        enddo      
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Allocation de linkValues")'
  call mpi_barrier(commWorld,iErr)
  
  allocate(linkValues(1:stride*linkVertSize))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Create coupling
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Creation du couplage testPiPj")'
  call mpi_barrier(commWorld,iErr)
  
  call cwipi_create_coupling_f(                  &
  &    couplingName="testPiPj"                  ,&
  &    couplingType=cwipi_cpl_parallel_with_part,&
  &    cplAppli=codeCoupledName                 ,&
  &    entitiesDim=2                            ,& !> Nature du couplage surfacique
  &    tolerance=1d-1                           ,& !> Tolerance geometrique 1d-1 par defaut
  &    meshT=cwipi_static_mesh                  ,&
  &    solvert=cwipi_solver_cell_vertex         ,&
  &    outputfreq=-1                            ,& !> Frequence du post-traitement
  &    outputfmt="Ensight Gold"                 ,&
  &    outputfmtopt="binary"                     )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> On se couple sur les sphères maillées
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Transmission à cwipi du maillage")'
  call mpi_barrier(commWorld,iErr)
  
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
  !> Definition du Triangle géométrique P2
  
  ! 3
  ! 6 5
  ! 1 4 2
  allocate( ijk(1:12) )  !> sommets
  ijk( 1: 2)=[0,0] !> 1
  ijk( 3: 4)=[2,0] !> 2
  ijk( 5: 6)=[0,2] !> 3
  ijk( 7: 8)=[1,0] !> 4
  ijk( 9:10)=[1,1] !> 5
  ijk(11:12)=[0,1] !> 6
  
  call cwipi_ho_ordering_from_IJK_set_f( & !> NEW Cwipi
  &   couplingName ="testPiPj"          ,&
  &   tElt         = CWIPI_FACE_TRIAHO  ,&
  &   nNodes       = 6                  ,&
  &   IJK          = ijk                 )
  
  deallocate (ijk)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Transmission à cwipi des coordonnees de couplage
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Transmission à cwipi de linkVert")'
  call mpi_barrier(commWorld,iErr)
  
  call cwipi_set_points_to_locate_f( &
  &    couplingName="testPiPj"      ,&
  &    nPts  =linkVertSize          ,&
  &    coords=linkVert               )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Localisation par cwipi des coordonnees de couplage
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Localisation")'
  call mpi_barrier(commWorld,iErr)
  
  call cwipi_locate_f(couplingName="testPiPj")
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"cwipi_exchange_f")'
  call mpi_barrier(commWorld,iErr)
  
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
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Delete coupling")'
  call mpi_barrier(commWorld,iErr)
  call cwipi_delete_coupling_f("testPiPj")
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Controling Coupling Results")'
  call mpi_barrier(commWorld,iErr)
  
  deltaMax=-1d0
  j=0
  k=0
  do iRank=0,sizeWorld-1
    if( iRank==rankWorld )then
      print '(/3x,"compOrder =",i2," controling linkValues - linkVertSize=",i2,t120,"@rkw",i3)',compOrder,linkVertSize,rankWorld
      do iVert=1,linkVertSize
        delta=norm2(linkVert(j+1:j+3)-linkValues(k+1:k+3)) !+( real(rankWorld,kind=8)-linkValues(k+4) )**2
        if( deltaMax<delta )deltaMax=delta
        if( 1d-08<=delta )then
          print'(6x,"iVert=",i2,1x,"linkVert=",3(e22.15,1x),"linkValues=",3(e22.15,1x),"Delta=",e22.15)',&
          & ivert,linkVert(j+1:j+3),linkValues(k+1:k+3),delta
        endif
        j=j+3
        k=k+4
      enddo
      delta=deltaMax
    endif
    call mpi_barrier(commWorld,iErr)
  enddo
  
  call mpi_allreduce(delta,deltaMax,1,mpi_real8,mpi_max,commWorld,iErr)
  
  if( rankWorld==0 )then
    print '(/"deltaMax=",e22.15/)',deltaMax
  endif
  
  if( deltaMax<1d-08 )then
    if( rankWorld==0 )print '(/"SUCCESS: fortran_surf_TriaPi_PiPj"/)'
  else
    if( rankWorld==0 )print '(/"FAILED: fortran_surf_TriaPi_PiPj"/)'
    stop
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(vertx,tetra,trian)
  deallocate(linkVert,linkValues)
  deallocate(vertx,cellsIdx,cells,myValues)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call cwipi_finalize_f()
  call mpi_finalize(iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
end program fortran_surf_TriaP2_PiPj
