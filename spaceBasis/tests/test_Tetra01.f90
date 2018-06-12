
subroutine tetraMaillageVisuNew()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Cette procedure sert à construire les maillages de visu pour le tetra d'ordre élevé
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !use modDeterminant
  use baseSimplex3D
  use table_tet_mesh
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer            :: ord,iOrd
  real(8), pointer   :: uvw  (:,:)
  integer, pointer   :: tetra(:,:)
  integer            :: iVert,iCell
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ord=15
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do iOrd=1,ord
    if(   0<=iOrd .and. iOrd<  10 )write(*,'("TetraP00",i1)')iOrd
    if(  10<=iOrd .and. iOrd< 100 )write(*,'("TetraP0" ,i2)')iOrd
    if( 100<=iOrd .and. iOrd<1000 )write(*,'("TetraP" ,i3)')iOrd
    
    call nodes3D(ord=iOrd,uvw=uvw,display=.false.)
    call driverTetMesh(ord=iOrd,node_xyz=uvw,tetra_node=tetra) !> c'est plus joli en construisant la connectivité avec les points reguliers
    call nodes3Dopt(ord=iOrd,uvw=uvw,display=.false. )
    call saveTetMesh(ord=iOrd,node_xyz=uvw,tetra_node=tetra)
    deallocate(uvw,tetra)
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Le tetra P2 geometrique
  
  write(*,'("TetraGeoP2.mesh")')
  allocate(uvw(1:3,1:10))
  uvw(1:3,01)=[0d0,0d0,0d0]
  uvw(1:3,02)=[1d0,0d0,0d0]
  uvw(1:3,03)=[0d0,1d0,0d0]
  uvw(1:3,04)=[0d0,0d0,1d0]
  uvw(1:3,05)=[5d-1,0d0 ,0d0]
  uvw(1:3,06)=[5d-1,5d-1,0d0]
  uvw(1:3,07)=[0d0 ,5d-1,0d0]
  uvw(1:3,08)=[0d0 ,0d0 ,5d-1]
  uvw(1:3,09)=[5d-1,0d0 ,5d-1]
  uvw(1:3,10)=[0d0 ,5d-1,5d-1]
  call driverTetMesh(ord=iOrd,node_xyz=uvw,tetra_node=tetra)
  
  open(unit=100,file="TetraGeoP2.mesh",action='write')
  write(100,'("MeshVersionFormatted 2")')
  write(100,'(/"Dimension")')
  write(100,'( "3")')

  write(100,'(/"Vertices")')
  write(100,'( "10")')
  do iVert=1,10
    write(100,'(3(f6.2,1x),1x,i1)')uvw(1:3,iVert),0
  enddo
  
  write(100,'(/"Tetrahedra")')
  write(100,'(i1)')size(tetra,2)
  do iCell=1,size(tetra,2)
    write(100,'(4(i2,1x),1x,i1)')tetra(1:4,iCell),iCell
  enddo
  
  write(100,'(/"End")')
  
  deallocate(uvw,tetra)
  
  close(100)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine tetraMaillageVisuNew


program main
  
  call tetraMaillageVisuNew()
  
end program main
