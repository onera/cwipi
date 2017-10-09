

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
  integer            :: ord,iOrd,ad,nVert
  real(8), pointer   :: uvw  (:,:)
  integer, pointer   :: tetra(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ord=15
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do iOrd=1,ord
    if(  0<=iOrd .and. iOrd< 10 )write(fileName,'("TetraP0",i1")')iOrd
    if( 10<=iOrd .and. iOrd<100 )write(fileName,'("TetraP" ,i2")')iOrd
    
    call nodes3D(ord=iOrd,uvw=uvw,display=.true.)
    call driverTetMesh(ord=iOrd,node_xyz=uvw,tetra_node=tetra) !> c'est plus joli en construisant la connectivité avec les points reguliers
    call nodes3Dopt(ord=iOrd,uvw=uvw,display=.true. )
    call saveTetMesh(ord=iOrd,node_xyz=uvw,tetra_node=tetra)
    deallocate(uvw,tetra)
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine tetraMaillageVisuNew


program main
  
  !> Test Tetra
  call tetraMaillageVisuNew()
  
end program main