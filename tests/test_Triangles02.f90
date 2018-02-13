

subroutine triangle_02()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use baseSimplex2D
  use baseSimplex3D
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !
  integer          :: iOrd,order
  real(8), pointer :: uv(:,:)
  character(80)    :: fileName
  real(8)          :: node_xy(1:2,1:3) !> Triangle unite
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  order=10
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  node_xy(1:2,1)=[0,0]
  node_xy(1:2,2)=[1,0]
  node_xy(1:2,3)=[0,1]
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '(/"Triangles en 2D")'
  do iOrd=1,order
    call nodes2D   (ord=iOrd,uvw=uv,display=.false.)
    call nodes2Dopt(ord=iOrd,uvw=uv,display=.false.)
    
    if(  0<=iOrd .and. iOrd< 10 )write(fileName,'("triangle2D_P0",i1,".eps")')iOrd
    if( 10<=iOrd .and. iOrd<100 )write(fileName,'("triangle2D_P" ,i2,".eps")')iOrd
    print '(3x,"writing File: ",a)',trim(fileName)
    
    call trianglePointsPlot(   &
    &    file_name=fileName   ,&
    &    node_xy=node_xy      ,&
    &    node_show=0          ,&
    &    point_num=size(uv,2) ,&
    &    point_xy=uv          ,&
    &    point_show=2          ) !> point_show=2, shows the points and number them
    
    deallocate(uv)
    
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '(/"Triangles en 3D (pas identiques aux Triangles en 2D)")'
  do iOrd=1,order
    call nodes3Dopt_2D(ord=order,uvw=uv,display=.false.)
    
    if(  0<=iOrd .and. iOrd< 10 )write(fileName,'("triangle3D_P0",i1,".eps")')iOrd
    if( 10<=iOrd .and. iOrd<100 )write(fileName,'("triangle3D_P" ,i2,".eps")')iOrd
    print '(3x,"writing File: ",a)',trim(fileName)
    
    call trianglePointsPlot(   &
    &    file_name=fileName   ,&
    &    node_xy=node_xy      ,&
    &    node_show=0          ,&
    &    point_num=size(uv,2) ,&
    &    point_xy=uv          ,&
    &    point_show=2          ) !> point_show=2, shows the points and number them
    
    deallocate(uv)
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine triangle_02



program test_triangles
  
  !> Test Triangles
 !call triangle_00()
  call triangle_02()
    
end program test_triangles