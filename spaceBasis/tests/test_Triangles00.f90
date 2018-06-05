
subroutine triangle_00()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use baseSimplex2D
  use modDeterminant
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !
  integer            :: i,j,order,nVert,ad
  real(8), pointer   :: vand(:,:)
  real(8), pointer   :: drVand  (:,:),dsVand  (:,:)
  real(8), pointer   :: drMatrix(:,:),dsMatrix(:,:)
  real(8), pointer   :: mass(:,:)
  real(8), pointer   :: psi (:,:)
  real(8), pointer   :: xyout(:,:),lxout(:,:),drLxout(:,:),dsLxout(:,:)
  real(8), pointer   :: leb(:,:)
  real(8), pointer   :: uvw(:,:),rs(:,:),a(:),b(:)
  real(8), parameter :: eps=1d-15
  real(8), pointer   :: eigv(:)
  character(80)      :: fileName
  real(8)            :: node_xy(1:2,1:3) !> Triangle
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'("Calcul des bases polynômiales 2D (Triangle)")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  order=15
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Points d''interpolation")')
  call nodes2D   (ord=order,uvw=uvw,display=.true.)
  call nodes2Dopt(ord=order,uvw=uvw,display=.true.)
  call nodes2Duv2rs(uv=uvw,rs=rs   ,display=.true.) !> rs(1:2,:)=2d0*uv(1:2,:)-1d0
  call nodes2Drs2ab(rs=rs,a=a,b=b  ,display=.true.) !> a=2 (1+r)/(1-s)-1 && b=s
  
  nVert=size(uvw,2)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(  0<=order .and. order< 10 )write(fileName,'("pointInterpolationP0",i1,".eps")')order
  if( 10<=order .and. order<100 )write(fileName,'("pointInterpolationP" ,i2,".eps")')order
  print '("writing File: ",a)',trim(fileName)
  node_xy(1:2,1)=[-1,-1]
  node_xy(1:2,2)=[ 1,-1]
  node_xy(1:2,3)=[-1, 1]
  
  call trianglePointsPlot(   &
  &    file_name=fileName   ,&
  &    node_xy=node_xy      ,&
  &    node_show=0          ,&
  &    point_num=nVert      ,&
  &    point_xy=rs          ,&
  &    point_show=2          ) !> point_show=2, shows the points and number them
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call vandermonde2D(ord=order,a=a,b=b,vand=vand)
  if( order<10 )then
    call display(title="Vandermonde = [P_j(xi_i)]",mat=vand)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Polynomes d'interpolation (on teste avec les points d'interpolation)
  
  allocate(lxOut(1:(order+1)*(order+2)/2,1:nVert))
  call lagrange2Dv(ord=order,vand=vand,a=a,b=b,lx=lxOut,transpose=.true.) !> true pour affichage
  if( order<10 )then
    call display(title="Test avec l(uvw)",mat=lxOut)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse (\intt ai aj dV)
  call massMatrix(vand=vand,mass=mass)
  if( order<10 )then
    call display(title="Mass Matrix",mat=mass)
   !call mathematica(title="Mass Matrix",mat=mass)
  endif
  print '(/"det(mass)=",e22.15/)',DMGT(eps=eps,n=size(mass,1),A=mass)
  call eigenVectors(mat=mass,w=eigv,display=.true.)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Derivees bases fonctionnelles en (a,b)
  call gradVandermonde2D(ord=order,a=a,b=b,drVand=drVand,dsVand=dsVand)
  call derive1D(vand=vand,dVand=drVand,dMat=drMatrix)
  call derive1D(vand=vand,dVand=dsVand,dMat=dsMatrix)
  
  if( order<10 )then
    call display(title="drVand Matrix",mat=drVand)
    call display(title="dsVand Matrix",mat=dsVand)
    call display(title="drMatrix=drVand.vand^{-1}",mat=drMatrix)
    call display(title="dsMatrix=dsVand.vand^{-1}",mat=dsMatrix)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(a,b,rs)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Points solution
  call readXYout2D(xyout=xyout)
  call nodes2Drs2ab(rs=xyout,a=a,b=b,display=.false.)
  nVert=size(xyOut,2)
  
  call simplex2D(ord=order,a=a,b=b,mode=psi,transpose=.false.)
  
  allocate(  lxOut(1:nVert,1:(order+1)*(order+2)/2))
  call lagrange2Dv(ord=order,vand=vand,a=a,b=b,lx=lxout,transpose=.false.)
  
  allocate(drLxout(1:nVert,1:(order+1)*(order+2)/2))
  call dLagrange1Dv(dMat=drMatrix,lx=lxout,dlx=drLxout,transpose=.false.)
  
  allocate(dsLxout(1:nVert,1:(order+1)*(order+2)/2))
  call dLagrange1Dv(dMat=dsMatrix,lx=lxout,dlx=dsLxout,transpose=.false.)
  
  if( order<6 )then
    call writeSolOut2D(title="simplex2D"   ,solOut=psi    )
    call writeSolOut2D(title="lagrange2D"  ,solOut=lxout  )
    call writeSolOut2D(title="drLagrange2D",solOut=drLxout)
    call writeSolOut2D(title="dsLagrange2D",solOut=dsLxout)
  endif
  
  allocate(leb(1:nVert,1))
  call lebesgue(lx=lxout,l=leb,transpose=.false.)
  call writeSolOut2D(title="lebesgue2D",solOut=leb)
  
  deallocate(psi)
  deallocate(a,b)
  deallocate(xyout)
  deallocate(lxout,drLxout,dsLxout)
  deallocate(leb)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine triangle_00


program test_triangles
  
  !> Test Triangles
  call triangle_00()
  
end program test_triangles