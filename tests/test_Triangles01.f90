

subroutine triangle_01()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use baseSimplex2D
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !
  integer          :: i,j,order,nVert,ad,iOrd
  real(8), pointer :: vand(:,:)
  real(8), pointer :: drVand  (:,:),dsVand  (:,:)
  real(8), pointer :: drMatrix(:,:),dsMatrix(:,:)
  real(8), pointer :: mass(:,:)
  real(8), pointer :: mode(:,:)
  real(8), pointer :: xyout(:,:),lxout(:,:),drLxOut(:,:),dsLxOut(:,:)
  real(8), pointer :: uvw(:,:),rs(:,:),a(:),b(:)
  real(8), pointer :: eigv(:)  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  order=10
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  do iOrd=1,order
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call nodes2D   (ord=iOrd,uvw=uvw,display=.false.)
   !call nodes2Dopt(ord=iOrd,uvw=uvw,display=.false.)
    call nodes2Duv2rs(uv=uvw,rs=rs  ,display=.false.) !> rs(1:2,:)=2d0*uv(1:2,:)-1d0
    call nodes2Drs2ab(rs=rs,a=a,b=b ,display=.false.) !> a=2 (1+r)/(1-s)-1 && b=s
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call vandermonde2D(ord=iOrd,a=a,b=b,vand=vand)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Matrice de masse
    call massMatrix(vand=vand,mass=mass)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call eigenVectors(mat=mass,w=eigv,display=.false.)
    print '("iOrd=",i2,2x,"cond(massEqu)=",e22.15)',iOrd,maxval(eigv)/minval(eigv)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw,a,b,rs,mass,vand,eigv)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call nodes2D   (ord=iOrd,uvw=uvw,display=.false.)
    call nodes2Dopt(ord=iOrd,uvw=uvw,display=.false.)
    call nodes2Duv2rs(uv=uvw,rs=rs  ,display=.false.)  !> rs(1:2,:)=2d0*uv(1:2,:)-1d0
    call nodes2Drs2ab(rs=rs,a=a,b=b ,display=.false.)  !> a=2 (1+r)/(1-s)-1 && b=s
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call vandermonde2D(ord=iOrd,a=a,b=b,vand=vand)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Matrice de masse
    call massMatrix(vand=vand,mass=mass)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call eigenVectors(mat=mass,w=eigv,display=.false.)
    print '("iOrd=",i2,2x,"cond(massOpt)=",e22.15)',iOrd,maxval(eigv)/minval(eigv)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw,a,b,rs,mass,vand,eigv)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    print '()'
    
  enddo
  
  return
end subroutine triangle_01



program test_triangles
  
  !> Test Triangles
 !call triangle_00()
  call triangle_01()
    
end program test_triangles