
subroutine edge_01()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use baseSimplex1D
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer          :: i,j,order,nVert,ad,iOrd
  real(8), pointer :: vand(:,:)
  real(8), pointer :: drVand  (:,:),dsVand  (:,:)
  real(8), pointer :: drMatrix(:,:),dsMatrix(:,:)
  real(8), pointer :: mass(:,:)
  real(8), pointer :: mode(:,:)
  real(8), pointer :: xyout(:,:),lxout(:,:),drLxout(:,:),dsLxout(:,:)
  real(8), pointer :: r(:)
  real(8), pointer :: eigv(:)
  real(8), parameter :: eps=1d-15
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  order=10
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  do iOrd=1,order
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call nodes1D   (ord=iOrd,uvw=r,display=.false.)
   !call nodes1Dopt(ord=iOrd,uvw=r,display=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call vandermonde1D(ord=iOrd,a=r,vand=vand)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Matrice de masse
    call massMatrix(vand=vand,mass=mass)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call eigenVectors(mat=mass,w=eigv,display=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '("iOrd=",i2,2x,"cond(massEqu)=",e22.15,2x,"det(vand)=",e22.15)'&
    & ,iOrd,maxval(eigv)/minval(eigv),DMGT(eps=eps,n=iOrd+1,A=vand)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(r,mass,vand,eigv)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call nodes1Dopt(ord=iOrd,uvw=r,display=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call vandermonde1D(ord=iOrd,a=r,vand=vand)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Matrice de masse
    call massMatrix(vand=vand,mass=mass)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call eigenVectors(mat=mass,w=eigv,display=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '("iOrd=",i2,2x,"cond(massOpt)=",e22.15,2x,"det(vand)=",e22.15)'&
    & ,iOrd,maxval(eigv)/minval(eigv),DMGT(eps=eps,n=iOrd+1,A=vand)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(r,mass,vand,eigv)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    print '()'
    
  enddo
  
  return
end subroutine edge_01


program main
  
  !> Test segments
  call edge_01()
  
end program main