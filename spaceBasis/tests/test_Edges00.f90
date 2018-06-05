subroutine edge_00()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use baseSimplex1D
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !
  integer            :: i,j,order,nVert,ad
  real(8), pointer   :: vand(:,:),dVand(:,:),jf(:),dr(:,:)
  real(8), pointer   :: mass(:,:)
  real(8), pointer   :: r(:)=>null()
  real(8), pointer   :: xout(:)
  real(8), pointer   :: lxout(:,:),dlxout(:,:)
  real(8), pointer   :: leb(:,:)
  real(8), parameter :: eps=1d-15
  real(8), pointer   :: eigv(:)  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'("Calcul des bases polynômiales Lagrange 1D")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  order=10
  !write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> points d'interpolation dans [-1,1]
 !call nodes1D   (ord=order,uvw=r,display=.true.)
  call nodes1Dopt(ord=order,uvw=r,display=.true.)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de Vandermonde 1D
  call vandermonde1D(ord=order,a=r,vand=vand)
  call display(title="Vandermonde = [P_j(xi_i)]",mat=vand)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse
  call massMatrix(vand=vand,mass=mass)
  call display(title="Mass=Inverse[Vand.Transpose[Vand]]",mat=mass)
  !print '("det(mass)=",f22.15)',deter(n=size(mass,1),A=mass)
  print '(/"det(mass)=",e22.15/)',DMGT(eps=eps,n=size(mass,1),A=mass)
  call eigenVectors(mat=mass,w=eigv,display=.true.)
  deallocate(eigv)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Points solution
  nVert=05 ; allocate(xout(1:nVert)) ; xout(1:nVert)=[( (-1d0+2d0*real(i-1,kind=8)/real(nVert-1,kind=8)), i=1,nVert)]
 !nVert=01 ; allocate(xout(1:nVert)) ; xout(1:nVert)=-5d-1
  print '(/"size(xOut)=",i3)',size(xout)
  print '("xOut(",i2,")=",f12.5)',(i,xOut(i),i=1,nVert)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Polynomes d'interpolation
  allocate(lxOut(1:order+1,1:nVert))
  call lagrange1Dv(ord=order,vand=vand,x=xOut,lx=lxOut,transpose=.true.)
  call display(title="lxOut",mat=lxOut)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Dérivée Polynomes d'interpolation
  call gradVandermonde1D(ord=order,a=r,dVand=dVand)
  call display(title="drVand Matrix",mat=dVand)
  
  call derive1D(vand=vand,dVand=dVand,dMat=dr)
  call display(title="Dr Matrix",mat=dr)
  
  allocate(dlXout(1:order+1,1:nVert))
  call dLagrange1Dv(dMat=dr,lx=lxOut,dlx=dlXout,transpose=.true.) ! <= true pour affichage
  call display(title="d lx=Dr^T lx",mat=dlXout)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  allocate(leb(1,1:nVert))
  call lebesgue(lx=lxOut,l=leb,transpose=.true.)
  print '("i=",i2," xout=",f12.5," lebesgue=",f12.5)',(i,xOut(i),leb(1,i),i=1,nVert)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '("Cleanning Memory")'
  deallocate(xOut,lxOut,dlXout,leb)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
  
end subroutine edge_00


program main
  
  !> Test segments
  call edge_00()
  
end program main