
subroutine quadTest()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use space_LDLt      , only: factorise
  use baseSimplex1D
  use baseSimplexTools, only: massMatrix,display,compactForm
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer            :: ord
  integer            :: i,j,k,l,n,row0,col0
  integer            :: nMod
  real(8), pointer   :: vand(:,:),dVand(:,:)
  real(8), pointer   :: r(:)=>null()
  real(8), pointer   :: mass(:)
  real(8), pointer   :: massL2(:,:)
  real(8), pointer   :: massQ4(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ord=10
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> points d'interpolation dans [-1,1]
 !call nodes1D   (ord=ord,uvw=r,display=.true.)
  call nodes1Dopt(ord=ord,uvw=r,display=.true.)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de Vandermonde 1D
  call vandermonde1D(ord=ord,a=r,vand=vand)
  call display(title="Vandermonde = [P_j(xi_i)]",mat=vand)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse compacte
  nMod=ord+1
  allocate(mass(1:nMod*(nMod+1)/2)) ; mass(:)=0d0 
  call massMatrix(vand=vand,mass=mass)
  
  call display(title="Mass=",vec=mass)
  deallocate(mass)
  !> Factorisation de la matrice
  !call factorise(n=ord+1,mat=mass)
  !call display(title="Factorized(Mass)=",vec=mass)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  nMod=ord+1
  allocate(massL2(1:(nMod+1),1:(nMod+1))) ; massL2(:,:)=0d0 
  call massMatrix(vand=vand,mass=massL2)
  call display(title="MassL2=",mat=massL2)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  allocate(mass(1:nMod*(nMod+1)/2)) ; mass(:)=0d0 
  call compactForm(mat0=massL2,mat1=mass)
  call display(title="MassL2=",vec=mass)
  deallocate(mass)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Produit tensoriel
  nMod=ord+1
  allocate(massQ4(1:nMod*nMod,1:nMod*nMod)) ; massQ4(:,:)=0d0
  do i=1,nMod
    row0=(i-1)*nMod
    do j=1,nMod
      col0=(j-1)*nMod
      do k=1,nMod
        do l=1,nMod
         !print '("row0+k=",i3,2x,"col0+l=",i3)',row0+k,col0+l
          massQ4(row0+k,col0+l)=massL2(i,j)*massL2(k,l) !> Tenseur massQ4(i,j,k,l)
        enddo
      enddo
    enddo
  enddo
  
  deallocate(massL2)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call display(title="MassQ4=",mat=massQ4)
  nMod=ord+1
  allocate(mass(1:nMod*nMod*(nMod*nMod+1)/2)) ; mass(:)=0d0
  call compactForm(mat0=massQ4,mat1=mass)
  call display(title="CompactedForm(MassQ4)=",vec=mass)
  
  call factorise(n=nMod*nMod,mat=mass)
 !call factorise(n=nMod*nMod,mat=mass)
  call display(title="Factorized(MassQ4)=",vec=mass)
  deallocate(massQ4)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine quadTest


program main
    
  !> Test Quad
  call quadTest()
  
end program main