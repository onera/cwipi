

subroutine pyramTestMass()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use baseSimplexTools
  use basePyramid
  use pyramidRule, only: P5_gauss
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer              :: i,j,iGauss,np,nPt
  integer              :: dim
  integer              :: ord,nMod
  integer              :: nNod
  real(8), pointer     :: uvw(:,:),w(:)
!  real(8), pointer     :: x(:),y(:),z(:),w(:)
  
  real(8), pointer     :: a(:),b(:),c(:)
  real(8), pointer     :: li(:,:),psi(:,:)
  real(8), pointer     :: vand(:,:)
  real(8), pointer     :: mass0(:,:),mass1(:,:),mass2(:,:),imas(:,:)
  real(8), pointer     :: mas1(:),mas2(:)
  real(8), pointer     :: eigv(:)
  real(8), parameter   :: eps=1d-12
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Calcul Matrice de Masse")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ord=10
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !call pyramidNodes   (ord=ord, uvw=uvw, display=.false.)  !> Points réguliers
  call pyramidNodesOpt(ord=ord, uvw=uvw, display=.false.)  !> Points optimises
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
  call pyramidVandermonde3D(ord=ord,a=a,b=b,c=c,vand=vand)
 !call display    (title="    Vandermonde Matrix",mat=vand)
 !call mathematica(title="    Vandermonde Matrix",mat=vand)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(a,b,c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/,">>> Test0 :  Base orthonormale")')
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Liste des points de Gauss
  call P5_gauss(    &
  &    order=ord   ,&
  &    nGauss=nPt  ,&
  &    uvwGauss=uvw,&
  &    pGauss=w     )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des fonctions de base orthonormales aux points uvw
  !> Transpose=.true. => psi(np,nPt)
  call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
  call pyramidBasePi(ord=ord,a=a,b=b,c=c,mode=psi,transpose=.true.) 
  deallocate(a,b,c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse = \int_0^1 \int_{-1+w}^{1-w} \int_{-1+w}^{1-w} psi_i(u,v,w) psi_j(u,v,w) du dv dw
  np=size(psi,1)
  allocate( mass0(np,np)) ; mass0(:,:)=0d0
  do iGauss=1,nPt
    do j=1,np
      do i=1,np
        mass0(i,j)=mass0(i,j)+psi(i,iGauss)*psi(j,iGauss)*w(iGauss)
      enddo
    enddo
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(psi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  call displaySparce(title="    Mass Matrix 0",mat=mass0,tol=1d-14)
 !print '(/4x,"det(mass0)=",e22.15/)',DMGT(eps=eps,n=size(mass0,1),A=mass0)
 !call eigenVectors(mat=mass0,w=eigv,display=.true.)
 !deallocate(eigv)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  write(*,'(/"<<< Test0 :  Base orthogonale")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/,">>> Test1 :  Mass = (vand x vand^t )^{-1}")')
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse
  call massMatrix(vand=vand,mass=mass1)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call display(title="    Mass Matrix 1",mat=mass1)
  print '(/4x,"det(mass1)=",e22.15/)',DMGT(eps=eps,n=size(mass1,1),A=mass1)
  call eigenVectors(mat=mass1,w=eigv,display=.true.)
  deallocate(eigv)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  write(*,'("<<< Test1 : Mass= (vand x vand^t )^{-1}")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/">>> Test2 : Mass=\int ai aj dV avec quadratures de Gauss")')
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !> Liste des points de Gauss
  call P5_gauss(    &
  &    order=ord   ,&
  &    nGauss=nPt  ,&
  &    uvwGauss=uvw,&
  &    pGauss=w     )
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !> Evaluation des fonctions de Lagrange aux points de Gauss uvw
  !> Transpose=.true. => li(np,nPt)
  nMod=(ord+1)*(ord+2)*(2*ord+3)/6 ; nNod=size(a)
  allocate(li(1:nMod,1:nNod)) ; li(:,:)=0d0
  
  call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
  call pyramidLagrange3Dv(ord=ord,vand=vand,a=a,b=b,c=c,lx=li,transpose=.true.)  !> li(uvw) = Inverse[Transpose[Vand]].Psi[uvw]
  deallocate(a,b,c)
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse
  np=size(li,1)
  allocate( mass2(1:np,1:np)) ; mass2(1:np,1:np)=0d0
  do iGauss=1,nPt
    do j=1,np
      do i=1,np
        mass2(i,j)=mass2(i,j)+li(i,iGauss)*li(j,iGauss)*w(iGauss)
      enddo
    enddo
  enddo
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(li)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  call display(title="    Mass Matrix 2",mat=mass2)
  print '(/4x,"det(mass2)=",e22.15/)',DMGT(eps=eps,n=size(mass1,1),A=mass2)
  call eigenVectors(mat=mass2,w=eigv,display=.true.)
  deallocate(eigv)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  call displaySparce(title="    Mass2-Mass1 (éléments <1d-15)",mat=mass2-mass1,tol=1d-15)
 !call display(title="    Mass2-mass1 ",mat=mass2-mass1)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  write(*,'("<<< Test2 : Mass=\int ai aj dV"/)')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/,">>> Test3 :  Mass= (vand x vand^t )^{-1} Format condensé")')
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse
  nMod=(ord+1)*(ord+2)*(2*ord+3)/6
  allocate(mas1(1:nMod*(nMod+1)/2)) ; mas1(:)=0d0
  
  call massMatrix(vand=vand,mass=mas1)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( ord<3 )call display(title="    Mass Matrix 1 Format condensé",vec=mas1)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  write(*,'("<<< Test3 : Mass= (vand x vand^t )^{-1} Format condensé")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/">>> Test4 : Mass=\int ai aj dV   format condensé")')
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  dim=(ord+1)*(ord+2)*(2*ord+3)/6 ; dim=dim*(dim+1)/2
  allocate(mas2(1:dim))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramideMassMatrix(ord=ord,vand=vand,mass=mas2)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  if( ord<3 )call display(title="    Mass Matrix 2 Format condensé",vec=mas2)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !call displaySparce(title="    Mass2 - Mass1 Format condensé  (éléments <1d-15)",vec=mas2-mas1,tol=1d-15)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  write(*,'(/"<<< Test 4 : Mass=\int ai aj dV Format condensé"/)')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramTestMass


program main
  
  
  !> Test Matrice de Masse des pyramides
  call pyramTestMass()
  
end program main