subroutine jacobiTest()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use baseSimplexTools
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !
  real(8), pointer :: jac(:),dJac(:)
  real(8), pointer :: jf (:),dJf (:)
  real(8), pointer :: u(:)
  real(8)          :: alpha,beta
  integer          :: iVert,nVert,iOrd,ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  print '(/"Polynomes de Jacobi"/)'
  
  nVert=10 ; allocate(u(1:nVert+1))
  u=[( (-1d0+real(2*iVert,kind=8)/real(nVert,kind=8)), iVert=0,nVert )]
  
  ! J(u_i,iOrd)
  alpha=3d0 ; beta=1d0 ; ord=2
  
  do iOrd=0,ord
    
    call jacobiP (alpha=alpha,beta=beta,n=iOrd,u=u,jf= jac)
    call dJacobiP(alpha=alpha,beta=beta,n=iOrd,u=u,jf=djac)
    
    call jacobi (u=u,alpha=alpha,beta=beta,n=iOrd,jf= jf)
    call dJacobi(u=u,alpha=alpha,beta=beta,n=iOrd,jf=dJf)
    
    print '(/"J^{alpha=",f5.2,", beta=",f5.2,"}_{",i2,"}(u)"/)',alpha,beta,iOrd
   !print '(3x,"u=",e22.15,3x,"J=",e22.15,3x,"dJ=",e22.15)',(u(iVert),jac(iVert,iOrd),djac(iVert,iOrd), iVert=0,nVert)
    print '(3x,"u=",f6.2,3x,"J=",e22.15,1x,e22.15)',(u(iVert),jac(iVert)-jf(iVert),dJac(iVert)-dJf(iVert), iVert=1,nVert+1)
    
    deallocate(jac,dJac)
    deallocate(jf ,dJf )
  enddo
  
  deallocate(u)
  
  return
end subroutine jacobiTest

subroutine testQuadratureGL()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use baseSimplexTools, only: gaussLobattoQuadratures,gaussLegendreQuadratures
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  real(8), pointer   :: xGL (:),wGL (:)
  real(8), pointer   :: xGLL(:),wGLL(:)
  real(8), pointer   :: xGJ (:),wGJ (:)
  integer            :: i,order
  real(8)            :: int
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '(/"Gauss Legendre Quadatures")'
  order=4
  call gaussLegendreQuadratures(ord=order,xGL=xGL,wGL=wGL)
  print '("i=",i2," xGL=",f12.5," wGL=",f12.5)',(i,xGL(i),wGL(i),i=1,size(xGL))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  int=0d0
  do i=1,order+1
    int=int+wGL(i)*func( xGL(i) )
  enddo
  
  print '("\int_{-1}^{+1} 2 - 3x + x^2 + 2x^3 - 6x^4 + 8x^5- 19x^6 + x^7 dx = ",e22.15)',int
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '(/"Gauss Legendre-Lobatto Quadatures")'
  order=4
  call gaussLobattoQuadratures(ord=order,xGL=xGLL,wGL=wGLL)
  print '("i=",i2," xGL=",f12.5," wGL=",f12.5)',(i,xGLL(i),wGLL(i),i=1,size(xGLL))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  int=0d0
  do i=1,order+1
    int=int+wGLL(i)*func( xGLL(i) )
  enddo
  
  print '("\int_{-1}^{+1} 2 - 3x + x^2 + 2x^3 - 6x^4 + 8x^5- 19x^6 + x^7 dx = ",e22.15)',int
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(xGL ,wGL )
  deallocate(xGLL,wGLL)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
  
contains
  
  function func(x) result(fx)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8) :: x,fx
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    fx=2d0-3d0*x+x**2+2d0*x**3-6d0*x**4+8d0*x**5-19d0*x**6+x**7
    
    return
  end function func
  
end subroutine testQuadratureGL


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
  write(*,'("Calcul des bases polynômiales 1D")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)order
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
  real(8)            :: node_xy(2,3) !> Triangle
  real(8), pointer   :: node_uv(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'("Calcul des bases polynômiales 2D")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Points d''interpolation")')
  call nodes2D   (ord=order,uvw=uvw,display=.true.)
  call nodes2Dopt(ord=order,uvw=uvw,display=.true.)
  call nodes2Duv2rs(uv=uvw,rs=rs ,display=.true.) !> rs(1:2,:)=2d0*uv(1:2,:)-1d0
  call nodes2Drs2ab(rs=rs,a=a,b=b,display=.true.) !> a=2 (1+r)/(1-s)-1 && b=s
  
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
  !> Polynomes d'interpolation (on teste avec les points d'interpolation
  
  allocate(lxOut(1:(order+1)*(order+2)/2,1:nVert))
  call lagrange2Dv(ord=order,vand=vand,a=a,b=b,lx=lxOut,transpose=.true.) !> true pour affichage
  if( order<10 )then
    call display(title="Test avec l(uvw)",mat=lxOut)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse
  call massMatrix(vand=vand,mass=mass)
  if( order<10 )then
    call display(title="Mass Matrix",mat=mass)
   !call mathematica(title="Mass Matrix",mat=mass)
  endif
  print '(/"det(mass)=",e22.15/)',DMGT(eps=eps,n=size(mass,1),A=mass)
  call eigenVectors(mat=mass,w=eigv,display=.true.)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Derivees base fonctionnelle
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
  write(*,'(/"Order: ")',advance='no') ; read(*,*)order
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
  write(*,'(/"Order: ")',advance='no') ; read(*,*)order
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


subroutine tetraTest()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use baseSimplex2D
  use baseSimplex3D
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer              :: i,j,order,nMod,nNod,ad,np,nPt,cpt
  real(8), pointer     :: vand(:,:),dVand(:,:),jf(:,:),dr(:,:)
  real(8), pointer     :: drVand  (:,:),dsVand  (:,:),dtVand  (:,:)
  real(8), pointer     :: drMatrix(:,:),dsMatrix(:,:),dtMatrix(:,:)
  real(8), pointer     :: mass(:,:)
  real(8), pointer     :: psi(:,:)
  real(8), pointer     :: xyzOut(:,:),lxOut(:,:),drLxOut(:,:),dsLxOut(:,:),dtLxOut(:,:),leb(:,:)
  real(8), pointer     :: uvw(:,:),rst(:,:),a(:),b(:),c(:)
  real(8), parameter   :: eps=1d-12
  real(8), pointer     :: eigv(:)
  real(8), pointer     :: fi(:)
  real(8)              :: f0
  character(3)         :: sfx
  real(8), pointer     :: xGLL(:)
  real(8), pointer     :: uv(:,:)
  integer, allocatable :: conec(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'("Calcul des bases polynômiales Tetra")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(   1<=order .and. order<  10 ) write(sfx,'("00",i1)')order
  if(  10<=order .and. order< 100 ) write(sfx,'("0" ,i2)')order
  if( 100<=order .and. order<1000 ) write(sfx,'(     i3)')order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Tétraèdre points d''interpolation")')
  call nodes3D   (ord=order,uvw=uvw,display=.false.)
  call nodes3Dopt(ord=order,uvw=uvw,display=.true. )
  
  call nodes3Duvw2rst(uvw=uvw,rst=rst) !> rst(1:3,:)=2d0*uvw(1:3,:)-1d0
 !write(*,'(/"Tetra (rst):")')
 !print '("rst(1:3,",i2,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,rst(1:3,ad),ad=1,size(rst,2))
  
  call nodes3Drst2abc(rst=rst,a=a,b=b,c=c)
 !write(*,'(/"Tetra (abc):")')
 !print '("a,b,c(",i2,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,a(ad),b(ad),c(ad),ad=1,size(a))
  
  nNod=size(uvw,2)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Extraction des faces et vérification des coordonnées barycentriques
  if( 0==1 )then
    call gaussLegendreLobatto(ord=order,xGLL=xGLL)
    print '(/"Point de Gauss Lobatto")'
    print '("ad=",i5,2x,"u=",f19.16)',(ad,5d-1*(xGLL(ad)+1d0),ad=1,size(xGLL))
    
    call nodes2D   (ord=order,uvw=uv,display=.false.)
    call nodes2Dopt(ord=order,uvw=uv,display=.true. )
    
    call trianglesConnectivity(ord=order,conec=conec)
    print '()'
   !do i=1,4
    do i=3,3
      print '("Triangle",i1)',i
      do j=1,size(conec,1)
        ad=conec(j,i)
        !print '("dg=",i6," uvw=",3(f19.16,1x))',conec(j,i),uvw(1:3,conec(j,i))
        !print '("dg=",i5,2x,"u=",f19.16,2x,"v=",f19.16,2x,"w=",f19.16)',ad,uvw(1:3,ad)
        print '("dg=",i5,2x,"du=",f19.16,2x,"dw=",f19.16)',ad,uvw(1,ad)-uv(1,j),uvw(3,ad)-uv(2,j)
      enddo
    enddo
    
    do j=1,size(conec,1)
      ad=conec(j,3)
      uv(1,j)=uvw(1,ad)
      uv(2,j)=uvw(3,ad)
      uv(3,j)=1d0-uv(1,j)-uv(2,j)
    enddo
    deallocate(conec,uv)
    
    call nodes3Dopt_2D(ord=order,uvw=uv,display=.true.)
    stop
    
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call writeMesh3D(ord=order,uvw=uvw)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call vandermonde3D(ord=order,a=a,b=b,c=c,vand=vand)
  if( order<3 )then
    call display(title="Vandermonde Matrix",mat=vand)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Polynomes d'interpolation (on teste avec les points d'interpolation) => Matrice Identité
  nMod=(order+1)*(order+2)*(order+3)/6
  allocate(lxOut(1:nMod,1:nNod))                  !> allocation pour true
  call lagrange3Dv(ord=order,vand=vand,a=a,b=b,c=c,lx=lxOut,transpose=.true.) !> true pour affichage
  if( order<3 )then
    call display(title="Test avec l(uvw)",mat=lxOut)
  else
    call displaySparce(title="Test avec l(uvw)",mat=lxOut,tol=1d-14)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse
  call massMatrix(vand=vand,mass=mass)
  if( order<3 )then
    call display(title="Mass Matrix",mat=mass)
  endif
  print '(/"det(mass)=",e22.15/)',DMGT(eps=eps,n=size(mass,1),A=mass)
  call eigenVectors(mat=mass,w=eigv,display=.true.)
  deallocate(eigv)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Derivees base fonctionnelle
  call gradVandermonde3D(ord=order,a=a,b=b,c=c,drVand=drVand,dsVand=dsVand,dtVand=dtVand)
  call derive1D(vand=vand,dVand=drVand,dMat=drMatrix) !> drMatrix = drVand.Inverse[vand]
  call derive1D(vand=vand,dVand=dsVand,dMat=dsMatrix) !> dsMatrix = dsVand.Inverse[vand]
  call derive1D(vand=vand,dVand=dtVand,dMat=dtMatrix) !> dtMatrix = dtVand.Inverse[vand]
  if( order<3 )then
    call display(title="drVand Matrix",mat=drVand)
    call display(title="dsVand Matrix",mat=dsVand)
    call display(title="dtVand Matrix",mat=dtVand)
    call display(title="drMatrix=drVand.vand^{-1}",mat=drMatrix)
    call display(title="dsMatrix=dsVand.vand^{-1}",mat=dsMatrix)
    call display(title="dsMatrix=dtVand.vand^{-1}",mat=dtMatrix)
  endif
  deallocate(drVand,dsVand,dtVand)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(a,b,c,rst)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des fonctions de Lagrange aux points xyzOut
  call readXYZout3D(xyzOut=xyzOut)
  
  nMod=(order+1)*(order+2)*(order+3)/6 ; nNod=size(xyzOut,2)
  
  call nodes3Drst2abc(rst=xyzOut,a=a,b=b,c=c)
  call simplex3D  (ord=order,a=a,b=b,c=c,mode=psi,transpose=.false.)           !> Psi(xyzOut)
  
  
  allocate(lxOut(1:nNod,1:nMod))
  call lagrange3Dv(ord=order,vand=vand,a=a,b=b,c=c,lx=lxOut,transpose=.false.)  !> lxOut= Inverse[Transpose[Vand]].Psi[xyzOut] lxOut(nPt,np)
  if( order<3 )then
    call writeSolOut3D(title="simplex3D"//sfx,solOut=psi)
    call writeSolOut3D(title="lagrange3D"//sfx  ,solOut=lxOut  )
  endif
  deallocate(psi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des dérivées des fonctions de Lagrange aux points xyzOut
  allocate(drLxOut(1:nNod,1:nMod))
  call dLagrange1Dv(dMat=drMatrix,lx=lxOut,dlx=drLxOut,transpose=.false.) !> drLxOut(1:nPt,1:np)= Transpose[drMatrix] lxOut
  allocate(dsLxOut(1:nNod,1:nMod))
  call dLagrange1Dv(dMat=dsMatrix,lx=lxOut,dlx=dsLxOut,transpose=.false.) !> dsLxOut= Transpose[dsMatrix] lxOut
  allocate(dtLxOut(1:nNod,1:nMod))
  call dLagrange1Dv(dMat=dtMatrix,lx=lxOut,dlx=dtLxOut,transpose=.false.) !> dtLxOut= Transpose[dtMatrix] lxOut
  
!  drLxOut=2d0*drLxOut ! car u'=2u-1
!  dsLxOut=2d0*dsLxOut ! car v'=2v-1
!  dtLxOut=2d0*dtLxOut ! car w'=2w-1
  
  if( order<3 )then
    call writeSolOut3D(title="drLagrange3D"//sfx,solOut=drLxOut)
    call writeSolOut3D(title="dsLagrange3D"//sfx,solOut=dsLxOut)
    call writeSolOut3D(title="dtLagrange3D"//sfx,solOut=dtLxOut)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions : preliminaires
  print '("size(lxOut)=",i6," x ",i3)',size(lxOut,1),size(lxOut,2)
  allocate(fi(1:nMod)) ; fi(1:nMod)=1d0
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions lxOut : f(xyzOut) = \sum_{i=1,nMod} lxOut(xyzOut,i) f_i
  cPt=0
  do i=1,nNod
    f0=0d0
    do j=1,nMod
      f0=f0+lxOut(i,j)*fi(j)
    enddo
    if( abs(f0-1d0)>eps )then
      print '("f(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(f)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions drLxOut : df/dr(xyzOut) = \sum_{i=1,nMod} drLxOut(xyzOut,i) f_i
  cPt=0
  do i=1,nNod
    f0=0d0
    do j=1,nMod
      f0=f0+drLxOut(i,j)*fi(j)
    enddo
    if( abs(f0)>eps )then
      print '("df/dr(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(df/dr)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions dsLxOut : df/ds(xyzOut) = \sum_{i=1,nMod} dsLxOut(xyzOut,i) f_i
  cPt=0
  do i=1,nNod
    f0=0d0
    do j=1,nMod
      f0=f0+dsLxOut(i,j)*fi(j)
    enddo
    if( abs(f0)>eps )then
      print '("df/ds(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(df/ds)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions dtLxOut : df/dt(xyzOut) = \sum_{i=1,nMod} dtLxOut(xyzOut,i) f_i
  cPt=0
  do i=1,nNod
    f0=0d0
    do j=1,nMod
      f0=f0+dtLxOut(i,j)*fi(j)
    enddo
    if( abs(f0)>eps )then
      print '("df/dt(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(df/dt)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(fi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  allocate(leb(1:nNod,1))
  call lebesgue(lx=lxout,l=leb,transpose=.false.)
  call writeSolOut3D(title="lebesgue3DP"//sfx,solOut=leb)
  
  deallocate(leb)
  deallocate(xyzOut)
  deallocate(a,b,c)
  !
  deallocate(lxOut)
  deallocate(drLxOut,dsLxOut,dtLxOut)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(vand)
  deallocate(drMatrix,dsMatrix,dtMatrix)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine tetraTest

subroutine tetraMaillageVisu()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Cette procedure sert à construire les maillages de visu pour le tetra d'ordre élevé
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use baseSimplex3D
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer            :: order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Construction maillage Tetra P_i")')
  write(*,'("Warning ghs3d is required")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Building Mesh for Tetra P_{ord}
  call writeMeshSkin3D(ord=order)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine tetraMaillageVisu

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
  real(8), pointer   :: uvw(:,:)
  real(8), pointer   :: uvw0(:,:)
  integer, pointer   :: tetra(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Construction maillage Tetra P_i")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !do iOrd=1,ord
  do iOrd=ord,ord
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call nodes3D(ord=iOrd,uvw=uvw,display=.true.)
    nVert=size(uvw,2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(uvw0(1:3,1:nVert))
    uvw0(1:3,:)=uvw(1:3,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call driverTetMesh(node_xyz=uvw0,tetra_node=tetra)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call nodes3Dopt(ord=iOrd,uvw=uvw,display=.true. )
    uvw0(1:3,:)=uvw(1:3,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call saveTetMesh(ord=iOrd, node_xyz=uvw0,tetra_node=tetra)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw,tetra,uvw0)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  return
end subroutine tetraMaillageVisuNew

subroutine testConnectivitiesOfSides()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use baseSimplex2D, only: edgesConnectivity
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer              :: i,j,order,nVert,ad
  integer, allocatable :: conec(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( order<10 )then ; write(*,'(/"Triangle P",i1)')order
  else               ; write(*,'(/"Triangle P",i2)')order
  endif
  
  call edgesConnectivity(ord=order,conec=conec)
  
  write(*,'(/"Local DOF / side")')
  do j=1,3
    write(*,'("sd",i1,": ",$)')j
    do i=1,order+1
      write(*,'(i4," ",$)')conec(i,j)
    enddo
    write(*,'()')
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(conec)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine testConnectivitiesOfSides


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
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
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


subroutine pyramBasis()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use basePyramid
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer              :: i,j,nPt,cpt
  integer              :: ord,nMod
  real(8), pointer     :: uv (:,:)
  real(8), pointer     :: uvw(:,:)
  integer              :: ad,iSide
  integer              :: iNod,nNod,iu,iv,iw
  integer, allocatable :: conec(:,:)
  integer, allocatable :: idx(:)
  real(8)              :: rot(3,3),xyz(1:3),cos_a,sin_a
  real(8)              :: alpha
  
  real(8), parameter   :: eps=1d-12
  real(8), pointer     :: a(:),b(:),c(:)
  real(8), pointer     :: vand(:,:),dVand(:,:) !,jf(:,:),dr(:,:)
  real(8), pointer     :: duPsi   (:,:),dvPsi  (:,:),dwPsi  (:,:)
  real(8), pointer     :: drMatrix(:,:),dsMatrix(:,:),dtMatrix(:,:)
  real(8), pointer     :: mass(:,:)
  real(8), pointer     :: xyzOut(:,:),lxOut(:,:),drLxOut(:,:),dsLxOut(:,:),dtLxOut(:,:),leb(:,:)
  real(8), pointer     :: mode(:,:)
  real(8), pointer     :: eigv(:)
  real(8), pointer     :: fi(:)
  real(8)              :: f0
  
  character(3)         :: sfx
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'("Calcul des bases polynômiales Pyram")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !call pyramidNodes   (ord=ord, uvw=uvw, display=.true.)  !> Points réguliers
  call pyramidNodesOpt(ord=ord, uvw=uvw, display=.true.)  !> Points optimises
 !call writeMesh3D    (ord=ord, uvw=uvw)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramidSkin3D(ord=ord, uvw=uvw, display=.true.)
  call pyramidMesh3D(ord=ord, uvw=uvw, display=.true.)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
 !write(*,'(/"Pyramid (abc):")')
 !print '("a,b,c(",i2,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,a(ad),b(ad),c(ad),ad=1,size(a))
  call pyramidVandermonde3D(ord=ord,a=a,b=b,c=c,vand=vand)
  if( ord<3 )then
    call display(title="Vandermonde Matrix",mat=vand)
  else
    call displaySparce(title="Test avec l(uvw)",mat=lxOut,tol=1d-14)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Polynomes d'interpolation (on teste avec les points d'interpolation) => Matrice Identité
  
  nMod=(ord+1)*(ord+2)*(2*ord+3)/6
  nNod=size(a)
  
  allocate(lxOut(1:nMod,1:nNod))
  call pyramidLagrange3Dv(ord=ord,vand=vand,a=a,b=b,c=c,lx=lxOut,transpose=.true.) !> true pour affichage
  if( ord<3 )then
    call display(title="Test avec l(uvw)",mat=lxOut)
  else
    call displaySparce(title="Test avec l(uvw)",mat=lxOut,tol=1d-12)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse
  call massMatrix(vand=vand,mass=mass)
  if( ord<3 )then
    call display(title="Mass Matrix",mat=mass)
  endif
 !print '(/"det(mass)=",e22.15/)',DMGT(eps=eps,n=size(mass,1),A=mass)
 !call eigenVectors(mat=mass,w=eigv,display=.true.)
 !deallocate(eigv)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Derivees base fonctionnelle
  call pyramidGradVandermonde3D(ord=ord,a=a,b=b,c=c,drMode=duPsi,dsMode=dvPsi,dtMode=dwPsi)
  call derive1D(vand=vand,dVand=duPsi,dMat=drMatrix) !> drMatrix = duPsi.Inverse[vand]
  call derive1D(vand=vand,dVand=dvPsi,dMat=dsMatrix) !> dsMatrix = dvPsi.Inverse[vand]
  call derive1D(vand=vand,dVand=dwPsi,dMat=dtMatrix) !> dtMatrix = dwPsi.Inverse[vand]
  if( ord<3 )then
    !call display(title="duPsi Matrix",mat=duPsi)
    !call display(title="dvPsi Matrix",mat=dvPsi)
    !call display(title="dwPsi Matrix",mat=dwPsi)
    call display(title="drMatrix=duPsi.vand^{-1}",mat=drMatrix)
    call display(title="dsMatrix=dvPsi.vand^{-1}",mat=dsMatrix)
    call display(title="dsMatrix=dwPsi.vand^{-1}",mat=dtMatrix)
  endif
  deallocate(duPsi,dvPsi,dwPsi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(a,b,c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramBasis


subroutine pyramLebesgue()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use basePyramid
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer              :: i,j,order,np,nPt,cpt
  integer              :: ord
  real(8), pointer     :: uv (:,:)
  real(8), pointer     :: uvw(:,:)
  integer              :: ad,iSide
  integer              :: iNod,nNod,iu,iv,iw
  integer, allocatable :: conec(:,:)
  integer, allocatable :: idx(:)
  real(8)              :: rot(3,3),xyz(1:3),cos_a,sin_a
  real(8)              :: alpha
  !>
  real(8), parameter   :: eps=1d-12
  real(8), pointer     :: a(:),b(:),c(:)
  real(8), pointer     :: vand(:,:),dVand(:,:) !,jf(:,:),dr(:,:)
  real(8), pointer     :: duPsi  (:,:),dvPsi  (:,:),dwPsi  (:,:)
  real(8), pointer     :: drMatrix(:,:),dsMatrix(:,:),dtMatrix(:,:)
  real(8), pointer     :: mass(:,:)
  real(8), pointer     :: xyzOut(:,:),lxOut(:,:),drLxOut(:,:),dsLxOut(:,:),dtLxOut(:,:),leb(:,:)
  real(8), pointer     :: mode(:,:)
  real(8), pointer     :: eigv(:)
  real(8), pointer     :: fi(:)
  real(8)              :: f0
  
  character(3)         :: sfx
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'("Calcul fonction de Lebesgue pour pyramides")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !call pyramidNodes   (ord=ord, uvw=uvw, display=.false.)  !> Points réguliers
  call pyramidNodesOpt(ord=ord, uvw=uvw, display=.false.)  !> Points optimises
 !call writeMesh3D    (ord=ord, uvw=uvw)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
 !write(*,'(/"Pyramid (abc):")')
 !print '("a,b,c(",i2,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,a(ad),b(ad),c(ad),ad=1,size(a))
  call pyramidVandermonde3D(ord=ord,a=a,b=b,c=c,vand=vand)
  if( ord<3 )then
    call display(title="Vandermonde Matrix",mat=vand)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Derivees base fonctionnelle
  call pyramidGradVandermonde3D(ord=ord,a=a,b=b,c=c,drMode=duPsi,dsMode=dvPsi,dtMode=dwPsi)
  call derive1D(vand=vand,dVand=duPsi,dMat=drMatrix) !> drMatrix = duPsi.Inverse[vand]
  call derive1D(vand=vand,dVand=dvPsi,dMat=dsMatrix) !> dsMatrix = dvPsi.Inverse[vand]
  call derive1D(vand=vand,dVand=dwPsi,dMat=dtMatrix) !> dtMatrix = dwPsi.Inverse[vand]
  if( ord<3 )then
    call display(title="duPsi Matrix",mat=duPsi)
    call display(title="dvPsi Matrix",mat=dvPsi)
    call display(title="dwPsi Matrix",mat=dwPsi)
    call display(title="drMatrix=duPsi.vand^{-1}",mat=drMatrix)
    call display(title="dsMatrix=dvPsi.vand^{-1}",mat=dsMatrix)
    call display(title="dsMatrix=dwPsi.vand^{-1}",mat=dtMatrix)
  endif
  deallocate(duPsi,dvPsi,dwPsi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(a,b,c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
  if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
  if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Liste des points xyzOut  
  call pyramidReadXYZout3D(xyzOut=xyzOut, display=.true.)
  call pyramiduvw2abc(uvw=xyzOut,a=a,b=b,c=c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des fonctions de Lagrange aux points xyzOut
  call pyramidBasePi(ord=ord,a=a,b=b,c=c,mode=mode,transpose=.false.)                !> Psi(xyzOut)
  call pyramidLagrange3Dv(ord=ord,vand=vand,a=a,b=b,c=c,lx=lxOut,transpose=.false.)  !> lxOut= Inverse[Transpose[Vand]].Psi[xyzOut] lxOut(nPt,np)
  if( ord<3 )then
    call pyramidWriteSolOut3D(title="simplex3D"//sfx,solOut=mode   )
    call pyramidWriteSolOut3D(title="lagrange3D"//sfx,solOut=lxOut  )
  endif
  deallocate(mode)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call lebesgue(lx=lxout,l=leb,transpose=.false.) !> leb = Sum abs(lxOut)
  call pyramidWriteSolOut3D(title="lebesgue3DP"//sfx,solOut=leb)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des dérivées des fonctions de Lagrange aux points xyzOut
  call dLagrange1Dv(dMat=drMatrix,lx=lxOut,dlx=drLxOut,transpose=.false.) !> drLxOut(1:nPt,1:np)= Transpose[drMatrix] lxOut
  call dLagrange1Dv(dMat=dsMatrix,lx=lxOut,dlx=dsLxOut,transpose=.false.) !> dsLxOut= Transpose[dsMatrix] lxOut
  call dLagrange1Dv(dMat=dtMatrix,lx=lxOut,dlx=dtLxOut,transpose=.false.) !> dtLxOut= Transpose[dtMatrix] lxOut
  
  if( ord<3 )then
    call pyramidWriteSolOut3D(title="drLagrange3D"//sfx,solOut=drLxOut)
    call pyramidWriteSolOut3D(title="dsLagrange3D"//sfx,solOut=dsLxOut)
    call pyramidWriteSolOut3D(title="dtLagrange3D"//sfx,solOut=dtLxOut)
  endif  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions : preliminaires
  print '("size(lxOut)=",i6," x ",i3)',size(lxOut,1),size(lxOut,2)
  nPt=size(lxOut,1)
  np =(ord+1)*(ord+2)*(2*ord+3)/6 !> = \sum_{k=1}^{ord+1} k^2  
  allocate(fi(1:np)) ; fi(1:np)=1d0
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions lxOut : f(xyzOut) = \sum_{i=1,np} lxOut(xyzOut,i) f_i
  print '("Test des fonctions lxOut : f(xyzOut) = \sum_{i=1,np} lxOut(xyzOut,i) f_i")'
  cPt=0
  do i=1,nPt
    f0=0d0
    do j=1,np
      f0=f0+lxOut(i,j)*fi(j)
    enddo
    if( abs(f0-1d0)>eps )then
      print '("f(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(f)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions drLxOut : df/dr(xyzOut) = \sum_{i=1,np} drLxOut(xyzOut,i) f_i
  print '("Test des fonctions drLxOut : df/dr(xyzOut) = \sum_{i=1,np} drLxOut(xyzOut,i) f_i")'
  cPt=0
  do i=1,nPt
    f0=0d0
    do j=1,np
      f0=f0+drLxOut(i,j)*fi(j)
    enddo
    if( abs(f0)>eps )then
      print '("df/dr(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(df/dr)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions dsLxOut : df/ds(xyzOut) = \sum_{i=1,np} dsLxOut(xyzOut,i) f_i
  print '("Test des fonctions dsLxOut : df/ds(xyzOut) = \sum_{i=1,np} dsLxOut(xyzOut,i) f_i")'
  cPt=0
  do i=1,nPt
    f0=0d0
    do j=1,np
      f0=f0+dsLxOut(i,j)*fi(j)
    enddo
    if( abs(f0)>eps )then
      print '("df/ds(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(df/ds)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions dtLxOut : df/dt(xyzOut) = \sum_{i=1,np} dtLxOut(xyzOut,i) f_i
  print '("Test des fonctions dtLxOut : df/dt(xyzOut) = \sum_{i=1,np} dtLxOut(xyzOut,i) f_i")'
  cPt=0
  do i=1,nPt
    f0=0d0
    do j=1,np
      f0=f0+dtLxOut(i,j)*fi(j)
    enddo
    if( abs(f0)>eps )then
      print '("df/dt(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(df/dt)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(fi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(xyzOut)
  deallocate(a,b,c)
  !
  deallocate(lxOut)
  deallocate(drLxOut,dsLxOut,dtLxOut)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(vand)
  deallocate(drMatrix,dsMatrix,dtMatrix)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramLebesgue

subroutine pyramMaillageVisu()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !use modDeterminant
  !use baseSimplex2D
  !use baseSimplex3D
  use basePyramid
  use table_tet_mesh
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer            :: ord,iOrd,ad
  real(8), pointer   :: uvw(:,:)
  integer, pointer   :: tetra(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Construction maillage Pyramid P_i")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do iOrd=1,ord
  !do iOrd=ord,ord
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call pyramidNodes   (ord=iOrd, uvw=uvw, display=.false.)
    call driverTetMesh  (node_xyz=uvw,tetra_node=tetra)
    deallocate(uvw)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call pyramidNodesOpt(ord=iOrd, uvw=uvw, display=.true.)  !> Points optimises
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call saveTetMesh(ord=iOrd, node_xyz=uvw,tetra_node=tetra)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !call pyramidSkin3D(ord=iOrd, uvw=uvw, display=.true.)
    !call pyramidMesh3D(ord=iOrd, uvw=uvw, display=.true.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw,tetra)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    print '("nn=",i10," ne=",i10," Pi=",i3)',size(uvw,2),size(tetra,2),iOrd
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramMaillageVisu

subroutine pyramTestQuadrature()
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !> Routine permettant de tester les
  !> bases fonctionnelles sur une liste
  !> de points
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use pyramidRule, only: P5_gauss
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer              :: i
  integer              :: iOrd,ord
  !>
  integer              :: n
  real(8), pointer     :: x(:),y(:),z(:),w(:)
  integer              :: power
  real(8), pointer     :: f(:)
  real(8)              :: s
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Contrôle des quadratures")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do iOrd=0,ord
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Liste des points xyzOut
   !write(*,'(/"Points de Gauss:")')
    call P5_gauss(   &
    &    order=iOrd ,&
    &    nGauss=n   ,&
    &    uGauss=x   ,&
    &    vGauss=y   ,&
    &    wGauss=z   ,&
    &    pGauss=w    )
    
    print '(4x,"3D Quadrature")'
    do i=1,n
      print '(8x,"i",i3,": uvw=",3(f12.5,1x),2x,"p=",f12.5)',i,x(i),y(i),z(i),w(i)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(f(1:n))
   !print '("f(x,y,z)= (x+1)^",i1," (y+1)^",i1," z^",i1)',n,n,n
    power=(2*iOrd+1)/3
    do i=1,n
      f(i)=(x(i)+1)**power *(y(i)+1)**power *(z(i))**power
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Test des quadratures
    s=sum([(f(i)*w(i), i=1,n)])
    if(                 power<10  )print '("\int_0^1 \int_{-(1-z)}^{+(1-z)} \int_{-(1-z)}^{+(1-z)} (x+1)^0",i1," (y+1)^0",i1," z^0",i1," dx dy dz = ",e22.15)',power,power,power,s
    if( 10<=power .and. power<100 )print '("\int_0^1 \int_{-(1-z)}^{+(1-z)} \int_{-(1-z)}^{+(1-z)} (x+1)^" ,i2," (y+1)^" ,i2," z^" ,i2," dx dy dz = ",e22.15)',power,power,power,s
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(x,y,z,w,f)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramTestQuadrature

subroutine pyramTestBasis()
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !> Routine permettant de tester les
  !> bases fonctionnelles sur une liste
  !> de points
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use basePyramid
  use pyramidRule, only: P5_gauss
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer              :: i,j,np,nPt
  integer              :: iOrd,ord
  real(8), pointer     :: uvw(:,:)=>null()
  !>
  real(8), pointer     :: x(:),y(:),z(:),w(:)
  real(8), pointer     :: a(:),b(:),c(:)
  real(8), pointer     :: vand(:,:),dVand(:,:)
  real(8), pointer     :: duPsi   (:,:),dvPsi   (:,:),dwPsi   (:,:)
  real(8), pointer     :: drMatrix(:,:),dsMatrix(:,:),dtMatrix(:,:)
  real(8), pointer     :: li(:,:),duLi(:,:),dvLi(:,:),dwLi(:,:)
  real(8), pointer     :: fi(:),dxfi(:),dyfi(:),dzfi(:)
  real(8), pointer     :: f (:),dxf (:),dyf (:),dzf (:)
  real(8)              :: f0,dxf0,dyf0,dzf0,delta,deltaMax
  real(8) , parameter  :: tol=1d-11
  integer              :: cpt0,cpt1,cpt2,cpt3,cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Contrôle de la base fonctionnelle pyramidale")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> CONSTRUCTION DES BASES
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Construction des bases")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !call pyramidNodes   (ord=ord, uvw=uvw, display=.false.)  !> Points réguliers
  call pyramidNodesOpt(ord=ord, uvw=uvw, display=.false.)  !> Points optimises
  call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> f(u,v,w),∂xf(u,v,w),∂yf(u,v,w),∂zf(u,v,w)
  nPt=size(uvw,2) ; allocate(fi(nPt),dxfi(nPt),dyfi(nPt),dzfi(nPt))
  call fxyz(xyz=uvw,f=fi,dxf=dxfi,dyf=dyfi,dzf=dzfi)
 !print '(3x,"ad=",i6,2x,"x=",f12.5,2x,"y=",f12.5,2x,"z=",f12.5,4x,"fi=",e22.15)',(i,uvw(1:3,i),fi(i),i=1,nPt)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Preparation calcul de la base
  call pyramidVandermonde3D(ord=ord,a=a,b=b,c=c,vand=vand)
 !call display(title="    vand Matrix",mat=vand)
  print '(/4x,"vand")'
  do i=1,nPt
    print '(4x,"i=",i4," vand =",14(e9.2,1x))',i,vand(i,:)
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>  Preparation calcul des dérivées de la base
  
  call pyramidGradLagrange3Dv(    &
  &    ord=ord,vand=vand         ,&
  &    a=a,b=b,c=c               ,&
  &    drPhi=drMatrix            ,& !> duai(np,nPt) ; nPt=size(u)  duai:= Inverse[Transpose[Vand]].duPsi[x]
  &    dsPhi=dsMatrix            ,& !> duai(np,nPt) ; nPt=size(u)  dvai:= Inverse[Transpose[Vand]].dvPsi[x]
  &    dtPhi=dtMatrix            ,& !> duai(np,nPt) ; nPt=size(u)  dwai:= Inverse[Transpose[Vand]].dwPsi[x]
  &    transpose=.true.           )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(a,b,c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> PyramideP1
  if( ord==1 )then
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(*,'(/">>> Test PyramideP1")')
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !call pyramidNodes   (ord=10, uvw=uvw, display=.false.)  !> Points réguliers
    call pyramidNodesOpt(ord=10, uvw=uvw, display=.false.)  !> Points optimises
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> f(u,v,w),∂xf(u,v,w),∂yf(u,v,w),∂zf(u,v,w)
    nPt=size(uvw,2) ; allocate(f(1:nPt),dxf(nPt),dyf(nPt),dzf(nPt))
    call fxyz(xyz=uvw,f=f,dxf=dxf,dyf=dyf,dzf=dzf)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Evaluation des fonctions de Lagrange aux points uvw
    !> Transpose=.true. => li(np,nPt)
    call pyramidBaseP1  (uvw=uvw, ai=li, transpose=.true.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Evaluation des dérivées des fonctions de Lagrange aux points uvw
    call pyramidGradBaseP1(uvw=uvw,duai=duLi,dvai=dvLi,dwai=dwLi,transpose=.true.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> comparaison des resultats (calcul direct et calcul interpolé)
    np =(ord+1)*(ord+2)*(2*ord+3)/6 !> = \sum_{k=1}^{ord+1} k^2
    print'(4x,"np x nPt =",i6," x ",i6)',np,nPt
    
    !> test li
    cpt0=0 ; deltaMax=0d0
    do j=1,nPt
      f0=0d0
      do i=1,np
        f0=f0+li(i,j)*fi(i)
      enddo
      delta=abs(f(j)-f0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt0=cpt0+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"   f(uvw)- ∑   ai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),f(j),f0,f(j)-f0
      endif
    enddo
    print '(4x,"erreur sur   f cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt0,nPt,deltaMax
    
    !> test duLi
    cpt1=0 ; deltaMax=0d0
    do j=1,nPt
      dxf0=0d0
      do i=1,np
        dxf0=dxf0+duLi(i,j)*fi(i)
      enddo
      delta=abs(dxf(j)-dxf0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt1=cpt1+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"∂u f(uvw)- ∑ ∂uai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),dxf(j),dxf0,dxf(j)-dxf0
      endif
    enddo
    print '(4x,"erreur sur ∂uf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt1,nPt,deltaMax
    
    !> test dvLi
    cpt2=0 ; deltaMax=0d0
    do j=1,nPt
      dyf0=0d0
      do i=1,np
        dyf0=dyf0+dvLi(i,j)*fi(i)
      enddo
      delta=abs(dyf(j)-dyf0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt2=cpt2+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"∂vf(uvw)- ∑ ∂vai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),dyf(j),dyf0,dyf(j)-dyf0
      endif
    enddo
    print '(4x,"erreur sur ∂vf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt2,nPt,deltaMax
    
    !> test dwLi
    cpt3=0 ; deltaMax=0d0
    do j=1,nPt
      dzf0=0d0
      do i=1,np
        dzf0=dzf0+dwLi(i,j)*fi(i)
      enddo
      delta=abs(dzf(j)-dzf0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt3=cpt3+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"∂wf(uvw)- ∑ ∂wai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),dzf(j),dzf0,dzf(j)-dzf0
      endif
    enddo
    print '(4x,"erreur sur ∂wf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt3,nPt,deltaMax
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw)
    deallocate(li)
    deallocate(duLi,dvLi,dwLi)
    deallocate(f,dxf,dyf,dzf)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cpt=cpt1+cpt2+cpt3
    if( .not.cpt==0 )then
      print '("xxx Arret sur test PyramideP1")'
      stop
    else
      write(*,'("<<< Test PyramideP1")')
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> TEST BASIC DES BASES
  
  do iOrd=0,10
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(*,'(/">>> Test basic ",i2," ord=",i2)')iOrd,ord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !call pyramidNodes   (ord=iOrd, uvw=uvw, display=.false.)  !> Points réguliers
    call pyramidNodesOpt(ord=iOrd, uvw=uvw, display=.false.)  !> Points optimises
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> f(u,v,w),∂xf(u,v,w),∂yf(u,v,w),∂zf(u,v,w)
    nPt=size(uvw,2) ; allocate(f(1:nPt),dxf(nPt),dyf(nPt),dzf(nPt))
    call fxyz(xyz=uvw,f=f,dxf=dxf,dyf=dyf,dzf=dzf)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Evaluation des fonctions de Lagrange aux points uvw
    !> Transpose=.true. => li(np,nPt)
    call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
    call pyramidLagrange3Dv(ord=ord,vand=vand,a=a,b=b,c=c,lx=li,transpose=.true.)  !> li(uvw) = Inverse[Transpose[Vand]].Psi[uvw]
    deallocate(a,b,c)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Evaluation des dérivées des fonctions de Lagrange aux points uvw
    call dLagrange1Dv(dMat=drMatrix,lx=li,dlx=duLi,transpose=.true.) !> duLi(1:np,1:nPt)=Transpose[drMatrix] li
    call dLagrange1Dv(dMat=dsMatrix,lx=li,dlx=dvLi,transpose=.true.) !> dvLi(1:np,1:nPt)=Transpose[dsMatrix] li
    call dLagrange1Dv(dMat=dtMatrix,lx=li,dlx=dwLi,transpose=.true.) !> dwLi(1:np,1:nPt)=Transpose[dtMatrix] li
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> comparaison des resultats (calcul direct et calcul interpolé)
    np =(ord+1)*(ord+2)*(2*ord+3)/6 !> = \sum_{k=1}^{ord+1} k^2
    print'(4x,"np x nPt =",i6," x ",i6)',np,nPt
    
    !> test li
    cpt0=0 ; deltaMax=0d0
    do j=1,nPt
      f0=0d0
      do i=1,np
        f0=f0+li(i,j)*fi(i)
      enddo
      delta=abs(f(j)-f0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt0=cpt0+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"   f(uvw)- ∑   ai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),f(j),f0,f(j)-f0
      endif
    enddo
    print '(4x,"erreur sur   f cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt0,nPt,deltaMax
    
    !> test duLi
    cpt1=0 ; deltaMax=0d0
    do j=1,nPt
      dxf0=0d0
      do i=1,np
        dxf0=dxf0+duLi(i,j)*fi(i)
      enddo
      delta=abs(dxf(j)-dxf0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt1=cpt1+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"∂u f(uvw)- ∑ ∂uai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),dxf(j),dxf0,dxf(j)-dxf0
      endif
    enddo
    print '(4x,"erreur sur ∂uf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt1,nPt,deltaMax
    
    !> test dvLi
    cpt2=0 ; deltaMax=0d0
    do j=1,nPt
      dyf0=0d0
      do i=1,np
        dyf0=dyf0+dvLi(i,j)*fi(i)
      enddo
      delta=abs(dyf(j)-dyf0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt2=cpt2+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"∂vf(uvw)- ∑ ∂vai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),dyf(j),dyf0,dyf(j)-dyf0
      endif
    enddo
    print '(4x,"erreur sur ∂vf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt2,nPt,deltaMax
    
    !> test dwLi
    cpt3=0 ; deltaMax=0d0
    do j=1,nPt
      dzf0=0d0
      do i=1,np
        dzf0=dzf0+dwLi(i,j)*fi(i)
      enddo
      delta=abs(dzf(j)-dzf0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt3=cpt3+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"∂wf(uvw)- ∑ ∂wai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),dzf(j),dzf0,dzf(j)-dzf0
      endif
    enddo
    print '(4x,"erreur sur ∂wf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt3,nPt,deltaMax
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw)
    deallocate(li)
    deallocate(duLi,dvLi,dwLi)
    deallocate(f,dxf,dyf,dzf)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cpt=cpt1+cpt2+cpt3
    if( .not.cpt==0 )then
      print '("xxx Arret sur test basic ",i2)',iOrd
      stop
    else
      write(*,'("<<< Test basic ",i2," ord=",i2)')iOrd,ord
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> TEST AVANCE DES BASES
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/">>> Test avancé ord=",i2)')ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Liste des points uvw
  call pyramidReadXYZout3D(xyzOut=uvw, display=.true.)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> f(u,v,w),∂xf(u,v,w),∂yf(u,v,w),∂zf(u,v,w)
  nPt=size(uvw,2) ; allocate(f(nPt),dxf(nPt),dyf(nPt),dzf(nPt))
  call fxyz(xyz=uvw,f=f,dxf=dxf,dyf=dyf,dzf=dzf)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des fonctions de Lagrange aux points uvw
  !> Transpose=.true. => li(np,nPt)
  call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
  call pyramidLagrange3Dv(ord=ord,vand=vand,a=a,b=b,c=c,lx=li,transpose=.true.)  !> li(uvw) = Inverse[Transpose[Vand]].Psi[uvw]
  deallocate(a,b,c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des dérivées des fonctions de Lagrange aux points uvw
  call dLagrange1Dv(dMat=drMatrix,lx=li,dlx=duLi,transpose=.true.) !> duLi(1:np,1:nPt)=Transpose[drMatrix] li
  call dLagrange1Dv(dMat=dsMatrix,lx=li,dlx=dvLi,transpose=.true.) !> dvLi(1:np,1:nPt)=Transpose[dsMatrix] li
  call dLagrange1Dv(dMat=dtMatrix,lx=li,dlx=dwLi,transpose=.true.) !> dwLi(1:np,1:nPt)=Transpose[dtMatrix] li
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> comparaison des resultats (calcul direct et calcul interpolé)
  np =(ord+1)*(ord+2)*(2*ord+3)/6 !> = \sum_{k=1}^{ord+1} k^2
  print'(4x,"np x nPt =",i6," x ",i6)',np,nPt
  
  !> test li
  cpt0=0 ; deltaMax=0d0
  do j=1,nPt
    f0=0d0
    do i=1,np
      f0=f0+li(i,j)*fi(i)
    enddo
    delta=abs(f(j)-f0) ; if( delta>deltaMax)deltaMax=delta
    if( delta>tol )then
      cpt0=cpt0+1
      print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"f(uvw)- ∑ ai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),f(j),f0,f(j)-f0
    endif
  enddo
  print '(4x,"erreur sur   f cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt0,nPt,deltaMax
  
  !> test duLi
  cpt1=0 ; deltaMax=0d0
  do j=1,nPt
    dxf0=0d0
    do i=1,np
      dxf0=dxf0+duLi(i,j)*fi(i)
    enddo
    delta=abs(dxf(j)-dxf0) ; if( delta>deltaMax)deltaMax=delta
    if( delta>tol )then
      cpt1=cpt1+1
      print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"   f(uvw)- ∑   ai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),f(j),f0,f(j)-f0
    endif
  enddo
  print '(4x,"erreur sur ∂uf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt1,nPt,deltaMax
  
  !> test dvLi
  cpt2=0 ; deltaMax=0d0
  do j=1,nPt
    dyf0=0d0
    do i=1,np
      dyf0=dyf0+dvLi(i,j)*fi(i)
    enddo
    delta=abs(dyf(j)-dyf0) ; if( delta>deltaMax)deltaMax=delta
    if( delta>tol )then
      cpt2=cpt2+1
      print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"∂vf(uvw)- ∑ ∂vai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),dyf(j),dyf0,dyf(j)-dyf0
    endif
  enddo
  print '(4x,"erreur sur ∂vf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt2,nPt,deltaMax
  
  !> test dwLi
  cpt3=0 ; deltaMax=0d0
  do j=1,nPt
    dzf0=0d0
    do i=1,np
      dzf0=dzf0+dwLi(i,j)*fi(i)
    enddo
    delta=abs(dzf(j)-dzf0) ; if( delta>deltaMax)deltaMax=delta
    if( delta>tol )then
      cpt3=cpt3+1
      print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"∂wf(uvw)- ∑ ∂wai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),dzf(j),dzf0,dzf(j)-dzf0
    endif
  enddo
  print '(4x,"erreur sur ∂wf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt3,nPt,deltaMax
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  cpt=cpt0+cpt1+cpt2+cpt3
  if( .not.cpt==0 )then
    print '("xxx Erreur sur test avancé ord=",i2)',ord
    stop
  else
    write(*,'("<<< Test avancé ord=",i2)')ord
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(li,duLi,dvLi,dwLi)
  deallocate(f,dxf,dyf,dzf)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> TEST SUR POINTS DE GAUSS DES BASES
  
  do iOrd=ord,ord
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(*,'(/">>> Test sur points de Gauss ord=",i2)')ord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Liste des points xyzOut
   !write(*,'(/"Points de Gauss:")')
    call P5_gauss(   &
    &    order=iOrd ,&
    &    nGauss=nPt ,&
    &    uGauss=x   ,&
    &    vGauss=y   ,&
    &    wGauss=z   ,&
    &    pGauss=w    )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(uvw(1:3,nPt))
    do i=1,nPt
      uvw(1:3,i)=[x(i),y(i),z(i)]
    enddo
    
    print '(/4x,"uvw")'
    do i=1,nPt
      print '(4x,"i=",i4," uvw =",3(f9.6,1x))',i,uvw(1:3,i)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(f(nPt),dxf(nPt),dyf(nPt),dzf(nPt))
    call fxyz(xyz=uvw,f=f,dxf=dxf,dyf=dyf,dzf=dzf)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Evaluation des fonctions de Lagrange aux points uvw
    !> Transpose=.true. => li(np,nPt)
    call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
    call pyramidLagrange3Dv(ord=ord,vand=vand,a=a,b=b,c=c,lx=li,transpose=.true.)  !> li(uvw) = Inverse[Transpose[Vand]].Psi[uvw]
    deallocate(a,b,c)
    
    print '(/4x,"li(:,uvw)")'
    do i=1,nPt
      print '(4x,"i=",i4," li  =",14(f9.6,1x))',i,li(:,i)
    enddo
    
    open(unit=100,file="BasisP5_ai.dat",action='write')
    write(100,'(4x,"ai")')
    do i=1,nPt
      write(100,'(4x,"i=",i3,2x,14(e19.12,1x))')i,li(:,i)
    enddo
    close(100)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Evaluation des dérivées des fonctions de Lagrange aux points uvw
    call dLagrange1Dv(dMat=drMatrix,lx=li,dlx=duLi,transpose=.true.) !> duLi(1:np,1:nPt)=Transpose[drMatrix] li
    call dLagrange1Dv(dMat=dsMatrix,lx=li,dlx=dvLi,transpose=.true.) !> dvLi(1:np,1:nPt)=Transpose[dsMatrix] li
    call dLagrange1Dv(dMat=dtMatrix,lx=li,dlx=dwLi,transpose=.true.) !> dwLi(1:np,1:nPt)=Transpose[dtMatrix] li
    
    print '(/4x,"∂uli(:,uvw)")'
    do i=1,nPt
      print '(4x,"i=",i4," ∂uli=",14(f9.6,1x))',i,duLi(:,i)
    enddo
    print '(/4x,"∂vli(:,uvw)")'
    do i=1,nPt
      print '(4x,"i=",i4," ∂vli=",14(f9.6,1x))',i,dvLi(:,i)
    enddo
    print '(/4x,"∂wli(:,uvw)")'
    do i=1,nPt
      print '(4x,"i=",i4," ∂wli=",14(f9.6,1x))',i,dwLi(:,i)
    enddo
    
    open(unit=100,file="BasisP5_grad_ai.dat",action='write')
    write(100,'(4x,"duai")')
    do i=1,nPt
      write(100,'(4x,"i=",i3,2x,14(e19.12,1x))')i,duLi(:,i)
    enddo
    write(100,'(/4x,"dvai")')
    do i=1,nPt
      write(100,'(4x,"i=",i3,2x,14(e19.12,1x))')i,dvLi(:,i)
    enddo
    write(100,'(/4x,"dwai")')
    do i=1,nPt
      write(100,'(4x,"i=",i3,2x,14(e19.12,1x))'),i,dwLi(:,i)
    enddo
    close(100)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> comparaison des resultats (calcul direct et calcul interpolé)
    np =(ord+1)*(ord+2)*(2*ord+3)/6 !> = \sum_{k=1}^{ord+1} k^2
    print'(/4x,"np x nPt =",i6," x ",i6)',np,nPt
    
    !> test li
    cpt0=0 ; deltaMax=0d0
    do j=1,nPt
      f0=0d0
      do i=1,np
        f0=f0+li(i,j)*fi(i)
      enddo
      delta=abs(f(j)-f0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt0=cpt0+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"f(uvw)- ∑ ai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),f(j),f0,f(j)-f0
      endif
    enddo
    print '(4x,"erreur sur   f cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt0,nPt,deltaMax
    
    !> test duLi
    cpt1=0 ; deltaMax=0d0
    do j=1,nPt
      dxf0=0d0
      do i=1,np
        dxf0=dxf0+duLi(i,j)*fi(i)
      enddo
      delta=abs(dxf(j)-dxf0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt1=cpt1+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"   f(uvw)- ∑   ai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),f(j),f0,f(j)-f0
      endif
    enddo
    print '(4x,"erreur sur ∂uf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt1,nPt,deltaMax
    
    !> test dvLi
    cpt2=0 ; deltaMax=0d0
    do j=1,nPt
      dyf0=0d0
      do i=1,np
        dyf0=dyf0+dvLi(i,j)*fi(i)
      enddo
      delta=abs(dyf(j)-dyf0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt2=cpt2+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"∂vf(uvw)- ∑ ∂vai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),dyf(j),dyf0,dyf(j)-dyf0
      endif
    enddo
    print '(4x,"erreur sur ∂vf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt2,nPt,deltaMax
    
    !> test dwLi
    cpt3=0 ; deltaMax=0d0
    do j=1,nPt
      dzf0=0d0
      do i=1,np
        dzf0=dzf0+dwLi(i,j)*fi(i)
      enddo
      delta=abs(dzf(j)-dzf0) ; if( delta>deltaMax)deltaMax=delta
      if( delta>tol )then
        cpt3=cpt3+1
        print '(4x,"ad=",i6,2x,"uvw=",3(f12.5,1x),"∂wf(uvw)- ∑ ∂wai fi= ",e22.15," - ",e22.15,"  =  ",e22.15)',j,uvw(1:3,j),dzf(j),dzf0,dzf(j)-dzf0
      endif
    enddo
    print '(4x,"erreur sur ∂wf cpt=",i6,"/",i6,3x,"deltaMax=",e22.15)',cpt3,nPt,deltaMax
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cpt=cpt0+cpt1+cpt2+cpt3
    if( .not.cpt==0 )then
      print '("xxx Erreur Test sur points de Gauss ord=",i2)',ord
      stop
    else
      write(*,'("<<< Test sur points de Gauss ord=",i2)')ord
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw)
    deallocate(x,y,z,w,f)
    deallocate(li,duLi,dvLi,dwLi)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(vand)
  deallocate(drMatrix,dsMatrix,dtMatrix)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Contrôles Réussis"/)')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
  
contains
  
  subroutine fxyz(xyz,f,dxf,dyf,dzf)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8) :: xyz(:,:)
    real(8) ::   f(:)
    real(8) :: dxf(:)
    real(8) :: dyf(:)
    real(8) :: dzf(:)
    !>
    real(8) :: x,y,z
    integer :: i,n
    integer :: power
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !print '("fxyz ord=",i2)',ord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    n=size(xyz,2) ! print '("fxyz n=",i6)',n
    
    do i=1,n
      x=xyz(1,i)
      y=xyz(2,i)
      z=xyz(3,i)
      
      !> ord0
      f  (i)=3
      dxf(i)=0
      dyf(i)=0
      dzf(i)=0
      
      if( ord>=1 )then
        f  (i)=1.5*x+2.5*y-3.5*z +  f(i)
        dxf(i)=1.5               +dxf(i)
        dyf(i)=2.5               +dyf(i)
        dzf(i)=-3.5              +dzf(i)
      endif
      
      if( ord>=2 )then
        f  (i)= 7*x**2+4*y**2-5*z**2+x*y-x*z+3*y*z +  f(i)
        dxf(i)=14*x+y-z                            +dxf(i)
        dyf(i)=8*y+x+3*z                           +dyf(i)
        dzf(i)=-10*z-x+3*y                         +dzf(i)
      endif
      
      if( ord>=3 )then
        f  (i)= x**3+2.5*y**3-1.5*z**3+x**2*y+x*y*z +  f(i)
        dxf(i)=3*x**2+2*x*y+y*z                     +dxf(i)
        dyf(i)=7.5*y**2+x**2+x*z                    +dyf(i)
        dzf(i)=-4.5*z**2+x*y                        +dzf(i)
      endif
      
      if( ord>=4 )then
        f  (i)= 3*x**4+7*y**4-8*z**4-x**2*y**2+3*x*y**2*z +  f(i)
        dxf(i)=12*x**3-2*x*y**2+3*y**2*z                  +dxf(i)
        dyf(i)=28*y**3-2*x**2*y+6*x*y*z                   +dyf(i)
        dzf(i)=-32*z**3+3*x*y**2                          +dzf(i)
      endif
      
      if( ord>=5 )then
        f  (i)=3*x**5+7*x**2*y**2*z+z**5 +  f(i)
        dxf(i)=15*x**4+14*x*y**2*z       +dxf(i)
        dyf(i)=14*x**2*y*z               +dyf(i)
        dzf(i)=7*x**2*y**2+5*z**4        +dzf(i)
      endif
      
      if( ord>=6 )then
        f  (i)=3*x**6+7*x**2*y**2*z**2+x*y*z**4 +  f(i)
        dxf(i)=18*x**5+14*x*y**2*z**2+y*z**4    +dxf(i)
        dyf(i)=14*x**2*y*z**2+X*z**4            +dyf(i)
        dzf(i)=14*x**2*y**2*z+4*x*y*z**3        +dzf(i)
      endif
      
    enddo
    
    if( ord>=7 )then
      print '("warning ord>=6 not tested")'
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return  
  end subroutine fxyz
  
end subroutine pyramTestBasis

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
  integer              :: ord
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
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !call pyramidNodes   (ord=ord, uvw=uvw, display=.false.)  !> Points réguliers
  call pyramidNodesOpt(ord=ord, uvw=uvw, display=.false.)  !> Points optimises
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
  call pyramidVandermonde3D(ord=ord,a=a,b=b,c=c,vand=vand)
  call display    (title="    Vandermonde Matrix",mat=vand)
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
  call displaySparce(title="    Mass Matrix 0",mat=mass0,tol=1d-15)
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
  
  !> Connectivite
  !call testConnectivitiesOfSides()
  !call jacobiTest()
  !stop
  
  !> Test des quadratures
  !call testQuadratureGL()
  
  !> Test segments
 !call edge_00()
 !call edge_01()
  
  !> Test Triangles
 !call triangle_00()
 !call triangle_01()
  
  !> Test Tetra
  !call tetraTest()
  !!call tetraMaillageVisu() !> maillages de visu pour le tetra d'ordre élevé
  !call tetraMaillageVisuNew() ; stop
  
  !> Test Quad
 !call quadTest()
  
  !> Test pyramids
  call pyramBasis()
  !call pyramLebesgue()
  !call pyramMaillageVisu() !> maillages de visu pour la pyramide d'ordre élevé
  !call pyramDegreesOverSides()
  stop
  
  
  !> Tests des quadratures pour pyramides
  call pyramTestQuadrature()
  
  !> Tests des bases pyramides et de leurs derivees
  call pyramTestBasis()
  
  !> Test Matrice de Masse des pyramides
  call pyramTestMass()
  
  
end program main