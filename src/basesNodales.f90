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
  int=0d0
  do i=1,order+1
    int=int+wGJ(i)*func( xGJ(i) )
  enddo
  
  print '("\int_{-1}^{+1} 2 - 3x + x^2 + 2x^3 - 6x^4 + 8x^5- 19x^6 + x^7 dx = ",e22.15)',int
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  deallocate(xGL,wGL)
  deallocate(xGLL,wGLL)
  deallocate(xGJ,wGJ)
  
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


subroutine test1D()
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
  print '(/"size(xout)=",i3)',size(xout)
  print '("xout(",i2,")=",f12.5)',(i,xout(i),i=1,nVert)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Polynomes d'interpolation
  call lagrange1Dv(ord=order,vand=vand,x=xout,lx=lxout,transpose=.true.) ! true pour affichage
  call display(title="lxout",mat=lxout)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Dérivée Polynomes d'interpolation
  call gradVandermonde1D(ord=order,a=r,dVand=dVand)
  call display(title="drVand Matrix",mat=dVand)
  
  call derive1D(vand=vand,dVand=dVand,dMat=dr)
  call display(title="Dr Matrix",mat=dr)
  
  call dLagrange1Dv(dMat=dr,lx=lxout,dlx=dlxout,transpose=.true.) ! <= true pour affichage
  call display(title="d lx=Dr^T lx",mat=dlxout)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call lebesgue(lx=lxout,l=leb,transpose=.false.)  
  print '("i=",i2," xout=",f12.5," lebesgue=",f12.5)',(i,xout(i),leb(i,1),i=1,nVert)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  deallocate(xout,lxout,dlxout,leb)
  
  return
  
end subroutine test1D

subroutine test2D()
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
  real(8), pointer   :: mode(:,:)
  real(8), pointer   :: xyout(:,:),lxout(:,:),drLxout(:,:),dsLxout(:,:)
  real(8), pointer   :: leb(:,:)
  real(8), pointer   :: uvw(:,:),rs(:,:),a(:),b(:)
  real(8), parameter :: eps=1d-15
  real(8), pointer   :: eigv(:)
  character(80)      :: fileName
  real(8)            :: node_xy(2,3) ! Triangle
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
  call nodes2Duv2rs(uv=uvw,rs=rs ,display=.true.) ! rs(1:2,:)=2d0*uv(1:2,:)-1d0
  call nodes2Drs2ab(rs=rs,a=a,b=b,display=.true.) ! a=2 (1+r)/(1-s)-1 && b=s
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
  &    point_num=size(uvw,2),&
  &    point_xy=rs          ,&
  &    point_show=1          ) ! 2, show the points and number them
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call vandermonde2D(ord=order,a=a,b=b,vand=vand)
  if( order<10 )then
    call display(title="Vandermonde = [P_j(xi_i)]",mat=vand)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Polynomes d'interpolation (on teste avec les points d'interpolation
  call lagrange2Dv(ord=order,vand=vand,a=a,b=b,lx=lxout,transpose=.true.) ! true pour affichage
  if( order<10 )then
    call display(title="Test avec l(uvw)",mat=lxout)
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
  
  call simplex2D(ord=order,a=a,b=b,mode=mode,transpose=.false.)
  
  call lagrange2Dv(ord=order,vand=vand,a=a,b=b,lx=lxout,transpose=.false.)
  call dLagrange1Dv(dMat=drMatrix,lx=lxout,dlx=drLxout,transpose=.false.)
  call dLagrange1Dv(dMat=dsMatrix,lx=lxout,dlx=dsLxout,transpose=.false.)
  
  if( order<6 )then
    call writeSolOut2D(title="simplex2D"   ,solOut=mode   )
    call writeSolOut2D(title="lagrange2D"  ,solOut=lxout  )
    call writeSolOut2D(title="drLagrange2D",solOut=drLxout)
    call writeSolOut2D(title="dsLagrange2D",solOut=dsLxout)
  endif
  
  call lebesgue(lx=lxout,l=leb,transpose=.false.)
  call writeSolOut2D(title="lebesgue2D",solOut=leb)
  
  deallocate(mode)
  deallocate(a,b)
  deallocate(xyout)
  deallocate(lxout,drLxout,dsLxout)
  deallocate(leb)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine test2D

subroutine test1D_01()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use baseSimplex1D
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
end subroutine test1D_01

subroutine test2D_01()
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
    call nodes2Duv2rs(uv=uvw,rs=rs  ,display=.false.)  ! rs(1:2,:)=2d0*uv(1:2,:)-1d0
    call nodes2Drs2ab(rs=rs,a=a,b=b ,display=.false.) ! a=2 (1+r)/(1-s)-1 && b=s
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
    call nodes2Duv2rs(uv=uvw,rs=rs  ,display=.false.)  ! rs(1:2,:)=2d0*uv(1:2,:)-1d0
    call nodes2Drs2ab(rs=rs,a=a,b=b ,display=.false.) ! a=2 (1+r)/(1-s)-1 && b=s
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
end subroutine test2D_01


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
  integer              :: i,j,order,nVert,ad,np,nPt,cpt
  real(8), pointer     :: vand(:,:),dVand(:,:),jf(:,:),dr(:,:)
  real(8), pointer     :: drVand  (:,:),dsVand  (:,:),dtVand  (:,:)
  real(8), pointer     :: drMatrix(:,:),dsMatrix(:,:),dtMatrix(:,:)
  real(8), pointer     :: mass(:,:)
  real(8), pointer     :: mode(:,:)
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
  
  call nodes3Duvw2rst(uvw=uvw,rst=rst) ! rst(1:3,:)=2d0*uvw(1:3,:)-1d0
 !write(*,'(/"Tetra (rst):")')
 !print '("rst(1:3,",i2,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,rst(1:3,ad),ad=1,size(rst,2))
  
  call nodes3Drst2abc(rst=rst,a=a,b=b,c=c)
 !write(*,'(/"Tetra (abc):")')
 !print '("a,b,c(",i2,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,a(ad),b(ad),c(ad),ad=1,size(a))
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
  call lagrange3Dv(ord=order,vand=vand,a=a,b=b,c=c,lx=lxOut,transpose=.true.) ! true pour affichage
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
  call nodes3Drst2abc(rst=xyzOut,a=a,b=b,c=c)
  
  call simplex3D  (ord=order,a=a,b=b,c=c,mode=mode,transpose=.false.)           !> Psi(xyzOut)
  call lagrange3Dv(ord=order,vand=vand,a=a,b=b,c=c,lx=lxOut,transpose=.false.)  !> lxOut= Inverse[Transpose[Vand]].Psi[xyzOut] lxOut(nPt,np)
  if( order<3 )then
    call writeSolOut3D(title="simplex3D"//sfx,solOut=mode   )
    call writeSolOut3D(title="lagrange3D"//sfx  ,solOut=lxOut  )
  endif
  deallocate(mode)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des dérivées des fonctions de Lagrange aux points xyzOut
  call dLagrange1Dv(dMat=drMatrix,lx=lxOut,dlx=drLxOut,transpose=.false.) !> drLxOut(1:nPt,1:np)= Transpose[drMatrix] lxOut
  call dLagrange1Dv(dMat=dsMatrix,lx=lxOut,dlx=dsLxOut,transpose=.false.) !> dsLxOut= Transpose[dsMatrix] lxOut
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
  nPt=size(lxOut,1)
  np =(order+1)*(order+2)*(order+3)/6 ! =size(lxOut,2)
  allocate(fi(1:np)) ; fi(1:np)=1d0
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions lxOut : f(xyzOut) = \sum_{i=1,np} lxOut(xyzOut,i) f_i
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
  !> Test des fonctions drLxOut : df/dt(xyzOut) = \sum_{i=1,np} dtLxOut(xyzOut,i) f_i
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
  call lebesgue(lx=lxout,l=leb,transpose=.false.)
  call writeSolOut3D(title="lebesgue3DP"//sfx,solOut=leb)
  
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


subroutine testConnectivitiesOfSides()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use baseSimplex2D, only: edgesConnectivity
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !
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
  call massMatrix(vand=vand,mass=mass)
  call display(title="Mass=",vec=mass)
  deallocate(mass)
  !> Factorisation de la matrice
  !call factorise(n=ord+1,mat=mass)
  !call display(title="Factorized(Mass)=",vec=mass)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  call massMatrix(vand=vand,mass=massL2)
  call display(title="MassL2=",mat=massL2)
  
  n=size(massL2,1) ; allocate(mass(n*(n+1)/2)) ; mass(:)=0d0
  call compactForm(mat0=massL2,mat1=mass)
  call display(title="MassL2=",vec=mass)
  deallocate(mass)
  
  !> Produit tensoriel
  allocate(massQ4((ord+1)**2,(ord+1)**2)) ; massQ4(:,:)=0d0
  do i=1,ord+1
    row0=(i-1)*(ord+1)
    do j=1,ord+1
      col0=(j-1)*(ord+1)
      do k=1,ord+1
        do l=1,ord+1
         !print '("row0+k=",i3,2x,"col0+l=",i3)',row0+k,col0+l
          massQ4(row0+k,col0+l)=massL2(i,j)*massL2(k,l) ! Tenseur massQ4(i,j,k,l)
        enddo
      enddo
    enddo
  enddo
  
 !call display(title="MassQ4=",mat=massQ4)
  call compactForm(mat0=massQ4,mat1=mass)
  call display(title="MassQ4=",vec=mass)
  call factorise(n=(ord+1)**2,mat=mass)
  call display(title="Factorized(MassQ4)=",vec=mass)
  
  return
end subroutine quadTest


subroutine pyramBasis()
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
  
  real(8), parameter   :: eps=1d-12
  real(8), pointer     :: a(:),b(:),c(:)
  real(8), pointer     :: vand(:,:),dVand(:,:) !,jf(:,:),dr(:,:)
  real(8), pointer     :: drVand  (:,:),dsVand  (:,:),dtVand  (:,:)
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
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Polynomes d'interpolation (on teste avec les points d'interpolation) => Matrice Identité
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
  call pyramidGradVandermonde3D(ord=ord,a=a,b=b,c=c,drVand=drVand,dsVand=dsVand,dtVand=dtVand)
  call derive1D(vand=vand,dVand=drVand,dMat=drMatrix) !> drMatrix = drVand.Inverse[vand]
  call derive1D(vand=vand,dVand=dsVand,dMat=dsMatrix) !> dsMatrix = dsVand.Inverse[vand]
  call derive1D(vand=vand,dVand=dtVand,dMat=dtMatrix) !> dtMatrix = dtVand.Inverse[vand]
  if( ord<3 )then
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
  deallocate(a,b,c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramBasis

subroutine pyramDegreesOverSides()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use basePyramid
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer              :: ord
  integer, allocatable :: sides(:)
  integer              :: sidesIdx(1:6)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramidSides3D(ord=ord, sidesIdx=sidesIdx, sides=sides, display=.true.)
  deallocate(sides)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramDegreesOverSides

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
  real(8), pointer     :: drVand  (:,:),dsVand  (:,:),dtVand  (:,:)
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
  call pyramidGradVandermonde3D(ord=ord,a=a,b=b,c=c,drVand=drVand,dsVand=dsVand,dtVand=dtVand)
  call derive1D(vand=vand,dVand=drVand,dMat=drMatrix) !> drMatrix = drVand.Inverse[vand]
  call derive1D(vand=vand,dVand=dsVand,dMat=dsMatrix) !> dsMatrix = dsVand.Inverse[vand]
  call derive1D(vand=vand,dVand=dtVand,dMat=dtMatrix) !> dtMatrix = dtVand.Inverse[vand]
  if( ord<3 )then
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
  deallocate(a,b,c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
  if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
  if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des fonctions de Lagrange aux points xyzOut
  call pyramidReadXYZout3D(xyzOut=xyzOut, display=.true.)
  call pyramiduvw2abc(uvw=xyzOut,a=a,b=b,c=c)
  
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
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer            :: ord
  real(8), pointer   :: uvw(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Construction maillage Pyramid P_i")')
  write(*,'("Warning ghs3d is required")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !call pyramidNodes   (ord=ord, uvw=uvw, display=.true.)
  call pyramidNodesOpt(ord=ord, uvw=uvw, display=.true.)  !> Points optimises
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramidSkin3D(ord=ord, uvw=uvw, display=.true.)
  call pyramidMesh3D(ord=ord, uvw=uvw, display=.true.)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramMaillageVisu


program main
  
  !> Connectivite 
 !call testConnectivitiesOfSides()
 !call jacobiTest()
  
  !> Test des quadratures
 !call testQuadratureGL()
  
  !> Test segments
 !call test1D()
 !call test1D_01()
  
  !> Test Triangles
 !call test2D()
 !call test2D_01()
  
  !> Test Tetra
  !call tetraTest()
  !call tetraMaillageVisu() !> maillages de visu pour le tetra d'ordre élevé
  
  !> Test Quad
  !call quadTest()
  
  !> Test pyramids
 !call pyramBasis()
 !call pyramLebesgue()
 !call pyramMaillageVisu() !> maillages de visu pour la pyramide d'ordre élevé
 call pyramDegreesOverSides()
  
end program main