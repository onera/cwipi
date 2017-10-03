module baseSimplexTools
  
  implicit none
  
  real(8), parameter   :: pi  =3.141592653589793238462643_8
  real(8), parameter   :: sqr2=1.414213562373095048801689_8
  real(8), parameter   :: sqr3=1.732050807568877293527446_8
  real(8), parameter   :: sqr6=2.449489742783178098197284_8
  real(8), parameter   :: sqr8=2.828427124746190097603377_8
  
  interface display       ; module procedure displayMatrix       ; end interface
  interface display       ; module procedure displayVector       ; end interface
  interface displaySparce ; module procedure displaySparceMatrix ; end interface
  
  interface mathematica ; module procedure mathematicaMatrix ; end interface
  
  interface massMatrix  ; module procedure massMatrix0 ; end interface
  interface massMatrix  ; module procedure massMatrix1 ; end interface
  
  contains
  
  subroutine gaussLegendreQuadratures(ord,xGL,wGL)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Cas particulier de jacobiGQ pour a=0 et b=0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer         , intent(in)  :: ord
    real(8), pointer, intent(out) :: xGL(:)
    real(8), pointer, intent(out) :: wGL(:)
    !>
    integer                       :: i,np
    real(8)                       :: h1(1:ord+1)
    real(8)                       :: diag0(1:ord+1)
    real(8)                       :: diag1(1:ord+1)
    real(8)                       :: mat  (1:ord+1,1:ord+1)
    real(8)                       :: eivec(1:ord+1,1:ord+1)
    real(8)                       :: y    (1:ord+1)
    integer                       :: cnt  (1:ord+1)
    integer                       :: rc
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=ord+1 ; allocate(xGL(1:np),wGL(1:np))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( ord==0 )then
      
      xGL(1)=0d0 ; wGL(1)=2d0
      
    else
      
      h1(1:np)=[( (2d0*i), i=0,ord )]
      diag0(1:np)=0d0
      
      do i=1,ord
        diag1(i+1)=2d0/(h1(i)+2d0)*sqrt( i*i*i*i/(h1(i)+1d0)/(h1(i)+3d0) )
      enddo
      
      mat(1:np,1:np)=0d0
      do i=1,np
        mat(i,i) = 1d0
      enddo
      
      call tql2(np,np,diag0,diag1,mat,rc)
      
      xGL(1:np)=diag0(1:np)
      where( abs(xGL)<1d-15 )xGL=0d0
      
      do i=1,np
        wGL(i)=mat(1,i)*mat(1,i)*2d0*gamma(1d0)
       !write(*,'(e22.15,3x,10e22.15)')xGL(i),mat(0:ord,i)
       !write(*,'(3x,"xGL=",e22.15,1x,"wGL=",e22.15)')xGauss(i),wGL(i)
      enddo
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef debug
    print '(/"Gauss Legendre Quadatures")'
    print '("i=",i2," xGL=",f12.5," wGL=",f12.5)',(i,xGL(i),wGL(i),i=1,size(xGL))  
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine gaussLegendreQuadratures
  
  subroutine gaussLobattoQuadratures(ord,xGL,wGL)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer         , intent(in)  :: ord
    real(8), pointer, intent(out) :: xGL(:)
    real(8), pointer, intent(out) :: wGL(:)
    !>
    integer                       :: i,k
    integer                       :: np
    real(8)                       :: coef
    real(8), pointer              :: p0(:)
    real(8), pointer              :: p1(:)
    real(8), pointer              :: pn(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( ord==0 )then
      
      allocate(xGL(1)) ; xGL(1)=0d0
      allocate(wGL(1)) ; wGL(1)=2d0
      
    elseif( ord==1 )then
      
      allocate(xGL(1:2)) ; xGL(1:2)=[-1d0,1d0]
      allocate(wGL(1:2)) ; wGL(1:2)=[ 1d0,1d0]
      
    else
      
      call gaussLegendreLobatto(ord=ord,xGLL=xGL)
      
      np=ord+1 ; allocate(wGL(1:np))
      allocate(p0(1:np)) ; p0(1:np)=1d0
      allocate(p1(1:np)) ; p1(1:np)=xGL(1:np)
      allocate(pn(1:np))
      do k=1,np-2
        pn(1:np)=( real(2*k+1,kind=8)*xGL(1:np)*p1(1:np) &
        &         -real(  k  ,kind=8)          *p0(1:np) &
        &        )/real(k+1,kind=8)
        
        p0(1:np)=p1(1:np)
        p1(1:np)=pn(1:np)
      enddo
      deallocate(p0,p1)
      
      coef=2d0/real(ord*(ord+1),kind=8)
      do k=1,np
        wGl(k)=coef/pn(k)**2
      enddo
      deallocate(pn)
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef debug
    print '(/"Gauss Legendre-Lobatto Quadatures")'
    print '("i=",i2," xGL=",f12.5," wGL=",f12.5)',(i,xGL(i),wGL(i),i=1,size(xGL))
    !stop
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine gaussLobattoQuadratures
  
  
  subroutine jacobiGQ(a,b,ord,nGauss,xGL,wGL)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Quadratures de Gauss - Polynômes de Jacobi
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    implicit none
    real(8) , intent(in)  :: a
    real(8) , intent(in)  :: b
    integer , intent(in)  :: ord
    !
    integer , intent(out) :: nGauss
    real(8) , intent(out) :: xGL(0:)
    real(8) , intent(out) :: wGL(0:)
    !>
    integer              :: i
    real(8)              :: amb,apb,a2mb2
    real(8)              :: h1(0:ord)
    real(8)              :: diag0(0:ord)
    real(8)              :: diag1(0:ord)
    real(8)              :: mat  (0:ord,0:ord)
    real(8)              :: eivec(0:ord,0:ord)
    real(8)              :: y    (0:ord)
    integer              :: cnt  (0:ord)
    integer              :: rc
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    nGauss=ord+1
    !
    apb=a+b
    amb=a-b
    a2mb2=a*a-b*b
    
    if( ord==0 )then
      
      xGL(0)=amb/(apb+2d0) ; wGL(0)=2d0
      
    else
      
      h1=[( (2d0*i+apb), i=0,ord )]
      
      if( a2mb2==0d0 )then
        do i=0,ord
          diag0(i)=0d0
        enddo
      else
        do i=0,ord
          diag0(i)=-5d-1*a2mb2/(h1(i)+2)/h1(i)
        enddo
      endif
      
      do i=1,ord
        diag1(i)=2d0/(h1(i-1)+2d0)*sqrt( i*(i+apb)*(i+a)*(i+b)/(h1(i-1)+1d0)/(h1(i-1)+3d0) )
      enddo
      
      mat=0d0
      do i=0,ord
        mat(i,i) = 1d0
      enddo
      
      call tql2(ord+1,ord+1,diag0,diag1,mat,rc)
      
      
      xGL=diag0
      where( abs(xGL)<1d-15 )xGL=0d0
      
     !write(*,*)'rc=',rc
      do i=0,ord
        wGL(i)=mat(0,i)**2 * 2**(apb+1)/(apb+1)*gamma(a+1d0)*gamma(b+1d0)/gamma(apb+1d0)
       !write(*,'(e22.15,3x,10e22.15)')xGL(i),mat(0:ord,i)
       !write(*,'(3x,"xGL=",e22.15,1x,"wGL=",e22.15)')xGL(i),wGL(i)
      enddo
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine jacobiGQ
    
  subroutine displayMatrix(title,mat)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*) :: title
    real(8)       :: mat(:,:)
    !>
    integer       :: i,j
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '(/a,2x,"size=",i4," x",i4)',trim(title),size(mat,1),size(mat,2)
    if( size(mat,2)<11 )then
      print '(3x,$)'
      do i=1,size(mat,1)
        do j=1,size(mat,2)
          print '(1x,e22.15,$)',mat(i,j)
        enddo
        print '()'
      enddo
    else
      print '(3x,$)'
      do i=1,size(mat,1)
        do j=1,size(mat,2)
          print '(1x,e12.5,$)',mat(i,j)
        enddo
        print '()'
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine displayMatrix
  
  subroutine displaySparceMatrix(title,mat,tol)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*) :: title
    real(8)       :: mat(:,:)
    real(8)       :: tol
    !>
    integer       :: i,j
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '(/a,2x,"size="i3," x",i3)',trim(title),size(mat,1),size(mat,2)
    do i=1,size(mat,1)
      do j=1,size(mat,2)
        if( abs(mat(i,j))>tol )print '(1x,i6,1x,i10,2x,e22.15)',i,j,mat(i,j)
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine displaySparceMatrix
  
  subroutine mathematicaMatrix(title,mat)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*) :: title
    real(8)       :: mat(:,:)
    !>
    integer       :: i,j
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '(/a)',trim(title)
    print '("{",$)'
    do i=1,size(mat,1)
      print '("{",$)'
      do j=1,size(mat,2)
        if( j< size(mat,2) )print '(1x,f22.15,","$)',mat(i,j)
        if( j==size(mat,2) )print '(1x,f22.15,   $)',mat(i,j)
      enddo
      if( i< size(mat,1) )print '("},",$)'
      if( i==size(mat,1) )print '("}}",$)'
    enddo
    print '(";")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine mathematicaMatrix
  
  subroutine displayVector(title,vec)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*)  :: title
    real(8)       :: vec(:)
    !>
    integer       :: i
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '(/a,2x,"size=",i10)',trim(title),size(vec)
    print '(4x,e22.15)',(vec(i),i=1,size(vec))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine displayVector
  
  function gamma(x) result(ga)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Purpose: Compute the gamma function ‚(x)
    ! Input :  x  --- Argument of ‚(x)
    !               ( x is not equal to 0,-1,-2,˙˙˙ )
    ! Output:  GA --- ‚(x)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    implicit none
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8)            :: x
    real(8)            :: ga
    !>
    integer            :: k,m1,m
    real(8)            :: r,gr
    real(8), parameter :: pi=3.141592653589793238462643_8
    real(8)            :: z
    !>
    real(8), parameter :: g(26)=[ 1.0d0               ,&
    &                             0.5772156649015329d0,&
    &                            -0.6558780715202538d0,&
    &                            -0.420026350340952d-1,&
    &                             0.1665386113822915d0,&
    &                            -0.421977345555443d-1,&
    &                            -0.96219715278770d-2 ,&
    &                             0.72189432466630d-2 ,&
    &                            -0.11651675918591d-2 ,&
    &                            -0.2152416741149d-3  ,&
    &                             0.1280502823882d-3  ,&
    &                            -0.201348547807d-4   ,&
    &                            -0.12504934821d-5    ,&
    &                             0.11330272320d-5    ,&
    &                            -0.2056338417d-6     ,&
    &                             0.61160950d-8       ,&
    &                             0.50020075d-8       ,&
    &                            -0.11812746d-8       ,&
    &                             0.1043427d-9        ,&
    &                             0.77823d-11         ,&
    &                            -0.36968d-11         ,&
    &                             0.51d-12            ,&
    &                            -0.206d-13           ,&
    &                            -0.54d-14            ,&
    &                             0.14d-14            ,&
    &                             0.1d-15              ]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( x==int(x) )then
      if( x>0d0 )then
        ga=1d0
        m1=x-1
        do k=2,m1
          ga=ga*real(k,kind=8)
        enddo
      else
        ga=301d0
      endif
    else
      if( abs(x)>1d0 )then
        z=abs(x)
        m=int(z)
        r=1d0
        do k=1,m
          r=r*(z-k)
        enddo
        z=z-m
      else
        z=x
      endif
      gr=g(26)
      do k=25,1,-1
        gr=gr*z+g(k)
      enddo
      ga=1d0/(gr*z)
      if( dabs(x)>1d0 )then
        ga=ga*r
        if( x<0d0 )ga=-pi/(x*ga*dsin(pi*x))
      endif
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end function gamma
  
  subroutine eigenVectors(mat,w,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: mat(:,:)
    real(8), intent(out), pointer :: w(:)
    logical, intent(in)           :: display
    !
    integer                       :: ad
    integer                       :: n,lda
    real(8), pointer              :: eigVector(:,:)
    integer                       :: lWork
    integer, pointer              :: ipiv(:)
    real(8), pointer              :: work(:)
    integer                       :: iErr,iFail
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(mat,1) ; lWork=64*n ; allocate(work(lWork))
    allocate(w(n))
    lda=n
    allocate(eigVector(n,lda)) ; eigVector(:,:)=mat(:,:)
    call dsyev('V', 'U', n, eigVector, lda, w, work, lWork, iErr)
    deallocate(work)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      print '(4x,"eigv(",i4,")=",e22.15)',(ad,w(ad),ad=1,size(w))
      print '(/4x,"cond(mass)=",e22.15/)',maxval(w)/minval(w)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine eigenVectors
    
  subroutine jacobiP(alpha,beta,n,u,jf)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> IN  jf is not allocated
    !> OUT jf is     allocated
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8)           :: alpha
    real(8)           :: beta
    integer           :: n
    real(8)           :: u(:)
    real(8) , pointer :: jf(:)
    !-
    integer           :: k
    integer           :: nVert
    real(8)           :: anew,bnew,aold
    real(8)           :: h1
    real(8)           :: apb,amb,a2mb2
    real(8)           :: gamma0,gamma1,gammai
    real(8) , pointer :: jf0(:), jf1(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nVert=size(u) ; allocate(jf(1:nVert))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    apb=alpha+beta
    gamma0=2d0**(apb+1d0)*gamma(alpha+1d0)*gamma(beta+1d0) /( (apb+1d0)*gamma(apb+1d0) )
    if( n==0 )then
      jf(:)=1d0/sqrt(gamma0)
    else if( n==1 )then
      gamma1=(alpha+1d0)*(beta+1d0)/(apb+3d0)*gamma0
      amb=alpha-beta
      
      jf(:)=5d-1/sqrt(gamma1)*(amb+(apb+2d0)*u(:))
    else
      allocate(jf0(1:nVert),jf1(1:nVert))
      gamma1=(alpha+1d0)*(beta+1d0)/(apb+3d0)*gamma0
      amb=alpha-beta
      
      jf0(:)=1d0/sqrt(gamma0)
      jf1(:)=5d-1/sqrt(gamma1)*(amb+(apb+2d0)*u(:))
      
      a2mb2=alpha*alpha-beta*beta
      aold=2d0/(apb+2d0)*sqrt((alpha+1d0)*(beta+1d0)/(apb+3d0))
      do k=2,n
        h1=2d0*(k-1)+apb
        anew=2d0/(h1+2d0)*sqrt( k*(k+apb)*(k+alpha)*(k+beta)/((h1+1d0)*(h1+3d0)) )
        bnew=-a2mb2/h1/(h1+2d0)
        jf(:)=1d0/anew*(-aold*jf0(:)+(u(:)-bnew)*jf1(:) )
        jf0(:)=jf1(:)
        jf1(:)=jf (:)
        aold =anew
      enddo
      deallocate(jf0,jf1)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine jacobiP
  
  subroutine dJacobiP(alpha,beta,n,u,jf)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> IN  jf is not allocated
    !> OUT jf is     allocated
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8)           :: alpha
    real(8)           :: beta
    integer           :: n
    real(8)           :: u(:)
    real(8) , pointer :: jf(:)
    !--
    integer           :: nVert
    integer           :: i
    real(8) , pointer :: jac(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nVert=size(u) ; allocate(jf(1:nVert))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( n==0 )then
      jf(:)=0d0
    else
      call JacobiP(alpha=alpha+1d0,beta=beta+1d0,n=n-1,u=u,jf=jac)
      jf(:)=sqrt(n*(n+alpha+beta+1d0))*jac(:)
      deallocate(jac)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine dJacobiP
  
  
  subroutine jacobiPtab(alpha,beta,n,u,jf)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8)           :: alpha
    real(8)           :: beta
    integer           :: n
    real(8)           :: u (  :)
    real(8) , pointer :: jf(:,:)
    !-
    integer           :: k
    integer           :: nVert
    real(8)           :: anew,bnew,aold
    real(8)           :: h1
    real(8)           :: apb,amb,a2mb2
    real(8)           :: gamma0,gamma1,gammai
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    nVert=size(u) ; allocate(jf(1:nVert,0:n))
    !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    apb=alpha+beta
    gamma0=2**(apb+1d0)*gamma(alpha+1d0)*gamma(beta+1d0) /( (apb+1d0)*gamma(apb+1d0) )
    if( n==0 )then
      jf(:,0)=1d0/sqrt(gamma0)
    else if( n==1 )then
      gamma1=(alpha+1d0)*(beta+1d0)/(apb+3d0)*gamma0
      amb=alpha-beta
      !
      jf(:,0)=1d0/sqrt(gamma0)
      jf(:,1)=5d-1/sqrt(gamma1)*(amb+(apb+2d0)*u(:))
    else
      gamma1=(alpha+1d0)*(beta+1d0)/(apb+3d0)*gamma0
      amb=alpha-beta
      !
      jf(:,0)=1d0/sqrt(gamma0)
      jf(:,1)=5d-1/sqrt(gamma1)*(amb+(apb+2d0)*u(:))
      !
      a2mb2=alpha*alpha-beta*beta
      aold=2d0/(apb+2d0)*sqrt((alpha+1d0)*(beta+1d0)/(apb+3d0))
      do k=2,n
        h1=2d0*(k-1)+apb
        anew=2d0/(h1+2d0)*sqrt( k*(k+apb)*(k+alpha)*(k+beta)/((h1+1d0)*(h1+3d0)) )
        bnew=-a2mb2/h1/(h1+2d0)
        jf(:,k)=1d0/anew*(-aold*jf(:,k-2)+(u(:)-bnew)*jf(:,k-1) )
        aold =anew
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine jacobiPtab
  
  subroutine gradJacobiPtab(alpha,beta,n,u,jf)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8)           :: alpha
    real(8)           :: beta
    integer           :: n
    real(8)           :: u(   :)
    real(8) , pointer :: jf(:,:)
    !--
    integer           :: nVert
    integer           :: i
    real(8) , pointer :: j(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nVert=size(u) ; allocate(jf(1:nVert,0:n))
    !
    if( n==0 )then
      jf(:,0)=0d0
    else
      jf(:,0)=0d0
      call jacobiPtab(alpha=alpha+1d0,beta=beta+1d0,n=n-1,u=u,jf=j)
      do i=1,n
        jf(:,i)=sqrt(i*(i+alpha+beta+1d0))*j(:,i-1)
      enddo
      deallocate(j)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine gradJacobiPtab  
    
  subroutine jacobi(u,alpha,beta,n,jf)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)           :: u(:)      ! [in]
    integer, intent(in)           :: n         ! [in]
    real(8), intent(in)           :: alpha     ! [in]
    real(8), intent(in)           :: beta      ! [in]
    real(8), intent(out), pointer :: jf(:)     ! [out]
    !---
    real(8) , pointer             :: jf0(:), jf1(:)
    integer                       :: i=0,k=0
    integer                       :: len
    real(8)                       :: apb,amb,a2mb2
    real(8)                       :: a1,a2,a3,a4
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> preliminaires
    len=size(u) ; allocate(jf(len))
    apb=alpha+beta
    amb=alpha-beta
    a2mb2=alpha*alpha-beta*beta
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    select case(n)
    case(0)
      jf(:)=sqrt( gamma(apb+2d0)/(2d0**(apb+1d0)*gamma(alpha+1d0)*gamma(beta+1d0)) )
    case(1 )
      jf(:)=sqrt( gamma(apb+2d0)/(2d0**(apb+1d0)*gamma(alpha+1d0)*gamma(beta+1d0)) )
      jf(:)=5d-1*jf(:)*sqrt((apb+3d0)/((alpha+1d0)*(beta+1d0)) )*((apb+2d0)*u(:)+amb)
    case(2:)
      
      allocate(jf0(len)) ; jf0(:)=sqrt( gamma(apb+2d0)/(2d0**(apb+1d0)*gamma(alpha+1d0)*gamma(beta+1d0)) )
      allocate(jf1(len)) ; jf1(:)=5d-1*jf0(:)*sqrt((apb+3d0)/((alpha+1d0)*(beta+1d0)) )*((apb+2d0)*u(:)+amb)
      
      do k=1,n-1
        jf (:)= ( (u(:)-b(k))*jf1(:)-a(k)*jf0(:) )/a(k+1)
        jf0(:)=jf1(:)
        jf1(:)=jf (:)
      enddo
      
      deallocate(jf0,jf1)
      
    end select
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
    
  contains
    
    function a(n) result (an)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      integer, intent(in) :: n
      real(8)             :: an
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      an= 2d0/(2d0*real(n,kind=8)+apb)         &
      &                                        &
      &  *sqrt(                                &
      &          real(n,kind=8)                &
      &        *(real(n,kind=8)+apb)           &
      &        *(real(n,kind=8)+alpha)         &
      &        *(real(n,kind=8)+beta)          &
      &        /(                              &
      &           (2d0*real(n,kind=8)+apb-1d0) &
      &          *(2d0*real(n,kind=8)+apb+1d0) &
      &          )                             &
      &        )
      return
    end function
    
    function b(n) result (bn)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      integer, intent(in) :: n
      real(8)             :: bn
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      bn=a2mb2/((2d0*real(n,kind=8)+apb)*(2d0*real(n,kind=8)+apb+2d0))
      return
    end function
    
  end subroutine jacobi
  
  subroutine dJacobi(n,alpha,beta,u, jf)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)           :: u(:)      ! [in]
    integer, intent(in)           :: n         ! [in]
    real(8), intent(in)           :: alpha     ! [in]
    real(8), intent(in)           :: beta      ! [in]
    real(8), intent(out), pointer :: jf(:)     ! [out]
    !
    integer                       :: len
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> preliminaires
    len=size(u) ; allocate(jf(len))   
    if( n==0 )then
      jf(1:len)=0d0
    else
      call jacobi(n=n-1,alpha=alpha+1d0,beta=beta+1d0,u=u, jf=jf) ! J_{n-1}^{a+1,b+1}(u)
      jf(1:len)=sqrt( real(n,kind=8)*(real(n+1,kind=8)+alpha+beta) )*jf(1:len)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine dJacobi
  
  
  subroutine gaussLegendreLobatto(ord,xGLL)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: xGLL(:)
    !---
    real(8), pointer              :: vand(:,:)
    integer                       :: i,j,k,np
    real(8), pointer              :: x0(:)
    real(8), parameter            :: eps=1d-15
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=ord+1
    allocate(xGLL(1:np))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( ord==0 )then
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      xGLL(1)=0d0
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
    else
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      allocate(vand(1:np,1:np)) ; vand(1:np,1:np)=0d0
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> Use the Chebyshev-Gauss-Lobatto nodes as the first guess
      xGLL(1:np)=[(cos(pi*real(i-1,kind=8)/real(ord,kind=8)),i=np,1,-1)]
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> Newton-Raphson iteration
      allocate(x0(1:np)) ; x0(1:np)=2d0
      do while( maxval(abs(xGLL(1:np)-x0(1:np)))>eps )
        x0(1:np)=xGLL(1:np)
        
        !> Polynomes de Legendre :
        !> p_{0}(x)= 1 && p_{1}(x)= x pour k>= 2 : (k+1) p_{k+1} = (2k+1) p_{k} x - k p_{k-1}
        vand(1:np,1)=1d0
        vand(1:np,2)=xGLL(1:np)
        
        do k=1,np-2
          vand(1:np,k+2)=( real(2*k+1,kind=8)*xGLL(1:np)*vand(1:np,k+1) &
          &               -real(  k  ,kind=8)           *vand(1:np,k  ) &
          &              )/real(k+1,kind=8)
        enddo
        
        xGLL(1:np)= x0(1:np)                            &
        &          -( xGLL(1:np)       *vand(1:np,np  ) &
        &            -                  vand(1:np,np-1) &
        &           )/( real(np,kind=8)*vand(1:np,np  ) )
        
      enddo
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      deallocate(vand,x0)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine gaussLegendreLobatto
  
  subroutine warpFactor(ord,xnodes,xout, warp)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! l_ord^i(xout) = \Prod_{j=0;j/=i}^{ord} (xout-xeq_j)/(xeq_i-xeq_j) est un polynôme de Lagrange sur grille equidistante
    ! warp   (xout) =    \sum_{i=0}^{ord} (xnodes_i-xeq_i) l_ord^i(xout)
    !
    ! warp'  (xout) = warp(xout)/(1-xout^2)   if xout|=-1 || xout|=+1
    !               = warp(xout)              if xout==-1 || xout==+1
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: xnodes(:)
    real(8), intent(in) , pointer :: xout  (:)
    real(8), intent(out), pointer :: warp(:)
    !--
    integer                       :: i,j
    integer                       :: iVert,nVert
    real(8), pointer              :: d   (:)
    real(8), pointer              :: xeq (:)
    real(8)                       :: var
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !> Dimensionnement
    nVert=size(xout) ! print '(/"nVert=",i6)',nVert
    
    !> Grille points equidistants dans [-1,1]
    allocate(xeq(1:ord+1)) ; xeq(1:ord+1)=[( (-1d0+2d0*real(i-1,kind=8)/real(ord,kind=8)), i=1,ord+1)]
    
    !> warp(xout) = \sum_{i=1}^{ord+1} (xnodes_i-xeq_i) \Prod_{j=1;j/=i}^{ord+1} (xout-xeq_j)/(xeq_i-xeq_j)
    allocate(warp(1:nVert)) ; warp(1:nVert)=0d0
    allocate(d   (1:nVert))
    do i=1,ord+1
      var=xnodes(i)-xeq(i)
      d(1:nVert)=var
      do j=1,ord+1
        if( .not.i==j )then
          var=1d0/(xeq(i)-xeq(j))
          d(1:nVert)=d(1:nVert)*(xout(1:nVert)-xeq(j))*var
        endif
      enddo
      warp(1:nVert)=warp(1:nVert)+d(1:nVert)
    enddo
    
    !> warp'  (xout) = warp(xout)/(1-xout^2)   if xout|=-1 || xout|=+1
    !>               = warp(xout)              if xout==-1 || xout==+1
    where( xout(1:nVert)/=-1d0 .and. xout(1:nVert)/=+1d0 )
      warp(:)=warp(:)/(1-xout(:)*xout(:))
    end where
    
    deallocate(d,xeq)
    
    return
    
    !> Idem mais un peu plus rapide
    !> warp(xout) = \sum_{i=0}^{ord} (xnodes_i-xeq_i) \Prod_{j=0;j/=i}^{ord} (xout-xeq_j)/(xeq_i-xeq_j)
    do i=0,ord
      var=xnodes(i)-xeq(i)
      d(1:nVert-2)=var
      
      do j=1,ord-1
        if( .not.i==j )then
          var=1d0/(xeq(i)-xeq(j))
          d(1:nVert-2)=d(1:nVert-2) * (xout(1:nVert-2)-xeq(j)) * var
        endif
      enddo
      
      !> deflate end roots
      if( .not.i==  0 )d(1:nVert-2)=-d(1:nVert-2)/(xeq(i)-xeq(  0))
      if( .not.i==ord )d(1:nVert-2)= d(1:nVert-2)/(xeq(i)-xeq(ord))
      
      warp(1:nVert-2)=warp(1:nVert-2)+d(1:nVert-2)
    enddo
    
    deallocate(d,xeq)
    
    return
  end subroutine warpFactor

  function matrixIsSymetric(mat) result(test)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), pointer :: mat(:,:)
    logical          :: test
    !>
    integer          :: i,j,ad,n
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    test=.true.
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not. size(mat,1)==size(mat,2) )then
      test=.false.
      return
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(mat,1)
    do i=1,n
      do j=i,n
        if( .not.mat(i,j)==mat(j,i) )then
          test=.false.
          return
        endif
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end function matrixIsSymetric
  
  subroutine compactForm(mat0,mat1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> [a11 a12 a13]
    !> [    a22 a23]  => [a11 a12 a22 a23 a33]
    !> [        a33]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)   , pointer :: mat0(:,:)
    real(8), intent(inout), pointer :: mat1(:)
    !>
    integer                         :: i,j,ad,n
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( matrixIsSymetric(mat0) )then
      n=size(mat0,1)
      ad=0
      do i=1,n
        do j=i,n
          ad=ad+1
          mat1(ad)=mat0(i,j)
        enddo
      enddo
    else
      print '("Matrix is not Symetric => Not compacted")'
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine compactForm
  
  subroutine massMatrix0(vand,mass)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! mass = Inverse[Vand.Transpose[Vand]]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: vand(:,:)
    real(8), intent(out), pointer :: mass(:,:)
    !---
    real(8), pointer              :: ap(:)
    integer                       :: ad
    real(8)                       :: t0,t1
    integer                       :: i,j,k,n
    integer                       :: iErr,iFail
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef debug
    print '("baseSimplexTools:massMatrix0")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef perfo
    call cpu_time(t0)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( associated(vand) )then
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      n=size(vand,1)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> ap=vand.Transpose[vand]
      !> [a11 a12 a13]
      !> [    a22 a23]  => [a11 a12 a22 a23 a33]
      !> [        a33]
      allocate(ap(n*(n+1)/2)) ; ap(:)=0d0
      do i=1,n
        do j=i,n
          ad=i+(j-1)*j/2
          do k=1,n
            ap(ad)=ap(ad)+vand(i,k)*vand(j,k)
          enddo
        enddo
      enddo
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> Cholevski factorization
      call dpptrf('U',n,ap,iErr)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> inversion of Symmetric Matrix
      call dpptri('U',n,ap,iErr)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      allocate(mass(n,n))
      do i=1,n
        do j=i,n
          ad=i+(j-1)*j/2
          mass(i,j)=ap(ad)
          mass(j,i)=ap(ad)
        enddo
      enddo
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      deallocate(ap)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
    else !> vand is not associated => ord=0 => mass = \int_{-1}^{+1} dV = 2
      
      allocate(mass(1,1)) ; mass(1,1)=2d0
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef perfo
    call cpu_time(t1)
    print '("massMatrix Dt=",f12.9)',t1-t0
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine massMatrix0
  
  subroutine massMatrix1(vand, mass)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define massMatrix1 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! mass = Inverse[Vand.Transpose[Vand]]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)   , pointer :: vand(:,:)
    real(8), intent(inout), pointer :: mass(  :)
    !>
    integer                         :: ad
    real(8)                         :: t0,t1
    integer                         :: i,j,k,n
    integer                         :: iErr,iFail
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if massMatrix1==1
    print '(">>> baseSimplexTools:massMatrix1")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef perfo
    call cpu_time(t0)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( associated(vand) )then
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      n=size(vand,1)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> Inv(mass)=vand.Transpose[vand]
      !> [a11        ]
      !> [a21 a22    ]  => [a11 a21 a31 a22 a32 a33]
      !> [a31 a32 a33]
     !allocate(mass(1:n*(n+1)/2)) ; mass(:)=0d0
      mass(:)=0d0
      ad=0
      do j=1,n
        do i=j,n
          ad=i+(j-1)*(2*n-j)/2
          do k=1,n
            mass(ad)=mass(ad)+vand(i,k)*vand(j,k)
          enddo
        enddo
      enddo
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> Cholevski factorization
      call dpptrf('L',n,mass,iErr)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> inversion of Symmetric Matrix
      call dpptri('L',n,mass,iErr)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
    else !> ord=0 => vand is not associated
      
      allocate(mass(1:1)) ; mass(1:1)=2d0
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef perfo
    call cpu_time(t1)
    print '("massMatrix Dt=",f12.9)',t1-t0
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if massMatrix1==1
    print '("<<< baseSimplexTools:massMatrix1")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef massMatrix1
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine massMatrix1
  
  subroutine invMassMatrix(vand,invMass)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! invMass = Vand.Transpose[Vand]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: vand   (:,:)
    real(8), intent(out), pointer :: invMass(:,:)
    !---
    real(8), pointer              :: ap(:)
    integer                       :: ad
    real(8) :: t0,t1
    integer                       :: i,j,k,n
    integer                       :: iErr,iFail
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef perfo
    call cpu_time(t0)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(vand,1)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> mat=vand.Transpose[vand]
    allocate(invMass(n,n)) ; invMass(:,:)=0d0
    do i=1,n
      do j=1,n
        do k=1,n
          invMass(i,j)=invMass(i,j)+vand(i,k)*vand(j,k)
        enddo
      enddo
    enddo
    !call displayMatrix(title="inverseMassMatrix",mat=invMass)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef perfo
    call cpu_time(t1)
    print '("invMassMatrix Dt=",f12.9)',t1-t0
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine invMassMatrix  
  
  subroutine inverseMatrix(mat,invMat)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> invMat = Inverse[mat]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in ), pointer ::    mat(:,:)
    real(8), intent(out), pointer :: invMat(:,:)
    !---
    integer                       :: n
    integer                       :: lWork
    integer, pointer              :: ipiv(:)
    real(8), pointer              :: work(:)
    integer                       :: iErr,iFail
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(mat,1) ; allocate(ipiv(n))
    lWork=64*n     ; allocate(work(lWork))
    allocate(invMat(n,n)) ; invMat(:,:)=mat(:,:)
    call dgetrf(n,n,invMat(1,1),n,ipiv(1),iErr)
    call dgetri(n,invMat(1,1),n,ipiv(1),work(1),lWork,iErr)
    deallocate(ipiv,work)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine inverseMatrix
  
  
  subroutine derive1D(vand,dVand,dMat)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> dMat = dVand.Inverse[vand]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in ), pointer ::  vand(:,:)
    real(8), intent(in ), pointer :: dVand(:,:)
    real(8), intent(out), pointer :: dMat (:,:)
    !---
    real(8), pointer              :: iVand(:,:) ! inverse(vand)
    integer                       :: i,j,k,nMod
    integer                       :: lWork
    integer, pointer              :: ipiv(:)
    real(8), pointer              :: work(:)
    integer                       :: iErr,iFail
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> mat= dVand.Inverse[vand]
    nMod=size(vand,1) ; allocate(ipiv(1:nMod))
    allocate(iVand(nMod,nMod)) ; iVand(1:nMod,1:nMod)=vand(1:nMod,1:nMod)
    call dgetrf(nMod,nMod,iVand(1,1),nMod,ipiv(1),iErr)
    
    if( .not.iErr==0 )then
      print '("XXX baseSimplexTools:derive1D:dgetrf iErr=",i10)',iErr
      
      if( iErr> 0 )then
        print '(/4x,"iErr = ",i3," > 0, U(i,i) is exactly zero. The factorization")',iErr
        print '( 4x,"has been completed, but the factor U is exactly")'
        print '( 4x,"singular, and division by zero will occur if it is used")'
        print '( 4x,"to solve a system of equations.")'
        print '( 4x,"vand(",i3,",",i3,")=",e22.15)',iErr,iErr,vand(iErr,iErr)
        print '()'
      else
        print '(/4x,"iErr < 0:  if iErr = -i, the i-th argument had an illegal value")'
        print '()'
      endif
      
      call displayMatrix(title="vand",mat=vand)
      call displayMatrix(title="iVand",mat=iVand)
      
      print '(/4x,"Mathematica Matrix Vand")'
      call mathematicaMatrix(title="vand=",mat=vand)
      print '()'
      
      stop '("stop @ baseSimplexTools:derive1D:dgetrf unsuccessful exit")'
    endif
    
    
    lWork=64*nMod ; allocate(work(1:lWork))
    call dgetri(nMod,iVand(1,1),nMod,ipiv(1),work(1),lWork,iErr)    
    if( .not.iErr==0 )stop '("stop @ baseSimplexTools:derive1D:dgetri unsuccessful exit")'
    deallocate(ipiv,work)
   !call displayMatrix(title="Inverse Vandermonde Matrix",mat=iVand)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> dMat = dVand.Transpose[vand]
    allocate(dMat(1:nMod,1:nMod))
    dMat(1:nMod,1:nMod)=0d0
    do i=1,nMod
      do j=1,nMod
        do k=1,nMod
          dMat(i,j)=dMat(i,j)+dvand(i,k)*iVand(k,j)
        enddo
      enddo
    enddo
    deallocate(iVand)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine derive1D
  
  subroutine dLagrange1Dv(dMat,lx,dlx,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! dLagrange1D := Transpose[dMat_{ij}] l_j(r);
    ! transpose = .false. => dlx(1:nPt,1:np)
    ! transpose = .true.  => dlx(1:np,1:nPt)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)    , pointer :: dMat(:,:)
    real(8), intent(in)    , pointer :: lx  (:,:)
    real(8), intent(inout) , pointer :: dlx (:,:)
    logical, intent(in)              :: transpose
    !-
    integer                          :: i,j,k
    integer                          :: nMod,nNod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( transpose )then
      nMod=size(lx,1) ; nNod=size(lx,2)
     !allocate(dlx(1:nMod,1:n)) ; dlx(1:n,1:nMod)=0d0
      dlx(:,:)=0d0
      do i=1,nNod
        do j=1,nMod
          do k=1,nMod
            dlx(j,i)=dlx(j,i)+dMat(k,j)*lx(k,i)
          enddo
        enddo
      enddo
    else
      nMod=size(lx,2) ; nNod=size(lx,1)
     !allocate(dlx(1:n,1:nMod)) ; dlx(1:n,1:nMod)=0d0
      dlx(:,:)=0d0
      do i=1,nNod
        do j=1,nMod
          do k=1,nMod
            dlx(i,j)=dlx(i,j)+dMat(k,j)*lx(i,k)
          enddo
        enddo
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine dLagrange1Dv
  
  subroutine lebesgue(lx,l,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in ), pointer :: lx (:,:)
    real(8), intent(out), pointer :: l  (:,:)
    logical, intent(in)           :: transpose
    !--
    integer                       :: i,nNod,nMod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !print '(">>> subroutine lebesgue")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( transpose )then
      nMod=size(lx,1) ;nNod=size(lx,2)
     !allocate(l(1,1:nNod)) ; l(1,1:nNod)=0d0
      l(1,1:nNod)=0d0
      do i=1,nMod
        l(1,1:nNod)=l(1,1:nNod)+abs( lx(i,1:nNod) )
      enddo
    else
      nMod=size(lx,2) ; nNod=size(lx,1)
     !allocate(l(1:nNod,1)) ; l(1:nNod,1)=0d0
      l(1:nNod,1)=0d0
      do i=1,nMod
        l(1:nNod,1)=l(1:nNod,1)+abs( lx(1:nNod,i) )
      enddo    
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !print '("<<< subroutine lebesgue")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine lebesgue
  
  
SUBROUTINE TQL2(NM,N,D,E,Z,IER)
!-------------------------------------------------------------------------
!     QL METHOD TO DETERMINE THE EIGENVALUES AND EIGENVECTORS OF:
!
!       1)  A SYMMETRIC TRIDIAGONAL MATRIX.
!       2)  A FULL SYMMETRIC MATRIX AFTER A PREVIOUS CALL TO TRED2.
!
!     CALLING MODE:
!               CALL TQL2(NM,N,D,E,Z,IER)
!     INPUTSS:
!     NM  (I4)  1ST DIMENSION OF MATRICES A AND Z IN CALLING PROGRAM
!     N   (I4)  SIZE OF Z
!     D  (R*8)  MAIN DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     E  (R*8)  SUB-DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     Z  (R*8)  TABLE (NM,N) STORING THE UNITY MATRIX IF THE TRIDIAGONAL
!               MATRIX IS DEFINED BY D AND E, CASE #1.
!               FOR CASE #2, IT CONTAINS THE ELEMENTS OF THE TRANSFORMATION
!               MATRIX AFTER A CALL TO TRED2.
!     OUTPUTS:
!     D  (R*8)  EIGENVALUES
!     Z  (R*8)  EIGENVECTORS
!     IER (I4)  ERROR CODE = 0,  CONVERGENCE OK.
!                          = L,  NO CONVERGENCE FOR THE Lth EIGENVALUE
!
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.
!-------------------------------------------------------------------------
      REAL(8) ::  D(N),E(N),Z(NM,N),B,C,F,G,H,P,R,S,EPS,EPS1
      INTEGER :: I,J,K,L,M,N,NM,JM,IER
      DATA EPS /0.D0/,JM /30/
      IER = 0
      IF (N==1) GO TO 38
!
!     MACHINE EPSILON
!
      IF (EPS.NE.0.D0) GO TO 12
      EPS = 1.D0
   10 EPS = EPS/2.D0
      EPS1 = 1.D0+EPS
      IF (EPS1.GT.1.D0) GO TO 10
!
   12 DO 14 I = 2,N
   14 E(I-1) = E(I)
      E(N) = 0.D0
      F = 0.D0
      B = 0.D0
!
      DO 28 L = 1,N
      J = 0
      H = EPS*(ABS(D(L))+ABS(E(L)))
      IF (B.LT.H) B = H
!
!     SEEK SMALLEST ELEMENT OF SUBDIAGONAL
!
      DO 16 M = L,N
      IF (ABS(E(M)).LE.B) GO TO 18
   16 CONTINUE
   18 IF (M.EQ.L) GO TO 26

!     START ITERATION

   20 IF (J.EQ.JM) GO TO 36
      J = J+1

!     SHIFT

      G = D(L)
      P = (D(L+1)-G)/(2.D0*E(L))
      R = SQRT(P*P+1.D0)
      D(L) = E(L)/(P+SIGN(R,P))
      H = G-D(L)
      DO 22 I = L+1,N
   22 D(I) = D(I)-H
      F = F+H

!     QL TRANSFORMATION

      P = D(M)
      C = 1.D0
      S = 0.D0
      DO 24 I = M-1,L,-1
      G = C*E(I)
      H = C*P
      IF (ABS(P).GE.ABS(E(I))) THEN
        C = E(I)/P
        R = SQRT(C*C+1.D0)
        E(I+1) = S*P*R
        S = C/R
        C = 1.D0/R
      ELSE
        C = P/E(I)
        R = SQRT(C*C+1.D0)
        E(I+1) = S*E(I)*R
        S = 1.D0/R
        C = C*S
      ENDIF
      P = C*D(I)-S*G
      D(I+1) = H+S*(C*G+S*D(I))

!     ELEMENTS OF EIGENVECTORS

      DO 24 K = 1,N
        H = Z(K,I+1)
        Z(K,I+1) = S*Z(K,I)+C*H
        Z(K,I) = Z(K,I)*C-S*H
   24 CONTINUE
      E(L) = S*P
      D(L) = C*P
      IF (ABS(E(L)).GT.B) GO TO 20

!     CONVERGENCE

   26 D(L) = D(L)+F
   28 CONTINUE

!     SORT EIGENVALUES AND EIGENVECTORS
!     IN ASVENDING ORDER

      DO 34 L = 2,N
      I = L-1
      K = I
      P = D(I)
      DO 30 J = L,N
      IF (D(J).GE.P) GO TO 30
      K = J
      P = D(J)
   30 CONTINUE
      IF (K.EQ.I) GO TO 34
      D(K) = D(I)
      D(I) = P
      DO 32 J = 1,N
      P = Z(J,I)
      Z(J,I) = Z(J,K)
   32 Z(J,K) = P
   34 CONTINUE
      GO TO 38

!     NO CONVERGENCE

   36 IER = L
   38 RETURN
END SUBROUTINE TQL2



      SUBROUTINE TRED2(NM,N,A,D,E,Z)
!---------------------------------------------------------------------------
!     TRIDIAGONALIZATION OF A SYMMETRIC MATRIX BY ORTHOGONAL TRANSFORMATIONS
!     (ALGORITHM OF HOUSEHOLDER)
!     CALLING MODE:
!               CALL TRED2(NM,N,A,D,E,Z)
!     INPUTS:
!     NM  (I4)  1ST DIMENSION OF MATRICES A AND Z IN CALLING PROGRAM
!     N   (I4)  SIZE OF A
!     A  (R*8)  TABLE(NM,N) STORING THE COEFFICIENTS OF SYMMETRIC A MATRIX
!               (LOWER HALF), A IS NOT DESTROYED DURING THE PROCESS
!               IF Z MATRIX HAS NOT THE SAME ADDRESS.
!     OUTPUTS:
!     D  (R*8)  MAIN DIAGONAL (N) OF REDUCED TRIDIAGONAL MATRIX
!     E  (R*8)  SUB-DIAGONAL (N) OF REDUCED TRIDIAGONAL MATRIX
!     Z  (R*8)  TABLE (NM,N) STORING THE ELEMENTS OF THE ORTHOGONAL 
!               TRANSFORMATION MATRIX.
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.
!-----------------------------------------------------------------------
      INTEGER I,J,K,L,N,NM
      REAL *8 A(NM,N),D(N),E(N),Z(NM,N),F,G,H,HH,SCALE

!     LOWER HALF OF A PUT INTO Z

      DO 10 I = 1,N
      DO 10 J = 1,I
   10 Z(I,J) = A(I,J)
      IF (N.EQ.1) GO TO 32

!     N-2 STAGE OF TRANSFORMATION

      DO 30 I = N,2,-1
      L = I-1
      H = 0.

!     CONDITIONNING BY NORM OF A

      SCALE = 0.
      IF (L.LT.2) GO TO 14
      DO 12 K = 1,L
   12 SCALE = SCALE+ABS(Z(I,K))
      IF (SCALE.NE.0.) GO TO 16

   14 E(I) = Z(I,L)
      GO TO 28

   16 DO 18 K = 1,L
      Z(I,K) = Z(I,K)/SCALE
      H = H+Z(I,K)*Z(I,K)
   18 CONTINUE

      F = Z(I,L)
      G = -SIGN(SQRT(H),F)
      E(I) = SCALE*G
      H = H-F*G
      Z(I,L) = F-G
      F = 0.
      DO 24 J = 1,L
      Z(J,I) = Z(I,J)/H
      G = 0.

!     ELEMENT OF A*U
      DO 20 K = 1,J
   20 G = G+Z(J,K)*Z(I,K)
      IF (L.GE.J+1) THEN
      DO 22 K = J+1,L
   22 G = G+Z(K,J)*Z(I,K)

!     ELEMENT OF P = A*U/H

      END IF
      E(J) = G/H
      F = F+E(J)*Z(I,J)
   24 CONTINUE

!     ELEMENT OF K

      HH = F/(H+H)

!     REDUCED FORM OF A

      DO 26 J = 1,L
      F = Z(I,J)
      G = E(J)-HH*F
      E(J) = G
      DO 26 K = 1,J
      Z(J,K) = Z(J,K)-F*E(K)-G*Z(I,K)
   26 CONTINUE
!
   28 D(I) = H
   30 CONTINUE

!     END OF TRANSFORMATION

   32 D(1) = 0.
      E(1) = 0.

!     ACCUMULATE TRANSFORMATION MATRICES IN Z

      DO 40 I = 1,N
      L = I-1
      IF (D(I).NE.0.) THEN
      DO 36 J = 1,L
      G = 0.
      DO 34 K = 1,L
   34 G = G+Z(I,K)*Z(K,J)
      DO 36 K = 1,L
      Z(K,J) = Z(K,J)-G*Z(K,I)
   36 CONTINUE
      END IF
      D(I) = Z(I,I)
      Z(I,I) = 1.
      IF (L.LT.1) GO TO 40
      DO 38 J = 1,L
      Z(I,J) = 0.
      Z(J,I) = 0.
   38 CONTINUE
   40 CONTINUE

      RETURN
      END SUBROUTINE TRED2  
  
  
end module baseSimplexTools