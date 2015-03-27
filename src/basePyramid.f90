module basePyramid
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use ieee_arithmetic
  use iso_fortran_env
  use baseSimplexTools
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>  Pyramid ord
  !>  nn=(ord+1)*(ord+2)*(2*ord+3)/6
  !>  ne=3*ord*ord+2
  !>  ni=(ord-1)*(ord-2)*(2*ord-3)/6 !> sans les faces
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  implicit none
  
  interface pyramiduvw2abc ; module procedure pyramiduvw2abc_0 ; end interface
  interface pyramiduvw2abc ; module procedure pyramiduvw2abc_1 ; end interface
  interface pyramiduvw2abc ; module procedure pyramiduvw2abc_2 ; end interface
  
contains
  
  subroutine pyramiduvw2abc_0(u,v,w, a,b,c)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> a = u/(1-w) pour w=1, u=0
    !> b = v/(1-w) pour w=1  v=0
    !> c = w
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)            :: u   ,v   ,w
    real(8), intent(out) , pointer :: a(:),b(:),c(:)
    !>
    real(8)                        :: w0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(a(1),b(1),c(1))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( w==1d0 )then
      a(1)=0d0
      b(1)=0d0
      c(1)=1d0
    else
      w0=1d0/(1d0-w)
      a(1)=u*w0
      b(1)=v*w0
      c(1)=w
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramiduvw2abc_0
  
  subroutine pyramiduvw2abc_1(u,v,w, a,b,c)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> a = u/(1-w) pour w=1, u=0
    !> b = v/(1-w) pour w=1  v=0
    !> c = w
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)  , pointer :: u(:),v(:),w(:)
    real(8), intent(out) , pointer :: a(:),b(:),c(:)
    !>
    integer                        :: i,n
    real(8)                        :: w0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(u)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(a(1:n),b(1:n),c(1:n))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do i=1,n
      if( w(i)==1d0 )then
        a(i)=0d0
        b(i)=0d0
        c(i)=1d0
      else
        w0=1d0/(1d0-w(i))
        a(i)=u(i)*w0
        b(i)=v(i)*w0
        c(i)=w(i)
      endif
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramiduvw2abc_1
  
  subroutine pyramiduvw2abc_2(uvw,a,b,c)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> a = u/(1-w) pour w=1, u=0
    !> b = v/(1-w) pour w=1  v=0
    !> c = w
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)  , pointer :: uvw(:,:)
    real(8), intent(out) , pointer :: a(:),b(:),c(:)
    !>
    integer                        :: i,n
    real(8)                        :: w0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(uvw,2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(a(1:n),b(1:n),c(1:n))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do i=1,n
      if( uvw(3,i)==1d0 )then
        a(i)=0d0
        b(i)=0d0
        c(i)=1d0
      else
        w0=1d0/(1d0-uvw(3,i))
        a(i)=uvw(1,i)*w0
        b(i)=uvw(2,i)*w0
        c(i)=uvw(3,i)
      endif
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramiduvw2abc_2
  
  subroutine pyramidLagrange3Dv(ord,vand,a,b,c,lx,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define pyramidLagrange3Dv 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! lagrange3Dv := Inverse[Transpose[Vand]].Psi[x];
    ! transpose = .true.  => lx(1:ord+1,1:nPt)
    ! transpose = .false. => lx(1:nPt,1:ord+1)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)            :: ord
    real(8), intent(in)  , pointer :: vand(:,:)
    real(8), intent(in)  , pointer :: a(:),b(:),c(:)
    real(8), intent(out) , pointer :: lx  (:,:)
    logical, intent(in)            :: transpose
    !---
    integer                        :: i,j,k,nPt,np
    real(8)                        :: gamma(0:ord+1)
    integer                        :: iOrd
    real(8), pointer               :: mode1(:),mode2(:),mode3(:)
    real(8), pointer               :: Psi(:,:),mode(:,:)
    real(8), pointer               :: mat(:,:)
    integer                        :: lWork
    integer, pointer               :: ipiv(:)
    real(8), pointer               :: work(:)
    integer                        :: iErr
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidLagrange3Dv==1
    print '(">>> pyramidLagrange3Dv size(vand)=",i3,"x",i3)',size(vand,1),size(vand,2)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)*(2*ord+3)/6 ; nPt=size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidLagrange3Dv==1
    print '("    pyramidLagrange3Dv step1 np=",i10," nPt=",i10)',np,nPt
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !call pyramidBasePi(ord=ord,a=a,b=b,c=c,mode=psi,transpose=.false.) !> psi(abc,:)
    call pyramidBasePi(ord=ord,a=a,b=b,c=c,mode=psi,transpose=.false.) !> psi(abc,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidLagrange3Dv==1
    print '("    pyramidLagrange3Dv step2")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> mat=Transpose[Vand]
    allocate(mat(1:np,1:np))
    do i=1,np
      do j=1,np
        mat(i,j)=vand(j,i)
      enddo
    enddo
    
    !> mat=Inverse[ Transpose[Vand] ]
    lWork=64*(np) ; allocate(work(lWork),ipiv(np))
    call dgetrf(np,np,mat(1,1),np,ipiv(1),iErr)
    call dgetri(np,mat(1,1),np,ipiv(1),work(1),lWork,iErr)
    deallocate(ipiv,work)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidLagrange3Dv==1
    print '("    pyramidLagrange3Dv step3")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> lx = Inverse[Transpose[Vand]].Psi
    if( transpose )then
      allocate(lx(1:np,1:nPt)) ; lx(:,:)=0d0
      !> psi(1:nPt,1:np)
      
      do i=1,nPt
        do j=1,np
          do k=1,np
            lx(j,i)=lx(j,i)+mat(j,k)*Psi(i,k)  !> Attention de bien prendre Psi(i,k)
          enddo
        enddo
      enddo
      
    else
      
      allocate(lx(1:nPt,1:np)) ; lx(:,:)=0d0
      ! psi(1:nPt,1:np)
      do i=1,nPt
        do j=1,np
          do k=1,np
            lx(i,j)=lx(i,j)+mat(j,k)*Psi(i,k)  !> Attention de bien prendre Psi(i,k)
          enddo
        enddo
      enddo
    endif
   !call displayMatrix(title="lx",mat=lx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidLagrange3Dv==1
    print '("    pyramidLagrange3Dv step4")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(mat)
    deallocate(psi)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidLagrange3Dv==1
    print '("<<< pyramidLagrange3Dv")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef pyramidLagrange3Dv
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidLagrange3Dv
  
  
  subroutine pyramidVandermonde3D(ord,a,b,c,vand)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:)
    real(8), intent(in) , pointer :: b(:)
    real(8), intent(in) , pointer :: c(:)
    real(8), intent(out), pointer :: vand(:,:)
    !>
    real(8)             , pointer :: weight(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(a) .or. .not.associated(b) .or. .not.associated(c) )then
      print '("3D Nodes a or/and b or/and c not associated")'
      stop "@pyramidVandermonde3D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Vandermonde matrix
   !call pyramidBasePi(ord=ord,a=a,b=b,c=c,mode=vand,transpose=.false.) !> vand(1:nPoint,1:np)
    call pyramidBasePi(ord=ord,a=a,b=b,c=c,mode=vand,transpose=.false.) !> vand(1:nPoint,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidVandermonde3D
  
  
  subroutine pyramidGradVandermonde3D(ord,a,b,c,drMode,dsMode,dtMode)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define pyramidGradVandermonde3D 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:),b(:),c(:)
    real(8), intent(out), pointer :: drMode(:,:),dsMode(:,:),dtMode(:,:)
    !>
    integer                       :: i,j
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(a) .or. .not.associated(b) .or. .not.associated(c) )then
      print '("3D Nodes a or/and b or/and c not associated")'
      stop "@pyramidGradVandermonde3D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidGradVandermonde3D==1
    print '(">>> pyramidGradVandermonde3D")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Grad Vandermonde matrix
    call pyramidGradBasePi( &
    &    ord=ord           ,&
    &    a=a,b=b,c=c       ,&
    &    drMode=drMode     ,&
    &    dsMode=dsMode     ,&
    &    dtMode=dtMode     ,&
    &    transpose=.false.  )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidGradVandermonde3D==1
    print '(/"    pyramidGradVandermonde3D: ∂Psi/∂x(1:n ,1:nP",i1,")")',ord
    do i=1,size(drMode,2)
      print '(4x,"i=",i3,2x,"∂Psi/∂x=",14(f9.6,1x))',i,drMode(i,:)
    enddo
    print '(/"    pyramidGradVandermonde3D: ∂Psi/∂y(1:n ,1:nP",i1,")")',ord
    do i=1,size(dsMode,2)
      print '(4x,"i=",i3,2x,"∂Psi/∂z=",14(f9.6,1x))',i,dsMode(i,:)
    enddo
    print '(/"    pyramidGradVandermonde3D: ∂Psi/∂z(1:n ,1:nP",i1,")")',ord
    do i=1,size(dtMode,2)
      print '(4x,"i=",i3,2x,"∂Psi/∂z=",14(f9.6,1x))',i,dtMode(i,:)
    enddo
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidGradVandermonde3D==1
    print '("<<< pyramidGradVandermonde3D")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef pyramidGradVandermonde3D
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidGradVandermonde3D
  
  subroutine pyramidBaseP1(uvw, ai, transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: uvw(:,:)
    real(8), intent(out), pointer :: ai(:,:)
    logical, intent(in)           :: transpose
    !>
    integer                       :: i,nn
    real(8)                       :: w0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => ai(1:np,1:nPt)
    !> Transpose = False => ai(1:nPt,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nn=size(uvw,2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(      transpose )allocate( ai(1:5,size(uvw,2)) )
    if( .not.transpose )allocate( ai(size(uvw,2),1:5) )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> vertex functions
    if( transpose )then
      do i=1,nn
        if( uvw(3,i)==1d0 )then
          ai(1,i)=0d0
          ai(2,i)=0d0
          ai(3,i)=0d0 !> retournement 3<->4
          ai(4,i)=0d0 !> retournement 3<->4
          ai(5,i)=1d0
        else
          w0=uvw(1,i)*uvw(2,i)/(1d0-uvw(3,i))
          ai(1,i) = 0.25d0*(1d0-uvw(1,i)-uvw(2,i)-uvw(3,i)+w0)
          ai(2,i) = 0.25d0*(1d0+uvw(1,i)-uvw(2,i)-uvw(3,i)-w0)
          ai(4,i) = 0.25d0*(1d0+uvw(1,i)+uvw(2,i)-uvw(3,i)+w0) !> retournement 3<->4
          ai(3,i) = 0.25d0*(1d0-uvw(1,i)+uvw(2,i)-uvw(3,i)-w0) !> retournement 3<->4
          ai(5,i) = uvw(3,i)
        endif
      enddo
    else
      do i=1,nn
        if( uvw(3,i)==1d0 )then
          ai(i,1)=0d0
          ai(i,2)=0d0
          ai(i,3)=0d0 !> retournement 3<->4
          ai(i,4)=0d0 !> retournement 3<->4
          ai(i,5)=1d0
        else
          w0=uvw(1,i)*uvw(2,i)/(1d0-uvw(3,i))
          ai(i,1) = 0.25d0*(1d0-uvw(1,i)-uvw(2,i)-uvw(3,i)+w0)
          ai(i,2) = 0.25d0*(1d0+uvw(1,i)-uvw(2,i)-uvw(3,i)-w0)
          ai(i,4) = 0.25d0*(1d0+uvw(1,i)+uvw(2,i)-uvw(3,i)+w0) !> retournement 3<->4
          ai(i,3) = 0.25d0*(1d0-uvw(1,i)+uvw(2,i)-uvw(3,i)-w0) !> retournement 3<->4
          ai(i,5) = uvw(3,i)
        endif
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( transpose )then
      if( 0==1 )then
        do i=1,size(ai,2)
          print '(i6,1x,"uvw=",3(f12.9,1x),3x,"ai=",5(f12.9,1x))',i,uvw(1:3,i),ai(1:5,i)
        enddo
      endif
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidBaseP1
  
  subroutine pyramidGradBaseP1(uvw,duai,dvai,dwai,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: uvw(:,:)
    real(8), intent(out), pointer :: duai(:,:),dvai(:,:),dwai(:,:)
    logical, intent(in)           :: transpose
    !>
    integer                       :: i,nn
    real(8)                       :: u0,v0,w0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => ai(1:np,1:nPt)
    !> Transpose = False => ai(1:nPt,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nn=size(uvw,2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(      transpose )allocate( duai(1:5,size(uvw,2)),dvai(1:5,size(uvw,2)),dwai(1:5,size(uvw,2)) )
    if( .not.transpose )allocate( duai(size(uvw,2),1:5),dvai(size(uvw,2),1:5),dwai(size(uvw,2),1:5) )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> vertex functions
    if( transpose )then
      
      do i=1,nn
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if( uvw(3,i)==1d0 )then
          u0=0d0
          v0=0d0
          w0=0d0 !????
        else
          w0=1d0/(1d0-uvw(3,i))      !> w0=1d0/(1d0-w)
          u0=uvw(1,i)*w0             !> u0=u  /(1d0-w)
          v0=uvw(2,i)*w0             !> v0=  v/(1d0-w)
          w0=uvw(1,i)*uvw(2,i)*w0*w0 !> w0=u*v/(1d0-w)**2
        endif
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        duai(1,i)=.25d0*(-1d0+v0)
        duai(2,i)=.25d0*( 1d0-v0)
        duai(4,i)=.25d0*( 1d0+v0) !> retournement 3<->4
        duai(3,i)=.25d0*(-1d0-v0) !> retournement 3<->4
        duai(5,i)= 0d0
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        dvai(1,i)=.25d0*(-1d0+u0)
        dvai(2,i)=.25d0*(-1d0-u0)
        dvai(4,i)=.25d0*( 1d0+u0) !> retournement 3<->4
        dvai(3,i)=.25d0*( 1d0-u0) !> retournement 3<->4
        dvai(5,i)= 0d0
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        dwai(1,i)=.25d0*(-1d0+w0)
        dwai(2,i)=.25d0*(-1d0-w0)
        dwai(4,i)=.25d0*(-1d0+w0) !> retournement 3<->4
        dwai(3,i)=.25d0*(-1d0-w0) !> retournement 3<->4
        dwai(5,i)= 1d0
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
      enddo
    else
      do i=1,nn
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if( uvw(3,i)==1d0 )then
          u0=0d0
          v0=0d0
          w0=0d0 !????
        else
          w0=1d0/(1d0-uvw(3,i))      !> w0=1d0/(1d0-w)
          u0=uvw(1,i)*w0             !> u0=u  /(1d0-w)
          v0=uvw(2,i)*w0             !> v0=  v/(1d0-w)
          w0=uvw(1,i)*uvw(2,i)*w0*w0 !> w0=u*v/(1d0-w)**2
        endif
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        duai(i,1)=.25d0*(-1d0+v0)
        duai(i,2)=.25d0*( 1d0-v0)
        duai(i,4)=.25d0*( 1d0+v0) !> retournement 3<->4
        duai(i,3)=.25d0*(-1d0-v0) !> retournement 3<->4
        duai(i,5)= 0d0
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        dvai(i,1)=.25d0*(-1d0+u0)
        dvai(i,2)=.25d0*(-1d0-u0)
        dvai(i,4)=.25d0*( 1d0+u0) !> retournement 3<->4
        dvai(i,3)=.25d0*( 1d0-u0) !> retournement 3<->4
        dvai(i,5)= 0d0
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        dwai(i,1)=.25d0*(-1d0+w0)
        dwai(i,2)=.25d0*(-1d0-w0)
        dwai(i,4)=.25d0*(-1d0+w0) !> retournement 3<->4
        dwai(i,3)=.25d0*(-1d0-w0) !> retournement 3<->4
        dwai(i,5)= 1d0
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidGradBaseP1
  
  subroutine pyramidH6BaseP1(uvw, ai, transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Hexa [-1,1] x [-1,1] [0,1] dégénéré
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: uvw(:,:)
    real(8), intent(out), pointer :: ai(:,:)
    logical, intent(in)           :: transpose
    !>
    integer                       :: i,nn
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> b1=1/4 (1-u) (1-v) (1-w)
    !> b2=1/4 (1+u) (1-v) (1-w)
    !> b3=1/4 (1+u) (1+v) (1-w)
    !> b4=1/4 (1-u) (1+v) (1-w)
    !> b5=                   w
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => ai(1:np,1:nPt)
    !> Transpose = False => ai(1:nPt,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nn=size(uvw,2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> vertex functions
    if( transpose )then
      allocate( ai(1:5,size(uvw,2)) )
      do i=1,nn
        ai(1,i) = 0.25d0*(1d0-uvw(1,i))*(1d0-uvw(2,i))*(1d0-uvw(3,i))
        ai(2,i) = 0.25d0*(1d0+uvw(1,i))*(1d0-uvw(2,i))*(1d0-uvw(3,i))
        ai(4,i) = 0.25d0*(1d0+uvw(1,i))*(1d0+uvw(2,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        ai(3,i) = 0.25d0*(1d0-uvw(1,i))*(1d0+uvw(2,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        ai(5,i) =                                           uvw(3,i)
      enddo
    else
      allocate( ai(size(uvw,2),1:5) )
      do i=1,nn
        ai(i,1) = 0.25d0*(1d0-uvw(1,i))*(1d0-uvw(2,i))*(1d0-uvw(3,i))
        ai(i,2) = 0.25d0*(1d0+uvw(1,i))*(1d0-uvw(2,i))*(1d0-uvw(3,i))
        ai(i,4) = 0.25d0*(1d0+uvw(1,i))*(1d0+uvw(2,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        ai(i,3) = 0.25d0*(1d0-uvw(1,i))*(1d0+uvw(2,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        ai(i,5) =                                           uvw(3,i)
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidH6BaseP1
  
  subroutine pyramidH6GradBaseP1(uvw,duai,dvai,dwai,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Hexa [-1,1] x [-1,1] [0,1] dégénéré
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: uvw(:,:)
    real(8), intent(out), pointer :: duai(:,:),dvai(:,:),dwai(:,:)
    logical, intent(in)           :: transpose
    !>
    integer                       :: i,nn
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> b1=1/4 (1-u) (1-v) (1-w)
    !> b2=1/4 (1+u) (1-v) (1-w)
    !> b3=1/4 (1+u) (1+v) (1-w)
    !> b4=1/4 (1-u) (1+v) (1-w)
    !> b5=                   w
    !>
    !> dub1=-1/4 (1-v) (1-w) ; dvb1=-1/4 (1-u) (1-w) ; dwb1=-1/4 (1-u) (1-v)
    !> dub2= 1/4 (1-v) (1-w) ; dvb2=-1/4 (1+u) (1-w) ; dwb2=-1/4 (1+u) (1-v)
    !> dub3= 1/4 (1+v) (1-w) ; dvb3= 1/4 (1-u) (1-w) ; dwb3=-1/4 (1+u) (1+v)
    !> dub4=-1/4 (1+v) (1-w) ; dvb4= 1/4 (1-u) (1-w) ; dwb4=-1/4 (1-u) (1+v)
    !> dub5= 0               ; dvb5= 0               ; dwb5= 1
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => ai(1:np,1:nPt)
    !> Transpose = False => ai(1:nPt,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nn=size(uvw,2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> vertex functions
    if( transpose )then
      allocate( duai(1:5,size(uvw,2)) )
      allocate( dvai(1:5,size(uvw,2)) )
      allocate( dwai(1:5,size(uvw,2)) )
      do i=1,nn
        
        duai(1,i) =-0.25d0*(1d0-uvw(2,i))*(1d0-uvw(3,i))
        duai(2,i) = 0.25d0*(1d0-uvw(2,i))*(1d0-uvw(3,i))
        duai(4,i) = 0.25d0*(1d0+uvw(2,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        duai(3,i) =-0.25d0*(1d0+uvw(2,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        duai(5,i) = 0.00d0
        
        dvai(1,i) =-0.25d0*(1d0-uvw(1,i))*(1d0-uvw(3,i))
        dvai(2,i) =-0.25d0*(1d0+uvw(1,i))*(1d0-uvw(3,i))
        dvai(4,i) = 0.25d0*(1d0+uvw(1,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        dvai(3,i) = 0.25d0*(1d0-uvw(1,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        dvai(5,i) = 0.00d0
        
        dwai(1,i) =-0.25d0*(1d0-uvw(1,i))*(1d0-uvw(2,i))
        dwai(2,i) =-0.25d0*(1d0+uvw(1,i))*(1d0-uvw(2,i))
        dwai(4,i) =-0.25d0*(1d0+uvw(1,i))*(1d0+uvw(2,i)) !> retournement 3<->4
        dwai(3,i) =-0.25d0*(1d0-uvw(1,i))*(1d0+uvw(2,i)) !> retournement 3<->4
        dwai(5,i) = 1.00d0
      enddo
    else
      allocate( duai(size(uvw,2),1:5) )
      allocate( dvai(size(uvw,2),1:5) )
      allocate( dwai(size(uvw,2),1:5) )
      do i=1,nn
        duai(i,1) =-0.25d0*(1d0-uvw(2,i))*(1d0-uvw(3,i))
        duai(i,2) = 0.25d0*(1d0-uvw(2,i))*(1d0-uvw(3,i))
        duai(i,4) = 0.25d0*(1d0+uvw(2,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        duai(i,3) =-0.25d0*(1d0+uvw(2,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        duai(i,5) = 0.00d0
        
        dvai(i,1) =-0.25d0*(1d0-uvw(1,i))*(1d0-uvw(3,i))
        dvai(i,2) =-0.25d0*(1d0+uvw(1,i))*(1d0-uvw(3,i))
        dvai(i,4) = 0.25d0*(1d0+uvw(1,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        dvai(i,3) = 0.25d0*(1d0-uvw(1,i))*(1d0-uvw(3,i)) !> retournement 3<->4
        dvai(i,5) = 0.00d0
        
        dwai(i,1) =-0.25d0*(1d0-uvw(1,i))*(1d0-uvw(2,i))
        dwai(i,2) =-0.25d0*(1d0+uvw(1,i))*(1d0-uvw(2,i))
        dwai(i,4) =-0.25d0*(1d0+uvw(1,i))*(1d0+uvw(2,i)) !> retournement 3<->4
        dwai(i,3) =-0.25d0*(1d0-uvw(1,i))*(1d0+uvw(2,i)) !> retournement 3<->4
        dwai(i,5) = 1.00d0
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidH6GradBaseP1
  
  subroutine pyramidBasePiOrthogonal(ord,a,b,c,mode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define pyramidBasePiOrthogonal 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Psi(a,b,c) = P_i^{0,0}(a)  P_j^{0,0}(b)  (1-c)^max(i,j)  P_k^{2max(i,j)+2,0}(2c-1)
    !> avec {a=u/(1-w) ; b=v/(1/w) ; c=w}
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => mode(1:np,1:n)
    !> Transpose = False => mode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use ieee_arithmetic
    use baseSimplexTools, only: jacobiP
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8)             , pointer :: a(:),b(:),c(:)
    real(8), intent(out), pointer :: mode(:,:)
    logical, intent(in)           :: transpose
    !>
    real(8)                       :: alpha,beta
    integer                       :: iu,iv,iw,iM,iMod,n,np
    real(8), pointer              :: P_iu(:),P_iv(:),P_iw(:)
    !>
    real(8), pointer              :: weight(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidBasePiOrthogonal==1
    print '(">>> pyramidBasePiOrthogonal transpose=",l)',transpose
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)*(2*ord+3)/6
    n =size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidBasePiOrthogonal==1
    print '(4x,"base orthogonale Psi(a,b,c) = P_i^{0,0}(a) P_j^{0,0}(b) (1-c)^max(i,j) P_k^{2max(i,j)+2,0}(2c-1)")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Base orthogonale
    if( .not.transpose )then
      
      allocate(mode(1:n,1:np))
      
      iMod=0
      do iu=0,ord ; do iv=0,ord
        iM=max(iu,iv)
        do iw=0,ord-iM
          iMod=iMod+1
          
          alpha=0d0                 ; beta=0d0
          call jacobiP(n=iu,alpha=alpha, beta=beta,u=a(1:n)      ,jf=P_iu)   !> J_iu^{0             ,0} (a)
          
          alpha=0d0                 ; beta=0d0
          call jacobiP(n=iv,alpha=alpha, beta=beta,u=b(1:n)      ,jf=P_iv)   !> J_iv^{0             ,0} (b)
          
          alpha=real(2*iM+2,kind=8) ; beta=0d0
          call jacobiP(n=iw,alpha=alpha,beta=beta,u=2d0*c(1:n)-1d0,jf=P_iw)  !> J_iw^{2*max(iu,iv)+2,0} (2c-1)
          
          !> Psi(a,b,c) = P_i^{0,0}(a)   P_j^{0,0}(b)  (1-c)^max(i,j)   P_k^{2max(i,j)+2,0}(2c-1)
          mode(1:n,iMod)= P_iu(1:n) * P_iv(1:n) * (1d0-c(1:n))**iM * P_iw(1:n)
          
        enddo
      enddo ; enddo
      
      deallocate(P_iu,P_iv,P_iw)
      
    else
      
      allocate(mode(1:np,1:n))
      
      iMod=0
      do iu=0,ord ; do iv=0,ord
        iM=max(iu,iv)
        do iw=0,ord-iM
          iMod=iMod+1
          
          alpha=0d0                 ; beta=0d0
          call jacobiP(n=iu,alpha=alpha, beta=beta,u=a(1:n)      ,jf=P_iu)   !> J_iu^{0             ,0} (a)
          
          alpha=0d0                 ; beta=0d0
          call jacobiP(n=iv,alpha=alpha, beta=beta,u=b(1:n)      ,jf=P_iv)   !> J_iv^{0             ,0} (b)
          
          alpha=real(2*iM+2,kind=8) ; beta=0d0
          call jacobiP(n=iw,alpha=alpha,beta=beta,u=2d0*c(1:n)-1d0,jf=P_iw)  !> J_iw^{2*max(iu,iv)+2,0} (2c-1)
          
          !> Psi(a,b,c) = P_i^{0,0}(a)   P_j^{0,0}(b)  (1-c)^max(i,j)   P_k^{2max(i,j)+2,0}(2c-1)
          mode(iMod,1:n)= P_iu(1:n) * P_iv(1:n) * (1d0-c(1:n))**iM * P_iw(1:n)
          
        enddo
      enddo ; enddo
      
      deallocate(P_iu,P_iv,P_iw)
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidBasePiOrthogonal==1
    print '("<<< pyramidBasePiOrthogonal")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef pyramidBasePiOrthogonal
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidBasePiOrthogonal
  
  subroutine pyramidNormalization(ord,weight)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define pyramidNormalization 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use pyramidRule, only: P5_gauss
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)            :: ord
    real(8), intent(out) , pointer :: weight(:)
    !>
    integer                        :: i,np
    integer                        :: iGauss,nGauss
    real(8), pointer               :: uGauss(:),vGauss(:),wGauss(:),pGauss(:)
    real(8), pointer               :: a(:),b(:),c(:)
    real(8), pointer               :: psi(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidNormalization==1
    print '(">>> pyramidNormalization")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call P5_gauss(                    &
    &    order=ord                   ,&
    &    nGauss=nGauss               ,&
    &    uGauss=uGauss               ,&
    &    vGauss=vGauss               ,&
    &    wGauss=wGauss               ,&
    &    pGauss=pGauss                )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call pyramiduvw2abc(              &
    &    u=uGauss                    ,&
    &    v=vGauss                    ,&
    &    w=wGauss                    ,&
    &    a=a                         ,&
    &    b=b                         ,&
    &    c=c                          )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call pyramidBasePiOrthogonal(     &
    &    ord=ord                     ,&
    &    a=a                         ,&
    &    b=b                         ,&
    &    c=c                         ,&
    &    mode=psi                    ,&
    &    transpose=.false.            ) !> psi(abc,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(a,b,c)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=size(psi,2) ; allocate(weight(np)) ; weight(1:np)=0d0
    
    do i=1,np
      do iGauss=1,nGauss
        weight(i)=weight(i)+psi(iGauss,i)*psi(iGauss,i)*pGauss(iGauss)
      enddo
    enddo
    
    do i=1,np
      weight(i)=1d0/sqrt(weight(i))
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uGauss,vGauss,wGauss,pGauss)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidNormalization==1
    print '(4x,"weight(",i6,")=",e22.15)',(i,weight(i),i=1,np)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if pyramidNormalization==1
    print '("<<< pyramidNormalization")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef pyramidNormalization
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidNormalization
  
  subroutine pyramidBasePi(ord,a,b,c,mode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define debug 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Psi' = Psi/ int Psi Psi dV
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => mode(1:np,1:n)
    !> Transpose = False => mode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use ieee_arithmetic
    use baseSimplexTools, only: jacobiP
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8)             , pointer :: a(:),b(:),c(:)
    real(8), intent(out), pointer :: mode(:,:)
    logical, intent(in)           :: transpose
    !>
    integer                       :: i,n,np
    real(8), pointer              :: weight(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if debug==1
    print '(">>> pyramidBasePi transpose=",l)',transpose
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call pyramidBasePiOrthogonal( &
    &    ord=ord                 ,&
    &    a=a                     ,&
    &    b=b                     ,&
    &    c=c                     ,&
    &    mode=mode               ,&
    &    transpose=transpose      )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Vandermonde matrix
    call pyramidNormalization(ord=ord,weight=weight)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)*(2*ord+3)/6
    n =size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if debug==1
    print '(4x,"base orthonormale Psi'' = Psi/sqrt(int Psi^2 dV)")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Base orthonormale
    if( .not.transpose )then !> mode(1:n,1:np)
      do i=1,np
        mode(1:n,i)=mode(1:n,i)*weight(i)
      enddo
    else                     !> mode(1:np,1:n)
      do i=1,np
        mode(i,1:n)=mode(i,1:n)*weight(i)
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(weight)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if debug==1
    print '("<<< pyramidBasePi")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef debug
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidBasePi
  
  
  subroutine pyramidGradBasePi(ord,a,b,c,drMode,dsMode,dtMode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define debug 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Psi' = Psi/ sqrt( int Psi^2 dV )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => mode(1:np,1:n)
    !> Transpose = False => mode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use ieee_arithmetic
    use baseSimplexTools, only: jacobiP
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8)             , pointer :: a(:),b(:),c(:)
    real(8), intent(out), pointer :: drMode(:,:)
    real(8), intent(out), pointer :: dsMode(:,:)
    real(8), intent(out), pointer :: dtMode(:,:)
    logical, intent(in)           :: transpose
    !>
    integer                       :: i,n,np
    real(8), pointer              :: weight(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if debug==1
    print '(">>> pyramidGradBasePi transpose=",l)',transpose
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Base orthogonale
    call pyramidGradBasePiOrthogonal(  &
    &    ord=ord                      ,&
    &    a=a                          ,&
    &    b=b                          ,&
    &    c=c                          ,&
    &    drMode=drMode                ,&
    &    dsMode=dsMode                ,&
    &    dtMode=dtMode                ,&
    &    transpose=transpose           )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Normalization coefficient weight=sqrt( 1/int(Psi^2 dV) )
    call pyramidNormalization(ord=ord,weight=weight)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)*(2*ord+3)/6
    n =size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if debug==1
    print '(4x,"base orthonormale Psi'' = Psi/sqrt(int Psi^2 dV)")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Base orthonormale
    if( .not.transpose )then !> mode(1:n,1:np)
      do i=1,np
        drMode(1:n,i)=drMode(1:n,i)*weight(i)
        dsMode(1:n,i)=dsMode(1:n,i)*weight(i)
        dtMode(1:n,i)=dtMode(1:n,i)*weight(i)
      enddo
    else                     !> mode(1:np,1:n)
      do i=1,np
        drMode(i,1:n)=drMode(i,1:n)*weight(i)
        dsMode(i,1:n)=dsMode(i,1:n)*weight(i)
        dtMode(i,1:n)=dtMode(i,1:n)*weight(i)
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(weight)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if debug==1
    print '("<<< pyramidGradBasePi")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef debug
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidGradBasePi
  
  subroutine pyramidGradBasePiOrthogonal(ord,a,b,c,drMode,dsMode,dtMode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define debug 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Psi(a,b,c) = P_i^{0,0}(a)  P_j^{0,0}(b) (1-c)^max(i,j)  P_k^{2 max(i,j)+2,0}(2c-1)
    !> avec {a=x/(1-z) ; b=y/(1-z) ; c=z}
    
    !> (∂Psi/∂x) = (∂Psi/∂a) (∂a/∂x) + (∂Psi/∂b) (∂b/∂x) + (∂Psi/∂c) (∂c/∂x)
    !> avec : (∂a/∂x)=1/(1-z)=1/(1-c) ; (∂b/∂x)=0 ; (∂c/∂x)=0
    !> soit : (∂Psi/∂x) =  (∂a/∂x) (∂Psi/∂a)
    !>                  =  1/(1-c) (∂Psi/∂a)
    !>                  =  (∂P_i^{0,0}/∂a)(a) P_j^{0,0}(b) P_k^{2 max(i,j)+2,0}(2c-1)  (1-c)^(max(i,j)-1)
    
    !> (∂Psi/∂y) = (∂Psi/∂a) (∂a/∂y) + (∂Psi/∂b) (∂b/∂y) + (∂Psi/∂c) (∂c/∂y)
    !> avec : (∂a/∂y)=0 ; (∂b/∂y)=1/(1-z)=1/(1-c) ; (∂c/∂y)=0
    !> soit : (∂Psi/∂x) =  (∂b/∂y) (∂Psi/∂b)
    !>                  =  1/(1-c) (∂Psi/∂b)
    !>                  =  1/(1-c) P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) (1-c)^max(i,j) P_k^{2 max(i,j)+2,0}(2c-1)
    !>                  =  P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) P_k^{2 max(i,j)+2,0}(2c-1)  (1-c)^(max(i,j)-1)
    
    !> (∂Psi/∂z) = (∂Psi/∂a) (∂a/∂z) + (∂Psi/∂b) (∂b/∂z) + (∂Psi/∂c) (∂c/∂z)
    !> avec : (∂a/∂z)=∂( x/(1-z) )/∂z=x/(1-z)^2= a/(1-c)
    !>        (∂b/∂z)=∂( y/(1-z) )/∂z=y/(1-z)^2= b/(1-c)
    !>        (∂c/∂z)= 1
    
    !> (∂a/∂z) (∂Psi/∂a) =  a/(1-c)  (∂P_i^{0,0}/∂a)(a) P_j^{0,0}(b) (1-c)^max(i,j) P_k^{2 max(i,j)+2,0}(2c-1)
    !>                   =  a        (∂P_i^{0,0}/∂a)(a) P_j^{0,0}(b) P_k^{2 max(i,j)+2,0}(2c-1)  (1-c)^(max(i,j)-1) 
    
    !> (∂b/∂z) (∂Psi/∂b) =  b/(1-c)  P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) (1-c)^max(i,j) P_k^{2 max(i,j)+2,0}(2c-1)
    !>                   =  b        P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) P_k^{2 max(i,j)+2,0}(2c-1)  (1-c)^(max(i,j)-1)
    
    !> (∂c/∂z) (∂Psi/∂c) = (∂Psi/∂c)
    !>                   =   P_i^{0,0}(a) P_j^{0,0}(b)   (∂(1-c)^max(i,j)/∂c)    P_k^{2 max(i,j)+2,0}     (2c-1)
    !>                      +P_i^{0,0}(a) P_j^{0,0}(b)     (1-c)^max(i,j)    2 (∂P_k^{2 max(i,j)+2,0}/∂2c)(2c-1)
    !>                   =  -max(i,j) P_i^{0,0}(a) P_j^{0,0}(b)   P_k^{2 max(i,j)+2,0}       (2c-1) (1-c)^(max(i,j)-1)
    !>                      +2        P_i^{0,0}(a) P_j^{0,0}(b) (∂P_k^{2 max(i,j)+2,0}/∂2c-1)(2c-1) (1-c)^(max(i,j)  )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => mode(1:np,1:n)
    !> Transpose = False => mode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use ieee_arithmetic
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8)             , pointer :: a(:),b(:),c(:)
    real(8), intent(out), pointer :: drMode(:,:)
    real(8), intent(out), pointer :: dsMode(:,:)
    real(8), intent(out), pointer :: dtMode(:,:)
    logical, intent(in)           :: transpose
    !
    real(8)                       :: alpha,beta
    integer                       :: i,n,ad
    integer                       :: iu,iv,iw,iM,iMod
    integer                       :: np
    real(8), pointer              :: Pi(:),dPi(:)
    real(8), pointer              :: Pj(:),dPj(:)
    real(8), pointer              :: Pk(:),dPk(:)
    real(8), pointer              :: tmp(:)
    real(8)                       :: coef
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if debug==1
    print '(">>> pyramidGradBasePiOrthogonal transpose=",l)',transpose
    do i=1,n
      print '(4x,"pyramidGradBasePiOrthogonal a,b,c=",3(e22.15,1x),3x,i3,"/",i3)',a(i),b(i),c(i),i,n
    enddo
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)*(2*ord+3)/6
    n =size(a) ! print '("size(a),size(b),size(c)=",3(i10,1x))',size(a),size(b),size(c)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(tmp(1:n))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.transpose )then
      
      allocate(drMode(1:n,1:np),dsMode(1:n,1:np),dtMode(1:n,1:np))
      
      iMod=0
      do iu=0,ord ; do iv=0,ord
        iM=max(iu,iv)
        do iw=0,ord-iM
          iMod=iMod+1 ! print '(/"pyramidGradBasePiOrthogonal: iMod=",i3)',iMod
          
          alpha=0d0                 ; beta=0d0
          call  jacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n)        ,jf= Pi)
          call dJacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n)        ,jf=dPi)
          
          alpha=0d0                 ; beta=0d0
          call  jacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n)        ,jf= Pj)
          call dJacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n)        ,jf=dPj)
          
          alpha=real(2*iM+2,kind=8) ; beta=0d0
          call  jacobiP(n=iw,alpha=alpha,beta=beta,u=2d0*c(1:n)-1d0,jf= Pk)
          call dJacobiP(n=iw,alpha=alpha,beta=beta,u=2d0*c(1:n)-1d0,jf=dPk)
          
          !> tmp(1:n)=(1-c)^(max(i,j)-1)
          do i=1,n
            if( c(i)==1d0 )then
              if    ( iM-1<0 )then
                tmp(i)=1d25 ; print '("pyramidGradBasePiOrthogonal iM=",i2,"<1")',iM
              elseif( iM-1==0 )then
                tmp(i)=1d0
              elseif( iM-1>0 )then
                tmp(i)=0d0
              endif
            else
              tmp(i)=(1d0-c(i))**(iM-1)
            endif
          enddo
          
          !> (∂Psi/∂x) = (∂P_i^{0,0}/∂a)(a)   P_j^{0,0}    (b)  P_k^{2 max(i,j)+2,0}(2c-1)  (1-c)^(max(i,j)-1)
          drMode(1:n,iMod)=dPi(1:n)* Pj(1:n)*Pk(1:n)*tmp(1:n)
          
          !> (∂Psi/∂y) =   P_i^{0,0}    (a) (∂P_j^{0,0}/∂b)(b)  P_k^{2*max(i,j)+2,0}(2c-1)  (1-c)^(max(i,j)-1)
          dsMode(1:n,iMod)= Pi(1:n)*dPj(1:n)*Pk(1:n)*tmp(1:n)
          
          !> (∂Psi/∂z) = (∂a/∂z) (∂Psi/∂a) + (∂b/∂z) (∂Psi/∂b) + (∂c/∂z) (∂2c-1/∂c) (∂Psi/∂2c-1)
          !>           =  a        (∂P_i^{0,0}/∂a)(a)   P_j^{0,0}    (b)   P_k^{2*max(i,j)+2,0}    (2c-1)  (1-c)**(max(i,j)-1)
          !>             +b          P_i^{0,0}    (a) (∂P_j^{0,0}/∂b)(b)   P_k^{2*max(i,j)+2,0}    (2c-1)  (1-c)**(max(i,j)-1)
          !>             -max(i,j)*  P_i^{0,0}    (a)   P_j^{0,0}    (b)   P_k^{2*max(i,j)+2,0}    (2c-1)  (1-c)**(max(i,j)-1)
          !>             +2       *  P_i^{0,0}    (a)   P_j^{0,0}    (b) (∂P_k^{2*max(i,j)+2,0}/∂c)(2c-1)  (1-c)**(max(i,j)  )
          
          dtMode(1:n,iMod)=( a(1:n)          * dPi(1:n) *  Pj(1:n) *  Pk(1:n)                      &
          &                 +b(1:n)          *  Pi(1:n) * dPj(1:n) *  Pk(1:n)                      &
          &                 -real(iM,kind=8) *  Pi(1:n) *  Pj(1:n) *  Pk(1:n) ) * tmp(1:n)         &
          &                 +2d0             *  Pi(1:n) *  Pj(1:n) * dPk(1:n)   * (1d0-c(1:n))**iM  
          
        enddo
      enddo ; enddo
      
    else
      
      allocate(drMode(1:np,1:n),dsMode(1:np,1:n),dtMode(1:np,1:n))
      
      iMod=0
      do iu=0,ord ; do iv=0,ord
        iM=max(iu,iv)
        do iw=0,ord-iM
          iMod=iMod+1 ! print '(/"pyramidGradBasePiOrthogonal: iMod=",i3)',iMod
          
          alpha=0d0                 ; beta=0d0
          call  jacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n)        ,jf= Pi)
          call dJacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n)        ,jf=dPi)
          
          alpha=0d0                 ; beta=0d0
          call  jacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n)        ,jf= Pj)
          call dJacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n)        ,jf=dPj)
          
          alpha=real(2*iM+2,kind=8) ; beta=0d0
          call  jacobiP(n=iw,alpha=alpha,beta=beta,u=2d0*c(1:n)-1d0,jf= Pk)
          call dJacobiP(n=iw,alpha=alpha,beta=beta,u=2d0*c(1:n)-1d0,jf=dPk)
          
          !> tmp(1:n)=(1-c)^(max(i,j)-1)
          do i=1,n
            if( c(i)==1d0 )then
              if    ( iM-1<0 )then
                tmp(i)=1d25 ; print '("pyramidGradBasePiOrthogonal iM=",i2,"<1")',iM
              elseif( iM-1==0 )then
                tmp(i)=1d0
              elseif( iM-1>0 )then
                tmp(i)=0d0
              endif
            else
              tmp(i)=(1d0-c(i))**(iM-1)
            endif
          enddo
          
          !> (∂Psi/∂x) = (∂P_i^{0,0}/∂a)(a)   P_j^{0,0}    (b)  P_k^{2 max(i,j)+2,0}(2c-1)  (1-c)^(max(i,j)-1)
          drMode(iMod,1:n)=dPi(1:n)* Pj(1:n)*Pk(1:n)*tmp(1:n)
          
          !> (∂Psi/∂y) =   P_i^{0,0}    (a) (∂P_j^{0,0}/∂b)(b)  P_k^{2*max(i,j)+2,0}(2c-1)  (1-c)^(max(i,j)-1)
          dsMode(iMod,1:n)= Pi(1:n)*dPj(1:n)*Pk(1:n)*tmp(1:n)
          
          !> (∂Psi/∂z) = (∂a/∂z) (∂Psi/∂a) + (∂b/∂z) (∂Psi/∂b) + (∂c/∂z) (∂2c-1/∂c) (∂Psi/∂2c-1)
          !>           =  a        (∂P_i^{0,0}/∂a)(a)   P_j^{0,0}    (b)   P_k^{2*max(i,j)+2,0}    (2c-1)  (1-c)**(max(i,j)-1)
          !>             +b          P_i^{0,0}    (a) (∂P_j^{0,0}/∂b)(b)   P_k^{2*max(i,j)+2,0}    (2c-1)  (1-c)**(max(i,j)-1)
          !>             -max(i,j)*  P_i^{0,0}    (a)   P_j^{0,0}    (b)   P_k^{2*max(i,j)+2,0}    (2c-1)  (1-c)**(max(i,j)-1)
          !>             +2       *  P_i^{0,0}    (a)   P_j^{0,0}    (b) (∂P_k^{2*max(i,j)+2,0}/∂c)(2c-1)  (1-c)**(max(i,j)  )
          
          dtMode(iMod,1:n)=( a(1:n)          * dPi(1:n) *  Pj(1:n) *  Pk(1:n)                      &
          &                 +b(1:n)          *  Pi(1:n) * dPj(1:n) *  Pk(1:n)                      &
          &                 -real(iM,kind=8) *  Pi(1:n) *  Pj(1:n) *  Pk(1:n) ) * tmp(1:n)         &
          &                 +2d0             *  Pi(1:n) *  Pj(1:n) * dPk(1:n)   * (1d0-c(1:n))**iM  
          
        enddo
      enddo ; enddo
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(tmp)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if debug==1
    print '("<<< pyramidGradBasePiOrthogonal")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef debug
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidGradBasePiOrthogonal
  
  subroutine pyramideMassMatrix(ord, vand, mass)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define pyramidMass 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use baseSimplexTools
    use pyramidRule, only: P5_gauss
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer              :: ord
    real(8), pointer     :: vand(:,:)
    real(8), pointer     :: mass(:)
    !>
    integer              :: i,j,iGauss,np,nPt
    integer              :: ad,cpt,dCpt
    real(8), pointer     :: uvw(:,:)
    real(8), pointer     :: x(:),y(:),z(:),w(:)
    real(8), pointer     :: a(:),b(:),c(:)
    real(8), pointer     :: li(:,:)
    real(8), parameter   :: eps=1d-12
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#if pyramidMass==1
    write(*,'(/">>> pyramideMassMatrix : Mass=\int ai aj dV   compact format ord=",i2)')ord
#endif
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !> Liste des points de Gauss
    call P5_gauss(   &
    &    order=ord  ,&
    &    nGauss=nPt ,&
    &    uGauss=x   ,&
    &    vGauss=y   ,&
    &    wGauss=z   ,&
    &    pGauss=w    )
    
    allocate(uvw(1:3,nPt))
    do i=1,nPt
      uvw(1:3,i)=[x(i),y(i),z(i)]
    enddo
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !> Evaluation des fonctions de Lagrange aux points uvw
    !> Transpose=.true. => li(np,nPt)
    call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
    call pyramidLagrange3Dv(ord=ord,vand=vand,a=a,b=b,c=c,lx=li,transpose=.true.)  !> li(uvw) = Inverse[Transpose[Vand]].Psi[uvw]
    deallocate(a,b,c)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    np=size(li,1)
    !allocate( mass(np*(np+1)/2)) ; mass(:)=0d0
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    mass(:)=0d0
    do iGauss=1,nPt
      ad=0 ; cpt=0 ; dCpt=0
      do i=1,np
        do j=1,np
          cpt=cpt+1
          if( cpt>dCpt )then
            ad=ad+1 ; mass(ad)=mass(ad)+li(i,iGauss)*li(j,iGauss)*w(iGauss)
          endif
        enddo
        cpt=0 ; dCpt=dCpt+1
      enddo
    enddo
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw)
    deallocate(li)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#if pyramidMass==1
    call display(title="    Mass Matrix Compact Format",vec=mass)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#if pyramidMass==1
    write(*,'("<<< pyramideMassMatrix : Mass=\int ai aj dV   compact format ord=",i2)')ord
#endif
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#undef pyramidMass
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    return
  end subroutine pyramideMassMatrix
  
  subroutine pyramidNodes(ord, uvw, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity pyramid
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:,:)
    logical, intent(in)           :: display
    !---
    integer                       :: iu,iv,iw,iOrd
    integer                       :: iNod,nNod
    real(8)                       :: a
    real(8)                       :: sub(0:ord)
    
    integer                       :: iCel,nPyr,nTet
    integer, allocatable          :: pyram(:,:)
    integer, allocatable          :: tetra(:,:)
    logical                       :: mesh
    logical                       :: tf
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(/"Building Pyramid Equidistant Nodes")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nNod=(ord+1)*(ord+2)*(2*ord+3)/6 !> = \sum_{k=1}^{ord+1} k^2
    if( display )print '(3x,"nDeg=",i6)',nNod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(uvw(3,nNod))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iNod=0
    do iw=0,ord ; print '("iw=",i3)',iw
      
      !>>>>>> square [-a,+a]x[-a,+a] @ level iw
      if( ord==0 )then
        a=0d0
      else
        a=1d0-real(iw,kind=8)/real(ord,kind=8) 
      endif
     !if( display )print '("a=",f12.9)',a
      !<<<<<<
      
      !>>>>>> subdivision of [-a,+a]
      sub(0)=-a
      do iu=1,ord-iw-1
        sub(iu)=-a+(2d0*a)*real(iu,kind=8)/real(ord-iw,kind=8)
      enddo
      sub(ord-iw)=+a
     !if( display )print '("sub=",15(f12.9,1x))',sub(0:ord-iw)
      !<<<<<<
      
      !>>>>>> 2D Grid (sub(iu),sub(iv))
      iOrd=ord-iw
      do iv=0,iOrd
        do iu=0,iOrd
          iNod=iNod+1
          uvw(1:3,iNod)=[ sub(iu), sub(iv), real(iw,kind=8)/real(ord,kind=8)  ]
        enddo
      enddo  
      !<<<<<<
      
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      print '(3x,"Vertices")'
      do iNod=1,nNod
        print '(6x,"uvw(",i6,")=",3(f12.9,1x))',iNod,uvw(1:3,iNod)
      enddo
      print '(3x,"end Vertices")'
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '("end Building Pyramid Equidistant Nodes")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidNodes
  
  subroutine pyramidNodesOpt(ord, uvw, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity pyramid
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)             :: ord
    real(8), intent(inout), pointer :: uvw(:,:) !> Pyramide optimized points
    logical, intent(in)             :: display
    !>
    real(8), pointer                :: uv (:,:) !> Triangle optimized points
    real(8)                         :: xyz(3)
    integer                         :: iu,iv,iw,ad
    integer                         :: iNod,iNod0,iNod1
    integer                         :: jNod,jNod0,jNod1
    integer                         :: nNod,nNodi,nNode
    real(8)                         :: sd2(3,0:ord)
    real(8)                         :: alphaU,betaU
    real(8)                         :: alphaV,betaV
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(/"Building Pyramid Optimized Nodes")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nNod=(ord+1)*(ord+2)*(2*ord+3)/6 !> = \sum_{k=1}^{ord+1} k^2
    if( display )print '(/3x,"nDeg=",i6)',nNod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(uvw(1:3,1:nNod))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( ord==0 )then
      
      uvw(1:3,1)=[0d0, 0d0, 0.2d0] !> barycentre de la pyramide
      
      nNode=1
      nNodi=0
      
    elseif( ord==1 )then
      
      uvw(1:3,1)=[-1d0,-1d0, 0d0]
      uvw(1:3,2)=[ 1d0,-1d0, 0d0]
      uvw(1:3,4)=[ 1d0, 1d0, 0d0] !> Retournement 3<->4
      uvw(1:3,3)=[-1d0, 1d0, 0d0] !> Retournement 3<->4
      uvw(1:3,5)=[ 0d0, 0d0, 1d0]
      
      nNode=5
      nNodi=0
      
    else
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call pyramidSide2NodesOpt(ord=ord, uv=uv, display=.false.) !> face triangle (necessaire)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if( display )then
        print '(3x,"3D Triangle Optimized Nodes")'
        iNod=0
        do iw=0,ord
          print '(6x,"level:",i3)',iw
          do iv=0,ord-iw
            iNod=iNod+1
            print '(9x,"uv(",i6,")=",3(f12.9,1x))',iNod,uv(1:3,iNod)
          enddo
        enddo
        print '(3x,"end")'
      endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> Traitement des faces
      if( display )print '(3x,"Nodes optimization over   pyramid")'
      !> iw=0      -> Side1
      !> iv=0      -> Side2
      !> iu=ord-iw -> Side3
      !> iv=ord-iw -> Side4
      !> iu=0      -> Side5
      nNode=0
      iNod=0 ; jNod=0
      do iw=0,ord
        
        sd2(1:3,0:ord-iw)=uv(1:3,jNod +1:jNod+ord-iw +1)
       !print '(/"sd2=",12(f12.5,2x))',sd2(1,0:ord-iw)
       !print '( 4x    ,12(f12.5,2x))',sd2(2,0:ord-iw)
       !print '( 4x    ,12(f12.5,2x))',sd2(3,0:ord-iw)
        
        do iv=0,ord-iw
          jNod=jNod+1
          do iu=0,ord-iw
            iNod=iNod+1
            
            if    ( iw==0 )then !> side1
              nNode=nNode+1
              uvw(1:3,iNod)=[sd2(1,iu),sd2(1,iv),0d0]
            elseif( iv==0 )then !> side2
              nNode=nNode+1
              ad=iu
              uvw(1:3,iNod)=sd2(1:3,ad)
            elseif( iu==ord-iw )then !> side3
              !> side2 -> side3 rotation +pi/2 (zz')
              !> [ 0 -1  0] [x]   [-y]
              !> [+1  0  0] [y] = [ x]
              !> [ 0  0 +1] [z]   [ z]
              nNode=nNode+1
              ad=iv
              uvw(1:3,iNod)=[-sd2(2,ad),+sd2(1,ad),sd2(3,ad)]
            elseif( iv==ord-iw )then !> side4
              !> side2 -> side4 rotation +pi (zz')
              !> [-1  0  0] [x]   [-x]
              !> [ 0 -1  0] [y] = [-y]
              !> [ 0  0 +1] [z]   [ z]
              nNode=nNode+1
              ad=ord-iw-iu
              uvw(1:3,iNod)=[-sd2(1,ad),-sd2(2,ad),sd2(3,ad)]
            elseif( iu==0 )then !> side5
              !> side2 -> side5 rotation 3pi/2 (zz')
              !> [ 0 +1  0] [x]   [+y]
              !> [-1  0  0] [y] = [-x]
              !> [ 0  0 +1] [z]   [ z]
              nNode=nNode+1
              ad=ord-iw-iv
              uvw(1:3,iNod)=[sd2(2,ad),-sd2(1,ad),sd2(3,ad)]
            endif
            
          enddo
        enddo
      enddo
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      deallocate(uv)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> Traitement hors faces
      if( display )print '(3x,"Nodes optimization inside pyramid")'
      
      nNodi=(ord-1)*(ord-2)*(2*ord-3)/6 !> internal nodes
      nodesInside: if( .not.nNodi==0 )then
        
        iNod =0
        do iw=0,ord
          if( display )print '(6x,"level:",i3)',iw
          do iv=0,ord-iw
            do iu=0,ord-iw
              iNod=iNod+1
              
              if(.not.(iw==0.or.iv==0.or.iv==ord-iw.or.iu==0 .or.iu==ord-iw) )then
                
                iNod0=pyramidIdx(ord=ord, iu=0     ,iv=iv    ,iw=iw)
                iNod1=pyramidIdx(ord=ord, iu=ord-iw,iv=iv    ,iw=iw)
                jNod0=pyramidIdx(ord=ord, iu=iu    ,iv=0     ,iw=iw)
                jNod1=pyramidIdx(ord=ord, iu=iu    ,iv=ord-iw,iw=iw)
                
                !if( display )then
                !  print '(/9x,"iNod=",i4,3x,"(iu,iv,iw)=",3(i4,1x))',iNod,iu,iv,iw
                !  print '( 9x,"iNod0-iNod1",i4,1x,i4)',iNod0,iNod1
                !  print '( 9x,"jNod0-jNod1",i4,1x,i4)',jNod0,jNod1
                !endif
                
                alphaU=real(iu,kind=8)/real(ord-iw,kind=8) ; betaU=1d0-alphaU
                alphaV=real(iv,kind=8)/real(ord-iw,kind=8) ; betaV=1d0-alphaV
                
                xyz(1:3)=( alphaU*uvw(1:3,iNod1)+betaU*uvw(1:3,iNod0) &
                &         +alphaV*uvw(1:3,jNod1)+betaV*uvw(1:3,jNod0) )/2d0
                
                uvw(1:3,iNod)=xyz(1:3)
                
              endif
              
            enddo
            iNod0=iNod1
          enddo
        enddo
        
      endif nodesInside
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      print '(3x,"Vertices")'
      iNod=0
      do iw=0,ord
        print '(6x,"level: ",i3)',iw
        do iv=0,ord-iw ; do iu=0,ord-iw
          iNod=iNod+1
          print '(6x,"uvw(",i6,")=",3(f12.9,1x))',iNod,uvw(1:3,iNod)
        enddo ; enddo
      enddo
      print '(3x,"end")'
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(3x,"nn=",i4,3x,"ne=",i4,3x,"ni=",i4)',nNod,nNode,nNodi
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '("end Building Pyramid Optimized Nodes")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidNodesOpt
  
  subroutine pyramidSide2NodesOpt(ord,uv,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Coordonnées des noeuds optimises
    !> sur la face2 de la pyramide (iv=0)
    !> nd 1 2 5
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use baseSimplex3D, only: nodes3D,nodes3DOpt,trianglesConnectivity,permutation
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uv (:,:)
    logical, intent(in)           :: display
    !>
    real(8), pointer     :: uvw(:,:)
    integer              :: ad,iSide
    integer              :: iNod,nNod,iu,iv,iw
    integer, allocatable :: conec(:,:)
    integer, allocatable :: idx(:)
    real(8)              :: rot(3,3),xyz(1:3),cos_a,sin_a
    real(8)              :: alpha
    character(3)         :: sfx
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(/"Pyramid Side2 (iv=0) Optimized Nodes")' !> confere routine tetraTest pour detail
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Points optimises sur une face triangle
    call nodes3D   (ord=ord,uvw=uvw,display=.false.) !> nodes3D sur tetra
    call nodes3Dopt(ord=ord,uvw=uvw,display=.false.) !> nodes3D sur tetra (optimisation) 
    call trianglesConnectivity(ord=ord,conec=conec)  !> Extraction des faces du tetra
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> On fait tourner la face pour mettre l'indice 1 en bas à gauche
    nNod=size(conec,1) ; allocate(idx(nNod))
    call permutation(order=ord, move=-1, flip=.false., dg=idx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Coordonnées de la face triangulaire 3
    allocate(uv(3,nNod))
    do iNod=1,nNod
      ad=conec(idx(iNod),3) !> idx sur face3 du triangle
      uv(1,iNod)=uvw(1,ad)
      uv(2,iNod)=uvw(3,ad)
      uv(3,iNod)=1d0-uv(1,iNod)-uv(2,iNod)
    enddo
    deallocate(uvw,conec,idx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Des transformations sont realisees pour plaquer sur Side2 (iu=0)
    
    !> Rotation autour de (zz') angle -3pi/4
    cos_a=-sqrt(2d0)/2d0 ; sin_a=-sqrt(2d0)/2d0
    rot(1,1:3) = [cos_a,-sin_a, 0d0]
    rot(2,1:3) = [sin_a, cos_a, 0d0]
    rot(3,1:3) = [0d0  , 0d0  , 1d0]
    do iNod=1,nNod
      xyz(1:3)=uv(1:3,iNod)
      uv(1:3,iNod)= matmul(rot(1:3,1:3),xyz(1:3)) !> rotation -3pi/4 autour de (zz')
    enddo
    
    !> Homothétie
    alpha=2d0/abs( uv(1,ord+1)-uv(1,1) )
    uv(1:2,1:nNod)=alpha*uv(1:2,1:nNod) !> on touche pas à uv(3,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Ecriture du triangle TriangleOptPi.mesh (Side2 de la pyramide)
    if( display )then
      
      if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
      if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
      if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
      
      write(*,'(3x,"writing: ",a)')"TriangleOptP"//sfx//".mesh"
      open(unit=10,file="TriangleOptP"//sfx//".mesh",action='write')
      write(10,'( "MeshVersionFormatted 1")' )
      write(10,'(/"Dimension"/,"3")' )
      write(10,'(/"Vertices"/,i3)' )nNod
      iNod=0
      do iw=0,ord
        do iv=0,ord-iw
          iNod=iNod+1
          write(10,'(3(e22.15,1x),i3)'),uv(1,iNod),uv(2,iNod),uv(3,iNod),0
        enddo
      enddo
      write(10,'(/"Triangles"/,i6)' )ord*ord
      iNod=0
      do iw=0,ord-1
        do iv=0,ord-iw-1
          iNod=iNod+1
          write(10,'(3(i3,1x),3x,i3)' )iNod,iNod+1,iNod+ord-iw+1, 0
          if( .not.iv==ord-iw-1 )then
            write(10,'(3(i3,1x),3x,i3)' )iNod+1,iNod+ord-iw+2,iNod+ord-iw+1, 0
          endif
        enddo
        iNod=iNod+1
      enddo  
      write(10,'(/"End")')
      close(10)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidSide2NodesOpt
    
  subroutine pyramidMesh3D(ord,uvw,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define fortran 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use M_libmesh6_api
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:,:)
    logical, intent(in)           :: display
    !>
    integer                       :: iu,iv,iw
    integer                       :: iNod ,nNod, iCel
    logical                       :: test
    real(4)                       :: dist
    character(3)                  :: sfx
    !> libmesh
    character(256)                :: name
    
    integer                       :: ins,ver,res,geo
    integer, allocatable          :: TypTab(:)
    real(4)                       :: xyz(3)
    integer                       :: nFld,kind(1)
    integer                       :: iNod0,nNod0
    real(8), allocatable          :: uvw0(:,:)
    integer, allocatable          :: mark(:)
    integer                       :: nTetr
    integer, allocatable          :: tetr(:,:)
    integer                       :: nTria
    integer, allocatable          :: tria(:,:)
    !> Reordering
    real(8), parameter            :: eps=1d-10
    integer, allocatable          :: indx(:)
    integer                       :: iNod0Min,iNodMin
    real(8)                       :: d2,d2Min
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(/"Building Pyramid Volumic Mesh")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
    if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
    if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#if fortran==1
    open(unit=150,file='nodesPyramid_P'//sfx//'.f90',action='write',status='unknown')
    !write(150,'("    select case(Pi)")')
#endif
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Ecriture
    nNod=0
    do iw=0,ord
      do iv=0,ord-iw
        do iu=0,ord-iw
          test= (iu==0.or.iu==ord-iw).or.( iv==0.or.iv==ord-iw).or.iw==0
          if( .not.test )nNod=nNod+1
        enddo
      enddo
    enddo
    
    if( .not. nNod==0 )then
      !> mesh
      name="nodes3DP"//sfx//".mesh" ; if( display )print '(3x,"Building and Writing: ",a)',trim(name)
      open(unit=10,file=trim(name),action='write')
      write(10,'( "MeshVersionFormatted 1")' )
      write(10,'(/"Dimension"/,"3")' )
      write(10,'(/"Vertices"/,i10)' )nNod
      iNod=0
      do iw=0,ord
        do iv=0,ord-iw
          do iu=0,ord-iw
            iNod=iNod+1
            test= (iu==0.or.iu==ord-iw).or.( iv==0.or.iv==ord-iw).or.iw==0
            if( .not.test )write(10,'(3(e22.15,1x),i3)')uvw(1:3,iNod),0
          enddo
        enddo
      enddo
      write(10,'(/"End")')
      close(10)
      
      !> sol
      dist=1e0/real(ord,kind=4)
      name="nodes3DP"//sfx//".sol" ;  if( display )print '(3x,"Building and Writing: ",a)',trim(name)
      open(unit=10,file=trim(name),action='write')
      write(10,'( "MeshVersionFormatted 1")' )
      write(10,'(/"Dimension"/,"3")' )
      write(10,'(/"SolAtVertices"/,i10)' )nNod
      write(10,'("1",1x,"1")' ) !> 1 field 1=scalar
      write(10,'(e22.15)')((dist),iNod=1,nNod)
      write(10,'(/"End")')
      close(10)
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(3x,"Building Pyramid Volumic Mesh with ghs3d: ",a)',"PyramidP"//sfx//".mesh"
    if( nNod==0 )then
      call system("ghs3d -O 1  -exit 3 -in PyramidSkinP"//sfx//".mesh                          -out PyramidP"//sfx//".mesh > ghs3d.log")
    else
      call system("ghs3d -O 1  -exit 3 -in PyramidSkinP"//sfx//".mesh -force nodes3DP"//sfx//" -out PyramidP"//sfx//".mesh > ghs3d.log")
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(3x,"rm -f nodes3DP",a,".mesh")',sfx
    if( display )print '(3x,"rm -f nodes3DP",a,".sol")',sfx
    if( display )print '(3x,"rm -f PyramidSkinP",a,".mesh")',sfx
    call system("rm -f nodes3DP"//sfx//".mesh")
    call system("rm -f nodes3DP"//sfx//".sol")
    call system("rm -f PyramidSkinP"//sfx//".mesh")
    call system("rm -f PyramidP"//sfx//".sol")
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    name="PyramidP"//sfx//".mesh"
    if( display )print '(3x,"Reading Pyramid Volumic Mesh: ",a)',trim(name)
    ins=GmfOpenMeshF77(trim(name),GmfRead,ver,geo) ! print '(3x,"ins=",i3)',ins
    nNod0=GmfStatKwdF77(ins,GmfVertices,ver,0,TypTab) ; if( display )print '(6x,"nNod0=",i10)',nNod0
    allocate(uvw0(geo+1,nNod0),mark(nNod0))  ! il faut mettre geo+1 pour ensuite reoordonner les noeuds
    res=GmfGotoKwdF77(ins,GmfVertices)
    select case(ver)
    case(1) !> real(4)
      do iNod0=1,nNod0
        call GmfGetVertex3dr4(ins,xyz(1),xyz(2),xyz(3),mark(iNod0))
        uvw0(1:3,iNod0)=xyz(1:3)
        uvw0(4,iNod0)=1d0-uvw0(1,iNod0)-uvw0(2,iNod0)-uvw0(3,iNod0)
      enddo
    case(2) !> real(8)
      do iNod0=1,nNod0
        call GmfGetVertex3dr8(ins,uvw0(1,iNod0),uvw0(2,iNod0),uvw0(3,iNod0),mark(iNod0))
        uvw0(4,iNod0)=1d0-uvw0(1,iNod0)-uvw0(2,iNod0)-uvw0(3,iNod0)
      enddo
    end select
    
    nTetr=GmfStatKwdF77(ins,GmfTetrahedra,0,0,TypTab) ; if( display )print '(6x,"nTetr=",i10)',nTetr
    allocate(tetr(5,nTetr))
    res=GmfGotoKwdF77(ins,GmfTetrahedra)
    do iCel=1,nTetr
      call GmfGetTetrahedron(ins,tetr(1,iCel),tetr(2,iCel),tetr(3,iCel),tetr(4,iCel),tetr(5,iCel))
    enddo
    
    nTria=GmfStatKwdF77(ins,GmfTriangles,0,0,TypTab) ; if( display )print '(6x,"nTria=",i10)',nTria
    allocate(tria(4,nTria))
    res=GmfGotoKwdF77(ins,GmfTriangles)
    do iCel=1,nTria
      call GmfGetTriangle(ins,tria(1,iCel),tria(2,iCel),tria(3,iCel),tria(4,iCel))
    enddo
    
    res=GmfCloseMeshF77(ins)
    if( display )print '(3x,"end Reading")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(3x,"Reordering Nodes")'
    
    nNod =size(uvw ,2)
    nNod0=size(uvw0,2)
    
    allocate(indx(nNod0)) ; indx(1:nNod0)=0
    do iNod=1,nNod
      iNod0=1 ; d2Min=1d50 ; iNodMin=1
      loop1 : do
        
        d2= (uvw0(1,iNod0)-uvw(1,iNod))**2 &
        &  +(uvw0(2,iNod0)-uvw(2,iNod))**2 &
        &  +(uvw0(3,iNod0)-uvw(3,iNod))**2
        
        if( d2<d2Min )then
          iNod0Min=iNod0 ; d2Min=d2
        endif
        
        if( d2<eps )then
          indx(iNod0)=iNod
          exit loop1
        endif
        
        iNod0=iNod0+1
        if( iNod0>nNod0 )then
          print '("uvw (",i6,")=",3(f12.5,1x),1x,"d2Min=",e12.5)',iNod    ,uvw (1:3,iNod    ),d2Min
          print '("uvw0(",i6,")=",3(f12.5,1x),1x,"d2Min=",e12.5)',iNod0Min,uvw0(1:3,iNod0Min),d2Min
          stop"problem @ writeMeshSkin3D"
        endif
      enddo loop1
    enddo
    
    if( display )print '(3x,"end Reordering")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Test pour verifier que tous les noeuds sont en correspondance
    if( .not.count( indx(:)==0 )==0 )then
      print '("count(indx==0)=",i10)',count(indx(:)==0)
      stop "problem @ pyramidMesh3D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(3x,"rm -f PyramidP",a,".mesh")',sfx
    call system("rm -f PyramidP"//sfx//".mesh")
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Reecriture de TetraPi.mesh
    name="PyramidP"//sfx//".mesh" ; if( display )print '(3x,"Writing Reordered Pyramid Volumic Mesh: ",a)',trim(name)
    ins=GmfOpenMeshF77(trim(name),GmfWrite,1,3) ; if( display )print '(6x,"nNod =",i10)',nNod
    
    res=GmfSetKwdF77(ins,GmfVertices,nNod,0,TypTab)
    do iNod=1,nNod
      call GmfSetVertex3dr4(ins                     ,&
      &                     real(uvw(1,iNod),kind=4),&
      &                     real(uvw(2,iNod),kind=4),&
      &                     real(uvw(3,iNod),kind=4),&
      &                     i3=0                     )
    enddo
    
    res=GmfSetKwdF77(ins,GmfTetrahedra,nTetr,0,TypTab) ; if( display )print '(6x,"nTetr=",i10)',nTetr
    do iCel=1,nTetr
      call GmfSetTetrahedron(ins               ,&
      &                      indx(tetr(1,iCel)),&
      &                      indx(tetr(2,iCel)),&
      &                      indx(tetr(3,iCel)),&
      &                      indx(tetr(4,iCel)),&
      &                           tetr(5,iCel)  )
    enddo
    
    res=GmfSetKwdF77(ins,GmfTriangles,nTria,0,TypTab) ;  if( display )print '(6x,"nTria=",i10)',nTria
    do iCel=1,nTria
      call GmfSetTriangle(ins               ,&
      &                  indx(tria(1,iCel)),&
      &                  indx(tria(2,iCel)),&
      &                  indx(tria(3,iCel)),&
      &                       tria(4,iCel)  )
    enddo
    
    res=GmfCloseMeshF77(ins)
    if( display )print '(3x,"end")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if fortran==1
    print '(/"Write Fortran Source Vertices")'
    
    write(150,'("    case(",i3,")")')ord
    !write(150,'("      ")')
    do iCel=1,nTetr
      write(150,'("      nd(1:4,",i4,")=[",i4,",",i4,",",i4,",",i4,"]")')&
      &                      iCel              ,&
      &                      indx(tetr(1,iCel)),&
      &                      indx(tetr(2,iCel)),&
      &                      indx(tetr(3,iCel)),&
      &                      indx(tetr(4,iCel))
    enddo
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#if fortran==1
    !write(150,'("    end select")')
    close(150)
#endif
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '("end Building Pyramid Volumic Mesh")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidMesh3D
  
  subroutine pyramidSkin3D(ord, uvw, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Cette routine construit le maillage surfacique
    !> d'une pyramide Pord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:,:)
    logical, intent(in)           :: display
    !---
    integer                       :: iu,iv,iw
    integer                       :: iNod,nNod
    integer                       :: iCel,nCel
    integer                       :: nQuadr,nTrian
    integer, allocatable          :: quadr(:,:)
    integer, allocatable          :: trian(:,:)
    character(3)                  :: sfx
    character(256)                :: name
    logical, allocatable          :: keep(:)
    integer, allocatable          :: indx(:)
    integer                       :: iNod0,nNod0
    real(8), pointer              :: uvw0(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(/"Building Pyramid Surfasic Mesh")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
    if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
    if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nNod=size(uvw,2)
    allocate(quadr(5,  ord*ord))
    allocate(trian(4,4*ord*ord))
    nQuadr=0
    nTrian=0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> iw=0b -> side1
    if( display )print '(3x,"Side1")'
    do iw=0,0
      do iv=0,ord-iw-1
        do iu=0,ord-iw-1
          nQuadr=nQuadr+1
          quadr(1:5,nQuadr)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
          &                  1                                           ]
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> iv=0 -> Side2
    if( display )print '(3x,"Side2")'
    do iw=0,ord-1
      do iv=0,0
        do iu=0,ord-iw-1
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  2                                           ]
          if( .not.iu==ord-iw-1 )then
            nTrian=nTrian+1
            trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
            &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw+1),&
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv , iw=iw+1),&
            &                  2                                           ]
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> iu=ord-iw -> Side3
    if( display )print '(3x,"Side3")'
    do iw=0,ord-1
      do iv=0,ord-iw-1
        do iu=ord-iw-1,ord-1-iw
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  3                                           ]
          if( .not.iv==ord-iw-1 )then
            nTrian=nTrian+1
            trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw+1),&
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
            &                  3                                           ]
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> iv=ord-iw -> Side4
    if( display )print '(3x,"Side4")'
    do iw=0,ord-1
      do iv=ord-iw-1,ord-iw-1
        do iu=0,ord-iw-1
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  4                                           ]
          if( .not.iu==ord-iw-1 )then
            nTrian=nTrian+1
            trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
            &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw+1),&
            &                  4                                           ]
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> iu=0      -> side5
    if( display )print '(3x,"Side5")'
    do iw=0,ord-1
      do iv=0,ord-iw-1
        do iu=0,0
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
          &                  5                                           ]
         !print '("Trian(",i3,")=",4(I3,1x))',nTrian,trian(1:4,nTrian)
          if( .not.iv==ord-iw-1 )then
            nTrian=nTrian+1
            trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),& 
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw+1),&
            &                  5                                           ]
           !print '("Trian(",i3,")=",4(I3,1x))',nTrian,trian(1:4,nTrian)
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Removing Extra vertices
    if( display ) print '(3x,"Removing Extra vertices")'
    
    allocate(keep(nNod)) ; keep(:)=.false.
    do iCel=1,nQuadr
      do iNod=1,4
        keep( quadr(iNod,iCel) )=.true.
      enddo
    enddo
    do iCel=1,nTrian
      do iNod=1,3
        keep( trian(iNod,iCel) )=.true.
      enddo
    enddo
    
    allocate(indx(nNod) )
    iNod0=0
    do iNod=1,nNod
      if( keep(iNod) )then
        iNod0=iNod0+1
        indx(iNod)=iNod0
      endif
    enddo
    nNod0=iNod0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Ecriture
    name="PyramidSkinP"//sfx//".mesh"       
    open(unit=10,file=trim(name),action='write')
    write(10,'( "MeshVersionFormatted 1")' )
    write(10,'(/"Dimension"/,"3")')
    write(10,'(/"Vertices"/,i10)' )nNod0
    do iNod=1,nNod
      if( keep(iNod) )write(10,'(3(e22.15,1x),i3)')uvw(1:3,iNod),0
    enddo
    write(10,'(/"Quadrilaterals"/,i10)' )nQuadr
    do iCel=1,nQuadr
      write(10,'(4(i6,1x),i3)')indx(quadr(1:4,iCel)),quadr(5,iCel)
    enddo
    write(10,'(/"Triangles"/,i10)' )nTrian
    do iCel=1,nTrian
      write(10,'(3(i6,1x),i3)')indx(trian(1:3,iCel)),trian(4,iCel)
    enddo
    write(10,'(/"End")')
    close(10)
    
    if( display )then
      print '(3x,"Writing ",a)',trim(name)
      print '(6x,"Vertices:       ",i10)',nNod0
      print '(6x,"Quadrilaterals: ",i10)',nQuadr
      print '(6x,"Triangles:      ",i10)',nTrian
      print '(3x,"end Writing")'
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display ) print '(3x,"Cleanning Memory")'
    deallocate(quadr,trian)
    deallocate(keep,indx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '("end Building Pyramid Surfasic Mesh")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidSkin3D
  
  function pyramidIdx(ord, iu,iv,iw) result(idx)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer :: ord
    integer :: iu,iv,iw
    integer :: idx
    !>
    integer :: k
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    idx=0
    do k=0,iw-1
      idx=idx+(ord+1-k)**2
    enddo
    
   !print '("iu,iv,iw=",3(i3,2x),"idx=",i3)',iu,iv,iw,idx
    idx=iu+iv*(ord+1-iw)+idx +1
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end function pyramidIdx
  
  subroutine pyramidSides3D(ord, sidesIdx, sides, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! input: ord=polynomial order of interpolant
    ! output: sidesIdx,sides (sides degrees)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>        y
    !>        ^
    !>    4x--|---x3
    !>     |\ S4 /|
    !>     | \| / |
    !>     |S5x-S3--->x
    !>     | /  \ |
    !>     |/ S2 \|
    !>    1x------x2
    !>
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use baseSimplex3D, only:  permutation
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer             , intent(in)  :: ord
    logical             , intent(in)  :: display
    integer             , intent(out) :: sidesIdx(6)
    integer, allocatable, intent(out) :: sides(:)
    !---
    integer                           :: iu,iv,iw
    integer                           :: iNod,iSide,iDg
    integer                           :: ijkToDeg(0:ord,0:ord,0:ord)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    stop "stop @ pyramidSides3D probleme avec les face"
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(/"Extract Pyramid List of DOF per Side")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(sides((ord+1)*(ord+1)+4*(ord+1)*(ord+2)/2))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sidesIdx(1)=0
    sidesIdx(2)=sidesIdx(1)+(ord+1)*(ord+1)   !> square
    sidesIdx(3)=sidesIdx(2)+(ord+1)*(ord+2)/2 !> triangle
    sidesIdx(4)=sidesIdx(3)+(ord+1)*(ord+2)/2 !> triangle
    sidesIdx(5)=sidesIdx(4)+(ord+1)*(ord+2)/2 !> triangle
    sidesIdx(6)=sidesIdx(5)+(ord+1)*(ord+2)/2 !> triangle
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iDg=0
    do iw=0,ord ; do iv=0,ord-iw ; do iu=0,ord-iw
      iDg=iDg+1 ; ijkToDeg(iu,iv,iw)=iDg
    enddo ; enddo ; enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Side1 square iw=0  (retournée)
    iw=0 ; do iu=0,ord ; do iv=0,ord
      sidesIdx(1)=sidesIdx(1)+1
      sides( sidesIdx(1) )=ijkToDeg(iu,iv,iw)
    enddo ; enddo
    
    !> Side2 triangle iv=0
    do iw=0,ord ; iv=0 ; do iu=0,ord-iw
      sidesIdx(2)=sidesIdx(2)+1
      sides( sidesIdx(2) )=ijkToDeg(iu,iv,iw)
    enddo ; enddo
    
    !> Side3 triangle iu=ord-iw
    do iw=0,ord ; do iv=0,ord-iw ; iu=ord-iw
      sidesIdx(3)=sidesIdx(3)+1
      sides( sidesIdx(3) )=ijkToDeg(iu,iv,iw)
    enddo ; enddo
    
    !> Side4 triangle iv=ord-iw (retournée)
    do iw=0,ord ; iv=ord-iw ; do iu=ord-iw,0,-1
      sidesIdx(4)=sidesIdx(4)+1
      sides( sidesIdx(4) )=ijkToDeg(iu,iv,iw)
    enddo ; enddo
    
    !> Side5 triangle iu=0  (retournée)
    do iw=0,ord ; do iv=ord-iw,0,-1 ; iu=0
      sidesIdx(5)=sidesIdx(5)+1
      sides( sidesIdx(5) )=ijkToDeg(iu,iv,iw)
    enddo ; enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      do iSide=1,5
        print '(3x,"Side",i1,": ",$)',iSide
        do iNod=sidesIdx(iSide)+1,sidesIdx(iSide+1)
          print '(i5,1x,$)',sides(iNod)
        enddo
        print '()'
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '("end Extract Pyramid List of DOF per Side")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidSides3D
  
  subroutine pyramidReadXYZout3D(xyzOut,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use M_libmesh6_api
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(out), pointer :: xyzOut (:,:)
    logical, intent(in)           :: display
    !-
    integer                       :: i,j,nVert
    character(256)                :: name
    integer                       :: iNod0,nNod0
    integer                       :: mark
    integer                       :: ins,ver,res,geo
    integer, allocatable          :: TypTab(:)
    real(4)                       :: xyz(3)
    integer                       :: nFld,kind(1)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    name="Pyramids.meshb"
    if( display )print '(3x,"Reading Pyramid Volumic Mesh: ",a)',trim(name)
    ins=GmfOpenMeshF77(trim(name),GmfRead,ver,geo) ! print '(3x,"ins=",i3)',ins
    nNod0=GmfStatKwdF77(ins,GmfVertices,ver,0,TypTab) ; if( display )print '(6x,"nNod0=",i10)',nNod0
    allocate(xyzOut(geo,nNod0))
    res=GmfGotoKwdF77(ins,GmfVertices)
    select case(ver)
    case(1) !> real(4)
      do iNod0=1,nNod0
        call GmfGetVertex3dr4(ins,xyz(1),xyz(2),xyz(3),mark)
        xyzOut(1:3,iNod0)=xyz(1:3)
      enddo
    case(2) !> real(8)
      do iNod0=1,nNod0
        call GmfGetVertex3dr8(ins,xyzOut(1,iNod0),xyzOut(2,iNod0),xyzOut(3,iNod0),mark)
      enddo
    end select
    
    res=GmfCloseMeshF77(ins)
    if( display )print '(3x,"end Reading")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidReadXYZout3D
  
  subroutine pyramidWriteSolOut3D(title,solOut)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use M_libmesh6_api
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*)                 :: title
    real(8), intent(in), pointer :: solOut(:,:)
    !-
    integer                      :: iOrd,nSol
    integer                      :: i,j
    character(3)                 :: sfx
    integer , parameter          :: iFile=100
    
    character(256)               :: name
    integer                      :: ins,ver,res,geo
    integer                      :: nFld,kind(1)
    real(8) , allocatable        :: sol(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
#if 0==1
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '(/"writing solOut")'
    nSol=size(solOut,1) ; allocate(sol(nSol))
    do iOrd=1,size(solOut,2)
      
      if(   1<=iOrd .and. iOrd<  10 ) write(sfx,'("00",i1)')iOrd
      if(  10<=iOrd .and. iOrd< 100 ) write(sfx,'("0" ,i2)')iOrd
      if( 100<=iOrd .and. iOrd<1000 ) write(sfx,'(     i3)')iOrd
      print '(3x,"saving ",a,"_",a)',title,sfx
      
      call system("ln -fs Pyramids.meshb "// title // "_" //sfx//".meshb")
      
      name=title//"_"//sfx//".sol" ; print '(/"Writing: ",a)',trim(name)
      open(unit=iFile,file=trim(name),action='write')
      write(iFile,'("MeshVersionFormatted 2"/)')
      write(iFile,'("Dimension 3"/)')
      write(iFile,'("SolAtVertices")')
      write(iFile,*)size(solOut,1)
      write(iFile,'("1 1")')
      do i=1,size(solOut,1) ! print '("i=",i6,"/",i6)',i,size(solOut,1)
        write(iFile,*)solOut(i,iOrd)
      enddo
      write(iFile,'(/"End")')
      close(iFile)
    enddo
    print '("end writing solOut")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#else
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '(/"writing solOut")'
    nSol=size(solOut,1) ; allocate(sol(nSol))
    do iOrd=1,size(solOut,2)
      
      if(   1<=iOrd .and. iOrd<  10 ) write(sfx,'("00",i1)')iOrd
      if(  10<=iOrd .and. iOrd< 100 ) write(sfx,'("0" ,i2)')iOrd
      if( 100<=iOrd .and. iOrd<1000 ) write(sfx,'(     i3)')iOrd
      print '(3x,"saving ",a,"_",a)',title,sfx
      
      call system("ln -fs Pyramids.meshb "// title // "_" //sfx//".meshb")
      geo=3
      nFld=1 ; kind(1)=1 ; sol(1:nSol)=solOut(1:nSol,iOrd)
      name=title//"_"//sfx//".solb" ; print '(/"Writing: ",a)',trim(name)
      ver=2
      ins=GmfOpenMeshF77(trim(name),GmfWrite,ver,geo) ; print '(3x,"nSolu=",i10)',nSol
      res=GmfSetKwdF77(ins,GmfSolAtVertices,nSol,nFld,kind(1:1))
      do i=1,nSol ! print '("i=",i6,2x,"sol=",e22.15)',i,sol(i)
        call gmfSetSolAtVertexR8(MshIdx=ins,SolTab=sol(i))
      enddo
      res=GmfCloseMeshF77(ins)
      
    enddo
    deallocate(sol)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#endif
    
    return
  end subroutine pyramidWriteSolOut3D
  
end module basePyramid