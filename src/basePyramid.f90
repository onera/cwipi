module basePyramid
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use baseSimplexTools
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>  Pyramid ord
  !>  nn=(ord+1)*(ord+2)*(2*ord+3)/6
  !>  ne=3*ord*ord+2
  !>  ni=(ord-1)*(ord-2)*(2*ord-3)/6
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  implicit none
  
contains
  
  
  subroutine pyramiduvw2abc(uvw,a,b,c)
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
    real(8), parameter             :: tol=1d-12
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(uvw,2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(a(1:n),b(1:n),c(1:n))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do i=1,n
      if( uvw(3,i)==1d0 )then
        a(i)=0d0
        b(i)=0d0
        c(i)=1d0
      else
        a(i)=uvw(1,i)/(1d0-uvw(3,i) +tol)
        b(i)=uvw(2,i)/(1d0-uvw(3,i) +tol)
        c(i)=uvw(3,i)
      endif
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramiduvw2abc
  
  subroutine pyramidLagrange3Dv(ord,vand,a,b,c,lx,transpose)
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
    np=(ord+1)*(ord+2)*(2*ord+3)/6 ; nPt=size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call pyramidBasePi(ord=ord,a=a,b=b,c=c,mode=psi,transpose=.false.)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> mat=Transpose[Vand]
    allocate(mat(np,np))
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
    !> lx = Inverse[Transpose[Vand]].Psi
    if( transpose )then
      allocate(lx(1:np,1:nPt)) ; lx(:,:)=0d0
      !> psi(1:nPt,1:np)
      
      do i=1,nPt
        do j=1,np
          do k=1,np
            lx(j,i)=lx(j,i)+mat(j,k)*Psi(i,k)  ! Attention de bien prendre Psi(i,k)
          enddo
        enddo
      enddo
      
    else
      
      allocate(lx(1:nPt,1:np)) ; lx(:,:)=0d0
      ! psi(1:nPt,1:np)
      do i=1,nPt
        do j=1,np
          do k=1,np
            lx(i,j)=lx(i,j)+mat(j,k)*Psi(i,k)  ! Attention de bien prendre Psi(i,k)
          enddo
        enddo
      enddo
    endif
   !call displayMatrix(title="lx",mat=lx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(mat)
    deallocate(psi)
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
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(a) .or. .not.associated(b) .or. .not.associated(c) )then
      print '("3D Nodes a or/and b or/and c not associated")'
      stop "@pyramidVandermonde3D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Vandermonde matrix
    call pyramidBasePi(ord=ord,a=a,b=b,c=c,mode=vand,transpose=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidVandermonde3D
  
  subroutine pyramidGradVandermonde3D(ord,a,b,c,drVand,dsVand,dtVand)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:),b(:),c(:)
    real(8), intent(out), pointer :: drVand(:,:),dsVand(:,:),dtVand(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(a) .or. .not.associated(b) .or. .not.associated(c) )then
      print '("3D Nodes a or/and b or/and c not associated")'
      stop "@pyramidGradVandermonde3D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Grad Vandermonde matrix
    call pyramidGradBasePi(  &
    &    ord=ord            ,&
    &    a=a,b=b,c=c        ,&
    &    drMode=drVand      ,&
    &    dsMode=dsVand      ,&
    &    dtMode=dtVand      ,&
    &    transpose=.false.   )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidGradVandermonde3D
  
  
  
  subroutine pyramidBaseP1(uvw, mode, transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: uvw(:,:)
    real(8), intent(out), pointer :: mode(:,:)
    logical, intent(in)           :: transpose
    !>
    integer                       :: i
    real(8), parameter            :: tol=1d-16
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => mode(1:np,1:n)
    !> Transpose = False => mode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(      transpose )allocate( mode(5,size(uvw,2)) )
    if( .not.transpose )allocate( mode(size(uvw,2),4) )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> vertex functions
    if( transpose )then
      mode(1,:) = .25d0*(1d0-uvw(1,:)-uvw(2,:)-uvw(3,:)+uvw(1,:)*uvw(2,:)/(1d0-uvw(3,:) + tol))
      mode(2,:) = .25d0*(1d0+uvw(1,:)-uvw(2,:)-uvw(3,:)-uvw(1,:)*uvw(2,:)/(1d0-uvw(3,:) + tol))
      mode(3,:) = .25d0*(1d0+uvw(1,:)+uvw(2,:)-uvw(3,:)+uvw(1,:)*uvw(2,:)/(1d0-uvw(3,:) + tol))
      mode(4,:) = .25d0*(1d0-uvw(1,:)+uvw(2,:)-uvw(3,:)-uvw(1,:)*uvw(2,:)/(1d0-uvw(3,:) + tol))
      mode(5,:) = uvw(3,:)
    else
      mode(:,1) = .25d0*(1d0-uvw(1,:)-uvw(2,:)-uvw(3,:)+uvw(1,:)*uvw(2,:)/(1d0-uvw(3,:) + tol))
      mode(:,2) = .25d0*(1d0+uvw(1,:)-uvw(2,:)-uvw(3,:)-uvw(1,:)*uvw(2,:)/(1d0-uvw(3,:) + tol))
      mode(:,3) = .25d0*(1d0+uvw(1,:)+uvw(2,:)-uvw(3,:)+uvw(1,:)*uvw(2,:)/(1d0-uvw(3,:) + tol))
      mode(:,4) = .25d0*(1d0-uvw(1,:)+uvw(2,:)-uvw(3,:)-uvw(1,:)*uvw(2,:)/(1d0-uvw(3,:) + tol))
      mode(:,5) = uvw(3,:)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( transpose )then
      if( 0==0 )then
        do i=1,size(mode,2)
          print '(i3,1x,"uvw=",3(f12.9,1x),1x,"func=",5(f12.9,1x))',i,uvw(1:3,i),mode(1:5,i)
        enddo
      endif
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidBaseP1
  
  subroutine pyramidBasePi(ord,a,b,c,mode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Psi(a,b,c) = P_i^{0,0}(a) P_j^{0,0}(b) P_k^{2*max(i,j)+2,0}(2*c-1) * (1-c)**max(i,j) 
    !> avec {a=u/(1-w) b=v/(1/w) et c=w}
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use baseSimplexTools, only: jacobiP
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:)
    real(8), intent(in) , pointer :: b(:)
    real(8), intent(in) , pointer :: c(:)
    real(8), intent(out), pointer :: mode(:,:)
    logical, intent(in)           :: transpose
    !>
    integer                       :: iu,iv,iw,iM,iNod,n,np
    real(8)                       :: alpha,beta
    real(8), pointer              :: mode1(:),mode2(:),mode3(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => mode(1:np,1:n)
    !> Transpose = False => mode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)*(2*ord+3)/6
    n =size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.transpose )then
      
      allocate(mode(1:n,1:np))
      
      iNod=0
      do iw=0,ord
        do iv=0,ord-iw
          do iu=0,ord-iw
            iNod=iNod+1
            
            iM=max(iu,iv)
            
            alpha=0d0                 ; beta=0d0
            call jacobiP(n=iu,alpha=alpha, beta=beta,u=a(1:n)    ,jf=mode1) ! J_iu^{0             ,0} ( x/(1-z) )
            
            alpha=0d0                 ; beta=0d0
            call jacobiP(n=iv,alpha=alpha, beta=beta,u=b(1:n)    ,jf=mode2) ! J_iv^{0             ,0} ( y/(1-z) )
            
            alpha=real(2*iM+2,kind=8) ; beta=0d0
            call jacobiP(n=iw,alpha=alpha,beta=beta,u=2*c(1:n)-1,jf=mode3) ! J_iw^{2*max(iu,iv)+2,0} ( 2*z-1   )
            
            mode(1:n,iNod)= mode1(1:n)                  &  !>   J_iu^{0             ,0}(x/(1-z))
            &              *mode2(1:n)*(1d0-c(1:n))**iM &  !>  *J_iv^{0             ,0}(y/(1-z))  (1-z)^max(iu,iv)
            &              *mode3(1:n)                     !>  *J_iw^{2*max(iu,iv)+2,0}(2*z-1  )  (2*z-1)
            
          enddo
        enddo
      enddo
      
      deallocate(mode1,mode2,mode3)
      
    else
      
      allocate(mode(1:np,1:n))
      
      iNod=0
      do iw=0,ord
        do iv=0,ord-iw
          do iu=0,ord-iw
            iNod=iNod+1
            
            iM=max(iu,iv)
            call jacobiP(n=iu,alpha=0d0             ,beta=0d0,u=a(1:n),jf=mode1) ! J_iu^{0             ,0}(x/(1-z))
            call jacobiP(n=iv,alpha=0d0             ,beta=0d0,u=b(1:n),jf=mode2) ! J_iv^{0             ,0}(y/(1-z))
            call jacobiP(n=iw,alpha=2d0*real(iM)+2d0,beta=0d0,u=c(1:n),jf=mode3) ! J_iw^{2*max(iu,iv)+2,0}(2*z-1  )
            
            mode(iNod,1:n)= mode1(1:n)                  &  !>   J_iu^{0             ,0}(x/(1-z))
            &              *mode2(1:n)*(1d0-c(1:n))**iM &  !>  *J_iv^{0             ,0}(y/(1-z))  (1-z)^max(iu,iv)
            &              *mode3(1:n)*(2*c(1:n)-1)        !>  *J_iw^{2*max(iu,iv)+2,0}(2*z-1  )  (2*z-1)
            
          enddo
        enddo
      enddo
      
      deallocate(mode1,mode2,mode3)
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidBasePi
  
  
  subroutine pyramidGradBasePi(ord,a,b,c,drMode,dsMode,dtMode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Psi(a,b,c) = P_i^{0,0}(a) P_j^{0,0}(b) (1-c)**max(i,j) P_k^{2*max(i,j)+2,0}(2*c+1)
    !> avec {a=x/(1-z) ; b=y/(1-z) ; c=z}
    
    !> (∂Psi/∂x) = (∂a/∂x) (∂Psi/∂a) + (∂b/∂x) (∂Psi/∂b) + 2 (∂c/∂x) (∂Psi/∂c)
    !> avec  (∂a/∂x)=1/(1-z) et (∂b/∂x)=(∂c/∂x)=0
    !> soit (∂Psi/∂x) =  1/(1-c)  (∂Psi/∂a)
    !>                =  (∂P_i^{0,0}/∂a)(a) P_j^{0,0}(b) (1-c)**(max(i,j)-1) P_k^{2*max(i,j)+2,0}(2*c+1)
    
    !> (∂Psi/∂y) = (∂a/∂y) (∂Psi/∂a) + (∂b/∂y) (∂Psi/∂b) + (∂c/∂y) (∂Psi/∂c)
    !> avec  (∂a/∂y)=0 ; (∂b/∂y)=1/(1-z)=1/(1-c) ; (∂c/∂y)=0
    !> soit (∂Psi/∂x) =  (∂b/∂y) (∂Psi/∂b)
    !>                =  1/(1-c) (∂Psi/∂b)
    !>                =  1/(1-c) P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) (1-c)**max(i,j) P_k^{2*max(i,j)+2,0}(2*c+1)
    !>                =  P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) (1-c)**(max(i,j)-1) P_k^{2*max(i,j)+2,0}(2*c+1)
    
    !> (∂Psi/∂z) = (∂a/∂z) (∂Psi/∂a) + (∂b/∂z) (∂Psi/∂b) + (∂c/∂z) (∂Psi/∂c)
    !> avec  (∂a/∂z)=-x/(1-z)^2=-a/(1-z)= -a/(1-c)
    !>       (∂b/∂y)=-y/(1-z)^2=-b/(1-z)= -b/(1-c)
    !>       (∂c/∂z)= 1
    
    !> (∂a/∂z) (∂Psi/∂a) = -a/(1-c)  (∂P_i^{0,0}/∂a)(a) P_j^{0,0}(b) (1-c)**max(i,j) P_k^{2*max(i,j)+2,0}(2*c+1)
    !>                   = -a  (∂P_i^{0,0}/∂a)(a) P_j^{0,0}(b) (1-c)**(max(i,j)-1) P_k^{2*max(i,j)+2,0}(2*c+1)
    !>                   = -a  (∂P_i^{0,0}/∂a)(a) P_j^{0,0}(b) P_k^{2*max(i,j)+2,0}(2*c+1)  (1-c)**(max(i,j)-1)
    
    !> (∂b/∂z) (∂Psi/∂b) = -b/(1-c)  P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) (1-c)**max(i,j) P_k^{2*max(i,j)+2,0}(2*c+1)
    !>                   = -b P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) (1-c)**(max(i,j)-1) P_k^{2*max(i,j)+2,0}(2*c+1)
    !>                   = -b P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) P_k^{2*max(i,j)+2,0}(2*c+1)  (1-c)**(max(i,j)-1) 
    
    !> (∂c/∂z) (∂Psi/∂c) = (∂Psi/∂c)
    !>                   =   P_i^{0,0}(a) P_j^{0,0}(b) (∂(1-c)**max(i,j)/∂c)   P_k^{2*max(i,j)+2,0}    (2*c+1)
    !>                      +P_i^{0,0}(a) P_j^{0,0}(b)   (1-c)**max(i,j)     (∂P_k^{2*max(i,j)+2,0}/∂c)(2*c+1)  ∂(2*c+1)/∂c
    !>                   =   max(i,j)*P_i^{0,0}(a) P_j^{0,0}(b)   P_k^{2*max(i,j)+2,0}    (2*c+1)  *(1-c)**(max(i,j)-1)
    !>                      +2       *P_i^{0,0}(a) P_j^{0,0}(b) (∂P_k^{2*max(i,j)+2,0}/∂c)(2*c+1)  *(1-c)**(max(i,j)  )
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
    real(8), intent(in) , pointer :: a(:),b(:),c(:)
    real(8), intent(out), pointer :: drMode(:,:)
    real(8), intent(out), pointer :: dsMode(:,:)
    real(8), intent(out), pointer :: dtMode(:,:)
    logical, intent(in)           :: transpose
    !
    real(8)                       :: alpha,beta
    integer                       :: i,n,ad
    integer                       :: iu,iv,iw,iM,iNod
    integer                       :: np
    real(8), pointer              :: fa(:),dfa(:)
    real(8), pointer              :: fb(:),dfb(:)
    real(8), pointer              :: fc(:),dfc(:)
    real(8), pointer              :: tmp(:),tmp1(:)
    real(8)                       :: coef
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '("pyramidGradBasePi")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)*(2*ord+3)/6
    n =size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.transpose )then
      
      allocate(tmp(1:n),tmp1(1:n))
      allocate(drMode(1:n,1:np),dsMode(1:n,1:np),dtMode(1:n,1:np))
      
      iNod=0
      do iw=0,ord ; do iv=0,ord-iw ; do iu=0,ord-iw
        iNod=iNod+1 !  print '(/"pyramidGradBasePi: iNod=",i3)',iNod
        
        iM=max(iu,iv)
        
        alpha=0d0                 ; beta=0d0
        call  jacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n)        ,jf= fa)
        call dJacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n)        ,jf=dfa)
        
        alpha=0d0                 ; beta=0d0
        call  jacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n)        ,jf= fb)
        call dJacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n)        ,jf=dfb)
        
        alpha=real(2*iM+2,kind=8) ; beta=0d0
        call  jacobiP(n=iw,alpha=alpha,beta=beta,u=2d0*c(1:n)-1d0,jf= fc)
        call dJacobiP(n=iw,alpha=alpha,beta=beta,u=2d0*c(1:n)-1d0,jf=dfc) ; dfc=2d0*dfc
        
        
        !> tmp1(1:n)=(1-c)**(max(i,j)-1)
        tmp1(1:n)=(1d0-c(1:n))**(iM-1)
        
        !> (∂Psi/∂x) = (∂P_i^{0,0}/∂a)(a) P_j^{0,0}(b) P_k^{2*max(i,j)+2,0}(2*c+1) (1-c)**(max(i,j)-1)
        tmp(1:n)=dfa(1:n)*fb(1:n)*fc(1:n)
        do i=1,n
          if( .not.tmp(i)==0d0 )tmp(i)=tmp(i)*tmp1(i)
        enddo
        
        !do i=1,n
        !  if( c(i)==1d0 .and. iM<1 )then
        !    print '(/"pyramidGradBasePi: iNod=",i4)',iNod
        !    print '( 3x,"drMode",i4,","i4,")=",f12.5)',i,iNod,drMode(i,iNod)
        !    print '( 3x,"(1-c(",i4,")^(",i4,")=",f12.5)',i,iM-1,tmp1(i)
        !   !print '( 6x,"(1-c(",i2,")^(",i3,")=",f12.5)',(i,iM-1,tmp1(i),i=1,n)
        !  endif
        !enddo
        
        
        if( iM>0 ) tmp(1:n)=tmp(1:n)*tmp1(1:n)
        drMode(1:n,iNod)=tmp(1:n)
        
        
        !> (∂Psi/∂x) =  P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) P_k^{2*max(i,j)+2,0}(2*c+1)  (1-c)**(max(i,j)-1)
        tmp(1:n)=fa(1:n)*dfb(1:n)*fc(1:n)
        do i=1,n
          if( .not.tmp(i)==0d0 )tmp(i)=tmp(i)*tmp1(i)
        enddo
        dsMode(1:n,iNod)=tmp(1:n)
        
        !> (∂Psi/∂z) = (∂a/∂z) (∂Psi/∂a) + (∂b/∂z) (∂Psi/∂b) + (∂c/∂z) (∂Psi/∂c)
        !>           = -a        (∂P_i^{0,0}/∂a)(a)   P_j^{0,0}    (b)   P_k^{2*max(i,j)+2,0}    (2*c+1)  *(1-c)**(max(i,j)-1)
        !>             -b          P_i^{0,0}    (a) (∂P_j^{0,0}/∂b)(b)   P_k^{2*max(i,j)+2,0}    (2*c+1)  *(1-c)**(max(i,j)-1)
        !>             +2       *  P_i^{0,0}    (a)   P_j^{0,0}    (b) (∂P_k^{2*max(i,j)+2,0}/∂c)(2*c+1)  *(1-c)**(max(i,j)  )
        !>             +max(i,j)*  P_i^{0,0}    (a)   P_j^{0,0}    (b)   P_k^{2*max(i,j)+2,0}    (2*c+1)  *(1-c)**(max(i,j)-1)
        
        
        !> dtMode(1:n,iNod) += (∂a/∂z) (∂Psi/∂a)
        tmp(1:n)=-a(1:n)*dfa(1:n)*fb(1:n)*fc(1:n)
        do i=1,n
          if( .not.tmp(i)==0d0 )tmp(i)=tmp(i)*tmp1(i)
        enddo
        dtMode(1:n,iNod)=tmp(1:n)
        
        !> dtMode(1:n,iNod) += +(∂b/∂z) (∂Psi/∂b)
        tmp(1:n)=-b(1:n)*fa(1:n)*dfb(1:n)*fc(1:n)
        do i=1,n
          if( .not.tmp(i)==0d0 )tmp(i)=tmp(i)*tmp1(i)
        enddo
        dtMode(1:n,iNod)=dtMode(1:n,iNod)+tmp(1:n)
        
        !> dtMode(1:n,iNod) += +(∂c/∂z) (∂Psi/∂c)
        tmp(1:n)=2d0*fa(1:n)*fb(1:n)*dfc(1:n)
        tmp(1:n)=tmp(1:n)*(1d0-c(1:n))**iM
        dtMode(1:n,iNod)=dtMode(1:n,iNod)+tmp(1:n)
        
        tmp(1:n)=real(iM,kind=8)*fa(1:n)*fb(1:n)*fc(1:n)
        if( iM>0 ) tmp(1:n)=tmp(1:n)*tmp1(1:n)
        dtMode(1:n,iNod)=dtMode(1:n,iNod)+tmp(1:n)
        
      enddo ; enddo ; enddo
      deallocate(tmp,tmp1)
      
    else
      
      allocate(tmp(1:n),tmp1(1:n))
      allocate(drMode(1:np,1:n),dsMode(1:np,1:n),dtMode(1:np,1:n))
      
      iNod=0
      do iw=0,ord ; do iv=0,ord-iw ; do iu=0,ord-iw
        iNod=iNod+1 !  print '(/"pyramidGradBasePi: iNod=",i3)',iNod
        
        iM=max(iu,iv)
        
        alpha=0d0                 ; beta=0d0
        call  jacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n)        ,jf= fa)
        call dJacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n)        ,jf=dfa)
        
        alpha=0d0                 ; beta=0d0
        call  jacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n)        ,jf= fb)
        call dJacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n)        ,jf=dfb)
        
        alpha=real(2*iM+2,kind=8) ; beta=0d0
        call  jacobiP(n=iw,alpha=alpha,beta=beta,u=2d0*c(1:n)-1d0,jf= fc)
        call dJacobiP(n=iw,alpha=alpha,beta=beta,u=2d0*c(1:n)-1d0,jf=dfc) ; dfc=2d0*dfc
        
        
        !> tmp1(1:n)=(1-c)**(max(i,j)-1)
        tmp1(1:n)=(1d0-c(1:n))**(iM-1)
        
        !> (∂Psi/∂x) = (∂P_i^{0,0}/∂a)(a) P_j^{0,0}(b) P_k^{2*max(i,j)+2,0}(2*c+1) (1-c)**(max(i,j)-1)
        tmp(1:n)=dfa(1:n)*fb(1:n)*fc(1:n)
        do i=1,n
          if( .not.tmp(i)==0d0 )tmp(i)=tmp(i)*tmp1(i)
        enddo
        
        !do i=1,n
        !  if( c(i)==1d0 .and. iM<1 )then
        !    print '(/"pyramidGradBasePi: iNod=",i4)',iNod
        !    print '( 3x,"drMode",i4,","i4,")=",f12.5)',i,iNod,drMode(i,iNod)
        !    print '( 3x,"(1-c(",i4,")^(",i4,")=",f12.5)',i,iM-1,tmp1(i)
        !   !print '( 6x,"(1-c(",i2,")^(",i3,")=",f12.5)',(i,iM-1,tmp1(i),i=1,n)
        !  endif
        !enddo
        
        
        if( iM>0 ) tmp(1:n)=tmp(1:n)*tmp1(1:n)
        drMode(iNod,1:n)=tmp(1:n)
        
        
        !> (∂Psi/∂x) =  P_i^{0,0}(a) (∂P_j^{0,0}/∂b)(b) P_k^{2*max(i,j)+2,0}(2*c+1)  (1-c)**(max(i,j)-1)
        tmp(1:n)=fa(1:n)*dfb(1:n)*fc(1:n)
        do i=1,n
          if( .not.tmp(i)==0d0 )tmp(i)=tmp(i)*tmp1(i)
        enddo
        dsMode(iNod,1:n)=tmp(1:n)
        
        !> (∂Psi/∂z) = (∂a/∂z) (∂Psi/∂a) + (∂b/∂z) (∂Psi/∂b) + (∂c/∂z) (∂Psi/∂c)
        !>           = -a        (∂P_i^{0,0}/∂a)(a)   P_j^{0,0}    (b)   P_k^{2*max(i,j)+2,0}    (2*c+1)  *(1-c)**(max(i,j)-1)
        !>             -b          P_i^{0,0}    (a) (∂P_j^{0,0}/∂b)(b)   P_k^{2*max(i,j)+2,0}    (2*c+1)  *(1-c)**(max(i,j)-1)
        !>             +2       *  P_i^{0,0}    (a)   P_j^{0,0}    (b) (∂P_k^{2*max(i,j)+2,0}/∂c)(2*c+1)  *(1-c)**(max(i,j)  )
        !>             +max(i,j)*  P_i^{0,0}    (a)   P_j^{0,0}    (b)   P_k^{2*max(i,j)+2,0}    (2*c+1)  *(1-c)**(max(i,j)-1)
        
        
        !> dtMode(1:n,iNod) += (∂a/∂z) (∂Psi/∂a)
        tmp(1:n)=-a(1:n)*dfa(1:n)*fb(1:n)*fc(1:n)
        do i=1,n
          if( .not.tmp(i)==0d0 )tmp(i)=tmp(i)*tmp1(i)
        enddo
        dtMode(iNod,1:n)=tmp(1:n)
        
        !> dtMode(1:n,iNod) += +(∂b/∂z) (∂Psi/∂b)
        tmp(1:n)=-b(1:n)*fa(1:n)*dfb(1:n)*fc(1:n)
        do i=1,n
          if( .not.tmp(i)==0d0 )tmp(i)=tmp(i)*tmp1(i)
        enddo
        dtMode(iNod,1:n)=dtMode(iNod,1:n)+tmp(1:n)
        
        !> dtMode(1:n,iNod) += +(∂c/∂z) (∂Psi/∂c)
        tmp(1:n)=2d0*fa(1:n)*fb(1:n)*dfc(1:n)
        tmp(1:n)=tmp(1:n)*(1d0-c(1:n))**iM
        dtMode(iNod,1:n)=dtMode(iNod,1:n)+tmp(1:n)
        
        tmp(1:n)=real(iM,kind=8)*fa(1:n)*fb(1:n)*fc(1:n)
        if( iM>0 ) tmp(1:n)=tmp(1:n)*tmp1(1:n)
        dtMode(iNod,1:n)=dtMode(iNod,1:n)+tmp(1:n)
        
      enddo ; enddo ; enddo
      deallocate(tmp,tmp1)
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '("end pyramidGradBasePi")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidGradBasePi
  
  
  
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
    do iw=0,ord
      
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
  
  subroutine pyramidNodesOpt(ord, uvw, uv, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity pyramid
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)             :: ord
    real(8), intent(in)   , pointer :: uv (:,:) !> Triangle optimized points
    real(8), intent(inout), pointer :: uvw(:,:) !> Pyramide optimized points
    logical, intent(in)             :: display
    !>
    real(8)                         :: xyz(3)
    real(8), pointer                :: dis(:,:) !> Displacement
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
    nNod=0
    do iu=1,ord+1
      nNod=nNod+iu*iu
    enddo
    if( display )print '(/3x,"nDeg=",i6)',nNod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(uvw) )allocate(uvw(3,nNod))
    allocate(dis(3,nNod)) ; dis(1:3,1:nNod)=0d0
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
            dis(1:3,iNod)=[sd2(1,iu),sd2(1,iv),0d0]-uvw(1:3,iNod)
            uvw(1:3,iNod)=[sd2(1,iu),sd2(1,iv),0d0]
          elseif( iv==0 )then !> side2
            nNode=nNode+1
            ad=iu
            dis(1:3,iNod)=sd2(1:3,ad)-uvw(1:3,iNod)
            uvw(1:3,iNod)=sd2(1:3,ad)
          elseif( iu==ord-iw )then !> side3
            !> side2 -> side3 rotation +pi/2 (zz')
            !> [ 0 -1  0] [x]   [-y]
            !> [+1  0  0] [y] = [ x]
            !> [ 0  0 +1] [z]   [ z]
            nNode=nNode+1
            ad=iv
            dis(1:3,iNod)=[-sd2(2,ad),+sd2(1,ad),sd2(3,ad)]-uvw(1:3,iNod)
            uvw(1:3,iNod)=[-sd2(2,ad),+sd2(1,ad),sd2(3,ad)]
          elseif( iv==ord-iw )then !> side4
            !> side2 -> side4 rotation +pi (zz')
            !> [-1  0  0] [x]   [-x]
            !> [ 0 -1  0] [y] = [-y]
            !> [ 0  0 +1] [z]   [ z]
            nNode=nNode+1
            ad=ord-iw-iu
            dis(1:3,iNod)=[-sd2(1,ad),-sd2(2,ad),sd2(3,ad)]-uvw(1:3,iNod)
            uvw(1:3,iNod)=[-sd2(1,ad),-sd2(2,ad),sd2(3,ad)]
          elseif( iu==0 )then !> side5
            !> side2 -> side5 rotation 3pi/2 (zz')
            !> [ 0 +1  0] [x]   [+y]
            !> [-1  0  0] [y] = [-x]
            !> [ 0  0 +1] [z]   [ z]
            nNode=nNode+1
            ad=ord-iw-iv
            dis(1:3,iNod)=[sd2(2,ad),-sd2(1,ad),sd2(3,ad)]-uvw(1:3,iNod)
            uvw(1:3,iNod)=[sd2(2,ad),-sd2(1,ad),sd2(3,ad)]
          endif
          
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Traitement hors faces
    if( display )print '(3x,"Nodes optimization inside pyramid")'
    
    nNodi=0
    iNod=0
    do iw=0,ord
      if( display )print '(6x,"level:",i3)',iw
      do iv=0,ord-iw
        do iu=0,ord-iw
          iNod=iNod+1
          
          if(.not.(iw==0.or.iv==0.or.iv==ord-iw.or.iu==0 .or.iu==ord-iw) )then
            nNodi=nNodi+1
            
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
            
            dis(1:3,iNod)=xyz(1:3)-uvw(1:3,iNod)
            uvw(1:3,iNod)=xyz(1:3)
            
          endif
          
        enddo
        iNod0=iNod1
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !if( display )then
    if( 0==1 )then
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
   !if( display )then
    if( 0==1 )then
      print '(3x,"displacement")'
      iNod=0
      do iw=0,ord
        print '(6x,"level: ",i3)',iw
        do iv=0,ord-iw ; do iu=0,ord-iw
          iNod=iNod+1
          print '(6x,"dis(",i6,")=",3(f12.9,1x))',iNod,dis(1:3,iNod)
        enddo ; enddo
      enddo
      print '(3x,"end")'
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(3x,"nn=",i4,3x,"ne=",i4,3x,"ni=",i4)',nNod,nNode,nNodi
   !if( display )print '(3x,"nn=",i4,3x,"ne=",i4,3x,"ni=",i4)',(ord+1)*(ord+2)*(2*ord+3)/6,3*ord*ord+2,(ord-1)*(ord-2)*(2*ord-3)/6
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
  
  subroutine pyramidSides3D(ord, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity pyramid
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
    integer, intent(in)    :: ord
    logical, intent(in)    :: display
    !---
    integer                :: iu,iv,iw
    integer                :: iNod,iSide
    integer                :: side((ord+1)*(ord+1)+4*(ord+1)*(ord+2)/2)
    integer                :: sideIdx(6)
    integer                :: idx((ord+1)*(ord+2)/2)
    integer                :: nod((ord+1)*(ord+2)/2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(/"Extract Pyramid List of DOF per Side")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sideIdx(1)=0
    sideIdx(2)=sideIdx(1)+(ord+1)*(ord+1)   !> square
    sideIdx(3)=sideIdx(2)+(ord+1)*(ord+2)/2 !> triangle
    sideIdx(4)=sideIdx(3)+(ord+1)*(ord+2)/2 !> triangle
    sideIdx(5)=sideIdx(4)+(ord+1)*(ord+2)/2 !> triangle
    sideIdx(6)=sideIdx(5)+(ord+1)*(ord+2)/2 !> triangle
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iNod=0
    do iw=0,ord ; do iv=0,ord-iw ; do iu=0,ord-iw
      iNod=iNod+1
      
      !> side1 iw=0
      if( iw==0 )then
        sideIdx(1)=sideIdx(1)+1
        side( sideIdx(1) )=iNod
      endif
      
      !> side2 iv=0
      if( iv==0 )then
        sideIdx(2)=sideIdx(2)+1
        side( sideIdx(2) )=iNod
      endif
      
      !> side3 iu=ord
      if( iu==ord-iw )then
        sideIdx(3)=sideIdx(3)+1
        side( sideIdx(3) )=iNod
      endif
      
      !> side4 iv=ord
      if( iv==ord-iw )then
        sideIdx(4)=sideIdx(4)+1
        side( sideIdx(4) )=iNod
      endif
      
      !> side5 iu=0
      if( iu==0 )then
        sideIdx(5)=sideIdx(5)+1
        side( sideIdx(5) )=iNod
      endif
    enddo ; enddo ; enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sideIdx(1)=0
    sideIdx(2)=sideIdx(1)+(ord+1)*(ord+1)
    sideIdx(3)=sideIdx(2)+(ord+1)*(ord+2)/2
    sideIdx(4)=sideIdx(3)+(ord+1)*(ord+2)/2
    sideIdx(5)=sideIdx(4)+(ord+1)*(ord+2)/2
    sideIdx(6)=sideIdx(5)+(ord+1)*(ord+2)/2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Les faces 4 et 5 doivent être retournées
    call permutation(order=ord, move=+1, flip=.true., dg=idx)
    
    do iSide=4,5
      nod(1:(ord+1)*(ord+2)/2)=side(sideIdx(iSide)+1:sideIdx(iSide+1) )
      do iNod=1,(ord+1)*(ord+2)/2
        side( sideIdx(iSide)+iNod )=nod( idx(iNod) )
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      do iSide=1,5
        print '(3x,"Side",i1,": ",$)',iSide
        do iNod=sideIdx(iSide)+1,sideIdx(iSide+1)
          print '(i5,1x,$)',side(iNod)
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
  
  subroutine pyramidReadXYZout3D(xyzOut)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(out), pointer :: xyzOut (:,:)
    !-
    integer :: i,j,nVert
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    open(unit=10,file='PyramidP100.mesh',status='old',action='read')
    do i=1,5 ; read(10,*) ; enddo
    read(10,*)nVert
    allocate(xyzOut(3,nVert))
    do i=1,nVert
      read(10,*)xyzOut(1:3,i)
    enddo
    close(10)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidReadXYZout3D
  
  subroutine pyramidWriteSolOut3D(title,solOut)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    character(*)                 :: title
    real(8), intent(in), pointer :: solOut(:,:)
    !-
    integer                      :: iOrd
    integer                      :: i,j
    character(3)                 :: sfx
    integer , parameter          :: iFile=100
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '(/"writing solOut")'
    do iOrd=1,size(solOut,2)
      
      if(   1<=iOrd .and. iOrd<  10 ) write(sfx,'("00",i1)')iOrd
      if(  10<=iOrd .and. iOrd< 100 ) write(sfx,'("0" ,i2)')iOrd
      if( 100<=iOrd .and. iOrd<1000 ) write(sfx,'(     i3)')iOrd
      print '(3x,"saving ",a,"_",a)',title,sfx
      
      call system("ln -fs PyramidP100.mesh "// title // "_" //sfx//".mesh")
      
      open(unit=iFile,file=title//"_"//sfx//".sol",status='unknown',action='write')
      write(iFile,'("MeshVersionFormatted 2"/)')
      write(iFile,'("Dimension 3"/)')
      write(iFile,'("SolAtVertices")')
      write(iFile,*)size(solOut,1)
      write(iFile,'("1 1")')
      do i=1,size(solOut,1)
        write(iFile,*)solOut(i,iOrd)
      enddo
      write(iFile,'(/"End")')
      close(iFile)
      
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidWriteSolOut3D
  
end module basePyramid