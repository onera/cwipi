module baseSimplex1D
  use baseSimplexTools
  implicit none
  contains
  
  subroutine nodes1D(ord, uvw, display)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity triangle
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:)
    logical, intent(in)           :: display
    !---
    integer                       :: i
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Setting uvw
    allocate(uvw(1:ord+1))
    uvw(1:ord+1)=[( (-1d0+2d0*real(i-1,kind=8)/real(ord,kind=8)), i=1,ord+1)]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Points d''interpolation")')
      print '("uvw(",i2,")=",f22.15)',(i,uvw(i),i=1,ord+1)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes1D
  
  subroutine nodes1Dopt(ord, uvw, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:)
    logical, intent(in)           :: display
    !
    integer                       :: i
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call gaussLegendreLobatto(ord=ord,xGLL=uvw)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Points d''interpolation optimisés")')
      print '("uvw(",i2,")=",f22.15)',(i,uvw(i),i=1,ord+1)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes1Dopt
  
  subroutine simplex1D(ord,r,mode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> psi_i(r) = P_i^{0,0}(r) , i=0,ord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => mode(1:np,1:n)
    !> Transpose = False => mode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: r(:)
    real(8), intent(out), pointer :: mode(:,:)
    logical, intent(in)           :: transpose
    !---
    integer                       :: k,n,np
    real(8), pointer              :: jf(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(r) )then
      print '("r not associated")'
      stop "@simplex1D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=ord+1 ; n=size(r)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    if( .not.transpose )then
      allocate(mode(1:n,1:np))
      do k=0,np-1
        call jacobiP(n=k,alpha=0d0,beta=0d0,u=r(1:n),jf=jf)
        mode(1:n,k+1)=jf(1:n)
        deallocate(jf)
      enddo
    else
      allocate(mode(1:np,1:n))
      do k=0,np-1
        call jacobiP(n=k,alpha=0d0,beta=0d0,u=r(1:n),jf=jf)
        mode(k+1,1:n)=jf(1:n)
        deallocate(jf)
      enddo
    endif
    
    return
  end subroutine simplex1D
  
  subroutine gradSimplex1D(ord,r,drMode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => drMode(1:np,1:n)
    !> Transpose = False => drMode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: r(:)
    real(8), intent(out), pointer :: drMode(:,:)
    logical, intent(in)           :: transpose
    !---
    integer                       :: k,n,np
    real(8), pointer              :: jf(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(r) )then
      print '("r not associated")'
      stop "@simplex1D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=ord+1 ; n=size(r)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.transpose )then
      allocate(drMode(1:n,1:np))
      do k=0,np-1
        call djacobiP(n=k,alpha=0d0,beta=0d0,u=r(1:n),jf=jf)
        drMode(1:n,k+1)=jf(1:n)
        deallocate(jf)
      enddo
    else
      allocate(drMode(1:np,1:n))
      do k=0,np-1
        call djacobiP(n=k,alpha=0d0,beta=0d0,u=r(1:n),jf=jf)
        drMode(k+1,1:n)=jf(1:n)
        deallocate(jf)
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine gradSimplex1D
  
  subroutine vandermonde1D(ord,a,vand)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:)
    real(8), intent(out), pointer :: vand(:,:)
    !---
    integer                       :: np
    real(8)                       :: gamma(0:ord+1)
    real(8), pointer              :: jf(:)
    integer                       :: i,j,k
    real(8), pointer              :: x0(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if( .not.associated(a) )then
      print '("ord=",i2)',ord
      print '("Interpolation Points not associated")'
      stop "@vandermonde1D"
    endif
    
    !> Gauss-Legendre Vandermonde matrix
    call simplex1D(ord=ord,r=a,mode=vand,transpose=.false.)
    
    return
  end subroutine vandermonde1D
  
  subroutine gradVandermonde1D(ord,a,dVand)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:)
    real(8), intent(out), pointer :: dVand(:,:)
    !---
    integer                       :: i,j,k
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(a) )then
      print '("ord=",i2)',ord
      print '("Interpolation Points not associated")'
      stop "@vandermonde1D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call gradSimplex1D(ord=ord,r=a,drMode=dVand,transpose=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine gradVandermonde1D
  
  subroutine lagrange1D(ord,uvw,xout,li)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: xout  (  :)
    real(8), intent(in) , pointer :: uvw   (  :)
    real(8), intent(out), pointer :: li    (:,:)
    !--
    integer :: i,j,nVert
    real(8) :: var
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nVert=size(xout)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(li(1:nVert,1:ord)) ; li(1:nVert,1:ord)=1d0
    do i=1,ord+1
      do j=1,ord+1
        if( .not.i==j )then
          var= (uvw(i)-uvw(j))
          var=1d0/var
          li(1:nVert,i)= li(1:nVert,i)          &
          &             *(xout(1:nVert)-uvw(j)) &
          &             *var
        endif
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine lagrange1D
  
  subroutine lagrange1Dv(ord,vand,x,lx,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! lagrange1D := Inverse[Transpose[Vand]].Psi[x];
    ! transpose = .true.  => lx(1:ord+1,1:nPt)
    ! transpose = .false. => lx(1:nPt,1:ord+1)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)            :: ord
    real(8), intent(in)  , pointer :: vand(:,:)
    real(8), intent(in)  , pointer :: x   (:)
    real(8), intent(out) , pointer :: lx  (:,:)
    logical, intent(in)            :: transpose
    !---
    integer                        :: i,j,k,nPt
    real(8)                        :: gamma(0:ord+1)
    integer                        :: iOrd,np
    real(8), pointer               :: Psi(:,:)
    real(8), pointer               :: mat(:,:)
    integer                        :: lWork
    integer, pointer               :: ipiv(:)
    real(8), pointer               :: work(:)
    integer                        :: iErr
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    np=ord+1 ; nPt=size(x)
    
    gamma(0:np) = [(sqrt(real(2*iOrd+1,kind=8)/2d0), iOrd=0,np)]
    
    !> Polynomes de Legendre normalisés :
    !> p_{0}= gamma(0) && p_{1}= gamma(1) x
    !> pour k>= 2 : gamma(k+1) (k+1) p_{k+1} = (2k+1) gamma(k) p_{k} x - k gamma(k-1) p_{k-1}
    allocate(psi(0:np-1,1:nPt))
    psi(0,1:nPt)=gamma(0)
    if( ord>=1 )then
      psi(1,1:nPt)=gamma(1)*x(1:nPt)
      do iOrd=1,ord-1
        psi(iOrd+1,1:nPt)=gamma(iOrd+1)/real(iOrd+1,kind=8)                                 &
        &                 *( real(2*iOrd+1,kind=8)/gamma(iOrd  )*psi(iOrd,1:nPt  )*x(1:nPt) &
        &                   -real(  iOrd  ,kind=8)/gamma(iOrd-1)*psi(iOrd-1,1:nPt)          )
      enddo
    endif
   !call displayMatrix(title="Psi",mat=psi)
    
    !> mat=Inverse[Transpose[Vand]]
    allocate(mat(np,np))
    do i=1,np
      do j=1,np
        mat(i,j)=vand(j,i)
      enddo
    enddo
    
    lWork=64*(np) ; allocate(work(lWork),ipiv(np))
    call dgetrf(np,np,mat(1,1),np,ipiv(1),iErr)
    call dgetri(np,mat(1,1),np,ipiv(1),work(1),lWork,iErr)
    deallocate(ipiv,work)
   !call displayMatrix(title="Inverse[Transpose[Vand]]",mat=mat)
    
    !> lx = Inverse[Transpose[Vand]].Psi
    if( .not.transpose )then
      allocate(lx(1:nPt,1:np)) ; lx(:,:)=0d0
      do i=1,nPt
        do j=1,np
          do k=1,np
            lx(i,j)=lx(i,j)+mat(j,k)*Psi(k-1,i)  ! Attention de bien prendre Psi(k-1,i)
          enddo
        enddo
      enddo
    else
      allocate(lx(1:np,1:nPt)) ; lx(:,:)=0d0
      do i=1,nPt
        do j=1,np
          do k=1,np
            lx(j,i)=lx(j,i)+mat(j,k)*Psi(k-1,i)  ! Attention de bien prendre Psi(k-1,i)
          enddo
        enddo
      enddo
    endif
   !call displayMatrix(title="lx",mat=lx)
    
    deallocate(mat)
    deallocate(psi)
    
    return
  end subroutine lagrange1Dv
    
end module baseSimplex1D