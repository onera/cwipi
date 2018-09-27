module baseSimplex1D
  use baseSimplexTools
  implicit none
  contains
  
  subroutine setL2MeshIJK(meshOrder,ijk)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)    :: meshOrder
    integer, intent(inout) :: ijk(:)
    !>
    integer                :: iMod,nMod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> 01 03 03 05 02
    nMod=(meshOrder+1)
    
    ijk(1)=0
    ijk(2)=meshOrder
    do iMod=3,nMod
      ijk(iMod)=iMod-2
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    
    return
  end subroutine setL2MeshIJK
  
  subroutine setL2BasisEqui(ord,ijk,uvw,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Calcul des base d'ordere ord pour des triangles dont les points
    ! d'interpolation sont equidistants.
    ! Numerotation des sommets suivant ijk
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use baseSimplexTools, only: monomialProduct
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)    :: ord 
    integer, intent(in)    :: ijk(:)       !> ijk(1:nMod)
    real(8), intent(in)    :: uvw(:)       !> uvw(1:nNod)
    real(8), intent(inout) :: ai(:,:)      !> ai(1:nMod,1:nNod)
    !>
    real(8), pointer       :: u(:),v(:)
    integer                :: iMod,nMod    !> nMod=(ord+1)*(ord+2)/2
    integer                :: iNod,nNod    !> nNod=size(u) 
    integer                :: iu,iv,iw
    integer                :: i
    real(8), pointer       :: fu(:),fv(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> bases de Lagrange
    !print '(/"calcul des bases Triangle P",i1)',ord
    
    nNod=size(uvw)
    allocate( u(1:nNod), v(1:nNod))
    allocate(fu(1:nNod),fv(1:nNod))
    do iNod=1,nNod
      u(iNod)=(uvw(iNod)-1d0)*5d-1 !> u \in [-1,+1] => u \in [0,1]
      v(iNod)=1d0-u(iNod)
    enddo
    
    !print '("size(ai)=",i2,"x",i2)',size(ai,1),size(ai,2)
    nMod=(ord+1)
    do iMod=1,nMod
      iu=ijk(iMod)
      iv=1d0-iu
      
      call monomialProduct(ord=ord,n=iu,u=u, fn=fu)
      call monomialProduct(ord=ord,n=iv,u=v, fn=fv)
      
      !print '(3x,"iMod=",i3,3x,"iu=",i3," iv=",i3)',iMod,iu,iv
      do iNod=1,nNod
        ai(iMod,iNod)=fu(iNod)*fv(iNod)
      enddo
      
    enddo
    
    deallocate(u,v,fu,fv)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine setL2BasisEqui  
  
  
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
    !> Transpose = True  => mode(1:nMod,1:nNod)
    !> Transpose = False => mode(1:nNod,1:nMod)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: r(:)
    real(8), intent(out), pointer :: mode(:,:)
    logical, intent(in)           :: transpose
    !>
    integer                       :: nMod,nNod
    integer                       :: k
    real(8), pointer              :: jf(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(r) )then
      print '("r not associated")'
      stop "@simplex1D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nMod=ord+1 ; nNod=size(r)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( transpose )then
      allocate(mode(1:nMod,1:nNod))
      do k=0,nMod-1
        call jacobiP(n=k,alpha=0d0,beta=0d0,u=r(1:nNod),jf=jf)
        mode(k+1,1:nNod)=jf(1:nNod)
        deallocate(jf)
      enddo
    else
      allocate(mode(1:nNod,1:nMod))
      do k=0,nMod-1
        call jacobiP(n=k,alpha=0d0,beta=0d0,u=r(1:nNod),jf=jf)
        mode(1:nNod,k+1)=jf(1:nNod)
        deallocate(jf)
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine simplex1D
  
  subroutine gradSimplex1D(ord,r,drMode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => drMode(1:nMod,1:nNod)
    !> Transpose = False => drMode(1:nNod,1:nMod)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: r(:)
    real(8), intent(out), pointer :: drMode(:,:)
    logical, intent(in)           :: transpose
    !>
    integer                       :: k,nMod,nNod
    real(8), pointer              :: jf(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(r) )then
      print '("r not associated")'
      stop "@simplex1D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nMod=ord+1 ; nNod=size(r)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.transpose )then
      allocate(drMode(1:nNod,1:nMod))
      do k=0,nMod-1
        call djacobiP(n=k,alpha=0d0,beta=0d0,u=r(1:nNod),jf=jf)
        drMode(1:nNod,k+1)=jf(1:nNod)
        deallocate(jf)
      enddo
    else
      allocate(drMode(1:nMod,1:nNod))
      do k=0,nMod-1
        call djacobiP(n=k,alpha=0d0,beta=0d0,u=r(1:nNod),jf=jf)
        drMode(k+1,1:nNod)=jf(1:nNod)
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
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(a) )then
      print '("ord=",i2)',ord
      print '("Interpolation Points not associated")'
      stop "@vandermonde1D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> V(i,j) = Phi_j(xi_i) 
    !> vand(1:nNod,1:nMod) => Transpose=.false.
    call simplex1D(ord=ord,r=a,mode=vand,transpose=.false.) !> (nNod,nMod)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
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
   !allocate(li(1:nVert,1:ord)) ; li(1:nVert,1:ord)=1d0
    li(:,:)=1d0
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
    integer, intent(in)             :: ord
    real(8), intent(in)   , pointer :: vand(:,:)
    real(8), intent(in)   , pointer :: x   (:)
    real(8), intent(inout), pointer :: lx  (:,:)
    logical, intent(in)             :: transpose
    !---
    integer                         :: i,j,k,nPt
    real(8)                         :: gamma(0:ord+1)
    integer                         :: iOrd,np
    real(8), pointer                :: psi(:,:)
    real(8), pointer                :: mat(:,:)
    integer                         :: lWork
    integer, pointer                :: ipiv(:)
    real(8), pointer                :: work(:)
    integer                         :: iErr
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=ord+1 ; nPt=size(x)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    gamma(0:np) = [(sqrt(real(2*iOrd+1,kind=8)/2d0), iOrd=0,np)]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
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
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    !> lx = Inverse[Transpose[Vand]].psi
    if( .not.transpose )then
     !allocate(lx(1:nPt,1:np)) ; lx(:,:)=0d0
      lx(:,:)=0d0
      do i=1,nPt
        do j=1,np
          do k=1,np
            lx(i,j)=lx(i,j)+mat(j,k)*psi(k-1,i)  ! Attention de bien prendre psi(k-1,i)
          enddo
        enddo
      enddo
    else
     !allocate(lx(1:np,1:nPt)) ; lx(:,:)=0d0
      lx(:,:)=0d0
      do i=1,nPt
        do j=1,np
          do k=1,np
            lx(j,i)=lx(j,i)+mat(j,k)*psi(k-1,i)  ! Attention de bien prendre psi(k-1,i)
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
  end subroutine lagrange1Dv
    
end module baseSimplex1D