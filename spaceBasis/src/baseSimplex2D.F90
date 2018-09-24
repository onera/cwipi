
module baseSimplex2D
  use baseSimplexTools
  implicit none 
  
  interface nodes2Duv2rs ; module procedure nodes2Duv2rs_0 ; end interface
  interface nodes2Duv2rs ; module procedure nodes2Duv2rs_1 ; end interface
  
  interface nodes2Duv2ab ; module procedure nodes2Duv2ab   ; end interface
  interface nodes2Duv2ab ; module procedure nodes2Duv2ab_0 ; end interface
  interface nodes2Duv2ab ; module procedure nodes2Duv2ab_1 ; end interface
  
  interface lagrange2Dv  ; module procedure lagrange2Dv_1  ; end interface
  interface lagrange2Dv  ; module procedure lagrange2Dv_2  ; end interface
  
  contains

  subroutine free_double_c(array, s_array)  BIND(C, name="SNB_free_double")
    use, intrinsic :: ISO_C_BINDING 
    implicit none 
    type (C_PTR), value    :: array
    integer (C_INT), value :: s_array

    real(8), pointer       :: array_f  !(:)

    integer                :: s_array_f

    s_array_f = s_array
    
   !call c_f_pointer (array, array_f, (/s_array_f/) )
    call c_f_pointer (array, array_f)
    
    deallocate(array_f)
    
  end subroutine free_double_c
    
  subroutine nodes2D_c(ord, uvw, display)  BIND(C, name="SNB_nodes2D") 
    use, intrinsic :: ISO_C_BINDING 
    implicit none 
    integer (C_INT), value :: ord 
    type (C_PTR)           :: uvw
    integer (C_INT), value :: display

    real(8), pointer :: uvw_f(:,:)
    integer          :: ord_f
    logical          :: display_f

    ord_f = ord
    display_f = (display == 0)

    call nodes2D (ord_f, uvw_f, display_f)

    uvw = c_loc (uvw_f)

  end subroutine nodes2D_c

  Subroutine nodes2D(ord, uvw,display)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity triangle
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    logical, intent(in)           :: display
    real(8), pointer :: uvw(:,:)
    !---
    integer                       :: iu,iv,iw,ad
    integer                       :: m,n
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Total number of nodes
    n=(ord+1)*(ord+2)/2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Create equidistributed nodes on unity triangle
    allocate(uvw(1:3,1:n))
    if( ord==0 )then
      uvw(1:3,1)=[1d0/3d0,1d0/3d0,1d0/3d0]
    elseif( ord==1 )then
      uvw(1:3,1)=[0d0,0d0,1d0]
      uvw(1:3,2)=[1d0,0d0,0d0]
      uvw(1:3,3)=[0d0,1d0,0d0]
    elseif( ord==2 )then
      uvw(1:3,1)=[0.0d0, 0.0d0, 1.0d0]
      uvw(1:3,2)=[0.5d0, 0.0d0, 0.5d0]
      uvw(1:3,3)=[1.0d0, 0.0d0, 0.0d0]
      uvw(1:3,4)=[0.0d0, 0.5d0, 0.5d0]
      uvw(1:3,5)=[0.5d0, 0.5d0, 0.0d0]
      uvw(1:3,6)=[0.0d0, 1.0d0, 0.0d0]
    else
      do iu=0,ord
        do iv=0,ord-iu
          do iw=0,ord-iu-iv
            ad=iu+iv*(ord+1)-(iv*(iv-1))/2 +1 !> Rangement façon space            
            uvw(1:3,ad)=[real(iu,kind=8)/real(ord,kind=8),& !> u
            &            real(iv,kind=8)/real(ord,kind=8),& !> v
            &            real(iw,kind=8)/real(ord,kind=8) ] !> w
          enddo
        enddo
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Triangle unité initial:")')
      print '("ad=",i5,2x,"u=",f19.16,2x,"v=",f19.16,2x,"w=",f19.16)',(ad,uvw(1:3,ad),ad=1,n)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes2D
  
  subroutine nodes2Dopt(ord,uvw,display)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity triangle
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)             :: ord
    real(8), intent(inout), pointer :: uvw(:,:)
    logical, intent(in)             :: display
    !>
    real(8), pointer                :: x(:),y(:)
    !--
    real(8)                         :: alpha
    integer                         :: iu,iv,iw,ad
    integer                         :: m,n
    real(8), pointer                :: xGLL(:)
    real(8), pointer                :: xout(:)
    real(8), pointer                :: l1(:),blend1(:),warpFactor1(:)
    real(8), pointer                :: l2(:),blend2(:),warpFactor2(:)
    real(8), pointer                :: l3(:),blend3(:),warpFactor3(:)
    !---
    real(8), parameter   :: alpOpt(15)=[0.0000d0, 0.0000d0, 1.4152d0, 0.1001d0, 0.2751d0,&
    &                                   0.9800d0, 1.0999d0, 1.2832d0, 1.3648d0, 1.4773d0,&
    &                                   1.4959d0, 1.5743d0, 1.5770d0, 1.6223d0, 1.6258d0 ]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Total number of nodes
    n=(ord+1)*(ord+2)/2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( ord>2 )then
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> Create barycentric coordinate l1,l2,l3
      allocate(l1(1:n),l2(1:n),l3(1:n))
      l1(1:n)=uvw(2,1:n) ! v     =  (s+1)/2
      l2(1:n)=uvw(3,1:n) ! 1-u-v = -(r+s)/2
      l3(1:n)=uvw(1,1:n) ! u    =   (r+1)/2
      !write(*,'(/"Triangle unité initial:")')
      !print '("ad=",i2,2x,"u=",f12.5,2x,"v=",f12.5,2x,"w=",f12.5)',(ad,l3(ad),l1(ad),l2(ad),ad=1,n)
      !print '("ad=",i2,2x,"l1=",f12.5)',(ad,l1(ad),ad=1,n)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> Passage triangle unité -> triangle equilateral
      !>   triangle unité (0,0)-(1,0)-(0,1)
      !>   triangle equilateral (-1,-1/sqrt(3))-(1,-1/sqrt(3))-(0,2/sqrt(3)) 
      allocate(x(1:n)) ; x(1:n)=(-l2(1:n)+l3(1:n)            )           ! -w+u
      allocate(y(1:n)) ; y(1:n)=(-l2(1:n)-l3(1:n)+2d0*l1(1:n))/sqr3      ! -w-u+2/sqrt(3) v
      !write(*,'(/"Triangle equilateral initial:")')
      !print '("ad=",i2,2x,"x=",f12.5,2x,"y=",f12.5)',(ad,x(ad),y(ad),ad=1,n)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> alpha
      if( ord<16 )then
        alpha=alpOpt(ord)
      else
        alpha=5d0/3d0
      endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> compute blending function at each node for each edge
      allocate(blend1(1:n)) ; blend1(1:n)=4d0*l2(1:n)*l3(1:n)
      allocate(blend2(1:n)) ; blend2(1:n)=4d0*l1(1:n)*l3(1:n)
      allocate(blend3(1:n)) ; blend3(1:n)=4d0*l1(1:n)*l2(1:n)
      !write(*,'()')
      !print '("ad=",i2,2x,"blends=",3(f12.5,1x))',(ad,blend1(ad),blend2(ad),blend3(ad),ad=1,n)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> compute Gauss-Lobatto-Legendre node distribution
      call gaussLegendreLobatto(ord=ord,xGLL=xGLL)
      !write(*,'()')
      !print '("xGLL(",i2,")=",f22.15)',(ad,xGLL(ad),ad=1,ord+1)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> amount of warp for each node, for each edge
      allocate(xout(1:n))
      xout(1:n)=l3(1:n)-l2(1:n) ; call warpFactor(ord=ord,xnodes=xGLL,xout=xout,warp=warpFactor1)
      xout(1:n)=l1(1:n)-l3(1:n) ; call warpFactor(ord=ord,xnodes=xGLL,xout=xout,warp=warpFactor2)
      xout(1:n)=l2(1:n)-l1(1:n) ; call warpFactor(ord=ord,xnodes=xGLL,xout=xout,warp=warpFactor3)
      deallocate(xout,xGLL)
      !write(*,'()')
      !write(*,'("ad=",i2,2x,"xout=",3(f12.5,1x))') (ad,l3(ad)-l2(ad),l1(ad)-l3(ad),l2(ad)-l1(ad),ad=1,n)
      !write(*,'()')
      !write(*,'("ad=",i2,2x,"warpf=",3(f12.5,1x))') (ad,warpFactor1(ad),warpFactor2(ad),warpFactor3(ad),ad=1,n)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> combine blend & warp
      warpFactor1(1:n)=blend1(1:n)*warpFactor1(1:n)*(1d0+(alpha*l1(1:n))**2)
      warpFactor2(1:n)=blend2(1:n)*warpFactor2(1:n)*(1d0+(alpha*l2(1:n))**2)
      warpFactor3(1:n)=blend3(1:n)*warpFactor3(1:n)*(1d0+(alpha*l3(1:n))**2)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> accumulate deformations associated with each edge
      x(1:n)=x(1:n) + 1d0*warpFactor1(1:n) + cos(2d0*pi/3d0)*warpFactor2(1:n) + cos(4d0*pi/3d0)*warpFactor3(1:n)
      y(1:n)=y(1:n) + 0d0*warpFactor1(1:n) + sin(2d0*pi/3d0)*warpFactor2(1:n) + sin(4d0*pi/3d0)*warpFactor3(1:n)
      !write(*,'(/"Triangle equilateral final:")')
      !print '("ad=",i2,2x,"x=",f12.5,2x,"y=",f12.5)',(ad,x(ad),y(ad),ad=1,n)
      !print '("ad=",i2,2x,"x=",f12.5)',(ad,x(ad),ad=1,n)
      !print '("ad=",i2,2x,"y=",f12.5)',(ad,y(ad),ad=1,n)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> Passage triangle equilateral -> triangle unité
      l1(1:n)=(            sqr3*y(1:n)+1d0)/3d0
      l2(1:n)=(-3d0*x(1:n)-sqr3*y(1:n)+2d0)/6d0
      l3(1:n)=( 3d0*x(1:n)-sqr3*y(1:n)+2d0)/6d0
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      uvw(1,1:n)=l3(1:n)
      uvw(2,1:n)=l1(1:n)
      uvw(3,1:n)=l2(1:n)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      deallocate(l1,blend1,warpFactor1)
      deallocate(l2,blend2,warpFactor2)
      deallocate(l3,blend3,warpFactor3)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Triangle unité optimisé")')
      !print '("uvw(1:3,",i4,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,uvw(:,ad),ad=1,n)
      print '("ad=",i5,2x,"u=",f19.16,2x,"v=",f19.16,2x,"w=",f19.16)',(ad,uvw(1:3,ad),ad=1,n)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes2Dopt
  
  subroutine nodes2Duv2rs_0(uv,rs,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Triangle (0,0)-(1,0)-(0,1) To Triangle (-1,-1)-(1,-1)-(1,1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: uv(:,:)
    real(8), intent(out), pointer :: rs(:,:)
    logical, intent(in)           :: display
    !
    integer                       :: i
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(rs(1:2,1:size(uv,2)))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    rs(1:2,:)=2d0*uv(1:2,:)-1d0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Triangle (uv2rs):")')
      print '("r,s(",i4,")=",f19.16,2x,f19.16)',(i,rs(1:2,i),i=1,size(uv,2))
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes2Duv2rs_0
  
  subroutine nodes2Duv2rs_1(u,v,rs,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Triangle (0,0)-(1,0)-(0,1) To Triangle (-1,-1)-(1,-1)-(1,1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: u(:),v(:)
    real(8), intent(out), pointer :: rs(:,:)
    logical, intent(in)           :: display
    !
    integer                       :: i
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(rs(1:2,1:size(u)))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do i=1,size(u)
      rs(1:2,i)=2d0*[u(i),v(i)]-1d0
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Triangle (uv2rs):")')
      print '("r,s(",i4,")=",f19.16,2x,f19.16)',(i,rs(1:2,i),i=1,size(u))
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes2Duv2rs_1
  

  subroutine nodes2Duv2ab_c(nVtx, uv, a, b, display)  BIND(C, name="SNB_nodes2Duv2ab") 
    use, intrinsic :: ISO_C_BINDING 
    implicit none 

    integer (C_INT), value :: nVtx 
    type (C_PTR),    value :: uv
    type (C_PTR)           :: a
    type (C_PTR)           :: b
    integer (C_INT), value :: display

    real(8), pointer   :: uv_f(:,:)
    real(8), pointer   :: a_f(:),b_f(:)
    logical            :: display_f

    integer          :: n

    display_f = (display == 1)
    call c_f_pointer (uv, uv_f, (/2,nVtx/) )
    
    call nodes2Duv2ab (uv_f, a_f, b_f, display_f)

    a = c_loc (a_f)
    b = c_loc (b_f)
    
  end subroutine nodes2Duv2ab_c
  
  subroutine nodes2Duv2ab(u,v,a,b,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Triangle (0,0)-(1,0)-(0,1) To Triangle (-1,-1)-(1,-1)-(1,1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)           :: u,v
    real(8), intent(out), pointer :: a(:),b(:)
    logical, intent(in)           :: display
    !
    integer                       :: i,n
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=1 ; allocate(a(1),b(1))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    a(1)=2d0*u-1d0
    b(1)=2d0*v-1d0
    
    if( b(1)==1d0 )then ! Hesthaven p174
      a(1)=-1d0
    else
      a(1)=2d0*(1d0+a(1))/(1d0-b(1))-1d0
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Triangle (uv2ab):")')
      print '("a,b=",f12.5,2x,f12.5)',a(1),b(1)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes2Duv2ab
  
  subroutine nodes2Duv2ab_0(uv,a,b,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Triangle (0,0)-(1,0)-(0,1) To Triangle (-1,-1)-(1,-1)-(1,1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)           :: uv(:,:)
    real(8), intent(out), pointer :: a(:),b(:)
    logical, intent(in)           :: display
    !
    integer                       :: i,n
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(uv,2) ; allocate(a(1:n),b(1:n))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do i=1,n
      a(i)=2d0*uv(1,i)-1d0
      b(i)=2d0*uv(2,i)-1d0
    enddo
    
    do i=1,n
      if( b(i)==1d0 )then ! Hesthaven p174
        a(i)=-1d0
      else
        a(i)=2d0*(1d0+a(i))/(1d0-b(i))-1d0
      endif
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Triangle (uv2ab):")')
      print '("a,b(",i4,")=",f12.5,2x,f12.5)',(i,a(i),b(i),i=1,size(uv,2))
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes2Duv2ab_0
  
  subroutine nodes2Duv2ab_1(u,v,a,b,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Triangle (0,0)-(1,0)-(0,1) To Triangle (-1,-1)-(1,-1)-(1,1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in)           :: u(:),v(:)
    real(8), intent(out), pointer :: a(:),b(:)
    logical, intent(in)           :: display
    !
    integer                       :: i,n
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(u) ; allocate(a(1:n),b(1:n))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do i=1,n
      a(i)=2d0*u(i)-1d0
      b(i)=2d0*v(i)-1d0      
    enddo
        
    do i=1,n
      if( b(i)==1d0 )then ! Hesthaven p174
        a(i)=-1d0
      else
        a(i)=2d0*(1d0+a(i))/(1d0-b(i))-1d0
      endif
    enddo    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Triangle (uv2ab):")')
      print '("a,b(",i4,")=",f12.5,2x,f12.5)',(i,a(i),b(i),i=1,size(u))
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes2Duv2ab_1
  
  subroutine nodes2Drs2ab(rs,a,b,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: rs(:,:)
    real(8), intent(out), pointer :: a(:)
    real(8), intent(out), pointer :: b(:)
    logical, intent(in)           :: display
    !
    integer                       :: i,n
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(rs,2) ; allocate(a(1:n),b(1:n))
    
    do i=1,n
      if( rs(2,i)==1d0 )then ! Hesthaven p174
        a(i)=-1d0
      else
        a(i)=2d0*(1d0+rs(1,i))/(1d0-rs(2,i))-1d0
      endif
    enddo
    
    do i=1,n
      b(i)=rs(2,i)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Triangle (rs2ab):")')
      print '("a,b(",i4,")=",f12.5,2x,f12.5)',(i,a(i),b(i),i=1,size(a))
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    return
  end subroutine nodes2Drs2ab
  
  subroutine simplex2D(ord,a,b,mode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> psi_{ij}(a,b) = sqr2 P_i^{0,0}(a) P_j^{2i+1,0}(b) (1-b)**i  avec {a=2(1+r)/(1-s)-1 et b=s}
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => mode(1:np,1:n)
    !> Transpose = False => mode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
   !real(8), intent(in) , pointer :: a(:),b(:)
    real(8), intent(in)           :: a(:),b(:)
    real(8), intent(out), pointer :: mode(:,:)
    logical, intent(in)           :: transpose
    !---
    integer                       :: i,n,ad
    integer                       :: iu,iv,iw
    integer                       :: np
    real(8), pointer              :: mode1(:),mode2(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( size(a)==0.or.size(b)==0 )then
      print '("a or/and b not associated")'
      stop "@simplex2D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)/2 ; n=size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.transpose )then
      
      allocate(mode(1:n,1:np))
      do iu=0,ord
        call jacobiP(n=iu,alpha=0d0,beta=0d0,u=a(1:n),jf=mode1)
        do iv=0,ord-iu
          ad=iu+iv*(ord+1)-(iv*(iv-1))/2+1 ! Rangement façon space
          call jacobiP(n=iv,alpha=2d0*real(iu,kind=8)+1d0,beta=0d0,u=b(1:n),jf=mode2)
          mode(1:n,ad)=sqr2*mode1(1:n)*mode2(1:n)*(1d0-b(1:n))**iu
          deallocate(mode2)
        enddo
        deallocate(mode1)
      enddo
      
    else
      
      allocate(mode(1:np,1:n))
      do iu=0,ord
        call jacobiP(n=iu,alpha=0d0      ,beta=0d0,u=a(:),jf=mode1)
        do iv=0,ord-iu
          ad=iu+iv*(ord+1)-(iv*(iv-1))/2+1 ! Rangement façon space
          call jacobiP(n=iv,alpha=2d0*real(iu,kind=8)+1d0,beta=0d0,u=b(:),jf=mode2)
          mode(ad,1:n)=sqr2*mode1(1:n)*mode2(1:n)*(1d0-b(1:n))**iu
          deallocate(mode2)
        enddo
        deallocate(mode1)
      enddo
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine simplex2D
  
  subroutine gradSimplex2D(ord,a,b,drMode,dsMode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Psi(a,b) = sqr2 P_i^{0,0}(a) P_j^{2i+1,0}(b) (1-b)**i  avec {a=2(1+r)/(1-s)-1 et b=s}
    !>
    !> r-derivative
    !> dPsi/dr = da/dr ∂a Psi + db/dr ∂b Psi  avec {da/dr=2/(1-s)=2(1-b) et db/dr=0}
    !>         = (2/(1-b)) ∂a Psi
    !>         = (2/(1-b)) sqr2 ∂a P_i^{0,0}(a) P_j^{2i+1,0}(b) (1-b)**i
    !>         = 2 sqr2 ∂a P_i^{0,0}(a) P_j^{2i+1,0}(b) (1-b)**(i-1)  (dans space)
    !>
    !> dPsi/dr = 2^{i+1/2) dP_i^{0,0}(a)/da P_j^{2i+1,0}(b) ((1-b)/2)**(i-1)  (Hesthaven identique space)
    !
    !
    !> s-derivative
    !> dPsi/ds = da/ds dPsi/da + db/ds dPsi/db  avec {da/ds= 2 (1+r)/(1-s)^2 et db/ds=1}
    !> 
    !> On remarque que si s/=1
    !> a=2(1+r)/(1-s)-1 et b=s => a=2(1+r)/(1-b)-1 => (1+r)=(1+a)(1-b)/2
    !> da/ds= 2 (1+r)/(1-s)^2 = 2 (1+a)(1-b)/2/(1-b)^2 = (1+a)/(1-b)
    
    !>dPsi/ds = (1+a)/(1-b) dPsi/da + dPsi/db
    !>dPsi/ds = (1+a)/(1-b) sqr2 ∂a[P_i^{0,0}(a)] P_j^{2i+1,0}(b) (1-b)**i
    !>         + sqr2 P_i^{0,0}(a) ∂b[ P_j^{2i+1,0}(b) (1-b)**i]
    !
    !>dPsi/ds =  sqr2 (1+a)/(1-b) ∂a[P_i^{0,0}(a)]   P_j^{2i+1,0}(b)    (1-b)**i
    !>         + sqr2               P_i^{0,0}(a)  ∂b[P_j^{2i+1,0}(b)]   (1-b)**i
    !>         + sqr2               P_i^{0,0}(a)     P_j^{2i+1,0}(b) ∂b[(1-b)**i]
    !
    !>dPsi/ds =  sqr2 ∂a[P_i^{0,0}(a)]   P_j^{2i+1,0}(b)  (1+a) (1-b)**(i-1)
    !>         - sqr2 i  P_i^{0,0}(a)    P_j^{2i+1,0}(b)        (1-b)**(i-1)
    !>         + sqr2    P_i^{0,0}(a) ∂b[P_j^{2i+1,0}(b)]       (1-b)**(i  )
    !
    !>dPsi/ds =[  ∂a[P_i^{0,0}(a)]   P_j^{2i+1,0}(b)  (1+a)
    !>          -    P_i^{0,0}(a)    P_j^{2i+1,0}(b)   i      ] sqr2 (1-b)**(i-1)
    !>         +     P_i^{0,0}(a) ∂b[P_j^{2i+1,0}(b) (1-b)**iu  sqr2
    !
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => drMode(1:np,1:n) && dsMode(1:np,1:n)
    !> Transpose = False => drMode(1:n,1:np) && dsMode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:),b(:)
    real(8), intent(out), pointer :: drMode(:,:),dsMode(:,:)
    logical, intent(in)           :: transpose
    !---
    integer                       :: i,ad
    integer                       :: nMod,nNod
    integer                       :: iu,iv,iw
    real(8), pointer              :: fa(:),dfa(:)
    real(8), pointer              :: fb(:),dfb(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nMod=(ord+1)*(ord+2)/2 ; nNod=size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.transpose )then
      
      allocate(drMode(1:nNod,1:nMod),dsMode(1:nNod,1:nMod))
      
      do iu=0,ord
        call  jacobiP(n=iu,alpha=0d0,beta=0d0,u=a(1:nNod),jf= fa)
        call djacobiP(n=iu,alpha=0d0,beta=0d0,u=a(1:nNod),jf=dfa)
        do iv=0,ord-iu
          ad=iu+iv*(ord+1)-(iv*(iv-1))/2 +1 ! Rangement façon space
          call  jacobiP(n=iv,alpha=2d0*real(iu,kind=8)+1d0,beta=0d0,u=b(1:nNod),jf= fb)
          call djacobiP(n=iv,alpha=2d0*real(iu,kind=8)+1d0,beta=0d0,u=b(1:nNod),jf=dfb)
          
          !> drMode
          if( iu==0 )then
            drMode(1:nNod,ad)=(2d0*sqr2)*dfa(1:nNod)* fb(1:nNod)
          else
            drMode(1:nNod,ad)=(2d0*sqr2)*dfa(1:nNod)* fb(1:nNod) *(1d0-b(1:nNod))**(iu-1)
          endif
          
          !> dsMode
          !dsMode(1:nNod,ad)= dfa(1:nNod)*fb(1:nNod)*(1d0+a(1:nNod))*(1d0-b(1:nNod))**(iu-1)
          
          if( iu==0 )then
            dsMode(1:nNod,ad)= dfa(1:nNod)* fb(1:nNod)*(1d0+a(1:nNod)) &
            &                 + fa(1:nNod)*dfb(1:nNod)
          elseif( iu==1 )then
            dsMode(1:nNod,ad)= dfa(1:nNod)* fb(1:nNod)*(1d0+a(1:nNod))                &
            &                 - fa(1:nNod)* fb(1:nNod)*real(iu,kind=8)                &
            &                 + fa(1:nNod)*dfb(1:nNod)                 *(1d0-b(1:nNod))
          else
            dsMode(1:nNod,ad)= dfa(1:nNod)* fb(1:nNod)*(1d0+a(1:nNod)) *(1d0-b(1:nNod))**(iu-1) &
            &                 - fa(1:nNod)* fb(1:nNod)*real(iu,kind=8) *(1d0-b(1:nNod))**(iu-1) &
            &                 + fa(1:nNod)*dfb(1:nNod)                 *(1d0-b(1:nNod))**iu
          endif
          
          
!            dsMode(1:nNod,ad)=( dsMode(1:nNod,ad)                  &
!            &               -fa(1:nNod)*fb(1:nNod)*real(iu,kind=8) &
!            &              )*(1d0-b(1:nNod))**(iu-1)
!          else
!            dsMode(1:nNod,ad)=( dsMode(1:nNod,ad)                  &
!            &               -fa(1:nNod)*fb(1:nNod)*real(iu,kind=8) &
!            &              )*(1d0-b(1:nNod))**(iu-1)
!          endif
!          
!          dsMode(1:nNod,ad)=dsMode(1:nNod,ad)+fa(1:nNod)*dfb(1:nNod)*(1d0-b(1:nNod))**iu
          dsMode(1:nNod,ad)=sqr2*dsMode(1:nNod,ad)
          
          deallocate(fb,dfb)
        enddo
        deallocate(fa,dfa)
      enddo
    else
      
      allocate(drMode(1:nMod,1:nNod),dsMode(1:nMod,1:nNod))
      stop "stop@gradSimplex2D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine gradSimplex2D


  subroutine vandermonde2D_c(order, a, b, vand)  BIND(C, name="SNB_vandermonde2D") 
    use, intrinsic :: ISO_C_BINDING 
    implicit none 

    integer (C_INT), value :: order ! int en C
    type (C_PTR),    value :: a     ! void * en C
    type (C_PTR),    value :: b     ! void * en C
    type (C_PTR)           :: vand  ! void ** en C

    real(8), pointer       :: a_f(:), b_f(:), vand_f(:,:)
    logical                :: display_f

    integer                :: ord_f
    integer                :: n

    n=(order+1)*(order+2)/2

    ord_f = order
    
    call c_f_pointer (a, a_f, (/n/))
    call c_f_pointer (b, b_f, (/n/))
    
    call vandermonde2D(ord_f, a_f, b_f, vand_f)

    vand = c_loc (vand_f)
    
  end subroutine vandermonde2D_c
  
  subroutine vandermonde2D(ord,a,b,vand)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>     [Phi_0(xi_1) ... Phi_{n-1}(xi_1)]
    !> V = [                               ] 
    !>     [Phi_0(xi_n) ... Phi_{n-1}(xi_n)]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:)
    real(8), intent(in) , pointer :: b(:)
    real(8), intent(out), pointer :: vand(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(a) .or. .not.associated(b) )then
      print '("2D Nodes a or/and b not associated")'
      stop "@vandermonde2D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> V(i,j) = Phi_j(xi_i) 
    !> vand(1:nNod,1:nMod) => Transpose=.false.
    call simplex2D(        &
    &    ord=ord          ,&
    &    a=a,b=b          ,&
    &    mode=vand        ,& !> vand(1:Nod,1:nMod) allocated in simplex2D
    &    transpose=.false. )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine vandermonde2D
  
  subroutine gradVandermonde2D(ord,a,b,drVand,dsVand)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>          [drPhi_0(xi_1) ... drPhi_{n-1}(xi_1)]
    !> drVand = [                                   ] 
    !>          [drPhi_0(xi_n) ... drPhi_{n-1}(xi_n)]
    !>
    !>          [dsPhi_0(xi_1) ... dsPhi_{n-1}(xi_1)]
    !> dsVand = [                                   ] 
    !>          [dsPhi_0(xi_n) ... dsPhi_{n-1}(xi_n)]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:),b(:)
    real(8), intent(out), pointer :: drVand(:,:),dsVand(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(a) .or. .not.associated(b) )then
      print '("2D Nodes a or/and b not associated")'
      stop "@gradVandermonde2D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> drVand(i,j) = drPhi_j(xi_i) => transpose=.false.
    !> dsVand(i,j) = dsPhi_j(xi_i) => transpose=.false. 
    call gradSimplex2D(    &
    &    ord=ord          ,&
    &    a=a,b=b          ,&
    &    drMode=drVand    ,& !> drVand(1:Nod,1:nMod) allocated in gradSimplex2D
    &    dsMode=dsVand    ,& !> drVand(1:Nod,1:nMod) allocated in gradSimplex2D
    &    transpose=.false. )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine gradVandermonde2D

  subroutine lagrange2Dv_1_c (order, nVtx, vand, a, b, lx, transpose)  BIND(C, name="SNB_lagrange2Dv") 
    use, intrinsic :: ISO_C_BINDING 
    implicit none 

    integer (C_INT), value :: order 
    integer (C_INT), value :: nVtx 
    type (C_PTR),    value :: vand
    type (C_PTR),    value :: a
    type (C_PTR),    value :: b
    type (C_PTR),    value :: lx
    integer (C_INT), value :: transpose

    integer            :: order_f
    integer            :: nVtx_f
    real(8), pointer   :: vand_f(:,:)
    real(8), pointer   :: lx_f(:,:)
    real(8), pointer   :: a_f(:),b_f(:)
    logical            :: transpose_f

    integer          :: n

    n=(order+1)*(order+2)/2

    order_f = order
    nVtx_f = nVtx

    transpose_f = transpose
    
    call c_f_pointer (vand, vand_f, (/n,n/))
    call c_f_pointer (a, a_f, (/nVtx_f/))
    call c_f_pointer (b, b_f, (/nVtx_f/))
    
    if (transpose_f) then
      call c_f_pointer (lx, lx_f, (/n, nVtx_f/))
    else
      call c_f_pointer (lx, lx_f, (/nVtx_f, n/))
    endif

    call lagrange2Dv_1(order_f, vand_f, a_f, b_f, lx_f, transpose_f)
    
  end subroutine lagrange2Dv_1_c

  
  subroutine lagrange2Dv_1(ord,vand,a,b,lx,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define lagrange2Dv_1 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !! procedure moins directe mais plus performante
    !! car a,b et vand ne sont pas recalculés à chaque appel
    ! calcule les nMod fonctions de base pour l'ensemble des nNod points (a,b)
    ! lagrange2D := Inverse[Transpose[Vand]].Psi[x]
    ! transpose = .true.  => lx(1:nMod,1:nNod)
    ! transpose = .false. => lx(1:nNod,1:nMod)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)              :: ord
    real(8), intent(in)    , pointer :: vand(:,:)
    real(8), intent(in)    , pointer :: a(:),b(:)
    real(8), intent(inout) , pointer :: lx  (:,:)
    logical, intent(in)              :: transpose
    !---
    integer                          :: nMod,nNod
    integer                          :: iu,iv,ad,i,j,k
    real(8)                          :: gamma(0:ord+1)
    integer                          :: iOrd
    real(8), pointer                 :: psi(:,:)
    real(8), pointer                 :: mat(:,:)
    integer                          :: lWork
    integer, pointer                 :: ipiv(:)
    real(8), pointer                 :: work(:)
    integer                          :: iErr
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if lagrange2Dv_1==1
    print '(">>> baseSimplex2D:lagrange2Dv_1")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nMod=(ord+1)*(ord+2)/2 ; nNod=size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Computing psi(a,b) 
    !> transpose=.false. pour les perfos
    call simplex2D(ord=ord,a=a,b=b,mode=psi,transpose=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if lagrange2Dv_1==1
    print '("    baseSimplex2D:lagrange2Dv_1: nMod=",i10," nNod=",i10)',nMod,nNod
    print '("    baseSimplex2D:lagrange2Dv_1: size(psi) =",i10," x ",i10)',size(psi,1),size(psi,2)
    print '("    baseSimplex2D:lagrange2Dv_1: size(vand)=",i10," x ",i10)',size(vand,1),size(vand,2)
    print '("    baseSimplex2D:lagrange2Dv_1: step1 OK")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> mat=Inverse[Transpose[Vand]]
    allocate(mat(1:nMod,1:nMod))
    do i=1,nMod
      do j=1,nMod
        mat(i,j)=vand(j,i)
      enddo
    enddo
    
    lWork=64*(nMod) ; allocate(work(1:lWork),ipiv(1:nMod))
    call dgetrf(nMod,nMod,mat(1,1),nMod,ipiv(1),iErr)
    call dgetri(nMod,mat(1,1),nMod,ipiv(1),work(1),lWork,iErr)
    deallocate(ipiv,work)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if lagrange2Dv_1==1
    print '(">>> baseSimplex2D:lagrange2Dv_1: step2 OK")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> lx = Inverse[Transpose[Vand]].psi = mat.psi
    if( transpose )then
     !allocate(lx(1:nMod,1:nNod)) ; lx(:,:)=0d0
      lx(:,:)=0d0
      do i=1,nNod
        do j=1,nMod ; do k=1,nMod
          lx(j,i)=lx(j,i)+mat(j,k)*psi(i,k)  !> Attention de bien prendre psi(i,k)
        enddo ; enddo
      enddo
    else
     !allocate(lx(1:nNod,1:nMod)) ; lx(:,:)=0d0
      lx(:,:)=0d0
      do i=1,nNod
        do j=1,nMod ; do k=1,nMod
          lx(i,j)=lx(i,j)+mat(j,k)*psi(i,k)  !> Attention de bien prendre psi(i,k)
        enddo ; enddo
      enddo
    endif
   !call displayMatrix(title="lx",mat=lx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if lagrange2Dv_1==1
    print '(">>> baseSimplex2D:lagrange2Dv_1: step3 OK")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(mat)
    deallocate(psi)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if lagrange2Dv_1==1
    print '("<<< baseSimplex2D:lagrange2Dv_1")'
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef lagrange2Dv_1
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine lagrange2Dv_1
  
  
  subroutine lagrange2Dv_2(ord,uvwOut,lagrange,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !! procedure plus directe mais moins performante
    !! car uv et vand doivent être recalculé à chaque appel
    ! calcule les nMod fonctions de base pour l'ensemble des nNod points (a,b)
    ! lagrange2D := Inverse[Transpose[Vand]].Psi[x]
    ! transpose = .true.  => lx(1:nMod,1:nNod)
    ! transpose = .false. => lx(1:nNod,1:nMod)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer         , intent(in)  :: ord
    real(8)         , intent(in)  :: uvwOut(:,:)
    real(8), pointer, intent(out) :: lagrange(:,:)
    logical         , intent(in)  :: transpose
    !>
    integer                       :: nMod,nNod
    real(8), pointer              :: uvw(:,:),a(:),b(:),c(:)
    real(8), pointer              :: vand(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nMod=(ord+1)*(ord+2)/2
    nNod=size(uvwOut,2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( transpose )then
      allocate(lagrange(1:nMod,1:nNod))
    else
      allocate(lagrange(1:nNod,1:nMod))
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( ord==1 )then
      if( transpose )then
        lagrange(1,1:nNod)=1d0-uvwOut(1,1:nNod)-uvwOut(2,1:nNod)
        lagrange(2,1:nNod)=    uvwOut(1,1:nNod)
        lagrange(3,1:nNod)=                     uvwOut(2,1:nNod)
      else
        lagrange(1:nNod,1)=1d0-uvwOut(1,1:nNod)-uvwOut(2,1:nNod)
        lagrange(1:nNod,2)=    uvwOut(1,1:nNod)
        lagrange(1:nNod,3)=                     uvwOut(2,1:nNod)
      endif
    else
      !> Calcul de Vand(:,:)
      call nodes2D    (ord=ord,uvw=uvw,display=.false.)
      call nodes2Dopt (ord=ord,uvw=uvw,display=.false.)
      call nodes2Duv2ab(uv=uvw,a=a,b=b,display=.false.)
      !> calcul de la matrice de Vandermonde
      call vandermonde2D(ord=ord,a=a,b=b,vand=vand)
      !> Calcul des polonômes de Lagrange d'ordre ord en uvwOut
      call nodes2Duv2ab(uv=uvwOut,a=a,b=b,display=.false.)
      call lagrange2Dv_1(ord=ord,vand=vand,a=a,b=b,lx=lagrange,transpose=.true.)  !> lagrange= Inverse[Transpose[Vand]].Psi[xyzOut] lxOut(nPt,np)
      !> Nettoyage memoire
      deallocate(uvw,a,b,vand)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    return
  end subroutine lagrange2Dv_2
  
  
  
  subroutine readXYout2D(xyout)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(out), pointer :: xyout (:,:)
    !-
    integer :: i,j,nVert
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    open(unit=10,file='Mesh2D/Triangles.mesh',status='old',action='read')
    do i=1,8 ; read(10,*) ; enddo
    read(10,*)nVert
    allocate(xyout(2,nVert))
    do i=1,nVert
      read(10,*)xyout(1:2,i)
    enddo
    close(10)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine readXYout2D
  
  subroutine writeSolOut2D(title,solOut)
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
      
      call system("ln -fs Mesh2D/Triangles.mesh "// title // "_" //sfx//".mesh")
      
      open(unit=iFile,file=title//"_"//sfx//".sol",status='unknown',action='write')
      write(iFile,'("MeshVersionFormatted 2"/)')
      write(iFile,'("Dimension 2"/)')
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
  end subroutine writeSolOut2D
  
  subroutine edgesConnectivity(ord,conec)!,edges,uvw)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer,              intent(in)  :: ord
    integer, allocatable, intent(out) :: conec(:,:)
    !>
    integer               :: iu,iv,iw,ad1,ad2,ad3,dg
    integer               :: i,j
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(conec(ord+1,3))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> sideI (iw=0)
    ad1=0 ; ad2=ord+2 ; ad3=0
    dg=0
    iw=0
    do iv=0,ord
      do iu=0,ord
        do iw=0,ord
          if( iu+iv+iw==ord )then
            dg=dg+1
            !> SideI
            if( iw==0 )then
              ad1=ad1+1 ; conec(ad1,1)=dg
            endif
            !> SideII
            if( iu==0 )then
              ad2=ad2-1 ; conec(ad2,2)=dg
            endif
            !> SideIII
            if( iv==0 )then
              ad3=ad3+1 ; conec(ad3,3)=dg
            endif
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef debug
    if( ord<10 )then
      write(*,'(/"Triangle P",i1)')ord
    else
      write(*,'(/"Triangle P",i2)')ord
    endif
    write(*,'(/"Local DOF / side")')
    do j=1,3
      write(*,'(/"sd",i1,": ",$)')j
      do i=1,size(conec,1)
        write(*,'(i4," ",$)')conec(i,j)
      enddo
      write(*,'()')
    enddo
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine edgesConnectivity
  
  
  subroutine trianglePointsPlot(file_name, node_xy, node_show, point_num, point_xy, point_show)
    !*****************************************************************************80
    !
    !! TRIANGLE_POINTS_PLOT plots a triangle and some points.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 October 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) FILE_NAME, the name of the output file.
    !
    !    Input, real(8) NODE_XY(2,3), the coordinates of the nodes
    !    of the triangle.
    !
    !    Input, integer NODE_SHOW,
    !   -1, do not show the triangle, or the nodes.
    !    0, show the triangle, do not show the nodes;
    !    1, show the triangle and the nodes;
    !    2, show the triangle, the nodes and number them.
    !
    !    Input, integer :: POINT_NUM, the number of points.
    !
    !    Input, real(8) POINT_XY(2,POINT_NUM), the coordinates of the
    !    points.
    !
    !    Input, integer :: POINT_SHOW,
    !    0, do not show the points;
    !    1, show the points;
    !    2, show the points and number them.
    !
    character(80), intent(in) :: file_name
    real(8)      , intent(in) :: node_xy(1:2,1:3)
    integer      , intent(in) :: node_show
    integer      , intent(in) :: point_num
    real(8)      , intent(in) :: point_xy(1:3,1:point_num)
    integer      , intent(in) :: point_show
    
    integer, parameter :: node_num = 3
    
    integer            :: circle_size
    integer            :: delta
    integer            :: e
    integer            :: file_unit
    integer            :: i
    integer            :: ios
    integer            :: node
    integer            :: point
    character(40)      :: string
    real(8)            :: x_min,x_max
    integer            :: x_ps
    integer            :: x_ps_max = 576
    integer            :: x_ps_max_clip = 594
    integer            :: x_ps_min = 36
    integer            :: x_ps_min_clip = 18
    real(8)            :: x_scale
    real(8)            :: y_min,y_max
    integer            :: y_ps
    integer            :: y_ps_max = 666
    integer            :: y_ps_max_clip = 684
    integer            :: y_ps_min = 126
    integer            :: y_ps_min_clip = 108
    real (8)           ::  y_scale
  !
  !  We need to do some figuring here, so that we can determine
  !  the range of the data, and hence the height and width
  !  of the piece of paper.
  !
    x_max = max ( maxval (  node_xy(1,1:node_num) ), &
                  maxval ( point_xy(1,1:point_num) ) )
    x_min = min ( minval (  node_xy(1,1:node_num) ), &
                  minval ( point_xy(1,1:point_num) ) )
    x_scale = x_max - x_min
  
    x_max = x_max + 0.05D+00 * x_scale
    x_min = x_min - 0.05D+00 * x_scale
    x_scale = x_max - x_min
  
    y_max = max ( maxval (  node_xy(2,1:node_num) ), &
                  maxval ( point_xy(2,1:point_num) ) )
    y_min = min ( minval (  node_xy(2,1:node_num) ), &
                  minval ( point_xy(2,1:point_num) ) )
    y_scale = y_max - y_min
  
    y_max = y_max + 0.05D+00 * y_scale
    y_min = y_min - 0.05D+00 * y_scale
    y_scale = y_max - y_min
  
    if ( x_scale < y_scale ) then
  
      delta = nint ( real ( x_ps_max - x_ps_min, kind = 8 ) &
        * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )
  
      x_ps_max = x_ps_max - delta
      x_ps_min = x_ps_min + delta
  
      x_ps_max_clip = x_ps_max_clip - delta
      x_ps_min_clip = x_ps_min_clip + delta
  
      x_scale = y_scale
  
    else if ( y_scale < x_scale ) then
  
      delta = nint ( real ( y_ps_max - y_ps_min, kind = 8 ) &
        * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )
  
      y_ps_max      = y_ps_max - delta
      y_ps_min      = y_ps_min + delta
  
      y_ps_max_clip = y_ps_max_clip - delta
      y_ps_min_clip = y_ps_min_clip + delta
  
      y_scale = x_scale
  
    end if
  
    call get_unit ( file_unit )
  
    open ( unit = file_unit, file = file_name, status = 'replace',iostat = ios )
    
    if( .not.ios== 0 )then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGLE_POINTS_PLOT - Fatal error!'
      write ( *, '(a)' ) '  Can not open output file.'
      return
    end if
    
    write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
    write ( file_unit, '(a)' ) '%%Creator: triangulation_order3_plot.f90'
    write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
    write ( file_unit, '(a)' ) '%%Pages: 1'
    write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ', &
      x_ps_min, y_ps_min, x_ps_max, y_ps_max
    write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
    write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
    write ( file_unit, '(a)' ) '%%EndComments'
    write ( file_unit, '(a)' ) '%%BeginProlog'
    write ( file_unit, '(a)' ) '/inch {72 mul} def'
    write ( file_unit, '(a)' ) '%%EndProlog'
    write ( file_unit, '(a)' ) '%%Page: 1 1'
    write ( file_unit, '(a)' ) 'save'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw a gray border around the page.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) 'newpath'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
    write ( file_unit, '(a)' ) 'stroke'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to black.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the font and its size.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.50 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Print a title.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  210  702  moveto'
    write ( file_unit, '(a)' ) '%  (Triangulation)  show'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) 'newpath'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
      x_ps_min_clip, y_ps_min_clip, ' moveto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
      x_ps_max_clip, y_ps_min_clip, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
      x_ps_max_clip, y_ps_max_clip, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
      x_ps_min_clip, y_ps_max_clip, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
      x_ps_min_clip, y_ps_min_clip, ' lineto'
    write ( file_unit, '(a)' ) 'clip newpath'
  !
  !  Draw the nodes.
  !
    if ( 1 <= node_show ) then
  
      circle_size = 5
  
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
      write ( file_unit, '(a)' ) '%'
  
      do node = 1, 3
  
        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
          + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                   - x_min ) )
  
        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
          + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                   - y_min ) )
        
        write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
          circle_size, '0 360 arc closepath fill'
        
      end do
      
    end if
  !
  !  Label the nodes.
  !
    if ( 2 <= node_show ) then
  
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Label the nodes:'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the RGB color to darker blue.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
      write ( file_unit, '(a)' ) '/Times-Roman findfont'
      write ( file_unit, '(a)' ) '0.20 inch scalefont'
      write ( file_unit, '(a)' ) 'setfont'
      write ( file_unit, '(a)' ) '%'
  
      do node = 1, node_num
  
        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
          + (       + node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                   - x_min ) )
  
        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
          + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                   - y_min ) )
  
        write ( string, '(i4)' ) node
        string = adjustl ( string )
  
        write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
          ' moveto (' // trim ( string ) // ') show'
  
      end do
  
    end if
  !
  !  Draw the points.
  !
    if ( point_num <= 200 ) then
      circle_size = 5
    else if ( point_num <= 500 ) then
      circle_size = 4
    else if ( point_num <= 1000 ) then
      circle_size = 3
    else if ( point_num <= 5000 ) then
      circle_size = 2
    else
      circle_size = 1
    end if
  
    if ( 1 <= point_show ) then
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Draw filled dots at the points.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the RGB color to green.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.150  0.750  0.000 setrgbcolor'
      write ( file_unit, '(a)' ) '%'
  
      do point = 1, point_num
  
        x_ps = int ( &
          ( ( x_max - point_xy(1,point)         ) &
          * real ( x_ps_min, kind = 8 )   &
          + (         point_xy(1,point) - x_min ) &
          * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                     - x_min ) )
  
        y_ps = int ( &
          ( ( y_max - point_xy(2,point)         ) &
          * real ( y_ps_min, kind = 8 )   &
          + (         point_xy(2,point) - y_min ) &
          * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                     - y_min ) )
  
        write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
          circle_size, '0 360 arc closepath fill'
  
      end do
  
    end if
  !
  !  Label the points.
  !
    if ( 2 <= point_show ) then
  
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Label the point:'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the RGB color to darker green.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.250  0.850  0.000 setrgbcolor'
      write ( file_unit, '(a)' ) '/Times-Roman findfont'
      write ( file_unit, '(a)' ) '0.20 inch scalefont'
      write ( file_unit, '(a)' ) 'setfont'
      write ( file_unit, '(a)' ) '%'
  
      do point = 1, point_num
  
        x_ps = int ( &
          ( ( x_max - point_xy(1,point)         ) &
          * real ( x_ps_min, kind = 8 )   &
          + (       + point_xy(1,point) - x_min ) &
          * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                     - x_min ) )
  
        y_ps = int ( &
          ( ( y_max - point_xy(2,point)         ) &
          * real ( y_ps_min, kind = 8 )   &
          + (         point_xy(2,point) - y_min ) &
          * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                     - y_min ) )
  
        write ( string, '(i4)' ) point
        string = adjustl ( string )
  
        write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
          ' moveto (' // trim ( string ) // ') show'
  
      end do
  
    end if
  !
  !  Draw the triangle.
  !
    if ( 0 <= node_show ) then
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Draw the triangle.'
      write ( file_unit, '(a)' ) '%'
  
      write ( file_unit, '(a)' ) 'newpath'
  
      do i = 1, 4
  
        node = i4_wrap( i, 1, 3 )
  
        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
          + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                   - x_min ) )
  
        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
          + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                   - y_min ) )
  
        if ( i == 1 ) then
          write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' moveto'
        else
          write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'
        end if
  
      end do
  
      write ( file_unit, '(a)' ) 'stroke'
  
    end if
  
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) 'restore  showpage'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  End of page.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%%Trailer'
    write ( file_unit, '(a)' ) '%%EOF'
    close ( unit = file_unit )
  
    return
  end subroutine trianglePointsPlot
  
  subroutine get_unit ( iunit )
  
  !*****************************************************************************80
  !
  !! GET_UNIT returns a free FORTRAN unit number.
  !
  !  Discussion:
  !
  !    A "free" FORTRAN unit number is an integer between 1 and 99 which
  !    is not currently associated with an I/O device.  A free FORTRAN unit
  !    number is needed in order to open a file with the OPEN command.
  !
  !    If IUNIT = 0, then no free FORTRAN unit could be found, although
  !    all 99 units were checked (except for units 5, 6 and 9, which
  !    are commonly reserved for console I/O).
  !
  !    Otherwise, IUNIT is an integer between 1 and 99, representing a
  !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
  !    are special, and will never return those values.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 September 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, integer :: IUNIT, the free unit number.
  !
    integer :: i
    integer :: ios
    integer :: iunit
    logical :: lopen
  
    iunit = 0
  
    do i = 1, 99
  
      if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
  
        inquire ( unit = i, opened = lopen, iostat = ios )
  
        if ( ios == 0 ) then
          if ( .not. lopen ) then
            iunit = i
            return
          end if
        end if
  
      end if
  
    end do
  
    return
  end subroutine get_unit
  
  
  integer function i4_wrap ( ival, ilo, ihi )
  
  !*****************************************************************************80
  !
  !! I4_WRAP forces an I4 to lie between given limits by wrapping.
  !
  !  Example:
  !
  !    ILO = 4, IHI = 8
  !
  !    I  Value
  !
  !    -2     8
  !    -1     4
  !     0     5
  !     1     6
  !     2     7
  !     3     8
  !     4     4
  !     5     5
  !     6     6
  !     7     7
  !     8     8
  !     9     4
  !    10     5
  !    11     6
  !    12     7
  !    13     8
  !    14     4
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 August 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer :: IVAL, an integer value.
  !
  !    Input, integer :: ILO, IHI, the desired bounds for the value.
  !
  !    Output, integer :: I4_WRAP, a "wrapped" version of IVAL.
  !
    
    integer :: ihi
    integer :: ilo
    integer :: ival
    integer :: jhi
    integer :: jlo
    integer :: value
    integer :: wide
  
    jlo = min ( ilo, ihi )
    jhi = max ( ilo, ihi )
  
    wide = jhi - jlo + 1
  
    if ( wide == 1 ) then
      value = jlo
    else
      value = jlo + i4_modp( ival - jlo, wide )
    end if
  
    i4_wrap = value
  
    return
  end function i4_wrap
  
  integer function i4_modp ( i, j )
  
  !*****************************************************************************80
  !
  !! I4_MODP returns the nonnegative remainder of I4 division.
  !
  !  Discussion:
  !
  !    If
  !      NREM = I4_MODP ( I, J )
  !      NMULT = ( I - NREM ) / J
  !    then
  !      I = J * NMULT + NREM
  !    where NREM is always nonnegative.
  !
  !    The MOD function computes a result with the same sign as the
  !    quantity being divided.  Thus, suppose you had an angle A,
  !    and you wanted to ensure that it was between 0 and 360.
  !    Then mod(A,360) would do, if A was positive, but if A
  !    was negative, your result would be between -360 and 0.
  !
  !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
  !
  !  Example:
  !
  !        I     J     MOD I4_MODP    Factorization
  !
  !      107    50       7       7    107 =  2 *  50 + 7
  !      107   -50       7       7    107 = -2 * -50 + 7
  !     -107    50      -7      43   -107 = -3 *  50 + 43
  !     -107   -50      -7      43   -107 =  3 * -50 + 43
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 March 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer :: I, the number to be divided.
  !
  !    Input, integer :: J, the number that divides I.
  !
  !    Output, integer :: I4_MODP, the nonnegative remainder when I is
  !    divided by J.
  !
    implicit none
  
    integer :: i
    integer :: j
    integer :: value
  
    if ( j == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4_MODP - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
      stop
    end if
  
    value = mod ( i, j )
  
    if ( value < 0 ) then
      value = value + abs ( j )
    end if
  
    i4_modp = value
  
    return
  end function i4_modp
  
  
  
end module baseSimplex2D
