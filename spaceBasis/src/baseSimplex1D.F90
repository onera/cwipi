module baseSimplex1D
  use baseSimplexTools
  implicit none
  contains
  
  subroutine setL2MeshIJK(meshOrder,ijk)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)    :: meshOrder
    integer, intent(inout) :: ijk(:,:)
    !>
    integer                :: iMod,nMod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> 01 03 04 05 02
    nMod=(meshOrder+1)
    
    ijk(1,1)=0
    ijk(1,2)=meshOrder
    do iMod=3,nMod
      ijk(1,iMod)=iMod-2
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    
    return
  end subroutine setL2MeshIJK
  
  subroutine setQ4MeshIJK(meshOrder,ij)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)    :: meshOrder
    integer, intent(inout) :: ij(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Quad (2D/3D)
    !>
    !>   04 03
    !>   01 02
    !> 2D
    !>   f1 : 01 02
    !>   f2 : 02 03
    !>   f3 : 03 04
    !>   f4 : 04 01
    !> 3D
    !>   f1 : 01 02 03 04
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Quad2 (2D/3D)
    !>
    !>   04 07 03
    !>   08 09 06
    !>   01 05 02
    !> 2D
    !>   f1 : 01 02 05
    !>   f2 : 02 03 06
    !>   f3 : 03 04 07
    !>   f4 : 04 01 08
    !> 3D
    !>   f1 : 01 02 03 04 05 06 07 08 09
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Quad3 (2D)
    !>
    !>   04 10 09 03
    !>   11 16 15 08
    !>   12 13 14 07
    !>   01 05 06 02
    !> 2D
    !>   f1 : 01 02 05 06
    !>   f2 : 02 03 07 08
    !>   f3 : 03 04 09 10
    !>   f4 : 04 01 11 12
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Quad4 (2D)
    !>
    !>   04 13 12 11 03
    !>   14 20 23 19 10
    !>   15 24 25 22 09
    !>   16 17 21 18 08
    !>   01 05 06 07 02
    !> 2D
    !>   f1 : 01 02 05 06 07
    !>   f2 : 02 03 08 09 10
    !>   f3 : 03 04 11 12 13
    !>   f4 : 04 01 14 15 16
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
    
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if    ( meshOrder==1 )then !> QuadQ1
      
      !>   04 03
      !>   01 02
      
      ij(1:2,01)=[0,0] !> 1
      ij(1:2,02)=[1,0] !> 2
      ij(1:2,03)=[1,1] !> 3
      ij(1:2,04)=[0,1] !> 4
      
    elseif( meshOrder==2 )then !> QuadQ2
      
      !>   04 07 03
      !>   08 09 06
      !>   01 05 02
      
      ij(1:2,01)=[0,0] !> 1
      ij(1:2,02)=[2,0] !> 2
      ij(1:2,03)=[2,2] !> 3
      ij(1:2,04)=[0,2] !> 4
      ij(1:2,05)=[1,0] !> 5
      ij(1:2,06)=[2,1] !> 6
      ij(1:2,07)=[1,2] !> 7
      ij(1:2,08)=[0,1] !> 8
      ij(1:2,09)=[1,1] !> 9
      
    elseif(meshOrder==3 )then !> QuadQ3
      
      !>   04 10 09 03
      !>   11 16 15 08
      !>   12 13 14 07
      !>   01 05 06 02
      
      ij(1:2,01)=[0,0] !> 01
      ij(1:2,02)=[3,0] !> 02
      ij(1:2,03)=[3,3] !> 03
      ij(1:2,04)=[0,3] !> 04
      
      ij(1:2,05)=[1,0] !> 05
      ij(1:2,06)=[2,0] !> 06
      ij(1:2,07)=[3,1] !> 07
      ij(1:2,08)=[3,2] !> 08
      ij(1:2,09)=[2,3] !> 09
      ij(1:2,10)=[1,3] !> 10
      ij(1:2,11)=[0,2] !> 11
      ij(1:2,12)=[0,1] !> 12
      ij(1:2,13)=[1,1] !> 13
      ij(1:2,14)=[2,1] !> 14
      ij(1:2,15)=[2,2] !> 15
      ij(1:2,16)=[1,2] !> 16
      
    else ; stop "setQ4MeshIJK meshOrder>3 not implemented"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    
    return
  end subroutine setQ4MeshIJK
  
  subroutine setL2BasisEqui(ord,ijk,u,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Calcul des base d'ordere ord pour des triangles dont les points
    ! d'interpolation sont equidistants.
    ! Numerotation des sommets suivant ijk
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !use baseSimplexTools, only: monomialProduct
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)    :: ord 
    integer, intent(in)    :: ijk(:,:)     !> ijk(1:1,1:nMod)
    real(8), intent(in)    :: u(:)         !> u(1:nNod)        \in [-1,1]
    real(8), intent(inout) :: ai(:,:)      !> ai(1:nMod,1:nNod)
    !>
    integer                :: iMod,nMod    !> nMod=(ord+1)*(ord+2)/2
    integer                :: iNod,nNod    !> nNod=size(u) 
    integer                :: i,j
    real(8), pointer       :: xi(:)
    real(8) :: var
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call nodes1D(ord=ord, uvw=xi, display=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> bases de Lagrange
    !> ai(1:nMod,1:nNod)=1d0
    nMod=ord+1
    nNod=size(u)
    ai(:,:)=1d0
    do i=1,nMod
      do j=1,nMod
        if( .not.i==j )then
          var=(xi(i)-xi(j))
          var=1d0/var
          
          iMod=ijk(1,i)+1
          do iNod=1,nNod
            ai(iMod,iNod)= ai(iMod,iNod)    &
            &              *(u(iNod)-xi(j)) &
            &              *var
          enddo
        endif
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(xi)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine setL2BasisEqui
  
  subroutine setQ4BasisEqui_u(ord,ijk,u,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Calcul des base d'ordere ord pour des triangles dont les points
    ! d'interpolation sont equidistants.
    ! Numerotation des sommets suivant ijk
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)    :: ord 
    integer, intent(in)    :: ijk(:,:)     !> ijk(1:2,1:(ord+1)*(ord+1))
    real(8), intent(in)    :: u(:)         !> u(1:nNod)                 \in [-1,1]
    real(8), intent(inout) :: ai(:,:)      !> ai(1:nMod,1:nNod)
    !>
    integer                :: iMod,nMod    !> nMod=(ord+1)
    integer                :: iNod,nNod    !> nNod=size(u) 
    integer, pointer       :: ij(:,:)      !> ij(1:1,1:ord+1)
    integer                :: i,j,k,l
    integer                :: iu,iv
    real(8), pointer       :: xi(:)
    real(8), pointer       :: lagrangeL2    (:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nMod=(ord+1)                                       !> Edge meshOrder (entree) = ord
    nNod=size(u)                                       !> Edge compOrder (sortie)
    
    allocate(ij(1:1,1:nMod)) 
    do iMod=1,ord+1                                    !> 01 02  ... meshOrder meshOrder+1 (rangement regulier)
      ij(1,iMod)=(iMod-1)
    enddo
    
    allocate(lagrangeL2(1:nMod,1:nNod))
    call setL2BasisEqui(ord=ord,ijk=ij,u=u,ai=lagrangeL2)
    
    !block
    !real(8) :: toto(11)
    !print '(/"setQ4BasisEqui")'
    ! print '(/"u=",*(e12.5,1x))',u(:)
    !print '(/"lagrangeL2")'
    !do iMod=1,nMod
    !  toto(1:nMod)=lagrangeL2(iMod,1:nNod)
    !  print '(3x,"ai(",i2,",1:nNod)=",*(e12.5,1x) )',iMod,toto(1:nMod)
    !enddo
    !print '("")'
    !end block
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iMod=0
    do j=1,nMod ; do i=1,nMod
      iMod=iMod+1 ! print '(3x,"iMod=",i3,3x,"i,j=",2(i3,1x))',iMod,ij(1:2,iMod)
      iu=ijk(1,iMod)+1
      iv=ijk(2,iMod)+1
      !>
      !print '(3x,"iMod=",i3,3x,"(i,j)=",2(i3,1x))',iMod,iu,iv          
      iNod=0
      do l=1,nNod ; do k=1,nNod
        iNod=iNod+1
        !print '(6x,"iNod=",i3,3x,"(u,v)=",2(e12.5,1x),3x,"aiL2=",2(e12.5,1x))',iNod,u(k),u(l),lagrangeL2(iu,k),lagrangeL2(iv,l)
        
        ai(iMod,iNod)= lagrangeL2(iu,k)*lagrangeL2(iv,l)
      enddo ; enddo
    enddo ; enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(ij)
    deallocate(lagrangeL2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          
    return
  end subroutine setQ4BasisEqui_u


  subroutine setQ4BasisEqui_uv_c(ord, n_vtx, ijk, u, v, ai)  BIND(C, name="SNB_setQ4BasisEqui_uv") 
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use, intrinsic :: ISO_C_BINDING 
    implicit none 
    integer (C_INT), value :: ord 
    integer (C_INT), value :: n_vtx 
    type (C_PTR),    value  :: ijk
    type (C_PTR),    value  :: u
    type (C_PTR),    value  :: v
    type (C_PTR),    value  :: ai
    !>
    integer, pointer :: ijk_f(:,:)
    real(8), pointer :: u_f(:)
    real(8), pointer :: v_f(:)
    real(8), pointer :: ai_f(:,:)
    integer          :: ord_f, n_vtx_f, n_mode
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ord_f   = ord
    n_vtx_f = n_vtx
    n_mode = (ord_f+1)*(ord_f+1)
    
    call c_f_pointer (ijk, ijk_f, [2     , n_mode ])
    call c_f_pointer (u, u_f, [n_vtx_f])
    call c_f_pointer (v, v_f, [n_vtx_f])
    call c_f_pointer (ai , ai_f , [n_mode, n_vtx_f])
    
    call setQ4BasisEqui_uv(ord_f, ijk_f, u_f, v_f, ai_f)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine setQ4BasisEqui_uv_c

  
  subroutine setQ4BasisEqui_uv(ord,ijk,u,v,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Calcul des base d'ordere ord pour des triangles dont les points
    ! d'interpolation sont equidistants.
    ! Numerotation des sommets suivant ijk
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)    :: ord 
    integer, intent(in)    :: ijk(:,:)     !> ijk(1:2,1:(ord+1)*(ord+1))
    real(8), intent(in)    :: u(:)         !> u(1:nNod)                 \in [-1,1]
    real(8), intent(in)    :: v(:)         !> v(1:nNod)                 \in [-1,1]
    real(8), intent(inout) :: ai(:,:)      !> ai(1:nMod,1:nNod)
    !>
    integer                :: iMod,nMod    !> nMod=(ord+1)
    integer                :: iNod,nNod    !> nNod=size(u) 
    integer, pointer       :: ij(:,:)      !> ij(1:1,1:ord+1)
    integer                :: i,j,k,l
    integer                :: iu,iv
    real(8), pointer       :: xi(:)
    real(8), pointer       :: lagrangeL2_u(:,:)
    real(8), pointer       :: lagrangeL2_v(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nMod=(ord+1)                                       !> Edge meshOrder (entree) = ord
    nNod=size(u) !=size(v)                             !> Edge compOrder (sortie)
    
    allocate(ij(1:1,1:nMod)) 
    do iMod=1,ord+1                                    !> 01 02  ... meshOrder meshOrder+1 (rangement regulier)
      ij(1,iMod)=(iMod-1)
    enddo
    
    allocate(lagrangeL2_u(1:nMod,1:nNod)) ; call setL2BasisEqui(ord=ord,ijk=ij,u=u,ai=lagrangeL2_u)
    allocate(lagrangeL2_v(1:nMod,1:nNod)) ; call setL2BasisEqui(ord=ord,ijk=ij,u=v,ai=lagrangeL2_v)    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iMod=0
    do j=1,nMod ; do i=1,nMod
      iMod=iMod+1 ! print '(3x,"iMod=",i3,3x,"i,j=",2(i3,1x))',iMod,ij(1:2,iMod)
      iu=ijk(1,iMod)+1
      iv=ijk(2,iMod)+1
      !>
      !print '(3x,"iMod=",i3,3x,"(i,j)=",2(i3,1x))',iMod,iu,iv          
      do iNod=1,nNod
        !print '(6x,"iNod=",i3,3x,"(u,v)=",2(e12.5,1x),3x,"aiL2=",2(e12.5,1x))',iNod,u(iNod),v(iNod),lagrangeL2_u(iu,iNod),lagrangeL2_v(iv,iNod)        
        ai(iMod,iNod)= lagrangeL2_u(iu,iNod)*lagrangeL2_v(iv,iNod)
      enddo
    enddo ; enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(ij)
    deallocate(lagrangeL2_u)
    deallocate(lagrangeL2_v)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          
    return
  end subroutine setQ4BasisEqui_uv
  
  subroutine nodes1D_c(ord, uvw, display)  BIND(C, name="SNB_nodes1D") 
    use, intrinsic :: ISO_C_BINDING 
    implicit none 
    integer (C_INT), value :: ord 
    type (C_PTR)           :: uvw
    integer (C_INT), value :: display

    real(8), pointer :: uvw_f(:)
    integer          :: ord_f
    logical          :: display_f

    ord_f = ord
    display_f = (display == 1)

    call nodes1D (ord_f, uvw_f, display_f)

    uvw = c_loc (uvw_f)

  end subroutine nodes1D_c
  
  
  subroutine nodes1D(ord, uvw, display)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity edge
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
  
  subroutine vandermonde1D_c(order, a, vand)  BIND(C, name="SNB_vandermonde1D") 
    use, intrinsic :: ISO_C_BINDING 
    implicit none 

    integer (C_INT), value :: order ! int en C
    type (C_PTR),    value :: a     ! void * en C
    type (C_PTR)           :: vand  ! void ** en C

    real(8), pointer       :: a_f(:), vand_f(:,:)
    logical                :: display_f

    integer                :: ord_f
    integer                :: n

    n=(order+1)

    ord_f = order
    
    call c_f_pointer (a, a_f, (/n/))
    
    call vandermonde1D(ord_f, a_f, vand_f)

    vand = c_loc (vand_f)
    
  end subroutine vandermonde1D_c

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
    real(8), intent(in) , pointer :: uvw   (  :)
    real(8), intent(out), pointer :: li    (:,:)
    real(8), intent(in) , pointer :: xout  (  :)
    !--
    integer :: i,j,nVert
    real(8) :: var
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nVert=size(xout)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! nMod=(ord+1)
   !allocate(li(1:nVert,1:nMod)) ; li(1:nVert,1:nMod)=1d0
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


  subroutine lagrange1Dv_c (order, nVtx, vand, x, lx, transpose)  BIND(C, name="SNB_lagrange1Dv") 
    use, intrinsic :: ISO_C_BINDING 
    implicit none 

    integer (C_INT), value :: order 
    integer (C_INT), value :: nVtx 
    type (C_PTR),    value :: vand
    type (C_PTR),    value :: x
    type (C_PTR),    value :: lx
    integer (C_INT), value :: transpose

    integer            :: order_f
    integer            :: nVtx_f
    real(8), pointer   :: vand_f(:,:)
    real(8), pointer   :: lx_f(:,:)
    real(8), pointer   :: x_f(:)
    logical            :: transpose_f

    integer          :: n

    n = order+1

    order_f = order
    nVtx_f = nVtx

    transpose_f = transpose
    
    call c_f_pointer (vand, vand_f, (/n,n/))
    call c_f_pointer (x, x_f, (/nVtx_f/))
    
    if (transpose_f) then
      call c_f_pointer (lx, lx_f, (/n, nVtx_f/))
    else
      call c_f_pointer (lx, lx_f, (/nVtx_f, n/))
    endif

    call lagrange1Dv(order_f, vand_f, x_f, lx_f, transpose_f)
    
  end subroutine lagrange1Dv_c


  
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
