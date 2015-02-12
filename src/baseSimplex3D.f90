module baseSimplex3D
  use baseSimplexTools
  implicit none
  
  interface nodes3Duvw2abc ; module procedure nodes3Duvw2abc   ; end interface
  interface nodes3Duvw2abc ; module procedure nodes3Duvw2abc_0 ; end interface
  interface nodes3Duvw2abc ; module procedure nodes3Duvw2abc_1 ; end interface
  
  interface trianConnectivity ; module procedure trianglesConnectivity  ; end interface
  
  contains
  
  subroutine nodes3DOpt_2D(ord,uvw,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:,:)
    logical, intent(in)           :: display
    !>
    integer                       :: i,ad,np
    real(8)             , pointer :: uvw0(:,:)
    integer, allocatable          :: conec(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !print '("nodes3DOpt_2D")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call nodes3D   (ord=ord,uvw=uvw0,display=.false.)
    call nodes3Dopt(ord=ord,uvw=uvw0,display=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call trianglesConnectivity(ord=ord,conec=conec)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=size(conec,1)
    allocate(uvw(3,np))
    
    do i=1,size(conec,1)
      ad=conec(i,3) ! Triangle3
      uvw(1,i)=uvw0(1,ad)
      uvw(2,i)=uvw0(3,ad)
      uvw(3,i)=1d0-uvw(1,i)-uvw(2,i)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw0,conec)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Triangle unité optimisé")')
      print '("ad=",i5,2x,"u=",f19.16,2x,"v=",f19.16,2x,"w=",f19.16)',(ad,uvw(1:3,ad),ad=1,np)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes3DOpt_2D
  
  
  subroutine nodes3D(ord, uvw, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity tetrahedron
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:,:)
    logical, intent(in)           :: display
    !---
    integer                       :: iu,iv,iw,ix,ad
    integer                       :: np
    integer, pointer              :: idx(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Total number of nodes
    np=(ord+1)*(ord+2)*(ord+3)/6
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call index3D(ord=ord,idx=idx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Create equidistributed nodes on unity triangle
    allocate(uvw(1:4,1:np))
    do ad=1,np
      iu=idx(1,ad) ; iv=idx(2,ad) ; iw=idx(3,ad) ; ix=ord-iu-iv-iw
      uvw(1:4,ad)=[real(iu,kind=8)/real(ord,kind=8),& ! u
      &            real(iv,kind=8)/real(ord,kind=8),& ! v
      &            real(iw,kind=8)/real(ord,kind=8),& ! w
      &            real(ix,kind=8)/real(ord,kind=8) ] ! x
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(idx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Tetraèdre unité initial:")')
      print '("ad=",i5,2x,"u=",f12.5,2x,"v=",f12.5,2x,"w=",f12.5,2x,"x=",f12.5)',(ad,uvw(1:4,ad),ad=1,np)
      !print '("ad=",i5,2x,"u=",f12.5,2x," v=",f12.5,2x," w=",f12.5," x=",f12.5)',(ad,2d0*uvw(1:4,ad)-1,ad=1,np)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes3D
  
  subroutine nodes3Dopt(ord,uvw,display)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity triangle
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)             :: ord
    real(8), intent(inout), pointer :: uvw(:,:)
    logical, intent(in)             :: display
    
    real(8), pointer                :: x(:)
    real(8), pointer                :: y(:)
    real(8), pointer                :: z(:)
    !--
    real(8), pointer                :: xGLL(:)
    real(8)                         :: alpha
    integer                         :: iu,iv,iw,ix,ad
    integer                         :: i,j,k,n
    real(8), pointer                :: l1(:)
    real(8), pointer                :: l2(:),blend2(:),warpFactor2(:)
    real(8), pointer                :: l3(:),blend3(:),warpFactor3(:)
    real(8), pointer                :: l4(:),blend4(:),warpFactor4(:)
    real(8)                         :: v1(1:3),v2(1:3),v3(1:3),v4(1:3)
    real(8)                         :: t1(1:4,1:4),t2(1:4,1:4)
    real(8), pointer                :: la(:),lb(:),lc(:),ld(:)
    integer                         :: iFace
    real(8), pointer                :: shiftX(:),shiftY(:),shiftZ(:)
    real(8), pointer                :: warp1(:),warp2(:),blend(:),denom(:)
    real(8)                         :: array(3,3)
    real(8), pointer                :: RHS(:,:),RST(:,:)
    integer                         :: lWork
    integer                         :: ipiv(  3)
    real(8)                         :: work(192)
    integer                         :: iErr,iFail
    !---
    real(8), parameter   :: tol=1d-10
    real(8), parameter   :: alpOpt(15)=[0.0000d0, 0.0000d0 ,0.0000d0, 0.1002d0, 1.13320d0,&
    &                                   1.5608d0, 1.3413d0 ,1.2577d0, 1.1603d0, 1.10153d0,&
    &                                   0.6080d0, 0.4523d0, 0.8856d0, 0.8717d0, 0.96550d0 ]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !call nodes3D(ord=ord,uvw=uvw,display=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Total number of nodes
    n=size(uvw,2) ! print '("n=",i6)',n
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> alpha
    if( ord<16 )then
      alpha=alpOpt(ord)
    else
      alpha=1d0
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Create barycentric coordinate l1,l2,l3,l4
    allocate(l1(1:n),l2(1:n),l3(1:n),l4(1:n))
    l1(1:n)=uvw(3,1:n) ! w
    l2(1:n)=uvw(2,1:n) ! v
    l3(1:n)=uvw(4,1:n) ! x
    l4(1:n)=uvw(1,1:n) ! u
   !print '(/"[L1 L2 L3 L4]")'
   !print '(4(f12.5,2x))',(l1(ad),l2(ad),l3(ad),l4(ad),ad=1,n)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> set vertices of equilateral tetrahedron
    v1(1:3)=[-1d0,-1d0/sqr3,-1d0/sqr6]
    v2(1:3)=[ 1d0,-1d0/sqr3,-1d0/sqr6]
    v3(1:3)=[ 0d0, 2d0/sqr3,-1d0/sqr6]
    v4(1:3)=[ 0d0, 0d0     , 3d0/sqr6]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if 0==1
    open (unit=10,file='NormalizedTetra.mesh',action='write')
    write(10,'("MeshVersionFormatted 1")')
    write(10,'(/"Dimension")')
    write(10,'("3")')
    write(10,'(/"Vertices")')
    write(10,'("4")')
    write(10,'(3(f12.5,1x),2x,i1)')v1(1:3),1
    write(10,'(3(f12.5,1x),2x,i1)')v2(1:3),1
    write(10,'(3(f12.5,1x),2x,i1)')v3(1:3),1
    write(10,'(3(f12.5,1x),2x,i1)')v4(1:3),1
    write(10,'(/"Triangles")')
    write(10,'("4")')
    write(10,'("2 3 4  1")')
    write(10,'("1 4 3  2")')
    write(10,'("1 2 4  3")')
    write(10,'("1 3 2  4")')
    write(10,'(/"Tetrahedra")')
    write(10,'("1")')
    write(10,'("1 2 3 4  0")')
    write(10,'(/"End")')
    close(10)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> orthogonal axis tangents on faces 1-4
    t1(1:3,1) = v2(1:3)-v1(1:3)
    t1(1:3,2) = v2(1:3)-v1(1:3)
    t1(1:3,3) = v3(1:3)-v2(1:3)
    t1(1:3,4) = v3(1:3)-v1(1:3)
    !-
    t2(1:3,1) = v3(1:3)-5d-1*(v1(1:3)+v2(1:3))
    t2(1:3,2) = v4(1:3)-5d-1*(v1(1:3)+v2(1:3))
    t2(1:3,3) = v4(1:3)-5d-1*(v2(1:3)+v3(1:3))
    t2(1:3,4) = v4(1:3)-5d-1*(v1(1:3)+v3(1:3))
    !
    !> normalize tangents
    do ad=1,4
      t1(1:3,ad)=t1(1:3,ad)/norm2(t1(1:3,ad))
      t2(1:3,ad)=t2(1:3,ad)/norm2(t2(1:3,ad))
    enddo
    
#if 0==1
    print '(/"t1_1=",3(f12.5,1x))',t1(1:3,1)
    print '( "t1_2=",3(f12.5,1x))',t1(1:3,2)
    print '( "t1_3=",3(f12.5,1x))',t1(1:3,3)
    print '( "t1_4=",3(f12.5,1x))',t1(1:3,4)
    print '(/"t2_1=",3(f12.5,1x))',t2(1:3,1)
    print '( "t2_2=",3(f12.5,1x))',t2(1:3,2)
    print '( "t2_3=",3(f12.5,1x))',t2(1:3,3)
    print '( "t2_4=",3(f12.5,1x))',t2(1:3,4)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Warp and blend for each face (accumulated in shiftXYZ)
    !
    ! [x]   [v1(1) v2(1) v3(1) v4(1)]
    ! [y] = [v1(2) v2(2) v3(2) v4(2)] x [13,14,12,l1]^t
    ! [z]   [v1(3) v2(3) v3(3) v4(3)]
    !
    allocate(x(1:n)) ; x(1:n)=l3(1:n)*v1(1)+l4(1:n)*v2(1)+l2(1:n)*v3(1)+l1(1:n)*v4(1) ! form undeformed coordinates
    allocate(y(1:n)) ; y(1:n)=l3(1:n)*v1(2)+l4(1:n)*v2(2)+l2(1:n)*v3(2)+l1(1:n)*v4(2) ! form undeformed coordinates
    allocate(z(1:n)) ; z(1:n)=l3(1:n)*v1(3)+l4(1:n)*v2(3)+l2(1:n)*v3(3)+l1(1:n)*v4(3) ! form undeformed coordinates
    
#if 0==1
    print '(/"[x y z]")'
    print '(3(f12.5,2x))',(x(ad),y(ad),z(ad),ad=1,n)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(shiftX(1:n)) ; shiftX(1:n)=0d0
    allocate(shiftY(1:n)) ; shiftY(1:n)=0d0
    allocate(shiftZ(1:n)) ; shiftZ(1:n)=0d0
    allocate(blend(1:n))
    allocate(denom(1:n))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call gaussLegendreLobatto(ord=ord,xGLL=xGLL)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do iFace=1,4
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      select case(iFace)
      case(1) ; la=>l1 ; lb=>l2 ; lc=>l3 ; ld=>l4
      case(2) ; la=>l2 ; lb=>l1 ; lc=>l3 ; ld=>l4
      case(3) ; la=>l3 ; lb=>l1 ; lc=>l4 ; ld=>l2
      case(4) ; la=>l4 ; lb=>l1 ; lc=>l3 ; ld=>l2
      end select
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> compute warp tangential to face
      call evalShift(ord=ord, alpha=alpha, l1=lb, l2=lc, l3=ld,xGLL=xGLL, dx=warp1, dy=warp2)
      
#if 0==1
      print'(/"[warp1 warp2]")'
      print '(2(f12.6,2x))',(warp1(ad),warp2(ad),ad=1,n)
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> compute volume blending
      blend(1:n)=lb(1:n)*lc(1:n)*ld(1:n)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> modify linear blend
      denom(1:n)= (lb(1:n)+5d-1*la(1:n)) &
      &          *(lc(1:n)+5d-1*la(1:n)) &
      &          *(ld(1:n)+5d-1*la(1:n))
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do i=1,n
        if( denom(i)>tol )blend(i)=(1d0+(alpha*la(i))**2)*blend(i)/denom(i)
      enddo
      
#if 0==1
      print'(/"[blend]")'
      print '(1(f12.4,2x))',(blend(ad),ad=1,n)
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> compute warp & blend
      shiftX(1:n)= shiftX(1:n) +(blend(1:n)*warp1(1:n))*t1(1,iFace) +(blend(1:n)*warp2(1:n))*t2(1,iFace) ! x
      shiftY(1:n)= shiftY(1:n) +(blend(1:n)*warp1(1:n))*t1(2,iFace) +(blend(1:n)*warp2(1:n))*t2(2,iFace) ! y
      shiftZ(1:n)= shiftZ(1:n) +(blend(1:n)*warp1(1:n))*t1(3,iFace) +(blend(1:n)*warp2(1:n))*t2(3,iFace) ! z
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !> fix face warp 
      !ids = find(La<tol & ( (Lb>tol) + (Lc>tol) + (Ld>tol) < 3));
      !  if( la(i)<tol .and. .not.(lb(i)>tol.and.lc(i)>tol.and.ld(i)>tol) )then
      !shift(ids,:) = warp1(ids)*t1(face,:) + warp2(ids)*t2(face,:);
      do i=1,n
        if( la(i)<tol .and. .not.(lb(i)>tol.and.lc(i)>tol.and.ld(i)>tol) )then
          shiftX(i) = warp1(i)*t1(1,iFace) + warp2(i)*t2(1,iFace)
          shiftY(i) = warp1(i)*t1(2,iFace) + warp2(i)*t2(2,iFace)
          shiftZ(i) = warp1(i)*t1(3,iFace) + warp2(i)*t2(3,iFace)
        endif
      enddo
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      deallocate(warp1,warp2)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(xGLL)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> XYZ( = XYZ + shift
    x(1:n)=x(1:n)+shiftX(1:n)
    y(1:n)=y(1:n)+shiftY(1:n)
    z(1:n)=z(1:n)+shiftZ(1:n)
    
#if 0==1
    print '("[x y z]")'
    print '(3(f12.4,2x))',(x(ad),y(ad),z(ad),ad=1,n)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(blend)
    deallocate(denom)
    deallocate(shiftX,shiftY,shiftZ)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if 0==1
    open(unit=10,file="nodes3D.dat",action='write',status='unknown')
    do ad=1,n
      write(10,'(i3,1x,3(f12.5,2x))')ad,x(ad),y(ad),z(ad)
    enddo
    close(10)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> back out right tet nodes
    
    !rhs = [X';Y';Z'] - 0.5*(v2'+v3'+v4'-v1')*ones(1,length(X));
    allocate(rhs(1:3,1:n)) ; rhs(1:3,1:n)=0d0
    rhs(1,1:n)=x(1:n)-5d-1*(v2(1)+v3(1)+v4(1)-v1(1))
    rhs(2,1:n)=y(1:n)-5d-1*(v2(2)+v3(2)+v4(2)-v1(2))
    rhs(3,1:n)=z(1:n)-5d-1*(v2(3)+v3(3)+v4(3)-v1(3))
    
#if 0==1
    print '("rhs")'
    print '(3(f12.4,2x))',(rhs(1:3,ad),ad=1,n)
#endif
    
    !> array = [0.5*(v2-v1)',0.5*(v3-v1)',0.5*(v4-v1)'];  definition colonne par colonne
    array(1:3,1) = 5d-1*(v2(1:3)-v1(1:3)) ! colonne 1
    array(1:3,2) = 5d-1*(v3(1:3)-v1(1:3)) ! colonne 2
    array(1:3,3) = 5d-1*(v4(1:3)-v1(1:3)) ! colonne 3
    
#if 0==1
   !print '("array")'
   !print '(3(f12.4,2x))',(array(ad,1:3),ad=1,3)
#endif
    
    lWork=64*3
    call dgetrf(3,3,array(1,1),3,ipiv(1)              ,iErr)
    call dgetri(3  ,array(1,1),3,ipiv(1),work(1),lWork,iErr)
    
#if 0==1
   !print '("array^-1")'
   !print '(3(f12.4,2x))',(array(ad,1:3),ad=1,3)
#endif
    
    !> rst = array\[rhs]; = array Inverse[rhs]
    allocate(rst(1:n,1:3)) ; rst(1:n,1:3)=0d0
    do i=1,n
      do j=1,3
        do k=1,3
          rst(i,j)=rst(i,j)+array(j,k)*rhs(k,i)
        enddo
      enddo
    enddo
    
#if 0==1
   print '("rst")'
   print '(3(f12.4,2x))',(rst(ad,1:3),ad=1,n)
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    uvw(1,1:n)= 5d-1*( rst(:,1)                  +1d0) ! l3(1:n)
    uvw(2,1:n)= 5d-1*(          rst(:,2)         +1d0) ! l2(1:n)
    uvw(3,1:n)= 5d-1*(                   rst(:,3)+1d0) ! l4(1:n)
    uvw(4,1:n)=-5d-1*( rst(:,1)+rst(:,2)+rst(:,3)+1d0) ! l1(1:n)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(rhs,rst)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Tétra unité optimisé")')
      print '("ad=",i5,2x,"u=",f19.16,2x,"v=",f19.16,2x,"w=",f19.16,2x,"x=",f19.16)',(ad,uvw(1:4,ad),ad=1,n)
      !print '("ad=",i5,2x,"u=",f12.5,2x,"v=",f12.5,2x,"w=",f12.5,2x,"x=",f12.5)',(ad,2d0*uvw(1:4,ad)-1d0,ad=1,n)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(l1)
    deallocate(l2)
    deallocate(l3)
    deallocate(l4)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes3Dopt
  
  subroutine evalShift(ord, alpha, l1, l2, l3, xGLL, dx,dy)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)            :: ord
    real(8), intent(in)            :: alpha
    real(8), intent(in) , pointer  :: l1(:)
    real(8), intent(in) , pointer  :: l2(:)
    real(8), intent(in) , pointer  :: l3(:)
    real(8), intent(in) , pointer  :: xGLL(:)
    real(8), intent(out), pointer  :: dx(:) ! doit être alloué puisqu'on incrémente dx
    real(8), intent(out), pointer  :: dy(:) ! doit être alloué puisqu'on incrémente dy
    
    integer                         :: ad,n
    real(8), pointer                :: xout(:)
    real(8), pointer                :: blend1(:),warpFactor1(:)
    real(8), pointer                :: blend2(:),warpFactor2(:)
    real(8), pointer                :: blend3(:),warpFactor3(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Total number of nodes
    n=(ord+1)*(ord+2)*(ord+3)/6
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> compute blending function at each node for each edge
    allocate(blend1(1:n)) ; blend1(1:n)=l2(1:n)*l3(1:n)
    allocate(blend2(1:n)) ; blend2(1:n)=l1(1:n)*l3(1:n)
    allocate(blend3(1:n)) ; blend3(1:n)=l1(1:n)*l2(1:n)
   !print'(/"[blend1 blend2 blend3]")'
   !print '(3(f12.6,2x))',(blend1(ad),blend2(ad),blend3(ad),ad=1,n)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> amount of warp for each node, for each edge
    allocate(xout(1:n))
    xout(1:n)=l3(1:n)-l2(1:n) ; call warpFactor(ord=ord,xnodes=xGLL,xout=xout,warp=warpFactor1) ; warpFactor1(1:n)=4d0*warpFactor1(1:n)
    xout(1:n)=l1(1:n)-l3(1:n) ; call warpFactor(ord=ord,xnodes=xGLL,xout=xout,warp=warpFactor2) ; warpFactor2(1:n)=4d0*warpFactor2(1:n)
    xout(1:n)=l2(1:n)-l1(1:n) ; call warpFactor(ord=ord,xnodes=xGLL,xout=xout,warp=warpFactor3) ; warpFactor3(1:n)=4d0*warpFactor3(1:n)
   !print'(/"[warpFactor1 warpFactor2 warpFactor3]")'
   !print '(3(f12.6,2x))',(warpFactor1(ad),warpFactor2(ad),warpFactor3(ad),ad=1,n)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> combine blend & warp
    warpFactor1(1:n)=blend1(1:n)*warpFactor1(1:n)*(1d0+(alpha*l1(1:n))**2)
    warpFactor2(1:n)=blend2(1:n)*warpFactor2(1:n)*(1d0+(alpha*l2(1:n))**2)
    warpFactor3(1:n)=blend3(1:n)*warpFactor3(1:n)*(1d0+(alpha*l3(1:n))**2)
   !print'(/"[warp1 warp2 warp3]")'
   !print '(3(f12.6,2x))',(warpFactor1(ad),warpFactor2(ad),warpFactor3(ad),ad=1,n)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> evaluate shift in equilateral triangle
    allocate(dx(1:n)) ; dx(1:n) = 1d0*warpFactor1(1:n) + cos(2d0*pi/3d0)*warpFactor2(1:n) + cos(4d0*pi/3d0)*warpFactor3(1:n)
    allocate(dy(1:n)) ; dy(1:n) = 0d0*warpFactor1(1:n) + sin(2d0*pi/3d0)*warpFactor2(1:n) + sin(4d0*pi/3d0)*warpFactor3(1:n)
   !print'(/"[dx dy]")'
   !print '(2(f12.6,2x))',(dx(ad),dy(ad),ad=1,n)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Un peu de ménage
    deallocate(xout)
    deallocate(blend1,blend2,blend3)
    deallocate(warpFactor1,warpFactor2,warpFactor3)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    return
  end subroutine evalShift
  
  subroutine nodes3Duvw2rst(uvw,rst)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Tetra (0,0,0)-(1,0,0)-(0,1,0) To Tetra (-1,-1,-1)-(1,-1,-1)-(1,1,-1)-(-1,-1,1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: uvw(:,:)
    real(8), intent(out), pointer :: rst(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(rst(3,size(uvw,2)))
    rst(1:3,:)=2d0*uvw(1:3,:)-1d0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes3Duvw2rst
  
  subroutine nodes3Drst2abc(rst,a,b,c)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! a= -2 (1+r)/(s+t) -1
    ! b=  2 (1+s)/(1-t) -1
    ! c=  t
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(in) , pointer :: rst(:,:)
    real(8), intent(out), pointer :: a(:)
    real(8), intent(out), pointer :: b(:)
    real(8), intent(out), pointer :: c(:)
    !
    integer                       :: i,n
    real(8)                       :: d
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(rst,2) ; allocate(a(1:n),b(1:n),c(1:n))
    
    do i=1,n
      d=rst(2,i)+rst(3,i)
      if( d==0d0 )then ! Hesthaven p174
        a(i)=-1d0
      else
        a(i)=-2d0*(1d0+rst(1,i))/d-1d0
      endif
      
      if( rst(3,i)==1d0 )then
        b(i)=-1d0
      else
        b(i)=2d0*(1d0+rst(2,i))/(1d0-rst(3,i))-1d0
      endif
      
      c(i)=rst(3,i)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes3Drst2abc
  
  
  subroutine nodes3Duvw2abc(u,v,w,a,b,c,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Tetra (0,0,0)-(1,0,0)-(0,1,0) To Tetra (-1,-1,-1)-(1,-1,-1)-(1,1,-1)-(-1,-1,1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Variables globales
    real(8), intent(in)           :: u,v,w
    real(8), intent(out), pointer :: a(:),b(:),c(:)
    logical, intent(in)           :: display
    !> Variables locales
    integer                       :: i,n
    real(8)                       :: d,r,s,t
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(a(1),b(1),c(1))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Tetra
    r=2d0*u-1d0
    s=2d0*v-1d0
    t=2d0*w-1d0
    
    !
    d=s+t
    if( d==0d0 )then ! Hesthaven p174
      a(1)=-1d0
    else
      a(1)=-2d0*(1d0+r)/d-1d0
    endif
    
    if( t==1d0 )then
      b(1)=-1d0
    else
      b(1)=2d0*(1d0+s)/(1d0-t)-1d0
    endif
    
    c(1)=t
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Tetra (uvw2abc):")')
      print '("a,b,c=",3(e22.15,2x))',a(1),b(1),c(1)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes3Duvw2abc
  
  subroutine nodes3Duvw2abc_0(uvw,a,b,c,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Tetra (0,0,0)-(1,0,0)-(0,1,0) To Tetra (-1,-1,-1)-(1,-1,-1)-(1,1,-1)-(-1,-1,1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Variables globales
    real(8), intent(in) , pointer :: uvw(:,:)
    real(8), intent(out), pointer :: a(:),b(:),c(:)
    logical, intent(in)           :: display
    !> Variables locales
    integer                       :: i,n
    real(8)                       :: d,r,s,t
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(uvw,2) ; allocate(a(1:n),b(1:n),c(1:n))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do i=1,n
      r=2d0*uvw(1,i)-1d0
      s=2d0*uvw(2,i)-1d0
      t=2d0*uvw(3,i)-1d0
      !
      d=s+t
      if( d==0d0 )then ! Hesthaven p174
        a(i)=-1d0
      else
        a(i)=-2d0*(1d0+r)/d-1d0
      endif
      
      d=1d0-t
      if( d==0d0 )then
        b(i)=-1d0
      else
        b(i)=2d0*(1d0+s)/d-1d0
      endif
      
      c(i)=t
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Tetra (uvw2abc):")')
      print '("a,b,c(",i4,")=",f15.12,2x,f15.12,2x,f15.12)',(i,a(i),b(i),c(i),i=1,size(uvw,2))
      print '()'
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes3Duvw2abc_0
  
  subroutine nodes3Duvw2abc_1(u,v,w,a,b,c,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Tetra (0,0,0)-(1,0,0)-(0,1,0) To Tetra (-1,-1,-1)-(1,-1,-1)-(1,1,-1)-(-1,-1,1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Variables globales
    real(8), intent(in) , pointer :: u(:),v(:),w(:)
    real(8), intent(out), pointer :: a(:),b(:),c(:)
    logical, intent(in)           :: display
    !> Variables locales
    integer                       :: i,n
    real(8)                       :: d,r,s,t
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n=size(u) ; allocate(a(1:n),b(1:n),c(1:n))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do i=1,n
      r=2d0*u(i)-1d0
      s=2d0*v(i)-1d0
      t=2d0*w(i)-1d0
      !
      d=s+t
      if( d==0d0 )then ! Hesthaven p174
        a(i)=-1d0
      else
        a(i)=-2d0*(1d0+r)/d-1d0
      endif
      
      if( t==1d0 )then
        b(i)=-1d0
      else
        b(i)=2d0*(1d0+s)/(1d0-t)-1d0
      endif
      
      c(i)=t
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      write(*,'(/"Tetra (uvw2abc):")')
      print '("a,b,c(",i4,")=",f12.5,2x,f12.5,2x,f12.5)',(i,a(i),b(i),c(i),i=1,size(u))
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine nodes3Duvw2abc_1
  
  
  subroutine index3D(ord,idx)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    integer, intent(out), pointer :: idx(:,:)
    !
    integer                       :: np
    integer                       :: iu,iv,iw,ix,ad
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)*(ord+3)/6
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(idx(1:3,np))
    ad=0
    do iw=0,ord
      do ix=ord,0,-1
        do iv=0,ord
          do iu=0,ord
            if( iu+iv+iw+ix==ord )then
              ad=ad+1
              idx(1:3,ad)=[iu,iv,iw]
            endif
          enddo
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Verif Matlab
!    print '("index3D: WARNING TESTING MATLAB")'
!    allocate(idx(1:3,np))
!    ad=0
!    do iw=0,ord
!      do iv=0,ord-iw
!        do iu=0,ord-iw-iv
!          ad=ad+1
!          idx(1:3,ad)=[iu,iv,iw]
!        enddo
!      enddo
!    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine index3D
  
  subroutine simplex3D(ord,a,b,c,mode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> psi_{ij}(a,b) = sqr8 P_i^{0,0}(a) P_j^{2i+1,0}(b) (1-b)**i P_k^{2i+2j+2,0}(c) (1-c)**(i+j)
    !> avec {a=-2(1+r)/(s+t)-1 et b=2(1+s)/(1-t)-1 et c=t}
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => mode(1:np,1:n)
    !> Transpose = False => mode(1:n,1:np)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:),b(:),c(:)
    real(8), intent(out), pointer :: mode(:,:)
    logical, intent(in)           :: transpose
    !-
    integer                       :: i,n,ad
    integer                       :: iu,iv,iw
    integer                       :: np
    real(8), pointer              :: mode1(:),mode2(:),mode3(:)
    integer, pointer              :: idx(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(a).or..not.associated(b) )then
      print '("a or/and b not associated")'
      stop "@simplex2D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)*(ord+3)/6 ; n=size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call index3D(ord=ord,idx=idx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.transpose )then
      
      allocate(mode(1:n,1:np))
      do ad=1,np
        iu=idx(1,ad) ; iv=idx(2,ad) ; iw=idx(3,ad)
        call jacobiP(n=iu,alpha=0d0                ,beta=0d0,u=a(1:n),jf=mode1)
        call jacobiP(n=iv,alpha=2d0*real(iu   )+1d0,beta=0d0,u=b(1:n),jf=mode2)
        call jacobiP(n=iw,alpha=2d0*real(iu+iv)+2d0,beta=0d0,u=c(1:n),jf=mode3)
        
        mode(1:n,ad)= 2d0*sqr2                         &
        &            *mode1(1:n)                       &
        &            *mode2(1:n)*(1d0-b(1:n))**(iu   ) &
        &            *mode3(1:n)*(1d0-c(1:n))**(iu+iv)
        
        deallocate(mode1,mode2,mode3)
      enddo
      
    else
      
      allocate(mode(1:np,1:n))
      do ad=1,np
        iu=idx(1,ad) ; iv=idx(2,ad) ; iw=idx(3,ad)
        call jacobiP(n=iu,alpha=0d0                ,beta=0d0,u=a(1:n),jf=mode1)
        call jacobiP(n=iv,alpha=2d0*real(iu   )+1d0,beta=0d0,u=b(1:n),jf=mode2)
        call jacobiP(n=iw,alpha=2d0*real(iu+iv)+2d0,beta=0d0,u=c(1:n),jf=mode3)
        
        mode(ad,1:n)= 2d0*sqr2                         &
        &            *mode1(1:n)                       &
        &            *mode2(1:n)*(1d0-b(1:n))**(iu   ) &
        &            *mode3(1:n)*(1d0-c(1:n))**(iu+iv)
        
        deallocate(mode1,mode2,mode3)
      enddo
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(idx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine simplex3D
  
  subroutine gradSimplex3D(ord,a,b,c,drMode,dsMode,dtMode,transpose)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> psi_{ijk}(a,b,c) = sqr8 P_i^{0,0}(a)  P_j^{2i+1,0}(b) (1-b)**i  P_k^{2i+2j+2,0}(c) (1-c)**(i+j)
    !> avec {a=-2(1+r)/(s+t)-1 ; b=2(1+s)/(1-t)-1 ; c=t}
    !
    !> (∂Psi/∂r) = (∂a/∂r) (∂Psi/∂a) + (∂b/∂r) (∂Psi/∂b) + (∂c/∂r) (∂Psi/∂c)
    !> avec : (∂a/∂r)=-2/(s+t) ; (∂b/∂r)=0 ; (∂c/∂r)=0
    !> soit : (∂Psi/∂r) = -2/(s+t) (∂Psi/∂a)
    !
    !> or : b=2(1+s)/(1-c)-1 => s=(b+1)(1-c)/2-1
    !>
    !> soit : -2/(s+c) = -2/( (b+1)(1-c)/2-1 +c)
    !>                 = -1/( (b+1)(1-c)-2+2c )
    !>                 = -1/( (b+1)(1-c)-2(1-c) )
    !>                 = -1/(1-c) 1/( b+1-2)
    !>                 = -1/(1-c) 1/(b-1)
    !>
    !> donc : (∂Psi/∂r) = -1/(1-c) 1/(b-1) (∂Psi/∂a)
    !
    !> (∂Psi/∂r) = -sqr8/(1-c) 1/(b-1) (∂P_i^{0,0}/∂a)(a)  P_j^{2i+1,0}(b) (1-b)**i  P_k^{2i+2j+2,0}(c) (1-c)**(i+j)
    !>           = -sqr8 (∂P_i^{0,0}/∂a)(a) P_j^{2i+1,0}(b) (1-b)**(i-1) P_k^{2i+2j+2,0}(c) (1-c)**(i+j-1)
    !
    
    !> (∂Psi/∂s) = (∂a/∂s) (∂Psi/∂a) + (∂b/∂s) (∂Psi/∂b) + (∂c/∂s) (∂Psi/∂c)
    !> a=-2(1+r)/(s+t)-1 => (∂a/∂s) = 2(1+r)/s^2
    !> b= 2(1+s)/(1-t)-1 => (∂b/∂s) = 2/(1-t) = 2/(1-c)
    !> c=t               => (∂c/∂s) = 0
    !
    !> (∂b/∂s) =  2(1+r)/s^2
    !> or : b=2(1+s)/(1-c)-1 => s=(b+1)(1-c)/2-1
    !>      a=-2(1+r)/(s+t)-1=-2(1+r)/(s+c)-1
    !> (1+r)=-(a+1)(s+c)/2
    !> r= -(a+1)(s+c)/2 -1
    !> r= -(a+1) ( (b+1)(1-c)-2+2c )/4 -1
    !> r= -(a+1) ( (b+1)(1-c)-2(1-c) )/4 -1
    !> r= -(a+1) (1-c) (b-1)/4 -1
    !> r+1 = -(a+1) (1-c) (b-1)/4
    
    !> (∂b/∂s) =  -2(a+1) (1-c) (b-1)/4/( (b+1)(1-c)/2-1 )^2
    !> (∂b/∂s) =   -(a+1) (1-c) (b-1)/( (b+1) (1-c)/2-1 )^2/2    
    
    !> ∂Psi/∂t = (∂a/∂t) (∂Psi/∂a) + (∂b/∂t) (∂Psi/∂b) + (∂c/∂t) (∂Psi/∂c)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Transpose = True  => mode(1:np,1:n)
    !> Transpose = False => mode(1:n,1:np)
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
    integer                       :: iu,iv,iw
    integer                       :: np
    integer, pointer              :: idx(:,:)
    real(8), pointer              :: fa(:),dfa(:)
    real(8), pointer              :: fb(:),dfb(:)
    real(8), pointer              :: fc(:),dfc(:)
    real(8), pointer              :: tmp(:)
    real(8)                       :: coef
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !print '("gradSimplex3D")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=(ord+1)*(ord+2)*(ord+3)/6 ; n=size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call index3D(ord=ord,idx=idx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.transpose )then
      
      allocate(tmp(1:n))
      allocate(drMode(1:n,1:np),dsMode(1:n,1:np),dtMode(1:n,1:np))
      do ad=1,np
        iu=idx(1,ad) ; iv=idx(2,ad) ; iw=idx(3,ad)
        !print '("iu=",i3," iv=",i3," iw=",i3)',iu,iv,iw
        
        alpha=0d0                      ; beta=0d0
        call  jacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n),jf= fa)
        call dJacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n),jf=dfa)
        
        alpha=real(2*iu+1,kind=8)      ; beta=0d0
        call  jacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n),jf= fb)
        call dJacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n),jf=dfb)
        
        alpha=real(2*(iu+iv)+2,kind=8) ; beta=0d0
        call  jacobiP(n=iw,alpha=alpha,beta=beta,u=c(1:n),jf= fc)
        call dJacobiP(n=iw,alpha=alpha,beta=beta,u=c(1:n),jf=dfc)
        
        !> drMode
        drMode(1:n,ad)=dfa(1:n)*fb(1:n)*fc(1:n)
        if( iu   >0 ) drMode(1:n,ad)=drMode(1:n,ad)* (5d-1*(1d0-b(1:n)))**(iu   -1)
        if( iu+iv>0 ) drMode(1:n,ad)=drMode(1:n,ad)* (5d-1*(1d0-c(1:n)))**(iu+iv-1)
        
        !> dsMode
        dsMode(1:n,ad)=5d-1*(1d0+a(1:n))*drMode(1:n,ad)
        tmp(1:n) = dfb(1:n)*(5d-1*(1d0-b(1:n)))**iu
        if( iu   >0 ) tmp(1:n)=tmp(1:n)+(-5d-1*iu)*(fb(1:n)*(5d-1*(1d0-b(1:n)))**(iu-1d0))
        if( iu+iv>0 ) tmp(1:n)=tmp(1:n)*(5d-1*(1d0-c(1:n)))**(iu+iv-1)
        tmp(1:n)=fa(1:n)*tmp(1:n)*fc(1:n)
        dsMode(1:n,ad) =dsMode(1:n,ad)+tmp(1:n)
        
        !> dtMode
        dtMode(1:n,ad)=5d-1*(1d0+a(1:n))*drMode(1:n,ad)+5d-1*(1d0+b(1:n))*tmp(1:n)
        tmp(1:n)=dfc(1:n) *(5d-1*(1d0-c(1:n)))**(iu+iv)
        if( iu+iv>0 )then
          tmp(1:n)=tmp(1:n)-(5d-1*(iu+iv))*(fc(1:n)*(5d-1*(1d0-c(1:n)))**(iu+iv-1))
        endif
        tmp(1:n)=fa(1:n)*fb(1:n)*tmp(1:n)
        tmp(1:n)=tmp(1:n)*(5d-1*(1d0-b(1:n)))**iu
        dtMode(1:n,ad)=dtMode(1:n,ad)+tmp(1:n)
        
        !> Normalize
        coef=(2d0**( real(2*iu+iv,kind=8)+1.5d0 ) )
        drMode(1:n,ad) = drMode(1:n,ad)*coef
        dsMode(1:n,ad) = dsMode(1:n,ad)*coef
        dtMode(1:n,ad) = dtMode(1:n,ad)*coef
      enddo
      deallocate(tmp)
      
    else
      
      allocate(tmp(1:n))
      allocate(drMode(1:np,1:n),dsMode(1:np,1:n),dtMode(1:np,1:n))
      do ad=1,np
        iu=idx(1,ad) ; iv=idx(2,ad) ; iw=idx(3,ad)
        
        alpha=0d0 ; beta=0d0
        call  jacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n),jf= fa)
        call dJacobiP(n=iu,alpha=alpha,beta=beta,u=a(1:n),jf=dfa)
        
        alpha=real(2*iu+1) ; beta=0d0
        call  jacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n),jf= fb)
        call dJacobiP(n=iv,alpha=alpha,beta=beta,u=b(1:n),jf=dfb)
        
        alpha=real(2*(iu+iv)+2) ; beta=0d0
        call  jacobiP(n=iw,alpha=alpha,beta=beta,u=c(1:n),jf= fc)
        call dJacobiP(n=iw,alpha=alpha,beta=beta,u=c(1:n),jf=dfc)
        
        !> drMode
        drMode(ad,1:n)=dfa(1:n)*fb(1:n)*fc(1:n)
        if( iu   >0 ) drMode(ad,1:n)=drMode(ad,1:n)* (5d-1*(1d0-b))**(iu   -1)
        if( iu+iv>0 ) drMode(ad,1:n)=drMode(ad,1:n)* (5d-1*(1d0-c))**(iu+iv-1)
        
        !> dsMode
        dsMode(ad,1:n)=5d-1*(1d0+a(1:n))*drMode(ad,1:n)
        tmp(1:n) = dfb(1:n)*(5d-1*(1d0-b(1:n)))**iu
        if( iu   >0 ) tmp(1:n)=tmp(1:n)+(-5d-1*iu)*(fb(1:n)*(5d-1*(1d0-b(1:n)))**(iu-1   ))
        if( iu+iv>0 ) tmp(1:n)=tmp(1:n)*(                   (5d-1*(1d0-c(1:n)))**(iu+iv-1))
        tmp(1:n)=fa(1:n)*tmp(1:n)*fc(1:n)
        dsMode(ad,1:n) =dsMode(ad,1:n)+tmp(1:n)
        
        !> dtMode
        dtMode(ad,1:n)=5d-1*(1d0+a(1:n))*drMode(ad,1:n)+5d-1*(1d0+b(1:n))*tmp(1:n)
        tmp(1:n)=dfc(1:n) *(5d-1*(1d0-c(1:n)))**(iu+iv)
        if(iu+iv>0)then
          tmp(1:n)=tmp(1:n)-(5d-1*real(iu+iv,kind=8))*(fc(1:n)*((5d-1*(1d0-c(1:n)))**(iu+iv-1)))
        endif
        tmp(1:n)=fa(1:n)*fb(1:n)*tmp(1:n)
        tmp(1:n)=tmp(1:n)*((5d-1*(1d0-b(1:n)))**iu)
        dtMode(ad,1:n)=dtMode(ad,1:n)+tmp(1:n)
        
        !> Normalize
        coef=(2d0**(2*iu+iv+1.5d0))
        drMode(ad,1:n) = drMode(ad,1:n)*coef
        dsMode(ad,1:n) = dsMode(ad,1:n)*coef
        dtMode(ad,1:n) = dtMode(ad,1:n)*coef
      enddo
      
      deallocate(tmp)
      
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(idx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !print '("end gradSimplex3D")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine gradSimplex3D
  
  subroutine lagrange3Dv(ord,vand,a,b,c,lx,transpose)
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
    np=(ord+1)*(ord+2)*(ord+3)/6 ; nPt=size(a)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call simplex3D(ord=ord,a=a,b=b,c=c,mode=psi,transpose=.false.) ! psi(1:nPt,1:np)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> mat=Transpose[Vand]
    allocate(mat(np,np))
    do i=1,np
      do j=1,np
        mat(i,j)=vand(j,i)
      enddo
    enddo
    
    !> mat=Inverse[mat]
    lWork=64*(np) ; allocate(work(lWork),ipiv(np))
    call dgetrf(np,np,mat(1,1),np,ipiv(1),iErr)
    call dgetri(np,mat(1,1),np,ipiv(1),work(1),lWork,iErr)
    deallocate(ipiv,work)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> lx = Inverse[Transpose[Vand]].Psi
    if( transpose )then
      allocate(lx(1:np,1:nPt)) ; lx(:,:)=0d0
      ! psi(1:nPt,1:np)
      
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
  end subroutine lagrange3Dv
  
  
  subroutine vandermonde3D(ord,a,b,c,vand)
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
      stop "@vandermonde3D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Vandermonde matrix
    call simplex3D(ord=ord,a=a,b=b,c=c,mode=vand,transpose=.false.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine vandermonde3D
  
  subroutine gradVandermonde3D(ord,a,b,c,drVand,dsVand,dtVand)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: a(:),b(:),c(:)
    real(8), intent(out), pointer :: drVand(:,:),dsVand(:,:),dtVand(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( .not.associated(a) .or. .not.associated(b) .or. .not.associated(c) )then
      print '("3D Nodes a or/and b or/and c not associated")'
      stop "@gradVandermonde3D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Grad Vandermonde matrix
    call gradSimplex3D(    &
    &    ord=ord          ,&
    &    a=a,b=b,c=c      ,&
    &    drMode=drVand    ,&
    &    dsMode=dsVand    ,&
    &    dtMode=dtVand    ,&
    &    transpose=.false. )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine gradVandermonde3D
  
  subroutine readXYZout3D(xyzOut)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8), intent(out), pointer :: xyzOut (:,:)
    !-
    integer :: i,j,nVert
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    open(unit=10,file='Tetras.mesh',status='old',action='read')
    do i=1,5 ; read(10,*) ; enddo
    read(10,*)nVert
    allocate(xyzOut(3,nVert))
    do i=1,nVert
      read(10,*)xyzOut(1:3,i)
    enddo
    close(10)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine readXYZout3D
  
  subroutine writeSolOut3D(title,solOut)
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
      
      call system("ln -fs Tetras.mesh "// title // "_" //sfx//".mesh")
      
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
  end subroutine writeSolOut3D
  
  
  subroutine writeMesh3D(ord,uvw)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)          :: ord
    real(8), intent(in), pointer :: uvw(:,:)
    !-
    character(3)                 :: sfx
    integer                      :: i
    integer , parameter          :: iFile=100
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print '(/"writing mesh3D")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
    if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
    if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
    print '(3x,"saving ",a,a)',"mesh3D_P",sfx
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    open(unit=iFile,name="mesh3D_P" // trim(sfx) // ".mesh",status='unknown',action='write')
    write(iFile,'("MeshVersionFormatted 2"/)')
    write(iFile,'("Dimension 3"/)')
    write(iFile,'("Vertices")')
    write(iFile,*)size(uvw,2)
    do i=1,size(uvw,2)
      write(iFile,'(3(f12.5,1x),i2)')uvw(1:3,i),0
    enddo
    write(iFile,'(/"End")')
    close(iFile)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine writeMesh3D
  
  subroutine writeMeshSkin3D(ord)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use M_libmesh6_api
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)          :: ord
    !>
    integer                      :: iOrd
    integer                      :: iVert ,nVert
    integer                      :: iVert0,nVert0,nVertVol
    integer                      :: iTria,nTria
    integer                      :: iTetr,nTetr
    character(3)                 :: sfx
    integer , parameter          :: iFile=100
    real(8) , pointer            :: vert(:,:),vert0(:,:)
    integer , allocatable        :: indx(:)
    integer                      :: iVert0Min,iVertMin
    real(8)                      :: d2,d2Min
    integer , pointer            :: mark(:)
    integer , pointer            :: tetr(:,:)
    integer , pointer            :: tria(:,:)
    real(8), parameter           :: eps=1d-6
    real(4)                      :: dist(1)
    !> libmesh
    character(256)               :: name
    integer                      :: ins,ver,res,geo
    integer , allocatable        :: TypTab(:)
    real(4)                      :: xyz(3)
    integer                      :: nFld,kind(1)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define fortran 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Point de depart : Tetra P1
    print '(/"writing TetraP1.mesh")'
    open(unit=iFile,name="TetraP1.mesh",status='unknown',action='write')
    write(iFile,'( "MeshVersionFormatted 1")')
    write(iFile,'( "Dimension")')
    write(iFile,'( "3")')
    write(iFile,'(/"Vertices")')
    write(iFile,'( "4")')
    write(iFile,'( "0 0 0  1")')
    write(iFile,'( "1 0 0  1")')
    write(iFile,'( "0 1 0  1")')
    write(iFile,'( "0 0 1  1")')
    write(iFile,'(/"Triangles")')
    write(iFile,'("4")')
    write(iFile,'( "2 3 4  1")')
    write(iFile,'( "1 4 3  2")')
    write(iFile,'( "1 2 4  3")')
    write(iFile,'( "1 3 2  4")')
    write(iFile,'(/"Tetrahedra")')
    write(iFile,'( "1")')
    write(iFile,'( "1 2 3 4  0")')
    write(iFile,'(/"End")')
    close(iFile)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#if fortran==1
    open(unit=150,file='nodesTetra.f90',action='write',status='unknown')
    write(150,'("    select case(Pi)")')
#endif
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !do iOrd=1,ord
    do iOrd=ord,ord
      !>>>>>>>>
      if(   1<=iOrd .and. iOrd<  10 ) write(sfx,'("00",i1)')iOrd
      if(  10<=iOrd .and. iOrd< 100 ) write(sfx,'("0" ,i2)')iOrd
      if( 100<=iOrd .and. iOrd<1000 ) write(sfx,'(     i3)')iOrd
      !<<<<<<<<
      
      !>>>>>>>>
      !> Ecriture DEFAULT.yams pour imposer Nbiter = iOrd
      print '(/"Writing DEFAULT.yams")'
      open(unit=iFile,name="DEFAULT.yams",status='unknown',action='write')
      write(iFile,'("Nbiter",i6)')iOrd
      close(iFile)
      !<<<<<<<<
      
      !>>>>>>>>
      !> yams2_V2 va creer la peau TetraSkinPi.mesh en utilisant TetraSkinP1.meh et DEFAULT.yams
      name="TetraSkinP"//sfx//".mesh" ; print '(/"Building ",a)',trim(name)
      call system("yams2_V2 -O 4 -f TetraP1.mesh "// trim(name) // " >> ./yams.log")
      ins=GmfOpenMeshF77(trim(name),GmfRead,ver,geo) ! print '(3x,"ins=",i3)',ins
      print '(3x,"nVert=",i10)',GmfStatKwdF77(ins,GmfVertices,ver,0,TypTab)
      print '(3x,"nTetr=",i10)',GmfStatKwdF77(ins,GmfTetrahedra,0,0,TypTab)
      print '(3x,"nTria=",i10)',GmfStatKwdF77(ins,GmfTriangles,0,0,TypTab)
      res=GmfCloseMeshF77(ins)
      !<<<<<<<<
      
      !>>>>>>>>
      !> Lecture de la peau TetraSkinPi.mesh
      print '(/"Reading: ",a)',trim(name)
      ins=GmfOpenMeshF77(trim(name),GmfRead,ver,geo) ! print '(3x,"ins=",i3)',ins
      nVert=GmfStatKwdF77(ins,GmfVertices,ver,0,TypTab) ; print '(3x,"nVert=",i10)',nVert
      allocate(vert(geo+1,nVert),mark(nVert))
      res=GmfGotoKwdF77(ins,GmfVertices)
      select case(ver)
      case(1) ! real(4)
        do iVert=1,nVert
          call GmfGetVertex3dr4(ins,xyz(1),xyz(2),xyz(3),mark(iVert))
          vert(1:3,iVert)=xyz(1:3)
          vert(4,iVert)=1d0-vert(1,iVert)-vert(2,iVert)-vert(3,iVert)
        enddo
      case(2) ! real(8)
        do iVert=1,nVert
          call GmfGetVertex3dr8(ins,vert(1,iVert),vert(2,iVert),vert(3,iVert),mark(iVert))
          vert(4,iVert)=1d0-vert(1,iVert)-vert(2,iVert)-vert(3,iVert)
        enddo
      end select
      nTria=GmfStatKwdF77(ins,GmfTriangles,0,0,TypTab) ; print '(3x,"nTria=",i10)',nTria
      allocate(tria(4,nTria))
      res=GmfGotoKwdF77(ins,GmfTriangles)
      do iTria=1,nTria
        call GmfGetTriangle(ins,tria(1,iTria),tria(2,iTria),tria(3,iTria),tria(4,iTria))
      enddo
      res=GmfCloseMeshF77(ins)
      !<<<<<<<<
      
      !>>>>>>>>
      !> List of 3D equidistant Nodes
      call nodes3D(ord=iOrd,uvw=vert0,display=.false.)
      nVert0=size(vert0,2)
      !<<<<<<<<
      
      !>>>>>>>>
      !> Mise en correspondance des Noeuds de ceux construits par nodes3D (vert0) avec ceux de TetraSkinPi.mesh (vert)
      allocate(indx(nVert0)) ; indx(1:nVert0)=0
      do iVert0=1,nVert0
        iVert=1
        d2Min=1d50 ; iVert0Min=1
        loop : do
          d2  = (vert0(1,iVert0)-vert(1,iVert))**2 &
          &    +(vert0(2,iVert0)-vert(2,iVert))**2 &
          &    +(vert0(3,iVert0)-vert(3,iVert))**2
          
          if( d2<d2Min )then
            d2Min=d2
            iVert0Min=iVert
          endif
          
          if( d2<eps )then
            indx(iVert0)=iVert
            exit loop
          endif
          iVert=iVert+1
          if( iVert>nVert )exit loop
        enddo loop
      enddo
      nVertVol=count(indx(:)==0)
      print '("nVertVol=",i10)',nVertVol
      !<<<<<<<<
      
      !>>>>>>>>
      !> Sauvegarde des noeuds non situés sur la peau
      if( .not.nVertVol==0 )then
        name="nodes3DP"//sfx//".mesh" ; print '(/"Writing: ",a)',trim(name)
        ver=1
        ins=GmfOpenMeshF77(trim(name),GmfWrite,ver,geo) ; print '(3x,"nVert=",i10)',nVertVol
        res=GmfSetKwdF77(ins,GmfVertices,nVertVol,0,TypTab)
        do iVert0=1,nVert0
          if( indx(iVert0)==0 )then
            call GmfSetVertex3dr4(ins                         ,&
            &                     real(vert0(1,iVert0),kind=4),&
            &                     real(vert0(2,iVert0),kind=4),&
            &                     real(vert0(3,iVert0),kind=4),&
            &                     i3=0                         )
          endif
        enddo
        res=GmfCloseMeshF77(ins)
        !<<<<<<<<
        !>>>>>>>>
        nFld=1 ; kind(1)=1 ; dist(1)=1e0/real(iOrd,kind=4)
        name="nodes3DP"//sfx//".sol" ; print '(/"Writing: ",a)',trim(name)
        ver=1
        ins=GmfOpenMeshF77(trim(name),GmfWrite,ver,geo) ; print '(3x,"nSolu=",i10)',nVertVol
        res=GmfSetKwdF77(ins,GmfSolAtVertices,nVertVol,nFld,kind(1:1))
        do iVert0=1,nVert0
          if( indx(iVert0)==0 )then
            call gmfSetSolAtVertexR4(ins,dist(1))
          endif
        enddo
        res=GmfCloseMeshF77(ins)
      endif
      nVert0=0
      nVert =0
      
      deallocate(indx,vert,vert0,tria)
      !<<<<<<<<
      
      !>>>>>>>>
      if( nVertVol==0 )then
        call system("ghs3d -O 1  -exit 3 -in TetraSkinP"//sfx//".mesh                          -out TetraP"//sfx//".mesh > ghs3d.log")
      else
        call system("ghs3d -O 1  -exit 3 -in TetraSkinP"//sfx//".mesh -force nodes3DP"//sfx//" -out TetraP"//sfx//".mesh > ghs3d.log")
        call system("rm -f nodes3DP"//sfx//".mesh")
        call system("rm -f nodes3DP"//sfx//".sol")
      endif
      call system("rm -f TetraSkinP"//sfx//".mesh")
      !<<<<<<<<
      
      !>>>>>>>>
      name="TetraP"//sfx//".mesh"
      print '(/"Reading: ",a)',trim(name)
      ins=GmfOpenMeshF77(trim(name),GmfRead,ver,geo) ! print '(3x,"ins=",i3)',ins
      nVert=GmfStatKwdF77(ins,GmfVertices,ver,0,TypTab) ; print '(3x,"nVert=",i10)',nVert
      allocate(vert(geo+1,nVert),mark(nVert))  ! il faut mettre geo+1 pour ensuite optimiser les noeuds
      res=GmfGotoKwdF77(ins,GmfVertices)
      select case(ver)
      case(1) ! real(4)
        do iVert=1,nVert
          call GmfGetVertex3dr4(ins,xyz(1),xyz(2),xyz(3),mark(iVert))
          vert(1:3,iVert)=xyz(1:3)
          vert(4,iVert)=1d0-vert(1,iVert)-vert(2,iVert)-vert(3,iVert)
        enddo
      case(2) ! real(8)
        do iVert=1,nVert
          call GmfGetVertex3dr8(ins,vert(1,iVert),vert(2,iVert),vert(3,iVert),mark(iVert))
          vert(4,iVert)=1d0-vert(1,iVert)-vert(2,iVert)-vert(3,iVert)
        enddo
      end select
      nTetr=GmfStatKwdF77(ins,GmfTetrahedra,0,0,TypTab) ; print '(3x,"nTetr=",i10)',nTetr
      allocate(tetr(5,nTetr))
      res=GmfGotoKwdF77(ins,GmfTetrahedra)
      do iTetr=1,nTetr
        call GmfGetTetrahedron(ins,tetr(1,iTetr),tetr(2,iTetr),tetr(3,iTetr),tetr(4,iTetr),tetr(5,iTetr))
      enddo
      nTria=GmfStatKwdF77(ins,GmfTriangles,0,0,TypTab) ; print '(3x,"nTria=",i10)',nTria
      allocate(tria(4,nTria))
      res=GmfGotoKwdF77(ins,GmfTriangles)
      do iTria=1,nTria
        call GmfGetTriangle(ins,tria(1,iTria),tria(2,iTria),tria(3,iTria),tria(4,iTria))
      enddo
      res=GmfCloseMeshF77(ins)
      !<<<<<<<<
      
      !>>>>>>>>
      print '(/"Reordering Nodes")'
      call nodes3D(ord=iOrd,uvw=vert0,display=.false.)
      nVert0=size(vert0,2)
      
      allocate(indx(nVert0)) ; indx(1:nVert0)=0
      do iVert=1,nVert
        iVert0=1 ; d2Min=1d50 ; iVertMin=1
        loop1 : do
          
          d2= (vert0(1,iVert0)-vert(1,iVert))**2 &
          &  +(vert0(2,iVert0)-vert(2,iVert))**2 &
          &  +(vert0(3,iVert0)-vert(3,iVert))**2
          
          if( d2<d2Min )then
            iVert0Min=iVert0 ; d2Min=d2
          endif
          
          if( d2<eps )then
            indx(iVert)=iVert0
            exit loop1
          endif
          
          iVert0=iVert0+1
          if( iVert0>nVert0 )then
            print '("vert (",i6,")=",3(f12.5,1x),1x,"d2Min=",e12.5)',iVert    ,vert (1:3,iVert    ),d2Min
            print '("vert0(",i6,")=",3(f12.5,1x),1x,"d2Min=",e12.5)',iVert0Min,vert0(1:3,iVert0Min),d2Min
            stop"problem @ writeMeshSkin3D"
          endif
        enddo loop1
      enddo
      !<<<<<<<<
      
      !>>>>>>>>
      !> Test pour verifier que tous les noeuds sont en correspondance
      if( .not.count( indx(:)==0 )==0 )then
        print '("count(indx==0)=",i10)',count(indx(:)==0)
        stop"problem @ writeMeshSkin3D"
      endif
      !<<<<<<<<
      
      !>>>>>>>>
      !> Reecriture de TetraPi.mesh
      name="TetraP"//sfx//".mesh" ; print '(/"OverWriting: ",a)',trim(name)
      ins=GmfOpenMeshF77(trim(name),GmfWrite,1,3) ; print '(3x,"nVert=",i10)',nVert0
      res=GmfSetKwdF77(ins,GmfVertices,nVert0,0,TypTab)
      do iVert0=1,nVert0
        call GmfSetVertex3dr4(ins                         ,&
        &                     real(vert0(1,iVert0),kind=4),&
        &                     real(vert0(2,iVert0),kind=4),&
        &                     real(vert0(3,iVert0),kind=4),&
        &                     i3=0                         )
      enddo
      res=GmfSetKwdF77(ins,GmfTetrahedra,nTetr,0,TypTab) ; print '(3x,"nTetr=",i10)',nTetr
      do iTetr=1,nTetr
        call GmfSetTetrahedron(ins                ,&
        &                      indx(tetr(1,iTetr)),&
        &                      indx(tetr(2,iTetr)),&
        &                      indx(tetr(3,iTetr)),&
        &                      indx(tetr(4,iTetr)),&
        &                           tetr(5,iTetr)  )
      enddo
      res=GmfSetKwdF77(ins,GmfTriangles,nTria,0,TypTab) ; print '(3x,"nTria=",i10)',nTria
      do iTria=1,nTria
        call GmfSetTriangle(ins               ,&
        &                  indx(tria(1,iTria)),&
        &                  indx(tria(2,iTria)),&
        &                  indx(tria(3,iTria)),&
        &                       tria(4,iTria)  )
      enddo
      res=GmfCloseMeshF77(ins)
      !<<<<<<<<
      
      !>>>>>>>>
#if fortran==1
      print '(/"Write Fortran Source Vertices")'
      
      write(150,'("    case(",i3,")")')iOrd
      !write(150,'("      ")')
      do iTetr=1,nTetr
        write(150,'("      nd(1:4,",i4,")=[",i4,",",i4,",",i4,",",i4,"]")')&
        &                      iTetr              ,&
        &                      indx(tetr(1,iTetr)),&
        &                      indx(tetr(2,iTetr)),&
        &                      indx(tetr(3,iTetr)),&
        &                      indx(tetr(4,iTetr))
      enddo
#endif
      !<<<<<<<<
      
      !>>>>>>>>
      print '(/"Moving Vertices")'
      call nodes3DOpt(ord=iOrd,uvw=vert0,display=.false.)
      !<<<<<<<<
      
      !>>>>>>>>
      name="TetraP"//sfx//"Opt.mesh" ; print '(/"Writing: ",a)',trim(name)
      ins=GmfOpenMeshF77(trim(name),GmfWrite,1,3) ; print '(3x,"nVert=",i10)',nVert0
      res=GmfSetKwdF77(ins,GmfVertices,nVert0,0,TypTab)
      do iVert0=1,nVert0
        call GmfSetVertex3dr4(ins                         ,&
        &                     real(vert0(1,iVert0),kind=4),&
        &                     real(vert0(2,iVert0),kind=4),&
        &                     real(vert0(3,iVert0),kind=4),&
        &                     i3=0                         )
      enddo
      res=GmfSetKwdF77(ins,GmfTetrahedra,nTetr,0,TypTab) ; print '(3x,"nTetr=",i10)',nTetr
      do iTetr=1,nTetr
        call GmfSetTetrahedron(ins                ,&
        &                      indx(tetr(1,iTetr)),&
        &                      indx(tetr(2,iTetr)),&
        &                      indx(tetr(3,iTetr)),&
        &                      indx(tetr(4,iTetr)),&
        &                           tetr(5,iTetr)  )
      enddo
      res=GmfSetKwdF77(ins,GmfTriangles,nTria,0,TypTab) ; print '(3x,"nTria=",i10)',nTria
      do iTria=1,nTria
        call GmfSetTriangle(ins               ,&
        &                  indx(tria(1,iTria)),&
        &                  indx(tria(2,iTria)),&
        &                  indx(tria(3,iTria)),&
        &                       tria(4,iTria)  )
      enddo
      res=GmfCloseMeshF77(ins)
      !<<<<<<<<
      
      !>>>>>>>>
      print '(/"Cleanning memory")'
      nVert=0 ; nVert0=0
      nTetr=0
      nTria=0
      deallocate(vert,vert0)
      deallocate(tetr)
      deallocate(tria)
      if( allocated(indx) )deallocate(indx)
      !<<<<<<<<
      
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#if fortran==1
    write(150,'("    end select")')
    close(150)
#endif
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef fortran
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine writeMeshSkin3D
  
  
  subroutine trianglesConnectivity(ord,conec)!,trian,uvw)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Cette procedure extrait les connec-
    !> tivités des 4 triangles entourant le
    !> tétra.
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define trianglesConnectivity 0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer,              intent(in)  :: ord
    integer, allocatable, intent(out) :: conec (:,:)
    !>
    real(8) , allocatable :: uvw   (:,:)
    integer , allocatable :: conec0(:,:)
    integer               :: iu,iv,iw,ix
    integer               :: i,j,k,l,n,s
    integer               :: n1=0,n2=0,n3=0,n4=0
    integer , allocatable :: dg(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(conec0((ord+1)*(ord+2)/2,4))
    allocate(dg    ((ord+1)*(ord+2)/2  ))
    
    s=0
    do iw=0,ord
      !
      do ix=ord,0,-1
        do iv=0,ord
          do iu=0,ord
            if( iu+iv+iw+ix==ord )then
              s=s+1
              if( ix==0 )then
                n1=iv+iw*(ord+1)+1-(iw*(iw-1))/2 ! numérotation à la space
                conec0(n1,1)=s
              endif
              if( iu==0 )then
                n2=iv+iw*(ord+1)+1-(iw*(iw-1))/2 ! numérotation à la space
                conec0(n2,2)=s
              endif
              if( iv==0 )then ; n3=n3+1
                n3=iw+iu*(ord+1)+1-(iu*(iu-1))/2 ! numérotation à la space
                conec0(n3,3)=s
              endif
              if( iw==0 )then
                n4=iu+iv*(ord+1)+1-(iv*(iv-1))/2 ! numérotation à la space
                conec0(n4,4)=s
              endif
            endif
          enddo
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(conec((ord+1)*(ord+2)/2,4))
    do i=1,4
      select case(i)
      case(1) ; call permutation(order=ord, move=0, flip=.false., dg=dg)
      case(2) ; call permutation(order=ord, move=0, flip=.true. , dg=dg)
      case(3) ; call permutation(order=ord, move=0, flip=.true. , dg=dg)
      case(4) ; call permutation(order=ord, move=0, flip=.true. , dg=dg)
      end select
      
      conec(:,i)=conec0(dg(:),i)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(conec0)
    deallocate(dg    )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if trianglesConnectivity==1
    if( ord<10 )then
      write(*,'(/"Lagrangian Tetrahedra P",i1)')ord
    else
      write(*,'(/"Lagrangian Tetrahedra P",i2)')ord
    endif
    write(*,'(/"Local DOF / side")')
    do j=1,4
      write(*,'(/"sd",i1,": ",$)')j
      do i=1,size(conec,1)
        write(*,'(i4," ",$)')conec(i,j)
      enddo
    enddo
    write(*,'()')
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#undef trianglesConnectivity
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine trianglesConnectivity
  
  
  subroutine permutation(order, move, flip, dg)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Cette procedure tourne une face triangulaire
    !> vers la gauche ou vers la droite et peut aussi
    !> la retourner
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer , intent(in)  :: order
    integer , intent(in)  :: move
    logical , intent(in)  :: flip
    integer               :: dg(:)
    !>
    integer               :: i,j,k, n
    integer               :: s
    integer , allocatable :: tab(:,:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! dimensionnement
    n=size(dg)  
    allocate(tab(0:order,0:order,0:order))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> definition du standard i j k
    s=0
    do j=0,order ; do i=0,order ; do k=0,order
      if( i+j+k==order )then
        s=s+1 ; tab(i,j,k)=s
      endif
    enddo ; enddo ; enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> application de la transformation
    dg(1:n)=0
    
    s=0
    do j=0,order
      do i=0,order
        do k=0,order
          if( i+j+k==order )then
            select case(move)
            case(-1) ; s=s+1
              select case(flip)
              case(.false.) ; dg(s)=tab(k,i,j)
              case(.true. ) ; dg(s)=tab(i,k,j)
              end select
            case( 0) ; s=s+1
              select case(flip)
              case(.false.) ; dg(s)=tab(i,j,k)
              case(.true. ) ; dg(s)=tab(j,i,k)
              end select
            case( 1) ; s=s+1
              select case(flip)
              case(.false.) ; dg(s)=tab(j,k,i)
              case(.true. ) ; dg(s)=tab(k,j,i)
              end select
            end select
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(tab)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine permutation
  
end module baseSimplex3D