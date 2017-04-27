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
  order=10
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)order
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



program main
  
  !> Connectivite
  call testConnectivitiesOfSides()
 !call jacobiTest()
  
end program main