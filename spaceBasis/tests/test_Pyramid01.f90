
subroutine pyramTestQuadrature()
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !> Routine permettant de tester les
  !> bases fonctionnelles sur une liste
  !> de points
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use pyramidRule, only: P5_gauss
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer              :: i
  integer              :: iOrd,ord
  !>
  integer              :: n
  real(8), pointer     :: x(:),y(:),z(:),w(:)
  integer              :: power
  real(8), pointer     :: f(:)
  real(8)              :: s
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Contrôle des quadratures")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ord=5
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do iOrd=0,ord
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Liste des points xyzOut
   !write(*,'(/"Points de Gauss:")')
    call P5_gauss(   &
    &    order=iOrd ,&
    &    nGauss=n   ,&
    &    uGauss=x   ,&
    &    vGauss=y   ,&
    &    wGauss=z   ,&
    &    pGauss=w    )
    
    print '(4x,"3D Quadrature")'
    do i=1,n
      print '(8x,"i",i6,": uvw=",3(f12.5,1x),2x,"p=",f12.5)',i,x(i),y(i),z(i),w(i)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(f(1:n))
   !print '("f(x,y,z)= (x+1)^",i1," (y+1)^",i1," z^",i1)',n,n,n
    power=(2*iOrd+1)/3
    do i=1,n
      f(i)=(x(i)+1)**power *(y(i)+1)**power *(z(i))**power
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Test des quadratures
    s=sum([(f(i)*w(i), i=1,n)])
    if(                 power<10  ) then
       print '("\int_0^1 \int_{-(1-z)}^{+(1-z)} \int_{-(1-z)}^{+(1-z)} (x+1)^0",i1," (y+1)^0",i1," z^0",i1," dx dy dz = ",e22.15)'&
              ,power,power,power,s
    endif
    
    if( 10<=power .and. power<100 ) then
       print '("\int_0^1 \int_{-(1-z)}^{+(1-z)} \int_{-(1-z)}^{+(1-z)} (x+1)^" ,i2," (y+1)^" ,i2," z^" ,i2," dx dy dz = ",e22.15)',&
              power,power,power,s
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(x,y,z,w,f)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramTestQuadrature



program main
  

  !> Tests des quadratures pour pyramides
  call pyramTestQuadrature()
  
  
end program main