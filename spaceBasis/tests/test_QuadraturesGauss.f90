subroutine testQuadratureGL()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use baseSimplexTools, only: gaussLobattoQuadratures,gaussLegendreQuadratures
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  real(8), pointer   :: xGL (:),wGL (:)
  real(8), pointer   :: xGLL(:),wGLL(:)
  real(8), pointer   :: xGJ (:),wGJ (:)
  integer            :: i,order
  real(8)            :: int
  real(8), parameter :: resultat=-332d0/105d0
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '(/"Gauss Legendre Quadatures")'
  order=4
  call gaussLegendreQuadratures(ord=order,xGL=xGL,wGL=wGL)
  print '("i=",i2," xGL=",f12.5," wGL=",f12.5)',(i,xGL(i),wGL(i),i=1,size(xGL))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  int=0d0
  do i=1,order+1
    int=int+wGL(i)*func( xGL(i) )
  enddo
  
  print '("\int_{-1}^{+1} 2 - 3x + x^2 + 2x^3 - 6x^4 + 8x^5- 19x^6 + x^7 dx = ",e22.15)',int
  print '("Erreur=",e22.15)',int-resultat
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '(/"Gauss Legendre-Lobatto Quadatures")'
  order=4
  call gaussLobattoQuadratures(ord=order,xGL=xGLL,wGL=wGLL)
  print '("i=",i2," xGL=",f12.5," wGL=",f12.5)',(i,xGLL(i),wGLL(i),i=1,size(xGLL))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  int=0d0
  do i=1,order+1
    int=int+wGLL(i)*func( xGLL(i) )
  enddo
  
  print '("\int_{-1}^{+1} 2 - 3x + x^2 + 2x^3 - 6x^4 + 8x^5- 19x^6 + x^7 dx = ",e22.15)',int
  print '("Erreur=",e22.15)',int-resultat
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(xGL ,wGL )
  deallocate(xGLL,wGLL)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
  
contains
  
  function func(x) result(fx)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8) :: x,fx
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    fx=2d0-3d0*x+x**2+2d0*x**3-6d0*x**4+8d0*x**5-19d0*x**6+x**7

    return
  end function func
  
end subroutine testQuadratureGL

program main
  implicit none
  call testQuadratureGL()
end program main

