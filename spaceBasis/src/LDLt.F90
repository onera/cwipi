module space_LDLt
  use omp_lib
  
  implicit none ; private
  
  interface Determinant        ; module procedure LDLt_determinant_r                 ; end interface
  interface Det                ; module procedure LDLt_determinant_r                 ; end interface
  public :: Determinant,Det
  
  interface Factorise          ; module procedure LDLt_factorise_r                   ; end interface
  public :: Factorise
  
  interface LDLt_test          ; module procedure LDLt_test_r                        ; end interface
  public :: LDLt_test
  
  interface Solve              ; module procedure LDLt_solve_r                       ; end interface
  interface Solve              ; module procedure LDLt_solve_c                       ; end interface
  public :: Solve
  
contains
  
  function LDLt_determinant_r(n,mat) result(det)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)  :: n
    real(8)              :: mat(:)
    real(8)              :: det
    !>
    integer              :: i
    integer              :: ad
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    det=1d0
    ad=1
    do i=0,n-1
      det=det*mat(ad)
     !write(*,'("ad=",i2)')ad
      ad=ad+n-i
    enddo
    det=1d0/det
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end function LDLt_determinant_r
  
  
  subroutine LDLt_factorise_r(n,mat)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! mat(  1)  
    ! mat(  2)  mat( n+1)
    ! mat(  3)  mat( n+2) mat(2n  )
    ! ...
    ! mat(n-1)  mat(2n-2) mat(3n-3)
    ! mat(n  )  mat(2n-1) mat(3n-2) ... mat(nn-n)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)  :: n
    real(8)              :: mat(:)
    !>
    integer              :: i,j,k,l,iRow,jRow
    integer              :: ki,kj
    real(8), allocatable :: tmp(:)
    integer              :: row0,row1,drow,drow1
    integer              :: row2,row3,drow2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !write(*,'("subroutine LDLt_factorise_r")')
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(tmp(1:size(mat)))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Factorisation LDLT
    row0=1 ; drow=n
    do i=1,n
      !> preliminaires
      call dcopy(drow,mat(row0),1,tmp(row0),1) ! tmp(row0:row0+drow)=mat(row0:row0+drow)
      ! operations on col i
      mat(row0)=1d0/mat(row0)
      drow=drow-1
      ! mat(row+1:row0+drow)=mat(row0)*mat(row0+1:row0+drow)
      do iRow=1,drow
        jRow=row0+iRow
        mat(jRow)=mat(row0)*mat(jRow)
      enddo
     !call dscal(drow,mat(row0),mat(row0+1),1) ! mat(row+1:row0+drow)=mat(row0)*mat(row0+1:row0+drow)
      ki=row0+1
      row0=row0+drow+1
      ! operations on col i+1 to n
      row1=row0 ; drow1=drow
      do k=i+1,n
        call daxpy(drow1,-tmp(ki),mat(ki),1,mat(row1),1) ! mat(row1:)=-tmp(ki)*mat(ki:)+mat(row1:)
        ki=ki+1
        row1=row1+drow1
        drow1=drow1-1
      enddo
     !print '("LDLt_factorise_r i=",i6," row0=",i," drow=",i6," size(mat)=",i6)',i,row0,drow,size(mat)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(tmp)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine LDLt_factorise_r
  
  
  subroutine LDLt_solve_c(n,m,mat,b)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Resolution du systeme lineaire
    !   A x = b
    ! La matrice A a ete factorise sous la forme
    ! LDLt et la solution du systeme lineaire est
    ! retournee dans le vecteur b
    !
    ! ker = nombre de second member
    ! n   = nombre d'equations et d'inconnues
    ! m   = l=n*(n+1)/2
    !
    ! mat(  1)  
    ! mat(  2)  mat( n+1)
    ! mat(  3)  mat( n+2) mat(2n  )
    ! ...
    ! mat(n-1)  mat(2n-2) mat(3n-3)
    ! mat(n  )  mat(2n-1) mat(3n-2) ... mat(nn-n)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !use space_config        , only: calc,ker
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in) :: n
    integer, intent(in) :: m
    real   (8)          :: mat(:)
    complex(8)          :: b(:)
    !>
    integer             :: stride
    integer             :: i,j,k,l
    integer             :: i0,i1
    integer             :: j0,j1
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    stride=3 ! ker(ob=calc)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> descente
    !
    ! y_i = x_i - \sum_{k=1}^{i-1} a_{ik} y_k
    ! y_i = x_i - a_{i1} y_1 - a_{i2} y_2 - a_{i1} y_3 - ...
    ! y_j = x_j - a_{j1} y_1 - a_{j2} y_2 - a_{j1} y_3 - ...
    !
    ! y(1)=x(1) ; do i=2,n ; y(i)=y(i)- a(i,1) x(1) ; enddo 
    ! y(2)=x(2) ; do i=3,n ; y(i)=y(i)- a(i,2) x(1) ; enddo
    ! y(3)=x(3) ; do i=4,n ; y(i)=y(i)- a(i,3) x(1) ; enddo
    !
    k=0
    do i=1,n-1
      k=k+1
      i0=stride*(i-1)
      do j=i+1,n
        k=k+1
        j0=stride*(j-1)
        do l=1,stride
          i1=i0+l
          j1=j0+l
          b(j1)=b(j1)-mat(k)*b(i1)
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> diagonale
    !> y_i = x_i/d_i
    k=1
    do i=0,n-1
      i0=stride*i
      do l=1,stride
        i1=i0+l
        b(i1)=b(i1)*mat(k)
      enddo
      k=k+n-i
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> remontee
    ! y_i = x_i - \sum_{k=i+1}^{n} a_{ki} x_k
    !
    k=n*(n+1)/2+1
    do i=n,1,-1
      i0=stride*i
      do j=n,i+1,-1
        k=k-1
        j0=stride*j
        do l=0,stride-1,1
          i1=i0-l
          j1=j0-l
          b(i1)=b(i1)-mat(k)*b(j1)
        enddo
      enddo
      k=k-1
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine LDLt_solve_c

  subroutine LDLt_solve_r(n,m,mat,b)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Resolution du systeme lineaire
    !   A x = b
    ! La matrice A a ete factorise sous la forme
    ! LDLt et la solution du systeme lineaire est
    ! retournee dans le vecteur b
    !
    ! ker = nombre de second member
    ! n   = nombre d'equations et d'inconnues
    ! m   = l=n*(n+1)/2
    !
    ! mat(  1)  
    ! mat(  2)  mat( n+1)
    ! mat(  3)  mat( n+2) mat(2n  )
    ! ...
    ! mat(n-1)  mat(2n-2) mat(3n-3)
    ! mat(n  )  mat(2n-1) mat(3n-2) ... mat(nn-n)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !use space_config        , only: calc,ker
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(8)             :: mat(:)
    real(8)             :: b  (:)
    !>
    integer             :: stride
    integer             :: i,j,k,l
    integer             :: i0,i1
    integer             :: j0,j1
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    stride=3 ! ker(ob=calc)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! descente
    !
    ! y_i = x_i - \sum_{k=1}^{i-1} a_{ik} y_k
    ! y_i = x_i - a_{i1} y_1 - a_{i2} y_2 - a_{i1} y_3 - ...
    ! y_j = x_j - a_{j1} y_1 - a_{j2} y_2 - a_{j1} y_3 - ...
    !
    ! y(1)=x(1) ; do i=2,n ; y(i)=y(i)-a(i,1) x(1) ; enddo 
    ! y(2)=x(2) ; do i=3,n ; y(i)=y(i)-a(i,2) x(1) ; enddo
    ! y(3)=x(3) ; do i=4,n ; y(i)=y(i)-a(i,3) x(1) ; enddo
    !
    k=0
    do i=1,n-1
      k=k+1
      i0=stride*(i-1)+1
      do j=i+1,n
        k=k+1
        j0=stride*(j-1)+1
        call daxpy(stride,-mat(k),b(i0),1,b(j0),1) ! b(j0:)=-mat(k)*b(i0:)+b(j0:)
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> diagonale
    !> y_i = x_i/d_i
    !
    k=1
    do i=0,n-1
      i0=stride*i+1
      call dscal(stride,mat(k),b(i0),1) ! b(i0:)=b(i0:)*mat(k)
      k=k+n-i
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! remontee
    ! y_i = x_i - \sum_{k=i+1}^{n} a_{ki} x_k
    !
    k=m+1
    do i=n,1,-1
      i0=stride*i
      do j=n,i+1,-1
        k=k-1
        j0=stride*j
        do l=0,stride-1,1
          i1=i0-l
          j1=j0-l
          b(i1)=b(i1)-mat(k)*b(j1)
        enddo
      enddo
      k=k-1
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine LDLt_solve_r

  subroutine LDLt_test_r
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Matrice compactee par colonne :
    !--------------------------------
    ! mat(  1)  
    ! mat(  2)  mat( n+1)
    ! mat(  3)  mat( n+2) mat(2n  )
    ! ...
    ! mat(n-1)  mat(2n-2) mat(3n-2)
    ! mat(n  )  mat(2n-1) mat(3n-3) ... mat((nn+n)/2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !use space_config        , only: calc,ker
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer              :: stride
    integer              :: n,m
    real(8)              :: mat(21)
    real(8), allocatable ::   b(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    stride=3 ! ker(ob=calc)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(*,*)'TEST 1'
    n=3
    write(*,*)'ker= ',stride
    write(*,*)'n  = ',n

    mat(1:3)=(/0.083333333333333333333, 0.041666666666666666667, 0.041666666666666666667/) 
    mat(4:5)=(/                         0.083333333333333333333, 0.041666666666666666667/)
    mat(6:6)=(/                                                  0.083333333333333333333/) 

    write(*,'(a)')'Matrice sysmetrique'
    write(*,'(    10f10.2,1x)')mat(1:3)
    write(*,'(10x,10f10.2,1x)')mat(4:5)
    write(*,'(20x,10f10.2,1x)')mat(6:6)
    !
    allocate(b(n*stride))
    !
    b(         1:  stride)= 1d1 ; b(2)= 2d1 ; b(3)= 3d1
    b(  stride+1:2*stride)=-2d0
    b(2*stride+1:3*stride)= 3d0
    !
    write(*,*)' '
    write(*,*)'initials'
    write(*,'(4(f10.2,1x))')b(         1:  stride)
    write(*,'(4(f10.2,1x))')b(  stride+1:2*stride)
    write(*,'(4(f10.2,1x))')b(2*stride+1:3*stride)
    write(*,*)' '
    call LDLt_factorise_r(n,mat)
    !
    write(*,'(a)')'Matrice Factorisee LDLT'
    write(*,'(    10f10.2,1x)')mat(1:3)
    write(*,'(10x,10f10.2,1x)')mat(4:5)
    write(*,'(20x,10f10.2,1x)')mat(6:6)
    !
    m=n*(n+1)/2
    call LDLt_solve_r(n,m,mat,b)
    !
    write(*,*)' '
    write(*,*)'solution'
    write(*,'(4(f10.2,1x))')b(         1:  stride)
    write(*,'(4(f10.2,1x))')b(  stride+1:2*stride)
    write(*,'(4(f10.2,1x))')b(2*stride+1:3*stride)
    write(*,*)' '
    !
    deallocate(b)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(*,*)' '
    write(*,*)'TEST 2'
    n=6

    write(*,*)'ker= ',stride
    write(*,*)'n  = ',n

    mat( 1: 6)=(/1d0, 2d0, 3d0, 4d0, 5d0, 6d0/)
    mat( 7:11)=(/     2d0, 3d0, 4d0, 5d0, 5d0/)
    mat(12:15)=(/          3d0, 4d0, 5d0, 4d0/)
    mat(16:18)=(/               4d0, 5d0, 3d0/)
    mat(19:20)=(/                    5d0, 2d0/)
    mat(21:21)=(/                         1d0/)
    !
    write(*,'(a)')'Matrice sysmetrique'
    write(*,'(    10f10.2,1x)')mat( 1: 6)
    write(*,'(10x,10f10.2,1x)')mat( 7:11)
    write(*,'(20x,10f10.2,1x)')mat(12:15)
    write(*,'(30x,10f10.2,1x)')mat(16:18)
    write(*,'(40x,10f10.2,1x)')mat(19:20)
    write(*,'(50x,10f10.2,1x)')mat(21:21)
    !
    allocate(b(n*stride))
    !
    b(         1:  stride)= 1d1 ; b(2)= 2d1 ; b(3)= 3d1
    b(  stride+1:2*stride)=-2d0
    b(2*stride+1:3*stride)= 3d0
    b(3*stride+1:4*stride)=-4d0
    b(4*stride+1:5*stride)= 5d0
    b(5*stride+1:6*stride)= 6d0
    !
    write(*,*)' '
    write(*,*)'initials'
    write(*,'(4(f10.2,1x))')b(         1:  stride)
    write(*,'(4(f10.2,1x))')b(  stride+1:2*stride)
    write(*,'(4(f10.2,1x))')b(2*stride+1:3*stride)
    write(*,'(4(f10.2,1x))')b(3*stride+1:4*stride)
    write(*,'(4(f10.2,1x))')b(4*stride+1:5*stride)
    write(*,'(4(f10.2,1x))')b(5*stride+1:6*stride)
    write(*,*)' '
    !
    call LDLt_factorise_r(n,mat)
    !
    write(*,'(a)')'Matrice Factorisee LDLT'
    write(*,'(    10f10.2,1x)')mat( 1: 6)
    write(*,'(10x,10f10.2,1x)')mat( 7:11)
    write(*,'(20x,10f10.2,1x)')mat(12:15)
    write(*,'(30x,10f10.2,1x)')mat(16:18)
    write(*,'(40x,10f10.2,1x)')mat(19:20)
    write(*,'(50x,10f10.2,1x)')mat(21:21)
    !
    m=n*(n+1)/2
    call LDLt_solve_r(n,m,mat,b)
    !
    write(*,*)' '
    write(*,*)'solution'
    write(*,'(4(f10.2,1x))')b(         1:  stride)
    write(*,'(4(f10.2,1x))')b(  stride+1:2*stride)
    write(*,'(4(f10.2,1x))')b(2*stride+1:3*stride)
    write(*,'(4(f10.2,1x))')b(3*stride+1:4*stride)
    write(*,'(4(f10.2,1x))')b(4*stride+1:5*stride)
    write(*,'(4(f10.2,1x))')b(5*stride+1:6*stride)
    write(*,*)' '
    !
    deallocate(b)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    stop
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine LDLt_test_r
  
end module space_LDLt
  
