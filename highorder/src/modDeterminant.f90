module modDeterminant

contains

!The subroutine TSRGT applies to input real square matrix A(n,n) the upper
!triangularization algorithm of Gauss method with full pivoting and keeps
!trace of successive transformations done in integer vectors KP and LP.
!-----------------------------------------------------------------------------
!  Input parameters:
!  eps        precision (real*8)
!  n          size of A matrix (integer)
!  A          pointer to input real square matrix (real*8)
!  Output parameters:
!  it         flag=1 if A matrix ok, =0 if A matrix is singular (integer)
!  C          pointer to table storing main diagonal elements and supra-
!             diagonal elements of upper triangular matrix and the multi-
!             plying coefficients used during triangularization process
!  KP         table storing informations concerning the column exchanges
!             during process (integer)
!  LP         table storing informations concerning the line exchanges
!             during process (integer)
!-----------------------------------------------------------------------------
!The table C is first initialized to A matrix, then receives at each step k
!of the triangularization process, usefull elements of A matrix at step k for
!k=1,2,...n.
!The variables po(real*8), lo and ko(integer) store respectively pivot at step k,
!its line number and its column number.
!------------------------------------------------------------------------------
Subroutine TSRGT(eps, n, A, it, C, Kp, Lp)
  real*8 eps
  integer n,it
  real*8 A(n,n), C(n,n)
  integer Kp(n),Lp(n) 
  real*8  po,t0
  C=A; it=1; k=1
  do while (it==1.and.k<n)
    po=C(k,k); lo=k; ko=k
    do i=k, n
      do j=k, n
        if (dabs(C(i,j))>dabs(po)) then
          po=C(i,j); lo=i; ko=j
        end if
      end do
    end do
    Lp(k)=lo; Kp(k)=ko
    if (dabs(po)<eps) then
      it=0
    else
      if (lo.ne.k) then
        do j=k, n
          t0=C(k,j); C(k,j)=C(lo,j); C(lo,j)=t0
        end do
      end if
      if (ko.ne.k) then
        do i=1, n
          t0=C(i,k); C(i,k)=C(i,ko); C(i,ko)=t0
        end do
      end if 
      do i=k+1, n
        C(i,k)=C(i,k)/po
        do j=k+1, n
          C(i,j)=C(i,j)-C(i,k)*C(k,j)
        end do 
      end do
      k=k+1
    end if
  end do
  if (it==1.and.dabs(C(n,n))<eps)  it=0
  return
End Subroutine TSRGT

!The function DMGT returns the determinant of a real square matrix
!A(n,n) by Gauss method with full pivoting.
!----------------------------------------------------------------------------
!  Input parameters:
!  eps        precision (real*8)
!  n          size of A matrix (integer)
!  A          pointer to input real square matrix
!  Output parameters:
!  None
!-----------------------------------------------------------------------------
!The procedure TSRGT is used to reduce A matrix to an upper triangular matrix.
!Output variables are it(integer), C(n,n), Kp(n) and Lp(n).
!If it=0, matrix A is singular, if it=1, matrix A is regular. Table C contains
!at location i,j (j>=i) the corresponding element of the upper triangular matrix.
!Tables Lp and Kp contain informations relative to exchanges of line or column
!that occured during the process. For instance, the element number k of Lp is
!an integer <> k if an exchange of line has been made at step k (k=1,2,...,n).
!The number of exchanges of lines and columns is stored in integer L. the
!determinant of A matrix is stored in d0 (real*8).
!-----------------------------------------------------------------------------
real*8 Function DMGT(eps, n, A)
integer n
real*8 eps, A(n,n)
real*8 d0

real*8, pointer :: C(:,:)
integer,pointer :: Kp(:), Lp(:)

!allocate local matrix C and vectors Kp, Lp
  allocate(C(n,n),STAT=ialloc)
  allocate(Kp(n),STAT=ialloc)
  allocate(Lp(n),STAT=ialloc)

  call TSRGT(eps,n,A,it,C,Kp,Lp)  !call triangularization subroutine
  if (it==0) then
    d0=0.d0  !matrix singular, det=0
  else       !matrix regular, det<>0
    d0=1.d0
    do k=1, n
	  d0=d0*C(k,k)
    end do
    l=0
    do k=1, n-1
      if (Lp(k).ne.k)  l=l+1
      if (Kp(k).ne.k)  l=l+1
    end do
    if (MOD(l,2).ne.0) d0=-d0  !l is odd
  end if
  DMGT=d0   !return determinant
  return
End Function DMGT


end module modDeterminant