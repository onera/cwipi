module basePyramid
  
  implicit none

contains


  subroutine pyramidEquiNodes3D(ord, uvw, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity pyramid
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:,:)
    logical, intent(in)           :: display
    !---
    integer                       :: i,j
    real(8)                       :: a
    real(8), pointer              :: r1D(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(r1D(0:ord))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do i=0,ord
      
      if( ord==0 )then
        a=0d0
      else
        a=1d0-real(i,kind=8)/real(ord,kind=8) 
      endif
      !if( display )print '("a=",f12.9)',a
      
      r1d(0)=-a
      do j=0,ord-i-1
        r1d(j)=-a+(2d0*a)*real(j,kind=8)/real(ord-i,kind=8)
      enddo
      r1d(ord-i)=+a
      
      if( display )print '("r1d=",15(f12.9,1x))',r1d(0:ord-i)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(r1D)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
!    x = []; y = []; z= [];
!    for level = 0:N
!    a = (1-t(level+1));
!    disp(a);
!    if level < N        
!        r1D = linspace(-a,a,N+1-level);
!    else
!        r1D = 0;
!    end
!    
!    [r s] = meshgrid(r1D);
!    x = [x; r(:)];
!    y = [y; s(:)];    
!    z = [z; t(level+1)*ones(size(r(:)))];
    
    return
  end subroutine pyramidEquiNodes3D
  
  


end module basePyramid