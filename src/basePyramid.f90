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
    integer                       :: iu,iv,iw
    integer                       :: np,iNod,iOrd
    real(8)                       :: a
    real(8)                       :: r1D(0:ord)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    np=0
    do iu=1,ord+1
      np=np+iu*iu
    enddo
    if( display )print '("np=",i6)',np
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(uvw(3,np))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iNod=0
    do iw=0,ord
      
      !> a : square [-a,+a]x[-a,+a] @ iw
      if( ord==0 )then
        a=0d0
      else
        a=1d0-real(iw,kind=8)/real(ord,kind=8) 
      endif
      !if( display )print '("a=",f12.9)',a
      
      !> r1d subdivision of [-a,+a]
      r1d(0)=-a
      do iu=1,ord-iw-1
        r1d(iu)=-a+(2d0*a)*real(iu,kind=8)/real(ord-iw,kind=8)
      enddo
      r1d(ord-iw)=+a
      !if( display )print '("r1d=",15(f12.9,1x))',r1d(0:ord-iw)
      
      
      do iu=0,ord-iw
        do iv=0,ord-iw
          iNod=iNod+1
          uvw(1:3,iNod)=[ r1d(iu), r1d(iv), real(iw,kind=8)/real(ord,kind=8)  ]
          if( display )print '("uvw(",i6,")=",3(f12.9,1x))',iNod,uvw(1:3,iNod)
        enddo
      enddo  
      
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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