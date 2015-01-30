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
    if( display )print '("nDeg=",i6)',np
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
      
      !> computation of nodes
      do iv=0,ord-iw
        do iu=0,ord-iw
          iNod=iNod+1
          uvw(1:3,iNod)=[ r1d(iu), r1d(iv), real(iw,kind=8)/real(ord,kind=8)  ]
        enddo
      enddo  
      
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      print '(/"Nodes coordinates")'
      do iNod=1,np
        print '("uvw(",i6,")=",3(f12.9,1x))',iNod,uvw(1:3,iNod)
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    return
  end subroutine pyramidEquiNodes3D
  
  subroutine pyramidSides3D(ord, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity pyramid
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    logical, intent(in)           :: display
    !---
    integer                       :: iu,iv,iw
    integer                       :: iNod,iSide
    integer                       :: side((ord+1)*(ord+1)+4*(ord+1)*(ord+2)/2 )
    integer                       :: sideIdx(6)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sideIdx(1)=0
    sideIdx(2)=sideIdx(1)+(ord+1)*(ord+1)   !> square
    sideIdx(3)=sideIdx(2)+(ord+1)*(ord+2)/2 !> triangle
    sideIdx(4)=sideIdx(3)+(ord+1)*(ord+2)/2 !> triangle
    sideIdx(5)=sideIdx(4)+(ord+1)*(ord+2)/2 !> triangle
    sideIdx(6)=sideIdx(5)+(ord+1)*(ord+2)/2 !> triangle
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iNod=0
    do iw=0,ord
      do iv=0,ord-iw
        do iu=0,ord-iw
          
          iNod=iNod+1
          
          !> side1 iw=0
          if( iw==0 )then
            sideIdx(1)=sideIdx(1)+1 ! print '("sideIdx(1)=",i2)',sideIdx(1)
            side( sideIdx(1) )=iNod
          endif
          
          !> side2 iu=0
          if( iu==0 )then   
            sideIdx(2)=sideIdx(2)+1
            side( sideIdx(2) )=iNod
          endif
          
          !> side3 iv=0
          if( iv==0 )then   
            sideIdx(3)=sideIdx(3)+1
            side( sideIdx(3) )=iNod
          endif
          
          !> side4 iu=ord
          if( iu==ord-iw )then   
            sideIdx(4)=sideIdx(4)+1
            side( sideIdx(4) )=iNod
          endif
          
          !> side5 iv=ord
          if( iv==ord-iw )then   
            sideIdx(5)=sideIdx(5)+1 ! print '("sideIdx(5)=",i2)',sideIdx(5)
            side( sideIdx(5) )=iNod
          endif
          
        enddo
      enddo      
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sideIdx(1)=0
    sideIdx(2)=sideIdx(1)+(ord+1)*(ord+1)
    sideIdx(3)=sideIdx(2)+(ord+1)*(ord+2)/2
    sideIdx(4)=sideIdx(3)+(ord+1)*(ord+2)/2
    sideIdx(5)=sideIdx(4)+(ord+1)*(ord+2)/2
    sideIdx(6)=sideIdx(5)+(ord+1)*(ord+2)/2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      print '(/"Degrees of freedom/side")'
      do iSide=1,5
        print '("Side",i1,": ",$)',iSide
        do iNod=sideIdx(iSide)+1,sideIdx(iSide+1)
          print '(i5,1x,$)',side(iNod)
        enddo
        print '()'
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidSides3D
  


end module basePyramid