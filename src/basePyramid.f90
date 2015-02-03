module basePyramid
  
  implicit none

contains

  subroutine pyramidMesh3D(ord,uvw)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use M_libmesh6_api
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:,:)
    !>
    integer                       :: iVert ,nVert
    real(4)                       :: dist(1)
    character(3)                  :: sfx
    !> libmesh
    character(256)                :: name
    integer                       :: ins,ver,res,geo
    integer , allocatable         :: TypTab(:)
    real(4)                       :: xyz(3)
    integer                       :: nFld,kind(1)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>
    if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
    if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
    if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
    !<<<<<<<<
  
    !>>>>>>>>
    name="nodes3DP"//sfx//".mesh" ; print '(/"Writing: ",a)',trim(name)
    ver=1
    nVert=size(uvw,2)
    ins=GmfOpenMeshF77(trim(name),GmfWrite,ver,geo) ; print '(3x,"nVert=",i10)',nVert
    res=GmfSetKwdF77(ins,GmfVertices,nVert,0,TypTab)
    do iVert=1,nVert
      call GmfSetVertex3dr4(ins                      ,&
      &                     real(uvw(1,iVert),kind=4),&
      &                     real(uvw(2,iVert),kind=4),&
      &                     real(uvw(3,iVert),kind=4),&
      &                     i3=0                      )
    enddo
    res=GmfCloseMeshF77(ins)
    !<<<<<<<<
    !>>>>>>>>
    nFld=1 ; kind(1)=1 ; dist(1)=1e0/real(ord,kind=4)
    name="nodes3DP"//sfx//".sol" ; print '(/"Writing: ",a)',trim(name)
    ver=1
    ins=GmfOpenMeshF77(trim(name),GmfWrite,ver,3) ; print '(3x,"nSolu=",i10)',nVert
    res=GmfSetKwdF77(ins,GmfSolAtVertices,nVert,nFld,kind(1:1))
    do iVert=1,nVert
      call gmfSetSolAtVertexR4(ins,dist(1))
    enddo
    res=GmfCloseMeshF77(ins)
    !<<<<<<<<

  
    return
  end subroutine pyramidMesh3D
  
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
    integer                       :: iu,iv,iw,iOrd
    integer                       :: iNod,nNod
    integer                       :: iCel,nCel
    real(8)                       :: a
    real(8)                       :: sub(0:ord)
    integer, allocatable          :: conec(:,:)
    logical                       :: mesh
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nNod=0
    do iu=1,ord+1
      nNod=nNod+iu*iu
    enddo
    if( display )print '("nDeg=",i6)',nNod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nCel=0
    do iw=ord,1,-1
      nCel=nCel+(iw*iw)+((iw-1)*(iw-1)) !> Pyramides droites + pyramides retournees
    enddo
    if( display )print '("nCel=",i6)',nCel
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(uvw(3,nNod))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iNod=0
    do iw=0,ord
      
      !>>>>>> square [-a,+a]x[-a,+a] @ level iw
      if( ord==0 )then
        a=0d0
      else
        a=1d0-real(iw,kind=8)/real(ord,kind=8) 
      endif
     !if( display )print '("a=",f12.9)',a
      !<<<<<<
      
      !>>>>>> subdivision of [-a,+a]
      sub(0)=-a
      do iu=1,ord-iw-1
        sub(iu)=-a+(2d0*a)*real(iu,kind=8)/real(ord-iw,kind=8)
      enddo
      sub(ord-iw)=+a
     !if( display )print '("sub=",15(f12.9,1x))',sub(0:ord-iw)
      !<<<<<<
      
      !>>>>>> 2D Grid (sub(iu),sub(iv))
      iOrd=ord-iw
      do iv=0,iOrd
        do iu=0,iOrd
          iNod=iNod+1
          uvw(1:3,iNod)=[ sub(iu), sub(iv), real(iw,kind=8)/real(ord,kind=8)  ]
        enddo
      enddo  
      !<<<<<<
      
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      print '(/"Nodes coordinates")'
      do iNod=1,nNod
        print '("uvw(",i6,")=",3(f12.9,1x))',iNod,uvw(1:3,iNod)
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    mesh=.false.
    if( mesh )then
      allocate(conec(6,nCel))
      nCel=0
      do iw=0,ord-1
       !print '("iw=",i3)',iw
        do iv=0,ord-1-iw
          do iu=0,ord-1-iw
            nCel=nCel+1
            conec(1:5,nCel)=[ idx(iu  ,iv  ,iw  ),& 
            &                 idx(iu+1,iv  ,iw  ),&
            &                 idx(iu+1,iv+1,iw  ),&
            &                 idx(iu  ,iv+1,iw  ),&
            &                 idx(iu  ,iv  ,iw+1) ]
           !print '(3x,"nCel=",i3,1x,"conec=",5(i3,1x))',nCel,conec(1:5,nCel)
           
           if( .not.iw==0 )then
             nCel=nCel+1
             conec(1:5,nCel)=[ idx(iu  ,iv  ,iw  ),& 
             &                 idx(iu+1,iv  ,iw  ),&
             &                 idx(iu+1,iv+1,iw  ),&
             &                 idx(iu  ,iv+1,iw  ),&
             &                 idx(iu+1,iv+1,iw-1) ]
           endif
           
          enddo
        enddo
      enddo
      
      open(unit=10,file="pymarid.mesh",action='write')
      write(10,'( "MeshVersionFormatted 2")' )
      
      write(10,'(/"Dimension")' )
      write(10,'( "3")' )
      
      write(10,'(/"Vertices")' )
      write(10,'(i10)')nNod
      do iNod=1,nNod
        write(10,'(3(e22.15,1x),i3)')uvw(1:3,iNod),0
      enddo
      
      write(10,'(/"Pyramids")' )
      write(10,'(i10)')nCel
      do iCel=1,nCel
        write(10,'(5(i6,1x),i3)')conec(1:5,iCel),0
      enddo
      
      write(10,'(/"End")')
      close(10)
       
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
    
    contains
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer function idx(iu,iv,iw)
      !---
      integer :: iu,iv,iw
      integer :: k
      !---
      idx=0
      do k=0,iw-1
        idx=idx+(ord+1 -k)**2
      enddo
      print '("iu,iv,iw=",3(i3,2x),"idx=",i3)',iu,iv,iw,idx
      idx=iu+iv*(ord+1-iw)+idx +1
      return
    end function
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
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