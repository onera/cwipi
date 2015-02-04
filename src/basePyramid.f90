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
    integer                       :: iu,iv,iw
    integer                       :: iNod ,nNod
    real(4)                       :: dist
    character(3)                  :: sfx
    !> libmesh
    character(256)                :: name
    integer                       :: ins,ver,res,geo
    integer , allocatable         :: TypTab(:)
    real(4)                       :: xyz(3)
    integer                       :: nFld,kind(1)
    logical                       :: test
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
    if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
    if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Ecriture
    nNod=0
    do iw=0,ord
      do iv=0,ord-iw
        do iu=0,ord-iw
          test= (iu==0.or.iu==ord-iw).or.( iv==0.or.iv==ord-iw).or.iw==0
          if( .not.test )nNod=nNod+1
        enddo
      enddo
    enddo
    
    if( .not. nNod==0 )then    
      !> mesh
      name="nodes3DP"//sfx//".mesh" ; print '(/"Writing: ",a)',trim(name)
      open(unit=10,file=trim(name),action='write')
      write(10,'( "MeshVersionFormatted 1")' )
      write(10,'(/"Dimension"/,"3")' )
      write(10,'(/"Vertices"/,i10)' )nNod
      iNod=0
      do iw=0,ord
        do iv=0,ord-iw
          do iu=0,ord-iw
            iNod=iNod+1
            test= (iu==0.or.iu==ord-iw).or.( iv==0.or.iv==ord-iw).or.iw==0
            if( .not.test )write(10,'(3(e22.15,1x),i3)')uvw(1:3,iNod),0
          enddo
        enddo
      enddo
      write(10,'(/"End")')
      close(10)
      
      !> sol
      dist=1e0/real(ord,kind=4)      
      name="nodes3DP"//sfx//".sol" ; print '( "Writing: ",a)',trim(name)
      open(unit=10,file=trim(name),action='write')
      write(10,'( "MeshVersionFormatted 1")' )
      write(10,'(/"Dimension"/,"3")' )
      write(10,'(/"SolAtVertices"/,i10)' )nNod
      write(10,'("1",1x,"1")' ) !> 1 field 1=scalar
      write(10,'(e22.15)')((dist),iNod=1,nNod)
      write(10,'(/"End")')
      close(10)
            
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    ! yams2_V2 -f -O 0 PyramidSkinP025.mesh
    !>>>>>>>>
    print '(/"Writing: ",a)',"PyramidP"//sfx//".mesh"
    if( nNod==0 )then
      call system("ghs3d -O 1  -exit 3 -in PyramidSkinP"//sfx//".mesh -out PyramidP"//sfx//".mesh > ghs3d.log")
    else
      !call system("ghs3d -O 1  -exit 3 -in PyramidSkinP"//sfx//".mesh -out PyramidP"//sfx//".mesh > ghs3d.log")
      call system("ghs3d -O 1  -exit 3 -in PyramidSkinP"//sfx//".mesh -force nodes3DP"//sfx//" -out PyramidP"//sfx//".mesh") ! > ghs3d.log")
      print '(/"Writing: ",a)',"PyramidP"//sfx//".mesh"
    endif
    
      !call system("yams2_V2 -f -O 10 PyramidSkinP"//sfx//".mesh > yams.log")
    !call system("rm -f nodes3DP"//sfx//".mesh")
    !call system("rm -f nodes3DP"//sfx//".sol")
!    call system("rm -f PyramidSkinP"//sfx//".d.mesh")
    !<<<<<<<<
      
    return
  end subroutine pyramidMesh3D
  
  subroutine pyramidSkin3D(ord, uvw)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Cette routine construit le maillage surfacique
    !> d'une pyramide Pord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:,:)
    !---
    integer                       :: iu,iv,iw
    integer                       :: iNod,nNod
    integer                       :: iCel,nCel
    integer                       :: nQuadr,nTrian
    integer, allocatable          :: quadr(:,:)
    integer, allocatable          :: trian(:,:)
    character(3)                  :: sfx
    character(256)                :: name
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
    if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
    if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nNod=size(uvw,2)
    allocate(quadr(5,  ord*ord))
    allocate(trian(4,4*ord*ord))
    nQuadr=0
    nTrian=0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Side1
    iw=0
    do iv=0,ord-1-iw
      do iu=0,ord-1-iw
        nQuadr=nQuadr+1
        quadr(1:5,nQuadr)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),& 
        &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
        &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
        &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
        &                  1                                           ]
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Side2 
    do iw=0,ord-1
      do iv=0,ord-1-iw
        do iu=0,0
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
          &                  2                                           ]
          if( .not. iv==ord-1-iw )then
            nTrian=nTrian+1
            trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),& 
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw+1),&
            &                  2                                           ]
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Side3
    do iw=0,ord-1
      do iv=0,0
        do iu=0,ord-1-iw
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),& 
          &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  3                                           ]
          if( .not. iu==ord-1-iw )then
            nTrian=nTrian+1
            trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
            &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw+1),& 
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv , iw=iw+1),&
            &                  3                                           ]
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Side4 
    do iw=0,ord-1
      do iv=0,ord-1-iw
        do iu=ord-1-iw,ord-1-iw
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  4                                           ]            
          if( .not.iv==ord-1-iw )then
            nTrian=nTrian+1
            trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw+1),&
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
            &                  4                                           ]
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Side5
    do iw=0,ord-1
      do iv=ord-1-iw,ord-1-iw
        do iu=0,ord-1-iw
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  5                                           ]
          if( .not.iu==ord-1-iw )then
            nTrian=nTrian+1
            trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
            &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw+1),&
            &                  5                                           ]
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Ecriture
    name="PyramidSkinP"//sfx//".mesh" ; print '(/"Building ",a)',trim(name)
    open(unit=10,file=trim(name),action='write')
    write(10,'( "MeshVersionFormatted 1")' )
    write(10,'(/"Dimension"/,"3")' )
    write(10,'(/"Vertices"/,i10)' )nNod
    write(10,'(3(e22.15,1x),i3)')((uvw(1:3,iNod),0),iNod=1,nNod)
    write(10,'(/"Quadrilaterals"/,i10)' )nQuadr
    write(10,'(4(i6,1x),i3)')(quadr(1:5,iCel),iCel=1,nQuadr)
    write(10,'(/"Triangles"/,i10)' )nTrian
    write(10,'(3(i6,1x),i3)')(trian(1:4,iCel),iCel=1,nTrian)
    write(10,'(/"End")')
    close(10)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
    deallocate(quadr,trian)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidSkin3D
  
  
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
    real(8)                       :: a
    real(8)                       :: sub(0:ord)
    
    integer                       :: iCel,nPyr,nTet
    integer, allocatable          :: pyram(:,:)
    integer, allocatable          :: tetra(:,:)
    logical                       :: mesh
    logical                       :: tf
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nNod=0
    do iu=1,ord+1
      nNod=nNod+iu*iu
    enddo
    if( display )print '("nDeg=",i6)',nNod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nPyr=0
    do iw=ord,1,-1
      nPyr=nPyr+(iw*iw)+((iw-1)*(iw-1)) !> Pyramides droites + pyramides retournees
    enddo
    if( display )print '("nPyr=",i6)',nPyr
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
    mesh=.true.
    if( mesh )then
      allocate(pyram(6,nPyr))
      allocate(tetra(6,1000))
      nPyr=0
      nTet=0
      tf=.false.
      do iw=0,ord-1
        do iv=0,ord-1-iw
          print '("tf=",l)',tf
          do iu=0,ord-1-iw
          
            if( tf )then
              nTet=nTet+1
              tetra(1:5,nTet)=[ pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
              &                 pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),&
              &                 pyramidIdx(ord=ord,iu=iu-2,iv=iv  ,iw=iw+1),&
              &                 pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
              &                 0                                           ]
            endif
          
          
            nPyr=nPyr+1
            pyram(1:6,nPyr)=[ pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),& 
            &                 pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
            &                 pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
            &                 pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
            &                 pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
            &                 0                                           ]
           !print '(3x,"nPyr=",i3,1x,"conec=",5(i3,1x))',nPyr,conec(1:5,nPyr)
            
            if( .not.iw==0 )then
              nPyr=nPyr+1
              pyram(1:6,nPyr)=[ pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),& 
              &                 pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
              &                 pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
              &                 pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
              &                 pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw-1),&
              &                 0                                           ]
            endif
            
            if( .not.tf .and. .not.iu==ord-1-iw )then
              nTet=nTet+1
              tetra(1:5,nTet)=[ pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
              &                 pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
              &                 pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
              &                 pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw+1),&
              &                 0                                           ]
            endif
            
          enddo
          tf=.not.tf
        enddo
      enddo
      
      open(unit=10,file="Pyramid.mesh",action='write')
      write(10,'( "MeshVersionFormatted 1")' )
      write(10,'(/"Dimension")' )
      write(10,'( "3")' )
      write(10,'(/"Vertices")' )
      write(10,'(i10)')nNod
      do iNod=1,nNod
        write(10,'(3(e22.15,1x),i3)')uvw(1:3,iNod),0
      enddo
      write(10,'(/"Pyramids")' )
      write(10,'(i10)')nPyr
      do iCel=1,nPyr
        write(10,'(5(i6,1x),i3)')pyram(1:6,iCel)
      enddo
      write(10,'(/"Tetrahedra")' )
      write(10,'(i10)')nTet
      do iCel=1,nTet
        write(10,'(5(i6,1x),i3)')tetra(1:5,iCel)
      enddo
      write(10,'(/"End")')
      close(10)
       
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return    
  end subroutine pyramidEquiNodes3D
  
  function pyramidIdx(ord, iu,iv,iw) result(idx)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer :: ord
    integer :: iu,iv,iw
    integer :: idx
    !>
    integer :: k
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    idx=0
    do k=0,iw-1
      idx=idx+(ord+1-k)**2
    enddo
    !print '("iu,iv,iw=",3(i3,2x),"idx=",i3)',iu,iv,iw,idx
    idx=iu+iv*(ord+1-iw)+idx +1
    
    return
  end function pyramidIdx 
  
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