module basePyramid
  
  implicit none

contains
  
  
  subroutine pyramidNodes(ord, uvw, display)
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
    if( display )print '(/"Building Pyramid Equidistant Nodes")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nNod=0
    do iu=1,ord+1
      nNod=nNod+iu*iu
    enddo
    if( display )print '(3x,"nDeg=",i6)',nNod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nPyr=0
    do iw=ord,1,-1
      nPyr=nPyr+(iw*iw)+((iw-1)*(iw-1)) !> Pyramides droites + pyramides retournees
    enddo
    if( display )print '(3x,"nPyr=",i6)',nPyr
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
      print '(3x,"Vertices")'
      do iNod=1,nNod
        print '(6x,"uvw(",i6,")=",3(f12.9,1x))',iNod,uvw(1:3,iNod)
      enddo
      print '(3x,"end Vertices")'
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    mesh=.false.
    if( mesh )then
      allocate(pyram(6,nPyr))
      allocate(tetra(6,1000))
      nPyr=0
      nTet=0
      tf=.false.
      do iw=0,ord-1
        do iv=0,ord-iw-1
          print '("tf=",l)',tf
          do iu=0,ord-iw-1
          
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
            
            if( .not.iw==0 )then
              nPyr=nPyr+1
              pyram(1:6,nPyr)=[ pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),& 
              &                 pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
              &                 pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
              &                 pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
              &                 pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw-1),&
              &                 0                                           ]
            endif
            
            if( .not.tf .and. .not.iu==ord-iw-1 )then
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
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '("end Building Pyramid Equidistant Nodes")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidNodes
  
  
  subroutine pyramidNodesOpt(ord, uvw, uv, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! input: ord=polynomial order of interpolant
    ! output: uvw(:,:) node coordinates in unity pyramid
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(in) , pointer :: uv (:,:) !> Triangle optimized points
    real(8), intent(out), pointer :: uvw(:,:) !> Tetra    optimized points
    logical, intent(in)           :: display
    !>
    integer                       :: iu,iv,iw
    integer                       :: iNod,jNod,nNod
    real(8)                       :: subV(0:ord)
    real(8)                       :: subW(0:ord)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(/"Building Pyramid Optimized Nodes")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      print '(3x,"3D Triangle Optimized Nodes")'
      iNod=0
      do iw=0,ord
        print '(6x,"level:",i3)',iw
        do iv=0,ord-iw
          iNod=iNod+1
          print '(9x,"uv(",i6,")=",2(f12.9,1x))',iNod,uv(1:2,iNod)
        enddo
      enddo
      print '(3x,"end")'
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nNod=0
    do iu=1,ord+1
      nNod=nNod+iu*iu
    enddo
    if( display )print '(3x,"nDeg=",i6)',nNod
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(uvw(3,nNod))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iNod=0 ; jNod=0
    
   !subW(0:ord)=uv(1,1:ord+1)
    do iw=0,ord
      
      subV(0:ord-iw)=2d0*uv(1,jNod +1:jNod+ord-iw +1)-uv(1,jNod+ord-iw +1)
      subW(0:ord-iw)=    uv(2,jNod +1:jNod+ord-iw +1)
      
      do iv=0,ord-iw
        jNod=jNod+1
        do iu=0,ord-iw
        iNod=iNod+1 !> tensorisation
        uvw(1:3,iNod)=[subV(iu),subV(iv),  (subW(iu)+subW(iv))/2 ]
       !uvw(1:3,iNod)=[subV(iu),subV(iv),subW(iw)]
      enddo ; enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )then
      print '(3x,"Vertices")'
      iNod=0
      do iw=0,ord
        print '(6x,"level: ",i3)',iw
        do iv=0,ord-iw ; do iu=0,ord-iw
          iNod=iNod+1
          print '(6x,"uvw(",i6,")=",3(f12.9,1x))',iNod,uvw(1:3,iNod)
        enddo ; enddo
      enddo
      print '(3x,"end")'
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '("end Building Pyramid Optimized Nodes")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidNodesOpt
  
  subroutine pyramidMesh3D(ord,uvw,display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use M_libmesh6_api
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:,:)
    logical, intent(in)           :: display
    !>
    integer                       :: iu,iv,iw
    integer                       :: iNod ,nNod, iCel
    logical                       :: test
    real(4)                       :: dist
    character(3)                  :: sfx
    !> libmesh
    character(256)                :: name
    
    integer                       :: ins,ver,res,geo
    integer, allocatable          :: TypTab(:)
    real(4)                       :: xyz(3)
    integer                       :: nFld,kind(1)
    integer                       :: iNod0,nNod0
    real(8), allocatable          :: uvw0(:,:)
    integer, allocatable          :: mark(:)
    integer                       :: nTetr
    integer, allocatable          :: tetr(:,:)
    integer                       :: nTria
    integer, allocatable          :: tria(:,:)
    !> Reordering
    real(8), parameter            :: eps=1d-6
    integer , allocatable         :: indx(:)
    integer                       :: iNod0Min,iNodMin
    real(8)                       :: d2,d2Min
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(/"Building Pyramid Volumic Mesh")'
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
      name="nodes3DP"//sfx//".mesh" ; if( display )print '(3x,"Building and Writing: ",a)',trim(name)
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
      name="nodes3DP"//sfx//".sol" ;  if( display )print '(3x,"Building and Writing: ",a)',trim(name)
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
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(3x,"Building Pyramid Volumic Mesh with ghs3d: ",a)',"PyramidP"//sfx//".mesh"
    if( nNod==0 )then
      call system("ghs3d -O 1  -exit 3 -in PyramidSkinP"//sfx//".mesh                          -out PyramidP"//sfx//".mesh > ghs3d.log")
    else
      call system("ghs3d -O 1  -exit 3 -in PyramidSkinP"//sfx//".mesh -force nodes3DP"//sfx//" -out PyramidP"//sfx//".mesh > ghs3d.log")
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(3x,"rm -f nodes3DP",a,".mesh")',sfx
    if( display )print '(3x,"rm -f nodes3DP",a,".sol")',sfx
    if( display )print '(3x,"rm -f PyramidSkinP",a,".mesh")',sfx
    call system("rm -f nodes3DP"//sfx//".mesh")
    call system("rm -f nodes3DP"//sfx//".sol")
    call system("rm -f PyramidSkinP"//sfx//".mesh")
    call system("rm -f PyramidP"//sfx//".sol")
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    name="PyramidP"//sfx//".mesh"
    if( display )print '(3x,"Reading Pyramid Volumic Mesh: ",a)',trim(name)
    ins=GmfOpenMeshF77(trim(name),GmfRead,ver,geo) ! print '(3x,"ins=",i3)',ins
    nNod0=GmfStatKwdF77(ins,GmfVertices,ver,0,TypTab) ; if( display )print '(6x,"nNod0=",i10)',nNod0
    allocate(uvw0(geo+1,nNod0),mark(nNod0))  ! il faut mettre geo+1 pour ensuite reoordonner les noeuds
    res=GmfGotoKwdF77(ins,GmfVertices)
    select case(ver)
    case(1) !> real(4)
      do iNod0=1,nNod0
        call GmfGetVertex3dr4(ins,xyz(1),xyz(2),xyz(3),mark(iNod0))
        uvw0(1:3,iNod0)=xyz(1:3)
        uvw0(4,iNod0)=1d0-uvw0(1,iNod0)-uvw0(2,iNod0)-uvw0(3,iNod0)
      enddo
    case(2) !> real(8)
      do iNod0=1,nNod0
        call GmfGetVertex3dr8(ins,uvw0(1,iNod0),uvw0(2,iNod0),uvw0(3,iNod0),mark(iNod0))
        uvw0(4,iNod0)=1d0-uvw0(1,iNod0)-uvw0(2,iNod0)-uvw0(3,iNod0)
      enddo
    end select
    
    nTetr=GmfStatKwdF77(ins,GmfTetrahedra,0,0,TypTab) ; if( display )print '(6x,"nTetr=",i10)',nTetr
    allocate(tetr(5,nTetr))
    res=GmfGotoKwdF77(ins,GmfTetrahedra)
    do iCel=1,nTetr
      call GmfGetTetrahedron(ins,tetr(1,iCel),tetr(2,iCel),tetr(3,iCel),tetr(4,iCel),tetr(5,iCel))
    enddo
    
    nTria=GmfStatKwdF77(ins,GmfTriangles,0,0,TypTab) ; if( display )print '(6x,"nTria=",i10)',nTria
    allocate(tria(4,nTria))
    res=GmfGotoKwdF77(ins,GmfTriangles)
    do iCel=1,nTria
      call GmfGetTriangle(ins,tria(1,iCel),tria(2,iCel),tria(3,iCel),tria(4,iCel))
    enddo
    
    res=GmfCloseMeshF77(ins)
    if( display )print '(3x,"end Reading")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(3x,"Reordering Nodes")'
    
    nNod =size(uvw ,2)
    nNod0=size(uvw0,2)
    
    allocate(indx(nNod0)) ; indx(1:nNod0)=0
    do iNod=1,nNod
      iNod0=1 ; d2Min=1d50 ; iNodMin=1
      loop1 : do
        
        d2= (uvw0(1,iNod0)-uvw(1,iNod))**2 &
        &  +(uvw0(2,iNod0)-uvw(2,iNod))**2 &
        &  +(uvw0(3,iNod0)-uvw(3,iNod))**2
        
        if( d2<d2Min )then
          iNod0Min=iNod0 ; d2Min=d2
        endif
        
        if( d2<eps )then
          indx(iNod0)=iNod
          exit loop1
        endif
        
        iNod0=iNod0+1
        if( iNod0>nNod0 )then
          print '("uvw (",i6,")=",3(f12.5,1x),1x,"d2Min=",e12.5)',iNod    ,uvw (1:3,iNod    ),d2Min
          print '("uvw0(",i6,")=",3(f12.5,1x),1x,"d2Min=",e12.5)',iNod0Min,uvw0(1:3,iNod0Min),d2Min
          stop"problem @ writeMeshSkin3D"
        endif
      enddo loop1
    enddo
    
    if( display )print '(3x,"end Reordering")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Test pour verifier que tous les noeuds sont en correspondance
    if( .not.count( indx(:)==0 )==0 )then
      print '("count(indx==0)=",i10)',count(indx(:)==0)
      stop "problem @ pyramidMesh3D"
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(3x,"rm -f PyramidP",a,".mesh")',sfx
    call system("rm -f PyramidP"//sfx//".mesh")
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Reecriture de TetraPi.mesh
    name="PyramidP"//sfx//".mesh" ; if( display )print '(3x,"Writing Reordered Pyramid Volumic Mesh: ",a)',trim(name)
    ins=GmfOpenMeshF77(trim(name),GmfWrite,1,3) ; if( display )print '(6x,"nNod =",i10)',nNod
    
    res=GmfSetKwdF77(ins,GmfVertices,nNod,0,TypTab)
    do iNod=1,nNod
      call GmfSetVertex3dr4(ins                     ,&
      &                     real(uvw(1,iNod),kind=4),&
      &                     real(uvw(2,iNod),kind=4),&
      &                     real(uvw(3,iNod),kind=4),&
      &                     i3=0                     )
    enddo
    
    res=GmfSetKwdF77(ins,GmfTetrahedra,nTetr,0,TypTab) ; if( display )print '(6x,"nTetr=",i10)',nTetr
    do iCel=1,nTetr
      call GmfSetTetrahedron(ins               ,&
      &                      indx(tetr(1,iCel)),&
      &                      indx(tetr(2,iCel)),&
      &                      indx(tetr(3,iCel)),&
      &                      indx(tetr(4,iCel)),&
      &                           tetr(5,iCel)  )
    enddo
    
    res=GmfSetKwdF77(ins,GmfTriangles,nTria,0,TypTab) ;  if( display )print '(6x,"nTria=",i10)',nTria
    do iCel=1,nTria
      call GmfSetTriangle(ins               ,&
      &                  indx(tria(1,iCel)),&
      &                  indx(tria(2,iCel)),&
      &                  indx(tria(3,iCel)),&
      &                       tria(4,iCel)  )
    enddo
    
    res=GmfCloseMeshF77(ins)
    if( display )print '(3x,"end")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '("end Building Pyramid Volumic Mesh")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidMesh3D
  
  subroutine pyramidSkin3D(ord, uvw, display)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Cette routine construit le maillage surfacique
    !> d'une pyramide Pord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer, intent(in)           :: ord
    real(8), intent(out), pointer :: uvw(:,:)
    logical, intent(in)           :: display
    !---
    integer                       :: iu,iv,iw
    integer                       :: iNod,nNod
    integer                       :: iCel,nCel
    integer                       :: nQuadr,nTrian
    integer, allocatable          :: quadr(:,:)
    integer, allocatable          :: trian(:,:)
    character(3)                  :: sfx
    character(256)                :: name
    
    logical, allocatable          :: keep(:)
    integer, allocatable          :: indx(:)
    integer                       :: iNod0,nNod0
    real(8), pointer              :: uvw0(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '(/"Building Pyramid Surfasic Mesh")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
    if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
    if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nNod=size(uvw,2)
    allocate(quadr(1:5,  ord*ord))
    allocate(trian(1:4,4*ord*ord))
    nQuadr=0
    nTrian=0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Side1
    if( display )print '(3x,"Side1")'
    iw=0
    do iv=0,ord-iw-1
      do iu=0,ord-iw-1
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
    if( display )print '(3x,"Side2")'
    do iw=0,ord-1
      do iv=0,ord-iw-1
        do iu=0,0
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
          &                  2                                           ]
          print '("Trian(",i3,")=",4(I3,1x))',nTrian,trian(1:4,nTrian)
          if( .not. iv==ord-iw-1 )then
            nTrian=nTrian+1
            trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),& 
            &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw+1),&
            &                  2                                           ]
            print '("Trian(",i3,")=",4(I3,1x))',nTrian,trian(1:4,nTrian)
          endif
        enddo
      enddo
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Side3
    if( display )print '(3x,"Side3")'
    do iw=0,ord-1
      do iv=0,0
        do iu=0,ord-iw-1
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw  ),& 
          &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  3                                           ]
          if( .not. iu==ord-iw-1 )then
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
    if( display )print '(3x,"Side4")'
    do iw=0,ord-1
      do iv=0,ord-iw-1
        do iu=ord-iw-1,ord-1-iw
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv  ,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  4                                           ]            
          if( .not.iv==ord-iw-1 )then
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
      do iv=ord-iw-1,ord-iw-1
        do iu=0,ord-iw-1
          nTrian=nTrian+1
          trian(1:4,nTrian)=[pyramidIdx(ord=ord,iu=iu+1,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv+1,iw=iw  ),&
          &                  pyramidIdx(ord=ord,iu=iu  ,iv=iv  ,iw=iw+1),&
          &                  5                                           ]
          if( .not.iu==ord-iw-1 )then
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
    !> Removing Extra vertices
    if( display ) print '(3x,"Removing Extra vertices")'
    
    allocate(keep(nNod)) ; keep(:)=.false.
    do iCel=1,nQuadr
      do iNod=1,4
        keep( quadr(iNod,iCel) )=.true.
      enddo
    enddo
    do iCel=1,nTrian
      do iNod=1,3
        keep( trian(iNod,iCel) )=.true.
      enddo
    enddo
    
    allocate(indx(nNod) )
    iNod0=0
    do iNod=1,nNod
      if( keep(iNod) )then
        iNod0=iNod0+1
        indx(iNod)=iNod0
      endif
    enddo
    nNod0=iNod0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Ecriture
    name="PyramidSkinP"//sfx//".mesh"       
    open(unit=10,file=trim(name),action='write')
    write(10,'( "MeshVersionFormatted 1")' )
    write(10,'(/"Dimension"/,"3")')
    write(10,'(/"Vertices"/,i10)' )nNod0
    do iNod=1,nNod
      if( keep(iNod) )write(10,'(3(e22.15,1x),i3)')uvw(1:3,iNod),0        
    enddo
    write(10,'(/"Quadrilaterals"/,i10)' )nQuadr
    do iCel=1,nQuadr
      write(10,'(4(i6,1x),i3)')indx(quadr(1:4,iCel)),quadr(5,iCel)
    enddo
    write(10,'(/"Triangles"/,i10)' )nTrian
    do iCel=1,nTrian
      write(10,'(3(i6,1x),i3)')indx(trian(1:3,iCel)),trian(4,iCel)
    enddo
    write(10,'(/"End")')
    close(10)
    
    if( display )then
      print '(3x,"Writing ",a)',trim(name)
      print '(6x,"Vertices:       ",i10)',nNod0
      print '(6x,"Quadrilaterals: ",i10)',nQuadr
      print '(6x,"Triangles:      ",i10)',nTrian
      print '(3x,"end Writing")'
    endif       
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display ) print '(3x,"Cleanning Memory")'
    deallocate(quadr,trian)
    deallocate(keep,indx)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '("end Building Pyramid Surfasic Mesh")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidSkin3D
  
  
  function pyramidIdx(ord, iu,iv,iw) result(idx)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer :: ord
    integer :: iu,iv,iw
    integer :: idx
    !>
    integer :: k
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    idx=0
    do k=0,iw-1
      idx=idx+(ord+1-k)**2
    enddo
    !print '("iu,iv,iw=",3(i3,2x),"idx=",i3)',iu,iv,iw,idx
    idx=iu+iv*(ord+1-iw)+idx +1
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
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
    if( display )print '(/"Extract Pyramid List of DOF per Side")'
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
      do iSide=1,5
        print '(3x,"Side",i1,": ",$)',iSide
        do iNod=sideIdx(iSide)+1,sideIdx(iSide+1)
          print '(i5,1x,$)',side(iNod)
        enddo
        print '()'
      enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( display )print '("end Extract Pyramid List of DOF per Side")'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine pyramidSides3D
  
  
end module basePyramid