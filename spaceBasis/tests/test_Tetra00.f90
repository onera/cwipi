

subroutine tetraTest()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use baseSimplex2D
  use baseSimplex3D
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer              :: i,j,order,nMod,nNod,ad,np,nPt,cpt
  real(8), pointer     :: vand(:,:),dVand(:,:),jf(:,:),dr(:,:)
  real(8), pointer     :: drVand  (:,:),dsVand  (:,:),dtVand  (:,:)
  real(8), pointer     :: drMatrix(:,:),dsMatrix(:,:),dtMatrix(:,:)
  real(8), pointer     :: mass(:,:)
  real(8), pointer     :: psi(:,:)
  real(8), pointer     :: xyzOut(:,:),lxOut(:,:),drLxOut(:,:),dsLxOut(:,:),dtLxOut(:,:),leb(:,:)
  real(8), pointer     :: uvw(:,:),rst(:,:),a(:),b(:),c(:)
  real(8), parameter   :: eps=1d-12
  real(8), pointer     :: eigv(:)
  real(8), pointer     :: fi(:)
  real(8)              :: f0
  character(3)         :: sfx
  real(8), pointer     :: xGLL(:)
  real(8), pointer     :: uv(:,:)
  integer, allocatable :: conec(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'("Calcul des bases polynômiales Tetra")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  order=10
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(   1<=order .and. order<  10 ) write(sfx,'("00",i1)')order
  if(  10<=order .and. order< 100 ) write(sfx,'("0" ,i2)')order
  if( 100<=order .and. order<1000 ) write(sfx,'(     i3)')order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Tétraèdre points d''interpolation")')
  call nodes3D   (ord=order,uvw=uvw,display=.false.)
  call nodes3Dopt(ord=order,uvw=uvw,display=.true. )
  
  call nodes3Duvw2rst(uvw=uvw,rst=rst) !> rst(1:3,:)=2d0*uvw(1:3,:)-1d0
 !write(*,'(/"Tetra (rst):")')
 !print '("rst(1:3,",i2,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,rst(1:3,ad),ad=1,size(rst,2))
  
  call nodes3Drst2abc(rst=rst,a=a,b=b,c=c)
 !write(*,'(/"Tetra (abc):")')
 !print '("a,b,c(",i2,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,a(ad),b(ad),c(ad),ad=1,size(a))
  
  nNod=size(uvw,2)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Extraction des faces et vérification des coordonnées barycentriques
  if( 0==1 )then
    call gaussLegendreLobatto(ord=order,xGLL=xGLL)
    print '(/"Point de Gauss Lobatto")'
    print '("ad=",i5,2x,"u=",f19.16)',(ad,5d-1*(xGLL(ad)+1d0),ad=1,size(xGLL))
    
    call nodes2D   (ord=order,uvw=uv,display=.false.)
    call nodes2Dopt(ord=order,uvw=uv,display=.true. )
    
    call trianglesConnectivity(ord=order,conec=conec)
    print '()'
   !do i=1,4
    do i=3,3
      print '("Triangle",i1)',i
      do j=1,size(conec,1)
        ad=conec(j,i)
        !print '("dg=",i6," uvw=",3(f19.16,1x))',conec(j,i),uvw(1:3,conec(j,i))
        !print '("dg=",i5,2x,"u=",f19.16,2x,"v=",f19.16,2x,"w=",f19.16)',ad,uvw(1:3,ad)
        print '("dg=",i5,2x,"du=",f19.16,2x,"dw=",f19.16)',ad,uvw(1,ad)-uv(1,j),uvw(3,ad)-uv(2,j)
      enddo
    enddo
    
    do j=1,size(conec,1)
      ad=conec(j,3)
      uv(1,j)=uvw(1,ad)
      uv(2,j)=uvw(3,ad)
      uv(3,j)=1d0-uv(1,j)-uv(2,j)
    enddo
    deallocate(conec,uv)
    
    call nodes3Dopt_2D(ord=order,uvw=uv,display=.true.)
    stop
    
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call writeMesh3D(ord=order,uvw=uvw)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call vandermonde3D(ord=order,a=a,b=b,c=c,vand=vand)
  if( order<3 )then
    call display(title="Vandermonde Matrix",mat=vand)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Polynomes d'interpolation (on teste avec les points d'interpolation) => Matrice Identité
  nMod=(order+1)*(order+2)*(order+3)/6
  allocate(lxOut(1:nMod,1:nNod))                  !> allocation pour true
  call lagrange3Dv(ord=order,vand=vand,a=a,b=b,c=c,lx=lxOut,transpose=.true.) !> true pour affichage
  if( order<3 )then
    call display(title="Test avec l(uvw)",mat=lxOut)
  else
    call displaySparce(title="Test avec l(uvw)",mat=lxOut,tol=1d-14)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse
  call massMatrix(vand=vand,mass=mass)
  if( order<3 )then
    call display(title="Mass Matrix",mat=mass)
  endif
  print '(/"det(mass)=",e22.15/)',DMGT(eps=eps,n=size(mass,1),A=mass)
  call eigenVectors(mat=mass,w=eigv,display=.true.)
  deallocate(eigv)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Derivees base fonctionnelle
  call gradVandermonde3D(ord=order,a=a,b=b,c=c,drVand=drVand,dsVand=dsVand,dtVand=dtVand)
  call derive1D(vand=vand,dVand=drVand,dMat=drMatrix) !> drMatrix = drVand.Inverse[vand]
  call derive1D(vand=vand,dVand=dsVand,dMat=dsMatrix) !> dsMatrix = dsVand.Inverse[vand]
  call derive1D(vand=vand,dVand=dtVand,dMat=dtMatrix) !> dtMatrix = dtVand.Inverse[vand]
  if( order<3 )then
    call display(title="drVand Matrix",mat=drVand)
    call display(title="dsVand Matrix",mat=dsVand)
    call display(title="dtVand Matrix",mat=dtVand)
    call display(title="drMatrix=drVand.vand^{-1}",mat=drMatrix)
    call display(title="dsMatrix=dsVand.vand^{-1}",mat=dsMatrix)
    call display(title="dsMatrix=dtVand.vand^{-1}",mat=dtMatrix)
  endif
  deallocate(drVand,dsVand,dtVand)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(a,b,c,rst)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des fonctions de Lagrange aux points xyzOut
  call readXYZout3D(xyzOut=xyzOut)
  
  nMod=(order+1)*(order+2)*(order+3)/6 ; nNod=size(xyzOut,2)
  
  call nodes3Drst2abc(rst=xyzOut,a=a,b=b,c=c)
  call simplex3D  (ord=order,a=a,b=b,c=c,mode=psi,transpose=.false.)           !> Psi(xyzOut)
  
  
  allocate(lxOut(1:nNod,1:nMod))
  call lagrange3Dv(ord=order,vand=vand,a=a,b=b,c=c,lx=lxOut,transpose=.false.)  !> lxOut= Inverse[Transpose[Vand]].Psi[xyzOut] lxOut(nPt,np)
  if( order<3 )then
    call writeSolOut3D(title="simplex3D"//sfx,solOut=psi)
    call writeSolOut3D(title="lagrange3D"//sfx  ,solOut=lxOut  )
  endif
  deallocate(psi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des dérivées des fonctions de Lagrange aux points xyzOut
  allocate(drLxOut(1:nNod,1:nMod))
  call dLagrange1Dv(dMat=drMatrix,lx=lxOut,dlx=drLxOut,transpose=.false.) !> drLxOut(1:nPt,1:np)= Transpose[drMatrix] lxOut
  allocate(dsLxOut(1:nNod,1:nMod))
  call dLagrange1Dv(dMat=dsMatrix,lx=lxOut,dlx=dsLxOut,transpose=.false.) !> dsLxOut= Transpose[dsMatrix] lxOut
  allocate(dtLxOut(1:nNod,1:nMod))
  call dLagrange1Dv(dMat=dtMatrix,lx=lxOut,dlx=dtLxOut,transpose=.false.) !> dtLxOut= Transpose[dtMatrix] lxOut
  
!  drLxOut=2d0*drLxOut ! car u'=2u-1
!  dsLxOut=2d0*dsLxOut ! car v'=2v-1
!  dtLxOut=2d0*dtLxOut ! car w'=2w-1
  
  if( order<3 )then
    call writeSolOut3D(title="drLagrange3D"//sfx,solOut=drLxOut)
    call writeSolOut3D(title="dsLagrange3D"//sfx,solOut=dsLxOut)
    call writeSolOut3D(title="dtLagrange3D"//sfx,solOut=dtLxOut)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions : preliminaires
  print '("size(lxOut)=",i6," x ",i3)',size(lxOut,1),size(lxOut,2)
  allocate(fi(1:nMod)) ; fi(1:nMod)=1d0
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions lxOut : f(xyzOut) = \sum_{i=1,nMod} lxOut(xyzOut,i) f_i
  cPt=0
  do i=1,nNod
    f0=0d0
    do j=1,nMod
      f0=f0+lxOut(i,j)*fi(j)
    enddo
    if( abs(f0-1d0)>eps )then
      print '("f(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(f)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions drLxOut : df/dr(xyzOut) = \sum_{i=1,nMod} drLxOut(xyzOut,i) f_i
  cPt=0
  do i=1,nNod
    f0=0d0
    do j=1,nMod
      f0=f0+drLxOut(i,j)*fi(j)
    enddo
    if( abs(f0)>eps )then
      print '("df/dr(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(df/dr)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions dsLxOut : df/ds(xyzOut) = \sum_{i=1,nMod} dsLxOut(xyzOut,i) f_i
  cPt=0
  do i=1,nNod
    f0=0d0
    do j=1,nMod
      f0=f0+dsLxOut(i,j)*fi(j)
    enddo
    if( abs(f0)>eps )then
      print '("df/ds(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(df/ds)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions dtLxOut : df/dt(xyzOut) = \sum_{i=1,nMod} dtLxOut(xyzOut,i) f_i
  cPt=0
  do i=1,nNod
    f0=0d0
    do j=1,nMod
      f0=f0+dtLxOut(i,j)*fi(j)
    enddo
    if( abs(f0)>eps )then
      print '("df/dt(",i6,")=",e22.15)',i,(f0-1d0)
      cpt=cpt+1
    endif
  enddo
  if( .not.cpt==0 )print '("cpt(df/dt)=",i6)',cpt
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(fi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  allocate(leb(1:nNod,1))
  call lebesgue(lx=lxout,l=leb,transpose=.false.)
  call writeSolOut3D(title="lebesgue3DP"//sfx,solOut=leb)
  
  deallocate(leb)
  deallocate(xyzOut)
  deallocate(a,b,c)
  !
  deallocate(lxOut)
  deallocate(drLxOut,dsLxOut,dtLxOut)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(vand)
  deallocate(drMatrix,dsMatrix,dtMatrix)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine tetraTest

subroutine tetraMaillageVisu()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Cette procedure sert à construire les maillages de visu pour le tetra d'ordre élevé
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use baseSimplex3D
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer            :: order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Construction maillage Tetra P_i")')
  write(*,'("Warning ghs3d is required")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Building Mesh for Tetra P_{ord}
  call writeMeshSkin3D(ord=order)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine tetraMaillageVisu

subroutine tetraMaillageVisuNew()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Cette procedure sert à construire les maillages de visu pour le tetra d'ordre élevé
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !use modDeterminant
  use baseSimplex3D
  use table_tet_mesh
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer            :: ord,iOrd,ad,nVert
  real(8), pointer   :: uvw(:,:)
  real(8), pointer   :: uvw0(:,:)
  integer, pointer   :: tetra(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Construction maillage Tetra P_i")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !do iOrd=1,ord
  do iOrd=ord,ord
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call nodes3D(ord=iOrd,uvw=uvw,display=.true.)
    nVert=size(uvw,2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    allocate(uvw0(1:3,1:nVert))
    uvw0(1:3,:)=uvw(1:3,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call driverTetMesh(node_xyz=uvw0,tetra_node=tetra)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call nodes3Dopt(ord=iOrd,uvw=uvw,display=.true. )
    uvw0(1:3,:)=uvw(1:3,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call saveTetMesh(ord=iOrd, node_xyz=uvw0,tetra_node=tetra)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw,tetra,uvw0)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  return
end subroutine tetraMaillageVisuNew


program main
  
  !> Test Tetra
  call tetraTest()
 !call tetraMaillageVisu() !> maillages de visu pour le tetra d'ordre élevé
  call tetraMaillageVisuNew()
  
end program main