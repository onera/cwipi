
subroutine pyramBasis()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use basePyramid
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer              :: i,j,nPt,cpt
  integer              :: ord,nMod
  real(8), pointer     :: uv (:,:)
  real(8), pointer     :: uvw(:,:)
  integer              :: ad,iSide
  integer              :: iNod,nNod,iu,iv,iw
  integer, allocatable :: conec(:,:)
  integer, allocatable :: idx(:)
  real(8)              :: rot(3,3),xyz(1:3),cos_a,sin_a
  real(8)              :: alpha
  
  real(8), parameter   :: eps=1d-12
  real(8), pointer     :: a(:),b(:),c(:)
  real(8), pointer     :: vand(:,:),dVand(:,:) !,jf(:,:),dr(:,:)
  real(8), pointer     :: duPsi   (:,:),dvPsi  (:,:),dwPsi  (:,:)
  real(8), pointer     :: drMatrix(:,:),dsMatrix(:,:),dtMatrix(:,:)
  real(8), pointer     :: mass(:,:)
  real(8), pointer     :: xyzOut(:,:),lxOut(:,:),drLxOut(:,:),dsLxOut(:,:),dtLxOut(:,:),leb(:,:)
  real(8), pointer     :: mode(:,:)
  real(8), pointer     :: eigv(:)
  real(8), pointer     :: fi(:)
  real(8)              :: f0
  
  character(3)         :: sfx
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Calcul des bases polynômiales Pyram")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ord=10 !25
  !write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !call pyramidNodes   (ord=ord, uvw=uvw, display=.true.)  !> Points réguliers
  call pyramidNodesOpt(ord=ord, uvw=uvw, display=.true.)  !> Points optimises
 !call writeMesh3D    (ord=ord, uvw=uvw)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramidSkin3D(ord=ord, uvw=uvw, display=.true.)
  call pyramidMesh3D(ord=ord, uvw=uvw, display=.true.)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
 !write(*,'(/"Pyramid (abc):")')
 !print '("a,b,c(",i2,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,a(ad),b(ad),c(ad),ad=1,size(a))
  call pyramidVandermonde3D(ord=ord,a=a,b=b,c=c,vand=vand)
  if( ord<3 )then
    call display(title="Vandermonde Matrix",mat=vand)
  else
    call displaySparce(title="Test avec l(uvw)",mat=lxOut,tol=1d-14)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Polynomes d'interpolation (on teste avec les points d'interpolation) => Matrice Identité
  
  nMod=(ord+1)*(ord+2)*(2*ord+3)/6
  nNod=size(a)
  
  allocate(lxOut(1:nMod,1:nNod))
  call pyramidLagrange3Dv(ord=ord,vand=vand,a=a,b=b,c=c,lx=lxOut,transpose=.true.) !> true pour affichage
  if( ord<3 )then
    call display(title="Test avec l(uvw)",mat=lxOut)
  else
    call displaySparce(title="Test avec l(uvw)",mat=lxOut,tol=1d-12)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Matrice de masse
  call massMatrix(vand=vand,mass=mass)
  if( ord<3 )then
    call display(title="Mass Matrix",mat=mass)
  endif
 !print '(/"det(mass)=",e22.15/)',DMGT(eps=eps,n=size(mass,1),A=mass)
 !call eigenVectors(mat=mass,w=eigv,display=.true.)
 !deallocate(eigv)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Derivees base fonctionnelle
  call pyramidGradVandermonde3D(ord=ord,a=a,b=b,c=c,drMode=duPsi,dsMode=dvPsi,dtMode=dwPsi)
  call derive1D(vand=vand,dVand=duPsi,dMat=drMatrix) !> drMatrix = duPsi.Inverse[vand]
  call derive1D(vand=vand,dVand=dvPsi,dMat=dsMatrix) !> dsMatrix = dvPsi.Inverse[vand]
  call derive1D(vand=vand,dVand=dwPsi,dMat=dtMatrix) !> dtMatrix = dwPsi.Inverse[vand]
  if( ord<3 )then
    !call display(title="duPsi Matrix",mat=duPsi)
    !call display(title="dvPsi Matrix",mat=dvPsi)
    !call display(title="dwPsi Matrix",mat=dwPsi)
    call display(title="drMatrix=duPsi.vand^{-1}",mat=drMatrix)
    call display(title="dsMatrix=dvPsi.vand^{-1}",mat=dsMatrix)
    call display(title="dsMatrix=dwPsi.vand^{-1}",mat=dtMatrix)
  endif
  deallocate(duPsi,dvPsi,dwPsi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(a,b,c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramBasis

subroutine pyramLebesgue()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use basePyramid
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer              :: i,j,order,np,nPt,cpt
  integer              :: ord
  real(8), pointer     :: uv (:,:)
  real(8), pointer     :: uvw(:,:)
  integer              :: ad,iSide
  integer              :: iNod,nNod,iu,iv,iw
  integer, allocatable :: conec(:,:)
  integer, allocatable :: idx(:)
  real(8)              :: rot(3,3),xyz(1:3),cos_a,sin_a
  real(8)              :: alpha
  !>
  real(8), parameter   :: eps=1d-12
  real(8), pointer     :: a(:),b(:),c(:)
  real(8), pointer     :: vand(:,:),dVand(:,:) !,jf(:,:),dr(:,:)
  real(8), pointer     :: duPsi  (:,:),dvPsi  (:,:),dwPsi  (:,:)
  real(8), pointer     :: drMatrix(:,:),dsMatrix(:,:),dtMatrix(:,:)
  real(8), pointer     :: mass(:,:)
  real(8), pointer     :: xyzOut(:,:),lxOut(:,:),drLxOut(:,:),dsLxOut(:,:),dtLxOut(:,:),leb(:,:)
  real(8), pointer     :: mode(:,:)
  real(8), pointer     :: eigv(:)
  real(8), pointer     :: fi(:)
  real(8)              :: f0
  
  character(3)         :: sfx
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'("Calcul fonction de Lebesgue pour pyramides")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !call pyramidNodes   (ord=ord, uvw=uvw, display=.false.)  !> Points réguliers
  call pyramidNodesOpt(ord=ord, uvw=uvw, display=.false.)  !> Points optimises
 !call writeMesh3D    (ord=ord, uvw=uvw)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call pyramiduvw2abc(uvw=uvw,a=a,b=b,c=c)
 !write(*,'(/"Pyramid (abc):")')
 !print '("a,b,c(",i2,")=",f12.5,2x,f12.5,2x,f12.5)',(ad,a(ad),b(ad),c(ad),ad=1,size(a))
  call pyramidVandermonde3D(ord=ord,a=a,b=b,c=c,vand=vand)
  if( ord<3 )then
    call display(title="Vandermonde Matrix",mat=vand)
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Derivees base fonctionnelle
  call pyramidGradVandermonde3D(ord=ord,a=a,b=b,c=c,drMode=duPsi,dsMode=dvPsi,dtMode=dwPsi)
  call derive1D(vand=vand,dVand=duPsi,dMat=drMatrix) !> drMatrix = duPsi.Inverse[vand]
  call derive1D(vand=vand,dVand=dvPsi,dMat=dsMatrix) !> dsMatrix = dvPsi.Inverse[vand]
  call derive1D(vand=vand,dVand=dwPsi,dMat=dtMatrix) !> dtMatrix = dwPsi.Inverse[vand]
  if( ord<3 )then
    call display(title="duPsi Matrix",mat=duPsi)
    call display(title="dvPsi Matrix",mat=dvPsi)
    call display(title="dwPsi Matrix",mat=dwPsi)
    call display(title="drMatrix=duPsi.vand^{-1}",mat=drMatrix)
    call display(title="dsMatrix=dvPsi.vand^{-1}",mat=dsMatrix)
    call display(title="dsMatrix=dwPsi.vand^{-1}",mat=dtMatrix)
  endif
  deallocate(duPsi,dvPsi,dwPsi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(uvw)
  deallocate(a,b,c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(   1<=ord .and. ord<  10 ) write(sfx,'("00",i1)')ord
  if(  10<=ord .and. ord< 100 ) write(sfx,'("0" ,i2)')ord
  if( 100<=ord .and. ord<1000 ) write(sfx,'(     i3)')ord
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Liste des points xyzOut  
  call pyramidReadXYZout3D(xyzOut=xyzOut, display=.true.)
  call pyramiduvw2abc(uvw=xyzOut,a=a,b=b,c=c)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des fonctions de Lagrange aux points xyzOut
  call pyramidBasePi(ord=ord,a=a,b=b,c=c,mode=mode,transpose=.false.)                !> Psi(xyzOut)
  call pyramidLagrange3Dv(ord=ord,vand=vand,a=a,b=b,c=c,lx=lxOut,transpose=.false.)  !> lxOut= Inverse[Transpose[Vand]].Psi[xyzOut] lxOut(nPt,np)
  if( ord<3 )then
    call pyramidWriteSolOut3D(title="simplex3D"//sfx,solOut=mode   )
    call pyramidWriteSolOut3D(title="lagrange3D"//sfx,solOut=lxOut  )
  endif
  deallocate(mode)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call lebesgue(lx=lxout,l=leb,transpose=.false.) !> leb = Sum abs(lxOut)
  call pyramidWriteSolOut3D(title="lebesgue3DP"//sfx,solOut=leb)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Evaluation des dérivées des fonctions de Lagrange aux points xyzOut
  call dLagrange1Dv(dMat=drMatrix,lx=lxOut,dlx=drLxOut,transpose=.false.) !> drLxOut(1:nPt,1:np)= Transpose[drMatrix] lxOut
  call dLagrange1Dv(dMat=dsMatrix,lx=lxOut,dlx=dsLxOut,transpose=.false.) !> dsLxOut= Transpose[dsMatrix] lxOut
  call dLagrange1Dv(dMat=dtMatrix,lx=lxOut,dlx=dtLxOut,transpose=.false.) !> dtLxOut= Transpose[dtMatrix] lxOut
  
  if( ord<3 )then
    call pyramidWriteSolOut3D(title="drLagrange3D"//sfx,solOut=drLxOut)
    call pyramidWriteSolOut3D(title="dsLagrange3D"//sfx,solOut=dsLxOut)
    call pyramidWriteSolOut3D(title="dtLagrange3D"//sfx,solOut=dtLxOut)
  endif  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions : preliminaires
  print '("size(lxOut)=",i6," x ",i3)',size(lxOut,1),size(lxOut,2)
  nPt=size(lxOut,1)
  np =(ord+1)*(ord+2)*(2*ord+3)/6 !> = \sum_{k=1}^{ord+1} k^2  
  allocate(fi(1:np)) ; fi(1:np)=1d0
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Test des fonctions lxOut : f(xyzOut) = \sum_{i=1,np} lxOut(xyzOut,i) f_i
  print '("Test des fonctions lxOut : f(xyzOut) = \sum_{i=1,np} lxOut(xyzOut,i) f_i")'
  cPt=0
  do i=1,nPt
    f0=0d0
    do j=1,np
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
  !> Test des fonctions drLxOut : df/dr(xyzOut) = \sum_{i=1,np} drLxOut(xyzOut,i) f_i
  print '("Test des fonctions drLxOut : df/dr(xyzOut) = \sum_{i=1,np} drLxOut(xyzOut,i) f_i")'
  cPt=0
  do i=1,nPt
    f0=0d0
    do j=1,np
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
  !> Test des fonctions dsLxOut : df/ds(xyzOut) = \sum_{i=1,np} dsLxOut(xyzOut,i) f_i
  print '("Test des fonctions dsLxOut : df/ds(xyzOut) = \sum_{i=1,np} dsLxOut(xyzOut,i) f_i")'
  cPt=0
  do i=1,nPt
    f0=0d0
    do j=1,np
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
  !> Test des fonctions dtLxOut : df/dt(xyzOut) = \sum_{i=1,np} dtLxOut(xyzOut,i) f_i
  print '("Test des fonctions dtLxOut : df/dt(xyzOut) = \sum_{i=1,np} dtLxOut(xyzOut,i) f_i")'
  cPt=0
  do i=1,nPt
    f0=0d0
    do j=1,np
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
end subroutine pyramLebesgue

subroutine pyramMaillageVisu()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !use modDeterminant
  !use baseSimplex2D
  !use baseSimplex3D
  use basePyramid
  use table_tet_mesh,only: driverTetMesh,saveTetMesh
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer            :: ord,iOrd,ad
  real(8), pointer   :: uvw  (:,:)
  integer, pointer   :: tetra(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'(/"Construction maillage Pyramid P_i")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)ord
  ord=10
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do iOrd=1,ord
 !do iOrd=ord,ord
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call pyramidNodes   (ord=iOrd, uvw=uvw, display=.false.)
    call driverTetMesh  (ord=iOrd, node_xyz=uvw,tetra_node=tetra)
    deallocate(uvw)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call pyramidNodesOpt(ord=iOrd, uvw=uvw, display=.true.)  !> Points optimises
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call saveTetMesh(ord=iOrd, node_xyz=uvw,tetra_node=tetra)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !call pyramidSkin3D(ord=iOrd, uvw=uvw, display=.true.)
    !call pyramidMesh3D(ord=iOrd, uvw=uvw, display=.true.)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(uvw,tetra)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    print '("nn=",i10," ne=",i10," Pi=",i3)',size(uvw,2),size(tetra,2),iOrd
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine pyramMaillageVisu



program main
  
  !> Test pyramids
  print '("coucou0")'
  call pyramBasis()
  print '("coucou1")'
  call pyramLebesgue()
  print '("coucou2")'
  call pyramMaillageVisu() !> maillages de visu pour la pyramide d'ordre élevé
  print '("coucou3")'
  !call pyramDegreesOverSides()
  
end program main