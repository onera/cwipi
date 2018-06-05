subroutine edge_02()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use modDeterminant
  use baseSimplex1D
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  !
  integer            :: order
  integer            :: iMode,nMode
  real(8)            :: coef0
  real(8), pointer   :: coef(:)
  integer            :: iGaus,nGaus
  real(8), pointer   :: xGL (:),wGL (:),modeGL (:,:)
  integer            :: iVert,nVert
  real(8), pointer   :: xOut(:),fxOut(:),rFxOut(:),rFxOutD(:)
  real(8), pointer   :: mode(:,:),modeD(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  write(*,'("Reconstruction 1D d''une fonction discontinue")')
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  order=11
  print '(/"order=",i2)',order
 !write(*,'(/"Order: ")',advance='no') ; read(*,*)order
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '(/"Gauss Legendre Quadatures")'
  call gaussLegendreQuadratures(ord=order,xGL=xGL,wGL=wGL)
  print '("i=",i2," xGL=",f12.5," wGL=",f12.5)',(iGaus,xGL(iGaus),wGL(iGaus),iGaus=1,size(xGL))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '(/"Valeurs des polynômes de Legendre en xGL")'
  nMode=order+1
  nGaus=size(xGL)
  allocate(modeGL(1:nMode,1:nGaus))
  call simplex1D(ord=order,r=xGL,mode=modeGL,transpose=.true.) !> mode(1:nMod,1:nGaus)
  call display(title="modeGL",mat=modeGL)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '(/"Projection de func sur la base des polynômes de Legendre")'
  
  allocate(coef(1:nMode))
  coef(1:nMode)=0d0
  
  do iGaus=1,nGaus
    coef0=wGL(iGaus)*func(xGL(iGaus))
    do iMode=1,nMode
      coef(iMode)=coef(iMode)+coef0*modeGL(iMode,iGaus)
    enddo
  enddo  
  print '(/"Liste des coefficients")'
  print '("i=",i2," coef=",f12.5)',(iMode,coef(iMode),iMode=1,nMode)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Points de reconstruction
  nVert=1001 ; allocate(xOut(1:nVert))
  xOut(1:nVert)=[( (-1d0+2d0*real(iVert-1,kind=8)/real(nVert-1,kind=8)), iVert=1,nVert)]
  print '(/"Liste des points de reconstruction")'
  print '("size(xOut)=",i4)',size(xOut)
 !print '("xOut(",i2,")=",f12.5)',(iVert,xOut(iVert),iVert=1,nVert)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Calcul func(xOut)
  allocate(fxOut(1:nVert))
  do iVert=1,nVert
    fxOut(iVert)=func(xOut(iVert))
  enddo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print '(/"Reconstruction de func avec la base des polynômes de Legendre")'
  
  allocate(mode (1:nMode,1:nVert))
  allocate(modeD(1:nMode,1:nVert))
  call     simplex1D(ord=order,r=xOut,  mode=mode ,transpose=.true.) !> mode (1:nMod,1:nVert)
  call gradSimplex1D(ord=order,r=xOut,drMode=modeD,transpose=.true.) !> modeD(1:nMod,1:nVert)
 !call display(title="mode",mat=mode)
  
  allocate(rFxOut (1:nVert))
  allocate(rFxOutD(1:nVert))
  do iVert=1,nVert
    rFxOut (iVert)=0d0
    rFxOutD(iVert)=0d0
    do iMode=1,nMode
      rFxOut (iVert)=rFxOut (iVert)+coef(iMode)*mode (iMode,iVert)
      rFxOutD(iVert)=rFxOutD(iVert)+coef(iMode)*modeD(iMode,iVert)
    enddo
  enddo
  
  
 !print '("size(rFxOut)=",i3)',size(rFxOut)
 !print '("rFxOut(",i2,")=",f12.5)',(iVert,rFxOut(iVert),iVert=1,nVert)
  
  print '("Sauvegarde des résultats dans le fichier rebuild1D.csv")'
  open(unit=10,file="rebuild1D.csv",action='write')
  write(10,'(" ""x"" , ""fonction discontinue"" , ""fonction reconstruite"" , ""dérivée fonction reconstruite"" ")')
  do iVert=1,nVert
    write(10,'(*(e22.15,","))')xOut(iVert),fxOut(iVert),rFxOut(iVert),rFxOutD(iVert)
  enddo
  close(10)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(coef)
  deallocate(xOut,fxOut,rFxOut,rFxOutD)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
  
contains
  
  function func(x) result(fx)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    real(8)            :: x,fx
    real(8), parameter :: pi=3.141592653589793
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    if( x<0d0 )then
      fx=sin(pi*x)-1d-1
    else
      fx=sin(pi*x)+1d-1
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end function func
  
end subroutine edge_02


program main
  call edge_02()
end program main