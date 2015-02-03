! Generated automatically by libmesh6

module M_libmesh6_api

public

interface
  integer function gmfopenmeshf77(FilNam, mod, ver, dim)
  character(len=*) :: FilNam
  integer          :: mod, ver, dim
  end function gmfopenmeshf77
end interface

interface
  integer function gmfclosemeshf77(MshIdx)
  integer   :: MshIdx
  end function gmfclosemeshf77
end interface

interface
  integer function gmfstatkwdf77(MshIdx, KwdIdx, NmbTyp, SolSiz, TypTab)
  integer   :: MshIdx, NmbTyp, KwdIdx, SolSiz, TypTab(*)
  end function gmfstatkwdf77
end interface

interface
  integer function gmfsetkwdf77(MshIdx, KwdIdx, NmbLin, NmbTyp, TypTab)
  integer   :: MshIdx, KwdIdx, NmbLin, NmbTyp, TypTab(*)
  end function gmfsetkwdf77
end interface

interface
  integer function gmfgotokwdf77(MshIdx, KwdIdx)
  integer   :: MshIdx, KwdIdx
  end function gmfgotokwdf77
end interface

interface
  subroutine gmfgetvertex2dr4(MshIdx, r0, r1, i2)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  integer :: i2
  end subroutine gmfgetvertex2dr4
end interface

interface
  subroutine gmfgetvertex2dr8(MshIdx, r0, r1, i2)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  integer :: i2
  end subroutine gmfgetvertex2dr8
end interface

interface
  subroutine gmfgetvertex3dr4(MshIdx, r0, r1, r2, i3)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  real(4) :: r2
  integer :: i3
  end subroutine gmfgetvertex3dr4
end interface

interface
  subroutine gmfgetvertex3dr8(MshIdx, r0, r1, r2, i3)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  real(8) :: r2
  integer :: i3
  end subroutine gmfgetvertex3dr8
end interface

interface
  subroutine gmfsetvertex2dr4(MshIdx, r0, r1, i2)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  integer :: i2
  end subroutine gmfsetvertex2dr4
end interface

interface
  subroutine gmfsetvertex2dr8(MshIdx, r0, r1, i2)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  integer :: i2
  end subroutine gmfsetvertex2dr8
end interface

interface
  subroutine gmfsetvertex3dr4(MshIdx, r0, r1, r2, i3)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  real(4) :: r2
  integer :: i3
  end subroutine gmfsetvertex3dr4
end interface

interface
  subroutine gmfsetvertex3dr8(MshIdx, r0, r1, r2, i3)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  real(8) :: r2
  integer :: i3
  end subroutine gmfsetvertex3dr8
end interface

interface
  subroutine gmfgetedge(MshIdx, i0, i1, i2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  end subroutine gmfgetedge
end interface

interface
  subroutine gmfsetedge(MshIdx, i0, i1, i2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  end subroutine gmfsetedge
end interface

interface
  subroutine gmfgettriangle(MshIdx, i0, i1, i2, i3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  end subroutine gmfgettriangle
end interface

interface
  subroutine gmfsettriangle(MshIdx, i0, i1, i2, i3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  end subroutine gmfsettriangle
end interface

interface
  subroutine gmfgetquadrilateral(MshIdx, i0, i1, i2, i3, i4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  end subroutine gmfgetquadrilateral
end interface

interface
  subroutine gmfsetquadrilateral(MshIdx, i0, i1, i2, i3, i4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  end subroutine gmfsetquadrilateral
end interface

interface
  subroutine gmfgettetrahedron(MshIdx, i0, i1, i2, i3, i4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  end subroutine gmfgettetrahedron
end interface

interface
  subroutine gmfsettetrahedron(MshIdx, i0, i1, i2, i3, i4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  end subroutine gmfsettetrahedron
end interface

interface
  subroutine gmfgetprism(MshIdx, i0, i1, i2, i3, i4, i5, i6)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  end subroutine gmfgetprism
end interface

interface
  subroutine gmfsetprism(MshIdx, i0, i1, i2, i3, i4, i5, i6)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  end subroutine gmfsetprism
end interface

interface
  subroutine gmfgethexahedron(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7, i8)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  integer :: i8
  end subroutine gmfgethexahedron
end interface

interface
  subroutine gmfsethexahedron(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7, i8)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  integer :: i8
  end subroutine gmfsethexahedron
end interface

interface
  subroutine gmfgetcorner(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetcorner
end interface

interface
  subroutine gmfsetcorner(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetcorner
end interface

interface
  subroutine gmfgetridge(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetridge
end interface

interface
  subroutine gmfsetridge(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetridge
end interface

interface
  subroutine gmfgetrequiredvertex(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetrequiredvertex
end interface

interface
  subroutine gmfsetrequiredvertex(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetrequiredvertex
end interface

interface
  subroutine gmfgetrequirededge(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetrequirededge
end interface

interface
  subroutine gmfsetrequirededge(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetrequirededge
end interface

interface
  subroutine gmfgetrequiredtriangle(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetrequiredtriangle
end interface

interface
  subroutine gmfsetrequiredtriangle(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetrequiredtriangle
end interface

interface
  subroutine gmfgetrequiredquadrilateral(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetrequiredquadrilateral
end interface

interface
  subroutine gmfsetrequiredquadrilateral(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetrequiredquadrilateral
end interface

interface
  subroutine gmfgettangentatedgevertex(MshIdx, i0, i1, i2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  end subroutine gmfgettangentatedgevertex
end interface

interface
  subroutine gmfsettangentatedgevertex(MshIdx, i0, i1, i2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  end subroutine gmfsettangentatedgevertex
end interface

interface
  subroutine gmfgetnormalatvertex(MshIdx, i0, i1)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  end subroutine gmfgetnormalatvertex
end interface

interface
  subroutine gmfsetnormalatvertex(MshIdx, i0, i1)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  end subroutine gmfsetnormalatvertex
end interface

interface
  subroutine gmfgetnormalattrianglevertex(MshIdx, i0, i1, i2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  end subroutine gmfgetnormalattrianglevertex
end interface

interface
  subroutine gmfsetnormalattrianglevertex(MshIdx, i0, i1, i2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  end subroutine gmfsetnormalattrianglevertex
end interface

interface
  subroutine gmfgetnormalatquadrilateralvertex(MshIdx, i0, i1, i2, i3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  end subroutine gmfgetnormalatquadrilateralvertex
end interface

interface
  subroutine gmfsetnormalatquadrilateralvertex(MshIdx, i0, i1, i2, i3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  end subroutine gmfsetnormalatquadrilateralvertex
end interface

interface
  subroutine gmfgettrianglep2(MshIdx, i0, i1, i2, i3, i4, i5, i6)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  end subroutine gmfgettrianglep2
end interface

interface
  subroutine gmfsettrianglep2(MshIdx, i0, i1, i2, i3, i4, i5, i6)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  end subroutine gmfsettrianglep2
end interface

interface
  subroutine gmfgetedgep2(MshIdx, i0, i1, i2, i3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  end subroutine gmfgetedgep2
end interface

interface
  subroutine gmfsetedgep2(MshIdx, i0, i1, i2, i3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  end subroutine gmfsetedgep2
end interface

interface
  subroutine gmfgetsolatpyramidr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfgetsolatpyramidr4
end interface

interface
  subroutine gmfgetsolatpyramidr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfgetsolatpyramidr8
end interface

interface
  subroutine gmfsetsolatpyramidr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfsetsolatpyramidr4
end interface

interface
  subroutine gmfsetsolatpyramidr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfsetsolatpyramidr8
end interface

interface
  subroutine gmfgetquadrilateralq2(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  integer :: i8
  integer :: i9
  end subroutine gmfgetquadrilateralq2
end interface

interface
  subroutine gmfsetquadrilateralq2(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  integer :: i8
  integer :: i9
  end subroutine gmfsetquadrilateralq2
end interface

interface
  subroutine gmfgetisolatpyramid(MshIdx, i0, i1, i2, i3, i4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  end subroutine gmfgetisolatpyramid
end interface

interface
  subroutine gmfsetisolatpyramid(MshIdx, i0, i1, i2, i3, i4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  end subroutine gmfsetisolatpyramid
end interface

interface
  subroutine gmfgetsubdomainfromgeom(MshIdx, i0, i1, i2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  end subroutine gmfgetsubdomainfromgeom
end interface

interface
  subroutine gmfsetsubdomainfromgeom(MshIdx, i0, i1, i2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  end subroutine gmfsetsubdomainfromgeom
end interface

interface
  subroutine gmfgettetrahedronp2(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  integer :: i8
  integer :: i9
  integer :: i10
  end subroutine gmfgettetrahedronp2
end interface

interface
  subroutine gmfsettetrahedronp2(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  integer :: i8
  integer :: i9
  integer :: i10
  end subroutine gmfsettetrahedronp2
end interface

interface
  subroutine gmfgetfault_neartri(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetfault_neartri
end interface

interface
  subroutine gmfsetfault_neartri(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetfault_neartri
end interface

interface
  subroutine gmfgetfault_inter(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetfault_inter
end interface

interface
  subroutine gmfsetfault_inter(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetfault_inter
end interface

interface
  subroutine gmfgethexahedronq2(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23, i24, i25, i26, i27)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  integer :: i8
  integer :: i9
  integer :: i10
  integer :: i11
  integer :: i12
  integer :: i13
  integer :: i14
  integer :: i15
  integer :: i16
  integer :: i17
  integer :: i18
  integer :: i19
  integer :: i20
  integer :: i21
  integer :: i22
  integer :: i23
  integer :: i24
  integer :: i25
  integer :: i26
  integer :: i27
  end subroutine gmfgethexahedronq2
end interface

interface
  subroutine gmfsethexahedronq2(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23, i24, i25, i26, i27)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  integer :: i8
  integer :: i9
  integer :: i10
  integer :: i11
  integer :: i12
  integer :: i13
  integer :: i14
  integer :: i15
  integer :: i16
  integer :: i17
  integer :: i18
  integer :: i19
  integer :: i20
  integer :: i21
  integer :: i22
  integer :: i23
  integer :: i24
  integer :: i25
  integer :: i26
  integer :: i27
  end subroutine gmfsethexahedronq2
end interface

interface
  subroutine gmfgetextraverticesatedge(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetextraverticesatedge
end interface

interface
  subroutine gmfsetextraverticesatedge(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetextraverticesatedge
end interface

interface
  subroutine gmfgetextraverticesattriangle(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetextraverticesattriangle
end interface

interface
  subroutine gmfsetextraverticesattriangle(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetextraverticesattriangle
end interface

interface
  subroutine gmfgetextraverticesatquadrilateral(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetextraverticesatquadrilateral
end interface

interface
  subroutine gmfsetextraverticesatquadrilateral(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetextraverticesatquadrilateral
end interface

interface
  subroutine gmfgetextraverticesattetrahedron(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetextraverticesattetrahedron
end interface

interface
  subroutine gmfsetextraverticesattetrahedron(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetextraverticesattetrahedron
end interface

interface
  subroutine gmfgetextraverticesatprism(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetextraverticesatprism
end interface

interface
  subroutine gmfsetextraverticesatprism(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetextraverticesatprism
end interface

interface
  subroutine gmfgetextraverticesathexahedron(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetextraverticesathexahedron
end interface

interface
  subroutine gmfsetextraverticesathexahedron(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetextraverticesathexahedron
end interface

interface
  subroutine gmfgetvertexongeometricvertexr4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfgetvertexongeometricvertexr4
end interface

interface
  subroutine gmfgetvertexongeometricvertexr8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfgetvertexongeometricvertexr8
end interface

interface
  subroutine gmfsetvertexongeometricvertexr4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfsetvertexongeometricvertexr4
end interface

interface
  subroutine gmfsetvertexongeometricvertexr8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfsetvertexongeometricvertexr8
end interface

interface
  subroutine gmfgetvertexongeometricedger4(MshIdx, i0, i1, r2, r3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  real(4) :: r3
  end subroutine gmfgetvertexongeometricedger4
end interface

interface
  subroutine gmfgetvertexongeometricedger8(MshIdx, i0, i1, r2, r3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  real(8) :: r3
  end subroutine gmfgetvertexongeometricedger8
end interface

interface
  subroutine gmfsetvertexongeometricedger4(MshIdx, i0, i1, r2, r3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  real(4) :: r3
  end subroutine gmfsetvertexongeometricedger4
end interface

interface
  subroutine gmfsetvertexongeometricedger8(MshIdx, i0, i1, r2, r3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  real(8) :: r3
  end subroutine gmfsetvertexongeometricedger8
end interface

interface
  subroutine gmfgetvertexongeometrictriangler4(MshIdx, i0, i1, r2, r3, r4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  real(4) :: r3
  real(4) :: r4
  end subroutine gmfgetvertexongeometrictriangler4
end interface

interface
  subroutine gmfgetvertexongeometrictriangler8(MshIdx, i0, i1, r2, r3, r4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  real(8) :: r3
  real(8) :: r4
  end subroutine gmfgetvertexongeometrictriangler8
end interface

interface
  subroutine gmfsetvertexongeometrictriangler4(MshIdx, i0, i1, r2, r3, r4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  real(4) :: r3
  real(4) :: r4
  end subroutine gmfsetvertexongeometrictriangler4
end interface

interface
  subroutine gmfsetvertexongeometrictriangler8(MshIdx, i0, i1, r2, r3, r4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  real(8) :: r3
  real(8) :: r4
  end subroutine gmfsetvertexongeometrictriangler8
end interface

interface
  subroutine gmfgetvertexongeometricquadrilateralr4(MshIdx, i0, i1, r2, r3, r4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  real(4) :: r3
  real(4) :: r4
  end subroutine gmfgetvertexongeometricquadrilateralr4
end interface

interface
  subroutine gmfgetvertexongeometricquadrilateralr8(MshIdx, i0, i1, r2, r3, r4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  real(8) :: r3
  real(8) :: r4
  end subroutine gmfgetvertexongeometricquadrilateralr8
end interface

interface
  subroutine gmfsetvertexongeometricquadrilateralr4(MshIdx, i0, i1, r2, r3, r4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  real(4) :: r3
  real(4) :: r4
  end subroutine gmfsetvertexongeometricquadrilateralr4
end interface

interface
  subroutine gmfsetvertexongeometricquadrilateralr8(MshIdx, i0, i1, r2, r3, r4)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  real(8) :: r3
  real(8) :: r4
  end subroutine gmfsetvertexongeometricquadrilateralr8
end interface

interface
  subroutine gmfgetedgeongeometricedger4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfgetedgeongeometricedger4
end interface

interface
  subroutine gmfgetedgeongeometricedger8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfgetedgeongeometricedger8
end interface

interface
  subroutine gmfsetedgeongeometricedger4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfsetedgeongeometricedger4
end interface

interface
  subroutine gmfsetedgeongeometricedger8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfsetedgeongeometricedger8
end interface

interface
  subroutine gmfgetfault_freeedge(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetfault_freeedge
end interface

interface
  subroutine gmfsetfault_freeedge(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetfault_freeedge
end interface

interface
  subroutine gmfgetpolyhedron(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23, i24, i25, i26, i27, i28, i29, i30, i31, i32)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  integer :: i8
  integer :: i9
  integer :: i10
  integer :: i11
  integer :: i12
  integer :: i13
  integer :: i14
  integer :: i15
  integer :: i16
  integer :: i17
  integer :: i18
  integer :: i19
  integer :: i20
  integer :: i21
  integer :: i22
  integer :: i23
  integer :: i24
  integer :: i25
  integer :: i26
  integer :: i27
  integer :: i28
  integer :: i29
  integer :: i30
  integer :: i31
  integer :: i32
  end subroutine gmfgetpolyhedron
end interface

interface
  subroutine gmfsetpolyhedron(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23, i24, i25, i26, i27, i28, i29, i30, i31, i32)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  integer :: i8
  integer :: i9
  integer :: i10
  integer :: i11
  integer :: i12
  integer :: i13
  integer :: i14
  integer :: i15
  integer :: i16
  integer :: i17
  integer :: i18
  integer :: i19
  integer :: i20
  integer :: i21
  integer :: i22
  integer :: i23
  integer :: i24
  integer :: i25
  integer :: i26
  integer :: i27
  integer :: i28
  integer :: i29
  integer :: i30
  integer :: i31
  integer :: i32
  end subroutine gmfsetpolyhedron
end interface

interface
  subroutine gmfgetfault_overlap(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetfault_overlap
end interface

interface
  subroutine gmfsetfault_overlap(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetfault_overlap
end interface

interface
  subroutine gmfgetpyramid(MshIdx, i0, i1, i2, i3, i4, i5)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  end subroutine gmfgetpyramid
end interface

interface
  subroutine gmfsetpyramid(MshIdx, i0, i1, i2, i3, i4, i5)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  end subroutine gmfsetpyramid
end interface

interface
  subroutine gmfgetprivatetable(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetprivatetable
end interface

interface
  subroutine gmfsetprivatetable(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetprivatetable
end interface

interface
  subroutine gmfgetfault_badshape(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetfault_badshape
end interface

interface
  subroutine gmfsetfault_badshape(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetfault_badshape
end interface

interface
  subroutine gmfgettriangleongeometrictriangler4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfgettriangleongeometrictriangler4
end interface

interface
  subroutine gmfgettriangleongeometrictriangler8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfgettriangleongeometrictriangler8
end interface

interface
  subroutine gmfsettriangleongeometrictriangler4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfsettriangleongeometrictriangler4
end interface

interface
  subroutine gmfsettriangleongeometrictriangler8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfsettriangleongeometrictriangler8
end interface

interface
  subroutine gmfgettriangleongeometricquadrilateralr4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfgettriangleongeometricquadrilateralr4
end interface

interface
  subroutine gmfgettriangleongeometricquadrilateralr8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfgettriangleongeometricquadrilateralr8
end interface

interface
  subroutine gmfsettriangleongeometricquadrilateralr4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfsettriangleongeometricquadrilateralr4
end interface

interface
  subroutine gmfsettriangleongeometricquadrilateralr8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfsettriangleongeometricquadrilateralr8
end interface

interface
  subroutine gmfgetquadrilateralongeometrictriangler4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfgetquadrilateralongeometrictriangler4
end interface

interface
  subroutine gmfgetquadrilateralongeometrictriangler8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfgetquadrilateralongeometrictriangler8
end interface

interface
  subroutine gmfsetquadrilateralongeometrictriangler4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfsetquadrilateralongeometrictriangler4
end interface

interface
  subroutine gmfsetquadrilateralongeometrictriangler8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfsetquadrilateralongeometrictriangler8
end interface

interface
  subroutine gmfgetquadrilateralongeometricquadrilateralr4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfgetquadrilateralongeometricquadrilateralr4
end interface

interface
  subroutine gmfgetquadrilateralongeometricquadrilateralr8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfgetquadrilateralongeometricquadrilateralr8
end interface

interface
  subroutine gmfsetquadrilateralongeometricquadrilateralr4(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(4) :: r2
  end subroutine gmfsetquadrilateralongeometricquadrilateralr4
end interface

interface
  subroutine gmfsetquadrilateralongeometricquadrilateralr8(MshIdx, i0, i1, r2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  real(8) :: r2
  end subroutine gmfsetquadrilateralongeometricquadrilateralr8
end interface

interface
  subroutine gmfgettangent2dr4(MshIdx, r0, r1)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  end subroutine gmfgettangent2dr4
end interface

interface
  subroutine gmfgettangent2dr8(MshIdx, r0, r1)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  end subroutine gmfgettangent2dr8
end interface

interface
  subroutine gmfgettangent3dr4(MshIdx, r0, r1, r2)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  real(4) :: r2
  end subroutine gmfgettangent3dr4
end interface

interface
  subroutine gmfgettangent3dr8(MshIdx, r0, r1, r2)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  real(8) :: r2
  end subroutine gmfgettangent3dr8
end interface

interface
  subroutine gmfsettangent2dr4(MshIdx, r0, r1)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  end subroutine gmfsettangent2dr4
end interface

interface
  subroutine gmfsettangent2dr8(MshIdx, r0, r1)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  end subroutine gmfsettangent2dr8
end interface

interface
  subroutine gmfsettangent3dr4(MshIdx, r0, r1, r2)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  real(4) :: r2
  end subroutine gmfsettangent3dr4
end interface

interface
  subroutine gmfsettangent3dr8(MshIdx, r0, r1, r2)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  real(8) :: r2
  end subroutine gmfsettangent3dr8
end interface

interface
  subroutine gmfgetnormal2dr4(MshIdx, r0, r1)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  end subroutine gmfgetnormal2dr4
end interface

interface
  subroutine gmfgetnormal2dr8(MshIdx, r0, r1)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  end subroutine gmfgetnormal2dr8
end interface

interface
  subroutine gmfgetnormal3dr4(MshIdx, r0, r1, r2)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  real(4) :: r2
  end subroutine gmfgetnormal3dr4
end interface

interface
  subroutine gmfgetnormal3dr8(MshIdx, r0, r1, r2)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  real(8) :: r2
  end subroutine gmfgetnormal3dr8
end interface

interface
  subroutine gmfsetnormal2dr4(MshIdx, r0, r1)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  end subroutine gmfsetnormal2dr4
end interface

interface
  subroutine gmfsetnormal2dr8(MshIdx, r0, r1)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  end subroutine gmfsetnormal2dr8
end interface

interface
  subroutine gmfsetnormal3dr4(MshIdx, r0, r1, r2)
  integer :: MshIdx
  real(4) :: r0
  real(4) :: r1
  real(4) :: r2
  end subroutine gmfsetnormal3dr4
end interface

interface
  subroutine gmfsetnormal3dr8(MshIdx, r0, r1, r2)
  integer :: MshIdx
  real(8) :: r0
  real(8) :: r1
  real(8) :: r2
  end subroutine gmfsetnormal3dr8
end interface

interface
  subroutine gmfgettangentatvertex(MshIdx, i0, i1)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  end subroutine gmfgettangentatvertex
end interface

interface
  subroutine gmfsettangentatvertex(MshIdx, i0, i1)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  end subroutine gmfsettangentatvertex
end interface

interface
  subroutine gmfgetsolatvertexr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfgetsolatvertexr4
end interface

interface
  subroutine gmfgetsolatvertexr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfgetsolatvertexr8
end interface

interface
  subroutine gmfsetsolatvertexr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfsetsolatvertexr4
end interface

interface
  subroutine gmfsetsolatvertexr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfsetsolatvertexr8
end interface

interface
  subroutine gmfgetsolatedger4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfgetsolatedger4
end interface

interface
  subroutine gmfgetsolatedger8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfgetsolatedger8
end interface

interface
  subroutine gmfsetsolatedger4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfsetsolatedger4
end interface

interface
  subroutine gmfsetsolatedger8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfsetsolatedger8
end interface

interface
  subroutine gmfgetsolattriangler4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfgetsolattriangler4
end interface

interface
  subroutine gmfgetsolattriangler8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfgetsolattriangler8
end interface

interface
  subroutine gmfsetsolattriangler4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfsetsolattriangler4
end interface

interface
  subroutine gmfsetsolattriangler8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfsetsolattriangler8
end interface

interface
  subroutine gmfgetsolatquadrilateralr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfgetsolatquadrilateralr4
end interface

interface
  subroutine gmfgetsolatquadrilateralr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfgetsolatquadrilateralr8
end interface

interface
  subroutine gmfsetsolatquadrilateralr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfsetsolatquadrilateralr4
end interface

interface
  subroutine gmfsetsolatquadrilateralr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfsetsolatquadrilateralr8
end interface

interface
  subroutine gmfgetsolattetrahedronr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfgetsolattetrahedronr4
end interface

interface
  subroutine gmfgetsolattetrahedronr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfgetsolattetrahedronr8
end interface

interface
  subroutine gmfsetsolattetrahedronr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfsetsolattetrahedronr4
end interface

interface
  subroutine gmfsetsolattetrahedronr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfsetsolattetrahedronr8
end interface

interface
  subroutine gmfgetsolatprismr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfgetsolatprismr4
end interface

interface
  subroutine gmfgetsolatprismr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfgetsolatprismr8
end interface

interface
  subroutine gmfsetsolatprismr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfsetsolatprismr4
end interface

interface
  subroutine gmfsetsolatprismr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfsetsolatprismr8
end interface

interface
  subroutine gmfgetsolathexahedronr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfgetsolathexahedronr4
end interface

interface
  subroutine gmfgetsolathexahedronr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfgetsolathexahedronr8
end interface

interface
  subroutine gmfsetsolathexahedronr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfsetsolathexahedronr4
end interface

interface
  subroutine gmfsetsolathexahedronr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfsetsolathexahedronr8
end interface

interface
  subroutine gmfgetdsolatvertexr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfgetdsolatvertexr4
end interface

interface
  subroutine gmfgetdsolatvertexr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfgetdsolatvertexr8
end interface

interface
  subroutine gmfsetdsolatvertexr4(MshIdx, SolTab)
  integer :: MshIdx
  real(4) :: SolTab(*)
  end subroutine gmfsetdsolatvertexr4
end interface

interface
  subroutine gmfsetdsolatvertexr8(MshIdx, SolTab)
  integer :: MshIdx
  real(8) :: SolTab(*)
  end subroutine gmfsetdsolatvertexr8
end interface

interface
  subroutine gmfgetisolatvertex(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetisolatvertex
end interface

interface
  subroutine gmfsetisolatvertex(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetisolatvertex
end interface

interface
  subroutine gmfgetisolatedge(MshIdx, i0, i1)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  end subroutine gmfgetisolatedge
end interface

interface
  subroutine gmfsetisolatedge(MshIdx, i0, i1)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  end subroutine gmfsetisolatedge
end interface

interface
  subroutine gmfgetisolattriangle(MshIdx, i0, i1, i2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  end subroutine gmfgetisolattriangle
end interface

interface
  subroutine gmfsetisolattriangle(MshIdx, i0, i1, i2)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  end subroutine gmfsetisolattriangle
end interface

interface
  subroutine gmfgetisolatquadrilateral(MshIdx, i0, i1, i2, i3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  end subroutine gmfgetisolatquadrilateral
end interface

interface
  subroutine gmfsetisolatquadrilateral(MshIdx, i0, i1, i2, i3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  end subroutine gmfsetisolatquadrilateral
end interface

interface
  subroutine gmfgetisolattetrahedron(MshIdx, i0, i1, i2, i3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  end subroutine gmfgetisolattetrahedron
end interface

interface
  subroutine gmfsetisolattetrahedron(MshIdx, i0, i1, i2, i3)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  end subroutine gmfsetisolattetrahedron
end interface

interface
  subroutine gmfgetisolatprism(MshIdx, i0, i1, i2, i3, i4, i5)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  end subroutine gmfgetisolatprism
end interface

interface
  subroutine gmfsetisolatprism(MshIdx, i0, i1, i2, i3, i4, i5)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  end subroutine gmfsetisolatprism
end interface

interface
  subroutine gmfgetisolathexahedron(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  end subroutine gmfgetisolathexahedron
end interface

interface
  subroutine gmfsetisolathexahedron(MshIdx, i0, i1, i2, i3, i4, i5, i6, i7)
  integer :: MshIdx
  integer :: i0
  integer :: i1
  integer :: i2
  integer :: i3
  integer :: i4
  integer :: i5
  integer :: i6
  integer :: i7
  end subroutine gmfsetisolathexahedron
end interface

interface
  subroutine gmfgetfault_smalltri(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetfault_smalltri
end interface

interface
  subroutine gmfsetfault_smalltri(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetfault_smalltri
end interface

interface
  subroutine gmfgetcoarsehexahedron(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfgetcoarsehexahedron
end interface

interface
  subroutine gmfsetcoarsehexahedron(MshIdx, i0)
  integer :: MshIdx
  integer :: i0
  end subroutine gmfsetcoarsehexahedron
end interface

interface
  subroutine gmfgetcomment(MshIdx)
  integer :: MshIdx
  end subroutine gmfgetcomment
end interface

interface
  subroutine gmfsetcomment(MshIdx)
  integer :: MshIdx
  end subroutine gmfsetcomment
end interface

integer,parameter :: GmfMaxTyp=1000
integer,parameter :: GmfMaxKwd=81
integer,parameter :: GmfRead=1
integer,parameter :: GmfWrite=2
integer,parameter :: GmfSca=1
integer,parameter :: GmfVec=2
integer,parameter :: GmfSymMat=3
integer,parameter :: GmfMat=4

integer,parameter :: GmfMeshVersionFormatted=1
integer,parameter :: GmfDimension=3
integer,parameter :: GmfVertices=4
integer,parameter :: GmfEdges=5
integer,parameter :: GmfTriangles=6
integer,parameter :: GmfQuadrilaterals=7
integer,parameter :: GmfTetrahedra=8
integer,parameter :: GmfPrisms=9
integer,parameter :: GmfHexahedra=10
integer,parameter :: GmfIterationsAll=11
integer,parameter :: GmfTimesAll=12
integer,parameter :: GmfCorners=13
integer,parameter :: GmfRidges=14
integer,parameter :: GmfRequiredVertices=15
integer,parameter :: GmfRequiredEdges=16
integer,parameter :: GmfRequiredTriangles=17
integer,parameter :: GmfRequiredQuadrilaterals=18
integer,parameter :: GmfTangentAtEdgeVertices=19
integer,parameter :: GmfNormalAtVertices=20
integer,parameter :: GmfNormalAtTriangleVertices=21
integer,parameter :: GmfNormalAtQuadrilateralVertices=22
integer,parameter :: GmfAngleOfCornerBound=23
integer,parameter :: GmfTrianglesP2=24
integer,parameter :: GmfEdgesP2=25
integer,parameter :: GmfSolAtPyramids=26
integer,parameter :: GmfQuadrilateralsQ2=27
integer,parameter :: GmfISolAtPyramids=28
integer,parameter :: GmfSubDomainFromGeom=29
integer,parameter :: GmfTetrahedraP2=30
integer,parameter :: GmfFault_NearTri=31
integer,parameter :: GmfFault_Inter=32
integer,parameter :: GmfHexahedraQ2=33
integer,parameter :: GmfExtraVerticesAtEdges=34
integer,parameter :: GmfExtraVerticesAtTriangles=35
integer,parameter :: GmfExtraVerticesAtQuadrilaterals=36
integer,parameter :: GmfExtraVerticesAtTetrahedra=37
integer,parameter :: GmfExtraVerticesAtPrisms=38
integer,parameter :: GmfExtraVerticesAtHexahedra=39
integer,parameter :: GmfVerticesOnGeometricVertices=40
integer,parameter :: GmfVerticesOnGeometricEdges=41
integer,parameter :: GmfVerticesOnGeometricTriangles=42
integer,parameter :: GmfVerticesOnGeometricQuadrilaterals=43
integer,parameter :: GmfEdgesOnGeometricEdges=44
integer,parameter :: GmfFault_FreeEdge=45
integer,parameter :: GmfPolyhedra=46
integer,parameter :: GmfPolygons=47
integer,parameter :: GmfFault_Overlap=48
integer,parameter :: GmfPyramids=49
integer,parameter :: GmfBoundingBox=50
integer,parameter :: GmfBody=51
integer,parameter :: GmfPrivateTable=52
integer,parameter :: GmfFault_BadShape=53
integer,parameter :: GmfEnd=54
integer,parameter :: GmfTrianglesOnGeometricTriangles=55
integer,parameter :: GmfTrianglesOnGeometricQuadrilaterals=56
integer,parameter :: GmfQuadrilateralsOnGeometricTriangles=57
integer,parameter :: GmfQuadrilateralsOnGeometricQuadrilaterals=58
integer,parameter :: GmfTangents=59
integer,parameter :: GmfNormals=60
integer,parameter :: GmfTangentAtVertices=61
integer,parameter :: GmfSolAtVertices=62
integer,parameter :: GmfSolAtEdges=63
integer,parameter :: GmfSolAtTriangles=64
integer,parameter :: GmfSolAtQuadrilaterals=65
integer,parameter :: GmfSolAtTetrahedra=66
integer,parameter :: GmfSolAtPrisms=67
integer,parameter :: GmfSolAtHexahedra=68
integer,parameter :: GmfDSolAtVertices=69
integer,parameter :: GmfISolAtVertices=70
integer,parameter :: GmfISolAtEdges=71
integer,parameter :: GmfISolAtTriangles=72
integer,parameter :: GmfISolAtQuadrilaterals=73
integer,parameter :: GmfISolAtTetrahedra=74
integer,parameter :: GmfISolAtPrisms=75
integer,parameter :: GmfISolAtHexahedra=76
integer,parameter :: GmfIterations=77
integer,parameter :: GmfTime=78
integer,parameter :: GmfFault_SmallTri=79
integer,parameter :: GmfCoarseHexahedra=80
integer,parameter :: GmfComments=81

end module M_libmesh6_api

