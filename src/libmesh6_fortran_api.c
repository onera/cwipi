#include "libmesh6.h"
#include <string.h>
char *GmfKwdFmt[ GmfMaxKwd + 1 ][4];


/*----------------------------------------------------------*/
/* Fortran 77 API											*/
/*----------------------------------------------------------*/

#if defined(F77_NO_UNDER_SCORE)
#define call(x) x
#else
#define call(x) x ## _
#endif


/*----------------------------------------------------------*/
/* Fortran 77 API											*/
/*----------------------------------------------------------*/

int call(gmfopenmeshf77)(char *FilNam, int *mod, int *ver, int *dim, int StrSiz)
{
	int i;
	char TmpNam[ GmfStrSiz ];

	for(i=0;i<StrSiz;i++)
		TmpNam[i] = FilNam[i];

	TmpNam[ StrSiz ] = 0;

	if(*mod == GmfRead)
		return(GmfOpenMesh(TmpNam, *mod, ver, dim));
	else
		return(GmfOpenMesh(TmpNam, *mod, *ver, *dim));
}

int call(gmfclosemeshf77)(int *idx)
{
	return(GmfCloseMesh(*idx));
}

int call(gmfstatkwdf77)(int *MshIdx, int *KwdIdx, int *NmbTyp, int *SolSiz, int *TypTab)
{
	if(!strcmp(GmfKwdFmt[ *KwdIdx ][3], "sr"))
		return(GmfStatKwd(*MshIdx, *KwdIdx, NmbTyp, SolSiz, TypTab));
	else
		return(GmfStatKwd(*MshIdx, *KwdIdx));
}

int call(gmfgotokwdf77)(int *MshIdx, int *KwdIdx)
{
	return(GmfGotoKwd(*MshIdx, *KwdIdx));
}

int call(gmfsetkwdf77)(int *MshIdx, int *KwdIdx, int *NmbLin, int *NmbTyp, int *TypTab)
{
	if(!strcmp(GmfKwdFmt[ *KwdIdx ][3], "sr"))
		return(GmfSetKwd(*MshIdx, *KwdIdx, *NmbLin, *NmbTyp, TypTab));
	else if(strlen(GmfKwdFmt[ *KwdIdx ][2]))
		return(GmfSetKwd(*MshIdx, *KwdIdx, *NmbLin));
	else
		return(GmfSetKwd(*MshIdx, *KwdIdx));
}
/* Generated automatically by libmesh6 */

void call(gmfgetvertex2dr4)(int *MshIdx, float *r0, float *r1, int *i2)
{
	GmfGetLin(*MshIdx, GmfVertices, r0, r1, i2);
}

void call(gmfgetvertex2dr8)(int *MshIdx, double *r0, double *r1, int *i2)
{
	GmfGetLin(*MshIdx, GmfVertices, r0, r1, i2);
}

void call(gmfgetvertex3dr4)(int *MshIdx, float *r0, float *r1, float *r2, int *i3)
{
	GmfGetLin(*MshIdx, GmfVertices, r0, r1, r2, i3);
}

void call(gmfgetvertex3dr8)(int *MshIdx, double *r0, double *r1, double *r2, int *i3)
{
	GmfGetLin(*MshIdx, GmfVertices, r0, r1, r2, i3);
}

void call(gmfsetvertex2dr4)(int *MshIdx, float *r0, float *r1, int *i2)
{
	GmfSetLin(*MshIdx, GmfVertices, *r0, *r1, *i2);
}

void call(gmfsetvertex2dr8)(int *MshIdx, double *r0, double *r1, int *i2)
{
	GmfSetLin(*MshIdx, GmfVertices, *r0, *r1, *i2);
}

void call(gmfsetvertex3dr4)(int *MshIdx, float *r0, float *r1, float *r2, int *i3)
{
	GmfSetLin(*MshIdx, GmfVertices, *r0, *r1, *r2, *i3);
}

void call(gmfsetvertex3dr8)(int *MshIdx, double *r0, double *r1, double *r2, int *i3)
{
	GmfSetLin(*MshIdx, GmfVertices, *r0, *r1, *r2, *i3);
}

void call(gmfgetedge)(int *MshIdx, int *i0, int *i1, int *i2)
{
	GmfGetLin(*MshIdx, GmfEdges, i0, i1, i2);
}

void call(gmfsetedge)(int *MshIdx, int *i0, int *i1, int *i2)
{
	GmfSetLin(*MshIdx, GmfEdges, *i0, *i1, *i2);
}

void call(gmfgettriangle)(int *MshIdx, int *i0, int *i1, int *i2, int *i3)
{
	GmfGetLin(*MshIdx, GmfTriangles, i0, i1, i2, i3);
}

void call(gmfsettriangle)(int *MshIdx, int *i0, int *i1, int *i2, int *i3)
{
	GmfSetLin(*MshIdx, GmfTriangles, *i0, *i1, *i2, *i3);
}

void call(gmfgetquadrilateral)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4)
{
	GmfGetLin(*MshIdx, GmfQuadrilaterals, i0, i1, i2, i3, i4);
}

void call(gmfsetquadrilateral)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4)
{
	GmfSetLin(*MshIdx, GmfQuadrilaterals, *i0, *i1, *i2, *i3, *i4);
}

void call(gmfgettetrahedron)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4)
{
	GmfGetLin(*MshIdx, GmfTetrahedra, i0, i1, i2, i3, i4);
}

void call(gmfsettetrahedron)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4)
{
	GmfSetLin(*MshIdx, GmfTetrahedra, *i0, *i1, *i2, *i3, *i4);
}

void call(gmfgetprism)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6)
{
	GmfGetLin(*MshIdx, GmfPrisms, i0, i1, i2, i3, i4, i5, i6);
}

void call(gmfsetprism)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6)
{
	GmfSetLin(*MshIdx, GmfPrisms, *i0, *i1, *i2, *i3, *i4, *i5, *i6);
}

void call(gmfgethexahedron)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8)
{
	GmfGetLin(*MshIdx, GmfHexahedra, i0, i1, i2, i3, i4, i5, i6, i7, i8);
}

void call(gmfsethexahedron)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8)
{
	GmfSetLin(*MshIdx, GmfHexahedra, *i0, *i1, *i2, *i3, *i4, *i5, *i6, *i7, *i8);
}

void call(gmfgetcorner)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfCorners, i0);
}

void call(gmfsetcorner)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfCorners, *i0);
}

void call(gmfgetridge)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfRidges, i0);
}

void call(gmfsetridge)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfRidges, *i0);
}

void call(gmfgetrequiredvertex)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfRequiredVertices, i0);
}

void call(gmfsetrequiredvertex)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfRequiredVertices, *i0);
}

void call(gmfgetrequirededge)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfRequiredEdges, i0);
}

void call(gmfsetrequirededge)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfRequiredEdges, *i0);
}

void call(gmfgetrequiredtriangle)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfRequiredTriangles, i0);
}

void call(gmfsetrequiredtriangle)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfRequiredTriangles, *i0);
}

void call(gmfgetrequiredquadrilateral)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfRequiredQuadrilaterals, i0);
}

void call(gmfsetrequiredquadrilateral)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfRequiredQuadrilaterals, *i0);
}

void call(gmfgettangentatedgevertex)(int *MshIdx, int *i0, int *i1, int *i2)
{
	GmfGetLin(*MshIdx, GmfTangentAtEdgeVertices, i0, i1, i2);
}

void call(gmfsettangentatedgevertex)(int *MshIdx, int *i0, int *i1, int *i2)
{
	GmfSetLin(*MshIdx, GmfTangentAtEdgeVertices, *i0, *i1, *i2);
}

void call(gmfgetnormalatvertex)(int *MshIdx, int *i0, int *i1)
{
	GmfGetLin(*MshIdx, GmfNormalAtVertices, i0, i1);
}

void call(gmfsetnormalatvertex)(int *MshIdx, int *i0, int *i1)
{
	GmfSetLin(*MshIdx, GmfNormalAtVertices, *i0, *i1);
}

void call(gmfgetnormalattrianglevertex)(int *MshIdx, int *i0, int *i1, int *i2)
{
	GmfGetLin(*MshIdx, GmfNormalAtTriangleVertices, i0, i1, i2);
}

void call(gmfsetnormalattrianglevertex)(int *MshIdx, int *i0, int *i1, int *i2)
{
	GmfSetLin(*MshIdx, GmfNormalAtTriangleVertices, *i0, *i1, *i2);
}

void call(gmfgetnormalatquadrilateralvertex)(int *MshIdx, int *i0, int *i1, int *i2, int *i3)
{
	GmfGetLin(*MshIdx, GmfNormalAtQuadrilateralVertices, i0, i1, i2, i3);
}

void call(gmfsetnormalatquadrilateralvertex)(int *MshIdx, int *i0, int *i1, int *i2, int *i3)
{
	GmfSetLin(*MshIdx, GmfNormalAtQuadrilateralVertices, *i0, *i1, *i2, *i3);
}

void call(gmfgettrianglep2)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6)
{
	GmfGetLin(*MshIdx, GmfTrianglesP2, i0, i1, i2, i3, i4, i5, i6);
}

void call(gmfsettrianglep2)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6)
{
	GmfSetLin(*MshIdx, GmfTrianglesP2, *i0, *i1, *i2, *i3, *i4, *i5, *i6);
}

void call(gmfgetedgep2)(int *MshIdx, int *i0, int *i1, int *i2, int *i3)
{
	GmfGetLin(*MshIdx, GmfEdgesP2, i0, i1, i2, i3);
}

void call(gmfsetedgep2)(int *MshIdx, int *i0, int *i1, int *i2, int *i3)
{
	GmfSetLin(*MshIdx, GmfEdgesP2, *i0, *i1, *i2, *i3);
}

void call(gmfgetsolatpyramidr4)(int *MshIdx, float *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtPyramids, SolTab);
}

void call(gmfgetsolatpyramidr8)(int *MshIdx, double *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtPyramids, SolTab);
}

void call(gmfsetsolatpyramidr4)(int *MshIdx, float *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtPyramids, SolTab);
}

void call(gmfsetsolatpyramidr8)(int *MshIdx, double *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtPyramids, SolTab);
}

void call(gmfgetquadrilateralq2)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8, int *i9)
{
	GmfGetLin(*MshIdx, GmfQuadrilateralsQ2, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9);
}

void call(gmfsetquadrilateralq2)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8, int *i9)
{
	GmfSetLin(*MshIdx, GmfQuadrilateralsQ2, *i0, *i1, *i2, *i3, *i4, *i5, *i6, *i7, *i8, *i9);
}

void call(gmfgetisolatpyramid)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4)
{
	GmfGetLin(*MshIdx, GmfISolAtPyramids, i0, i1, i2, i3, i4);
}

void call(gmfsetisolatpyramid)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4)
{
	GmfSetLin(*MshIdx, GmfISolAtPyramids, *i0, *i1, *i2, *i3, *i4);
}

void call(gmfgetsubdomainfromgeom)(int *MshIdx, int *i0, int *i1, int *i2)
{
	GmfGetLin(*MshIdx, GmfSubDomainFromGeom, i0, i1, i2);
}

void call(gmfsetsubdomainfromgeom)(int *MshIdx, int *i0, int *i1, int *i2)
{
	GmfSetLin(*MshIdx, GmfSubDomainFromGeom, *i0, *i1, *i2);
}

void call(gmfgettetrahedronp2)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8, int *i9, int *i10)
{
	GmfGetLin(*MshIdx, GmfTetrahedraP2, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10);
}

void call(gmfsettetrahedronp2)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8, int *i9, int *i10)
{
	GmfSetLin(*MshIdx, GmfTetrahedraP2, *i0, *i1, *i2, *i3, *i4, *i5, *i6, *i7, *i8, *i9, *i10);
}

void call(gmfgetfault_neartri)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfFault_NearTri, i0);
}

void call(gmfsetfault_neartri)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfFault_NearTri, *i0);
}

void call(gmfgetfault_inter)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfFault_Inter, i0);
}

void call(gmfsetfault_inter)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfFault_Inter, *i0);
}

void call(gmfgethexahedronq2)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8, int *i9, int *i10, int *i11, int *i12, int *i13, int *i14, int *i15, int *i16, int *i17, int *i18, int *i19, int *i20, int *i21, int *i22, int *i23, int *i24, int *i25, int *i26, int *i27)
{
	GmfGetLin(*MshIdx, GmfHexahedraQ2, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23, i24, i25, i26, i27);
}

void call(gmfsethexahedronq2)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8, int *i9, int *i10, int *i11, int *i12, int *i13, int *i14, int *i15, int *i16, int *i17, int *i18, int *i19, int *i20, int *i21, int *i22, int *i23, int *i24, int *i25, int *i26, int *i27)
{
	GmfSetLin(*MshIdx, GmfHexahedraQ2, *i0, *i1, *i2, *i3, *i4, *i5, *i6, *i7, *i8, *i9, *i10, *i11, *i12, *i13, *i14, *i15, *i16, *i17, *i18, *i19, *i20, *i21, *i22, *i23, *i24, *i25, *i26, *i27);
}

void call(gmfgetextraverticesatedge)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfExtraVerticesAtEdges, i0);
}

void call(gmfsetextraverticesatedge)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfExtraVerticesAtEdges, *i0);
}

void call(gmfgetextraverticesattriangle)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfExtraVerticesAtTriangles, i0);
}

void call(gmfsetextraverticesattriangle)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfExtraVerticesAtTriangles, *i0);
}

void call(gmfgetextraverticesatquadrilateral)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfExtraVerticesAtQuadrilaterals, i0);
}

void call(gmfsetextraverticesatquadrilateral)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfExtraVerticesAtQuadrilaterals, *i0);
}

void call(gmfgetextraverticesattetrahedron)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfExtraVerticesAtTetrahedra, i0);
}

void call(gmfsetextraverticesattetrahedron)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfExtraVerticesAtTetrahedra, *i0);
}

void call(gmfgetextraverticesatprism)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfExtraVerticesAtPrisms, i0);
}

void call(gmfsetextraverticesatprism)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfExtraVerticesAtPrisms, *i0);
}

void call(gmfgetextraverticesathexahedron)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfExtraVerticesAtHexahedra, i0);
}

void call(gmfsetextraverticesathexahedron)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfExtraVerticesAtHexahedra, *i0);
}

void call(gmfgetvertexongeometricvertexr4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfGetLin(*MshIdx, GmfVerticesOnGeometricVertices, i0, i1, r2);
}

void call(gmfgetvertexongeometricvertexr8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfGetLin(*MshIdx, GmfVerticesOnGeometricVertices, i0, i1, r2);
}

void call(gmfsetvertexongeometricvertexr4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfSetLin(*MshIdx, GmfVerticesOnGeometricVertices, *i0, *i1, *r2);
}

void call(gmfsetvertexongeometricvertexr8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfSetLin(*MshIdx, GmfVerticesOnGeometricVertices, *i0, *i1, *r2);
}

void call(gmfgetvertexongeometricedger4)(int *MshIdx, int *i0, int *i1, float *r2, float *r3)
{
	GmfGetLin(*MshIdx, GmfVerticesOnGeometricEdges, i0, i1, r2, r3);
}

void call(gmfgetvertexongeometricedger8)(int *MshIdx, int *i0, int *i1, double *r2, double *r3)
{
	GmfGetLin(*MshIdx, GmfVerticesOnGeometricEdges, i0, i1, r2, r3);
}

void call(gmfsetvertexongeometricedger4)(int *MshIdx, int *i0, int *i1, float *r2, float *r3)
{
	GmfSetLin(*MshIdx, GmfVerticesOnGeometricEdges, *i0, *i1, *r2, *r3);
}

void call(gmfsetvertexongeometricedger8)(int *MshIdx, int *i0, int *i1, double *r2, double *r3)
{
	GmfSetLin(*MshIdx, GmfVerticesOnGeometricEdges, *i0, *i1, *r2, *r3);
}

void call(gmfgetvertexongeometrictriangler4)(int *MshIdx, int *i0, int *i1, float *r2, float *r3, float *r4)
{
	GmfGetLin(*MshIdx, GmfVerticesOnGeometricTriangles, i0, i1, r2, r3, r4);
}

void call(gmfgetvertexongeometrictriangler8)(int *MshIdx, int *i0, int *i1, double *r2, double *r3, double *r4)
{
	GmfGetLin(*MshIdx, GmfVerticesOnGeometricTriangles, i0, i1, r2, r3, r4);
}

void call(gmfsetvertexongeometrictriangler4)(int *MshIdx, int *i0, int *i1, float *r2, float *r3, float *r4)
{
	GmfSetLin(*MshIdx, GmfVerticesOnGeometricTriangles, *i0, *i1, *r2, *r3, *r4);
}

void call(gmfsetvertexongeometrictriangler8)(int *MshIdx, int *i0, int *i1, double *r2, double *r3, double *r4)
{
	GmfSetLin(*MshIdx, GmfVerticesOnGeometricTriangles, *i0, *i1, *r2, *r3, *r4);
}

void call(gmfgetvertexongeometricquadrilateralr4)(int *MshIdx, int *i0, int *i1, float *r2, float *r3, float *r4)
{
	GmfGetLin(*MshIdx, GmfVerticesOnGeometricQuadrilaterals, i0, i1, r2, r3, r4);
}

void call(gmfgetvertexongeometricquadrilateralr8)(int *MshIdx, int *i0, int *i1, double *r2, double *r3, double *r4)
{
	GmfGetLin(*MshIdx, GmfVerticesOnGeometricQuadrilaterals, i0, i1, r2, r3, r4);
}

void call(gmfsetvertexongeometricquadrilateralr4)(int *MshIdx, int *i0, int *i1, float *r2, float *r3, float *r4)
{
	GmfSetLin(*MshIdx, GmfVerticesOnGeometricQuadrilaterals, *i0, *i1, *r2, *r3, *r4);
}

void call(gmfsetvertexongeometricquadrilateralr8)(int *MshIdx, int *i0, int *i1, double *r2, double *r3, double *r4)
{
	GmfSetLin(*MshIdx, GmfVerticesOnGeometricQuadrilaterals, *i0, *i1, *r2, *r3, *r4);
}

void call(gmfgetedgeongeometricedger4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfGetLin(*MshIdx, GmfEdgesOnGeometricEdges, i0, i1, r2);
}

void call(gmfgetedgeongeometricedger8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfGetLin(*MshIdx, GmfEdgesOnGeometricEdges, i0, i1, r2);
}

void call(gmfsetedgeongeometricedger4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfSetLin(*MshIdx, GmfEdgesOnGeometricEdges, *i0, *i1, *r2);
}

void call(gmfsetedgeongeometricedger8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfSetLin(*MshIdx, GmfEdgesOnGeometricEdges, *i0, *i1, *r2);
}

void call(gmfgetfault_freeedge)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfFault_FreeEdge, i0);
}

void call(gmfsetfault_freeedge)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfFault_FreeEdge, *i0);
}

void call(gmfgetpolyhedron)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8, int *i9, int *i10, int *i11, int *i12, int *i13, int *i14, int *i15, int *i16, int *i17, int *i18, int *i19, int *i20, int *i21, int *i22, int *i23, int *i24, int *i25, int *i26, int *i27, int *i28, int *i29, int *i30, int *i31, int *i32)
{
	GmfGetLin(*MshIdx, GmfPolyhedra, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23, i24, i25, i26, i27, i28, i29, i30, i31, i32);
}

void call(gmfsetpolyhedron)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8, int *i9, int *i10, int *i11, int *i12, int *i13, int *i14, int *i15, int *i16, int *i17, int *i18, int *i19, int *i20, int *i21, int *i22, int *i23, int *i24, int *i25, int *i26, int *i27, int *i28, int *i29, int *i30, int *i31, int *i32)
{
	GmfSetLin(*MshIdx, GmfPolyhedra, *i0, *i1, *i2, *i3, *i4, *i5, *i6, *i7, *i8, *i9, *i10, *i11, *i12, *i13, *i14, *i15, *i16, *i17, *i18, *i19, *i20, *i21, *i22, *i23, *i24, *i25, *i26, *i27, *i28, *i29, *i30, *i31, *i32);
}

void call(gmfgetfault_overlap)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfFault_Overlap, i0);
}

void call(gmfsetfault_overlap)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfFault_Overlap, *i0);
}

void call(gmfgetpyramid)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5)
{
	GmfGetLin(*MshIdx, GmfPyramids, i0, i1, i2, i3, i4, i5);
}

void call(gmfsetpyramid)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5)
{
	GmfSetLin(*MshIdx, GmfPyramids, *i0, *i1, *i2, *i3, *i4, *i5);
}

void call(gmfgetprivatetable)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfPrivateTable, i0);
}

void call(gmfsetprivatetable)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfPrivateTable, *i0);
}

void call(gmfgetfault_badshape)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfFault_BadShape, i0);
}

void call(gmfsetfault_badshape)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfFault_BadShape, *i0);
}

void call(gmfgettriangleongeometrictriangler4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfGetLin(*MshIdx, GmfTrianglesOnGeometricTriangles, i0, i1, r2);
}

void call(gmfgettriangleongeometrictriangler8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfGetLin(*MshIdx, GmfTrianglesOnGeometricTriangles, i0, i1, r2);
}

void call(gmfsettriangleongeometrictriangler4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfSetLin(*MshIdx, GmfTrianglesOnGeometricTriangles, *i0, *i1, *r2);
}

void call(gmfsettriangleongeometrictriangler8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfSetLin(*MshIdx, GmfTrianglesOnGeometricTriangles, *i0, *i1, *r2);
}

void call(gmfgettriangleongeometricquadrilateralr4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfGetLin(*MshIdx, GmfTrianglesOnGeometricQuadrilaterals, i0, i1, r2);
}

void call(gmfgettriangleongeometricquadrilateralr8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfGetLin(*MshIdx, GmfTrianglesOnGeometricQuadrilaterals, i0, i1, r2);
}

void call(gmfsettriangleongeometricquadrilateralr4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfSetLin(*MshIdx, GmfTrianglesOnGeometricQuadrilaterals, *i0, *i1, *r2);
}

void call(gmfsettriangleongeometricquadrilateralr8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfSetLin(*MshIdx, GmfTrianglesOnGeometricQuadrilaterals, *i0, *i1, *r2);
}

void call(gmfgetquadrilateralongeometrictriangler4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfGetLin(*MshIdx, GmfQuadrilateralsOnGeometricTriangles, i0, i1, r2);
}

void call(gmfgetquadrilateralongeometrictriangler8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfGetLin(*MshIdx, GmfQuadrilateralsOnGeometricTriangles, i0, i1, r2);
}

void call(gmfsetquadrilateralongeometrictriangler4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfSetLin(*MshIdx, GmfQuadrilateralsOnGeometricTriangles, *i0, *i1, *r2);
}

void call(gmfsetquadrilateralongeometrictriangler8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfSetLin(*MshIdx, GmfQuadrilateralsOnGeometricTriangles, *i0, *i1, *r2);
}

void call(gmfgetquadrilateralongeometricquadrilateralr4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfGetLin(*MshIdx, GmfQuadrilateralsOnGeometricQuadrilaterals, i0, i1, r2);
}

void call(gmfgetquadrilateralongeometricquadrilateralr8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfGetLin(*MshIdx, GmfQuadrilateralsOnGeometricQuadrilaterals, i0, i1, r2);
}

void call(gmfsetquadrilateralongeometricquadrilateralr4)(int *MshIdx, int *i0, int *i1, float *r2)
{
	GmfSetLin(*MshIdx, GmfQuadrilateralsOnGeometricQuadrilaterals, *i0, *i1, *r2);
}

void call(gmfsetquadrilateralongeometricquadrilateralr8)(int *MshIdx, int *i0, int *i1, double *r2)
{
	GmfSetLin(*MshIdx, GmfQuadrilateralsOnGeometricQuadrilaterals, *i0, *i1, *r2);
}

void call(gmfgettangent2dr4)(int *MshIdx, float *r0, float *r1)
{
	GmfGetLin(*MshIdx, GmfTangents, r0, r1);
}

void call(gmfgettangent2dr8)(int *MshIdx, double *r0, double *r1)
{
	GmfGetLin(*MshIdx, GmfTangents, r0, r1);
}

void call(gmfgettangent3dr4)(int *MshIdx, float *r0, float *r1, float *r2)
{
	GmfGetLin(*MshIdx, GmfTangents, r0, r1, r2);
}

void call(gmfgettangent3dr8)(int *MshIdx, double *r0, double *r1, double *r2)
{
	GmfGetLin(*MshIdx, GmfTangents, r0, r1, r2);
}

void call(gmfsettangent2dr4)(int *MshIdx, float *r0, float *r1)
{
	GmfSetLin(*MshIdx, GmfTangents, *r0, *r1);
}

void call(gmfsettangent2dr8)(int *MshIdx, double *r0, double *r1)
{
	GmfSetLin(*MshIdx, GmfTangents, *r0, *r1);
}

void call(gmfsettangent3dr4)(int *MshIdx, float *r0, float *r1, float *r2)
{
	GmfSetLin(*MshIdx, GmfTangents, *r0, *r1, *r2);
}

void call(gmfsettangent3dr8)(int *MshIdx, double *r0, double *r1, double *r2)
{
	GmfSetLin(*MshIdx, GmfTangents, *r0, *r1, *r2);
}

void call(gmfgetnormal2dr4)(int *MshIdx, float *r0, float *r1)
{
	GmfGetLin(*MshIdx, GmfNormals, r0, r1);
}

void call(gmfgetnormal2dr8)(int *MshIdx, double *r0, double *r1)
{
	GmfGetLin(*MshIdx, GmfNormals, r0, r1);
}

void call(gmfgetnormal3dr4)(int *MshIdx, float *r0, float *r1, float *r2)
{
	GmfGetLin(*MshIdx, GmfNormals, r0, r1, r2);
}

void call(gmfgetnormal3dr8)(int *MshIdx, double *r0, double *r1, double *r2)
{
	GmfGetLin(*MshIdx, GmfNormals, r0, r1, r2);
}

void call(gmfsetnormal2dr4)(int *MshIdx, float *r0, float *r1)
{
	GmfSetLin(*MshIdx, GmfNormals, *r0, *r1);
}

void call(gmfsetnormal2dr8)(int *MshIdx, double *r0, double *r1)
{
	GmfSetLin(*MshIdx, GmfNormals, *r0, *r1);
}

void call(gmfsetnormal3dr4)(int *MshIdx, float *r0, float *r1, float *r2)
{
	GmfSetLin(*MshIdx, GmfNormals, *r0, *r1, *r2);
}

void call(gmfsetnormal3dr8)(int *MshIdx, double *r0, double *r1, double *r2)
{
	GmfSetLin(*MshIdx, GmfNormals, *r0, *r1, *r2);
}

void call(gmfgettangentatvertex)(int *MshIdx, int *i0, int *i1)
{
	GmfGetLin(*MshIdx, GmfTangentAtVertices, i0, i1);
}

void call(gmfsettangentatvertex)(int *MshIdx, int *i0, int *i1)
{
	GmfSetLin(*MshIdx, GmfTangentAtVertices, *i0, *i1);
}

void call(gmfgetsolatvertexr4)(int *MshIdx, float *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtVertices, SolTab);
}

void call(gmfgetsolatvertexr8)(int *MshIdx, double *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtVertices, SolTab);
}

void call(gmfsetsolatvertexr4)(int *MshIdx, float *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtVertices, SolTab);
}

void call(gmfsetsolatvertexr8)(int *MshIdx, double *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtVertices, SolTab);
}

void call(gmfgetsolatedger4)(int *MshIdx, float *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtEdges, SolTab);
}

void call(gmfgetsolatedger8)(int *MshIdx, double *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtEdges, SolTab);
}

void call(gmfsetsolatedger4)(int *MshIdx, float *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtEdges, SolTab);
}

void call(gmfsetsolatedger8)(int *MshIdx, double *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtEdges, SolTab);
}

void call(gmfgetsolattriangler4)(int *MshIdx, float *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtTriangles, SolTab);
}

void call(gmfgetsolattriangler8)(int *MshIdx, double *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtTriangles, SolTab);
}

void call(gmfsetsolattriangler4)(int *MshIdx, float *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtTriangles, SolTab);
}

void call(gmfsetsolattriangler8)(int *MshIdx, double *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtTriangles, SolTab);
}

void call(gmfgetsolatquadrilateralr4)(int *MshIdx, float *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtQuadrilaterals, SolTab);
}

void call(gmfgetsolatquadrilateralr8)(int *MshIdx, double *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtQuadrilaterals, SolTab);
}

void call(gmfsetsolatquadrilateralr4)(int *MshIdx, float *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtQuadrilaterals, SolTab);
}

void call(gmfsetsolatquadrilateralr8)(int *MshIdx, double *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtQuadrilaterals, SolTab);
}

void call(gmfgetsolattetrahedronr4)(int *MshIdx, float *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtTetrahedra, SolTab);
}

void call(gmfgetsolattetrahedronr8)(int *MshIdx, double *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtTetrahedra, SolTab);
}

void call(gmfsetsolattetrahedronr4)(int *MshIdx, float *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtTetrahedra, SolTab);
}

void call(gmfsetsolattetrahedronr8)(int *MshIdx, double *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtTetrahedra, SolTab);
}

void call(gmfgetsolatprismr4)(int *MshIdx, float *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtPrisms, SolTab);
}

void call(gmfgetsolatprismr8)(int *MshIdx, double *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtPrisms, SolTab);
}

void call(gmfsetsolatprismr4)(int *MshIdx, float *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtPrisms, SolTab);
}

void call(gmfsetsolatprismr8)(int *MshIdx, double *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtPrisms, SolTab);
}

void call(gmfgetsolathexahedronr4)(int *MshIdx, float *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtHexahedra, SolTab);
}

void call(gmfgetsolathexahedronr8)(int *MshIdx, double *SolTab)
{
	GmfGetLin(*MshIdx, GmfSolAtHexahedra, SolTab);
}

void call(gmfsetsolathexahedronr4)(int *MshIdx, float *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtHexahedra, SolTab);
}

void call(gmfsetsolathexahedronr8)(int *MshIdx, double *SolTab)
{
	GmfSetLin(*MshIdx, GmfSolAtHexahedra, SolTab);
}

void call(gmfgetdsolatvertexr4)(int *MshIdx, float *SolTab)
{
	GmfGetLin(*MshIdx, GmfDSolAtVertices, SolTab);
}

void call(gmfgetdsolatvertexr8)(int *MshIdx, double *SolTab)
{
	GmfGetLin(*MshIdx, GmfDSolAtVertices, SolTab);
}

void call(gmfsetdsolatvertexr4)(int *MshIdx, float *SolTab)
{
	GmfSetLin(*MshIdx, GmfDSolAtVertices, SolTab);
}

void call(gmfsetdsolatvertexr8)(int *MshIdx, double *SolTab)
{
	GmfSetLin(*MshIdx, GmfDSolAtVertices, SolTab);
}

void call(gmfgetisolatvertex)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfISolAtVertices, i0);
}

void call(gmfsetisolatvertex)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfISolAtVertices, *i0);
}

void call(gmfgetisolatedge)(int *MshIdx, int *i0, int *i1)
{
	GmfGetLin(*MshIdx, GmfISolAtEdges, i0, i1);
}

void call(gmfsetisolatedge)(int *MshIdx, int *i0, int *i1)
{
	GmfSetLin(*MshIdx, GmfISolAtEdges, *i0, *i1);
}

void call(gmfgetisolattriangle)(int *MshIdx, int *i0, int *i1, int *i2)
{
	GmfGetLin(*MshIdx, GmfISolAtTriangles, i0, i1, i2);
}

void call(gmfsetisolattriangle)(int *MshIdx, int *i0, int *i1, int *i2)
{
	GmfSetLin(*MshIdx, GmfISolAtTriangles, *i0, *i1, *i2);
}

void call(gmfgetisolatquadrilateral)(int *MshIdx, int *i0, int *i1, int *i2, int *i3)
{
	GmfGetLin(*MshIdx, GmfISolAtQuadrilaterals, i0, i1, i2, i3);
}

void call(gmfsetisolatquadrilateral)(int *MshIdx, int *i0, int *i1, int *i2, int *i3)
{
	GmfSetLin(*MshIdx, GmfISolAtQuadrilaterals, *i0, *i1, *i2, *i3);
}

void call(gmfgetisolattetrahedron)(int *MshIdx, int *i0, int *i1, int *i2, int *i3)
{
	GmfGetLin(*MshIdx, GmfISolAtTetrahedra, i0, i1, i2, i3);
}

void call(gmfsetisolattetrahedron)(int *MshIdx, int *i0, int *i1, int *i2, int *i3)
{
	GmfSetLin(*MshIdx, GmfISolAtTetrahedra, *i0, *i1, *i2, *i3);
}

void call(gmfgetisolatprism)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5)
{
	GmfGetLin(*MshIdx, GmfISolAtPrisms, i0, i1, i2, i3, i4, i5);
}

void call(gmfsetisolatprism)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5)
{
	GmfSetLin(*MshIdx, GmfISolAtPrisms, *i0, *i1, *i2, *i3, *i4, *i5);
}

void call(gmfgetisolathexahedron)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7)
{
	GmfGetLin(*MshIdx, GmfISolAtHexahedra, i0, i1, i2, i3, i4, i5, i6, i7);
}

void call(gmfsetisolathexahedron)(int *MshIdx, int *i0, int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7)
{
	GmfSetLin(*MshIdx, GmfISolAtHexahedra, *i0, *i1, *i2, *i3, *i4, *i5, *i6, *i7);
}

void call(gmfgetfault_smalltri)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfFault_SmallTri, i0);
}

void call(gmfsetfault_smalltri)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfFault_SmallTri, *i0);
}

void call(gmfgetcoarsehexahedron)(int *MshIdx, int *i0)
{
	GmfGetLin(*MshIdx, GmfCoarseHexahedra, i0);
}

void call(gmfsetcoarsehexahedron)(int *MshIdx, int *i0)
{
	GmfSetLin(*MshIdx, GmfCoarseHexahedra, *i0);
}

void call(gmfgetcomment)(int *MshIdx)
{
	GmfGetLin(*MshIdx, GmfComments);
}

void call(gmfsetcomment)(int *MshIdx)
{
	GmfSetLin(*MshIdx, GmfComments);
}

