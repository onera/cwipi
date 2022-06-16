/*
  This file is part of the CWIPI library.

  Copyright (C) 2011-2017  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more detailstr_options.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/license_width/>.
*/


/**
 * \cond
 */


#include "pdm_timer.h"
#include "pdm.h"
#include "pdm_config.h"
#include "time.h"
#include "surfMeshGenerator.hxx"
#include "surfMeshGeneratorDB.hxx"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <pdm_geom_elem.h>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace std;

const double epsilon = numeric_limits<float>().epsilon();
const numeric_limits<double> DOUBLE;
const double MIN = DOUBLE.min();
const double MAX = DOUBLE.max();

struct Point { double x, y; };

struct Edge {
    Point a, b;

    CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
    bool operator()(const Point& p) const
    {
        if (a.y > b.y) return Edge{ b, a }(p);
        if (p.y == a.y || p.y == b.y) return operator()({ p.x, p.y + epsilon });
        if (p.y > b.y || p.y < a.y || p.x > max(a.x, b.x)) return false;
        if (p.x < min(a.x, b.x)) return true;
        auto blue = abs((long) (a.x - p.x)) > MIN ? (p.y - a.y) / (p.x - a.x) : MAX;
        auto red = abs((long) (a.x - b.x)) > MIN ? (b.y - a.y) / (b.x - a.x) : MAX;
        return blue >= red;
    }
    CWP_GCC_SUPPRESS_WARNING_POP
};

struct Figure {
    const string  name;
    vector<Edge> edges;

    bool contains(Point& p)
    {
        auto c = 0;
        for (auto e : edges) if (e(p)) c++;
        return c % 2 != 0;
    }

};





namespace cwipi {



surfMeshGenerator::surfMeshGenerator()
                  :_nx(0),_ny(0),_color(0),_prop(1.0),_width(0.0),_randomVar(1.0)
{

}


void surfMeshGenerator::init(int nx, int ny, int nPart, MPI_Comm* comm, double prop, double width, double randomVar)
{

  _nx = nx;
  _ny = ny;
  _prop = prop;
  _width = width;
  _comm = _interfComm = PDM_MPI_mpi_2_pdm_mpi_comm(&comm);
  _randomVar = randomVar;
  _nPart = nPart;

  _nVtx    .resize (_nPart);
  _coords  .resize (_nPart);
  _vtxGnum .resize (_nPart);

  _nPoly              .resize (_nPart);
  _eltsConnecPolyIndex.resize (_nPart);
  _eltsConnecPoly     .resize (_nPart);
  _eltsGnumPoly       .resize (_nPart);
  _surfVecPoly   .resize (_nPart,NULL);
  _centerPoly    .resize (_nPart,NULL);
  _charLengthPoly.resize (_nPart,NULL);
  _isDegPoly     .resize (_nPart,NULL);

  _nTri         .resize (_nPart);
  _eltsConnecTri.resize (_nPart);
  _eltsGnumTri  .resize (_nPart);
  _surfVecTri   .resize (_nPart);
  _centerTri    .resize (_nPart);
  _charLengthTri.resize (_nPart);
  _isDegTri     .resize (_nPart);


  _nQuad         .resize (_nPart);
  _eltsConnecQuad.resize (_nPart);
  _eltsGnumQuad  .resize (_nPart);
  _surfVecQuad   .resize (_nPart);
  _centerQuad    .resize (_nPart);
  _charLengthQuad.resize (_nPart);
  _isDegQuad     .resize (_nPart);


  _nFace      .resize (_nPart);
  _nEdge      .resize (_nPart);
  _faceEdgeIdx.resize (_nPart);
  _faceEdge   .resize (_nPart);
  _edgeVtxIdx .resize (_nPart);
  _edgeVtx    .resize (_nPart);
  _faceLNToGN .resize (_nPart);
  _nElts   .resize (_nPart);
  _eltsConnecIndex.resize (_nPart);
  _eltsConnec     .resize (_nPart);
  _eltsGnum       .resize (_nPart);

  _specialFieldTri.resize(_nPart);
  _specialFieldQuad.resize(_nPart);
  _specialFieldPoly.resize(_nPart);

  /* Interface communicator
   * ---------------------- */

  MPI_Comm_rank(*comm, &_rank);
  MPI_Comm_size(*comm, &_commSize);

  if (_rank < _prop * _commSize) {
    _color = 1;
  }
  if (_rank == 0) {
    _color = 1;
  }

  MPI_Comm interfComm = MPI_COMM_NULL;
  MPI_Comm_split(*comm, _color, _rank, &interfComm);
  MPI_Comm_size(interfComm, &_interfCommSize);
  _interfComm = PDM_MPI_mpi_2_pdm_mpi_comm(&interfComm);

  _xmin = -_width/2.;
  _xmax =  _width/2.;
  _ymin = -_width/2.;
  _ymax =  _width/2.;

}


surfMeshGenerator::~surfMeshGenerator() {

}


double surfMeshGenerator::_inBox(double x, double y, double x1, double y1 , double x2, double y2) {

  double xr = x/_width;
  double yr = y/_width;

  if( x1< xr && xr<x2 && yr<y2 && y1<yr)
    return 1.;
  else
    return 0.5;
}


Point _transform0 (Point p){
  return { 3.0 * ((p.x - 5.0)/10.0 + 0.295), 3.0 * ((p.y - 5.0)/10.0 + 0.0244) } ;
}


Point _transform (Point p){
  Point center = {-165.,-26.};
  Point size = {7.,7.};
  return { 0.8 * ((p.x - center.x)/size.x ) -0.183,  -0.8 * ((p.y - center.y)/size.y ) - 0.0414 } ;
}



double surfMeshGenerator::_inBox2(double x,double y) {

  double xr = x/_width;
  double yr = y/_width;

  Figure square = { "Square",
       {  {{0.0, 0.0}, {0.2, 0.0}}, {{0.2, 0.0}, {0.2, 0.2}}, {{0.2, 0.2}, {0.0, 0.2}}, {{0.0, 0.2}, {0.0, 0.0}} }
  };


std::vector<Point> path = {
{-166.22656,-26.4375},{-166.220597128,-26.6790365708},{-166.19357045,-26.9201755029},{-166.146217463,-27.1573567309},{-166.079288711,-27.389159316},
{-165.992940712,-27.6154946081},{-165.888212905,-27.833561281},{-165.765799333,-28.0419162912},{-165.62593096,-28.2397014179},{-165.469951007,-28.4244658927},
{-165.29811,-28.595306},{-165.10881387,-28.7628296233},{-164.904684394,-28.9111590685},{-164.687625132,-29.0395521856},{-164.459877506,-29.1468927885},
{-164.223626805,-29.2320648381},{-163.980511843,-29.2940401983},{-163.731691903,-29.3316020717},{-163.480084342,-29.3431032592},{-163.22914985,-29.3272449178},
{-162.98038,-29.282822},{-162.865343868,-29.2476699695},{-162.75231701,-29.2055055526},{-162.642475303,-29.1569071568},{-162.535974442,-29.1021654297},
{-162.431771388,-29.0408092591},{-162.332466826,-28.9745516278},{-162.236404804,-28.9025177792},{-162.143578485,-28.8246580343},{-162.055964837,-28.7427741346},
{-161.97266,-28.65625},{-161.926360529,-28.7223008755},{-161.880547106,-28.7894053503},{-161.835320567,-28.8567603178},{-161.790375705,-28.9241420914},
{-161.744876547,-28.992124386},{-161.699587454,-29.0588938685},{-161.653395481,-29.1254352186},{-161.60585664,-29.191691871},{-161.557606422,-29.2561464207},
{-161.50781,-29.319339},{-161.462640405,-29.3354845807},{-161.415639267,-29.3363017009},{-161.378455871,-29.3089206952},{-161.360919197,-29.2653780165},
{-161.35430791,-29.2197468753},{-161.352664132,-29.1719764988},{-161.353592942,-29.1255826877},{-161.355500696,-29.0790809238},{-161.357006357,-29.0308159508},
{-161.35551,-28.984198},{-161.35586447,-28.8014388541},{-161.355830273,-28.6175136017},{-161.355555458,-28.4339844983},{-161.355191113,-28.2508107361},
{-161.354886464,-28.0679368528},{-161.354789469,-27.884925056},{-161.355051634,-27.7010252033},{-161.355820882,-27.517856793},{-161.357257325,-27.333628186},
{-161.35951,-27.149421},{-161.379865498,-27.1040305718},{-161.415202513,-27.0694847934},{-161.463123475,-27.0564748559},{-161.511422956,-27.0654493774},
{-161.554232917,-27.0887746265},{-161.590679806,-27.122395948},{-161.618307055,-27.1630013956},{-161.63547985,-27.2090329095},{-161.639751868,-27.2580695506},
{-161.62993,-27.305781},{-161.671780695,-27.5199680661},{-161.731325714,-27.7316806121},{-161.807858054,-27.9368410182},{-161.901902761,-28.1346811413},
{-162.01422618,-28.323412918},{-162.144396687,-28.4988966619},{-162.294340113,-28.6599056032},{-162.462054921,-28.8006693021},{-162.647151957,-28.9181416531},
{-162.84697,-29.008662},{-163.082326751,-29.0692397986},{-163.323294426,-29.0930914193},{-163.566143856,-29.0817773153},{-163.80490967,-29.0375623613},
{-164.03476895,-28.9635556608},{-164.254665729,-28.8615656194},{-164.460300415,-28.7342248987},{-164.650007911,-28.5823668901},{-164.819745972,-28.4078614475},
{-164.96471,-28.21377},{-165.09120349,-27.9735644512},{-165.195157731,-27.7226700602},{-165.27667927,-27.4647371229},{-165.337396604,-27.1995042154},
{-165.377889236,-26.9308282622},{-165.399476785,-26.6610153236},{-165.40356399,-26.388628767},{-165.391208471,-26.1176247272},{-165.363370151,-25.8469268189},
{-165.32121,-25.57972},{-165.28084362,-25.4127343518},{-165.232680034,-25.2470435946},{-165.175871662,-25.0838906109},{-165.109650105,-24.9248499274},
{-165.03324047,-24.7714687438},{-164.94475098,-24.6235919572},{-164.844055494,-24.4844794931},{-164.730677371,-24.3561924916},{-164.604071924,-24.2404409813},
{-164.46531,-24.13981},{-164.313204307,-24.0417560622},{-164.151617732,-23.9572766421},{-163.983393276,-23.8879025837},{-163.810054456,-23.8342497835},
{-163.633087293,-23.7969387754},{-163.451814693,-23.7764668735},{-163.27097657,-23.7740110297},{-163.089382624,-23.7902860766},{-162.912071185,-23.825594358},
{-162.73843,-23.880734},{-162.551579501,-23.9626839975},{-162.375474462,-24.0663627213},{-162.213096213,-24.1890696805},{-162.065733563,-24.3287694393},
{-161.933942773,-24.4844306804},{-161.819438793,-24.6541459373},{-161.724774329,-24.8340531274},{-161.650334353,-25.0237770187},{-161.597860642,-25.220437219},
{-161.56851,-25.422639},{-161.574863935,-25.4687686609},{-161.57366386,-25.5143246731},{-161.566791824,-25.5592766867},{-161.554525895,-25.6031754993},
{-161.535609985,-25.6457992586},{-161.507876882,-25.6830029342},{-161.468353285,-25.7042874169},{-161.424454253,-25.694589486},{-161.38991829,-25.6647567782},
{-161.36182,-25.626213},{-161.364248842,-25.4026219036},{-161.395417122,-25.1813174464},{-161.453688058,-24.9644457068},{-161.536437659,-24.7560288132},
{-161.640620029,-24.5589551161},{-161.76495526,-24.3729232034},{-161.906859788,-24.2004531321},{-162.065705008,-24.0418560959},{-162.239583274,-23.899420935},
{-162.42548,-23.776054},{-162.634928614,-23.6782030879},{-162.854215425,-23.6048056849},{-163.079401204,-23.5557098759},{-163.309952821,-23.5297293291},
{-163.541752266,-23.5261939951},{-163.772071501,-23.5438583743},{-164.001361518,-23.581953182},{-164.225492566,-23.6392551289},{-164.444517634,-23.7154562915},
{-164.65496,-23.809344},{-164.910990499,-23.9588799819},{-165.150204544,-24.1342917983},{-165.3699668,-24.3326262857},{-165.568521254,-24.5514337567},
{-165.745089271,-24.789775856},{-165.896647699,-25.0440337321},{-166.021380029,-25.311819654},{-166.117713975,-25.5916316275},{-166.183291016,-25.879731855},
{-166.21651,-26.174572},{-166.218492488,-26.2012310059},{-166.22014095,-26.2270714681},{-166.221571436,-26.2534064501},{-166.222760143,-26.279348918},{-166.223800755,-26.3066029717},
{-166.224639268,-26.3332669881},{-166.225314169,-26.3593950481},{-166.225883442,-26.3863008709},{-166.226331393,-26.4119332545},{-166.22671,-26.437503},
{-166.22656,-26.4375},{-166.22656,-26.4375},{-166.22656,-26.4375},{-166.22656,-26.4375},{-166.22656,-26.4375},{-166.22656,-26.4375},{-166.22656,-26.4375},
{-166.22656,-26.4375},{-166.22656,-26.4375},{-166.22656,-26.4375} };

  std::vector<Point> path0 = {
  {1.19,4.77},{1.45,5.46},{2.05,5.69},{2.39,5.58},{2.50,5.47},{2.64,5.66},{2.69,5.70},{2.69,4.99},{2.67,4.93},
  {2.62,4.97},{2.52,5.31},{2.26,5.58},{1.95,5.61},{1.59,5.37},{1.43,4.86},{1.47,4.33},{1.76,3.92},{2.07,3.82},
  {2.41,3.95},{2.59,4.20},{2.63,4.44},{2.67,4.47},{2.70,4.42},{2.50,3.96},{2.29,3.79},{2.04,3.74},{1.64,3.85},
  {1.36,4.11},{1.21,4.39},{1.19,4.77} };

  std::vector<Edge> path2;
  for (size_t i=0; i<path.size()-1; i++){
     Edge e = { {1.19,4.77},{1.45,5.46} };

     //path[0].x = 3.0 * ( (path[0].x - 5.0)/10.0 + 0.295 ) ;
     //path[0].y = 3.0 * ( (path[0].y - 5.0)/10.0 + 0.0244 ) ;
     path[0] = _transform (path[0]);


     if(i+1 != path.size()){
       //path[i+1].x = 3.0 * ( (path[i+1].x - 5.0)/10.0 + 0.295 ) ;
       //path[i+1].y = 3.0 * ( (path[i+1].y - 5.0)/10.0 + 0.0244 ) ;
       path[i+1] = _transform (path[i+1]);
       e = {path[i],path[i+1]};
     }
     else {
       e = {path[i],path[0]};
     }
     //printf("edge %3.2f,%3.2f  %3.2f,%3.2f\n",e.a.x,e.a.y,e.b.x,e.b.y);
     path2.push_back(e);
  }

  Figure cell = { "C",path2};

  Point p = {xr,yr};


  if(cell.contains(p))
    return p.y;
  else
    return 0.0;
}


double surfMeshGenerator::_inCircle(double x, double y, double R) {

  double xr = x/_width;
  double yr = y/_width;

  if(xr*xr + yr*yr < R*R)
    return 1.;
  else
    return 0.;
}


double surfMeshGenerator::_motif(double x, double y) {


   std::vector<double> box_lim = {0.1,0.1, 0.4,0.4};
   double result = 0.0;

   result = _inBox2(x,y);//,  x1,y1, x2,y2 );

/*
   for(int i=0; i<nBox; i++){
     double x1 = box_lim[4*i];
     double y1 = box_lim[4*i+1];
     double x2 = box_lim[4*i+2];
     double y2 = box_lim[4*i+3];
     result = _inBox2(x,y);//,  x1,y1, x2,y2 );
   }
  */
   return result;
}

double* surfMeshGenerator::specialFieldTriGet(int i_part) {

  _specialFieldTri[i_part] = (double*)malloc(sizeof(double) * _nTri[i_part]);
  _surfVecTri     [i_part] = (double*)malloc(sizeof(double) * 3 * _nTri[i_part]);
  _centerTri      [i_part] = (double*)malloc(sizeof(double) * 3 * _nTri[i_part]);
  _charLengthTri  [i_part] = (double*)malloc(sizeof(double) * _nTri[i_part]);
  _isDegTri       [i_part] = (int*)   malloc(sizeof(int) * _nTri[i_part]);

  PDM_geom_elem_tria_properties(_nTri[i_part],
                                _eltsConnecTri[i_part],
                                _coords[i_part],
                                _surfVecTri[i_part],
                                _centerTri[i_part],
                                _charLengthTri[i_part],
                                _isDegTri[i_part]);

  for(int i_tri=0; i_tri < _nTri[i_part]; i_tri++){

    double xc = _centerTri[i_part][3*i_tri]  ;
    double yc = _centerTri[i_part][3*i_tri+1];

    _specialFieldTri[i_part][i_tri] = _motif(xc,yc);

  }

  return _specialFieldTri[i_part];
}


double* surfMeshGenerator::specialFieldQuadGet(int i_part) {

  _specialFieldQuad[i_part] = (double*)malloc(sizeof(double) * _nQuad[i_part]);
  _surfVecQuad     [i_part] = (double*)malloc(sizeof(double) * 3 * _nQuad[i_part]);
  _centerQuad      [i_part] = (double*)malloc(sizeof(double) * 3 * _nQuad[i_part]);
  _charLengthQuad  [i_part] = (double*)malloc(sizeof(double) * _nQuad[i_part]);
  _isDegQuad       [i_part] = (int*)   malloc(sizeof(int) * _nQuad[i_part]);

  PDM_geom_elem_quad_properties(_nQuad[i_part],
                                _eltsConnecQuad[i_part],
                                _coords[i_part],
                                _surfVecQuad[i_part],
                                _centerQuad[i_part],
                                _charLengthQuad[i_part],
                                _isDegQuad[i_part]);

  for(int i_Quad=0; i_Quad < _nQuad[i_part]; i_Quad++){

    double xc = _centerQuad[i_part][3*i_Quad]  ;
    double yc = _centerQuad[i_part][3*i_Quad+1];

    _specialFieldQuad[i_part][i_Quad] = _motif(xc,yc);

  }

  return _specialFieldQuad[i_part];
}


double* surfMeshGenerator::specialFieldPolyGet(int i_part) {

  _specialFieldPoly[i_part] = (double*)malloc(sizeof(double) * _nPoly[i_part]);
  _surfVecPoly     [i_part] = (double*)malloc(sizeof(double) * 3 * _nPoly[i_part]);
  _centerPoly      [i_part] = (double*)malloc(sizeof(double) * 3 * _nPoly[i_part]);
  _charLengthPoly  [i_part] = (double*)malloc(sizeof(double) * _nPoly[i_part]);
  _isDegPoly       [i_part] = (int*)   malloc(sizeof(int) * _nPoly[i_part]);

  PDM_geom_elem_polygon_properties(_nPoly[i_part],
                                   _eltsConnecPolyIndex[i_part],
                                   _eltsConnecPoly[i_part],
                                   _coords[i_part],
                                   _surfVecPoly[i_part],
                                   _centerPoly[i_part],
                                   _charLengthPoly[i_part],
                                   _isDegPoly[i_part]);

  for(int i_poly=0; i_poly < _nPoly[i_part]; i_poly++){

    double xc = _centerPoly[i_part][3*i_poly]  ;
    double yc = _centerPoly[i_part][3*i_poly+1];

    _specialFieldPoly[i_part][i_poly] = _motif(xc,yc);

  }
  return _specialFieldPoly[i_part];
}




void surfMeshGenerator::computeMesh() {

  /* Define mesh in interface communicator
   * ------------------------------------- */

  if (_color == 1) {

    const int haveRandom = 1;
    const int initRandom = time(NULL) * _randomVar;
    PDM_g_num_t  nGFace;
    PDM_g_num_t  nGVtx;
    PDM_g_num_t  nGEdge;
    int         d_nVtx;
    double     *dVtxCoord;
    int         dNFace;
    int        *dFaceVtxIdx;
    PDM_g_num_t *dFaceVtx;
    PDM_g_num_t *dFaceEdge;
    int         dNEdge;
    PDM_g_num_t *dEdgeVtx;
    PDM_g_num_t *dEdgeFace;
    int         nEdgeGroup;
    int        *dEdgeGroupIdx;
    PDM_g_num_t *dEdgeGroup;

    PDM_poly_surf_gen (_interfComm,
                       _xmin,
                       _xmax,
                       _ymin,
                       _ymax,
                       haveRandom,
                       initRandom,
                       _nx,
                       _ny,
                       &nGFace,
                       &nGVtx,
                       &nGEdge,
                       &d_nVtx,
                       &dVtxCoord,
                       &dNFace,
                       &dFaceVtxIdx,
                       &dFaceVtx,
                       &dFaceEdge,
                       &dNEdge,
                       &dEdgeVtx,
                       &dEdgeFace,
                       &nEdgeGroup,
                       &dEdgeGroupIdx,
                       &dEdgeGroup);

    PDM_part_t *ppartId;
#ifdef PDM_HAVE_PARMETIS
    PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
    PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#endif
#endif
    int have_dCellPart = 0;

    int *dCellPart = (int *) malloc(dNFace*sizeof(int));
    int *renum_properties_cell = NULL;
    int *renum_properties_face = NULL;
    int nPropertyCell = 0;
    int nPropertyFace = 0;

    int *dEdgeVtxIdx = (int*)malloc (sizeof(int) * (dNEdge+1));
    dEdgeVtxIdx[0] = 0;
    for (int i = 0; i < dNEdge; i++) {
      dEdgeVtxIdx[i+1] = dEdgeVtxIdx[i] + 2;
    }

    ppartId = PDM_part_create (
                     _interfComm,
                     method,
                     "PDM_PART_RENUM_CELL_NONE",
                     "PDM_PART_RENUM_FACE_NONE",
                     nPropertyCell,
                     renum_properties_cell,
                     nPropertyFace,
                     renum_properties_face,
                     _nPart,
                     dNFace,
                     dNEdge,
                     d_nVtx,
                     nEdgeGroup,
                     NULL,
                     NULL,
                     NULL,
                     NULL,
                     have_dCellPart,
                     dCellPart,
                     dEdgeFace,
                     dEdgeVtxIdx,
                     dEdgeVtx,
                     NULL,
                     dVtxCoord,
                     NULL,
                     dEdgeGroupIdx,
                     dEdgeGroup);

    free (dCellPart);


    PDM_Mesh_nodal_t *id_mn = PDM_Mesh_nodal_create(_nPart, _interfComm);
    for(int i_part =0;i_part<_nPart;i_part++){

      int nEdgePartBound;
      int nVtx1;
      int nProc;
      int sFaceEdge;
      int sEdgeVtx;
      int sEdgeGroup;
      int nEdgeGroup2;
      int nTPart;

      PDM_part_part_dim_get (ppartId,
                           i_part,
                           &(_nFace[i_part]),
                           &(_nEdge[i_part]),
                           &nEdgePartBound,
                           &nVtx1,
                           &nProc,
                           &nTPart,
                           &sFaceEdge,
                           &sEdgeVtx,
                           &sEdgeGroup,
                           &nEdgeGroup2);

      int          *faceTag;
      int          *edgeTag;
      int          *edgeFace;
      PDM_g_num_t *edgeLNToGN;
      int          *edgePartBoundProcIdx;
      int          *edgePartBoundPartIdx;
      int          *edgePartBound;
      int          *vtxTag;
      double       *vtx;
      PDM_g_num_t *vtxLNToGN;
      int          *edgeGroupIdx;
      int          *edgeGroup;
      PDM_g_num_t *edgeGroupLNToGN;

      PDM_part_part_val_get (ppartId,
                           i_part,
                           &faceTag,
                           &(_faceEdgeIdx[i_part]),
                           &(_faceEdge[i_part]),
                           &(_faceLNToGN[i_part]),
                           &edgeTag,
                           &edgeFace,
                           &(_edgeVtxIdx[i_part]),
                           &(_edgeVtx[i_part]),
                           &edgeLNToGN,
                           &edgePartBoundProcIdx,
                           &edgePartBoundPartIdx,
                           &edgePartBound,
                           &vtxTag,
                           &vtx,
                           &vtxLNToGN,
                           &edgeGroupIdx,
                           &edgeGroup,
                           &edgeGroupLNToGN);

      _nVtx[i_part] = nVtx1;
      _coords[i_part] = (double*)malloc (sizeof(double) * 3 * nVtx1);
      for (int i = 0; i < 3 * nVtx1; i++) {
        _coords[i_part][i] = vtx[i];
      }

      _vtxGnum[i_part] = (CWP_g_num_t*)malloc (sizeof(CWP_g_num_t) * nVtx1);
      for (int i = 0; i <  nVtx1; i++) {
        _vtxGnum[i_part][i] = (CWP_g_num_t) vtxLNToGN[i];
      }

      int *edgeVtxN = (int*)malloc(sizeof(int) * _nEdge[i_part]);
      for (int i = 0; i < _nEdge[i_part]; i++) {
        edgeVtxN[i] = _edgeVtxIdx[i_part][i+1] - _edgeVtxIdx[i_part][i];
      }

      int *faceEdgeN = (int*)malloc(sizeof(int) * _nFace[i_part]);
      for (int i = 0; i < _nFace[i_part]; i++) {
        faceEdgeN[i] = _faceEdgeIdx[i_part][i+1] - _faceEdgeIdx[i_part][i];
      }


      PDM_Mesh_nodal_coord_set(id_mn,
                               i_part,
                               nVtx1,
                               _coords[i_part],
                               _vtxGnum[i_part],
                               PDM_OWNERSHIP_USER);


      PDM_Mesh_nodal_cell2d_celledge_add (id_mn,
                                        i_part,
                                        _nFace[i_part],
                                        _nEdge[i_part],
                                        _edgeVtxIdx[i_part],
                                        edgeVtxN,
                                        _edgeVtx[i_part],
                                        _faceEdgeIdx[i_part],
                                        faceEdgeN,
                                        _faceEdge[i_part],
                                        _faceLNToGN[i_part],
                                        PDM_OWNERSHIP_USER);

      _nElts[i_part] = _nFace[i_part];
    }


    int n_block = PDM_Mesh_nodal_n_blocks_get (id_mn);
    assert (n_block == 3);

    int *block_ids = PDM_Mesh_nodal_blocks_id_get (id_mn);

    assert (PDM_Mesh_nodal_block_type_get (id_mn, block_ids[0]) == PDM_MESH_NODAL_TRIA3);
    assert (PDM_Mesh_nodal_block_type_get (id_mn, block_ids[1]) == PDM_MESH_NODAL_QUAD4);
    assert (PDM_Mesh_nodal_block_type_get (id_mn, block_ids[2]) == PDM_MESH_NODAL_POLY_2D);



    for (int i = 0; i < n_block; i++) {
      int id_block = block_ids[i];

      PDM_Mesh_nodal_elt_t t_block = PDM_Mesh_nodal_block_type_get (id_mn, block_ids[i]);

      if (t_block == PDM_MESH_NODAL_TRIA3){
        for(int i_part =0;i_part<_nPart;i_part++){
          _nTri[i_part] = PDM_Mesh_nodal_block_n_elt_get (id_mn, id_block, i_part);
          PDM_Mesh_nodal_block_std_get (id_mn, id_block, i_part, &_eltsConnecTri[i_part]);
          _eltsGnumTri[i_part] = PDM_Mesh_nodal_g_num_get (id_mn, id_block, i_part);

          /*
          for(int i =0; i<_nTri[i_part]; i++){
            for (int j=0;j<3;j++){
             printf("eltsConnecTri[%i] %i rank %i nVtx %i _nTri %i color %i n_block %i\n",3*i+j,_eltsConnecTri[i_part][3*i+j],_rank,_nVtx[i_part],_nTri[i_part],_color,n_block);
            }
          }
          */
        }
      }
      else if(t_block == PDM_MESH_NODAL_QUAD4) {
        for(int i_part =0;i_part<_nPart;i_part++){
          _nQuad[i_part] = PDM_Mesh_nodal_block_n_elt_get (id_mn, id_block, i_part);
          PDM_Mesh_nodal_block_std_get (id_mn, id_block, i_part, &_eltsConnecQuad[i_part]);
          _eltsGnumQuad[i_part] = PDM_Mesh_nodal_g_num_get (id_mn, id_block, i_part);
          /*
          for(int i =0; i<_nTri;i++){
            for (int j=0;j<3;j++){
             printf("eltsConnecTri[%i] %i rank %i nVtx %i _nTri %i color %i n_block %i\n",3*i+j,_eltsConnecTri[3*i+j],_rank,_nVtx,_nTri,_color,n_block);
            }
          }
          */
        }

      }
      else if(t_block == PDM_MESH_NODAL_POLY_2D){
        for(int i_part =0;i_part<_nPart;i_part++){
          _nPoly[i_part]  = PDM_Mesh_nodal_block_n_elt_get (id_mn, id_block, i_part);
          PDM_Mesh_nodal_block_poly2d_get (id_mn, id_block, i_part, &_eltsConnecPolyIndex[i_part] , &_eltsConnecPoly[i_part] );

          if(_nPoly[i_part] == 0){
            _eltsConnecPolyIndex[i_part] = (int*) malloc(sizeof(int));
            _eltsConnecPolyIndex[i_part][0]=0;
          }

          _eltsGnumPoly[i_part]  = PDM_Mesh_nodal_g_num_get (id_mn, id_block, i_part);
         /* for(int i =0; i<_nPoly;i++){
            for (int j=_eltsConnecPolyIndex[i];j<_eltsConnecPolyIndex[i+1];j++){
              printf("_eltsConnecPoly[%i] %i rank %i nVtx %i _nPoly %i color %i n_block %i\n",j,_eltsConnecPoly[j],_rank,_nVtx,_nPoly,_color,n_block);
            }
          }
         */
        }
      }
    }

    for(int i_part =0;i_part<_nPart;i_part++){
      assert(_nElts[i_part] == _nTri[i_part] + _nQuad[i_part] + _nPoly[i_part]);

      _eltsConnecIndex[i_part] = (int*) malloc( sizeof(int) * (_nElts[i_part]+1) );
      _eltsConnecIndex[i_part][0] = 0;

      int idx =0;
      for(int i =0; i<_nTri[i_part];i++){
        _eltsConnecIndex[i_part][idx+1] = _eltsConnecIndex[i_part][idx] + 3;
        idx++;
      }

      for(int i =0; i< _nQuad[i_part];i++){
        _eltsConnecIndex[i_part][idx+1] = _eltsConnecIndex[i_part][idx] + 4;
        idx++;
      }

      for(int i =0; i< _nPoly[i_part]; i++){
        _eltsConnecIndex[i_part][idx+1] = _eltsConnecIndex[i_part][idx] + (_eltsConnecPolyIndex[i_part][i+1]-_eltsConnecPolyIndex[i_part][i]);
        idx++;
      }

      _nElts[i_part] = _nPoly[i_part] + _nTri[i_part] + _nQuad[i_part];

      _eltsConnec[i_part] = (int*) malloc(sizeof(int) * _eltsConnecIndex[i_part][_nElts[i_part]] );
      memcpy( _eltsConnec[i_part], _eltsConnecTri[i_part], sizeof(int)*3*_nTri[i_part] );
      memcpy( &(_eltsConnec[i_part][3*_nTri[i_part]]), _eltsConnecQuad[i_part], sizeof(int)*4*_nQuad[i_part] );

      if(_nPoly[i_part]>0)
        memcpy( &(_eltsConnec[i_part][ 3*_nTri[i_part]+4*_nQuad[i_part] ]), _eltsConnecPoly[i_part], sizeof(int)*_eltsConnecPolyIndex[i_part][ _nPoly[i_part] ] );
      /*
      for(int i =0; i<_nElts;i++){
        for (int j=_eltsConnecIndex[i];j<_eltsConnecIndex[i+1];j++){
          printf("_eltsConnec[%i] %i rank %i Elts %i nVtx %i _nElts %i color %i n_block %i\n",j,_eltsConnec[j],_rank,i,_nVtx,_nElts,_color,n_block);
        }
      }
     */

   }//end for i_part

     //PDM_Mesh_nodal_free (id_mn);
    //PDM_part_free (ppartId);

   /* free (dEdgeVtxIdx);
    free (edgeVtxN);
    free (faceEdgeN);
*/
    //PDM_MPI_Comm_free(&_interfComm);
  }
  else {

    for(int i_part =0;i_part<_nPart;i_part++){
      _nVtx[i_part] = 0;
      _nElts[i_part] = 0;
      _coords[i_part] = (double*)malloc (sizeof(double) * 3 * _nVtx[i_part]);
      _vtxGnum[i_part] = (CWP_g_num_t*)malloc (sizeof(CWP_g_num_t) * _nVtx[i_part]);

      _eltsConnecTri[i_part] = (int*)malloc(sizeof(int) * (_nElts[i_part]));
      _eltsConnecQuad[i_part] = (int*)malloc(sizeof(int) * (_nElts[i_part]));
      _eltsConnecPolyIndex[i_part] = (int*)malloc(sizeof(int) * (_nElts[i_part]+1));
      _eltsConnecPolyIndex[i_part][0] = 0;
      _eltsConnecPoly[i_part] = (int*)malloc(sizeof(int) * _eltsConnecPolyIndex[i_part][_nElts[i_part]]);
      _eltsGnumTri[i_part] = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t) * _nElts[i_part]);
      _eltsGnumQuad[i_part] = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t) * _nElts[i_part]);
      _eltsGnumPoly[i_part] = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t) * _nElts[i_part]);
    }
  }





}

}//end namespace cwipi


/**
 * \endcond
 */
