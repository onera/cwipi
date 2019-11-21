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
 
    bool operator()(const Point& p) const
    {
        if (a.y > b.y) return Edge{ b, a }(p);
        if (p.y == a.y || p.y == b.y) return operator()({ p.x, p.y + epsilon });
        if (p.y > b.y || p.y < a.y || p.x > max(a.x, b.x)) return false;
        if (p.x < min(a.x, b.x)) return true;
        auto blue = abs(a.x - p.x) > MIN ? (p.y - a.y) / (p.x - a.x) : MAX;
        auto red = abs(a.x - b.x) > MIN ? (b.y - a.y) / (b.x - a.x) : MAX;
        return blue >= red;
    }
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
                  :_nx(0),_ny(0),_prop(1.0),_color(0),_width(0.0),_randomVar(1.0)
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


Point transform (Point p){
  return { 3.0 * (p.x - 5.0)/10.0 + 0.295, 3.0 * (p.y - 5.0)/10.0 + 0.244 } ; 
}

double surfMeshGenerator::_inBox2(double x,double y) {

  double xr = x/_width;
  double yr = y/_width;

  Figure square = { "Square",
       {  {{0.0, 0.0}, {0.2, 0.0}}, {{0.2, 0.0}, {0.2, 0.2}}, {{0.2, 0.2}, {0.0, 0.2}}, {{0.0, 0.2}, {0.0, 0.0}} }
  };
  
  std::vector<Point> path = {
  {1.19,4.77},{1.45,5.46},{2.05,5.69},{2.39,5.58},{2.50,5.47},{2.64,5.66},{2.69,5.70},{2.69,4.99},{2.67,4.93},
  {2.62,4.97},{2.52,5.31},{2.26,5.58},{1.95,5.61},{1.59,5.37},{1.43,4.86},{1.47,4.33},{1.76,3.92},{2.07,3.82},
  {2.41,3.95},{2.59,4.20},{2.63,4.44},{2.67,4.47},{2.70,4.42},{2.50,3.96},{2.29,3.79},{2.04,3.74},{1.64,3.85},
  {1.36,4.11},{1.21,4.39},{1.19,4.77} };
  
  std::vector<Edge> path2;
  for (int i=0; i<path.size()-1; i++){
     Point p = {0.0,0.0};
     Edge e = { {1.19,4.77},{1.45,5.46} };
     
     path[0].x = 3.0 * ( (path[0].x - 5.0)/10.0 + 0.295 ) ; 
     path[0].y = 3.0 * ( (path[0].y - 5.0)/10.0 + 0.0244 ) ; 
     
     if(i+1 != path.size()){
       path[i+1].x = 3.0 * ( (path[i+1].x - 5.0)/10.0 + 0.295 ) ; 
       path[i+1].y = 3.0 * ( (path[i+1].y - 5.0)/10.0 + 0.0244 ) ; 
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
   int nBox = box_lim.size()/4;
   
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
    double zc = _centerTri[i_part][3*i_tri+2];
    
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
    double zc = _centerQuad[i_part][3*i_Quad+2];
    
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
    double zc = _centerPoly[i_part][3*i_poly+2];
    
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

    int ppartId;
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

    PDM_part_create (&ppartId,
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

    
    int id_mn = PDM_Mesh_nodal_create(_nPart, _interfComm);
    for(int i_part =0;i_part<_nPart;i_part++){
    
      int nFace;
      int nEdge;
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
                           &nFace,
                           &nEdge,
                           &nEdgePartBound,
                           &nVtx1,
                           &nProc,
                           &nTPart,
                           &sFaceEdge,
                           &sEdgeVtx,
                           &sEdgeGroup,
                           &nEdgeGroup2);
      
      int          *faceTag;
      int          *faceEdgeIdx;
      int          *faceEdge;
      PDM_g_num_t *faceLNToGN;
      int          *edgeTag;
      int          *edgeFace;
      int          *edgeVtxIdx;
      int          *edgeVtx;
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
                           &faceEdgeIdx,
                           &faceEdge,
                           &faceLNToGN,
                           &edgeTag,
                           &edgeFace,
                           &edgeVtxIdx,
                           &edgeVtx,
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

      int *edgeVtxN = (int*)malloc(sizeof(int) * nEdge);
      for (int i = 0; i < nEdge; i++) {
        edgeVtxN[i] = edgeVtxIdx[i+1] - edgeVtxIdx[i];
      }

      int *faceEdgeN = (int*)malloc(sizeof(int) * nFace);
      for (int i = 0; i < nFace; i++) {
        faceEdgeN[i] = faceEdgeIdx[i+1] - faceEdgeIdx[i];
      }


      PDM_Mesh_nodal_coord_set(id_mn,
                               i_part,
                               nVtx1,
                               _coords[i_part],
                               _vtxGnum[i_part]);


      PDM_Mesh_nodal_cell2d_celledge_add (id_mn,
                                        i_part,
                                        nFace,
                                        nEdge,
                                        edgeVtxIdx,
                                        edgeVtxN,
                                        edgeVtx,
                                        faceEdgeIdx,
                                        faceEdgeN,
                                        faceEdge,
                                        faceLNToGN);

      _nElts[i_part] = nFace;
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
          _eltsGnumTri[i_part] = PDM_Mesh_nodal_block_g_num_get (id_mn, id_block, i_part);      
          /*
          for(int i =0; i<_nTri;i++){
            for (int j=0;j<3;j++){
             printf("eltsConnecTri[%i] %i rank %i nVtx %i _nTri %i color %i n_block %i\n",3*i+j,_eltsConnecTri[3*i+j],_rank,_nVtx,_nTri,_color,n_block);
            }
          }        
          */
        }
      }
      else if(t_block == PDM_MESH_NODAL_QUAD4) {
        for(int i_part =0;i_part<_nPart;i_part++){
          _nQuad[i_part] = PDM_Mesh_nodal_block_n_elt_get (id_mn, id_block, i_part);            
          PDM_Mesh_nodal_block_std_get (id_mn, id_block, i_part, &_eltsConnecQuad[i_part]);
          _eltsGnumQuad[i_part] = PDM_Mesh_nodal_block_g_num_get (id_mn, id_block, i_part);      
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
          
          _eltsGnumPoly[i_part]  = PDM_Mesh_nodal_block_g_num_get (id_mn, id_block, i_part);    
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




