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

  _nTri         .resize (_nPart);
  _eltsConnecTri.resize (_nPart);
  _eltsGnumTri  .resize (_nPart);
  
  _nQuad         .resize (_nPart);
  _eltsConnecQuad.resize (_nPart);
  _eltsGnumQuad  .resize (_nPart);

  _nElts   .resize (_nPart);  
  _eltsConnecIndex.resize (_nPart);
  _eltsConnec     .resize (_nPart);   
  _eltsGnum       .resize (_nPart);
  
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




