/*
  This file is part of the CWIPI library.

  Copyright (C) 2012  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <map>
#include <mesh.hxx>
#include <mpi.h>
#include <cstdlib>

#include <pdm_mesh_nodal.h>
#include <pdm_gnum.h>
#include <pdm_error.h>
#include <pdm_printf.h>
#include "cwp.h"
#include "factory.hpp"
#include "block.hxx"
#include "blockStd.hxx"
#include "blockCP.hxx"
#include "blockFP.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"

/**
 * \cond
 */

namespace cwipi {

  /**
   * \typedef FB
   *
   * \brief Block Factory
   *
   *  A block \ref Factory wich makes \ref Block
   *  class objects.
   *  The type of block objects build depends on the
   *  block type \ref CWP_Block_t .
   *
   */

  typedef Factory<Block, CWP_Block_t> FB;

  /**
    * \brief Mesh constructor
    *
    * Construct the CWIPI mesh by using paradigm nodal methods.
    *
    * \param [in] npart Number of mesh partitions.
    *
    */

  Mesh::Mesh
  (
    const MPI_Comm &localComm,
    Visu* visu,
    const int npart,
    const CWP_Dynamic_mesh_t displacement,
    Coupling *cpl
  )
  :
    _localComm(localComm),
    _nBlocks(0),
                  //_hoOrdering (NULL),
    _visu(visu),
    _displacement(displacement),
    _cpl(cpl),
    _faceEdgeMethod(0),
    _cellFaceMethod(0),
    _pdmNodal_handle_index(),
    _isVtxGnumComputed(false),
    _isEltGnumComputed(false)
  {

    _pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&localComm));
    // pdm_nodal building

    _npart                 = npart;
    _nVertex   .resize(npart,0);
    _nElts     .resize(npart,0);
    _connec_idx.resize(npart,NULL);
    _connec    .resize(npart,NULL);
    _gnum_elt  .resize(npart,NULL);
    _elt_centers  .resize(npart,NULL);

    _coords .resize(npart,NULL);
    _global_num_vtx .resize(npart,NULL);
    _global_num_elt .resize(npart,NULL);

    _nCells      .resize(npart, 0)  ;
    _cellFaceIdx .resize(npart, NULL);
    _cellFace    .resize(npart, NULL);
    _nFace       .resize(npart,0)   ;
    _faceEdgeIdx .resize(npart,NULL);
    _faceEdge    .resize(npart,NULL);
    _faceVtxIdx  .resize(npart,NULL);
    _faceVtx     .resize(npart,NULL);
    _nEdge       .resize(npart,0)   ;
    _edgeVtxIdx  .resize(npart,NULL);
    _edgeVtx     .resize(npart,NULL);

    _cellFaceNb  .resize(npart,NULL);
    _edgeVtxNb   .resize(npart,NULL);
    _faceEdgeNb  .resize(npart,NULL);
    _faceVtxNb   .resize(npart,NULL);

    _faceLNToGN  .resize(npart,NULL);
    _cellLNToGN  .resize(npart,NULL);

  }

  Mesh::~Mesh()
  {
#if defined(DEBUG) && 0
    cout << "destroying mesh of partition  : TODO" << endl;
#endif
    for (int i = 0; i < _npart; i++) {
      free (_connec_idx[i]);
      free (_connec[i]);
      free (_gnum_elt[i]);
    }

    for (int i = 0; i < (int) _blockDB.size(); i++) {
      delete _blockDB[i];
    }

  }



  int
  Mesh::_Mesh_nodal_block_std_type_size_get
  (
    CWP_Block_t type
  )
  {

    switch (type) {

      case CWP_BLOCK_NODE: return 1;
      break;

      case CWP_BLOCK_EDGE2: return 2;
      break;

      case CWP_BLOCK_FACE_TRIA3: return 3;
      break;

      case CWP_BLOCK_FACE_QUAD4: return 4;
      break;

      case CWP_BLOCK_CELL_TETRA4: return 4;
      break;

      case CWP_BLOCK_CELL_HEXA8: return 8;
      break;

      case CWP_BLOCK_CELL_PYRAM5: return 5;
      break;

      case CWP_BLOCK_CELL_PRISM6: return 6;
      break;

      default: return -1;
      PDM_error(__FILE__, __LINE__, 0, "This CWP_Block_t is not available as function argument.\n");

    }
  }

  void
  Mesh::eltCentersCompute
  (
    int i_part
  )
  {

    int n_elt_part = getPartNElts(i_part);

    if(_elt_centers[i_part] == NULL)
      _elt_centers[i_part] = (double*)malloc(sizeof(double)*3* n_elt_part);


    int ind_idx=0;
    for(int i=0;i<_nBlocks;i++){
      int n_elt = _blockDB[i]->NEltsGet()[i_part];

      const double* elt_centers_block = _blockDB[i]->eltCentersGet(i_part);

      for(int j=0;j<n_elt;j++){
        for(int k=0;k<3;k++){
          _elt_centers[i_part][3*ind_idx+k]=elt_centers_block[3*j+k];
        }
        ind_idx++;
      }//end loop j
    }//end loop on block
  }


  CWP_g_num_t*
  Mesh::GNumEltsGet
  (
    int i_part
  )
  {
    // if (_connec_idx[i_part]==NULL) {
    //   connecCompute(i_part);
    // }//end if NULL

    return PDM_Mesh_nodal_g_num_get_from_part (_pdmNodal_handle_index, i_part);
  }


  double*
  Mesh::eltCentersGet
  (
    int i_part
  )
  {
    if(_elt_centers[i_part]==NULL || _displacement != CWP_DYNAMIC_MESH_STATIC) {
      eltCentersCompute(i_part);
    }//end if NULL

    return _elt_centers[i_part];
  }


  void
  Mesh::coordSet
  (
    const int   i_part,
    const int   n_vtx,
    double      coords[],
    CWP_g_num_t global_num[]
  )
  {
    _coords[i_part] = coords;
    _nVertex[i_part]  = n_vtx;
    _global_num_vtx[i_part] = global_num;
  }



  /**************************************************************/


  void
  Mesh::updateBlockDB()
  {
    int n_block = PDM_Mesh_nodal_n_blocks_get (_pdmNodal_handle_index);
    int *block_ids = PDM_Mesh_nodal_blocks_id_get (_pdmNodal_handle_index);
    for (int i = 0; i < n_block; i++) {
      PDM_Mesh_nodal_elt_t t_block = PDM_Mesh_nodal_block_type_get (_pdmNodal_handle_index, block_ids[i]);

      if (t_block == PDM_MESH_NODAL_TRIA3){
        int block_id = blockAdd(CWP_BLOCK_FACE_TRIA3);
        for(int i_part =0;i_part<_npart;i_part++){
          int n_tri = PDM_Mesh_nodal_block_n_elt_get (_pdmNodal_handle_index,  block_ids[i], i_part);
          int* connec = NULL;
          CWP_g_num_t* gnum = NULL;
          PDM_Mesh_nodal_block_std_get (_pdmNodal_handle_index, block_ids[i], i_part, &connec);
          gnum = PDM_Mesh_nodal_g_num_get (_pdmNodal_handle_index, block_ids[i], i_part);

          stdBlockSet (i_part  ,
                       block_id,
                       n_tri  ,
                       connec ,
                       gnum);

        }
      }

      else if (t_block == PDM_MESH_NODAL_QUAD4) {
        int block_id = blockAdd(CWP_BLOCK_FACE_QUAD4);
        for (int i_part =0;i_part<_npart;i_part++){
          int n_quad = PDM_Mesh_nodal_block_n_elt_get (_pdmNodal_handle_index, block_ids[i], i_part);
          int* connec = NULL;
          CWP_g_num_t* gnum = NULL;
          PDM_Mesh_nodal_block_std_get (_pdmNodal_handle_index, block_ids[i], i_part, &connec);
          gnum = PDM_Mesh_nodal_g_num_get (_pdmNodal_handle_index, block_ids[i], i_part);

          stdBlockSet (i_part  ,
                       block_id,
                       n_quad  ,
                       connec ,
                       gnum   );

        }
      }

      else if(t_block == PDM_MESH_NODAL_POLY_2D){
        int block_id = blockAdd(CWP_BLOCK_FACE_POLY);
        for(int i_part =0;i_part<_npart;i_part++){
          int n_poly = PDM_Mesh_nodal_block_n_elt_get (_pdmNodal_handle_index, block_ids[i], i_part);
          int* connec = NULL;
          int* connec_idx = NULL;
          CWP_g_num_t* gnum = NULL;
          PDM_Mesh_nodal_block_poly2d_get (_pdmNodal_handle_index, block_ids[i], i_part, &connec_idx , &connec );
          gnum = PDM_Mesh_nodal_g_num_get (_pdmNodal_handle_index, block_ids[i], i_part);

          poly2DBlockSet (i_part  ,
                          block_id,
                          n_poly  ,
                          connec_idx,
                          connec ,
                          gnum   );
        }
      }

      else if (t_block == PDM_MESH_NODAL_PYRAMID5) {
        int block_id = blockAdd(CWP_BLOCK_CELL_PYRAM5);
        for (int i_part = 0 ; i_part < _npart ; i_part++) {
          int n_pyramid = PDM_Mesh_nodal_block_n_elt_get(_pdmNodal_handle_index, block_ids[i], i_part);
          int *connec = NULL;
          CWP_g_num_t *gnum;
          PDM_Mesh_nodal_block_std_get(_pdmNodal_handle_index, block_ids[i], i_part, &connec);
          gnum = PDM_Mesh_nodal_g_num_get(_pdmNodal_handle_index, block_ids[i], i_part);

          stdBlockSet(i_part, block_id, n_pyramid, connec, gnum);
        }
      }

      else if (t_block == PDM_MESH_NODAL_PRISM6) {
        int block_id = blockAdd(CWP_BLOCK_CELL_PRISM6);
        for (int i_part = 0 ; i_part < _npart ; i_part++) {
          int n_prism = PDM_Mesh_nodal_block_n_elt_get(_pdmNodal_handle_index, block_ids[i], i_part);
          int *connec = NULL;
          CWP_g_num_t *gnum;
          PDM_Mesh_nodal_block_std_get(_pdmNodal_handle_index, block_ids[i], i_part, &connec);
          gnum = PDM_Mesh_nodal_g_num_get(_pdmNodal_handle_index, block_ids[i], i_part);

          stdBlockSet(i_part, block_id, n_prism, connec, gnum);
        }
      }

      else if (t_block == PDM_MESH_NODAL_HEXA8) {
        int block_id = blockAdd(CWP_BLOCK_CELL_HEXA8);
        for (int i_part = 0 ; i_part < _npart ; i_part++) {
          int n_hexa = PDM_Mesh_nodal_block_n_elt_get(_pdmNodal_handle_index, block_ids[i], i_part);
          int *connec = NULL;
          CWP_g_num_t *gnum;
          PDM_Mesh_nodal_block_std_get(_pdmNodal_handle_index, block_ids[i], i_part, &connec);
          gnum = PDM_Mesh_nodal_g_num_get(_pdmNodal_handle_index, block_ids[i], i_part);

          stdBlockSet(i_part, block_id, n_hexa, connec, gnum);
        }
      }

      else if (t_block == PDM_MESH_NODAL_TETRA4) {
        int block_id = blockAdd(CWP_BLOCK_CELL_TETRA4);
        for (int i_part = 0 ; i_part < _npart ; i_part++) {
          int n_tetra = PDM_Mesh_nodal_block_n_elt_get(_pdmNodal_handle_index, block_ids[i], i_part);
          int *connec = NULL;
          CWP_g_num_t *gnum;
          PDM_Mesh_nodal_block_std_get(_pdmNodal_handle_index, block_ids[i], i_part, &connec);
          gnum = PDM_Mesh_nodal_g_num_get(_pdmNodal_handle_index, block_ids[i], i_part);

          stdBlockSet(i_part, block_id, n_tetra, connec, gnum);
        }
      }
       // Define all the other types
      else if(t_block == PDM_MESH_NODAL_POLY_3D){
        int block_id = blockAdd(CWP_BLOCK_CELL_POLY);
        for(int i_part =0;i_part<_npart;i_part++){
          int n_poly = PDM_Mesh_nodal_block_n_elt_get (_pdmNodal_handle_index, block_ids[i], i_part);
          PDM_l_num_t       n_face;
          PDM_l_num_t      *facvtx_idx;
          PDM_l_num_t      *facvtx;
          PDM_l_num_t      *cellfac_idx;
          PDM_l_num_t      *cellfa;
          CWP_g_num_t* gnum = NULL;

          PDM_Mesh_nodal_block_poly3d_get (_pdmNodal_handle_index, block_ids[i], i_part, &n_face, &facvtx_idx, &facvtx, &cellfac_idx, &cellfa);

          gnum = PDM_Mesh_nodal_g_num_get (_pdmNodal_handle_index, block_ids[i], i_part);

          poly3DBlockSet (i_part,
                          block_id,
                          n_poly,
                          n_face,
                          facvtx_idx,
                          facvtx,
                          cellfac_idx,
                          cellfa,
                          gnum);
        }
      }
    }
  }



  /**********************************************************************/

  void
  Mesh::geomFinalize
  (
  )
  {

    int unionRank;
    MPI_Comm_rank(_cpl->communicationGet()->unionCommGet(),&unionRank);

    int globalRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&globalRank);

    _pdmNodal_handle_index = PDM_Mesh_nodal_create (_npart,_pdm_localComm);

    if(gnumVtxRequired () ){
      _isVtxGnumComputed = true;

      PDM_gen_gnum_t *pdmGNum_handle_index = PDM_gnum_create(3, _npart, PDM_FALSE, 1e-3, _pdm_localComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

      for(int i_part=0;i_part<_npart;i_part++) {
        PDM_gnum_set_from_coords (pdmGNum_handle_index, i_part, _nVertex[i_part], _coords[i_part], NULL);
      }

      PDM_gnum_compute (pdmGNum_handle_index);

      for(int i_part=0;i_part<_npart;i_part++) {
        _global_num_vtx[i_part] =const_cast<CWP_g_num_t*>(PDM_gnum_get (pdmGNum_handle_index, i_part));
      }

      PDM_gnum_free (pdmGNum_handle_index);
    }

    for(int i_part=0;i_part<_npart;i_part++) {

      PDM_Mesh_nodal_coord_set(_pdmNodal_handle_index  ,
                               i_part                 ,
                               _nVertex       [i_part],
                               _coords        [i_part],
                               _global_num_vtx[i_part],
                               PDM_OWNERSHIP_USER);

      if(_visu->isCreated() && _displacement == CWP_DYNAMIC_MESH_STATIC) {
        _visu->GeomCoordSet(i_part,
        _nVertex       [i_part],
        _coords        [i_part],
        _global_num_vtx[i_part]);
      }

    }//loop i_part


    if (_faceEdgeMethod == 1) {

      int compute_gnum = 0;
      for (int i_part = 0; i_part < _npart; i_part++) {
        if (_faceLNToGN[i_part] != NULL) {
          compute_gnum = 1;
          break;
        }
      }

      if(compute_gnum) {
        _isEltGnumComputed = true;

        PDM_gen_gnum_t *pdmGNum_handle_index = PDM_gnum_create(3, _npart, PDM_FALSE, 1e-3, _pdm_localComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        double ** face_center = new double* [_npart];

        for (int i_part = 0; i_part < _npart; i_part++) {
          face_center[i_part] = new double[3*_nFace[i_part]];
          for (int j = 0; j < 3*_nFace[i_part]; j++) {
            face_center[i_part][j] = 0.;
          }

          for (int j = 0; j < _nFace[i_part]; j++) {
            int idx = _faceEdgeIdx[i_part][j];
            int nb = _faceEdgeNb[i_part][j];

            for (int k = idx; k < idx + nb; k++) {
              int iEdge = _faceEdge[i_part][k] - 1;
              for (int k1 = 0; k1 < 2; k1++) {
                int ivtx = _edgeVtx[i_part][2*iEdge + k1] - 1;
                for (int k2 = 0; k2 < 3; k2++) {
                  face_center[i_part][3*j+k2] += _coords[i_part][3*ivtx+k2];
                }
              }
            }

            for (int k2 = 0; k2 < 3; k2++) {
              face_center[i_part][3*j+k2] /= 2*nb;
            }
          }
        }

        for (int i_part = 0; i_part < _npart; i_part++) {
          PDM_gnum_set_from_coords (pdmGNum_handle_index, i_part, _nFace[i_part], face_center[i_part], NULL);
        }

        PDM_gnum_compute (pdmGNum_handle_index);

        for(int i_part=0;i_part<_npart;i_part++) {
          _faceLNToGN[i_part] = const_cast<CWP_g_num_t*>(PDM_gnum_get (pdmGNum_handle_index, i_part));
        }

        PDM_gnum_free (pdmGNum_handle_index);

        for (int i_part = 0; i_part < _npart; i_part++) {
          delete[] face_center[i_part];
        }

        delete[] face_center;

      }

      for (int i_part=0; i_part < _npart; i_part++) {

        PDM_Mesh_nodal_cell2d_celledge_add (_pdmNodal_handle_index,
                                            i_part,
                                            _nFace[i_part],
                                            _nEdge[i_part]    ,
                                            _edgeVtxIdx[i_part] ,
                                            _edgeVtxNb[i_part]  , //Number of vertices for each edge
                                            _edgeVtx[i_part]    ,
                                            _faceEdgeIdx[i_part],
                                            _faceEdgeNb[i_part] , //Number of edges for each faces
                                            _faceEdge[i_part]   ,
                                            _faceLNToGN[i_part],
                                            PDM_OWNERSHIP_USER);

      }//end i_part loop

      updateBlockDB();

    }

    else if(_cellFaceMethod == 1){

      int compute_gnum = 0;
      for (int i_part = 0; i_part < _npart; i_part++) {
        if (_faceLNToGN[i_part] != NULL) {
          compute_gnum = 1;
          break;
        }
      }

      if(compute_gnum) {

        _isEltGnumComputed = true;

        PDM_gen_gnum_t *pdmGNum_handle_index = PDM_gnum_create(3, _npart, PDM_FALSE, 1e-3, _pdm_localComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        double ** cell_center = new double* [_npart];

        for (int i_part = 0; i_part < _npart; i_part++) {
          cell_center[i_part] = new double[3*_nCells[i_part]];
          for (int j = 0; j < 3*_nCells[i_part]; j++) {
            cell_center[i_part][j] = 0.;
          }

          for (int j = 0; j < _nCells[i_part]; j++) {

            int weights = 0;
            int idx_face = _cellFaceIdx[i_part][j];
            int nb_face = _cellFaceNb[i_part][j];

            for (int j1 = idx_face; j1 < idx_face + nb_face; j1++) {

              int i_face = _cellFace[i_part][j1] - 1;

              int idx = _faceVtxIdx[i_part][i_face];
              int nb = _faceVtxNb[i_part][i_face];

              for (int k = idx; k < idx + nb; k++) {

                int i_vtx = _faceVtx[i_part][k] - 1;

                for (int k2 = 0; k2 < 3; k2++) {
                  cell_center[i_part][3*j+k2] += _coords[i_part][3*i_vtx+k2];
                }
                weights++;
              }
            }

            for (int k2 = 0; k2 < 3; k2++) {
              cell_center[i_part][3*j+k2] /= weights;
            }
          }
        }

        for (int i_part = 0; i_part < _npart; i_part++) {
          PDM_gnum_set_from_coords (pdmGNum_handle_index, i_part, _nCells[i_part], cell_center[i_part], NULL);
        }

        PDM_gnum_compute (pdmGNum_handle_index);

        for(int i_part=0;i_part<_npart;i_part++) {
          _cellLNToGN[i_part] = const_cast<CWP_g_num_t*>(PDM_gnum_get (pdmGNum_handle_index, i_part));
        }

        PDM_gnum_free (pdmGNum_handle_index);

        for (int i_part = 0; i_part < _npart; i_part++) {
          delete[] cell_center[i_part];
        }

        delete[] cell_center;

      }

      for(int i_part=0;i_part<_npart;i_part++){
        PDM_Mesh_nodal_cell3d_cellface_add(_pdmNodal_handle_index,
                                           i_part,
                                           _nCells[i_part],
                                           _nFace[i_part]    ,
                                           _faceVtxIdx[i_part],
                                           _faceVtxNb[i_part],
                                           _faceVtx[i_part],
                                           _cellFaceIdx[i_part],
                                           _cellFaceNb[i_part],
                                           _cellFace[i_part],
                                           _cellLNToGN[i_part],
                                           PDM_OWNERSHIP_USER);


      }//end i_part loop

      updateBlockDB();

    }

    else {

      for(int i_block=0;i_block<_nBlocks;i_block++){
        CWP_Block_t block_type = _blockDB[i_block]->blockTypeGet();
        PDM_Mesh_nodal_elt_t pdm_block_type;

        if (block_type == CWP_BLOCK_FACE_POLY) {
          pdm_block_type = PDM_MESH_NODAL_POLY_2D;
        }
        else if (block_type == CWP_BLOCK_CELL_POLY) {
          pdm_block_type = PDM_MESH_NODAL_POLY_3D;

        }
        else if (block_type == CWP_BLOCK_CELL_TETRA4) {
          pdm_block_type = PDM_MESH_NODAL_TETRA4;

        }
        else if (block_type == CWP_BLOCK_CELL_HEXA8) {
          pdm_block_type = PDM_MESH_NODAL_HEXA8;

        }
        else if (block_type == CWP_BLOCK_CELL_PRISM6) {
          pdm_block_type = PDM_MESH_NODAL_PRISM6;

        }
        else if (block_type == CWP_BLOCK_CELL_PYRAM5) {
          pdm_block_type = PDM_MESH_NODAL_PYRAMID5;

        }
        else if (block_type == CWP_BLOCK_FACE_QUAD4) {
          pdm_block_type = PDM_MESH_NODAL_QUAD4;

        }
        else if (block_type == CWP_BLOCK_FACE_TRIA3) {
          pdm_block_type = PDM_MESH_NODAL_TRIA3;

        }
        else if (block_type == CWP_BLOCK_EDGE2) {
          pdm_block_type = PDM_MESH_NODAL_BAR2;

        }
        else {
          PDM_error (__FILE__, __LINE__, 0, "unknown block type\n");
        }

        _blockDB[i_block]->blockIDPDMSet(PDM_Mesh_nodal_block_add (_pdmNodal_handle_index, pdm_block_type, PDM_OWNERSHIP_USER));

      } //end loop on block

      int compute_gnum = 0;

      for(int i_part = 0; i_part < _npart; i_part++){

        for(int i_block = 0; i_block < _nBlocks; i_block++){

          if (_blockDB[i_block]->GNumMeshGet(i_part) == NULL) {
            compute_gnum = 1;
            break;
          }
        }
      }

      if (compute_gnum) {

        _isEltGnumComputed = true;

        PDM_gen_gnum_t *pdmGNum_handle_index = PDM_gnum_create(3, _npart, PDM_FALSE, 1e-3, _pdm_localComm, PDM_OWNERSHIP_KEEP);

        double **cell_center = new double* [_npart];

        for(int i_part = 0; i_part < _npart; i_part++){
          cell_center[i_part] = new double[3*_nElts[i_part]];

          for(int i = 0; i < 3*_nElts[i_part]; i++){
            cell_center[i_part][i] = 0.;
          }

          int ielt = 0;
          for(int i_block = 0; i_block < _nBlocks; i_block++){

            int n_elt = _blockDB[i_block]->NEltsGet(i_part);

            CWP_Block_t block_type = _blockDB[i_block]->blockTypeGet();

            if (block_type == CWP_BLOCK_FACE_POLY) {
              BlockFP *block = dynamic_cast<BlockFP *>(_blockDB[i_block]);

              for (int j = 0; j < n_elt; j++) {
                int idx = block->ConnecIDXGet()[i_part][j];
                int nb  = block->ConnecIDXGet()[i_part][j+1] - idx;

                for (int k = idx; k < idx + nb; k++) {
                  int ivtx = block->ConnecGet()[i_part][k] - 1;
                  for (int k2 = 0; k2 < 3; k2++) {
                    cell_center[i_part][3*ielt+k2] += _coords[i_part][3*ivtx+k2];
                  }
                }

                for (int k2 = 0; k2 < 3; k2++) {
                  cell_center[i_part][3*ielt+k2] /= nb;
                }

                ielt++;
              }

            }

            else if (block_type == CWP_BLOCK_CELL_POLY) {
              BlockCP *block = dynamic_cast<BlockCP *>(_blockDB[i_block]);

              for (int j = 0; j < n_elt; j++) {

                int weights = 0;
                int idx_face = block->ConnecIDXGet()[i_part][j];
                int nb_face = block->ConnecIDXGet()[i_part][j+1] - idx_face;

                for (int j1 = idx_face; j1 < idx_face + nb_face; j1++) {

                  int i_face = std::abs(block->ConnecGet()[i_part][j1]) - 1;

                  int idx = block->ConnecFacesIDXGet()[i_part][i_face];
                  int nb  = block->ConnecFacesIDXGet()[i_part][i_face+1] - idx;

                  for (int k = idx; k < idx + nb; k++) {

                    int i_vtx = block->ConnecFacesGet()[i_part][k] - 1;

                    for (int k2 = 0; k2 < 3; k2++) {
                      cell_center[i_part][3*ielt+k2] += _coords[i_part][3*i_vtx+k2];
                    }
                    weights++;
                  }
                }

                for (int k2 = 0; k2 < 3; k2++) {
                  cell_center[i_part][3*ielt+k2] /= weights;
                }

                ielt++;
              }
            }

            else {
              BlockStd *block = dynamic_cast<BlockStd *>(_blockDB[i_block]);
              int n_vtx_elt = 0;

              if (block_type == CWP_BLOCK_CELL_TETRA4) {
                n_vtx_elt = 4;
              }
              else if (block_type == CWP_BLOCK_CELL_HEXA8) {
                n_vtx_elt = 8;
              }
              else if (block_type == CWP_BLOCK_CELL_PRISM6) {
                n_vtx_elt = 6;
              }
              else if (block_type == CWP_BLOCK_CELL_PYRAM5) {
                n_vtx_elt = 5;
              }
              else if (block_type == CWP_BLOCK_FACE_QUAD4) {
                n_vtx_elt = 4;
              }
              else if (block_type == CWP_BLOCK_FACE_TRIA3) {
                n_vtx_elt = 3;
              }
              else if (block_type == CWP_BLOCK_EDGE2) {
                n_vtx_elt = 2;
              }
              else {
                PDM_error (__FILE__, __LINE__, 0, "unknown block type\n");
              }

              for (int j = 0; j < n_elt; j++) {

                for (int k = 0; k < n_vtx_elt; k++) {

                  int i_vtx = block->ConnecGet()[i_part][j*n_vtx_elt + k] - 1;

                  for (int k2 = 0; k2 < 3; k2++) {
                    cell_center[i_part][3*ielt+k2] += _coords[i_part][3*i_vtx+k2];
                  }  
                }

                for (int k2 = 0; k2 < 3; k2++) {
                  cell_center[i_part][3*ielt+k2] /= n_vtx_elt;
                }
                ielt++;
              }

              block->ConnecGet()[i_part];

            }
          }

          PDM_gnum_set_from_coords (pdmGNum_handle_index, i_part, _nElts[i_part], cell_center[i_part], NULL);

        }

        PDM_gnum_compute (pdmGNum_handle_index);

        for (int i_part = 0; i_part < _npart; i_part++) {

          PDM_g_num_t *gnum = const_cast<PDM_g_num_t*>(PDM_gnum_get (pdmGNum_handle_index, i_part));

          for(int i_block = 0; i_block < _nBlocks; i_block++){

            int n_elt = _blockDB[i_block]->NEltsGet(i_part);
            CWP_g_num_t *block_gnum = (CWP_g_num_t *) malloc (sizeof(CWP_g_num_t) * n_elt);

            memcpy(block_gnum, gnum, n_elt * sizeof(CWP_g_num_t));

            _blockDB[i_block]->GNumMeshSet(i_part, block_gnum, PDM_OWNERSHIP_KEEP);

          }


        }

        PDM_gnum_free (pdmGNum_handle_index);

        for (int i_part = 0; i_part < _npart; i_part++) {
          delete[] cell_center[i_part];
        }

        delete[] cell_center;

      }

      for(int i_part = 0; i_part < _npart; i_part++){

        for(int i_block = 0; i_block < _nBlocks; i_block++){
          int n_elt = _blockDB[i_block]->NEltsGet(i_part);

          CWP_Block_t block_type = _blockDB[i_block]->blockTypeGet();

          int pdm_id_block = _blockDB[i_block]->blockIDPDMGet();

          if (block_type == CWP_BLOCK_FACE_POLY) {
            BlockFP *block = dynamic_cast<BlockFP *>(_blockDB[i_block]);

            PDM_Mesh_nodal_block_poly2d_set (_pdmNodal_handle_index,
                                             pdm_id_block,
                                             i_part,
                                             n_elt,
                                             block->ConnecIDXGet()[i_part],
                                             block->ConnecGet()[i_part],
                                             _blockDB[i_block]->GNumMeshGet(i_part),
                                             NULL);

          }

          else if (block_type == CWP_BLOCK_CELL_POLY) {
            BlockCP *block = dynamic_cast<BlockCP *>(_blockDB[i_block]);

            PDM_Mesh_nodal_block_poly3d_set (_pdmNodal_handle_index,
                                             pdm_id_block,
                                             i_part,
                                             n_elt,
                                             block->NFacesGet()[i_part],
                                             block->ConnecFacesIDXGet()[i_part],
                                             block->ConnecFacesGet()[i_part],
                                             block->ConnecIDXGet()[i_part],
                                             block->ConnecGet()[i_part],
                                             block->GNumMeshGet(i_part),
                                             NULL);

          }

          else {
            BlockStd *block = dynamic_cast<BlockStd *>(_blockDB[i_block]);

            PDM_Mesh_nodal_block_std_set (_pdmNodal_handle_index,
                                          pdm_id_block,
                                          i_part,
                                          n_elt,
                                          block->ConnecGet()[i_part],
                                          block->GNumMeshGet(i_part),
                                          NULL);
          }
        }
      }   //end loop on block
    }

    for (int i_block = 0 ; i_block < _nBlocks ; i_block++) {
      _blockDB[i_block]->geomFinalize();
    } //Loop on blockDB


    _nBlocks     = PDM_Mesh_nodal_n_blocks_get (_pdmNodal_handle_index);
    _blocks_id   = PDM_Mesh_nodal_blocks_id_get(_pdmNodal_handle_index);

    printf("n elts : %d\n ", PDM_Mesh_nodal_n_cell_get (_pdmNodal_handle_index, 0));

    if(_visu->isCreated() && _displacement == CWP_DYNAMIC_MESH_STATIC ) {
      _visu->GeomWrite(this);
    }

  }


  void
  Mesh::stdBlockSet
  (
    const int              i_part,
    const int              block_id,
    const int              n_elts,
    int                    connec[],
    CWP_g_num_t            global_num[]
  )
  {

    BlockStd *block = dynamic_cast <BlockStd *> (_blockDB [block_id]);
    block->blockSet(i_part,n_elts,connec,global_num);

    _nElts[i_part]  += n_elts;

  }

  /*************************************************/

  void
  Mesh::poly2DBlockSet
  (
    const int              i_part,
    const int              block_id,
    const int              n_elts,
    int                    connec_idx[],
    int                    connec[],
    CWP_g_num_t            global_num[]
  )
  {

    BlockFP *block = dynamic_cast <BlockFP *> (_blockDB [block_id]);

    block->blockSet(i_part,
                    n_elts,
                    connec_idx,
                    connec,
                    global_num);
    _nElts[i_part]  += n_elts;

  }

  /**********************************************************************/


  void
  Mesh::poly3DBlockSet
  (
    const int              i_part,
    const int              block_id,
    const int              n_elts,
    const int              n_faces,
    int                    connec_faces_idx[],
    int                    connec_faces[],
    int                    connec_cells_idx[],
    int                    connec_cells[],
    CWP_g_num_t            global_num[]
  )
  {

    BlockCP *block = dynamic_cast <BlockCP *> (_blockDB [block_id]);
    block->blockSet(i_part,n_elts,
                    n_faces,
                    connec_faces_idx,
                    connec_faces,
                    connec_cells_idx,
                    connec_cells,
                    global_num);

    if (_visu->isCreated() && _displacement == CWP_DYNAMIC_MESH_STATIC) {
      _visu->GeomBlockPoly3D (_id_visu[block_id],
                                i_part,
                                n_elts,
                                n_faces,
                                connec_faces_idx,
                                connec_faces,
                                connec_cells_idx,
                                connec_cells,
                                global_num);
    }

    _nElts[i_part]  +=  n_elts;

  }


  void
  Mesh::meshDel()
  {
    if (_pdmNodal_handle_index != NULL) {
      PDM_Mesh_nodal_free(_pdmNodal_handle_index);
    }
  }


  int
  Mesh::blockAdd
  (
    const CWP_Block_t      block_type
  )
  {
    Block *myBlock = FB::getInstance().CreateObject(block_type);
    myBlock->BlockAdd(block_type,this);
    int block_id_cwipi = _nBlocks;
    _blockDB.push_back (myBlock);
    myBlock->blockIDCWIPISet(block_id_cwipi);

    _nBlocks   = _blockDB.size();

    if(_visu->isCreated() && _displacement == CWP_DYNAMIC_MESH_STATIC) {
      int id_visu = _visu->GeomBlockAdd(block_type);
      _id_visu.insert(std::pair <int,int> (myBlock->blockIDCWIPIGet(),id_visu));
    }

    return myBlock->blockIDCWIPIGet();

  }

  void
  Mesh::fromCellFaceSet
  (
    const int   i_part,
    const int   n_cells,
    int         cell_face_idx[],
    int         cell_face[],
    int         n_faces,
    int         face_vtx_idx[],
    int         face_vtx[],
    CWP_g_num_t global_num[]
  )
  {
    _cellFaceMethod = 1;
    _cellLNToGN[i_part] = global_num;

    _faceVtxIdx[i_part] = face_vtx_idx;
    _cellFaceIdx[i_part] = cell_face_idx;
    _faceVtx[i_part] = face_vtx;
    _cellFace[i_part] = cell_face;

    _faceVtxNb[i_part] = (int *) malloc(sizeof(int) * n_faces);
    for (int i = 0; i < n_faces; i++) {
      _faceVtxNb[i_part][i] = face_vtx_idx[i + 1] - face_vtx_idx[i];
    }

    _cellFaceNb[i_part] = (int *) malloc(sizeof(int) * n_cells);
    for (int i = 0; i < n_cells; i++) {
      _cellFaceNb[i_part][i] = cell_face_idx[i + 1] - cell_face_idx[i];
    }

    _nFace[i_part] = n_faces;
    _nCells[i_part] = n_cells;

  }


  void
  Mesh::fromFacesEdgeSet
  (
    const int   i_part,
    const int   n_faces,
    int         face_edge_idx[],
    int         face_edge[],
    const int   n_edges,
    int         edge_vtx_idx[],
    int         edge_vtx[],
    CWP_g_num_t global_num[]
  )
  {
    _faceEdgeMethod = 1;

    _faceLNToGN[i_part] = global_num;
    _edgeVtxIdx[i_part] = edge_vtx_idx;
    _faceEdgeIdx[i_part] = face_edge_idx;
    _edgeVtx[i_part] = edge_vtx;
    _faceEdge[i_part] = face_edge;

    _edgeVtxNb[i_part] = (int *) malloc(sizeof(int) * n_edges);
    for (int i = 0 ; i < n_edges ; i++) {
      _edgeVtxNb[i_part][i] = edge_vtx_idx[i + 1] - edge_vtx_idx[i];
    }

    _faceEdgeNb[i_part] = (int *) malloc(sizeof(int) * n_faces);

    for (int i = 0; i < n_faces; i++) {
      _faceEdgeNb[i_part][i] = face_edge_idx[i + 1] - face_edge_idx[i];
    }

    _nEdge[i_part] = n_edges;
    _nFace[i_part] = n_faces;
  }


  int Mesh::getPartNElts(int id_part) const
  {
    CWP_Interface_t ed = _cpl->entitiesDimGet();

    if (ed == CWP_INTERFACE_SURFACE) {
      return PDM_Mesh_nodal_n_cell_get(_pdmNodal_handle_index, id_part);
    }
    else if (ed == CWP_INTERFACE_VOLUME) {
      return PDM_Mesh_nodal_n_cell_get(_pdmNodal_handle_index, id_part);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "getPartNElts : Element type is no taking account.\n");
    }

    return -1;
  }

}

/**
 * \endcond
 */
