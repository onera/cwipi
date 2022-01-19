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
#include <cmath>
#include "pdm_gnum_location.h"

#include "spatialInterpLocation.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"

/**
 * \cond
 */

namespace cwipi {
  SpatialInterpLocation::SpatialInterpLocation() = default;

  SpatialInterpLocation::~SpatialInterpLocation
  (
  )
  {

    for (int i = 0; i < _nPart; i++) {
      if (_tgt_distance[i] != NULL) {
        free (_tgt_distance[i]);
      }
      if (_tgt_projected[i] != NULL) {
        free (_tgt_projected[i]);
      }
      if (_tgt_closest_elt_gnum[i] != NULL) {
        free (_tgt_closest_elt_gnum[i]);
      }

      if (_elt_pts_inside_idx[i] != NULL) {
        free (_elt_pts_inside_idx[i]);
        free (_points_gnum[i]);
        free (_points_coords[i]);
        free (_points_dist2[i]);
        free (_points_projected_coords[i]);
      }

      if (_points_uvw[i] != NULL) {
        free (_points_uvw[i]);
      }      
    }
    delete[] _tgt_distance;
    delete[] _tgt_projected;
    delete[] _tgt_closest_elt_gnum;

    delete[] _elt_pts_inside_idx;
    delete[] _points_gnum;
    delete[] _points_coords;
    delete[] _points_uvw;
    delete[] _points_dist2;
    delete[] _points_projected_coords;
  }


  /**
    *
    * \brief SpatialInterp location Init.
    *
    */

  void 
  SpatialInterpLocation::init 
  (
    Coupling           *coupling, 
    CWP_Dof_location_t localCodeDofLOcation,
    CWP_Dof_location_t cplCodeDofLOcation,
    SpatialInterpExchDirection exchDirection 
  )
  {
    SpatialInterp::init (coupling, 
                         localCodeDofLOcation, 
                         cplCodeDofLOcation, 
                         exchDirection);

    _interpolation_time = CWP_SPATIAL_INTERP_AT_SEND;


    // Create nodal mesh
    _pdm_CplNodal = PDM_Mesh_nodal_create (_nPart, _pdmCplComm);

    if (!_coupledCodeProperties->localCodeIs()) {
      if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
        for (int i_part = 0 ; i_part < _nPart ; i_part++) {
          PDM_Mesh_nodal_coord_set(_pdm_CplNodal,
                                   i_part,
                                   _mesh->getPartNVertex(i_part),
                                   _mesh->getVertexCoords(i_part),
                                   _mesh->getVertexGNum(i_part));
        }

        for (int i_block = 0 ; i_block < _mesh->nBlockGet() ; i_block++) {
          int CWP_block_type = _mesh->blockTypeGet(i_block);
          PDM_Mesh_nodal_elt_t pdm_block_type;
          switch (CWP_block_type) {
            case CWP_BLOCK_NODE: pdm_block_type = PDM_MESH_NODAL_POINT;
              break;
            case CWP_BLOCK_EDGE2: pdm_block_type = PDM_MESH_NODAL_BAR2;
              break;
            case CWP_BLOCK_FACE_TRIA3: pdm_block_type = PDM_MESH_NODAL_TRIA3;
              break;
            case CWP_BLOCK_FACE_QUAD4: pdm_block_type = PDM_MESH_NODAL_QUAD4;
              break;
            case CWP_BLOCK_CELL_TETRA4: pdm_block_type = PDM_MESH_NODAL_TETRA4;
              break;
            case CWP_BLOCK_FACE_POLY: pdm_block_type = PDM_MESH_NODAL_POLY_2D;
              break;
            case CWP_BLOCK_CELL_HEXA8: pdm_block_type = PDM_MESH_NODAL_HEXA8;
              break;
            case CWP_BLOCK_CELL_PYRAM5: pdm_block_type = PDM_MESH_NODAL_PYRAMID5;
              break;
            case CWP_BLOCK_CELL_PRISM6: pdm_block_type = PDM_MESH_NODAL_PRISM6;
              break;
            case CWP_BLOCK_CELL_POLY: pdm_block_type = PDM_MESH_NODAL_POLY_3D;
              break;
            default:pdm_block_type = PDM_MESH_NODAL_POINT;
              PDM_error(__FILE__, __LINE__, 0, "No referenced CWP_Block_t.\n");
          }

          int _block_id_pdm = PDM_Mesh_nodal_block_add(_pdm_CplNodal,
                                                       PDM_FALSE,
                                                       pdm_block_type);

          for (int i_part = 0 ; i_part < _nPart ; i_part++) {
            if (CWP_block_type == CWP_BLOCK_FACE_POLY) {
              PDM_Mesh_nodal_block_poly2d_set(_pdm_CplNodal,
                                              _block_id_pdm,
                                              i_part,
                                              _mesh->getPartNElts(i_part),
                                              _mesh->getEltConnectivityIndex(i_block, i_part),
                                              _mesh->getEltConnectivity(i_block, i_part),
                                              _mesh->GNumEltsGet(i_part),
                                              NULL);
            }
            else if (CWP_block_type == CWP_BLOCK_CELL_POLY) {
//            PDM_Mesh_nodal_block_poly3d_set (_pdm_CplNodal,
//                                             _block_id_pdm,
//                                             i_part,
//                                             _mesh->getPartNElts(i_part),
//                                             _mesh->getPartNFaces(i_part),
//                                             _mesh->getEltConnectivityIndex(i_block, i_part),
//                                             _mesh->getEltConnectivity(i_block, i_part),
//                                             _connec_cells_idx[i_part],
//                                             _connec_cells[i_part],
//                                             _global_num[i_part],
//                                             NULL);
            }
            else {
              PDM_Mesh_nodal_block_std_set(_pdm_CplNodal,
                                           _block_id_pdm,
                                           i_part,
                                           _mesh->getPartNElts(i_part),
                                           _mesh->getEltConnectivity(i_block, i_part),
                                           _mesh->GNumEltsGet(i_part),
                                           NULL);
            }
          }
        }
      }

        // Receive
      else {
        for (int i_part = 0 ; i_part < _nPart ; i_part++) {
          PDM_Mesh_nodal_coord_set(_pdm_CplNodal,
                                   i_part,
                                   0,
                                   NULL,
                                   NULL);

          for (int i_block = 0 ; i_block < _mesh->nBlockGet() ; i_block++) {
            int CWP_block_type = _mesh->blockTypeGet(i_block);
            PDM_Mesh_nodal_elt_t pdm_block_type;
            switch (CWP_block_type) {
              case CWP_BLOCK_NODE: pdm_block_type = PDM_MESH_NODAL_POINT;
                break;
              case CWP_BLOCK_EDGE2: pdm_block_type = PDM_MESH_NODAL_BAR2;
                break;
              case CWP_BLOCK_FACE_TRIA3: pdm_block_type = PDM_MESH_NODAL_TRIA3;
                break;
              case CWP_BLOCK_FACE_QUAD4: pdm_block_type = PDM_MESH_NODAL_QUAD4;
                break;
              case CWP_BLOCK_CELL_TETRA4: pdm_block_type = PDM_MESH_NODAL_TETRA4;
                break;
              case CWP_BLOCK_FACE_POLY: pdm_block_type = PDM_MESH_NODAL_POLY_2D;
                break;
              case CWP_BLOCK_CELL_HEXA8: pdm_block_type = PDM_MESH_NODAL_HEXA8;
                break;
              case CWP_BLOCK_CELL_PYRAM5: pdm_block_type = PDM_MESH_NODAL_PYRAMID5;
                break;
              case CWP_BLOCK_CELL_PRISM6: pdm_block_type = PDM_MESH_NODAL_PRISM6;
                break;
              case CWP_BLOCK_CELL_POLY: pdm_block_type = PDM_MESH_NODAL_POLY_3D;
                break;
              default:pdm_block_type = PDM_MESH_NODAL_POINT;
                PDM_error(__FILE__, __LINE__, 0, "No referenced CWP_Block_t.\n");
            }

            int _block_id_pdm = PDM_Mesh_nodal_block_add(_pdm_CplNodal,
                                                         PDM_FALSE,
                                                         pdm_block_type);

            for (int i_part = 0 ; i_part < _nPart ; i_part++) {
              if (CWP_block_type == CWP_BLOCK_FACE_POLY) {
                PDM_Mesh_nodal_block_poly2d_set(_pdm_CplNodal,
                                                _block_id_pdm,
                                                i_part,
                                                0,
                                                fake_idx,
                                                NULL,
                                                NULL,
                                                NULL);
              }
              else if (CWP_block_type == CWP_BLOCK_CELL_POLY) {
//            PDM_Mesh_nodal_block_poly3d_set (_pdm_CplNodal,
//                                             _block_id_pdm,
//                                             i_part,
//                                             0,
//                                             0,
//                                             fake_idx,
//                                             NULL,
//                                             fake_idx,
//                                             NULL,
//                                             NULL,
//                                             NULL);
              }
              else {
                PDM_Mesh_nodal_block_std_set(_pdm_CplNodal,
                                             _block_id_pdm,
                                             i_part,
                                             0,
                                             NULL,
                                             NULL,
                                             NULL);
              }
            }
          }
        }
      }
    }
    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
          // Comme au-dessus

          for (int i_part = 0 ; i_part < _nPart ; i_part++) {
            PDM_Mesh_nodal_coord_set(_pdm_CplNodal,
                                     i_part,
                                     _mesh->getPartNVertex(i_part),
                                     _mesh->getVertexCoords(i_part),
                                     _mesh->getVertexGNum(i_part));
          }

          for (int i_block = 0 ; i_block < _mesh->nBlockGet() ; i_block++) {
            int CWP_block_type = _mesh->blockTypeGet(i_block);
            PDM_Mesh_nodal_elt_t pdm_block_type;
            switch (CWP_block_type) {
              case CWP_BLOCK_NODE: pdm_block_type = PDM_MESH_NODAL_POINT;
                break;
              case CWP_BLOCK_EDGE2: pdm_block_type = PDM_MESH_NODAL_BAR2;
                break;
              case CWP_BLOCK_FACE_TRIA3: pdm_block_type = PDM_MESH_NODAL_TRIA3;
                break;
              case CWP_BLOCK_FACE_QUAD4: pdm_block_type = PDM_MESH_NODAL_QUAD4;
                break;
              case CWP_BLOCK_CELL_TETRA4: pdm_block_type = PDM_MESH_NODAL_TETRA4;
                break;
              case CWP_BLOCK_FACE_POLY: pdm_block_type = PDM_MESH_NODAL_POLY_2D;
                break;
              case CWP_BLOCK_CELL_HEXA8: pdm_block_type = PDM_MESH_NODAL_HEXA8;
                break;
              case CWP_BLOCK_CELL_PYRAM5: pdm_block_type = PDM_MESH_NODAL_PYRAMID5;
                break;
              case CWP_BLOCK_CELL_PRISM6: pdm_block_type = PDM_MESH_NODAL_PRISM6;
                break;
              case CWP_BLOCK_CELL_POLY: pdm_block_type = PDM_MESH_NODAL_POLY_3D;
                break;
              default:pdm_block_type = PDM_MESH_NODAL_POINT;
                PDM_error(__FILE__, __LINE__, 0, "No referenced CWP_Block_t.\n");
            }

            int _block_id_pdm = PDM_Mesh_nodal_block_add(_pdm_CplNodal,
                                                         PDM_FALSE,
                                                         pdm_block_type);

            for (int i_part = 0 ; i_part < _nPart ; i_part++) {
              if (CWP_block_type == CWP_BLOCK_FACE_POLY) {
                PDM_Mesh_nodal_block_poly2d_set(_pdm_CplNodal,
                                                _block_id_pdm,
                                                i_part,
                                                _mesh->getPartNElts(i_part),
                                                _mesh->getEltConnectivityIndex(i_block, i_part),
                                                _mesh->getEltConnectivity(i_block, i_part),
                                                _mesh->GNumEltsGet(i_part),
                                                NULL);
              }
              else if (CWP_block_type == CWP_BLOCK_CELL_POLY) {
//            PDM_Mesh_nodal_block_poly3d_set (_pdm_CplNodal,
//                                             _block_id_pdm,
//                                             i_part,
//                                             _mesh->getPartNElts(i_part),
//                                             _mesh->getPartNFaces(i_part),
//                                             _mesh->getEltConnectivityIndex(i_block, i_part),
//                                             _mesh->getEltConnectivity(i_block, i_part),
//                                             _connec_cells_idx[i_part],
//                                             _connec_cells[i_part],
//                                             _global_num[i_part],
//                                             NULL);
              }
              else {
                PDM_Mesh_nodal_block_std_set(_pdm_CplNodal,
                                             _block_id_pdm,
                                             i_part,
                                             _mesh->getPartNElts(i_part),
                                             _mesh->getEltConnectivity(i_block, i_part),
                                             _mesh->GNumEltsGet(i_part),
                                             NULL);
              }
            }
          }
        }
        else { // Recv
          // Remplir avec les données du code couplé

          cwipi::Coupling &cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());
          cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();

          for (int i_part = 0 ; i_part < _nPart ; i_part++) {
            PDM_Mesh_nodal_coord_set(_pdm_CplNodal,
                                     i_part,
                                     cpl_mesh->getPartNVertex(i_part),
                                     cpl_mesh->getVertexCoords(i_part),
                                     cpl_mesh->getVertexGNum(i_part));
          }

          for (int i_block = 0 ; i_block < cpl_mesh->nBlockGet() ; i_block++) {
            int CWP_block_type = cpl_mesh->blockTypeGet(i_block);
            PDM_Mesh_nodal_elt_t pdm_block_type;
            switch (CWP_block_type) {
              case CWP_BLOCK_NODE: pdm_block_type = PDM_MESH_NODAL_POINT;
                break;
              case CWP_BLOCK_EDGE2: pdm_block_type = PDM_MESH_NODAL_BAR2;
                break;
              case CWP_BLOCK_FACE_TRIA3: pdm_block_type = PDM_MESH_NODAL_TRIA3;
                break;
              case CWP_BLOCK_FACE_QUAD4: pdm_block_type = PDM_MESH_NODAL_QUAD4;
                break;
              case CWP_BLOCK_CELL_TETRA4: pdm_block_type = PDM_MESH_NODAL_TETRA4;
                break;
              case CWP_BLOCK_FACE_POLY: pdm_block_type = PDM_MESH_NODAL_POLY_2D;
                break;
              case CWP_BLOCK_CELL_HEXA8: pdm_block_type = PDM_MESH_NODAL_HEXA8;
                break;
              case CWP_BLOCK_CELL_PYRAM5: pdm_block_type = PDM_MESH_NODAL_PYRAMID5;
                break;
              case CWP_BLOCK_CELL_PRISM6: pdm_block_type = PDM_MESH_NODAL_PRISM6;
                break;
              case CWP_BLOCK_CELL_POLY: pdm_block_type = PDM_MESH_NODAL_POLY_3D;
                break;
              default:pdm_block_type = PDM_MESH_NODAL_POINT;
                PDM_error(__FILE__, __LINE__, 0, "No referenced CWP_Block_t.\n");
            }

            int _block_id_pdm = PDM_Mesh_nodal_block_add(_pdm_CplNodal,
                                                         PDM_FALSE,
                                                         pdm_block_type);

            for (int i_part = 0 ; i_part < _nPart ; i_part++) {
              if (CWP_block_type == CWP_BLOCK_FACE_POLY) {
                PDM_Mesh_nodal_block_poly2d_set(_pdm_CplNodal,
                                                _block_id_pdm,
                                                i_part,
                                                cpl_mesh->getPartNElts(i_part),
                                                cpl_mesh->getEltConnectivityIndex(i_block, i_part),
                                                cpl_mesh->getEltConnectivity(i_block, i_part),
                                                cpl_mesh->GNumEltsGet(i_part),
                                                NULL);
              }
              else if (CWP_block_type == CWP_BLOCK_CELL_POLY) {
//            PDM_Mesh_nodal_block_poly3d_set (_pdm_CplNodal,
//                                             _block_id_pdm,
//                                             i_part,
//                                             cpl_mesh->getPartNElts(i_part),
//                                             cpl_mesh->getPartNFaces(i_part),
//                                             cpl_mesh->getEltConnectivityIndex(i_block, i_part),
//                                             cpl_mesh->getEltConnectivity(i_block, i_part),
//                                             _connec_cells_idx[i_part],
//                                             _connec_cells[i_part],
//                                             _global_num[i_part],
//                                             NULL);
              }
              else {
                PDM_Mesh_nodal_block_std_set(_pdm_CplNodal,
                                             _block_id_pdm,
                                             i_part,
                                             cpl_mesh->getPartNElts(i_part),
                                             cpl_mesh->getEltConnectivity(i_block, i_part),
                                             cpl_mesh->GNumEltsGet(i_part),
                                             NULL);
              }
            } // i_part
          } // i_block
        } // recv
      } // _localCodeProperties->idGet() < _coupledCodeProperties->idGet()
    }
    //
    // Data for PDM_part1_to_selected_part2_t

    if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
      for (int i_part = 0 ; i_part < _nPart ; i_part++) { 
       _src_gnum[i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
       _src_n_gnum[i_part] = _mesh->getPartNElts (i_part);
      }
    }

    else {
      for (int i_part = 0 ; i_part < _nPart ; i_part++) { 
        if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
          _tgt_gnum[i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
          _tgt_n_gnum[i_part] = _mesh->getPartNElts (i_part);

        }
        else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
          _tgt_gnum[i_part] = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
          _tgt_n_gnum[i_part] = _mesh->getPartNVertex (i_part);            
        }
        else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
          _tgt_gnum[i_part] = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
          _tgt_n_gnum[i_part] = _cpl->userTargetNGet (i_part);
        }
      }
    }

    //
    // Target properties
    
    _tgt_distance = new double* [_nPart];                 // Distance to the closest source element surface by partition
    _tgt_projected = new double* [_nPart];                // Projected point coordinates (on the closest source element surface)
    _tgt_closest_elt_gnum = new CWP_g_num_t* [_nPart];    // Closest source element global numbering

    //
    // Source properties

    _elt_pts_inside_idx = new int* [_nPart];
    _points_gnum = new CWP_g_num_t* [_nPart];
    _points_coords = new double* [_nPart];
    _points_uvw = new double* [_nPart];
    _points_dist2 = new double* [_nPart];
    _points_projected_coords = new double* [_nPart];

    for (int i_part = 0; i_part < _nPart; i_part++) {
      _tgt_distance[i_part] = NULL;
      _tgt_projected[i_part] = NULL;
      _tgt_closest_elt_gnum[i_part] = NULL;
      _elt_pts_inside_idx[i_part] = NULL;
      _points_gnum[i_part] = NULL;
      _points_coords[i_part] = NULL;
      _points_uvw[i_part] = NULL;
      _points_dist2[i_part] = NULL;
      _points_projected_coords[i_part] = NULL;
    }

  }


  void SpatialInterpLocation::weightsCompute() 
  {
    for (int i_part = 0; i_part < _nPart; i_part++) {
      if (_tgt_distance[i_part] != NULL) {
        free (_tgt_distance[i_part]);
        free (_tgt_projected[i_part]);
        free (_tgt_closest_elt_gnum[i_part]);
      }
      if (_elt_pts_inside_idx[i_part] != NULL) {
        free (_elt_pts_inside_idx[i_part]);
        free (_points_gnum[i_part]);
        free (_points_coords[i_part]);
        free (_points_uvw[i_part]);
        free (_points_dist2[i_part]);
        free (_points_projected_coords[i_part]);
      }

      if (_weights_idx[i_part] != NULL) {
        free (_weights_idx[i_part]);
        free (_weights[i_part]);
      }

      if (_computed_tgt[i_part] != NULL) {
        free (_computed_tgt[i_part]);   
      }

      if (_uncomputed_tgt[i_part] != NULL) {
        free (_uncomputed_tgt[i_part]);   
      }

      _n_elt_weights[i_part] = 0;
      _weights_idx[i_part] = NULL;
      _weights[i_part] = NULL;

      _n_computed_tgt[i_part] = 0;
      _computed_tgt[i_part] = NULL;

      _n_uncomputed_tgt[i_part] = 0;
      _uncomputed_tgt[i_part] = NULL;

      _tgt_distance[i_part] = NULL;
      _tgt_projected[i_part] = NULL;
      _tgt_closest_elt_gnum[i_part] = NULL;
      _elt_pts_inside_idx[i_part] = NULL;
      _points_gnum[i_part] = NULL;
      _points_coords[i_part] = NULL;
      _points_uvw[i_part] = NULL;
      _points_dist2[i_part] = NULL;
      _points_projected_coords[i_part] = NULL;

    }

    localization_init();

    localization_points_cloud_setting();
    
    localization_surface_setting();

    localization_compute();

    localization_get();

    localization_free();

    if (_ptsp != nullptr) {
      if (!_coupledCodeProperties->localCodeIs()) {
        PDM_part1_to_selected_part2_free (_ptsp);
        _ptsp = nullptr;
      }
      else {

        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
          PDM_part1_to_selected_part2_free (_ptsp);
          _ptsp = nullptr;
  
          SpatialInterpLocation *cpl_spatial_interp;

          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet(); 
            cpl_spatial_interp = 
              dynamic_cast <SpatialInterpLocation *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }

          else {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 
            cpl_spatial_interp = 
              dynamic_cast <SpatialInterpLocation *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }

          cpl_spatial_interp->_ptsp = NULL;
        }

      }
    }


    if (!_coupledCodeProperties->localCodeIs()) {
//      printf("_src_gnum %d :", _src_n_gnum[0]);
//      for (int i = 0; i < _src_n_gnum[0]; i++) {
//        printf(" %ld", _src_gnum[0][i]);
//      }
//      printf("\n");

//      printf("_tgt_gnum %d :", _tgt_n_gnum[0]);
//      for (int i = 0; i < _tgt_n_gnum[0]; i++) {
//        printf(" %ld", _tgt_gnum[0][i]);
//      }
//      printf("\n");

//      for (int i = 0; i < _src_n_gnum[0]; i++) {
//        printf("_elt_pts_inside_idx %d %d: ", _elt_pts_inside_idx[0][i], _elt_pts_inside_idx[0][i+1]);
//        for (int j = _elt_pts_inside_idx[0][i]; j < _elt_pts_inside_idx[0][i+1]; j++) {
//          printf(" %ld", _points_gnum[0][j]);
//        }
//        printf("\n");
//      }
//      printf("\n");
//      fflush(stdout);
      _ptsp = PDM_part1_to_selected_part2_create ((const PDM_g_num_t **)_src_gnum,
                                                  (const int *)_src_n_gnum,
                                                  _nPart,
                                                  (const PDM_g_num_t **)_tgt_gnum,
                                                  (const int *)_tgt_n_gnum,
                                                  _nPart,
                                                  (const int **)_elt_pts_inside_idx,
                                                  (const PDM_g_num_t **)_points_gnum,
                                                  _pdmCplComm);                         
    }
    else {

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        SpatialInterpLocation *cpl_spatial_interp;

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet(); 
          cpl_spatial_interp = 
            dynamic_cast <SpatialInterpLocation *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }

        else {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 
          cpl_spatial_interp = 
            dynamic_cast <SpatialInterpLocation *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }

        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
          _ptsp = PDM_part1_to_selected_part2_create ((const PDM_g_num_t **)_src_gnum,
                                                      (const int *)_src_n_gnum,
                                                      _nPart,
                                                      (const PDM_g_num_t **)cpl_spatial_interp->_tgt_gnum,
                                                      (const int *)cpl_spatial_interp->_tgt_n_gnum,
                                                      _cplNPart,
                                                      (const int **)_elt_pts_inside_idx,
                                                      (const PDM_g_num_t **)_points_gnum,
                                                      _pdmCplComm);                         
        }
        else {
          _ptsp = PDM_part1_to_selected_part2_create ((const PDM_g_num_t **)cpl_spatial_interp->_src_gnum,
                                                      (const int *)cpl_spatial_interp->_src_n_gnum,
                                                      _cplNPart,
                                                      (const PDM_g_num_t **) _tgt_gnum,
                                                      (const int *) _tgt_n_gnum,
                                                      _nPart,
                                                      (const int **)cpl_spatial_interp->_elt_pts_inside_idx,
                                                      (const PDM_g_num_t **)cpl_spatial_interp->_points_gnum,
                                                      _pdmCplComm);                         

        }

        cpl_spatial_interp->_ptsp = _ptsp;
      } 
    }
  }

  void SpatialInterpLocation::interpolate (Field *referenceField, double **buffer) 
  {
    int nComponent = referenceField->nComponentGet();
    CWP_Dof_location_t referenceFieldType = referenceField->locationGet();
    int dataTypeSize = referenceField->dataTypeSizeGet();
    CWP_Interpolation_t interpolationType = referenceField->interpolationTypeGet();


    // int  _nPart; !< Mesh partition number                                                    


    if (interpolationType == CWP_INTERPOLATION_USER) {
    // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
          // int         *part_elt_pts_inside_idx      = _elt_pts_inside_idx[i_part];
          // PDM_g_num_t *part_points_gnum             = _points_gnum[i_part];
          // double      *part_points_coords           = _points_coords[i_part];
          // double      *part_points_uvw              = _points_uvw[i_part];
          // double      *part_points_dist2            = _points_dist2[i_part];
          // double      *part_points_projected_coords = _points_projected_coords[i_part];
          // double      *referenceData                = (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_SOURCE);

          // int          part_n_elt_weights           = _n_elt_weights[i_part];
          // int         *part_weights_idx             = _weights_idx[i_part];
          // double      *part_weights                 = _weights[i_part];

          // int          part_n_elt                   = _mesh->getPartNElts(i_part);
    //       (*interpolationFunction)(CWP_INTERFACE_SURFACE, _n_vtx[i_part], _n_elt[i_part], n_tgt, coords, connecIdx, connec,
    //                                tgt_pts_projected_coords, tgt_pts_location, tgt_pts_dist, tgt_pts_bary_coords_idx, tgt_pts_bary_coords,
    //                                nComponent, referenceFieldType, referenceData, referenceFieldType, tmpData);
    //
    // }

    }

    else { 
      if (referenceFieldType == CWP_DOF_LOCATION_CELL_CENTER) {

        for (int i_part = 0; i_part < _nPart; i_part++) {
          int         *part_elt_pts_inside_idx      = _elt_pts_inside_idx[i_part];
          void        *referenceData                = referenceField->dataGet(i_part, CWP_FIELD_MAP_SOURCE);

          int          part_n_elt                   = _mesh->getPartNElts(i_part);

          char *p_local_buffer = (char *) buffer;
          char *p_referenceData = (char *) referenceData;

          for (int i = 0; i < part_n_elt; i++) {
            for (int j = part_elt_pts_inside_idx[i]; j < part_elt_pts_inside_idx[i+1]; j++) {
              memcpy((char *) p_local_buffer, p_referenceData, nComponent * dataTypeSize);
              p_local_buffer += (part_elt_pts_inside_idx[j+1] - part_elt_pts_inside_idx[j]) * nComponent * dataTypeSize;
            } 
            p_referenceData += nComponent * dataTypeSize;        
          }

        }
      }
      
      else if (referenceFieldType == CWP_DOF_LOCATION_NODE || referenceFieldType == CWP_DOF_LOCATION_USER) {

        for (int i_part = 0; i_part < _nPart; i_part++) {
          int         *part_elt_pts_inside_idx      = _elt_pts_inside_idx[i_part];
          double      *referenceData                = (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_SOURCE);

          int         *part_weights_idx             = _weights_idx[i_part];
          double      *part_weights                 = _weights[i_part];

          int          part_n_elt                   = _mesh->getPartNElts(i_part);

          int         *connec_idx                   = _mesh->connecIdxGet(i_part);
          int         *connec                       = _mesh->connecGet(i_part);

          double *local_buffer = (double *) *buffer;

          int ival = 0;
          for (int i = 0; i < part_n_elt; i++) {
            for (int j = part_elt_pts_inside_idx[i]; j < part_elt_pts_inside_idx[i+1]; j++) {
//              if (connec_idx[i+1] - connec_idx[i] != part_weights_idx[j + 1] - part_weights_idx[j]) {
//                printf("\ti, j, connecs, weights %d %d: %d %d %d / %d %d %d\n", i, j, connec_idx[i],       connec_idx[i + 1],       connec_idx[i + 1]       - connec_idx[i],
//                                                                                      part_weights_idx[j], part_weights_idx[j + 1], part_weights_idx[j + 1] - part_weights_idx[j]);
//              }
              for (int k1 = 0; k1 < nComponent; k1++) {
                local_buffer[ival] = 0;
                int k2 = connec_idx[i];
                assert(connec_idx[i+1] - connec_idx[i] == part_weights_idx[j + 1] - part_weights_idx[j]);
                for (int k = part_weights_idx[j]; k < part_weights_idx[j+1]; k++) {
                  int isom = connec[k2++] - 1;
//                  printf("\t\ti, k, weights isom %d %d %f %d\n", i, k, part_weights[k], isom);
                  local_buffer[ival] += part_weights[k] * referenceData[isom*nComponent+k1];
                }
                ival++;
              }

//              printf("coords dist projected_x gnum weight_idx local_buffer %f %f %f %f %f %ld %d %f\n", _points_coords[i_part][3 * j], _points_coords[i_part][3 * j + 1], _points_coords[i_part][3 * j + 2], _points_dist2[i_part][j], _points_projected_coords[i_part][3 * j], _points_gnum[i_part][j], part_weights_idx[j], local_buffer[ival - 1]);
            }
          }
        }
      }
    }
  }

  /**********************************************************
  ***********************************************************
  **                                                       **
  **            Localization object functions              **
  **                                                       **
  ***********************************************************
  **********************************************************/

  void SpatialInterpLocation::localization_init()
  {
    PDM_error(__FILE__, __LINE__, 0, "Unknown location method.\n");
  }

  void SpatialInterpLocation::localization_points_cloud_setting() 
  {
    PDM_error(__FILE__, __LINE__, 0, "Unknown location method.\n");
  }

  void SpatialInterpLocation::localization_surface_setting() 
  {
    PDM_error(__FILE__, __LINE__, 0, "Unknown location method.\n");
  }

  void SpatialInterpLocation::localization_compute() 
  {
    PDM_error(__FILE__, __LINE__, 0, "Unknown location method.\n");
  }

  void SpatialInterpLocation::localization_get() 
  {
    PDM_error(__FILE__, __LINE__, 0, "Unknown location method.\n");
  }

  void SpatialInterpLocation::localization_free() 
  {
    PDM_error(__FILE__, __LINE__, 0, "Unknown location method.\n");
  }



}

/**
 * \endcond
 */
