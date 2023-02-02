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

#include <spatialInterpIntersection.hxx>
#include "pdm_mesh_intersection.h"
#include "cwp_priv.h"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include "pdm_array.h"
#include "pdm_logging.h"

namespace cwipi {
  SpatialInterpIntersection::SpatialInterpIntersection() = default;

  SpatialInterpIntersection::~SpatialInterpIntersection
  (
   )
  {
    for (int i_part = 0; i_part < _nPart; i_part++) {
      if (_src_to_tgt_idx[i_part] != NULL) {
        free(_src_to_tgt_idx[i_part]);
      }
      if (_src_to_tgt_gnum[i_part] != NULL) {
        free(_src_to_tgt_gnum[i_part]);
      }
      if (_src_to_tgt_weight[i_part] != NULL) {
        free(_src_to_tgt_weight[i_part]);
      }

      if (_tgt_to_src_weight[i_part] != NULL) {
        free(_tgt_to_src_weight[i_part]);
      }
    }

    delete[] _src_to_tgt_idx;
    delete[] _src_to_tgt_gnum;
    delete[] _src_to_tgt_weight;
    delete[] _tgt_to_src_weight;
  }


  /**
     *
     * \brief SpatialInterpIntersection Init.
     *
     */

    void
    SpatialInterpIntersection::init
    (
     Coupling                   *coupling,
     CWP_Dof_location_t          localCodeDofLocation,
     CWP_Dof_location_t          cplCodeDofLocation,
     SpatialInterpExchDirection  exchDirection
     )
    {
      SpatialInterp::init (coupling,
                           localCodeDofLocation,
                           cplCodeDofLocation,
                           exchDirection);

      _interpolation_time = CWP_SPATIAL_INTERP_AT_RECV;

      //
      // Map nodal mesh
      _pdm_CplNodal = _mesh->getPdmNodalIndex();

      //
      // Data for PDM_part_to_part_t
      if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
        for (int i_part = 0 ; i_part < _nPart ; i_part++) {
          if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            _src_gnum  [i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
            _src_n_gnum[i_part] = _mesh->getPartNElts (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
            PDM_error(__FILE__, __LINE__, 0, "node-centered version of PDM_mesh_intersection not implemented yet\n");
            _src_gnum  [i_part] = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
            _src_n_gnum[i_part] = _mesh->getPartNVertex (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            PDM_error(__FILE__, __LINE__, 0, "Invalid dof location %d\n", (int) _localCodeDofLocation);
            _src_gnum  [i_part] = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
            _src_n_gnum[i_part] = _cpl->userTargetNGet (i_part);
          }
        }
      }
      else {
        for (int i_part = 0 ; i_part < _nPart ; i_part++) {
          if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            _tgt_gnum  [i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
            _tgt_n_gnum[i_part] = _mesh->getPartNElts (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
            PDM_error(__FILE__, __LINE__, 0, "node-centered version of PDM_mesh_intersection not implemented yet\n");
            _tgt_gnum  [i_part] = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
            _tgt_n_gnum[i_part] = _mesh->getPartNVertex (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            PDM_error(__FILE__, __LINE__, 0, "Invalid dof location %d\n", (int) _localCodeDofLocation);
            _tgt_gnum  [i_part] = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
            _tgt_n_gnum[i_part] = _cpl->userTargetNGet (i_part);
          }
        }
      }


      //
      // Target properties
      _src_to_tgt_idx    = new int*         [_nPart];
      _src_to_tgt_gnum   = new PDM_g_num_t* [_nPart];
      _src_to_tgt_weight = new double*      [_nPart];
      _tgt_to_src_weight = new double*      [_nPart];

      for (int i_part = 0; i_part < _nPart; i_part++) {
        _src_to_tgt_idx   [i_part] = NULL;
        _src_to_tgt_gnum  [i_part] = NULL;
        _src_to_tgt_weight[i_part] = NULL;
        _tgt_to_src_weight[i_part] = NULL;
      }
    }

  void SpatialInterpIntersection::weightsCompute() {

    if (!_coupledCodeProperties->localCodeIs() ||
        (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet())) {

      for (int i_part = 0; i_part < _nPart; i_part++) {

        if (_src_to_tgt_idx[i_part] != NULL) {
          free(_src_to_tgt_idx[i_part]);
        }
        if (_src_to_tgt_gnum[i_part] != NULL) {
          free(_src_to_tgt_gnum[i_part]);
        }
        if (_src_to_tgt_weight[i_part] != NULL) {
          free(_src_to_tgt_weight[i_part]);
        }
        if (_tgt_to_src_weight[i_part] != NULL) {
          free(_tgt_to_src_weight[i_part]);
        }
        _src_to_tgt_idx   [i_part] = NULL;
        _src_to_tgt_gnum  [i_part] = NULL;
        _src_to_tgt_weight[i_part] = NULL;
        _tgt_to_src_weight[i_part] = NULL;


        if (_weights_idx[i_part] != NULL) {
          free(_weights_idx[i_part]);
        }

        if (_weights[i_part] != NULL) {
          free(_weights[i_part]);
        }

        if (_computed_tgt[i_part] != NULL) {
          free(_computed_tgt[i_part]);
        }

        if (_uncomputed_tgt[i_part] != NULL) {
          free(_uncomputed_tgt[i_part]);
        }

        if (_involved_sources_tgt[i_part] != NULL) {
          free(_involved_sources_tgt[i_part]);
        }

        _n_elt_weights[i_part] = 0;
        _weights_idx  [i_part] = NULL;
        _weights      [i_part] = NULL;

        _n_computed_tgt[i_part] = 0;
        _computed_tgt  [i_part] = NULL;

        _n_uncomputed_tgt[i_part] = 0;
        _uncomputed_tgt  [i_part] = NULL;

        _n_involved_sources_tgt[i_part] = 0;
        _involved_sources_tgt  [i_part] = NULL;
      }

    }

    if (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
      SpatialInterpIntersection *cpl_spatial_interp;

      cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

      if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
        cpl_spatial_interp =
        dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
      }
      else {
        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
        cpl_spatial_interp =
        dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
      }

      for (int i_part = 0; i_part < _nPart; i_part++) {

        if (cpl_spatial_interp->_src_to_tgt_idx[i_part] != NULL) {
          free(cpl_spatial_interp->_src_to_tgt_idx[i_part]);
        }
        if (cpl_spatial_interp->_src_to_tgt_gnum[i_part] != NULL) {
          free(cpl_spatial_interp->_src_to_tgt_gnum[i_part]);
        }
        if (cpl_spatial_interp->_src_to_tgt_weight[i_part] != NULL) {
          free(cpl_spatial_interp->_src_to_tgt_weight[i_part]);
        }
        if (cpl_spatial_interp->_tgt_to_src_weight[i_part] != NULL) {
          free(cpl_spatial_interp->_tgt_to_src_weight[i_part]);
        }
        cpl_spatial_interp->_src_to_tgt_idx   [i_part] = NULL;
        cpl_spatial_interp->_src_to_tgt_gnum  [i_part] = NULL;
        cpl_spatial_interp->_src_to_tgt_weight[i_part] = NULL;
        cpl_spatial_interp->_tgt_to_src_weight[i_part] = NULL;


        if (cpl_spatial_interp->_weights_idx[i_part] != NULL) {
          free(cpl_spatial_interp->_weights_idx[i_part]);
        }

        if (cpl_spatial_interp->_weights[i_part] != NULL) {
          free(cpl_spatial_interp->_weights[i_part]);
        }

        if (cpl_spatial_interp->_computed_tgt[i_part] != NULL) {
          free(cpl_spatial_interp->_computed_tgt[i_part]);
        }

        if (cpl_spatial_interp->_uncomputed_tgt[i_part] != NULL) {
          free(cpl_spatial_interp->_uncomputed_tgt[i_part]);
        }

        if (cpl_spatial_interp->_involved_sources_tgt[i_part] != NULL) {
          free(cpl_spatial_interp->_involved_sources_tgt[i_part]);
        }

        cpl_spatial_interp->_n_elt_weights[i_part] = 0;
        cpl_spatial_interp->_weights_idx  [i_part] = NULL;
        cpl_spatial_interp->_weights      [i_part] = NULL;

        cpl_spatial_interp->_n_computed_tgt[i_part] = 0;
        cpl_spatial_interp->_computed_tgt  [i_part] = NULL;

        cpl_spatial_interp->_n_uncomputed_tgt[i_part] = 0;
        cpl_spatial_interp->_uncomputed_tgt  [i_part] = NULL;

        cpl_spatial_interp->_n_involved_sources_tgt[i_part] = 0;
        cpl_spatial_interp->_involved_sources_tgt  [i_part] = NULL;
      }
    }



    /* Set source (A) and target (B) meshes */

    // get dimension of local and coupled codes
    int local_mesh_dim = 0;
    int cpl_mesh_dim   = 0;

    CWP_Interface_t interf_dim = _cpl->entitiesDimGet();
    switch (interf_dim) {
      case CWP_INTERFACE_POINT:
      local_mesh_dim = 0;
      break;
      case CWP_INTERFACE_LINEAR:
      local_mesh_dim = 1;
      break;
      case CWP_INTERFACE_SURFACE:
      local_mesh_dim = 2;
      break;
      case CWP_INTERFACE_VOLUME:
      local_mesh_dim = 3;
      break;
      default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid interface dimension %d\n", local_mesh_dim);
    }

    cpl_mesh_dim = local_mesh_dim; // ?

    if (local_mesh_dim != cpl_mesh_dim) {
      PDM_error(__FILE__, __LINE__, 0,
                "Both meshes must have the same dimension (got %d / %d)\n",
                local_mesh_dim, cpl_mesh_dim);
    }

    int mesh_dim = local_mesh_dim;

    int n_part_src = 0;
    int n_part_tgt = 0;
    if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
      n_part_src = _nPart;
    }
    else {
      n_part_src = _cplNPart;
    }
    if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
      n_part_tgt = _nPart;
    }
    else {
      n_part_tgt = _cplNPart;
    }

    // A = source
    // B = target

    if (!_coupledCodeProperties->localCodeIs()) {

      _id_pdm = PDM_mesh_intersection_create(PDM_MESH_INTERSECTION_KIND_SOFT,
                                             mesh_dim,
                                             mesh_dim,
                                             0.,
                                             _pdmCplComm,
                                             PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

      // source mesh
      if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
        PDM_mesh_intersection_mesh_nodal_set(_id_pdm, 0, _pdm_CplNodal);
      }
      else {
        // empty mesh
        PDM_mesh_intersection_n_part_set2(_id_pdm, 0, n_part_src);
        // PDM_mesh_intersection_n_part_set(_id_pdm, 0, n_part_src);
        // for (int i = 0; i < n_part_src; i++) {
        //   PDM_mesh_intersection_part_set(_id_pdm,
        //                                  0,     // i_mesh
        //                                  i,     // i_part
        //                                  0,     // n_cell
        //                                  0,     // n_face
        //                                  0,     // n_edge
        //                                  0,     // n_vtx
        //                                  NULL,  // cell_face_idx
        //                                  NULL,  // cell_face
        //                                  NULL,  // face_edge_idx
        //                                  NULL,  // face_edge
        //                                  NULL,  // edge_vtx
        //                                  NULL,  // face_vtx_idx
        //                                  NULL,  // face_vtx
        //                                  NULL,  // cell_ln_to_gn
        //                                  NULL,  // face_ln_to_gn
        //                                  NULL,  // edge_ln_to_gn
        //                                  NULL,  // vtx_ln_to_gn
        //                                  NULL); // vtx_coord
        // }
      }

      // target mesh
      if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
        PDM_mesh_intersection_mesh_nodal_set(_id_pdm, 1, _pdm_CplNodal);
      }
      else {
        // empty mesh
        PDM_mesh_intersection_n_part_set2(_id_pdm, 1, n_part_tgt);
        // PDM_mesh_intersection_n_part_set(_id_pdm, 1, n_part_tgt);
        // for (int i = 0; i < n_part_tgt; i++) {
        //   PDM_mesh_intersection_part_set(_id_pdm,
        //                                  1,     // i_mesh
        //                                  i,     // i_part
        //                                  0,     // n_cell
        //                                  0,     // n_face
        //                                  0,     // n_edge
        //                                  0,     // n_vtx
        //                                  NULL,  // cell_face_idx
        //                                  NULL,  // cell_face
        //                                  NULL,  // face_edge_idx
        //                                  NULL,  // face_edge
        //                                  NULL,  // edge_vtx
        //                                  NULL,  // face_vtx_idx
        //                                  NULL,  // face_vtx
        //                                  NULL,  // cell_ln_to_gn
        //                                  NULL,  // face_ln_to_gn
        //                                  NULL,  // edge_ln_to_gn
        //                                  NULL,  // vtx_ln_to_gn
        //                                  NULL); // vtx_coord
        // }
      }
    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        _id_pdm = PDM_mesh_intersection_create(PDM_MESH_INTERSECTION_KIND_SOFT,
                                               mesh_dim,
                                               mesh_dim,
                                               0.,
                                               _pdmCplComm,
                                               PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        SpatialInterpIntersection *cpl_spatial_interp;
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
          cpl_spatial_interp =
          dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }
        else {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
          cpl_spatial_interp =
          dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }

        cpl_spatial_interp->_id_pdm = _id_pdm;

        cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();
        CWP_UNUSED(cpl_mesh);

        // source mesh
        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
          PDM_mesh_intersection_mesh_nodal_set(_id_pdm, 0, _pdm_CplNodal);
        }
        else {
          PDM_mesh_intersection_mesh_nodal_set(_id_pdm, 0, cpl_mesh->getPdmNodalIndex());
        }

        // target mesh
        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          PDM_mesh_intersection_mesh_nodal_set(_id_pdm, 1, _pdm_CplNodal);
        }
        else {
          PDM_mesh_intersection_mesh_nodal_set(_id_pdm, 1, cpl_mesh->getPdmNodalIndex());
        }
      }
    }


    /* Compute */
    if (!_coupledCodeProperties->localCodeIs() ||
        (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet())) {
      PDM_mesh_intersection_compute(_id_pdm);
      // PDM_mesh_intersection_dump_times(_id_pdm); // not implemented
    }

    /* Reset part_to_part object */
    if (_ptsp != nullptr) {
      if (!_coupledCodeProperties->localCodeIs()) {
        PDM_part_to_part_free(_ptsp);
        _ptsp = nullptr;
      }
      else {

        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
          PDM_part_to_part_free(_ptsp);
          _ptsp = nullptr;

          SpatialInterpIntersection *cpl_spatial_interp;

          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            cpl_spatial_interp =
            dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }
          else {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            cpl_spatial_interp =
            dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }

          cpl_spatial_interp->_ptsp = NULL;
        }

      }
    }

    /* Get PDM part_to_part object */
    if (_id_pdm != NULL) {
      PDM_mesh_intersection_part_to_part_get(_id_pdm,
                                             &_ptsp,
                                             PDM_OWNERSHIP_USER);
    }

    /* Create part_to_part object if null */
    if (_ptsp == NULL) {
      if (!_coupledCodeProperties->localCodeIs()) {
        _ptsp = PDM_part_to_part_create((const PDM_g_num_t **) _src_gnum,
                                        (const int          *) _src_n_gnum,
                                        _nPart,
                                        (const PDM_g_num_t **) _tgt_gnum,
                                        (const int          *) _tgt_n_gnum,
                                        _nPart,
                                        (const int         **) _src_to_tgt_idx,
                                        (const PDM_g_num_t **) _src_to_tgt_gnum,
                                        _pdmCplComm);
      }
      else {
        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
          SpatialInterpIntersection *cpl_spatial_interp;

          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            cpl_spatial_interp =
            dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }
          else {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            cpl_spatial_interp =
            dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }

          if (_ptsp == NULL) {
            if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
              _ptsp = PDM_part_to_part_create((const PDM_g_num_t **) _src_gnum,
                                              (const int          *) _src_n_gnum,
                                              _nPart,
                                              (const PDM_g_num_t **) cpl_spatial_interp->_tgt_gnum,
                                              (const int          *) cpl_spatial_interp->_tgt_n_gnum,
                                              _cplNPart,
                                              (const int         **) _src_to_tgt_idx,
                                              (const PDM_g_num_t **) _src_to_tgt_gnum,
                                              _pdmCplComm);
            }
            else {
              _ptsp = PDM_part_to_part_create((const PDM_g_num_t **) cpl_spatial_interp->_src_gnum,
                                              (const int          *) cpl_spatial_interp->_src_n_gnum,
                                              _cplNPart,
                                              (const PDM_g_num_t **) _tgt_gnum,
                                              (const int          *) _tgt_n_gnum,
                                              _nPart,
                                              (const int         **) cpl_spatial_interp->_src_to_tgt_idx,
                                              (const PDM_g_num_t **) cpl_spatial_interp->_src_to_tgt_gnum,
                                              _pdmCplComm);
            }
          }

          cpl_spatial_interp->_ptsp = _ptsp;
        }
      }
    }

    /* Get PDM results */
    if (!_coupledCodeProperties->localCodeIs()) {
      int  *n_ref_tgt = NULL;
      int **ref_tgt   = NULL;
      PDM_part_to_part_ref_lnum2_get(_ptsp, &n_ref_tgt, &ref_tgt);

      for (int i_part = 0; i_part < _nPart; i_part++) {

        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {

          PDM_mesh_intersection_result_from_a_get(_id_pdm,
                                                  i_part,
                                                  &(_src_to_tgt_idx   [i_part]),
                                                  &(_src_to_tgt_gnum  [i_part]),
                                                  &(_src_to_tgt_weight[i_part]));

          _n_involved_sources_tgt[i_part] = _src_n_gnum[i_part];
          _involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * _n_involved_sources_tgt[i_part]);

          int count = 0;
          for (int i = 0 ; i < _src_n_gnum[i_part] ; ++i) {
            if (_src_to_tgt_idx[i_part][i + 1] > _src_to_tgt_idx[i_part][i]) {
              _involved_sources_tgt[i_part][count] = i + 1;
              ++count;
            }
          }

          _n_involved_sources_tgt[i_part] = count;
          _involved_sources_tgt[i_part] = (int*) realloc(_involved_sources_tgt[i_part], sizeof(int) * count);
        }
        else {
          _src_to_tgt_idx[i_part] = (int *) malloc(sizeof(int));
          _src_to_tgt_idx[i_part][0] = 0;

          PDM_mesh_intersection_result_from_b_get(_id_pdm,
                                                  i_part,
                                                  &(_tgt_to_src_weight[i_part]));

          _n_computed_tgt[i_part] = n_ref_tgt[i_part];
          _computed_tgt[i_part] = (int *) malloc(sizeof(int) * _n_computed_tgt[i_part]);
          memcpy(_computed_tgt[i_part], ref_tgt[i_part], sizeof(int) * _n_computed_tgt[i_part]);
        }
      }
    }
    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        SpatialInterpIntersection *cpl_spatial_interp;

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
          cpl_spatial_interp = dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }
        else {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
          cpl_spatial_interp = dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }

        int  *n_ref_tgt = NULL;
        int **ref_tgt   = NULL;
        PDM_part_to_part_ref_lnum2_get(_ptsp, &n_ref_tgt, &ref_tgt);

        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
          for (int i_part = 0; i_part < _nPart; i_part++) {
            PDM_mesh_intersection_result_from_a_get(_id_pdm,
                                                    i_part,
                                                    &(_src_to_tgt_idx   [i_part]),
                                                    &(_src_to_tgt_gnum  [i_part]),
                                                    &(_src_to_tgt_weight[i_part]));

            _n_involved_sources_tgt[i_part] = _src_n_gnum[i_part];
            _involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * _n_involved_sources_tgt[i_part]);

            int count = 0;
            for (int i = 0 ; i < _src_n_gnum[i_part] ; ++i) {
              if (_src_to_tgt_idx[i_part][i + 1] > _src_to_tgt_idx[i_part][i]) {
                _involved_sources_tgt[i_part][count] = i + 1;
                ++count;
              }
            }

            _n_involved_sources_tgt[i_part] = count;
            _involved_sources_tgt[i_part] = (int*) realloc(_involved_sources_tgt[i_part], sizeof(int) * count);
          }

          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            cpl_spatial_interp->_src_to_tgt_idx[i_part] = (int *) malloc(sizeof(int));
            cpl_spatial_interp->_src_to_tgt_idx[i_part][0] = 0;

            PDM_mesh_intersection_result_from_b_get(_id_pdm,
                                                    i_part,
                                                    &(cpl_spatial_interp->_tgt_to_src_weight[i_part]));

            cpl_spatial_interp->_n_computed_tgt[i_part] = n_ref_tgt[i_part];
            cpl_spatial_interp->_computed_tgt[i_part] = (int *) malloc(sizeof(int) * cpl_spatial_interp->_n_computed_tgt[i_part]);
            memcpy(cpl_spatial_interp->_computed_tgt[i_part], ref_tgt[i_part], sizeof(int) * cpl_spatial_interp->_n_computed_tgt[i_part]);
          }
        }
        else {
          for (int i_part = 0; i_part < _nPart; i_part++) {
            _src_to_tgt_idx[i_part] = (int *) malloc(sizeof(int));
            _src_to_tgt_idx[i_part][0] = 0;

            PDM_mesh_intersection_result_from_b_get(_id_pdm,
                                                    i_part,
                                                    &(_tgt_to_src_weight[i_part]));

            _n_computed_tgt[i_part] = n_ref_tgt[i_part];
            _computed_tgt[i_part] = (int *) malloc(sizeof(int) * _n_computed_tgt[i_part]);
            memcpy(_computed_tgt[i_part], ref_tgt[i_part], sizeof(int) * _n_computed_tgt[i_part]);
          }

          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            PDM_mesh_intersection_result_from_a_get(_id_pdm,
                                                    i_part,
                                                    &(cpl_spatial_interp->_src_to_tgt_idx   [i_part]),
                                                    &(cpl_spatial_interp->_src_to_tgt_gnum  [i_part]),
                                                    &(cpl_spatial_interp->_src_to_tgt_weight[i_part]));

            cpl_spatial_interp->_n_involved_sources_tgt[i_part] = cpl_spatial_interp->_src_n_gnum[i_part];
            cpl_spatial_interp->_involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * cpl_spatial_interp->_n_involved_sources_tgt[i_part]);

            int count = 0;
            for (int i = 0 ; i < cpl_spatial_interp->_src_n_gnum[i_part] ; ++i) {
              if (cpl_spatial_interp->_src_to_tgt_idx[i_part][i + 1] > cpl_spatial_interp->_src_to_tgt_idx[i_part][i]) {
                cpl_spatial_interp->_involved_sources_tgt[i_part][count] = i + 1;
                ++count;
              }
            }

            cpl_spatial_interp->_n_involved_sources_tgt[i_part] = count;
            cpl_spatial_interp->_involved_sources_tgt[i_part] = (int*) realloc(cpl_spatial_interp->_involved_sources_tgt[i_part], sizeof(int) * count);
          }
        }
      }
    }


    /* Free PDM object */
    if (!_coupledCodeProperties->localCodeIs()) {
      PDM_mesh_intersection_free(_id_pdm);
      _id_pdm = nullptr;
    }
    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        PDM_mesh_intersection_free(_id_pdm);
        _id_pdm = nullptr;

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        SpatialInterpIntersection *cpl_spatial_interp;

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
          cpl_spatial_interp = dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }
        else {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
          cpl_spatial_interp = dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }
        cpl_spatial_interp->_id_pdm = nullptr;
      }
    }
  }

  void SpatialInterpIntersection::interpolate(Field *referenceField, double **buffer) {

    int nComponent = referenceField->nComponentGet();
    CWP_Dof_location_t referenceFieldType = referenceField->locationGet();
    CWP_Interpolation_t interpolationType = referenceField->interpolationTypeGet();
    const CWP_Field_storage_t storage = referenceField->storageTypeGet();

    if (interpolationType == CWP_INTERPOLATION_USER) {
      PDM_error(__FILE__, __LINE__, 0, "user interpolation not implemented yet");
    }

    else {

      int  *n_ref_tgt = NULL;
      int **ref_tgt   = NULL;
      PDM_part_to_part_ref_lnum2_get(_ptsp, &n_ref_tgt, &ref_tgt);

      int         **come_from_idx = NULL;
      PDM_g_num_t **come_from     = NULL;
      PDM_part_to_part_gnum1_come_from_get(_ptsp,
                                           &come_from_idx,
                                           &come_from);

      for (int i_part = 0; i_part < _nPart; i_part++) {

        double *referenceData = (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_TARGET);

        int n_tgt = 0;
        if (referenceFieldType == CWP_DOF_LOCATION_CELL_CENTER) {
          n_tgt = _mesh->getPartNElts(i_part);
        }
        else if (referenceFieldType == CWP_DOF_LOCATION_NODE) {
          n_tgt = _mesh->getPartNVertex(i_part);
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "user tgt not supported yet");
        }

        for (int j = 0; j < nComponent*n_tgt; j++) {
          referenceData[j] = 0;
        }

        for (int iref = 0; iref < n_ref_tgt[i_part]; iref++) {

          int i = ref_tgt[i_part][iref] - 1;

          for (int k = come_from_idx[i_part][iref]; k < come_from_idx[i_part][iref+1]; k++) {

            double w = _tgt_to_src_weight[i_part][k];

            if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
              for (int j = 0; j < nComponent; j++) {
                referenceData[n_tgt*j + i] += w*buffer[i_part][come_from_idx[i_part][n_ref_tgt[i_part]]*j + k];
              }
            }
            else {
              for (int j = 0; j < nComponent; j++) {
                referenceData[nComponent*i + j] += w*buffer[i_part][nComponent*k + j];
              }
            }

          }

        }

      }
    }
  }
};
