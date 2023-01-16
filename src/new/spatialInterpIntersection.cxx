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

namespace cwipi {
  SpatialInterpIntersection::SpatialInterpIntersection() = default;

  SpatialInterpIntersection::~SpatialInterpIntersection
  (
   )
  {
    for (int i_part = 0; i_part < _nPart; i_part++) {
      if (_tgt_to_src_idx[i_part] != NULL) {
        free(_tgt_to_src_idx[i_part]);
      }
      if (_tgt_to_src_gnum[i_part] != NULL) {
        free(_tgt_to_src_gnum[i_part]);
      }
      if (_tgt_to_src_weight[i_part] != NULL) {
        free(_tgt_to_src_weight[i_part]);
      }
    }

    delete[] _tgt_to_src_idx;
    delete[] _tgt_to_src_gnum;
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
      // Data for PDM_part_to_part_t
      if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
        for (int i_part = 0 ; i_part < _nPart ; i_part++) {
          if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            _src_gnum  [i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
            _src_n_gnum[i_part] = _mesh->getPartNElts (i_part);

          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
            _src_gnum  [i_part] = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
            _src_n_gnum[i_part] = _mesh->getPartNVertex (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            // ???
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
            _tgt_gnum  [i_part] = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
            _tgt_n_gnum[i_part] = _mesh->getPartNVertex (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            // ???
            _tgt_gnum  [i_part] = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
            _tgt_n_gnum[i_part] = _cpl->userTargetNGet (i_part);
          }
        }
      }

      //
      // Target properties
      _tgt_to_src_idx    = new int*         [_nPart];
      _tgt_to_src_gnum   = new PDM_g_num_t* [_nPart];
      _tgt_to_src_weight = new double*      [_nPart];

      for (int i_part = 0; i_part < _nPart; i_part++) {
        _tgt_to_src_idx   [i_part] = NULL;
        _tgt_to_src_gnum  [i_part] = NULL;
        _tgt_to_src_weight[i_part] = NULL;
      }
    }

  void SpatialInterpIntersection::weightsCompute() {

    if (!_coupledCodeProperties->localCodeIs() ||
        (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet())) {

      for (int i_part = 0; i_part < _nPart; i_part++) {

        if (_tgt_to_src_idx[i_part] != NULL) {
          free(_tgt_to_src_idx[i_part]);
        }
        if (_tgt_to_src_gnum[i_part] != NULL) {
          free(_tgt_to_src_gnum[i_part]);
        }
        if (_tgt_to_src_weight[i_part] != NULL) {
          free(_tgt_to_src_weight[i_part]);
        }
        _tgt_to_src_idx   [i_part] = NULL;
        _tgt_to_src_gnum  [i_part] = NULL;
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

        if (cpl_spatial_interp->_tgt_to_src_idx[i_part] != NULL) {
          free(cpl_spatial_interp->_tgt_to_src_idx[i_part]);
        }
        if (cpl_spatial_interp->_tgt_to_src_gnum[i_part] != NULL) {
          free(cpl_spatial_interp->_tgt_to_src_gnum[i_part]);
        }
        if (cpl_spatial_interp->_tgt_to_src_weight[i_part] != NULL) {
          free(cpl_spatial_interp->_tgt_to_src_weight[i_part]);
        }
        cpl_spatial_interp->_tgt_to_src_idx   [i_part] = NULL;
        cpl_spatial_interp->_tgt_to_src_gnum  [i_part] = NULL;
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



    /* Set source (B) and target (A) meshes */
    /* CWIPI stores a nodal mesh but pdm_mesh_intersection requires ngons... */

    //...


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

    /* Get PDM results */
    if (!_coupledCodeProperties->localCodeIs()) {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
          _tgt_to_src_idx[i_part] = (int *) malloc(sizeof(int));
          _tgt_to_src_idx[i_part][0] = 0;

          // involved sources...
          // computed tgt...
        }
        else {
          PDM_mesh_intersection_result_from_a_get(_id_pdm,
                                                  i_part,
                                                  &(_tgt_to_src_idx   [i_part]),
                                                  &(_tgt_to_src_gnum  [i_part]),
                                                  &(_tgt_to_src_weight[i_part]));

          // involved sources...
          // computed tgt...
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

        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
          for (int i_part = 0; i_part < _nPart; i_part++) {
            _tgt_to_src_idx[i_part] = (int *) malloc(sizeof(int));
            _tgt_to_src_idx[i_part][0] = 0;

            // involved sources...
            // computed tgt...
          }

          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            PDM_mesh_intersection_result_from_a_get(_id_pdm,
                                                    i_part,
                                                    &(cpl_spatial_interp->_tgt_to_src_idx   [i_part]),
                                                    &(cpl_spatial_interp->_tgt_to_src_gnum  [i_part]),
                                                    &(cpl_spatial_interp->_tgt_to_src_weight[i_part]));

            // involved sources...
            // computed tgt...
          }
        }
        else {
          for (int i_part = 0; i_part < _nPart; i_part++) {
            PDM_mesh_intersection_result_from_a_get(_id_pdm,
                                                    i_part,
                                                    &(_tgt_to_src_idx   [i_part]),
                                                    &(_tgt_to_src_gnum  [i_part]),
                                                    &(_tgt_to_src_weight[i_part]));

            // involved sources...
            // computed tgt...
          }

          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            cpl_spatial_interp->_tgt_to_src_idx[i_part] = (int *) malloc(sizeof(int));
            cpl_spatial_interp->_tgt_to_src_idx[i_part][0] = 0;

            // involved sources...
            // computed tgt...
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

    /* Create part_to_part object if null */
    if (_ptsp == NULL) {
      if (!_coupledCodeProperties->localCodeIs()) {
        _ptsp = PDM_part_to_part_create((const PDM_g_num_t **) _tgt_gnum,
                                        (const int          *) _tgt_n_gnum,
                                        _nPart,
                                        (const PDM_g_num_t **) _src_gnum,
                                        (const int          *) _src_n_gnum,
                                        _nPart,
                                        (const int         **) _tgt_to_src_idx,
                                        (const PDM_g_num_t **) _tgt_to_src_gnum,
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
              _ptsp = PDM_part_to_part_create((const PDM_g_num_t **) _tgt_gnum,
                                              (const int          *) _tgt_n_gnum,
                                              _nPart,
                                              (const PDM_g_num_t **) cpl_spatial_interp->_src_gnum,
                                              (const int          *) cpl_spatial_interp->_src_n_gnum,
                                              _cplNPart,
                                              (const int         **) _tgt_to_src_idx,
                                              (const PDM_g_num_t **) _tgt_to_src_gnum,
                                              _pdmCplComm);
            }
            else {
              _ptsp = PDM_part_to_part_create((const PDM_g_num_t **) cpl_spatial_interp->_tgt_gnum,
                                              (const int          *) cpl_spatial_interp->_tgt_n_gnum,
                                              _cplNPart,
                                              (const PDM_g_num_t **) _src_gnum,
                                              (const int          *) _src_n_gnum,
                                              _nPart,
                                              (const int         **) cpl_spatial_interp->_tgt_to_src_idx,
                                              (const PDM_g_num_t **) cpl_spatial_interp->_tgt_to_src_gnum,
                                              _pdmCplComm);
            }
          }

          cpl_spatial_interp->_ptsp = _ptsp;
        }
      }
    }
  }


  void SpatialInterpIntersection::issend(Field *referenceField) {

    if (referenceField->currentStepWasExchangedGet()) {
      PDM_error(__FILE__, __LINE__, 0,
                "The field has already been exchanged for the current time step "
                "(CWP_Time_update must be called before the exchange)\n");
    }

    if (!_coupledCodeProperties->localCodeIs()) {

      const int intId            = referenceField->fieldIDIntGet();
      const CWP_Type_t data_type = referenceField->dataTypeGet();
      CWP_UNUSED (data_type);
      const size_t s_data        = sizeof(double);
      const int stride           = referenceField->nComponentGet();
      const CWP_Field_storage_t storage = referenceField->storageTypeGet();
      PDM_stride_t pdm_storage;
      if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
        pdm_storage = PDM_STRIDE_CST_INTERLEAVED;
      } else {
        pdm_storage = PDM_STRIDE_CST_INTERLACED;
      }

      int  n_part_tgt;
      int  n_part_src;
      int *n_tgt = NULL;
      int *n_src = NULL;
      PDM_part_to_part_n_part_and_n_elt_get(_ptsp,
                                            &n_part_tgt,
                                            &n_part_src,
                                            &n_tgt,
                                            &n_src);

      int         **come_from_idx = NULL;
      PDM_g_num_t **come_from     = NULL;
      PDM_part_to_part_gnum1_come_from_get(_ptsp,
                                           &come_from_idx,
                                           &come_from);

      _send_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);
      for (int i = 0; i < _nPart; i++) {
        double *referenceData = (double *) referenceField->dataGet(i, CWP_FIELD_MAP_SOURCE);
        _send_buffer[intId][i] = (double *) malloc(sizeof(double) * stride * n_src[i]);

        if (storage == CWP_FIELD_STORAGE_INTERLACED) {
          for (int j = 0; j < n_src[i]; j++) {
            for (int k = come_from_idx[i][j]; k < come_from_idx[i][j+1]; k++) {
              for (int l = 0; l < stride; l++) {
                _send_buffer[intId][i][stride*k + l] = referenceData[stride*j + l];
              }
            }
          }
        }
        else {
          for (int l = 0; l < stride; l++) {
            for (int j = 0; j < n_src[i]; j++) {
              for (int k = come_from_idx[i][j]; k < come_from_idx[i][j+1]; k++) {
                _send_buffer[intId][i][come_from_idx[i][n_src[i]]*l + k] = referenceData[n_src[i]*l + j];
              }
            }
          }
        }
      }


      PDM_part_to_part_reverse_iexch(_ptsp,
                                     PDM_MPI_COMM_KIND_P2P,
                                     pdm_storage,
                                     PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                     stride,
                                     s_data,
                                     NULL,
                    (const void  **) _send_buffer[intId],
                                     NULL,
                    (      void ***) &_recv_buffer[intId],
                                     &(_send_request[intId]));
      _recv_request[intId] = _send_request[intId];
    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        const int intId            = referenceField->fieldIDIntGet();
        const CWP_Type_t data_type = referenceField->dataTypeGet();
        CWP_UNUSED(data_type);
        const size_t s_data        = sizeof(double);
        const int stride           = referenceField->nComponentGet();
        const CWP_Field_storage_t storage = referenceField->storageTypeGet();
        PDM_stride_t pdm_storage;
        if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
          pdm_storage = PDM_STRIDE_CST_INTERLEAVED;
        } else {
          pdm_storage = PDM_STRIDE_CST_INTERLACED;
        }

        int  n_part_tgt;
        int  n_part_src;
        int *n_tgt = NULL;
        int *n_src = NULL;
        PDM_part_to_part_n_part_and_n_elt_get(_ptsp,
                                              &n_part_tgt,
                                              &n_part_src,
                                              &n_tgt,
                                              &n_src);

        int         **come_from_idx = NULL;
        PDM_g_num_t **come_from     = NULL;
        PDM_part_to_part_gnum1_come_from_get(_ptsp,
                                             &come_from_idx,
                                             &come_from);

        _send_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);
        for (int i = 0; i < _nPart; i++) {
          double *referenceData = (double *) referenceField->dataGet(i, CWP_FIELD_MAP_SOURCE);
          _send_buffer[intId][i] = (double *) malloc(sizeof(double) * stride * n_src[i]);

          if (storage == CWP_FIELD_STORAGE_INTERLACED) {
            for (int j = 0; j < n_src[i]; j++) {
              for (int k = come_from_idx[i][j]; k < come_from_idx[i][j+1]; k++) {
                for (int l = 0; l < stride; l++) {
                  _send_buffer[intId][i][stride*k + l] = referenceData[stride*j + l];
                }
              }
            }
          }
          else {
            for (int l = 0; l < stride; l++) {
              for (int j = 0; j < n_src[i]; j++) {
                for (int k = come_from_idx[i][j]; k < come_from_idx[i][j+1]; k++) {
                  _send_buffer[intId][i][come_from_idx[i][n_src[i]]*l + k] = referenceData[n_src[i]*l + j];
                }
              }
            }
          }
        }


        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
        SpatialInterpIntersection *cpl_spatial_interp =
        dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];

        const int cpl_intId = cpl_referenceField->fieldIDIntGet();
        cpl_spatial_interp->_send_buffer[cpl_intId] = (double **) malloc(sizeof(double *) * _cplNPart);
        cpl_spatial_interp->_recv_buffer[cpl_intId] = (double **) malloc(sizeof(double *) * _cplNPart);

        for (int i = 0; i < _cplNPart; i++) {
          cpl_spatial_interp->_send_buffer[cpl_intId][i] = nullptr;
        }
        PDM_part_to_part_reverse_iexch(_ptsp,
                                       PDM_MPI_COMM_KIND_P2P,
                                       pdm_storage,
                                       PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                       stride,
                                       s_data,
                                       NULL,
                       (const void **) _send_buffer[intId],
                                       NULL,
                            (void ***) &cpl_spatial_interp->_recv_buffer[cpl_intId],
                                       &(_send_request[intId]));
        cpl_spatial_interp->_recv_request[cpl_intId] = _send_request[intId];
      }
    }
  }

  void SpatialInterpIntersection::waitIssend(Field* referenceField) {
    if (!_coupledCodeProperties->localCodeIs()) {

      const int intId = referenceField->fieldIDIntGet();

      PDM_part_to_part_reverse_iexch_wait (_ptsp, _send_request[intId]);

      if (_send_buffer[intId] != NULL) {
        for (int i = 0; i < _nPart; i++) {
          if (_send_buffer[intId][i] != NULL) {
            free (_send_buffer[intId][i]);
            _send_buffer[intId][i] = NULL;
          }
        }
        free (_send_buffer[intId]);
        _send_buffer[intId] = NULL;
      }

      if (_recv_buffer[intId] != NULL) {
        for (int i = 0; i < _nPart; i++) {
          if (_recv_buffer[intId][i] != NULL) {
            free (_recv_buffer[intId][i]);
            _recv_buffer[intId][i] = NULL;
          }
        }
        free (_recv_buffer[intId]);
        _recv_buffer[intId] = NULL;
      }

      // if(_visu->isCreated() && referenceField->visuStatusGet() == CWP_STATUS_ON) {
      //   _visu->WriterField(referenceField, nullptr, nullptr, CWP_FIELD_MAP_SOURCE);
      // }
    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        const int intId            = referenceField->fieldIDIntGet();

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
        SpatialInterpIntersection *cpl_spatial_interp =
        dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];

        const int cpl_intId = cpl_referenceField->fieldIDIntGet();

        PDM_part_to_part_reverse_iexch_wait(_ptsp, cpl_spatial_interp->_recv_request[cpl_intId]);

        for (int i = 0; i < _nPart; i++) {
          if (_send_buffer[intId] != NULL) {
            if (_send_buffer[intId][i] != NULL) {
              free (_send_buffer[intId][i]);
              _send_buffer[intId][i] = NULL;
            }
            free (_send_buffer[intId]);
            _send_buffer[intId] = NULL;
          }
          if (_recv_buffer[intId] != NULL) {
            if (_recv_buffer[intId][i] != NULL) {
              free (_recv_buffer[intId][i]);
              _recv_buffer[intId][i] = NULL;
            }
            free (_recv_buffer[intId]);
            _recv_buffer[intId] = NULL;
          }
        }

        int nComponent   = cpl_referenceField->nComponentGet();
        int dataTypeSize = cpl_referenceField->dataTypeSizeGet();

        int          *n_tgt;
        int         **tgt_to_src_idx;
        PDM_g_num_t **tgt_to_src;
        PDM_part_to_part_part1_to_part2_get(_ptsp,
                                            &n_tgt,
                                            &tgt_to_src_idx,
                                            &tgt_to_src);

        for (int i = 0; i < _cplNPart; i++) {
          double *referenceData  = (double *) cpl_referenceField->dataGet(i, CWP_FIELD_MAP_TARGET);
          for (int j = 0; j < n_tgt[i]; j++) {
            assert ((tgt_to_src[i][j+1] - tgt_to_src[i][j]) == 1); // why??
          }
          memcpy(referenceData, cpl_spatial_interp->_recv_buffer[cpl_intId][i], dataTypeSize * nComponent * n_tgt[i]);
        }

        for (int i = 0; i < _cplNPart; i++) {
          if (cpl_spatial_interp->_send_buffer[cpl_intId] != NULL) {
            if (cpl_spatial_interp->_send_buffer[cpl_intId][i] != NULL) {
              free (cpl_spatial_interp->_send_buffer[cpl_intId][i]);
              cpl_spatial_interp->_send_buffer[cpl_intId][i] = NULL;
            }
            free (cpl_spatial_interp->_send_buffer[cpl_intId]);
            cpl_spatial_interp->_send_buffer[cpl_intId] = NULL;
          }
          if (cpl_spatial_interp->_recv_buffer[cpl_intId] != NULL) {
            if (cpl_spatial_interp->_recv_buffer[cpl_intId][i] != NULL) {
              free (cpl_spatial_interp->_recv_buffer[cpl_intId][i]);
              cpl_spatial_interp->_recv_buffer[cpl_intId][i] = NULL;
            }
            free (cpl_spatial_interp->_recv_buffer[cpl_intId]);
            cpl_spatial_interp->_recv_buffer[cpl_intId] = NULL;
          }

          // if(cpl_spatial_interp->_visu->isCreated() && cpl_referenceField->visuStatusGet() == CWP_STATUS_ON) {
          //   cpl_spatial_interp->_visu->WriterField(cpl_referenceField, ptp2_n_ref_gnum2, ptp2_ref_gnum2, CWP_FIELD_MAP_TARGET);
          // }

        }
      }
    }
  }

  void SpatialInterpIntersection::irecv(Field *referenceField) {

    if (referenceField->currentStepWasExchangedGet()) {
      PDM_error(__FILE__, __LINE__, 0,
                "The field has already been exchanged for the current time step "
                "(CWP_Time_update must be called before the exchange)\n");
    }

    if (!_coupledCodeProperties->localCodeIs()) {

      const int intId            = referenceField->fieldIDIntGet();
      const CWP_Type_t data_type = referenceField->dataTypeGet();
      CWP_UNUSED(data_type);
      const size_t s_data        = sizeof(double);
      const int stride           = referenceField->nComponentGet();
      const CWP_Field_storage_t storage = referenceField->storageTypeGet();
      PDM_stride_t pdm_storage;
      if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
        pdm_storage = PDM_STRIDE_CST_INTERLEAVED;
      } else {
        pdm_storage = PDM_STRIDE_CST_INTERLACED;
      }

      _send_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);

      for (int i = 0; i < _nPart; i++) {
        _send_buffer[intId][i] = nullptr;
      }

      PDM_part_to_part_reverse_iexch(_ptsp,
                                     PDM_MPI_COMM_KIND_P2P,
                                     pdm_storage,
                                     PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                     stride,
                                     s_data,
                                     NULL,
                     (const void **) _send_buffer[intId],
                                     NULL,
                          (void ***) &_recv_buffer[intId],
                                     &(_send_request[intId]));
      _recv_request[intId] = _send_request[intId];
    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        const int intId            = referenceField->fieldIDIntGet();
        const CWP_Type_t data_type = referenceField->dataTypeGet();
        CWP_UNUSED(data_type);
        const size_t s_data        = sizeof(double);
        const int stride           = referenceField->nComponentGet();
        const CWP_Field_storage_t storage = referenceField->storageTypeGet();
        PDM_stride_t pdm_storage;
        if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
          pdm_storage = PDM_STRIDE_CST_INTERLEAVED;
        } else {
          pdm_storage = PDM_STRIDE_CST_INTERLACED;
        }

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
        SpatialInterpIntersection *cpl_spatial_interp =
        dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];

        const int cpl_intId = cpl_referenceField->fieldIDIntGet();


        int  cpl_n_part_tgt;
        int  cpl_n_part_src;
        int *cpl_n_tgt = NULL;
        int *cpl_n_src = NULL;
        PDM_part_to_part_n_part_and_n_elt_get(_ptsp,
                                              &cpl_n_part_tgt,
                                              &cpl_n_part_src,
                                              &cpl_n_tgt,
                                              &cpl_n_src);

        int         **cpl_come_from_idx = NULL;
        PDM_g_num_t **cpl_come_from     = NULL;
        PDM_part_to_part_gnum1_come_from_get(_ptsp,
                                             &cpl_come_from_idx,
                                             &cpl_come_from);


        cpl_spatial_interp->_send_buffer[cpl_intId] = (double **) malloc(sizeof(double *) * _cplNPart);

        for (int i = 0; i < _cplNPart; i++) {
          double *cpl_referenceData = (double *) cpl_referenceField->dataGet(i, CWP_FIELD_MAP_SOURCE);

          cpl_spatial_interp->_send_buffer[intId][i] = (double *) malloc(sizeof(double) * stride * cpl_n_src[i]);

          if (storage == CWP_FIELD_STORAGE_INTERLACED) {
            for (int j = 0; j < cpl_n_src[i]; j++) {
              for (int k = cpl_come_from_idx[i][j]; k < cpl_come_from_idx[i][j+1]; k++) {
                for (int l = 0; l < stride; l++) {
                  cpl_spatial_interp->_send_buffer[intId][i][stride*k + l] = cpl_referenceData[stride*j + l];
                }
              }
            }
          }
          else {
            for (int l = 0; l < stride; l++) {
              for (int j = 0; j < cpl_n_src[i]; j++) {
                for (int k = cpl_come_from_idx[i][j]; k < cpl_come_from_idx[i][j+1]; k++) {
                  cpl_spatial_interp->_send_buffer[intId][i][cpl_come_from_idx[i][cpl_n_src[i]]*l + k] = cpl_referenceData[cpl_n_src[i]*l + j];
                }
              }
            }
          }
        }

        PDM_part_to_part_reverse_iexch(_ptsp,
                                       PDM_MPI_COMM_KIND_P2P,
                                       pdm_storage,
                                       PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                       stride,
                                       s_data,
                                       NULL,
                       (const void **) cpl_spatial_interp->_send_buffer[cpl_intId],
                                       NULL,
                            (void ***) &_recv_buffer[intId],
                                       &(cpl_spatial_interp->_send_request[cpl_intId]));
        _recv_request[intId] = cpl_spatial_interp->_send_request[cpl_intId];
      }
    }
  }

  void SpatialInterpIntersection::waitIrecv(Field *referenceField) {

    if (!_coupledCodeProperties->localCodeIs()) {

      const int intId = referenceField->fieldIDIntGet();

      PDM_part_to_part_reverse_iexch_wait(_ptsp, _send_request[intId]);

      interpolate(referenceField, _recv_buffer[intId]);

      if (_send_buffer[intId] != NULL) {
        for (int i = 0; i < _nPart; i++) {
          if (_send_buffer[intId][i] != NULL) {
            free (_send_buffer[intId][i]);
            _send_buffer[intId][i] = NULL;
          }
        }
        free (_send_buffer[intId]);
        _send_buffer[intId] = NULL;
      }

      if (_recv_buffer[intId] != NULL) {
        for (int i = 0; i < _nPart; i++) {
          if (_recv_buffer[intId][i] != NULL) {
            free (_recv_buffer[intId][i]);
            _recv_buffer[intId][i] = NULL;
          }
        }
        free (_recv_buffer[intId]);
        _recv_buffer[intId] = NULL;
      }

      // if(_visu->isCreated() && referenceField->visuStatusGet() == CWP_STATUS_ON) {

      //   int  *ptp2_n_ref_gnum2;
      //   int **ptp2_ref_gnum2;
      //   PDM_part_to_part_ref_lnum2_get (_ptsp,
      //                                  &ptp2_n_ref_gnum2,
      //                                  &ptp2_ref_gnum2);

      //   _visu->WriterField(referenceField, ptp2_n_ref_gnum2, ptp2_ref_gnum2, CWP_FIELD_MAP_TARGET);

      // }
    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        const int intId = referenceField->fieldIDIntGet();

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.sendSpatialInterpGet();
        SpatialInterpIntersection *cpl_spatial_interp =
        dynamic_cast <SpatialInterpIntersection *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];

        const int cpl_intId = cpl_referenceField->fieldIDIntGet();

        PDM_part_to_part_reverse_iexch_wait(_ptsp, cpl_spatial_interp->_send_request[cpl_intId]);

        interpolate(referenceField, _recv_buffer[intId]);

        for (int i = 0; i < _cplNPart; i++) {
          if (cpl_spatial_interp->_send_buffer[intId] != NULL) {
            if (cpl_spatial_interp->_send_buffer[intId][i] != NULL) {
              free (cpl_spatial_interp->_send_buffer[intId][i]);
              cpl_spatial_interp->_send_buffer[intId][i] = NULL;
            }
            free (cpl_spatial_interp->_send_buffer[intId]);
            cpl_spatial_interp->_send_buffer[intId] = NULL;
          }
          if (cpl_spatial_interp->_recv_buffer[intId] != NULL) {
            if (cpl_spatial_interp->_recv_buffer[intId][i] != NULL) {
              free (cpl_spatial_interp->_recv_buffer[intId][i]);
              cpl_spatial_interp->_recv_buffer[intId][i] = NULL;
            }
            free (cpl_spatial_interp->_send_buffer[intId]);
            cpl_spatial_interp->_send_buffer[intId] = NULL;
          }
        }

        for (int i = 0; i < _nPart; i++) {
          if (_send_buffer[intId] != NULL) {
            if (_send_buffer[intId][i] != NULL) {
              free (_send_buffer[intId][i]);
              _send_buffer[intId][i] = NULL;
            }
            free (_send_buffer[intId]);
            _send_buffer[intId] = NULL;
          }
          if (_recv_buffer[intId] != NULL) {
            if (_recv_buffer[intId][i] != NULL) {
              free (_recv_buffer[intId][i]);
              _recv_buffer[intId][i] = NULL;
            }
            free (_recv_buffer[intId]);
            _recv_buffer[intId] = NULL;
          }
        }

        // if(_visu->isCreated() && referenceField->visuStatusGet() == CWP_STATUS_ON) {

        //   int  *ptp2_n_ref_gnum2;
        //   int **ptp2_ref_gnum2;
        //   PDM_part_to_part_ref_lnum2_get (_ptsp,
        //                                  &ptp2_n_ref_gnum2,
        //                                  &ptp2_ref_gnum2);
        //   _visu->WriterField(referenceField, ptp2_n_ref_gnum2, ptp2_ref_gnum2, CWP_FIELD_MAP_TARGET);
        // }

        // if(cpl_spatial_interp->_visu->isCreated() && cpl_referenceField->visuStatusGet() == CWP_STATUS_ON) {
        //   cpl_spatial_interp->_visu->WriterField(cpl_referenceField, nullptr, nullptr, CWP_FIELD_MAP_SOURCE);
        // }
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

        for (int i = 0; i < n_tgt; i++) {

          for (int k = _tgt_to_src_idx[i_part][i]; k < _tgt_to_src_idx[i_part][i+1]; k++) {

            double w = _tgt_to_src_weight[i_part][k];

            if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
              for (int j = 0; j < nComponent; j++) {
                referenceData[n_tgt*j + i] += w*buffer[i_part][_tgt_to_src_idx[i_part][n_tgt]*j + k];
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
