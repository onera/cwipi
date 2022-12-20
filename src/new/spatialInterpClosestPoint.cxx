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

#include "pdm_closest_points.h"
#include "spatialInterpClosestPoint.hxx"
#include "cwp_priv.h"
#include "coupling.hxx"
#include "coupling_i.hxx"

CWP_CLANG_SUPPRESS_WARNING("-Wunused-private-field")


namespace cwipi {
    //CWIPI_CLANG_SUPPRESS_WARNING("-Wunused-private-field")
    SpatialInterpClosestPoint::SpatialInterpClosestPoint() = default;

    // SpatialInterpClosestPoint::~SpatialInterpClosestPoint() = default;

    SpatialInterpClosestPoint::~SpatialInterpClosestPoint
    (
     )
    {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_src_to_tgt_idx[i_part] != NULL) {
          free(_src_to_tgt_idx[i_part]);
        }
        if (_src_to_tgt[i_part] != NULL) {
          free(_src_to_tgt[i_part]);
        }
        if (_src_to_tgt_dist[i_part] != NULL) {
          free(_src_to_tgt_dist[i_part]);
        }
      }


      delete[] _src_to_tgt_idx;
      delete[] _src_to_tgt;
      delete[] _src_to_tgt_dist;
    }


    /**
     *
     * \brief SpatialInterpClosestPoint Init.
     *
     */

    void
    SpatialInterpClosestPoint::init
    (
     Coupling                   *coupling,
     CWP_Dof_location_t          localCodeDofLOcation,
     CWP_Dof_location_t          cplCodeDofLOcation,
     SpatialInterpExchDirection  exchDirection
     )
    {
      SpatialInterp::init (coupling,
                           localCodeDofLOcation,
                           cplCodeDofLOcation,
                           exchDirection);

      _interpolation_time = CWP_SPATIAL_INTERP_AT_SEND;

      //
      // Data for PDM_part_to_part_t
      // ?

      //
      // Target properties
      // not really useful, do we keep this?

      //
      // Source properties
      _src_to_tgt_idx  = new int*         [_nPart];
      _src_to_tgt      = new PDM_g_num_t* [_nPart];
      _src_to_tgt_dist = new double*      [_nPart];

      for (int i_part = 0; i_part < _nPart; i_part++) {
        _src_to_tgt_idx [i_part] = NULL;
        _src_to_tgt     [i_part] = NULL;
        _src_to_tgt_dist[i_part] = NULL;
      }

    }

    void SpatialInterpClosestPoint::weightsCompute() {
//         // In case of withOutPart the user provided not null data only on the root rank (senderRank).
//         if (!_both_codes_are_local) {
//             if (_Texch_t == CWP_FIELD_EXCH_RECV && _pointsCloudLocation == CWP_DOF_LOCATION_USER) user_targets_gnum_compute();
//         }
//         else {
//             if (_Texch_t == CWP_FIELD_EXCH_SEND) {
//                 _spatial_interp_cpl->_Texch_t = CWP_FIELD_EXCH_RECV;
//                 if (_pointsCloudLocation == CWP_DOF_LOCATION_USER) _spatial_interp_cpl->user_targets_gnum_compute();
//             }
//         }

//         // Get informations about the local and the coupled meshes
//         info_mesh();

//         _id_pdm = PDM_closest_points_create(_pdm_cplComm, 5, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);
// //        cout << cplComm_rank << ": PDM_closest_point created with id " << _id_pdm << endl;

//         PDM_closest_points_n_part_cloud_set(_id_pdm, 1, _nb_part);
// //        cout << cplComm_rank << ": PDM_closest_point n_part cloud set" << endl;

//         for (int i_part = 0 ; i_part < _nb_part ; i_part++)
//             PDM_closest_points_src_cloud_set(_id_pdm, i_part, _n_user_targets[i_part], _coords_user_targets[i_part], _gnum_user_targets[i_part]);
// //        cout << cplComm_rank << ": PDM_closest_point src cloud set" << endl;

//         for (int i_part = 0 ; i_part < _nb_part ; i_part++)
//             PDM_closest_points_tgt_cloud_set(_id_pdm, i_part, _n_target[i_part], _coords_target[i_part], _gnum_target[i_part]);

//         PDM_closest_points_compute(_id_pdm);
// //        cout << cplComm_rank << ": PDM_closest_point computed" << endl;

//         PDM_closest_points_dump_times(_id_pdm);
// //        cout << cplComm_rank << ": PDM_closest_point times dumped" << _id_pdm << endl;

//         for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
//             PDM_closest_points_get(_id_pdm, i_part, &closest_src_gnum, &closest_src_dstance);
//         }
// //        cout << cplComm_rank << ": PDM_closest_point got " << _id_pdm << endl;

//         PDM_closest_points_free (_id_pdm);
// //        cout << cplComm_rank << ": PDM_closest_point freed " << _id_pdm << endl;

      if (!_coupledCodeProperties->localCodeIs() ||
          _localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        for (int i_part = 0; i_part < _nPart; i_part++) {
          if (_src_to_tgt_idx[i_part] != NULL) {
            free(_src_to_tgt_idx[i_part]);
          }
          if (_src_to_tgt[i_part] != NULL) {
            free(_src_to_tgt[i_part]);
          }
          if (_src_to_tgt_dist[i_part] != NULL) {
            free(_src_to_tgt_dist[i_part]);
          }
          _src_to_tgt_idx [i_part] = NULL;
          _src_to_tgt     [i_part] = NULL;
          _src_to_tgt_dist[i_part] = NULL;


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

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        SpatialInterpClosestPoint *cpl_spatial_interp;

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();

          cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }

        else {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();

          cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }

        for (int i_part = 0; i_part < _nPart; i_part++) {

          if (cpl_spatial_interp->_src_to_tgt_idx[i_part] != NULL) {
            free(cpl_spatial_interp->_src_to_tgt_idx[i_part]);
          }
          if (cpl_spatial_interp->_src_to_tgt[i_part] != NULL) {
            free(cpl_spatial_interp->_src_to_tgt[i_part]);
          }
          if (cpl_spatial_interp->_src_to_tgt_dist[i_part] != NULL) {
            free(cpl_spatial_interp->_src_to_tgt_dist[i_part]);
          }
          cpl_spatial_interp->_src_to_tgt_idx [i_part] = NULL;
          cpl_spatial_interp->_src_to_tgt     [i_part] = NULL;
          cpl_spatial_interp->_src_to_tgt_dist[i_part] = NULL;


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

      // localization_init();

      // localization_points_cloud_setting();

      // localization_surface_setting();

      // localization_compute();

      // Reset part_to_part object
      if (_ptsp != nullptr) {
        if (!_coupledCodeProperties->localCodeIs()) {
          PDM_part_to_part_free (_ptsp);
          _ptsp = nullptr;
        }
        else {

          if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
            PDM_part_to_part_free (_ptsp);
            _ptsp = nullptr;

            SpatialInterpClosestPoint *cpl_spatial_interp;

            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

            if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
              std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
              cpl_spatial_interp =
              dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
            }

            else {
              std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
              cpl_spatial_interp =
              dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
            }

            cpl_spatial_interp->_ptsp = NULL;
          }

        }
      }

      // ptp get...

      // localization_get();

      // localization_free();

      // create ptp if null

    }

    void SpatialInterpClosestPoint::interpolate(Field *referenceField, double **buffer) {
      CWP_UNUSED (referenceField);    
      CWP_UNUSED (buffer);    
    }




};
