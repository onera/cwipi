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
#include "pdm_array.h"
#include "pdm_logging.h"

CWP_CLANG_SUPPRESS_WARNING("-Wunused-private-field")


namespace cwipi {
    //CWIPI_CLANG_SUPPRESS_WARNING("-Wunused-private-field")
    SpatialInterpClosestPoint::SpatialInterpClosestPoint() = default;

    SpatialInterpClosestPoint::~SpatialInterpClosestPoint
    (
     )
    {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_closest_src_gnum[i_part] != NULL) {
          free(_closest_src_gnum[i_part]);
        }
        if (_closest_src_dist[i_part] != NULL) {
          free(_closest_src_dist[i_part]);
        }
      }

      delete[] _closest_src_gnum;
      delete[] _closest_src_dist;

      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_tgt_in_src_idx[i_part] != NULL) {
          free(_tgt_in_src_idx[i_part]);
        }
        if (_tgt_in_src_gnum[i_part] != NULL) {
          free(_tgt_in_src_gnum[i_part]);
        }
        if (_tgt_in_src_dist[i_part] != NULL) {
          free(_tgt_in_src_dist[i_part]);
        }
      }

      delete[] _tgt_in_src_idx;
      delete[] _tgt_in_src_gnum;
      delete[] _tgt_in_src_dist;
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
            _tgt_gnum  [i_part] = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
            _tgt_n_gnum[i_part] = _cpl->userTargetNGet (i_part);
          }
        }
      }

      //
      // Target properties
      _closest_src_gnum = new PDM_g_num_t* [_nPart];
      _closest_src_dist = new double*      [_nPart];

      for (int i_part = 0; i_part < _nPart; i_part++) {
        _closest_src_gnum[i_part] = NULL;
        _closest_src_dist[i_part] = NULL;
      }

      //
      // Source properties
      _tgt_in_src_idx  = new int*         [_nPart];
      _tgt_in_src_gnum = new PDM_g_num_t* [_nPart];
      _tgt_in_src_dist = new double*      [_nPart];

      for (int i_part = 0; i_part < _nPart; i_part++) {
        _tgt_in_src_idx [i_part] = NULL;
        _tgt_in_src_gnum[i_part] = NULL;
        _tgt_in_src_dist[i_part] = NULL;
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

//         for (int i_part = 0; i_part < _nb_part ; i_part++)
//             PDM_closest_points_src_cloud_set(_id_pdm, i_part, _n_user_targets[i_part], _coords_user_targets[i_part], _gnum_user_targets[i_part]);
// //        cout << cplComm_rank << ": PDM_closest_point src cloud set" << endl;

//         for (int i_part = 0; i_part < _nb_part ; i_part++)
//             PDM_closest_points_tgt_cloud_set(_id_pdm, i_part, _n_target[i_part], _coords_target[i_part], _gnum_target[i_part]);

//         PDM_closest_points_compute(_id_pdm);
// //        cout << cplComm_rank << ": PDM_closest_point computed" << endl;

//         PDM_closest_points_dump_times(_id_pdm);
// //        cout << cplComm_rank << ": PDM_closest_point times dumped" << _id_pdm << endl;

//         for (int i_part = 0; i_part < _nb_part ; i_part++) {
//             PDM_closest_points_get(_id_pdm, i_part, &closest_src_gnum, &closest_src_dstance);
//         }
// //        cout << cplComm_rank << ": PDM_closest_point got " << _id_pdm << endl;

//         PDM_closest_points_free (_id_pdm);
// //        cout << cplComm_rank << ": PDM_closest_point freed " << _id_pdm << endl;


      if (!_coupledCodeProperties->localCodeIs() ||
          (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet())) {

        for (int i_part = 0; i_part < _nPart; i_part++) {
          if (_closest_src_gnum[i_part] != NULL) {
            free(_closest_src_gnum[i_part]);
          }
          if (_closest_src_dist[i_part] != NULL) {
            free(_closest_src_dist[i_part]);
          }
          _closest_src_gnum[i_part] = NULL;
          _closest_src_dist[i_part] = NULL;


          if (_tgt_in_src_idx[i_part] != NULL) {
            free(_tgt_in_src_idx[i_part]);
          }
          if (_tgt_in_src_gnum[i_part] != NULL) {
            free(_tgt_in_src_gnum[i_part]);
          }
          if (_tgt_in_src_dist[i_part] != NULL) {
            free(_tgt_in_src_dist[i_part]);
          }
          _tgt_in_src_idx [i_part] = NULL;
          _tgt_in_src_gnum[i_part] = NULL;
          _tgt_in_src_dist[i_part] = NULL;


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

          if (cpl_spatial_interp->_closest_src_gnum[i_part] != NULL) {
            free(cpl_spatial_interp->_closest_src_gnum[i_part]);
          }
          if (cpl_spatial_interp->_closest_src_dist[i_part] != NULL) {
            free(cpl_spatial_interp->_closest_src_dist[i_part]);
          }
          cpl_spatial_interp->_closest_src_gnum[i_part] = NULL;
          cpl_spatial_interp->_closest_src_dist[i_part] = NULL;


          if (cpl_spatial_interp->_tgt_in_src_idx[i_part] != NULL) {
            free(cpl_spatial_interp->_tgt_in_src_idx[i_part]);
          }
          if (cpl_spatial_interp->_tgt_in_src_dist[i_part] != NULL) {
            free(cpl_spatial_interp->_tgt_in_src_dist[i_part]);
          }
          if (cpl_spatial_interp->_tgt_in_src_dist[i_part] != NULL) {
            free(cpl_spatial_interp->_tgt_in_src_dist[i_part]);
          }
          cpl_spatial_interp->_tgt_in_src_idx [i_part] = NULL;
          cpl_spatial_interp->_tgt_in_src_gnum[i_part] = NULL;
          cpl_spatial_interp->_tgt_in_src_dist[i_part] = NULL;


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
      int n_closest_pts = 0;
      if (!_coupledCodeProperties->localCodeIs() ||
          (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet())) {

        n_closest_pts = CWP_CLOSEST_POINTS_N_CLOSEST_PTS;
        std::map<std::string, int> prop = _cpl->SpatialInterpPropertiesIntGet();
        std::map<std::string, int>::iterator it;

        it = prop.find("n_closest_pts");
        if (it != prop.end()) {
          n_closest_pts = it->second;
        }

        _id_pdm = PDM_closest_points_create(_pdmCplComm,
                                            n_closest_pts,
                                            PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        // not 100% sure about this...
        int n_part_src = 0;
        int n_part_tgt = 0;
        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          n_part_src = _cplNPart;
        }
        else {
          n_part_src = _nPart;
        }

        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
          n_part_tgt = _cplNPart;
        }
        else {
          n_part_tgt = _nPart;
        }

        PDM_closest_points_n_part_cloud_set(_id_pdm,
                                            n_part_src,
                                            n_part_tgt);


        if (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
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

          cpl_spatial_interp->_id_pdm = _id_pdm;
        }
      }


      // localization_points_cloud_setting();
      // localization_surface_setting();
      if (!_coupledCodeProperties->localCodeIs()) {

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            const double      *src_coord = NULL;
            const PDM_g_num_t *src_g_num = NULL;
            int n_src = 0;

            PDM_closest_points_src_cloud_set(_id_pdm,
                                             i_part,
                                             n_src,
                             (double      *) src_coord,
                             (PDM_g_num_t *) src_g_num);
          }
        }
        else {
          for (int i_part = 0; i_part < _nPart; i_part++) {
            const double      *src_coord = NULL;
            const PDM_g_num_t *src_g_num = NULL;
            int n_src = 0;

            if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              src_g_num = (const PDM_g_num_t *) _mesh->GNumEltsGet  (i_part);
              src_coord =                       _mesh->eltCentersGet(i_part);
              n_src     =                       _mesh->getPartNElts (i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              src_g_num = (const PDM_g_num_t *) _mesh->getVertexGNum  (i_part);
              src_coord =                       _mesh->getVertexCoords(i_part);
              n_src     =                       _mesh->getPartNVertex (i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
              src_g_num = (const PDM_g_num_t *) _cpl->userTargetGNumGet  (i_part);
              src_coord =                       _cpl->userTargetCoordsGet(i_part);
              n_src     =                       _cpl->userTargetNGet     (i_part);
            }

            PDM_closest_points_src_cloud_set(_id_pdm,
                                             i_part,
                                             n_src,
                             (double      *) src_coord,
                             (PDM_g_num_t *) src_g_num);
          }
        }

        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            const double      *tgt_coord = NULL;
            const PDM_g_num_t *tgt_g_num = NULL;
            int n_tgt = 0;


            PDM_closest_points_tgt_cloud_set(_id_pdm,
                                             i_part,
                                             n_tgt,
                             (double      *) tgt_coord,
                             (PDM_g_num_t *) tgt_g_num);
          }
        }
        else {
          for (int i_part = 0; i_part < _nPart; i_part++) {
            const double      *tgt_coord = NULL;
            const PDM_g_num_t *tgt_g_num = NULL;
            int n_tgt = 0;

            if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              tgt_g_num = (const PDM_g_num_t *) _mesh->GNumEltsGet  (i_part);
              tgt_coord =                       _mesh->eltCentersGet(i_part);
              n_tgt     =                       _mesh->getPartNElts (i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              tgt_g_num = (const PDM_g_num_t *) _mesh->getVertexGNum  (i_part);
              tgt_coord =                       _mesh->getVertexCoords(i_part);
              n_tgt     =                       _mesh->getPartNVertex (i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
              tgt_g_num = (const PDM_g_num_t *) _cpl->userTargetGNumGet  (i_part);
              tgt_coord =                       _cpl->userTargetCoordsGet(i_part);
              n_tgt     =                       _cpl->userTargetNGet     (i_part);
            }

            PDM_closest_points_tgt_cloud_set(_id_pdm,
                                             i_part,
                                             n_tgt,
                             (double      *) tgt_coord,
                             (PDM_g_num_t *) tgt_g_num);
          }
        }
      }

      else {
        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());
          cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            for (int i_part = 0; i_part < _cplNPart; i_part++) {
              const double      *src_coord = NULL;
              const PDM_g_num_t *src_g_num = NULL;
              int n_src = 0;

              if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
                src_g_num = (const PDM_g_num_t *) cpl_mesh->GNumEltsGet  (i_part);
                src_coord =                       cpl_mesh->eltCentersGet(i_part);
                n_src     =                       cpl_mesh->getPartNElts (i_part);
              }
              else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
                src_g_num = (const PDM_g_num_t *) cpl_mesh->getVertexGNum  (i_part);
                src_coord =                       cpl_mesh->getVertexCoords(i_part);
                n_src     =                       cpl_mesh->getPartNVertex (i_part);
              }
              else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
                src_g_num = (const PDM_g_num_t *) _cpl->userTargetGNumGet  (i_part);
                src_coord =                       _cpl->userTargetCoordsGet(i_part);
                n_src     =                       _cpl->userTargetNGet     (i_part);
              }

              PDM_closest_points_src_cloud_set(_id_pdm,
                                               i_part,
                                               n_src,
                               (double      *) src_coord,
                               (PDM_g_num_t *) src_g_num);
            }
          }
          else {
            for (int i_part = 0; i_part < _nPart; i_part++) {
              const double      *src_coord = NULL;
              const PDM_g_num_t *src_g_num = NULL;
              int n_src = 0;

              if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
                src_g_num = (const PDM_g_num_t *) _mesh->GNumEltsGet  (i_part);
                src_coord =                       _mesh->eltCentersGet(i_part);
                n_src     =                       _mesh->getPartNElts (i_part);
              }
              else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
                src_g_num = (const PDM_g_num_t *) _mesh->getVertexGNum  (i_part);
                src_coord =                       _mesh->getVertexCoords(i_part);
                n_src     =                       _mesh->getPartNVertex (i_part);
              }
              else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
                src_g_num = (const PDM_g_num_t *) _cpl->userTargetGNumGet  (i_part);
                src_coord =                       _cpl->userTargetCoordsGet(i_part);
                n_src     =                       _cpl->userTargetNGet     (i_part);
              }

              PDM_closest_points_src_cloud_set(_id_pdm,
                                               i_part,
                                               n_src,
                               (double      *) src_coord,
                               (PDM_g_num_t *) src_g_num);
            }
          }

          if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
            for (int i_part = 0; i_part < _cplNPart; i_part++) {
              const double      *tgt_coord = NULL;
              const PDM_g_num_t *tgt_g_num = NULL;
              int n_tgt = 0;

              if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
                tgt_g_num = (const PDM_g_num_t *) cpl_mesh->GNumEltsGet  (i_part);
                tgt_coord =                       cpl_mesh->eltCentersGet(i_part);
                n_tgt     =                       cpl_mesh->getPartNElts (i_part);
              }
              else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
                tgt_g_num = (const PDM_g_num_t *) cpl_mesh->getVertexGNum  (i_part);
                tgt_coord =                       cpl_mesh->getVertexCoords(i_part);
                n_tgt     =                       cpl_mesh->getPartNVertex (i_part);
              }
              else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
                tgt_g_num = (const PDM_g_num_t *) _cpl->userTargetGNumGet  (i_part);
                tgt_coord =                       _cpl->userTargetCoordsGet(i_part);
                n_tgt     =                       _cpl->userTargetNGet     (i_part);
              }

              PDM_closest_points_tgt_cloud_set(_id_pdm,
                                               i_part,
                                               n_tgt,
                               (double      *) tgt_coord,
                               (PDM_g_num_t *) tgt_g_num);
            }
          }
          else {
            for (int i_part = 0; i_part < _nPart; i_part++) {
              const double      *tgt_coord = NULL;
              const PDM_g_num_t *tgt_g_num = NULL;
              int n_tgt = 0;

              if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
                tgt_g_num = (const PDM_g_num_t *) _mesh->GNumEltsGet  (i_part);
                tgt_coord =                       _mesh->eltCentersGet(i_part);
                n_tgt     =                       _mesh->getPartNElts (i_part);
              }
              else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
                tgt_g_num = (const PDM_g_num_t *) _mesh->getVertexGNum  (i_part);
                tgt_coord =                       _mesh->getVertexCoords(i_part);
                n_tgt     =                       _mesh->getPartNVertex (i_part);
              }
              else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
                tgt_g_num = (const PDM_g_num_t *) _cpl->userTargetGNumGet  (i_part);
                tgt_coord =                       _cpl->userTargetCoordsGet(i_part);
                n_tgt     =                       _cpl->userTargetNGet     (i_part);
              }

              PDM_closest_points_tgt_cloud_set(_id_pdm,
                                               i_part,
                                               n_tgt,
                               (double      *) tgt_coord,
                               (PDM_g_num_t *) tgt_g_num);
            }
          }

        }
      }


      // localization_compute();
      if (!_coupledCodeProperties->localCodeIs() ||
          (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet())) {
        PDM_closest_points_compute(_id_pdm);
        PDM_closest_points_dump_times(_id_pdm);
      }


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

      // ptp get
      if (_id_pdm != NULL) {
        PDM_closest_points_part_to_part_get(_id_pdm,
                                            &_ptsp,
                                            PDM_OWNERSHIP_USER);
      }

      // localization_get();
      if (!_coupledCodeProperties->localCodeIs()) {
        for (int i_part = 0; i_part < _nPart; i_part++) {
          if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
            PDM_closest_points_tgt_in_src_get(_id_pdm,
                                              i_part,
                                              &(_tgt_in_src_idx [i_part]),
                                              &(_tgt_in_src_gnum[i_part]));
            PDM_closest_points_tgt_in_src_dist_get(_id_pdm,
                                                   i_part,
                                                   &(_tgt_in_src_idx [i_part]),
                                                   &(_tgt_in_src_dist[i_part]));

            _n_involved_sources_tgt[i_part] = _src_n_gnum[i_part];
            _involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * _n_involved_sources_tgt[i_part]);

            int count = 0;
            for (int i = 0 ; i < _src_n_gnum[i_part] ; ++i) {
              if (_tgt_in_src_idx[i_part][i + 1] > _tgt_in_src_idx[i_part][i]) {
                _involved_sources_tgt[i_part][count] = i + 1;
                ++count;
              }
            }

            _n_involved_sources_tgt[i_part] = count;
            _involved_sources_tgt[i_part] = (int*) realloc(_involved_sources_tgt[i_part], sizeof(int) * count);

          }
          else {
            _tgt_in_src_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
            _tgt_in_src_idx[i_part][0] = 0;

            PDM_closest_points_get(_id_pdm,
                                   i_part,
                                   &(_closest_src_gnum[i_part]),
                                   &(_closest_src_dist[i_part]));

            _n_computed_tgt  [i_part] = PDM_closest_points_n_tgt_get(_id_pdm, i_part);
            _n_uncomputed_tgt[i_part] = 0;

            _computed_tgt[i_part] = (int *) malloc(sizeof(int) * _n_computed_tgt[i_part]);
            for (int i = 0; i < _n_computed_tgt[i_part]; i++) {
              _computed_tgt[i_part][i] = i + 1;
            }
          }

        }
      }
      else {
        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

          SpatialInterpClosestPoint *cpl_spatial_interp;

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            cpl_spatial_interp = dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }
          else {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            cpl_spatial_interp = dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }

          if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
            for (int i_part = 0; i_part < _nPart; i_part++) {
              PDM_closest_points_tgt_in_src_get(_id_pdm,
                                                i_part,
                                                &(_tgt_in_src_idx [i_part]),
                                                &(_tgt_in_src_gnum[i_part]));
              PDM_closest_points_tgt_in_src_dist_get(_id_pdm,
                                                     i_part,
                                                     &(_tgt_in_src_idx [i_part]),
                                                     &(_tgt_in_src_dist[i_part]));

              _n_involved_sources_tgt[i_part] = _src_n_gnum[i_part];
              _involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * _n_involved_sources_tgt[i_part]);

              int count = 0;
              for (int i = 0 ; i < _src_n_gnum[i_part] ; ++i) {
                if (_tgt_in_src_idx[i_part][i + 1] > _tgt_in_src_idx[i_part][i]) {
                  _involved_sources_tgt[i_part][count] = i + 1;
                  ++count;
                }
              }

              _n_involved_sources_tgt[i_part] = count;
              _involved_sources_tgt[i_part] = (int*) realloc(_involved_sources_tgt[i_part], sizeof(int) * count);
            }

            for (int i_part = 0; i_part < _cplNPart; i_part++) {
              cpl_spatial_interp->_tgt_in_src_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
              cpl_spatial_interp->_tgt_in_src_idx[i_part][0] = 0;

              PDM_closest_points_get(_id_pdm,
                                     i_part,
                                     &(cpl_spatial_interp->_closest_src_gnum[i_part]),
                                     &(cpl_spatial_interp->_closest_src_dist[i_part]));

              cpl_spatial_interp->_n_computed_tgt  [i_part] = PDM_closest_points_n_tgt_get(_id_pdm, i_part);
              cpl_spatial_interp->_n_uncomputed_tgt[i_part] = 0;

              cpl_spatial_interp->_computed_tgt[i_part] = (int *) malloc(sizeof(int) * cpl_spatial_interp->_n_computed_tgt[i_part]);
              for (int i = 0; i < cpl_spatial_interp->_n_computed_tgt[i_part]; i++) {
                cpl_spatial_interp->_computed_tgt[i_part][i] = i + 1;
              }
            }
          }
          else {
            for (int i_part = 0; i_part < _nPart; i_part++) {
              _tgt_in_src_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
              _tgt_in_src_idx[i_part][0] = 0;

              PDM_closest_points_get(_id_pdm,
                                     i_part,
                                     &(_closest_src_gnum[i_part]),
                                     &(_closest_src_dist[i_part]));

              _n_computed_tgt  [i_part] = PDM_closest_points_n_tgt_get(_id_pdm, i_part);
              _n_uncomputed_tgt[i_part] = 0;

              _computed_tgt[i_part] = (int *) malloc(sizeof(int) * _n_computed_tgt[i_part]);
              for (int i = 0; i < _n_computed_tgt[i_part]; i++) {
                _computed_tgt[i_part][i] = i + 1;
              }
            }

            for (int i_part = 0; i_part < _cplNPart; i_part++) {
              PDM_closest_points_tgt_in_src_get(_id_pdm,
                                                i_part,
                                                &(cpl_spatial_interp->_tgt_in_src_idx [i_part]),
                                                &(cpl_spatial_interp->_tgt_in_src_gnum[i_part]));
              PDM_closest_points_tgt_in_src_dist_get(_id_pdm,
                                                     i_part,
                                                     &(cpl_spatial_interp->_tgt_in_src_idx [i_part]),
                                                     &(cpl_spatial_interp->_tgt_in_src_dist[i_part]));

              cpl_spatial_interp->_n_involved_sources_tgt[i_part] = cpl_spatial_interp->_src_n_gnum[i_part];
              cpl_spatial_interp->_involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * cpl_spatial_interp->_n_involved_sources_tgt[i_part]);

              int count = 0;
              for (int i = 0 ; i < cpl_spatial_interp->_src_n_gnum[i_part] ; ++i) {
                if (cpl_spatial_interp->_tgt_in_src_idx[i_part][i + 1] > cpl_spatial_interp->_tgt_in_src_idx[i_part][i]) {
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

      // localization_free();
      if (!_coupledCodeProperties->localCodeIs()) {
        PDM_closest_points_free(_id_pdm);
        _id_pdm = nullptr;
      }
      else {
        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
          PDM_closest_points_free(_id_pdm);
          _id_pdm = nullptr;

          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

          SpatialInterpClosestPoint *cpl_spatial_interp;

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            cpl_spatial_interp = dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }
          else {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            cpl_spatial_interp = dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }
          cpl_spatial_interp->_id_pdm = nullptr;
        }
      }

      // create ptp if null
      if (_ptsp == NULL) {

        if (!_coupledCodeProperties->localCodeIs()) {
          _ptsp = PDM_part_to_part_create ((const PDM_g_num_t **) _src_gnum,
                                           (const int          *) _src_n_gnum,
                                           _nPart,
                                           (const PDM_g_num_t **) _tgt_gnum,
                                           (const int          *) _tgt_n_gnum,
                                           _nPart,
                                           (const int         **) _tgt_in_src_idx,
                                           (const PDM_g_num_t **) _tgt_in_src_gnum,
                                           _pdmCplComm);
        }
        else {
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

            if (_ptsp == NULL) {
              if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
                _ptsp = PDM_part_to_part_create ((const PDM_g_num_t **) _src_gnum,
                                                 (const int          *) _src_n_gnum,
                                                 _nPart,
                                                 (const PDM_g_num_t **) cpl_spatial_interp->_tgt_gnum,
                                                 (const int          *) cpl_spatial_interp->_tgt_n_gnum,
                                                 _cplNPart,
                                                 (const int         **) _tgt_in_src_idx,
                                                 (const PDM_g_num_t **) _tgt_in_src_gnum,
                                                 _pdmCplComm);
              }
              else {
                _ptsp = PDM_part_to_part_create ((const PDM_g_num_t **) cpl_spatial_interp->_src_gnum,
                                                 (const int          *) cpl_spatial_interp->_src_n_gnum,
                                                 _cplNPart,
                                                 (const PDM_g_num_t **) _tgt_gnum,
                                                 (const int          *) _tgt_n_gnum,
                                                 _nPart,
                                                 (const int         **) cpl_spatial_interp->_tgt_in_src_idx,
                                                 (const PDM_g_num_t **) cpl_spatial_interp->_tgt_in_src_gnum,
                                                 _pdmCplComm);
              }
            }

            cpl_spatial_interp->_ptsp = _ptsp;

          }
        }
      }

    }



    void SpatialInterpClosestPoint::interpolate(Field *referenceField, double **buffer) {

      int nComponent = referenceField->nComponentGet();
      CWP_Dof_location_t referenceFieldType = referenceField->locationGet();
      CWP_Interpolation_t interpolationType = referenceField->interpolationTypeGet();
      const CWP_Field_storage_t storage = referenceField->storageTypeGet();

      if (interpolationType == CWP_INTERPOLATION_USER) {
        PDM_error(__FILE__, __LINE__, 0, "user interpolation not implemented yet");
      }

      else {

        int n_closest_pts = CWP_CLOSEST_POINTS_N_CLOSEST_PTS;
        std::map<std::string, int> prop = _cpl->SpatialInterpPropertiesIntGet();
        std::map<std::string, int>::iterator it;

        it = prop.find("n_closest_pts");
        if (it != prop.end()) {
          n_closest_pts = it->second;
        }

        for (int i_part = 0; i_part < _nPart; i_part++) {

          double *referenceData = (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_TARGET);

          int n_pts = 0;
          if (referenceFieldType == CWP_DOF_LOCATION_CELL_CENTER) {
            n_pts = _mesh->getPartNElts(i_part);
          }
          else if (referenceFieldType == CWP_DOF_LOCATION_NODE) {
            n_pts = _mesh->getPartNVertex(i_part);
          }
          else {
            PDM_error(__FILE__, __LINE__, 0, "user tgt not supported yet");
          }


          // TO DO : least square interpolation

          if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {

            for (int i = 0; i < n_pts; i++) {
              double sum_w = 0;
              for (int k = 0; k < n_closest_pts; k++) {
                double w = 1./std::max(1e-16, sqrt(_closest_src_dist[i_part][n_closest_pts*i + k]));
                sum_w += w;
                for (int j = 0; j < nComponent; j++) {
                  referenceData[n_pts*j + i] += w*buffer[i_part][n_closest_pts*(n_pts*j + i) + k];
                }
              }

              for (int j = 0; j < nComponent; j++) {
                referenceData[n_pts*j + i] /= sum_w;
              }

              if (0) {
                log_trace("received:\n");
                for (int k = 0; k < n_closest_pts; k++) {
                  log_trace("from " PDM_FMT_G_NUM ", at dist %f : ",
                            _closest_src_gnum[i_part][n_closest_pts*i + k],
                            sqrt(_closest_src_dist[i_part][n_closest_pts*i + k]));
                  for (int j = 0; j < nComponent; j++) {
                    log_trace("%f ", buffer[i_part][n_closest_pts*(n_pts*j + i) + k]);
                  }
                  log_trace("\n");
                }
                for (int j = 0; j < nComponent; j++) {
                  log_trace("%f ", referenceData[n_pts*j + i]);
                }
                log_trace("\n");
              }
            }
          }

          else { // if (storage == CWP_FIELD_STORAGE_INTERLACED) {

            if (0) {
              log_trace(">>>\nreceived :\n");
              for (int i = 0; i < n_pts; i++) {
                for (int k = 0; k < n_closest_pts; k++) {
                  log_trace("from " PDM_FMT_G_NUM " : ", _closest_src_gnum[i_part][n_closest_pts*i + k]);
                  for (int j = 0; j < nComponent; j++) {
                    log_trace("%f ", buffer[i_part][nComponent*(n_closest_pts*i + k) + j]);
                  }
                  log_trace("\n");
                }
              }
              log_trace("<<<\n");
            }

            for (int i = 0; i < n_pts; i++) {
              for (int j = 0; j < nComponent; j++) {
                referenceData[nComponent*i + j] = 0;
              }

              double sum_w = 0;
              for (int k = 0; k < n_closest_pts; k++) {
                double w = 1./std::max(1e-16, sqrt(_closest_src_dist[i_part][n_closest_pts*i + k]));
                sum_w += w;
                for (int j = 0; j < nComponent; j++) {
                  referenceData[nComponent*i + j] += w*buffer[i_part][nComponent*(n_closest_pts*i + k) + j];
                }
              }

              for (int j = 0; j < nComponent; j++) {
                referenceData[nComponent*i + j] /= sum_w;
              }

              if (0) {
                log_trace("received:\n");
                for (int k = 0; k < n_closest_pts; k++) {
                  log_trace("from " PDM_FMT_G_NUM ", at dist %f : ",
                            _closest_src_gnum[i_part][n_closest_pts*i + k],
                            sqrt(_closest_src_dist[i_part][n_closest_pts*i + k]));
                  PDM_log_trace_array_double(&buffer[i_part][nComponent*(n_closest_pts*i + k)],
                                             nComponent,
                                             "");
                }
                PDM_log_trace_array_double(&referenceData[nComponent*i],
                                           nComponent,
                                           "interpolated : ");
              }
            }
          }

        }


      }


    }




};
