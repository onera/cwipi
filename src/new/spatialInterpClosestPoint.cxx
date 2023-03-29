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
#include "pdm_linear_algebra.h"

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
        // if (_weights[i_part] != NULL) {
        //   free(_weights[i_part]);
        // }
      }

      free ( _closest_src_gnum);
      // free ( _weights);

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

      free ( _tgt_in_src_idx);
      free ( _tgt_in_src_gnum);
      free ( _tgt_in_src_dist);

      if (_send_coord != NULL) {
        free(_send_coord);
        _send_coord = NULL;
      }

      if (_recv_coord != NULL) {
        int n_part_tgt = 0;
        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          n_part_tgt = _nPart;
        }
        else {
          n_part_tgt  = _cplNPart;
        }

        for (int i_part = 0; i_part < n_part_tgt; i_part++) {
          if (_recv_coord[i_part] != NULL) {
            free(_recv_coord[i_part]);
          }
        }
        free(_recv_coord);
        _recv_coord = NULL;
      }
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

      _coordinates_exchanged = 0;

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
      _closest_src_gnum = (PDM_g_num_t**) malloc (sizeof(PDM_g_num_t*) * _nPart);

      for (int i_part = 0; i_part < _nPart; i_part++) {
        _closest_src_gnum[i_part] = NULL;
        // _closest_src_dist[i_part] = NULL;
      }

      //
      // Source properties
      _tgt_in_src_idx  =  (int **) malloc (sizeof(int *) * (_nPart));
      _tgt_in_src_gnum = (PDM_g_num_t**) malloc (sizeof(PDM_g_num_t*) * _nPart);
      
      _tgt_in_src_dist =  (double **) malloc (sizeof(double *) * (_nPart));

      for (int i_part = 0; i_part < _nPart; i_part++) {
        _tgt_in_src_idx [i_part] = NULL;
        _tgt_in_src_gnum[i_part] = NULL;
        _tgt_in_src_dist[i_part] = NULL;
      }

      _send_coord = NULL;
      _recv_coord = NULL;
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
      _coordinates_exchanged = 0;

      if (!_coupledCodeProperties->localCodeIs() ||
          (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet())) {

        for (int i_part = 0; i_part < _nPart; i_part++) {
          if (_closest_src_gnum[i_part] != NULL) {
            free(_closest_src_gnum[i_part]);
          }
          if (_weights[i_part] != NULL) {
            free(_weights[i_part]);
          }
          _closest_src_gnum[i_part] = NULL;
          _weights[i_part] = NULL;


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


        CWP_Dynamic_mesh_t dyn_mesh = _cpl->DisplacementGet();

        if (dyn_mesh != CWP_DYNAMIC_MESH_STATIC) {
          if (_send_coord != NULL) {
            free(_send_coord);
            _send_coord = NULL;
          }
        }

        if (_recv_coord != NULL) {
          for (int i_part = 0; i_part < _nPart; i_part++) {
            if (_recv_coord[i_part] != NULL) {
              free(_recv_coord[i_part]);
            }
          }
          free(_recv_coord);
          _recv_coord = NULL;
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

        for (int i_part = 0; i_part < _cplNPart; i_part++) {

          if (cpl_spatial_interp->_closest_src_gnum[i_part] != NULL) {
            free(cpl_spatial_interp->_closest_src_gnum[i_part]);
          }
          if (cpl_spatial_interp->_weights[i_part] != NULL) {
            free(cpl_spatial_interp->_weights[i_part]);
          }
          cpl_spatial_interp->_closest_src_gnum[i_part] = NULL;
          cpl_spatial_interp->_weights[i_part] = NULL;


          if (cpl_spatial_interp->_tgt_in_src_idx[i_part] != NULL) {
            free(cpl_spatial_interp->_tgt_in_src_idx[i_part]);
          }
          if (cpl_spatial_interp->_tgt_in_src_gnum[i_part] != NULL) {
            free(cpl_spatial_interp->_tgt_in_src_gnum[i_part]);
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

        if (cpl_spatial_interp->_send_coord != NULL) {
          free(cpl_spatial_interp->_send_coord);
          cpl_spatial_interp->_send_coord = NULL;
        }

        if (cpl_spatial_interp->_recv_coord != NULL) {
          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            if (cpl_spatial_interp->_recv_coord[i_part] != NULL) {
              free(cpl_spatial_interp->_recv_coord[i_part]);
            }
          }
          free(cpl_spatial_interp->_recv_coord);
          cpl_spatial_interp->_recv_coord = NULL;
        }
      }

      /* Set source and target point clouds */
      int n_closest_pts = CWP_CLOSEST_POINTS_N_CLOSEST_PTS;
      std::map<std::string, int> prop = _cpl->SpatialInterpPropertiesIntGet();
      std::map<std::string, int>::iterator it;

      it = prop.find("n_closest_pts");
      if (it != prop.end()) {
        n_closest_pts = it->second;
      }


      if (!_coupledCodeProperties->localCodeIs()) {
        _id_pdm = PDM_closest_points_create(_pdmUnionComm,
                                            n_closest_pts,
                                            PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        int n_part_src = 0;
        int n_part_tgt = 0;
        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
          n_part_src = _nPart;
        }
        else {
          n_part_src  = _cplNPart;
        }
        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          n_part_tgt = _nPart;
        }
        else {
          n_part_tgt  = _cplNPart;
        }

        PDM_closest_points_n_part_cloud_set(_id_pdm,
                                            n_part_src,
                                            n_part_tgt);

        _send_coord = (const double **) malloc(sizeof(double *) * n_part_src);

        // source point cloud
        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
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

            _send_coord[i_part] = src_coord;
          }
        }
        else {
          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            PDM_closest_points_src_cloud_set(_id_pdm, i_part, 0, NULL, NULL);
            _send_coord[i_part] = NULL;
          }
        }

        // target point cloud
        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          for (int i_part = 0 ; i_part < _nPart; i_part++) {
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
        else {
          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            PDM_closest_points_tgt_cloud_set(_id_pdm, i_part, 0, NULL, NULL);
          }
        }

      }

      else {
        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
          _id_pdm = PDM_closest_points_create(_pdmUnionComm,
                                              n_closest_pts,
                                              PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

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

          int n_part_src = 0;
          int n_part_tgt = 0;
          if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
            n_part_src = _nPart;
          }
          else {
            n_part_src  = _cplNPart;
          }
          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            n_part_tgt = _nPart;
          }
          else {
            n_part_tgt  = _cplNPart;
          }

          PDM_closest_points_n_part_cloud_set(_id_pdm,
                                              n_part_src,
                                              n_part_tgt);

          _send_coord = (const double **) malloc(sizeof(double *) * _nPart);//n_part_src);
          cpl_spatial_interp->_send_coord = (const double **) malloc(sizeof(double *) * _cplNPart);

          // source point cloud
          if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
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

              _send_coord[i_part] = src_coord;
            }

            for (int i_part = 0; i_part < _cplNPart; i_part++) {
              cpl_spatial_interp->_send_coord[i_part] = NULL;
            }
          }
          else {

            for (int i_part = 0; i_part < _nPart; i_part++) {
              _send_coord[i_part] = NULL;
            }

            cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();

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

              cpl_spatial_interp->_send_coord[i_part] = src_coord;
            }
          }

          // target point cloud
          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
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
          else {

            cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();

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
        }
      }


      /* Compute */
      if (!_coupledCodeProperties->localCodeIs() ||
          (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet())) {
        PDM_closest_points_compute(_id_pdm);
        PDM_closest_points_dump_times(_id_pdm);
      }


      /* Reset part_to_part object */
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

      /* Get PDM part_to_part object */
      if (_id_pdm != NULL) {
        PDM_closest_points_part_to_part_get(_id_pdm,
                                            &_ptsp,
                                            PDM_OWNERSHIP_USER);

        // if (0) {
        //   int  *n_ref_gnum2;
        //   int **ref_gnum2;
        //   PDM_part_to_part_ref_lnum2_get (_ptsp,
        //                                   &n_ref_gnum2,
        //                                   &ref_gnum2);

        //   int          **gnum1_come_from_idx;
        //   PDM_g_num_t  **gnum1_come_from;
        //   PDM_part_to_part_gnum1_come_from_get (_ptsp,
        //                                         &gnum1_come_from_idx,
        //                                         &gnum1_come_from);

        //   int n_part1;
        //   int n_part2;
        //   PDM_part_to_part_n_part_get(_ptsp, &n_part1, &n_part2);

        //   for (int i = 0; i < n_part2; i++) {
        //     log_trace("n_ref_gnum2[%d] = %d, end idx = %d\n",
        //               i, n_ref_gnum2[i], gnum1_come_from_idx[i][n_ref_gnum2[i]]);
        //     for (int j = 0; j < n_ref_gnum2[i]; j++) {
        //       log_trace("  %d -> %d\n", j, gnum1_come_from_idx[i][j+1] - gnum1_come_from_idx[i][j]);
        //     }
        //   }
        // }
      }

      /* Get PDM results */
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
                                   &(_weights[i_part]));

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
                                     &(cpl_spatial_interp->_weights[i_part]));

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
                                     &(_weights[i_part]));

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

      /* Free PDM object */
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

      /* Create part_to_part object if null */

      if (!_coupledCodeProperties->localCodeIs()) {
        if (_ptsp == NULL) {
          _ptsp = PDM_part_to_part_create ((const PDM_g_num_t **) _src_gnum,
                                           (const int          *) _src_n_gnum,
                                           _nPart,
                                           (const PDM_g_num_t **) _tgt_gnum,
                                           (const int          *) _tgt_n_gnum,
                                           _nPart,
                                           (const int         **) _tgt_in_src_idx,
                                           (const PDM_g_num_t **) _tgt_in_src_gnum,
                                           _pdmUnionComm);
        }
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
                                               _pdmUnionComm);
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
                                               _pdmUnionComm);
            }
          }

          cpl_spatial_interp->_ptsp = _ptsp;

        }
      }
    }




    void SpatialInterpClosestPoint::issend(Field *referenceField) {

      if (!_coordinates_exchanged) {
        /* Send source points coordinates */
        if (!_coupledCodeProperties->localCodeIs()) {

          PDM_part_to_part_iexch(_ptsp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                                 1,
                                 3*sizeof(double),
                                 NULL,
                (const void  **) _send_coord,
                                 NULL,
                (      void ***) &_recv_coord,
                                 &(_send_coord_request));
          _recv_coord_request = _send_coord_request;
        }

        else {
          if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            SpatialInterpClosestPoint *cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

            PDM_part_to_part_iexch(_ptsp,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                                   1,
                                   3*sizeof(double),
                                   NULL,
                  (const void  **) _send_coord,
                                   NULL,
                  (      void ***) &cpl_spatial_interp->_recv_coord,
                                   &(_send_coord_request));
            cpl_spatial_interp->_recv_coord_request = _send_coord_request;
          }
        }
      }

      SpatialInterp::issend(referenceField);
    }


    void SpatialInterpClosestPoint::waitIssend(Field *referenceField) {

      if (!_coordinates_exchanged) {

        if (!_coupledCodeProperties->localCodeIs()) {
          PDM_part_to_part_iexch_wait(_ptsp, _send_coord_request);
        }

        else {
          if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            SpatialInterpClosestPoint *cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

            PDM_part_to_part_iexch_wait(_ptsp, cpl_spatial_interp->_recv_coord_request);
          }
        }

        _coordinates_exchanged = 1;
      }

      SpatialInterp::waitIssend(referenceField);
    }


    void SpatialInterpClosestPoint::irecv(Field *referenceField) {

      /* Receive source points coordinates */
      if (!_coordinates_exchanged) {

        if (!_coupledCodeProperties->localCodeIs()) {

          PDM_part_to_part_iexch(_ptsp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                                 1,
                                 3*sizeof(double),
                                 NULL,
                (const void  **) _send_coord,
                                 NULL,
                (      void ***) &_recv_coord,
                                 &(_send_coord_request));
          _recv_coord_request = _send_coord_request;
        }

        else {
          if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            SpatialInterpClosestPoint *cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);


            PDM_part_to_part_iexch(_ptsp,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                                   1,
                                   3*sizeof(double),
                                   NULL,
                  (const void  **) cpl_spatial_interp->_send_coord,
                                   NULL,
                  (      void ***) &_recv_coord,
                                   &(cpl_spatial_interp->_send_coord_request));
            _recv_coord_request = cpl_spatial_interp->_send_coord_request;
          }
        }

      }

      SpatialInterp::irecv(referenceField);
    }


    void SpatialInterpClosestPoint::waitIrecv(Field *referenceField) {

      if (!_coordinates_exchanged) {

        if (!_coupledCodeProperties->localCodeIs()) {
          PDM_part_to_part_iexch_wait(_ptsp, _send_coord_request);
        }

        else {
          if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            SpatialInterpClosestPoint *cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

            PDM_part_to_part_iexch_wait(_ptsp, cpl_spatial_interp->_send_coord_request);
          }
        }

        _coordinates_exchanged = 1;
      }

      SpatialInterp::waitIrecv(referenceField);
    }





    static void _interp_idw
    (
      const int     n_closest_pts,
      const int     stride,
      const double *src_value,
      const double *src_dist2,
            double *tgt_value
     )
    {
      const double eps_dist2 = 1e-30;
      for (int j = 0; j < stride; j++) {
        tgt_value[j] = 0.;
      }

      double sum_w = 0.;
      for (int i = 0; i < n_closest_pts; i++) {
        double w = 1./std::max(eps_dist2, src_dist2[i]);
        sum_w += w;
        for (int j = 0; j < stride; j++) {
          tgt_value[j] += w*src_value[stride*i + j];
        }
      }

      sum_w = 1./sum_w;
      for (int j = 0; j < stride; j++) {
        tgt_value[j] *= sum_w;
      }
    }


    /**
     * Solve the linear system Ax = b using Gaussian elimination,
     * where A is a n*n matrix and b, x are n*stride matrices
     * (Aij = A[n*i+j], bij = b[stride*i+j], xij = x[stride*i+j])
     *
     * return 1 if A is singular, 0 else
     */
    static int _linsolve
    (
      const int     n,
      const int     stride,
            double *A,
            double *b,
            double *x
     )
    {
      const double eps = 1e-15;

      if (0) {
        log_trace("A = \n");
        for (int ii = 0; ii < n; ii++) {
          for (int jj = 0; jj < n; jj++) {
            log_trace("%f ", A[n*ii+jj]);
          }
          log_trace("\n");
        }

        // log_trace("b = \n");
        // for (int i = 0; i < n; i++) {
        //   for (int j = 0; j < stride; j++) {
        //     log_trace("%f ", b[stride*i+j]);
        //   }
        //   log_trace("\n");
        // }
      }

      for (int i = 0; i < n; i++) {
        /* Seek best pivot */
        double amax = std::fabs(A[n*i+i]);
        int imax = i;
        for (int k = i+1; k < n; k++) {
          double aki = std::fabs(A[n*k+i]);
          if (aki > amax) {
            amax = aki;
            imax = k;
          }
        }

        if (amax <= eps) {
          /* matrix A is singular */
          // log_trace("A is singular\n");
          // for (int ii = 0; ii < n; ii++) {
          //   for (int jj = 0; jj < n; jj++) {
          //     log_trace("%f ", A[n*ii+jj]);
          //   }
          //   log_trace("\n");
          // }
          return 1;
        }

        /* Swap rows i and imax */
        if (i != imax) {
          for (int j = 0; j < n; j++) {
            double tmp = A[n*i+j];
            A[n*i   +j] = A[n*imax+j];
            A[n*imax+j] = tmp;
          }

          for (int j = 0; j < stride; j++) {
            double tmp = b[stride*i + j];
            b[stride*i    + j] = b[stride*imax + j];
            b[stride*imax + j] = tmp;
          }
        }

        /* Eliminate subdiagonal terms */
        double inv_amax = 1./A[n*i+i];

        for (int k = i+1; k < n; k++) {
          double r = A[n*k+i] * inv_amax;
          for (int j = i+1; j < n; j++) {
            A[n*k+j] -= r * A[n*i+j];
          }
          A[n*k+i] = 0.;

          for (int j = 0; j < stride; j++) {
            b[stride*k + j] -= r * b[stride*i + j];
          }
        }
      }

      /* Solve triangular system */
      memcpy(x, b, sizeof(double) * n * stride);

      for (int i = n-1; i >= 0; i--) {
        for (int j = i+1; j < n; j++) {
          for (int k = 0; k < stride; k++) {
            x[stride*i + k] -= x[stride*j + k] * A[n*i+j];
          }
        }

        double inv_ai = 1./A[n*i+i];
        for (int k = 0; k < stride; k++) {
          x[stride*i + k] *= inv_ai;
        }
      }

      if (0) {
        log_trace("x = \n");
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < stride; j++) {
            log_trace("%f ", x[stride*i+j]);
          }
          log_trace("\n");
        }
      }

      return 0;
    }


    static void _interp_least_squares
    (
      const int     n_closest_pts,
      const int     stride,
      const double *src_value,
      const double *src_coord,
      const double *src_dist2,
      const double *tgt_coord,
            double *tgt_value
     )
    {
      double A[4*4] = {0.};
      double b[4*stride] = {0.};

      for (int i = 0; i < n_closest_pts; i++) {
        // log_trace("i = %d / %d\n", i, n_closest_pts);
        double x = src_coord[3*i    ];
        double y = src_coord[3*i + 1];
        double z = src_coord[3*i + 2];

        A[4*0+0] += x * x;
        A[4*0+1] += x * y;
        A[4*0+2] += x * z;
        A[4*0+3] += x;

        A[4*1+1] += y * y;
        A[4*1+2] += y * z;
        A[4*1+3] += y;

        A[4*2+2] += z * z;
        A[4*2+3] += z;

        A[4*3+3] += 1.;

        for (int j = 0; j < stride; j++) {
          double f = src_value[stride*i + j];
          b[stride*0 + j] += x * f;
          b[stride*1 + j] += y * f;
          b[stride*2 + j] += z * f;
          b[stride*3 + j] += f;
        }
      }

      /* Symmetrize */
      A[4*1+0] = A[4*0+1];
      A[4*2+0] = A[4*0+2];
      A[4*3+0] = A[4*0+3];

      A[4*2+1] = A[4*1+2];
      A[4*3+1] = A[4*1+3];

      A[4*3+2] = A[4*2+3];

      double coeff[4*stride];
      if (_linsolve(4, stride, A, b, coeff) == 0) {
        for (int j = 0; j < stride; j++) {
          tgt_value[j] =
          coeff[stride*0 + j] * tgt_coord[0] +
          coeff[stride*1 + j] * tgt_coord[1] +
          coeff[stride*2 + j] * tgt_coord[2] +
          coeff[stride*3 + j];
        }
      }
      else {
        // what do we do if A is singular? SVD? IDW?
        _interp_idw(n_closest_pts,
                    stride,
                    src_value,
                    src_dist2,
                    tgt_value);
        // for (int j = 0; j < stride; j++) {
        //   tgt_value[j] = (j+1)*1000;
        // }
      }
    }


    static inline void _basis_vector
    (
     const int     dim,
     const int     degree,
     const double *x,
     const double *x0,
           double *b
     )
    {
      b[0] = 1.;

      for (int j = 0; j < degree; j++) {
        for (int i = 0; i < dim; i++) {
          double _x = x[i];
          if (x0 != NULL) {
            _x -= x0[i];
          }
          b[1+dim*j+i] = 1.;
          for (int k = 0; k <= j; k++) {
            b[1+dim*j+i] *= _x;
          }
        }
      }
    }


    static inline double _weight_function
    (
     const double d2
     )
    {
      const double eps2 = 1e-12;

      return 1. / (d2 + eps2);
    }

    static void _interp_weighted_least_squares
    (
      const int     degree,
      const int     n_closest_pts,
      const int     stride,
      const double *src_value,
      const double *src_coord,
      const double *src_dist2,
      const double *tgt_coord,
            double *tgt_value
     )
    {
      #define siz (1 + degree*3)

      double A[siz*siz] = {0.};
      double rhs[siz*stride] = {0.};
      double b[siz];

      for (int i = 0; i < n_closest_pts; i++) {
        double wi = _weight_function(src_dist2[i]);

        _basis_vector(3,
                      degree,
                      src_coord + 3*i,
                      tgt_coord,
                      b);

        for (int j = 0; j < siz; j++) {
          for (int k = 0; k < stride; k++) {
            rhs[stride*j+k] += wi * b[j] * src_value[stride*i+k];
          }

          for (int k = j; k < siz; k++) {
            A[siz*j+k] += wi * b[j] * b[k];
          }
          for (int k = 0; k < j; k++) {
            A[siz*j+k] = A[siz*k+j];
          }
        }
      }

      double c[siz*stride];
      int stat = PDM_linear_algebra_linsolve_svd(siz,
                                                 siz,
                                                 stride,
                                                 0.,
                                      (double *) A,
                                      (double *) rhs,
                                                 c);
      if (stat == 0) {

        for (int j = 0; j < stride; j++) {
          tgt_value[j] = c[j];
        }

      }
      else {
        // log_trace("singular matrix! ");
        // PDM_log_trace_array_double(tgt_coord, 3, "tgt_coord : ");
        // what do we do if A is singular? SVD? IDW?
        _interp_idw(n_closest_pts,
                    stride,
                    src_value,
                    src_dist2,
                    tgt_value);
      }

      #undef siz
    }




    void SpatialInterpClosestPoint::interpolate(Field *referenceField, double **buffer) {

      int nComponent = referenceField->nComponentGet();
      CWP_Dof_location_t referenceFieldType = referenceField->locationGet();
      CWP_Interpolation_t interpolationType = referenceField->interpolationTypeGet();
      const CWP_Field_storage_t storage = referenceField->storageTypeGet();

      if (interpolationType == CWP_INTERPOLATION_USER) {
        CWP_Interp_function_t interpolationFunction = referenceField->interpolationFunctionGet();

        if (interpolationFunction != NULL) {

          for (int i_part = 0 ; i_part < _nPart ; i_part++) {

            (*interpolationFunction) (_localCodeProperties->nameGet().c_str(),
                                      _cpl->IdGet().c_str(),
                                      referenceField->fieldIDGet().c_str(),
                                      i_part,
                                      _cpl->spatialInterpAlgoGet(),
                                      storage,
                           (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_SOURCE),
                           (double *) buffer[i_part]);
          }

        }

      }

      else {

        int n_closest_pts = CWP_CLOSEST_POINTS_N_CLOSEST_PTS;
        std::map<std::string, int> prop = _cpl->SpatialInterpPropertiesIntGet();
        std::map<std::string, int>::iterator it;

        it = prop.find("n_closest_pts");
        if (it != prop.end()) {
          n_closest_pts = it->second;
        }

        int use_idw_interpolation = 0;
        double       *src_coord = NULL;
        const double *tgt_coord = NULL;

        double *src_value = NULL;
        double *tgt_value = NULL;

        if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
          src_value = (double *) malloc(sizeof(double) * nComponent * n_closest_pts);
          tgt_value = (double *) malloc(sizeof(double) * nComponent);
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
            n_pts = _cpl->userTargetNGet(i_part);
            // PDM_error(__FILE__, __LINE__, 0, "user tgt not supported yet");
          }

          if (!use_idw_interpolation) {
            if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              tgt_coord = _mesh->eltCentersGet(i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              tgt_coord = _mesh->getVertexCoords(i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
              tgt_coord = _cpl->userTargetCoordsGet(i_part);
            }
          }

          src_coord = _recv_coord[i_part];

          for (int i = 0; i < n_pts; i++) {
            // log_trace("tgt pt %d / %d\n", i, n_pts);

            if (0) {
              for (int k = 0; k < n_closest_pts; k++) {
                log_trace("src point " PDM_FMT_G_NUM " : %f %f %f\n",
                          _closest_src_gnum[i_part][n_closest_pts*i + k],
                          src_coord[3*(n_closest_pts*i+k)  ],
                          src_coord[3*(n_closest_pts*i+k)+1],
                          src_coord[3*(n_closest_pts*i+k)+2]);
              }
            }

            if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
              // interlace src_value
              for (int k = 0; k < n_closest_pts; k++) {
                for (int j = 0; j < nComponent; j++) {
                  src_value[nComponent*k + j] = buffer[i_part][n_closest_pts*(n_pts*j + i) + k];
                }
              }
            }
            else {
              src_value = buffer[i_part] + nComponent*n_closest_pts*i;
              tgt_value = referenceData + nComponent*i;
            }

            if (use_idw_interpolation) {
              _interp_idw(n_closest_pts,
                          nComponent,
                          src_value,
                          _weights[i_part] + n_closest_pts*i,
                          tgt_value);
            }
            else {
              char *env_var = NULL;
              int use_wls = 0;
              env_var = getenv("USE_WLS");
              if (env_var != NULL) {
                use_wls = atoi(env_var);
              }
              if (use_wls) {
                int degree = 1;
                env_var = getenv("WLS_DEGREE");
                if (env_var != NULL) {
                  degree = atoi(env_var);
                }
                _interp_weighted_least_squares(degree,
                                               n_closest_pts,
                                               nComponent,
                                               src_value,
                                               src_coord + 3*n_closest_pts*i,
                                               _weights[i_part] + n_closest_pts*i,
                                               tgt_coord + 3*i,
                                               tgt_value);
              }
              else {
                _interp_least_squares(n_closest_pts,
                                      nComponent,
                                      src_value,
                                      src_coord + 3*n_closest_pts*i,
                                      _weights[i_part] + n_closest_pts*i,
                                      tgt_coord + 3*i,
                                      tgt_value);
              }


              if (0) {
                log_trace("recv :\n");
                for (int k = 0; k < n_closest_pts; k++) {
                  log_trace("  from " PDM_FMT_G_NUM " : ", _closest_src_gnum[i_part][n_closest_pts*i + k]);
                  PDM_log_trace_array_double(src_value + nComponent*k,
                                             nComponent,
                                             "");
                }
                PDM_log_trace_array_double(tgt_value,
                                           nComponent,
                                           "interpolated : ");
                PDM_log_trace_array_double(tgt_coord + 3*i,
                                           3,
                                           "tgt_coord    : ");
              }
            }

            if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
              // de-interlace tgt_value
              for (int j = 0; j < nComponent; j++) {
                referenceData[n_pts*j + i] = tgt_value[j];
              }
            }

          }
        }

        if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
          free(src_value);
          free(tgt_value);
        }


      }


    }


    double **
    SpatialInterpClosestPoint::closest_src_coord_get()
    {
      return _recv_coord;
    }

};
