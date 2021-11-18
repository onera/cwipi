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
#include <pdm_gnum_location.h>

#include <spatialInterpLocation.hxx>

/**
 * \cond
 */

namespace cwipi {
  SpatialInterpLocation::SpatialInterpLocation() = default;

  SpatialInterpLocation::~SpatialInterpLocation
  (
  )
  {

    delete[] _tgt_distance;
    delete[] _tgt_projected;
    delete[] _tgt_closest_elt_gnum;

    delete[] _elt_pts_inside_idx;
    delete[] _points_gnum;
    delete[] _points_coords;
    delete[] _points_uvw;
    delete[] _points_dist2;
    delete[] _points_projected_coords;

    printf("delete SpatialInterpLocation\n");

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

    //
    // Data for PDM_part1_to_selected_part2_t

    if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
      for (int i_part = 0 ; i_part < _nPart ; i_part++) { 
       _src_gnum[i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
       printf("_mesh->getPartNElts (i_part) : %d\n", _mesh->getPartNElts (i_part));
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
      printf("_src_gnum %d :", _src_n_gnum[0]);
      for (int i = 0; i < _src_n_gnum[0]; i++) {
        printf(" %ld", _src_gnum[0][i]);
      }
      printf("\n");

      printf("_tgt_gnum %d :", _tgt_n_gnum[0]);
      for (int i = 0; i < _tgt_n_gnum[0]; i++) {
        printf(" %ld", _tgt_gnum[0][i]);
      }
      printf("\n");

      printf("_points_gnum %d :", _src_n_gnum[0]);
      for (int i = 0; i < _src_n_gnum[0]; i++) {
        for (int j = _elt_pts_inside_idx[0][i]; j < _elt_pts_inside_idx[0][i+1]; j++) {
          printf(" %ld", _points_gnum[0][j]);
        }
        printf("\n");
      }
      printf("\n");
      fflush(stdout);
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

  void *SpatialInterpLocation::interpolate (Field *referenceField, double **buffer) 
  {
    int nComponent = referenceField->nComponentGet();
    CWP_Dof_location_t referenceFieldType = referenceField->locationGet();
    int dataTypeSize = referenceField->dataTypeSizeGet();
    CWP_Interpolation_t interpolationType = referenceField->interpolationTypeGet();

    


    // void *interpolatedData = referenceField->sendBufferGet();

    // // Allocate interpolatedData and weights
    // if (interpolatedData != NULL) {
    //   free(interpolatedData);
    // }

    // interpolatedData = malloc(dataTypeSize * nComponent * _n_tot_target_cpl);

    // if (_weights_src_idx == NULL) {
    //   _weights_src_idx = (int **) malloc(sizeof(int *) * _nPart);
    //   for (int i = 0 ; i < _nPart ; i++) _weights_src_idx[i] = NULL;
    //   _weights_src = (double **) malloc(sizeof(double *) * _nPart);
    //   for (int i = 0 ; i < _nPart ; i++) _weights_src[i] = NULL;
    // }

    // // Actually calculate weights with PDM_geom_elem_compute_polygon_barycentric_coordinates for each part
    // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
    //   bool weights_src_empty = false;
    //   int s_weights_src = -1;

    //   if (_weights_src_idx[i_part] == NULL) {
    //     weights_src_empty = true;
    //     _weights_src_idx[i_part] = (int *) malloc(sizeof(int) * (_targets_localization_idx_cpl[cplComm_size - 1][i_part + 1] + 1));
    //     _weights_src_idx[i_part][0] = 0;

    //     s_weights_src = 4 * _targets_localization_idx_cpl[cplComm_size - 1][i_part + 1];
    //     _weights_src[i_part] = (double *) malloc(sizeof(double) * s_weights_src);
    //   }

    //   void *referenceData = referenceField->dataGet(i_part);
    //   // For a cell center field : give the value of the located cell

    //   for (int i_proc = 0 ; i_proc < cplComm_size ; i_proc++) {
    //     // User interpolation
    //     if (interpolationType == CWP_INTERPOLATION_USER) {
    //       CWP_Interp_from_location_t interpolationFunction = referenceField->interpolationFunctionGet();
    //       int n_tgt = _targets_localization_idx_cpl[i_proc][i_part + 1] - _targets_localization_idx_cpl[i_proc][i_part];

    //       int *connecIdx = _mesh->connecIdxGet(i_part);
    //       int *connec = _mesh->connecGet(i_part);
    //       double *coords = _mesh->getVertexCoords(i_part);

    //       int *tgt_pts_location = (int *) malloc(sizeof(int) * n_tgt);
    //       int *tgt_pts_location_p1 = (int *) malloc(sizeof(int) * n_tgt);
    //       auto *tgt_pts_dist = (double *) malloc(sizeof(double) * n_tgt);
    //       auto *tgt_pts_projected_coords = (double *) malloc(3 * sizeof(double) * n_tgt);
    //       int *tgt_pts_bary_coords_idx = NULL;
    //       double *tgt_pts_bary_coords = NULL;

    //       int i = 0;
    //       for (int itarget = _targets_localization_idx_cpl[i_proc][i_part] ; itarget < _targets_localization_idx_cpl[i_proc][i_part + 1] ; itarget++) {
    //           tgt_pts_location[i] = _targets_localization_data_cpl[itarget].lnum;
    //           tgt_pts_location_p1[i] = _targets_localization_data_cpl[itarget].lnum + 1;
    //           tgt_pts_dist[i] = _targets_localization_data_cpl[itarget].distance;
    //           double x_target = _targets_localization_data_cpl[itarget].projectedX;
    //           double y_target = _targets_localization_data_cpl[itarget].projectedY;
    //           double z_target = _targets_localization_data_cpl[itarget].projectedZ;

    //           tgt_pts_projected_coords[3 * i] = x_target;
    //           tgt_pts_projected_coords[3 * i + 1] = y_target;
    //           tgt_pts_projected_coords[3 * i + 2] = z_target;
    //           i++;
    //       }

    //       PDM_geom_elem_compute_polygon_barycentric_coordinates(n_tgt, tgt_pts_location_p1, tgt_pts_projected_coords,
    //                                                             connecIdx, connec, coords, &tgt_pts_bary_coords_idx, &tgt_pts_bary_coords);

    //       void *tmpData = (char *) interpolatedData + dataTypeSize * nComponent * _targets_localization_idx_cpl[i_proc][i_part];

    //       (*interpolationFunction)(CWP_INTERFACE_SURFACE, _n_vtx[i_part], _n_elt[i_part], n_tgt, coords, connecIdx, connec,
    //                                tgt_pts_projected_coords, tgt_pts_location, tgt_pts_dist, tgt_pts_bary_coords_idx, tgt_pts_bary_coords,
    //                                nComponent, referenceFieldType, referenceData, referenceFieldType, tmpData);

    //       if (tgt_pts_location != NULL) free(tgt_pts_location);
    //       if (tgt_pts_location_p1 != NULL) free(tgt_pts_location_p1);
    //       if (tgt_pts_dist != NULL) free(tgt_pts_dist);
    //       if (tgt_pts_projected_coords != NULL) free(tgt_pts_projected_coords);
    //       if (tgt_pts_bary_coords_idx != NULL) free(tgt_pts_bary_coords_idx);
    //       if (tgt_pts_bary_coords != NULL) free(tgt_pts_bary_coords);
    //     }
    //     else { // No user interpolation
    //       if (referenceFieldType == CWP_DOF_LOCATION_CELL_CENTER) {
    //         for (int itarget = _targets_localization_idx_cpl[i_proc][i_part] ; itarget < _targets_localization_idx_cpl[i_proc][i_part + 1] ; itarget++) {
    //           // Index of the corresponding local reference Data.
    //           int iel = _targets_localization_data_cpl[itarget].lnum;
    //           // Index in the interpolated Data array
    //           int interpInd = itarget;
    //           for (int k = 0 ; k < nComponent ; k++) {
    //             memcpy((char *) interpolatedData + dataTypeSize * (nComponent * interpInd + k), (char *) referenceData + dataTypeSize * (nComponent * iel + k), dataTypeSize);
    //           }
    //         } // loop on itarget
    //       } // if referenceFieldType == CWP_DOF_LOCATION_CELL_CENTER
    //       else if (referenceFieldType == CWP_DOF_LOCATION_NODE) {
    //         int *connecIdx = _mesh->connecIdxGet(i_part);
    //         int *connec = _mesh->connecGet(i_part);
    //         double *coords = _mesh->getVertexCoords(i_part);

    //         for (int itarget = _targets_localization_idx_cpl[i_proc][i_part] ; itarget < _targets_localization_idx_cpl[i_proc][i_part + 1] ; itarget++) {
    //           // Index in the interpolated data array
    //           int interpInd = itarget;
    //           int iel = _targets_localization_data_cpl[itarget].lnum;
    //           int ielP1 = iel + 1;

    //           double value = 0.0;
    //           if (_targets_localization_data_cpl[itarget].distance != INFINITY) {
    //             double x_target = _targets_localization_data_cpl[itarget].projectedX;
    //             double y_target = _targets_localization_data_cpl[itarget].projectedY;
    //             double z_target = _targets_localization_data_cpl[itarget].projectedZ;

    //             double tgtCoords[3] = {x_target, y_target, z_target};
    //             double *barCoords;

    //             // TODO Should affect value to weights for any referenceFieldType and also for user interpolation
    //             if (weights_src_empty) {
    //               int *_barCoordsIndex = NULL;
    //               double *_barCoords = NULL;
    //               PDM_geom_elem_compute_polygon_barycentric_coordinates(1, &ielP1, tgtCoords, connecIdx, connec, coords, &_barCoordsIndex, &_barCoords);
    //               int n_elt = connecIdx[iel + 1] - connecIdx[iel];
    //               _weights_src_idx[i_part][itarget + 1] = _weights_src_idx[i_part][itarget] + n_elt;

    //               if (s_weights_src <= _weights_src_idx[i_part][itarget + 1]) {
    //                 s_weights_src *= 2;
    //                 _weights_src[i_part] = (double *) realloc((void *) (_weights_src[i_part]), sizeof(double) * s_weights_src);
    //               }

    //               for (int i = 0 ; i < n_elt ; i++) _weights_src[i_part][_weights_src_idx[i_part][itarget] + i] = _barCoords[i];

    //               free(_barCoordsIndex);
    //               free(_barCoords);
    //             }

    //             barCoords = &(_weights_src[i_part][_weights_src_idx[i_part][itarget]]);

    //             for (int k = 0 ; k < nComponent ; k++) {
    //               value = 0.0;
    //               int k1 = 0;
    //               for (int i_vtx = connecIdx[iel] ; i_vtx < connecIdx[iel + 1] ; i_vtx++) value += barCoords[k1++] * (*(double *) ((char *) referenceData + dataTypeSize * (nComponent * (connec[i_vtx] - 1) + k)));
    //               memcpy((char *) interpolatedData + dataTypeSize * (nComponent * interpInd + k), &value, dataTypeSize);
    //             }
    //           }
    //           else {
    //             for (int k = 0 ; k < nComponent ; k++) {
    //               value = 1000.0;
    //               memcpy((char *) interpolatedData + dataTypeSize * (nComponent * interpInd + k), &value, dataTypeSize);
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }
    // }

    // referenceField->sendBufferSet(interpolatedData);
    // return interpolatedData;
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

  /**********************************************************
  ***********************************************************
  **                                                       **
  **   Process, partition, num triplet location from       **
  **           global numbering functions                  **
  **                                                       **
  ***********************************************************
  **********************************************************/
  
  void SpatialInterpLocation::triplet_location_request() 
  {
    // _id_gnum_location = PDM_gnum_location_create(_nPart_cpl, _nPart, _pdmCplComm);

    // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
    //   PDM_gnum_location_requested_elements_set(_id_gnum_location, i_part, _n_target[i_part], &(_closest_elt_gnum[i_part][0]));
    // }

    // if (!_both_codes_are_local) {
    //   for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
    //     PDM_gnum_location_elements_set(_id_gnum_location, i_part, 0, NULL);
    //   }
    // }
    // else {
    //   for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
    //     Mesh *mesh_cpl = _spatial_interp_cpl->_mesh;
    //     CWP_g_num_t *gnum_elt_cpl = mesh_cpl->GNumEltsGet(i_part);
    //     int n_elt_cpl = mesh_cpl->getPartNElts(i_part);

    //     PDM_gnum_location_elements_set(_id_gnum_location, i_part, n_elt_cpl, gnum_elt_cpl);
    //   }
    // }
  }

  void SpatialInterpLocation::triplet_location_set() 
  {
    // _id_gnum_location = PDM_gnum_location_create(_nPart, _nPart_cpl, _pdmCplComm);

    // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
    //   CWP_g_num_t *gnum_elt = _mesh->GNumEltsGet(i_part);
    //   PDM_gnum_location_elements_set(_id_gnum_location, i_part, _n_elt[i_part], gnum_elt);
    // }

    // if (!_both_codes_are_local)
    //   for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) PDM_gnum_location_requested_elements_set(_id_gnum_location, i_part, 0, NULL);
    // else {
    //   for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
    //     int n_target_cpl = _spatial_interp_cpl->_n_target[i_part];
    //     PDM_gnum_location_requested_elements_set(_id_gnum_location, i_part, n_target_cpl, _spatial_interp_cpl->_closest_elt_gnum[i_part]);
    //   }
    // }
  }

  void SpatialInterpLocation::triplet_location_null_send() 
  {
    // _id_gnum_location = PDM_gnum_location_create(_nPart, _nPart_cpl, _pdmCplComm);

    // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
    //   PDM_gnum_location_elements_set(_id_gnum_location, i_part, 0, NULL);
    // }
    // for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
    //   PDM_gnum_location_requested_elements_set(_id_gnum_location, i_part, 0, NULL);
    // }
  }

  void SpatialInterpLocation::triplet_location_null_recv()
  {
    // _id_gnum_location = PDM_gnum_location_create(_nPart_cpl, _nPart, _pdmCplComm);

    // for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
    //   PDM_gnum_location_elements_set(_id_gnum_location, i_part, 0, NULL);
    // }
    // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
    //   PDM_gnum_location_requested_elements_set(_id_gnum_location, i_part, 0, NULL);
    // }
  }

  void SpatialInterpLocation::triplet_location_compute() const 
  {
    PDM_gnum_location_compute(_id_gnum_location);
  }

  void SpatialInterpLocation::triplet_location_get() 
  {
    // _target_proc_part_num_idx = (int **) malloc(sizeof(int *) * _nPart);
    // _target_proc_part_num = (int **) malloc(sizeof(int *) * _nPart);

    // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
    //   PDM_gnum_location_get(_id_gnum_location, i_part, &(_target_proc_part_num_idx[i_part]), &(_target_proc_part_num[i_part]));
    // }
  }

  void SpatialInterpLocation::triplet_location_get_cpl() 
  {
    _spatial_interp_cpl->_target_proc_part_num_idx = (int **) malloc(sizeof(int *) * _nPart_cpl);
    _spatial_interp_cpl->_target_proc_part_num = (int **) malloc(sizeof(int *) * _nPart_cpl);

    for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
      PDM_gnum_location_get(_id_gnum_location, i_part, &(_spatial_interp_cpl->_target_proc_part_num_idx[i_part]), &(_spatial_interp_cpl->_target_proc_part_num[i_part]));
    }
  }

  /**********************************************************
  ***********************************************************
  **            Communication tree array functions         **
  **                                                       **
  ***********************************************************
  **********************************************************/
  void SpatialInterpLocation::initialization_of_receving_communication_tree_array() 
  {
    if (_targets_localization_idx_cpl == NULL) {
      _targets_localization_idx_cpl = (int **) malloc(sizeof(int *) * cplComm_size);
      for (int i = 0 ; i < cplComm_size ; i++) {
        _targets_localization_idx_cpl[i] = NULL;
      }
    }

    for (int i = 0 ; i < cplComm_size ; i++) {
      if (_targets_localization_idx_cpl[i] == NULL) {
        _targets_localization_idx_cpl[i] = (int *) malloc(sizeof(int) * (1 + _nPart));
      }
      for (int i_part = 0 ; i_part < _nPart + 1 ; i_part++) {
        _targets_localization_idx_cpl[i][i_part] = 0;
      }
    }
  }

  void SpatialInterpLocation::filling_of_sending_communication_tree_array() 
  {
  //   _targets_localization_idx = (int **) malloc(sizeof(int *) * cplComm_size);
  //   for (int i_proc = 0 ; i_proc < cplComm_size ; i_proc++) _targets_localization_idx[i_proc] = NULL;
  //   _process_and_partition_count = (int **) malloc(sizeof(int *) * cplComm_size);

  //   for (int i_proc = 0 ; i_proc < cplComm_size ; i_proc++) {
  //     _targets_localization_idx[i_proc] = (int *) malloc(sizeof(int) * (1 + _nPart_cpl));
  //     _process_and_partition_count[i_proc] = (int *) malloc(sizeof(int) * (1 + _nPart_cpl));

  //     for (int i_part = 0 ; i_part < _nPart_cpl + 1 ; i_part++) {
  //       _targets_localization_idx[i_proc][i_part] = 0;
  //       _process_and_partition_count[i_proc][i_part] = 0;
  //     }
  //   }

  //   for (int i_part = 0 ; i_part < _nPart ; i_part++) {
  //     for (int k = 0 ; k < _n_target[i_part] ; k++) {
  //       int elt_proc = _target_proc_part_num[i_part][_target_proc_part_num_idx[i_part][k]];
  //       int elt_part = _target_proc_part_num[i_part][_target_proc_part_num_idx[i_part][k] + 1];

  //       _targets_localization_idx[elt_proc][elt_part]++;
  //       _process_and_partition_count[elt_proc][elt_part]++;
  //     }
  //   } //end i_part

  //   _transform_to_index(_targets_localization_idx, cplComm_size, _nPart_cpl);

  //   int **idx_proc = (int **) malloc(sizeof(int *) * cplComm_size);
  //   for (int i_proc = 0 ; i_proc < cplComm_size ; i_proc++) {
  //     idx_proc[i_proc] = (int *) malloc(sizeof(int) * _nPart_cpl);
  //     for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
  //       idx_proc[i_proc][i_part] = 0;
  //     }
  //   }

  //   if (_targets_localization_data != NULL) {
  //     free(_targets_localization_data);
  //   }
  //   _targets_localization_data = (target_data *) malloc(sizeof(target_data) * _targets_localization_idx[cplComm_size - 1][_nPart_cpl]);

  //   for (int i_part = 0 ; i_part < _nPart ; i_part++) {
  //     for (int k = 0 ; k < _n_target[i_part] ; k++) {
  //       int elt_proc = _target_proc_part_num[i_part][_target_proc_part_num_idx[i_part][k]];
  //       int elt_part = _target_proc_part_num[i_part][_target_proc_part_num_idx[i_part][k] + 1];
  //       int num = _target_proc_part_num[i_part][_target_proc_part_num_idx[i_part][k] + 2] - 1;

  //       int idx = _targets_localization_idx[elt_proc][elt_part] + idx_proc[elt_proc][elt_part];
  //       //Local Numbering
  //       _targets_localization_data[idx].lnum = num;
  //       //Coupled numbering
  //       _targets_localization_data[idx].origin_part = i_part;
  //       _targets_localization_data[idx].projectedX = _projected[i_part][3 * k];
  //       _targets_localization_data[idx].projectedY = _projected[i_part][3 * k + 1];
  //       _targets_localization_data[idx].projectedZ = _projected[i_part][3 * k + 2];
  //       _targets_localization_data[idx].distance = _distance[i_part][k];
  //       //Coupled process origin
  //       _targets_localization_data[idx].origin_proc = cplComm_rank;
  //       //Coupled origin partition
  //       _targets_localization_data[idx].l_num_origin = k;
  //       idx_proc[elt_proc][elt_part]++;
  //     }
  //   }

  //   for (int i_proc = 0 ; i_proc < cplComm_size ; i_proc++) {
  //     free(idx_proc[i_proc]);
  //   }

  //   for (int i_part = 0 ; i_part < _nPart ; i_part++) {
  //     free(_distance[i_part]);
  //     free(_projected[i_part]);
  //     free(_closest_elt_gnum[i_part]);
  //     free(_target_proc_part_num_idx[i_part]);
  //     free(_target_proc_part_num[i_part]);
  //   }

  //   free(_target_proc_part_num_idx);
  //   free(_target_proc_part_num);
  //   free(idx_proc);
  //   free(_distance);
  //   free(_projected);
  //   free(_closest_elt_gnum);
  }
}

/**
 * \endcond
 */
