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

  SpatialInterpLocation::~SpatialInterpLocation() = default;

  void SpatialInterpLocation::spatialInterpWeightsCompute(CWP_Field_exch_t Texch_t) 
  {
    _Texch_t = Texch_t;

    // In case of withOutPart the user provided not null data only on the root rank (senderRank).
    if (!_both_codes_are_local) {
        if (_Texch_t == CWP_FIELD_EXCH_RECV && _pointsCloudLocation == CWP_DOF_LOCATION_USER) user_targets_gnum_compute();
    }
    else {
      if (_Texch_t == CWP_FIELD_EXCH_SEND) {
          _spatial_interp_cpl->_Texch_t = CWP_FIELD_EXCH_RECV;
          if (_pointsCloudLocation == CWP_DOF_LOCATION_USER) _spatial_interp_cpl->user_targets_gnum_compute();
      }
    }

    // Get informations about the local and the coupled meshes
    info_mesh();

    if ((!_both_codes_are_local && _cpl->commTypeGet() == CWP_COMM_PAR_WITH_PART)
        || (!_both_codes_are_local && _cpl->commTypeGet() == CWP_COMM_PAR_WITHOUT_PART && cplComm_rank == _senderRank)) {
      // Localization
      // Surface and cloud points localization setting
      if (_Texch_t == CWP_FIELD_EXCH_SEND) localization_surface_setting();
      if (_Texch_t == CWP_FIELD_EXCH_RECV) localization_points_cloud_setting();
      // Localization compute, get and free
      localization_compute();
      if (_Texch_t == CWP_FIELD_EXCH_RECV) localization_get();
      localization_free();

      // Communication tree building
      // From a global number obtained the MPI rank and mesh partition of the element
      // Setting and request
      if (_Texch_t == CWP_FIELD_EXCH_RECV) triplet_location_request();
      if (_Texch_t == CWP_FIELD_EXCH_SEND) triplet_location_set();
      // Compute, get and free
      triplet_location_compute();
      if (_Texch_t == CWP_FIELD_EXCH_RECV) triplet_location_get();
      PDM_gnum_location_free(_id_gnum_location, 1);
      // Initialization
      if (_Texch_t == CWP_FIELD_EXCH_RECV) filling_of_sending_communication_tree_array();
      // targets_localization_idx_cpl allocation and init
      if (_Texch_t == CWP_FIELD_EXCH_SEND) initialization_of_receving_communication_tree_array();
      // Communication  of the communication tree index
      if (_Texch_t == CWP_FIELD_EXCH_RECV) data_index_communication_send_p2p();
      if (_Texch_t == CWP_FIELD_EXCH_SEND) data_index_communication_recv_p2p();

      // Communication of the communication tree
      // Preparation
      if (_Texch_t == CWP_FIELD_EXCH_RECV) prepare_data_communication_send();
      if (_Texch_t == CWP_FIELD_EXCH_SEND) prepare_data_communication_recv();
      // MPI asynchronous communication //
      if (_Texch_t == CWP_FIELD_EXCH_RECV) data_communication_send_p2p();
      if (_Texch_t == CWP_FIELD_EXCH_SEND) data_communication_recv_p2p();
      // MPI Wait
      if (_Texch_t == CWP_FIELD_EXCH_RECV) data_communication_wait_send();
      if (_Texch_t == CWP_FIELD_EXCH_SEND) data_communication_wait_recv();

    }
    else if (_both_codes_are_local && _cpl->commTypeGet() == CWP_COMM_PAR_WITH_PART) {
      if (_Texch_t == CWP_FIELD_EXCH_SEND) {
          _spatial_interp_cpl->_Texch_t = CWP_FIELD_EXCH_RECV;
        localization_surface_setting();
        localization_compute();
        localization_get_cpl();
        localization_free();

        triplet_location_set();
        triplet_location_compute();
        triplet_location_get_cpl();

        PDM_gnum_location_free(_id_gnum_location, 1);

        _spatial_interp_cpl->filling_of_sending_communication_tree_array();
        initialization_of_receving_communication_tree_array();

        both_index_communication_p2p();

        _spatial_interp_cpl->prepare_data_communication_send();
        prepare_data_communication_recv();
        both_data_communication_p2p();

        _spatial_interp_cpl->data_communication_wait_send();
        data_communication_wait_recv();
      }
    }
  }

  void *SpatialInterpLocation::interpolate (Field *referenceField) 
  {
    int nComponent = referenceField->nComponentGet();
    CWP_Dof_location_t referenceFieldType = referenceField->locationGet();
    int dataTypeSize = referenceField->dataTypeSizeGet();
    CWP_Interpolation_t interpolationType = referenceField->interpolationTypeGet();
    void *interpolatedData = referenceField->sendBufferGet();

    // Allocate interpolatedData and weights
    if (interpolatedData != NULL) {
      free(interpolatedData);
    }

    interpolatedData = malloc(dataTypeSize * nComponent * _n_tot_target_cpl);

    if (_weights_src_idx == NULL) {
      _weights_src_idx = (int **) malloc(sizeof(int *) * _nb_part);
      for (int i = 0 ; i < _nb_part ; i++) _weights_src_idx[i] = NULL;
      _weights_src = (double **) malloc(sizeof(double *) * _nb_part);
      for (int i = 0 ; i < _nb_part ; i++) _weights_src[i] = NULL;
    }

    // Actually calculate weights with PDM_geom_elem_compute_polygon_barycentric_coordinates for each part
    for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
      bool weights_src_empty = false;
      int s_weights_src = -1;

      if (_weights_src_idx[i_part] == NULL) {
        weights_src_empty = true;
        _weights_src_idx[i_part] = (int *) malloc(sizeof(int) * (_targets_localization_idx_cpl[cplComm_size - 1][i_part + 1] + 1));
        _weights_src_idx[i_part][0] = 0;

        s_weights_src = 4 * _targets_localization_idx_cpl[cplComm_size - 1][i_part + 1];
        _weights_src[i_part] = (double *) malloc(sizeof(double) * s_weights_src);
      }

      void *referenceData = referenceField->dataGet(i_part);
      // For a cell center field : give the value of the located cell

      for (int i_proc = 0 ; i_proc < cplComm_size ; i_proc++) {
        // User interpolation
        if (interpolationType == CWP_INTERPOLATION_USER) {
          CWP_Interp_from_location_t interpolationFunction = referenceField->interpolationFunctionGet();
          int n_tgt = _targets_localization_idx_cpl[i_proc][i_part + 1] - _targets_localization_idx_cpl[i_proc][i_part];

          int *connecIdx = _mesh->connecIdxGet(i_part);
          int *connec = _mesh->connecGet(i_part);
          double *coords = _mesh->getVertexCoords(i_part);

          int *tgt_pts_location = (int *) malloc(sizeof(int) * n_tgt);
          int *tgt_pts_location_p1 = (int *) malloc(sizeof(int) * n_tgt);
          auto *tgt_pts_dist = (double *) malloc(sizeof(double) * n_tgt);
          auto *tgt_pts_projected_coords = (double *) malloc(3 * sizeof(double) * n_tgt);
          int *tgt_pts_bary_coords_idx = NULL;
          double *tgt_pts_bary_coords = NULL;

          int i = 0;
          for (int itarget = _targets_localization_idx_cpl[i_proc][i_part] ; itarget < _targets_localization_idx_cpl[i_proc][i_part + 1] ; itarget++) {
              tgt_pts_location[i] = _targets_localization_data_cpl[itarget].lnum;
              tgt_pts_location_p1[i] = _targets_localization_data_cpl[itarget].lnum + 1;
              tgt_pts_dist[i] = _targets_localization_data_cpl[itarget].distance;
              double x_target = _targets_localization_data_cpl[itarget].projectedX;
              double y_target = _targets_localization_data_cpl[itarget].projectedY;
              double z_target = _targets_localization_data_cpl[itarget].projectedZ;

              tgt_pts_projected_coords[3 * i] = x_target;
              tgt_pts_projected_coords[3 * i + 1] = y_target;
              tgt_pts_projected_coords[3 * i + 2] = z_target;
              i++;
          }

          PDM_geom_elem_compute_polygon_barycentric_coordinates(n_tgt, tgt_pts_location_p1, tgt_pts_projected_coords,
                                                                connecIdx, connec, coords, &tgt_pts_bary_coords_idx, &tgt_pts_bary_coords);

          void *tmpData = (char *) interpolatedData + dataTypeSize * nComponent * _targets_localization_idx_cpl[i_proc][i_part];

          (*interpolationFunction)(CWP_INTERFACE_SURFACE, _n_vtx[i_part], _n_elt[i_part], n_tgt, coords, connecIdx, connec,
                                   tgt_pts_projected_coords, tgt_pts_location, tgt_pts_dist, tgt_pts_bary_coords_idx, tgt_pts_bary_coords,
                                   nComponent, referenceFieldType, referenceData, referenceFieldType, tmpData);

          if (tgt_pts_location != NULL) free(tgt_pts_location);
          if (tgt_pts_location_p1 != NULL) free(tgt_pts_location_p1);
          if (tgt_pts_dist != NULL) free(tgt_pts_dist);
          if (tgt_pts_projected_coords != NULL) free(tgt_pts_projected_coords);
          if (tgt_pts_bary_coords_idx != NULL) free(tgt_pts_bary_coords_idx);
          if (tgt_pts_bary_coords != NULL) free(tgt_pts_bary_coords);
        }
        else { // No user interpolation
          if (referenceFieldType == CWP_DOF_LOCATION_CELL_CENTER) {
            for (int itarget = _targets_localization_idx_cpl[i_proc][i_part] ; itarget < _targets_localization_idx_cpl[i_proc][i_part + 1] ; itarget++) {
              // Index of the corresponding local reference Data.
              int iel = _targets_localization_data_cpl[itarget].lnum;
              // Index in the interpolated Data array
              int interpInd = itarget;
              for (int k = 0 ; k < nComponent ; k++) {
                memcpy((char *) interpolatedData + dataTypeSize * (nComponent * interpInd + k), (char *) referenceData + dataTypeSize * (nComponent * iel + k), dataTypeSize);
              }
            } // loop on itarget
          } // if referenceFieldType == CWP_DOF_LOCATION_CELL_CENTER
          else if (referenceFieldType == CWP_DOF_LOCATION_NODE) {
            int *connecIdx = _mesh->connecIdxGet(i_part);
            int *connec = _mesh->connecGet(i_part);
            double *coords = _mesh->getVertexCoords(i_part);

            for (int itarget = _targets_localization_idx_cpl[i_proc][i_part] ; itarget < _targets_localization_idx_cpl[i_proc][i_part + 1] ; itarget++) {
              // Index in the interpolated data array
              int interpInd = itarget;
              int iel = _targets_localization_data_cpl[itarget].lnum;
              int ielP1 = iel + 1;

              double value = 0.0;
              if (_targets_localization_data_cpl[itarget].distance != INFINITY) {
                double x_target = _targets_localization_data_cpl[itarget].projectedX;
                double y_target = _targets_localization_data_cpl[itarget].projectedY;
                double z_target = _targets_localization_data_cpl[itarget].projectedZ;

                double tgtCoords[3] = {x_target, y_target, z_target};
                double *barCoords;

                // TODO Should affect value to weights for any referenceFieldType and also for user interpolation
                if (weights_src_empty) {
                  int *_barCoordsIndex = NULL;
                  double *_barCoords = NULL;
                  PDM_geom_elem_compute_polygon_barycentric_coordinates(1, &ielP1, tgtCoords, connecIdx, connec, coords, &_barCoordsIndex, &_barCoords);
                  int n_elt = connecIdx[iel + 1] - connecIdx[iel];
                  _weights_src_idx[i_part][itarget + 1] = _weights_src_idx[i_part][itarget] + n_elt;

                  if (s_weights_src <= _weights_src_idx[i_part][itarget + 1]) {
                    s_weights_src *= 2;
                    _weights_src[i_part] = (double *) realloc((void *) (_weights_src[i_part]), sizeof(double) * s_weights_src);
                  }

                  for (int i = 0 ; i < n_elt ; i++) _weights_src[i_part][_weights_src_idx[i_part][itarget] + i] = _barCoords[i];

                  free(_barCoordsIndex);
                  free(_barCoords);
                }

                barCoords = &(_weights_src[i_part][_weights_src_idx[i_part][itarget]]);

                for (int k = 0 ; k < nComponent ; k++) {
                  value = 0.0;
                  int k1 = 0;
                  for (int i_vtx = connecIdx[iel] ; i_vtx < connecIdx[iel + 1] ; i_vtx++) value += barCoords[k1++] * (*(double *) ((char *) referenceData + dataTypeSize * (nComponent * (connec[i_vtx] - 1) + k)));
                  memcpy((char *) interpolatedData + dataTypeSize * (nComponent * interpInd + k), &value, dataTypeSize);
                }
              }
              else {
                for (int k = 0 ; k < nComponent ; k++) {
                  value = 1000.0;
                  memcpy((char *) interpolatedData + dataTypeSize * (nComponent * interpInd + k), &value, dataTypeSize);
                }
              }
            }
          }
        }
      }
    }

    referenceField->sendBufferSet(interpolatedData);
    return interpolatedData;
  }

  void SpatialInterpLocation::init(Coupling *coupling, CWP_Dof_location_t pointsCloudLocation, bool slave) 
  {
    SpatialInterp::init(coupling, pointsCloudLocation, slave);

    interpolation_time = CWP_INTERP_AT_SEND;

    CouplingDB *cplDB = _cpl->couplingDBGet();
    string cplId = coupling->IdGet();
    if (_both_codes_are_local && !_slave) {
      Coupling coupling_cpl = cplDB->couplingGet(*_coupledCodeProperties, cplId);
      _spatial_interp_cpl = dynamic_cast<SpatialInterpLocation *>( coupling_cpl.spatialInterpGet(_pointsCloudLocation));
      _spatial_interp_cpl->_spatial_interp_cpl = this;
    }

    int tmp1 = -1000000; // For dummy communications
    if (!_both_codes_are_local) {
      if (_id < _id_cpl) {
        MPI_Bcast(&_nb_part_cpl, 1, MPI_INT, _senderRank, _cplComm);
        MPI_Bcast(&tmp1, 1, MPI_INT, _senderRank_cpl, _cplComm);
      }
      else {
        MPI_Bcast(&tmp1, 1, MPI_INT, _senderRank_cpl, _cplComm);
        MPI_Bcast(&_nb_part_cpl, 1, MPI_INT, _senderRank, _cplComm);
      }
    }
    else if (!_slave) {
      if (_id < _id_cpl) {
        MPI_Bcast(&_nb_part_cpl, 1, MPI_INT, _senderRank, _cplComm);
        MPI_Bcast(&(_spatial_interp_cpl->_nb_part_cpl), 1, MPI_INT, _senderRank_cpl, _cplComm);
      }
      else {
        MPI_Bcast(&(_spatial_interp_cpl->_nb_part_cpl), 1, MPI_INT, _senderRank_cpl, _cplComm);
        MPI_Bcast(&_nb_part_cpl, 1, MPI_INT, _senderRank, _cplComm);
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
  void SpatialInterpLocation::localization_points_cloud_setting() 
  {
    PDM_error(__FILE__, __LINE__, 0, "Unknown location method.\n");
  }

  void SpatialInterpLocation::localization_null_setting_send() 
  {
    PDM_error(__FILE__, __LINE__, 0, "Unknown location method.\n");
  }

  void SpatialInterpLocation::localization_null_setting_recv() 
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

  void SpatialInterpLocation::localization_get_cpl() 
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
    _id_gnum_location = PDM_gnum_location_create(_nb_part_cpl, _nb_part, _pdm_cplComm);

    for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
      PDM_gnum_location_requested_elements_set(_id_gnum_location, i_part, _n_target[i_part], &(_closest_elt_gnum[i_part][0]));
    }

    if (!_both_codes_are_local) {
      for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
        PDM_gnum_location_elements_set(_id_gnum_location, i_part, 0, NULL);
      }
    }
    else {
      for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
        Mesh *mesh_cpl = _spatial_interp_cpl->_mesh;
        CWP_g_num_t *gnum_elt_cpl = mesh_cpl->GNumEltsGet(i_part);
        int n_elt_cpl = mesh_cpl->getPartNElts(i_part);

        PDM_gnum_location_elements_set(_id_gnum_location, i_part, n_elt_cpl, gnum_elt_cpl);
      }
    }
  }

  void SpatialInterpLocation::triplet_location_set() 
  {
    _id_gnum_location = PDM_gnum_location_create(_nb_part, _nb_part_cpl, _pdm_cplComm);

    for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
      CWP_g_num_t *gnum_elt = _mesh->GNumEltsGet(i_part);
      PDM_gnum_location_elements_set(_id_gnum_location, i_part, _n_elt[i_part], gnum_elt);
    }

    if (!_both_codes_are_local)
      for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) PDM_gnum_location_requested_elements_set(_id_gnum_location, i_part, 0, NULL);
    else {
      for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
        int n_target_cpl = _spatial_interp_cpl->_n_target[i_part];
        PDM_gnum_location_requested_elements_set(_id_gnum_location, i_part, n_target_cpl, _spatial_interp_cpl->_closest_elt_gnum[i_part]);
      }
    }
  }

  void SpatialInterpLocation::triplet_location_null_send() 
  {
    _id_gnum_location = PDM_gnum_location_create(_nb_part, _nb_part_cpl, _pdm_cplComm);

    for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
      PDM_gnum_location_elements_set(_id_gnum_location, i_part, 0, NULL);
    }
    for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
      PDM_gnum_location_requested_elements_set(_id_gnum_location, i_part, 0, NULL);
    }
  }

  void SpatialInterpLocation::triplet_location_null_recv()
  {
    _id_gnum_location = PDM_gnum_location_create(_nb_part_cpl, _nb_part, _pdm_cplComm);

    for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
      PDM_gnum_location_elements_set(_id_gnum_location, i_part, 0, NULL);
    }
    for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
      PDM_gnum_location_requested_elements_set(_id_gnum_location, i_part, 0, NULL);
    }
  }

  void SpatialInterpLocation::triplet_location_compute() const 
  {
    PDM_gnum_location_compute(_id_gnum_location);
  }

  void SpatialInterpLocation::triplet_location_get() 
  {
    _target_proc_part_num_idx = (int **) malloc(sizeof(int *) * _nb_part);
    _target_proc_part_num = (int **) malloc(sizeof(int *) * _nb_part);

    for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
      PDM_gnum_location_get(_id_gnum_location, i_part, &(_target_proc_part_num_idx[i_part]), &(_target_proc_part_num[i_part]));
    }
  }

  void SpatialInterpLocation::triplet_location_get_cpl() 
  {
    _spatial_interp_cpl->_target_proc_part_num_idx = (int **) malloc(sizeof(int *) * _nb_part_cpl);
    _spatial_interp_cpl->_target_proc_part_num = (int **) malloc(sizeof(int *) * _nb_part_cpl);

    for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
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
        _targets_localization_idx_cpl[i] = (int *) malloc(sizeof(int) * (1 + _nb_part));
      }
      for (int i_part = 0 ; i_part < _nb_part + 1 ; i_part++) {
        _targets_localization_idx_cpl[i][i_part] = 0;
      }
    }
  }

  void SpatialInterpLocation::filling_of_sending_communication_tree_array() 
  {
    _targets_localization_idx = (int **) malloc(sizeof(int *) * cplComm_size);
    for (int i_proc = 0 ; i_proc < cplComm_size ; i_proc++) _targets_localization_idx[i_proc] = NULL;
    _process_and_partition_count = (int **) malloc(sizeof(int *) * cplComm_size);

    for (int i_proc = 0 ; i_proc < cplComm_size ; i_proc++) {
      _targets_localization_idx[i_proc] = (int *) malloc(sizeof(int) * (1 + _nb_part_cpl));
      _process_and_partition_count[i_proc] = (int *) malloc(sizeof(int) * (1 + _nb_part_cpl));

      for (int i_part = 0 ; i_part < _nb_part_cpl + 1 ; i_part++) {
        _targets_localization_idx[i_proc][i_part] = 0;
        _process_and_partition_count[i_proc][i_part] = 0;
      }
    }

    for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
      for (int k = 0 ; k < _n_target[i_part] ; k++) {
        int elt_proc = _target_proc_part_num[i_part][_target_proc_part_num_idx[i_part][k]];
        int elt_part = _target_proc_part_num[i_part][_target_proc_part_num_idx[i_part][k] + 1];

        _targets_localization_idx[elt_proc][elt_part]++;
        _process_and_partition_count[elt_proc][elt_part]++;
      }
    } //end i_part

    _transform_to_index(_targets_localization_idx, cplComm_size, _nb_part_cpl);

    int **idx_proc = (int **) malloc(sizeof(int *) * cplComm_size);
    for (int i_proc = 0 ; i_proc < cplComm_size ; i_proc++) {
      idx_proc[i_proc] = (int *) malloc(sizeof(int) * _nb_part_cpl);
      for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
        idx_proc[i_proc][i_part] = 0;
      }
    }

    if (_targets_localization_data != NULL) {
      free(_targets_localization_data);
    }
    _targets_localization_data = (target_data *) malloc(sizeof(target_data) * _targets_localization_idx[cplComm_size - 1][_nb_part_cpl]);

    for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
      for (int k = 0 ; k < _n_target[i_part] ; k++) {
        int elt_proc = _target_proc_part_num[i_part][_target_proc_part_num_idx[i_part][k]];
        int elt_part = _target_proc_part_num[i_part][_target_proc_part_num_idx[i_part][k] + 1];
        int num = _target_proc_part_num[i_part][_target_proc_part_num_idx[i_part][k] + 2] - 1;

        int idx = _targets_localization_idx[elt_proc][elt_part] + idx_proc[elt_proc][elt_part];
        //Local Numbering
        _targets_localization_data[idx].lnum = num;
        //Coupled numbering
        _targets_localization_data[idx].origin_part = i_part;
        _targets_localization_data[idx].projectedX = _projected[i_part][3 * k];
        _targets_localization_data[idx].projectedY = _projected[i_part][3 * k + 1];
        _targets_localization_data[idx].projectedZ = _projected[i_part][3 * k + 2];
        _targets_localization_data[idx].distance = _distance[i_part][k];
        //Coupled process origin
        _targets_localization_data[idx].origin_proc = cplComm_rank;
        //Coupled origin partition
        _targets_localization_data[idx].l_num_origin = k;
        idx_proc[elt_proc][elt_part]++;
      }
    }

    for (int i_proc = 0 ; i_proc < cplComm_size ; i_proc++) {
      free(idx_proc[i_proc]);
    }

    for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
      free(_distance[i_part]);
      free(_projected[i_part]);
      free(_closest_elt_gnum[i_part]);
      free(_target_proc_part_num_idx[i_part]);
      free(_target_proc_part_num[i_part]);
    }

    free(_target_proc_part_num_idx);
    free(_target_proc_part_num);
    free(idx_proc);
    free(_distance);
    free(_projected);
    free(_closest_elt_gnum);
  }
}

/**
 * \endcond
 */
