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

#include "spatialInterpLocationMeshLocation.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"

/**
 * \cond
 */

namespace cwipi {

  void SpatialInterpLocationMeshLocation::localization_init() {

    if (!_coupledCodeProperties->localCodeIs()) {

      _id_pdm = PDM_mesh_location_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm);

      PDM_mesh_location_method_set(_id_pdm, _location_method);
      PDM_mesh_location_tolerance_set(_id_pdm, _tolerance);
  
      if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

        PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _nPart);
        // PDM_mesh_location_mesh_global_data_set(_id_pdm, _nPart_cpl);
            //To be continued


      }

    }

    else {

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        _id_pdm = PDM_mesh_location_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm);

        SpatialInterpLocationMeshLocation *cpl_spatial_interp;

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet(); 

          cpl_spatial_interp = 
            dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

            //To be continued 

        }

        else {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 

          cpl_spatial_interp = 
            dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

            //To be continued

        }

        cpl_spatial_interp->_id_pdm = _id_pdm;

      }

    }
  }

  void SpatialInterpLocationMeshLocation::localization_points_cloud_setting() {

      if (!_coupledCodeProperties->localCodeIs()) {

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

              //To be continued


        }

      }

      else {

        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet(); 

            SpatialInterpLocationMeshLocation * cpl_spatial_interp_send = 
              dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

              //To be continued 

          }

          else {

            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 

            SpatialInterpLocationMeshLocation * cpl_spatial_interp_recv = 
              dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

              //To be continued

          }


        }

      }

        // if (_coupledCodeProperties->localCodeIs() && cp < dl) { 
        //     if 



        // _id_pdm = PDM_mesh_location_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm);

        // PDM_mesh_location_method_set(_id_pdm, _location_method);
        // PDM_mesh_location_tolerance_set(_id_pdm, _tolerance);
        // PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _nPart);
        // PDM_mesh_location_mesh_global_data_set(_id_pdm, _nPart_cpl);

        // for (int i_part = 0 ; i_part < _nPart ; i_part++) PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, _n_target[i_part], _coords_target[i_part], _gnum_target[i_part]);

        // if (!_both_codes_are_local) {
        //     CWP_Interface_t interf_dim_cpl = _cpl->entitiesDimGet();

        //     for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
        //         if (interf_dim_cpl == CWP_INTERFACE_SURFACE) PDM_mesh_location_part_set_2d(_id_pdm, i_part, 0, NULL, NULL, NULL, 0, NULL, NULL, NULL, 0, NULL, NULL);
        //         else if (interf_dim_cpl == CWP_INTERFACE_VOLUME) PDM_mesh_location_part_set(_id_pdm, i_part, 0, NULL, NULL, NULL, 0, NULL, NULL, NULL, 0, NULL, NULL);
        //     }
        // }
        // else {
        //     CWP_Interface_t interf_dim_cpl = _spatial_interp_cpl->_cpl->entitiesDimGet();
        //     Mesh *mesh_cpl = _spatial_interp_cpl->_mesh;

        //     for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
        //         int n_vtx = mesh_cpl->getPartNVertex(i_part);
        //         int n_face = mesh_cpl->getNFace(i_part);

        //         double *coords = mesh_cpl->getVertexCoords(i_part);

        //         CWP_g_num_t *vtx_gnum = mesh_cpl->getVertexGNum(i_part);
        //         CWP_g_num_t *elt_gnum = mesh_cpl->GNumEltsGet(i_part);

        //         if (interf_dim_cpl == CWP_INTERFACE_SURFACE) {
        //             int n_edge = mesh_cpl->getNEdge(i_part);

        //             int *face_edge_idx = mesh_cpl->getFaceEdgeIndex(i_part);
        //             int *face_edge = mesh_cpl->getFaceEdge(i_part);
        //             int *edge_vtx_idx = mesh_cpl->getEdgeVtxIndex(i_part);
        //             int *edge_vtx = mesh_cpl->getEdgeVtx(i_part);

        //             PDM_mesh_location_part_set_2d(_id_pdm, i_part, n_face, face_edge_idx, face_edge, elt_gnum, n_edge, edge_vtx_idx, edge_vtx, NULL, n_vtx, coords, vtx_gnum);
        //         }
        //         if (interf_dim_cpl == CWP_INTERFACE_VOLUME) {
        //             int n_cell = mesh_cpl->getNCell(i_part);

        //             int *cell_face_idx = mesh_cpl->getCellFaceIndex(i_part);
        //             int *cell_face = mesh_cpl->getCellFace(i_part);
        //             int *face_vtx_idx = mesh_cpl->getFaceVtxIndex(i_part);
        //             int *face_vtx = mesh_cpl->getFaceVtx(i_part);

        //             PDM_mesh_location_part_set(_id_pdm, i_part, n_cell, cell_face_idx, cell_face, elt_gnum, n_face, face_vtx_idx, face_vtx, NULL, n_vtx, coords, vtx_gnum);
        //         }
        //     }
        // }
    }

    void SpatialInterpLocationMeshLocation::localization_null_setting_send() {
        //  _id_pdm = PDM_mesh_location_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm);

        // PDM_mesh_location_method_set(_id_pdm, _location_method);
        // PDM_mesh_location_tolerance_set(_id_pdm, _tolerance);
        // PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _nPart_cpl);
        // PDM_mesh_location_mesh_global_data_set(_id_pdm, _nPart);

        // for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, 0, NULL, NULL);

        CWP_Interface_t interf_dim = _cpl->entitiesDimGet();

        // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
        //     if (interf_dim == CWP_INTERFACE_SURFACE) PDM_mesh_location_part_set_2d(_id_pdm, i_part, 0, NULL, NULL, NULL, 0, NULL, NULL, NULL, 0, NULL, NULL);
        //     else if (interf_dim == CWP_INTERFACE_VOLUME) PDM_mesh_location_part_set(_id_pdm, i_part, 0, NULL, NULL, NULL, 0, NULL, NULL, NULL, 0, NULL, NULL);
        // }
    }

    void SpatialInterpLocationMeshLocation::localization_null_setting_recv() {
        // _id_pdm = PDM_mesh_location_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm);

        // PDM_mesh_location_method_set(_id_pdm, _location_method);
        // PDM_mesh_location_tolerance_set(_id_pdm, _tolerance);
        // PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _nPart);
        // PDM_mesh_location_mesh_global_data_set(_id_pdm, _nPart_cpl);

        // for (int i_part = 0 ; i_part < _nPart ; i_part++) PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, 0, NULL, NULL);

        // CWP_Interface_t interf_dim_cpl = _spatial_interp_cpl->_cpl->entitiesDimGet();
        // for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
        //     if (interf_dim_cpl == CWP_INTERFACE_SURFACE) PDM_mesh_location_part_set_2d(_id_pdm, i_part, 0, NULL, NULL, NULL, 0, NULL, NULL, NULL, 0, NULL, NULL);
        //     else if (interf_dim_cpl == CWP_INTERFACE_VOLUME) PDM_mesh_location_part_set(_id_pdm, i_part, 0, NULL, NULL, NULL, 0, NULL, NULL, NULL, 0, NULL, NULL);
        // }
    }

    void SpatialInterpLocationMeshLocation::localization_surface_setting() {
        // _id_pdm = PDM_mesh_location_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm);

        // PDM_mesh_location_method_set(_id_pdm, _location_method);
        // PDM_mesh_location_tolerance_set(_id_pdm, _tolerance);
        // PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _nPart_cpl);
        // PDM_mesh_location_mesh_global_data_set(_id_pdm, _nPart);

        // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
        //     int n_vtx = _mesh->getPartNVertex(i_part);
        //     int n_face = _mesh->getNFace(i_part);

        //     double *coords = _mesh->getVertexCoords(i_part);

        //     CWP_g_num_t *vtx_gnum = _mesh->getVertexGNum(i_part);
        //     CWP_g_num_t *elt_gnum = _mesh->GNumEltsGet(i_part);

        //     CWP_Interface_t interf_dim = _cpl->entitiesDimGet();

        //     if (interf_dim == CWP_INTERFACE_SURFACE) {
        //         int n_edge = _mesh->getNEdge(i_part);

        //         int *face_edge_idx = _mesh->getFaceEdgeIndex(i_part);
        //         int *face_edge = _mesh->getFaceEdge(i_part);
        //         int *edge_vtx_idx = _mesh->getEdgeVtxIndex(i_part);
        //         int *edge_vtx = _mesh->getEdgeVtx(i_part);

        //         PDM_mesh_location_part_set_2d(_id_pdm, i_part, n_face, face_edge_idx, face_edge, elt_gnum, n_edge, edge_vtx_idx, edge_vtx, NULL, n_vtx, coords, vtx_gnum);
        //     }
        //     else if (interf_dim == CWP_INTERFACE_VOLUME) {
        //         int n_cell = _mesh->getNCell(i_part);

        //         int *cell_face_idx = _mesh->getCellFaceIndex(i_part);
        //         int *cell_face = _mesh->getCellFace(i_part);
        //         int *face_vtx_idx = _mesh->getFaceVtxIndex(i_part);
        //         int *face_vtx = _mesh->getFaceVtx(i_part);

        //         PDM_mesh_location_part_set(_id_pdm, i_part, n_cell, cell_face_idx, cell_face, elt_gnum, n_face, face_vtx_idx, face_vtx, NULL, n_vtx, coords, vtx_gnum);
        //     }
        // }

        // if (!_both_codes_are_local)
        //     for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, 0, NULL, NULL);
        // else {
        //     for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
        //         int n_target_cpl = _spatial_interp_cpl->_n_target[i_part];
        //         CWP_g_num_t *gnum_target_cpl = _spatial_interp_cpl->_gnum_target[i_part];
        //         double *coords_target_cpl = _spatial_interp_cpl->_coords_target[i_part];
        //         PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, n_target_cpl, coords_target_cpl, gnum_target_cpl);
        //     }
        // }
    }

    void SpatialInterpLocationMeshLocation::localization_compute() {
        PDM_mesh_location_compute(_id_pdm);
        PDM_mesh_location_dump_times(_id_pdm);
    }

    void SpatialInterpLocationMeshLocation::localization_get_cpl() {
//         _spatial_interp_cpl->_distance = (double **) malloc(sizeof(double *) * _nPart_cpl);
//         _spatial_interp_cpl->_projected = (double **) malloc(sizeof(double *) * _nPart_cpl);
//         _spatial_interp_cpl->_closest_elt_gnum = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t * ) * _nPart_cpl);

//         for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
//             int n_target_cpl = _spatial_interp_cpl->_n_target[i_part];

// //            PDM_mesh_location_get(_id_pdm, 0, i_part, &(_spatial_interp_cpl->_closest_elt_gnum[i_part]), &unused_weights_idx, &unused_weights, &(_spatial_interp_cpl->_projected[i_part]));
//             int n_unlocated = PDM_mesh_location_n_unlocated_get (_id_pdm, 0, i_part);

//             assert (n_unlocated == 0);

//             PDM_mesh_location_point_location_get (_id_pdm,
//                                                   0,
//                                                   i_part,
//                                                   &(_spatial_interp_cpl -> _closest_elt_gnum[i_part]),
//                                                   &(_spatial_interp_cpl -> _distance[i_part]),
//                                                   &(_spatial_interp_cpl -> _projected[i_part]));
//             _spatial_interp_cpl->_distance[i_part] = (double *) malloc(sizeof(double) * n_target_cpl);

//             for (int i = 0 ; i < n_target_cpl ; i++) {
//                 _spatial_interp_cpl->_distance[i_part][i] = 0.;

//                 if (_spatial_interp_cpl->_closest_elt_gnum[i_part][i] > CWP_g_num_t(_n_g_elt_cpl_over_part) ||
//                     _spatial_interp_cpl->_closest_elt_gnum[i_part][i] < CWP_g_num_t(1) ||
//                     _spatial_interp_cpl->_distance[i_part][i] > 0.1) {
//                     _spatial_interp_cpl->_closest_elt_gnum[i_part][i] = CWP_g_num_t(1);
//                     _spatial_interp_cpl->_distance[i_part][i] = INFINITY;
//                 }
//             }
//         }
    }

    void SpatialInterpLocationMeshLocation::localization_get() {
//         _distance = (double **) malloc(sizeof(double *) * _nPart);
//         _projected = (double **) malloc(sizeof(double *) * _nPart);
//         _closest_elt_gnum = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t * ) * _nPart);

//         for (int i_part = 0 ; i_part < _nPart ; i_part++) {
// //            PDM_mesh_location_get(_id_pdm, 0, i_part, &(_closest_elt_gnum[i_part]), &unused_weights_idx, &unused_weights, &(_projected[i_part]));
//             _distance[i_part] = (double *) malloc (sizeof(double) * _n_target[i_part]);

//             int n_unlocated = PDM_mesh_location_n_unlocated_get (_id_pdm, 0, i_part);

//             assert (n_unlocated == 0);

//             PDM_mesh_location_point_location_get (_id_pdm,
//                                                   0,
//                                                   i_part,
//                                                   &(_closest_elt_gnum[i_part]),
//                                                   &(_distance[i_part]),
//                                                   &(_projected[i_part]));

//             _distance[i_part] = (double *) malloc(sizeof(double) * _n_target[i_part]);
//             double *coords_target = _coords_target[i_part];

//             for (int i = 0 ; i < _n_target[i_part] ; i++) {
//                 _distance[i_part][i] = 0.;
//                 for (int j = 0 ; j < 3 ; j++) {
//                     double d = _projected[i_part][3 * i + j] - coords_target[3 * i + j];
//                     _distance[i_part][i] += d * d;
//                 }
//                 if (_closest_elt_gnum[i_part][i] > CWP_g_num_t(_n_g_elt_cpl_over_part) ||
//                     _closest_elt_gnum[i_part][i] < CWP_g_num_t(1) ||
//                     _distance[i_part][i] > 0.1) {
//                     _closest_elt_gnum[i_part][i] = CWP_g_num_t(1);
//                     _distance[i_part][i] = INFINITY;
//                 }
//             }
//         }
    }

    void SpatialInterpLocationMeshLocation::localization_free() {
        PDM_mesh_location_free(_id_pdm, 1);
    }
}

/**
 * \endcond
 */
