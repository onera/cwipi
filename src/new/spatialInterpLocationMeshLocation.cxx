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
#include "pdm_logging.h"

/**
 * \cond
 */

namespace cwipi {

  void 
  SpatialInterpLocationMeshLocation::localization_init
  (
  ) 
  {   

    if (!_coupledCodeProperties->localCodeIs()) {

      _id_pdm = PDM_mesh_location_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

      PDM_mesh_location_reverse_results_enable (_id_pdm);

      PDM_mesh_location_method_set(_id_pdm, _location_method);
      PDM_mesh_location_tolerance_set(_id_pdm, _tolerance);
  
      if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

        PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _nPart);

      }

      else {

        PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _cplNPart);

      }

    }

    else {

      // Attention : 
      //     - creation d'un objet unique pdm sur le plus petit des 2 id + Copie de l'id dans l'objet du code couple
      //     - rien a faire pour id le plus grand

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        printf("localization_init - 2.1\n");
        fflush(stdout);

        _id_pdm = PDM_mesh_location_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        SpatialInterpLocationMeshLocation *cpl_spatial_interp;

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          printf("localization_init - 2.2\n");
          fflush(stdout);

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet(); 

          cpl_spatial_interp = 
            dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }

        else {

          printf("localization_init - 2.3\n");
          fflush(stdout);

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 

          cpl_spatial_interp = 
            dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }

        cpl_spatial_interp->_id_pdm = _id_pdm;

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

//          PDM_mesh_location_mesh_global_data_set(_id_pdm, _cplNPart);
          PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _nPart);

          printf("localization_init - 2.4\n");
          fflush(stdout);

        }

        else {

//          PDM_mesh_location_mesh_global_data_set(_id_pdm, _nPart);
          PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _cplNPart);

          printf("localization_init - 2.5\n");
          fflush(stdout);

        }
      }
    }
  }


  SpatialInterpLocationMeshLocation::~SpatialInterpLocationMeshLocation
  (
  )
  {

  }




  void 
  SpatialInterpLocationMeshLocation::localization_points_cloud_setting
  (
  ) 
  {

    if (!_coupledCodeProperties->localCodeIs()) {

      if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

        for (int i_part = 0 ; i_part < _nPart ; i_part++) { 

          const double *part_coord = NULL;
          const PDM_g_num_t *part_gnum = NULL;
          int part_n = -1;

          if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            part_gnum = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
            part_coord = _mesh->eltCentersGet (i_part);
            part_n = _mesh->getPartNElts (i_part);

          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
            part_gnum = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
            part_coord = _mesh->getVertexCoords (i_part);
            part_n = _mesh->getPartNVertex (i_part);            
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            part_gnum = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
            part_coord = _cpl->userTargetCoordsGet (i_part);
            part_n = _cpl->userTargetNGet (i_part);
          }

          PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, part_n, (double *) part_coord, (PDM_g_num_t*) part_gnum);
        }
      }

      else {
        for (int i_part = 0 ; i_part < _nPart ; i_part++) { 
          PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, 0, NULL, NULL);
        }
      }
    }

    else {

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          for (int i_part = 0 ; i_part < _nPart ; i_part++) { 

            const double *part_coord = NULL;
            const PDM_g_num_t *part_gnum = NULL;
            int part_n = -1;

            if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              part_gnum = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
              part_coord = _mesh->eltCentersGet (i_part);
              part_n = _mesh->getPartNElts (i_part);

            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              part_gnum = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
              part_coord = _mesh->getVertexCoords (i_part);
              part_n = _mesh->getPartNVertex (i_part);            
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
              part_gnum = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
              part_coord = _cpl->userTargetCoordsGet (i_part);
              part_n = _cpl->userTargetNGet (i_part);
            }

            PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, part_n, (double *) part_coord, (PDM_g_num_t*) part_gnum);

          }

        }

        else {
  
          cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();

          for (int i_part = 0 ; i_part < _nPart ; i_part++) { 

            const double *part_coord = NULL;
            const PDM_g_num_t *part_gnum = NULL;
            int part_n = -1;

            if (_coupledCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              part_gnum = (const PDM_g_num_t *) cpl_mesh->GNumEltsGet (i_part);
              part_coord = cpl_mesh->eltCentersGet (i_part);
              part_n = cpl_mesh->getPartNElts (i_part);

            }
            else if (_coupledCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              part_gnum = (const PDM_g_num_t *) cpl_mesh->getVertexGNum (i_part);
              part_coord = cpl_mesh->getVertexCoords (i_part);
              part_n = cpl_mesh->getPartNVertex (i_part);            
            }
            else if (_coupledCodeDofLocation == CWP_DOF_LOCATION_USER) {
              part_gnum = (const PDM_g_num_t *) cpl_cpl.userTargetGNumGet (i_part);
              part_coord = cpl_cpl.userTargetCoordsGet (i_part);
              part_n = cpl_cpl.userTargetNGet (i_part);
            }

            PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, part_n, (double *) part_coord, (PDM_g_num_t*) part_gnum);
          }
        }
      }
    }
  }



  void 
  SpatialInterpLocationMeshLocation::localization_surface_setting
  (
  )
  {

    if (!_coupledCodeProperties->localCodeIs()) {
      if (_mesh->getNFace(0) == 0) {
        printf("No faces, using nodal\n");

        PDM_mesh_location_shared_nodal_mesh_set(_id_pdm, _pdm_CplNodal);

      }
      else {
        CWP_Interface_t interf_dim = _cpl->entitiesDimGet();
        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {

          for (int i_part = 0 ; i_part < _nPart ; i_part++) { 
            int n_vtx = _mesh->getPartNVertex(i_part);
            int n_face = _mesh->getNFace(i_part);
            double *coords = _mesh->getVertexCoords(i_part);
            CWP_g_num_t *vtx_gnum = _mesh->getVertexGNum(i_part);
            CWP_g_num_t *elt_gnum = _mesh->GNumEltsGet(i_part);

            if (interf_dim == CWP_INTERFACE_SURFACE) {
              int n_edge = _mesh->getNEdge(i_part);
              int *face_edge_idx = _mesh->getFaceEdgeIndex(i_part);
              int *face_edge = _mesh->getFaceEdge(i_part);
              int *edge_vtx_idx = _mesh->getEdgeVtxIndex(i_part);
              int *edge_vtx = _mesh->getEdgeVtx(i_part);

              PDM_mesh_location_part_set_2d (_id_pdm, 
                                             i_part, 
                                             n_face, 
                                             face_edge_idx,
                                             face_edge, 
                                             elt_gnum, 
                                             n_edge, 
                                             edge_vtx_idx, 
                                             edge_vtx, 
                                             NULL, 
                                             n_vtx, 
                                             coords, 
                                             vtx_gnum);
            }

            else if (interf_dim == CWP_INTERFACE_VOLUME) {

              int n_cell = _mesh->getNCell(i_part);
              int *cell_face_idx = _mesh->getCellFaceIndex(i_part);
              int *cell_face = _mesh->getCellFace(i_part);
              int *face_vtx_idx = _mesh->getFaceVtxIndex(i_part);
              int *face_vtx = _mesh->getFaceVtx(i_part);

              PDM_mesh_location_part_set (_id_pdm, 
                                          i_part, 
                                          n_cell, 
                                          cell_face_idx, 
                                          cell_face, 
                                          elt_gnum, 
                                          n_face, 
                                          face_vtx_idx, 
                                          face_vtx, 
                                          NULL, 
                                          n_vtx, 
                                          coords, 
                                          vtx_gnum);
            }
          }
        }

        else {
//          CWP_Interface_t interf_dim = _cpl->entitiesDimGet();
          for (int i_part = 0 ; i_part < _nPart ; i_part++) {
            if (interf_dim == CWP_INTERFACE_SURFACE) {
              PDM_mesh_location_part_set_2d(_id_pdm,
                                            i_part,
                                            0,
                                            fake_idx,
                                            NULL,
                                            NULL,
                                            0,
                                            fake_idx,
                                            NULL,
                                            NULL,
                                            0,
                                            NULL,
                                            NULL);
            }

            else if (interf_dim == CWP_INTERFACE_VOLUME) {
              PDM_mesh_location_part_set(_id_pdm,
                                         i_part,
                                         0,
                                         fake_idx,
                                         NULL,
                                         NULL,
                                         0,
                                         fake_idx,
                                         NULL,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL);
            }
          }
        }
      }
    }

    else {
      CWP_Interface_t interf_dim = _cpl->entitiesDimGet();

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        if (_mesh->getNFace(0) == 0) {
          printf("No faces, using nodal\n");
          PDM_mesh_location_shared_nodal_mesh_set(_id_pdm, _pdm_CplNodal);
        }
        else {
          if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {

            for (int i_part = 0 ; i_part < _nPart ; i_part++) {
              int n_vtx = _mesh->getPartNVertex(i_part);
              int n_face = _mesh->getNFace(i_part);
              double *coords = _mesh->getVertexCoords(i_part);
              CWP_g_num_t *vtx_gnum = _mesh->getVertexGNum(i_part);
              CWP_g_num_t *elt_gnum = _mesh->GNumEltsGet(i_part);

              if (interf_dim == CWP_INTERFACE_SURFACE) {
                int n_edge = _mesh->getNEdge(i_part);
                int *face_edge_idx = _mesh->getFaceEdgeIndex(i_part);
                int *face_edge = _mesh->getFaceEdge(i_part);
                int *edge_vtx_idx = _mesh->getEdgeVtxIndex(i_part);
                int *edge_vtx = _mesh->getEdgeVtx(i_part);

                PDM_mesh_location_part_set_2d(_id_pdm,
                                              i_part,
                                              n_face,
                                              face_edge_idx,
                                              face_edge,
                                              elt_gnum,
                                              n_edge,
                                              edge_vtx_idx,
                                              edge_vtx,
                                              NULL,
                                              n_vtx,
                                              coords,
                                              vtx_gnum);
              }

              else if (interf_dim == CWP_INTERFACE_VOLUME) {

                int n_cell = _mesh->getNCell(i_part);
                int *cell_face_idx = _mesh->getCellFaceIndex(i_part);
                int *cell_face = _mesh->getCellFace(i_part);
                int *face_vtx_idx = _mesh->getFaceVtxIndex(i_part);
                int *face_vtx = _mesh->getFaceVtx(i_part);

                PDM_mesh_location_part_set(_id_pdm,
                                           i_part,
                                           n_cell,
                                           cell_face_idx,
                                           cell_face,
                                           elt_gnum,
                                           n_face,
                                           face_vtx_idx,
                                           face_vtx,
                                           NULL,
                                           n_vtx,
                                           coords,
                                           vtx_gnum);
              }
            }
          }

          else {
            cwipi::Coupling &cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());
            cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();

            for (int i_part = 0 ; i_part < _cplNPart ; i_part++) {

              int n_vtx = cpl_mesh->getPartNVertex(i_part);
              int n_face = cpl_mesh->getNFace(i_part);
              double *coords = cpl_mesh->getVertexCoords(i_part);
              CWP_g_num_t *vtx_gnum = cpl_mesh->getVertexGNum(i_part);
              CWP_g_num_t *elt_gnum = cpl_mesh->GNumEltsGet(i_part);

              if (interf_dim == CWP_INTERFACE_SURFACE) {
                int n_edge = cpl_mesh->getNEdge(i_part);
                int *face_edge_idx = cpl_mesh->getFaceEdgeIndex(i_part);
                int *face_edge = cpl_mesh->getFaceEdge(i_part);
                int *edge_vtx_idx = cpl_mesh->getEdgeVtxIndex(i_part);
                int *edge_vtx = cpl_mesh->getEdgeVtx(i_part);
                PDM_mesh_location_part_set_2d(_id_pdm,
                                              i_part,
                                              n_face,
                                              face_edge_idx,
                                              face_edge,
                                              elt_gnum,
                                              n_edge,
                                              edge_vtx_idx,
                                              edge_vtx,
                                              NULL,
                                              n_vtx,
                                              coords,
                                              vtx_gnum);
              }

              else if (interf_dim == CWP_INTERFACE_VOLUME) {
                int n_cell = cpl_mesh->getNCell(i_part);
                int *cell_face_idx = cpl_mesh->getCellFaceIndex(i_part);
                int *cell_face = cpl_mesh->getCellFace(i_part);
                int *face_vtx_idx = cpl_mesh->getFaceVtxIndex(i_part);
                int *face_vtx = cpl_mesh->getFaceVtx(i_part);
                PDM_mesh_location_part_set(_id_pdm,
                                           i_part,
                                           n_cell,
                                           cell_face_idx,
                                           cell_face,
                                           elt_gnum,
                                           n_face,
                                           face_vtx_idx,
                                           face_vtx,
                                           NULL,
                                           n_vtx,
                                           coords,
                                           vtx_gnum);
              }
            }
          }
        }
      }
    }
  }


  void 
  SpatialInterpLocationMeshLocation::localization_compute
  (
  ) 
  {

    if (!_coupledCodeProperties->localCodeIs()) {

      PDM_mesh_location_compute(_id_pdm);
      PDM_mesh_location_dump_times(_id_pdm);
    }

    else { 
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        printf("localization_compute - 1.2\n");
        fflush(stdout);

        PDM_mesh_location_compute(_id_pdm);
        PDM_mesh_location_dump_times(_id_pdm);
      }
    }
  }


  void 
  SpatialInterpLocationMeshLocation::localization_get
  (
  ) 
  {

    if (!_coupledCodeProperties->localCodeIs()) {

      if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {

        for (int i_part = 0 ; i_part < _nPart ; i_part++) {


          PDM_mesh_location_cell_vertex_get (_id_pdm,
                                             i_part,
                                             &(_cell_vtx_idx[i_part]),
                                             &(_cell_vtx[i_part]));

          PDM_mesh_location_points_in_elt_get (_id_pdm,
                                               i_part,
                                               0,
                                               &(_elt_pts_inside_idx[i_part]),
                                               &(_points_gnum[i_part]),
                                               &(_points_coords[i_part]),
                                               &(_points_uvw[i_part]),
                                               &(_weights_idx[i_part]),
                                               &(_weights[i_part]),
                                               &(_points_dist2[i_part]),
                                               &(_points_projected_coords[i_part]));

          int n_elt = 0;
          CWP_Interface_t interf_dim = _cpl->entitiesDimGet();

          if (interf_dim == CWP_INTERFACE_SURFACE) {
            n_elt = _mesh->getNFace(i_part);
          }
          else if (interf_dim == CWP_INTERFACE_VOLUME) {
            n_elt = _mesh->getNCell(i_part);
          }

          _n_elt_weights[i_part] = _elt_pts_inside_idx[i_part][n_elt];

        }
      }

      else {

        for (int i_part = 0; i_part < _nPart; i_part++) {

          _n_computed_tgt[i_part] = PDM_mesh_location_n_located_get (_id_pdm,
                                                                     0, 
                                                                     i_part);

          _n_uncomputed_tgt[i_part] = PDM_mesh_location_n_unlocated_get (_id_pdm,
                                                                         0, 
                                                                         i_part);

          _computed_tgt[i_part] = PDM_mesh_location_located_get (_id_pdm,
                                                                 0, 
                                                                 i_part);

          _uncomputed_tgt[i_part] = PDM_mesh_location_unlocated_get (_id_pdm,
                                                                     0, 
                                                                     i_part);

          PDM_mesh_location_point_location_get (_id_pdm,
                                                0,
                                                i_part,
                                                &(_tgt_closest_elt_gnum[i_part]),
                                                &(_tgt_distance[i_part]),
                                                &(_tgt_projected[i_part]));
          // For pdm_part1_to_selected_part2

          _elt_pts_inside_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
          _elt_pts_inside_idx[i_part][0] = 0;

        }      
      }
    }

    else {

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        SpatialInterpLocationMeshLocation *cpl_spatial_interp;

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet(); 

          cpl_spatial_interp = dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }

        else {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 

          cpl_spatial_interp = dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }

        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {

          for (int i_part = 0 ; i_part < _nPart ; i_part++) {

            PDM_mesh_location_points_in_elt_get (_id_pdm,
                                                 i_part,
                                                 0,
                                                 &(_elt_pts_inside_idx[i_part]),
                                                 &(_points_gnum[i_part]),
                                                 &(_points_coords[i_part]),
                                                 &(_points_uvw[i_part]),
                                                 &(_weights_idx[i_part]),
                                                 &(_weights[i_part]),
                                                 &(_points_dist2[i_part]),
                                                 &(_points_projected_coords[i_part]));

            int n_elt = 0;
            CWP_Interface_t interf_dim = _cpl->entitiesDimGet();

            if (interf_dim == CWP_INTERFACE_SURFACE) {
              n_elt = _mesh->getNFace(i_part);
            }
            else if (interf_dim == CWP_INTERFACE_VOLUME) {
              n_elt = _mesh->getNCell(i_part);
            }

            _n_elt_weights[i_part] = _elt_pts_inside_idx[i_part][n_elt];

          }

          for (int i_part = 0; i_part < _cplNPart; i_part++) {

            cpl_spatial_interp->_n_computed_tgt[i_part] = PDM_mesh_location_n_located_get (_id_pdm,
                                                                                           0, 
                                                                                           i_part);

            cpl_spatial_interp->_n_uncomputed_tgt[i_part] = PDM_mesh_location_n_unlocated_get (_id_pdm,
                                                                                               0, 
                                                                                               i_part);

            cpl_spatial_interp->_computed_tgt[i_part] = PDM_mesh_location_located_get (_id_pdm,
                                                                                       0, 
                                                                                       i_part);

            cpl_spatial_interp->_uncomputed_tgt[i_part] = PDM_mesh_location_unlocated_get (_id_pdm,
                                                                                           0, 
                                                                                           i_part);

            PDM_mesh_location_point_location_get (_id_pdm,
                                                  0,
                                                  i_part,
                                                  &(cpl_spatial_interp->_tgt_closest_elt_gnum[i_part]),
                                                  &(cpl_spatial_interp->_tgt_distance[i_part]),
                                                  &(cpl_spatial_interp->_tgt_projected[i_part]));
            // For pdm_part1_to_selected_part2

            cpl_spatial_interp->_elt_pts_inside_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
            cpl_spatial_interp->_elt_pts_inside_idx[i_part][0] = 0; 

          }      

        }

        else {

          for (int i_part = 0; i_part < _nPart; i_part++) {

            _n_computed_tgt[i_part] = PDM_mesh_location_n_located_get (_id_pdm,
                                                                       0, 
                                                                       i_part);

            _n_uncomputed_tgt[i_part] = PDM_mesh_location_n_unlocated_get (_id_pdm,
                                                                           0, 
                                                                           i_part);

            _computed_tgt[i_part] = PDM_mesh_location_located_get (_id_pdm,
                                                                   0, 
                                                                   i_part);

            _uncomputed_tgt[i_part] = PDM_mesh_location_unlocated_get (_id_pdm,
                                                                       0, 
                                                                       i_part);

            PDM_mesh_location_point_location_get (_id_pdm,
                                                  0,
                                                  i_part,
                                                  &(_tgt_closest_elt_gnum[i_part]),
                                                  &(_tgt_distance[i_part]),
                                                  &(_tgt_projected[i_part]));
            // For pdm_part1_to_selected_part2

            _elt_pts_inside_idx[i_part] = (int*)  malloc(sizeof(int)); // Use malloc not new [] !
            _elt_pts_inside_idx[i_part][0] = 0;


          }

          for (int i_part = 0 ; i_part < _cplNPart ; i_part++) {

            PDM_mesh_location_points_in_elt_get (_id_pdm,
                                                 i_part,
                                                 0,
                                                 &(cpl_spatial_interp->_elt_pts_inside_idx[i_part]),
                                                 &(cpl_spatial_interp->_points_gnum[i_part]),
                                                 &(cpl_spatial_interp->_points_coords[i_part]),
                                                 &(cpl_spatial_interp->_points_uvw[i_part]),
                                                 &(cpl_spatial_interp->_weights_idx[i_part]),
                                                 &(cpl_spatial_interp->_weights[i_part]),
                                                 &(cpl_spatial_interp->_points_dist2[i_part]),
                                                 &(cpl_spatial_interp->_points_projected_coords[i_part]));

            int n_elt = 0;
            CWP_Interface_t interf_dim = cpl_cpl.entitiesDimGet();

            if (interf_dim == CWP_INTERFACE_SURFACE) {
              n_elt = cpl_spatial_interp->_mesh->getNFace(i_part);
            }
            else if (interf_dim == CWP_INTERFACE_VOLUME) {
              n_elt = cpl_spatial_interp->_mesh->getNCell(i_part);
            }

            cpl_spatial_interp->_n_elt_weights[i_part] = cpl_spatial_interp->_elt_pts_inside_idx[i_part][n_elt];

          }

        }
      }
    }
  }


  void SpatialInterpLocationMeshLocation::localization_free() {

    if (!_coupledCodeProperties->localCodeIs()) {
      PDM_mesh_location_free(_id_pdm);
      _id_pdm = nullptr;
    }

    else {

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        PDM_mesh_location_free(_id_pdm);
        _id_pdm = nullptr;

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        SpatialInterpLocationMeshLocation *cpl_spatial_interp;

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet(); 

          cpl_spatial_interp = dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }

        else {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 

          cpl_spatial_interp = dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }
        cpl_spatial_interp->_id_pdm = nullptr;
      }
    }
  }
}

/**
 * \endcond
 */
