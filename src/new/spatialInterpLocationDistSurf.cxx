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

#include <pdm_dist_cloud_surf.h>

#include "spatialInterpLocationDistSurf.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"

/**
 * \cond
 */

namespace cwipi {
    void
    SpatialInterpLocationDistSurf::localization_init
    (
    )

    {
        if (!_coupledCodeProperties->localCodeIs()) {
            _id_pdm = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

            if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
                PDM_dist_cloud_surf_surf_mesh_global_data_set(_id_pdm, _cplNPart); // TODO What should be in args 2 and 3 ?
                PDM_dist_cloud_surf_n_part_cloud_set(_id_pdm, 0, _nPart);
            }
            else {
                PDM_dist_cloud_surf_surf_mesh_global_data_set(_id_pdm, _nPart); // TODO What should be in args 2 and 3 ?
                PDM_dist_cloud_surf_n_part_cloud_set(_id_pdm, 0, _cplNPart);
            }
        }
        else {
            if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
                _id_pdm = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

                SpatialInterpLocationDistSurf *cpl_spatial_interp;

                cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

                if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
                    std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();

                    cpl_spatial_interp =
                            dynamic_cast <SpatialInterpLocationDistSurf *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
                }
                else {
                    std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();

                    cpl_spatial_interp =
                            dynamic_cast <SpatialInterpLocationDistSurf *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
                }

                cpl_spatial_interp->_id_pdm = _id_pdm;

                if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
                    PDM_dist_cloud_surf_surf_mesh_global_data_set(_id_pdm, _cplNPart); // TODO What should be in args 2 and 3 ?
                    PDM_dist_cloud_surf_n_part_cloud_set(_id_pdm, 0, _nPart);
                }
                else {
                    PDM_dist_cloud_surf_surf_mesh_global_data_set(_id_pdm, _nPart); // TODO What should be in args 2 and 3 ?
                    PDM_dist_cloud_surf_n_part_cloud_set(_id_pdm, 0, _cplNPart);
                }
            }
        }
    }

    void SpatialInterpLocationDistSurf::localization_points_cloud_setting
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

                    PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, part_n, (double *) part_coord, (PDM_g_num_t*) part_gnum);
                }
            }
            else {
                for (int i_part = 0 ; i_part < _nPart ; i_part++) {
                    PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, 0, NULL, NULL);
                }
            }
        }
        else {
            if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
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

                        PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, part_n, (double *) part_coord, (PDM_g_num_t*) part_gnum);
                    }
                }
                else {
                    cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

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

                        PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, part_n, (double *) part_coord, (PDM_g_num_t*) part_gnum);
                    }
                }
            }
        }

        // _id_pdm = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        // PDM_dist_cloud_surf_n_part_cloud_set(_id_pdm, 0, _nPart);
        // PDM_dist_cloud_surf_surf_mesh_global_data_set(_id_pdm, _n_g_elt_cpl_over_part, _n_g_vtx_cpl_over_part, _nPart_cpl);

        // for (int i_part = 0 ; i_part < _nPart ; i_part++) PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, _n_target[i_part], _coords_target[i_part], _gnum_target[i_part]);

        // if (!_both_codes_are_local) {
        //     for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
        //         int *connec_idx_null = (int *) malloc(sizeof(int));
        //         connec_idx_null[0] = 0;
        //         PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, 0, connec_idx_null, NULL, NULL, 0, NULL, NULL);
        //     }
        // }
        // else {
        //     for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
        //         Mesh *mesh_cpl = _spatial_interp_cpl->_mesh;
        //         int *connecIdx_cpl = mesh_cpl->connecIdxGet(i_part);
        //         int *connec_cpl = mesh_cpl->connecGet(i_part);

        //         int n_vtx_cpl = mesh_cpl->getPartNVertex(i_part);
        //         int n_elts_cpl = mesh_cpl->getPartNElts(i_part);
        //         double *coords_cpl = mesh_cpl->getVertexCoords(i_part);
        //         CWP_g_num_t *gnum_vtx_cpl = mesh_cpl->getVertexGNum(i_part);
        //         CWP_g_num_t *gnum_elt_cpl = mesh_cpl->GNumEltsGet(i_part);

        //         PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, n_elts_cpl, connecIdx_cpl, connec_cpl, gnum_elt_cpl, n_vtx_cpl, coords_cpl, gnum_vtx_cpl);
        //     }
        // }
    }

    void SpatialInterpLocationDistSurf::localization_surface_setting
    (
    )
    {
        // Ne fonctionne pas. Il faudra reconstruire connec et connecIdx (_cell_vtx et cell_vtx_id) commme pour la localisation (PDM_Mesh_location) 

        if (!_coupledCodeProperties->localCodeIs()) {
            if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
                for (int i_part = 0 ; i_part < _nPart ; i_part++) {
                    int n_vtx = _mesh->getPartNVertex(i_part);
                    int n_elts = _mesh->getPartNElts(i_part);
                    double *coords = _mesh->getVertexCoords(i_part);
                    CWP_g_num_t *vtx_gnum = _mesh->getVertexGNum(i_part);
                    CWP_g_num_t *elt_gnum = _mesh->GNumEltsGet(i_part);
                    int *connecIdx = NULL;
                    int *connec = NULL;
                    // int *connecIdx = _mesh->connecIdxGet(i_part);
                    // int *connec = _mesh->connecGet(i_part);

                    PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, n_elts, connecIdx, connec, elt_gnum, n_vtx, coords, vtx_gnum);
                }
            }
            else {
                for (int i_part = 0 ; i_part < _nPart ; i_part++) {
                    PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, 0, fake_idx, NULL, NULL, 0, NULL, NULL);
                }
            }
        }
        else {
            if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
                if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
                    for (int i_part = 0 ; i_part < _nPart ; i_part++) {
                        int n_vtx = _mesh->getPartNVertex(i_part);
                        int n_elts = _mesh->getPartNElts(i_part);
                        double *coords = _mesh->getVertexCoords(i_part);
                        CWP_g_num_t *vtx_gnum = _mesh->getVertexGNum(i_part);
                        CWP_g_num_t *elt_gnum = _mesh->GNumEltsGet(i_part);
                        int *connecIdx = NULL;
                        int *connec = NULL;
                        // int *connecIdx = _mesh->connecIdxGet(i_part);
                        // int *connec = _mesh->connecGet(i_part);

                        PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, n_elts, connecIdx, connec, elt_gnum, n_vtx, coords, vtx_gnum);
                    }
                }
                else {
                    cwipi::Coupling &cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

                    cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();

                    for (int i_part = 0 ; i_part < _cplNPart ; i_part++) {
                        int n_vtx = cpl_mesh->getPartNVertex(i_part);
                        int n_elts = cpl_mesh->getPartNElts(i_part);
                        double *coords = cpl_mesh->getVertexCoords(i_part);
                        CWP_g_num_t *vtx_gnum = cpl_mesh->getVertexGNum(i_part);
                        CWP_g_num_t *elt_gnum = cpl_mesh->GNumEltsGet(i_part);
                        int *connecIdx = NULL;
                        int *connec = NULL;
                        // int *connecIdx = cpl_mesh->connecIdxGet(i_part);
                        // int *connec = cpl_mesh->connecGet(i_part);

                        PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, n_elts, connecIdx, connec, elt_gnum, n_vtx, coords, vtx_gnum);
                    }
                }
            }
        }

        // _id_pdm = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdmCplComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        // PDM_dist_cloud_surf_n_part_cloud_set(_id_pdm, 0, _nPart_cpl);
        // PDM_dist_cloud_surf_surf_mesh_global_data_set(_id_pdm, _n_g_elt_over_part, _n_g_vtx_over_part, _nPart);

        // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
        //     int *connecIdx = _mesh->connecIdxGet(i_part);
        //     int *connec = _mesh->connecGet(i_part);

        //     int n_vtx = _mesh->getPartNVertex(i_part);
        //     int n_elts = _mesh->getPartNElts(i_part);
        //     double *coords = _mesh->getVertexCoords(i_part);
        //     CWP_g_num_t *gnum_vtx = _mesh->getVertexGNum(i_part);
        //     CWP_g_num_t *gnum_elt = _mesh->GNumEltsGet(i_part);

        //     PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, n_elts, connecIdx, connec, gnum_elt, n_vtx, coords, gnum_vtx);
        // }

        // if (!_both_codes_are_local)
        //     for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, 0, NULL, NULL);
        // else
        //     for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, _spatial_interp_cpl->_n_target[i_part], _spatial_interp_cpl->_coords_target[i_part], _spatial_interp_cpl->_gnum_target[i_part]);
    }

    void SpatialInterpLocationDistSurf::localization_compute
    (
    )
    {
        if (!_coupledCodeProperties->localCodeIs()) {
            PDM_dist_cloud_surf_compute(_id_pdm);
            PDM_dist_cloud_surf_dump_times(_id_pdm);
        }
        else {
            if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
                PDM_dist_cloud_surf_compute(_id_pdm);
                PDM_dist_cloud_surf_dump_times(_id_pdm);
            }
        }
    }

//    void SpatialInterpLocationDistSurf::localization_get_cpl() {
        // _spatial_interp_cpl->_distance = (double **) malloc(sizeof(double *) * _nPart_cpl);
        // _spatial_interp_cpl->_projected = (double **) malloc(sizeof(double *) * _nPart_cpl);
        // _spatial_interp_cpl->_closest_elt_gnum = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * _nPart_cpl);

        // for (int i_part = 0 ; i_part < _nPart_cpl ; i_part++) {
        //     int n_target_cpl = _spatial_interp_cpl->_n_target[i_part];
        //     PDM_dist_cloud_surf_get(_id_pdm, 0, i_part,
        //                             &(_spatial_interp_cpl->_distance[i_part]),
        //                             &(_spatial_interp_cpl->_projected[i_part]),
        //                             &(_spatial_interp_cpl->_closest_elt_gnum[i_part]));

        //     for (int i = 0 ; i < n_target_cpl ; i++) {
        //         if (_spatial_interp_cpl->_closest_elt_gnum[i_part][i] > CWP_g_num_t(_n_g_elt_over_part)
        //             || _spatial_interp_cpl->_closest_elt_gnum[i_part][i] < CWP_g_num_t(1)
        //             || _spatial_interp_cpl->_distance[i_part][i] > 0.01) {
        //             _spatial_interp_cpl->_closest_elt_gnum[i_part][i] = CWP_g_num_t(1);
        //             _spatial_interp_cpl->_distance[i_part][i] = INFINITY;
        //         }
        //     }
        // }
//    }

    void SpatialInterpLocationDistSurf::localization_get
    (
    )
    {
        if (!_coupledCodeProperties->localCodeIs()) {
            if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
                // TODO I don t find an equivalent for PDM_mesh_location_points_in_elt_get
                //  which would define _elt_pts_inside_idx, _points_gnum and other fields marked as source in SpatialInterpLocation
            }
            else {
                for (int i_part = 0 ; i_part < _nPart ; i_part++) {
                    // The method does not generate orphans
                    // TODO There are no orphans in the method but I don t know how to fill these
                    _n_computed_tgt[i_part] = _mesh->getPartNVertex(i_part); // ?

                    _n_uncomputed_tgt[i_part] = 0;

                    _computed_tgt[i_part] = NULL;

                    _uncomputed_tgt[i_part] = NULL;

                    PDM_dist_cloud_surf_get(_id_pdm, 0, i_part, &(_tgt_distance[i_part]), &(_tgt_projected[i_part]), &(_tgt_closest_elt_gnum[i_part]));

                    // For pdm_part1_to_selected_part2
                    _elt_pts_inside_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
                    _elt_pts_inside_idx[i_part][0] = 0;
                }
            }
        }
        else {
            if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
                cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

                SpatialInterpLocationDistSurf *cpl_spatial_interp;

                cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();

                if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
                    std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();

                    cpl_spatial_interp = dynamic_cast <SpatialInterpLocationDistSurf *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
                }
                else {
                    std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();

                    cpl_spatial_interp = dynamic_cast <SpatialInterpLocationDistSurf *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
                }

                if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
                    // TODO I don t find an equivalent for PDM_mesh_location_points_in_elt_get
                    //  which would define _elt_pts_inside_idx, _points_gnum and other fields marked as source in SpatialInterpLocation

                    for (int i_part = 0; i_part < _cplNPart; i_part++) {
                        cpl_spatial_interp->_n_computed_tgt[i_part] = cpl_mesh->getPartNVertex(i_part); // ?

                        cpl_spatial_interp->_n_uncomputed_tgt[i_part] = 0;

                        cpl_spatial_interp->_computed_tgt[i_part] = NULL;

                        cpl_spatial_interp->_uncomputed_tgt[i_part] = NULL;

                        PDM_dist_cloud_surf_get (_id_pdm, 0, i_part, &(cpl_spatial_interp->_tgt_distance[i_part]), &(cpl_spatial_interp->_tgt_projected[i_part]), &(cpl_spatial_interp->_tgt_closest_elt_gnum[i_part]));

                        // For pdm_part1_to_selected_part2
                        cpl_spatial_interp->_elt_pts_inside_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
                        cpl_spatial_interp->_elt_pts_inside_idx[i_part][0] = 0;
                    }
                }
                else {
                    for (int i_part = 0; i_part < _nPart; i_part++) {
                        // The method does not generate orphans
                        // TODO There are no orphans in the method but I don t know how to fill these
                        _n_computed_tgt[i_part] = _mesh->getPartNVertex(i_part); // ?

                        _n_uncomputed_tgt[i_part] = 0;

                        _computed_tgt[i_part] = NULL;

                        _uncomputed_tgt[i_part] = NULL;

                        PDM_dist_cloud_surf_get(_id_pdm, 0, i_part, &(_tgt_distance[i_part]), &(_tgt_projected[i_part]), &(_tgt_closest_elt_gnum[i_part]));

                        // For pdm_part1_to_selected_part2
                        _elt_pts_inside_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
                        _elt_pts_inside_idx[i_part][0] = 0;
                    }

                    // TODO I don t find an equivalent for PDM_mesh_location_points_in_elt_get
                    //  which would define _elt_pts_inside_idx, _points_gnum and other fields marked as source in SpatialInterpLocation
                }
            }
        }

        // _distance = (double **) malloc(sizeof(double *) * _nPart);
        // _projected = (double **) malloc(sizeof(double *) * _nPart);
        // _closest_elt_gnum = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * _nPart);

        // for (int i_part = 0 ; i_part < _nPart ; i_part++) {
        //     PDM_dist_cloud_surf_get(_id_pdm, 0, i_part, &(_distance[i_part]), &(_projected[i_part]), &(_closest_elt_gnum[i_part]));

        //     for (int i = 0 ; i < _n_target[i_part] ; i++) {
        //         if (_closest_elt_gnum[i_part][i] > CWP_g_num_t(_n_g_elt_cpl_over_part) ||
        //             _closest_elt_gnum[i_part][i] < CWP_g_num_t(1) ||
        //             _distance[i_part][i] > 0.1) {
        //             _closest_elt_gnum[i_part][i] = CWP_g_num_t(1);
        //             _distance[i_part][i] = INFINITY;
        //         }
        //     }
        // }
    }

    void SpatialInterpLocationDistSurf::localization_free() {
        if (!_coupledCodeProperties->localCodeIs()) {
            PDM_dist_cloud_surf_free(_id_pdm);
            _id_pdm = nullptr;
        }
        else {
            if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
                PDM_dist_cloud_surf_free(_id_pdm);
                _id_pdm = nullptr;

                cwipi::Coupling &cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

                SpatialInterpLocationDistSurf *cpl_spatial_interp;

                if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
                    std::map<std::pair<CWP_Dof_location_t, CWP_Dof_location_t>, SpatialInterp *> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();

                    cpl_spatial_interp = dynamic_cast <SpatialInterpLocationDistSurf *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
                }
                else {
                    std::map<std::pair<CWP_Dof_location_t, CWP_Dof_location_t>, SpatialInterp *> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();

                    cpl_spatial_interp = dynamic_cast <SpatialInterpLocationDistSurf *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
                }
                cpl_spatial_interp->_id_pdm = nullptr;
            }
        }
    }
}

/**
 * \endcond
 */
