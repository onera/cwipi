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

#include <spatialInterpLocationDistSurf.hxx>

/**
 * \cond
 */

namespace cwipi {
    void SpatialInterpLocationDistSurf::localization_points_cloud_setting() {
        _id_pdm = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        PDM_dist_cloud_surf_n_part_cloud_set(_id_pdm, 0, _nb_part);
        PDM_dist_cloud_surf_surf_mesh_global_data_set(_id_pdm, _n_g_elt_cpl_over_part, _n_g_vtx_cpl_over_part, _nb_part_cpl);

        for (int i_part = 0 ; i_part < _nb_part ; i_part++) PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, _n_target[i_part], _coords_target[i_part], _gnum_target[i_part]);

        if (!_both_codes_are_local) {
            for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
                int *connec_idx_null = (int *) malloc(sizeof(int));
                connec_idx_null[0] = 0;
                PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, 0, connec_idx_null, NULL, NULL, 0, NULL, NULL);
            }
        }
        else {
            for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
                Mesh *mesh_cpl = _spatial_interp_cpl->_mesh;
                int *connecIdx_cpl = mesh_cpl->connecIdxGet(i_part);
                int *connec_cpl = mesh_cpl->connecGet(i_part);

                int n_vtx_cpl = mesh_cpl->getPartNVertex(i_part);
                int n_elts_cpl = mesh_cpl->getPartNElts(i_part);
                double *coords_cpl = mesh_cpl->getVertexCoords(i_part);
                CWP_g_num_t *gnum_vtx_cpl = mesh_cpl->getVertexGNum(i_part);
                CWP_g_num_t *gnum_elt_cpl = mesh_cpl->GNumEltsGet(i_part);

                PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, n_elts_cpl, connecIdx_cpl, connec_cpl, gnum_elt_cpl, n_vtx_cpl, coords_cpl, gnum_vtx_cpl);
            }
        }
    }

    void SpatialInterpLocationDistSurf::localization_null_setting_send() {
        _id_pdm = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);
        PDM_dist_cloud_surf_n_part_cloud_set(_id_pdm, 0, _nb_part_cpl);
        PDM_dist_cloud_surf_surf_mesh_global_data_set(_id_pdm, _n_g_elt_over_part, _n_g_vtx_over_part, _nb_part);

        for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, 0, NULL, NULL);
        for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
            int *connec_idx_null = (int *) malloc(sizeof(int));
            connec_idx_null[0] = 0;
            PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, 0, connec_idx_null, NULL, NULL, 0, NULL, NULL);
        }
    }

    void SpatialInterpLocationDistSurf::localization_null_setting_recv() {
        _id_pdm = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        PDM_dist_cloud_surf_n_part_cloud_set(_id_pdm, 0, _nb_part);
        PDM_dist_cloud_surf_surf_mesh_global_data_set(_id_pdm, _n_g_elt_cpl_over_part, _n_g_vtx_cpl_over_part, _nb_part_cpl);

        for (int i_part = 0 ; i_part < _nb_part ; i_part++) PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, 0, NULL, NULL);
        for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
            int *connec_idx_null = (int *) malloc(sizeof(int));
            connec_idx_null[0] = 0;
            PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, 0, connec_idx_null, NULL, NULL, 0, NULL, NULL);
        }
    }

    void SpatialInterpLocationDistSurf::localization_surface_setting() {
        _id_pdm = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        PDM_dist_cloud_surf_n_part_cloud_set(_id_pdm, 0, _nb_part_cpl);
        PDM_dist_cloud_surf_surf_mesh_global_data_set(_id_pdm, _n_g_elt_over_part, _n_g_vtx_over_part, _nb_part);

        for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
            int *connecIdx = _mesh->connecIdxGet(i_part);
            int *connec = _mesh->connecGet(i_part);

            int n_vtx = _mesh->getPartNVertex(i_part);
            int n_elts = _mesh->getPartNElts(i_part);
            double *coords = _mesh->getVertexCoords(i_part);
            CWP_g_num_t *gnum_vtx = _mesh->getVertexGNum(i_part);
            CWP_g_num_t *gnum_elt = _mesh->GNumEltsGet(i_part);

            PDM_dist_cloud_surf_surf_mesh_part_set(_id_pdm, i_part, n_elts, connecIdx, connec, gnum_elt, n_vtx, coords, gnum_vtx);
        }

        if (!_both_codes_are_local)
            for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, 0, NULL, NULL);
        else
            for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) PDM_dist_cloud_surf_cloud_set(_id_pdm, 0, i_part, _spatial_interp_cpl->_n_target[i_part], _spatial_interp_cpl->_coords_target[i_part], _spatial_interp_cpl->_gnum_target[i_part]);
    }

    void SpatialInterpLocationDistSurf::localization_compute() {
        PDM_dist_cloud_surf_compute(_id_pdm);
        PDM_dist_cloud_surf_dump_times(_id_pdm);
    }

    void SpatialInterpLocationDistSurf::localization_get_cpl() {
        _spatial_interp_cpl->_distance = (double **) malloc(sizeof(double *) * _nb_part_cpl);
        _spatial_interp_cpl->_projected = (double **) malloc(sizeof(double *) * _nb_part_cpl);
        _spatial_interp_cpl->_closest_elt_gnum = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * _nb_part_cpl);

        for (int i_part = 0 ; i_part < _nb_part_cpl ; i_part++) {
            int n_target_cpl = _spatial_interp_cpl->_n_target[i_part];
            PDM_dist_cloud_surf_get(_id_pdm, 0, i_part,
                                    &(_spatial_interp_cpl->_distance[i_part]),
                                    &(_spatial_interp_cpl->_projected[i_part]),
                                    &(_spatial_interp_cpl->_closest_elt_gnum[i_part]));

            for (int i = 0 ; i < n_target_cpl ; i++) {
                if (_spatial_interp_cpl->_closest_elt_gnum[i_part][i] > CWP_g_num_t(_n_g_elt_over_part)
                    || _spatial_interp_cpl->_closest_elt_gnum[i_part][i] < CWP_g_num_t(1)
                    || _spatial_interp_cpl->_distance[i_part][i] > 0.01) {
                    _spatial_interp_cpl->_closest_elt_gnum[i_part][i] = CWP_g_num_t(1);
                    _spatial_interp_cpl->_distance[i_part][i] = INFINITY;
                }
            }
        }
    }

    void SpatialInterpLocationDistSurf::localization_get() {
        _distance = (double **) malloc(sizeof(double *) * _nb_part);
        _projected = (double **) malloc(sizeof(double *) * _nb_part);
        _closest_elt_gnum = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * _nb_part);

        for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
            PDM_dist_cloud_surf_get(_id_pdm, 0, i_part, &(_distance[i_part]), &(_projected[i_part]), &(_closest_elt_gnum[i_part]));

            for (int i = 0 ; i < _n_target[i_part] ; i++) {
                if (_closest_elt_gnum[i_part][i] > CWP_g_num_t(_n_g_elt_cpl_over_part) ||
                    _closest_elt_gnum[i_part][i] < CWP_g_num_t(1) ||
                    _distance[i_part][i] > 0.1) {
                    _closest_elt_gnum[i_part][i] = CWP_g_num_t(1);
                    _distance[i_part][i] = INFINITY;
                }
            }
        }
    }

    void SpatialInterpLocationDistSurf::localization_free() {
        PDM_dist_cloud_surf_free(_id_pdm);
    }
}

/**
 * \endcond
 */
