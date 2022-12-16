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
#include "pdm_gnum_location.h"

#include "spatialInterpLocation.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"

#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
extern "C" {
  void PROCF(callfortlocinterpfct, CALLFORTLOCINTERPFCT)
  (
    int    *interface_type,
    char   *code_name,
    int    *src_n_block,
    int    *src_blocks_type,
    int    *src_i_part,
    int    *src_n_vtx,
    double *src_vtx_coords,
    int    *src_vtx_global_num,
    int    *src_n_elts,
    int    *src_id_block,
    int    *src_elt_in_block,
    int    *src_elt_vtx_idx,
    int    *src_elt_vtx,
    int    *src_elts_global_num,
    int    *tgt_n_pts,
    int    *tgt_pts_elt_idx,
    double *tgt_pts_coords,
    double *tgt_pts_dist,
    double *tgt_pts_uvw,
    int    *tgt_pts_weights_idx,
    double *tgt_pts_weights,
    int    *stride,
    int    *src_field_dof_location,
    void   *src_field,
    void   *tgt_field,
    void   *ptFortranLocInterpolationFct
  );
}
#endif

/**
 * \cond
 */

namespace cwipi {
  SpatialInterpLocation::SpatialInterpLocation() = default;

  SpatialInterpLocation::~SpatialInterpLocation
  (
  )
  {

    for (int i = 0; i < _nPart; i++) {
      if (_tgt_distance[i] != NULL) {
        free (_tgt_distance[i]);
      }
      if (_tgt_projected[i] != NULL) {
        free (_tgt_projected[i]);
      }
      if (_tgt_closest_elt_gnum[i] != NULL) {
        free (_tgt_closest_elt_gnum[i]);
      }

      if (_elt_pts_inside_idx[i] != NULL) {
        free (_elt_pts_inside_idx[i]);
        free (_points_gnum[i]);
        free (_points_coords[i]);
        free (_points_dist2[i]);
        free (_points_projected_coords[i]);
      }

      if (_cell_vtx[i] != NULL) {
        free (_cell_vtx[i]);
        free (_cell_vtx_idx[i]);
      }

      if (_points_uvw[i] != NULL) {
        free (_points_uvw[i]);
      }      
    }

    delete[] _cell_vtx;
    delete[] _cell_vtx_idx;

    delete[] _tgt_distance;
    delete[] _tgt_projected;
    delete[] _tgt_closest_elt_gnum;

    delete[] _elt_pts_inside_idx;
    delete[] _points_gnum;
    delete[] _points_coords;
    delete[] _points_uvw;
    delete[] _points_dist2;
    delete[] _points_projected_coords;


    //PDM_Mesh_nodal_free (_pdm_CplNodal);
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
    // Map nodal mesh

    _pdm_CplNodal = NULL;

    // if (!_coupledCodeProperties->localCodeIs()) {

    if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
      _pdm_CplNodal = _mesh->getPdmNodalIndex();
    }
    else {
      _pdm_CplNodal = NULL; 
    }

    // }

    // else {
    //   if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
    //     if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
    //       _pdm_CplNodal = _mesh->getPdmNodalIndex();
    //     }
    //     else {
    //       cwipi::Coupling &cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());
    //       cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();
    //       _pdm_CplNodal = cpl_mesh->getPdmNodalIndex();

    //       _pdm_CplNodal = NULL; 
    //     }
    //   }
    // }

    if (_pdm_CplNodal != NULL) {
      printf("SpatialInterpLocation::init Mesh_nodal_n_blocks :%d\n", PDM_Mesh_nodal_n_blocks_get(_pdm_CplNodal));
      fflush(stdout);
    }

    //
    // Data for PDM_part_to_part_t

    if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
      for (int i_part = 0 ; i_part < _nPart ; i_part++) { 
       _src_gnum[i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
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
    _cell_vtx_idx = new int* [_nPart];
    _cell_vtx = new int* [_nPart];

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
      _cell_vtx_idx[i_part] = NULL;
      _cell_vtx[i_part] = NULL;
    }

  }


  void SpatialInterpLocation::weightsCompute() 
  {
    if (!_coupledCodeProperties->localCodeIs()) {

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

        if (_cell_vtx[i_part] != NULL) {
          free (_cell_vtx[i_part]);
          free (_cell_vtx_idx[i_part]);
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

        if (_involved_sources_tgt[i_part] != NULL) {
          free(_involved_sources_tgt[i_part]);
        }

        _n_elt_weights[i_part] = 0;
        _weights_idx[i_part] = NULL;
        _weights[i_part] = NULL;

        _n_computed_tgt[i_part] = 0;
        _computed_tgt[i_part] = NULL;

        _n_uncomputed_tgt[i_part] = 0;
        _uncomputed_tgt[i_part] = NULL;

        _n_involved_sources_tgt[i_part] = 0;
        _involved_sources_tgt[i_part] = NULL;

        _tgt_distance[i_part] = NULL;
        _tgt_projected[i_part] = NULL;
        _tgt_closest_elt_gnum[i_part] = NULL;
        _elt_pts_inside_idx[i_part] = NULL;
        _points_gnum[i_part] = NULL;
        _points_coords[i_part] = NULL;
        _points_uvw[i_part] = NULL;
        _points_dist2[i_part] = NULL;
        _points_projected_coords[i_part] = NULL;
        _cell_vtx[i_part] = NULL;
        _cell_vtx_idx[i_part] = NULL;

      }
    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

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

          if (_cell_vtx[i_part] != NULL) {
            free (_cell_vtx[i_part]);
            free (_cell_vtx_idx[i_part]);
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

          if (_involved_sources_tgt[i_part] != NULL) {
            free(_involved_sources_tgt[i_part]);
          }

          _n_elt_weights[i_part] = 0;
          _weights_idx[i_part] = NULL;
          _weights[i_part] = NULL;

          _n_computed_tgt[i_part] = 0;
          _computed_tgt[i_part] = NULL;

          _n_uncomputed_tgt[i_part] = 0;
          _uncomputed_tgt[i_part] = NULL;

          _n_involved_sources_tgt[i_part] = 0;
          _involved_sources_tgt[i_part] = NULL;

          _tgt_distance[i_part] = NULL;
          _tgt_projected[i_part] = NULL;
          _tgt_closest_elt_gnum[i_part] = NULL;
          _elt_pts_inside_idx[i_part] = NULL;
          _points_gnum[i_part] = NULL;
          _points_coords[i_part] = NULL;
          _points_uvw[i_part] = NULL;
          _points_dist2[i_part] = NULL;
          _points_projected_coords[i_part] = NULL;
          _cell_vtx[i_part] = NULL;
          _cell_vtx_idx[i_part] = NULL;

        }


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

        for (int i_part = 0; i_part < _nPart; i_part++) {
          if (cpl_spatial_interp->_tgt_distance[i_part] != NULL) {
            free (cpl_spatial_interp->_tgt_distance[i_part]);
            free (cpl_spatial_interp->_tgt_projected[i_part]);
            free (cpl_spatial_interp->_tgt_closest_elt_gnum[i_part]);
          }
          if (cpl_spatial_interp->_elt_pts_inside_idx[i_part] != NULL) {
            free (cpl_spatial_interp->_elt_pts_inside_idx[i_part]);
            free (cpl_spatial_interp->_points_gnum[i_part]);
            free (cpl_spatial_interp->_points_coords[i_part]);
            free (cpl_spatial_interp->_points_uvw[i_part]);
            free (cpl_spatial_interp->_points_dist2[i_part]);
            free (cpl_spatial_interp->_points_projected_coords[i_part]);
          }

          if (cpl_spatial_interp->_cell_vtx[i_part] != NULL) {
            free (cpl_spatial_interp->_cell_vtx[i_part]);
            free (cpl_spatial_interp->_cell_vtx_idx[i_part]);
          }

          if (cpl_spatial_interp->_weights_idx[i_part] != NULL) {
            free (cpl_spatial_interp->_weights_idx[i_part]);
            free (cpl_spatial_interp->_weights[i_part]);
          }

          if (cpl_spatial_interp->_computed_tgt[i_part] != NULL) {
            free (cpl_spatial_interp->_computed_tgt[i_part]);   
          }

          if (cpl_spatial_interp->_uncomputed_tgt[i_part] != NULL) {
            free (cpl_spatial_interp->_uncomputed_tgt[i_part]);
          }

          if (cpl_spatial_interp->_involved_sources_tgt[i_part] != NULL) {
            free(cpl_spatial_interp->_involved_sources_tgt[i_part]);
          }

          cpl_spatial_interp->_n_elt_weights[i_part] = 0;
          cpl_spatial_interp->_weights_idx[i_part] = NULL;
          cpl_spatial_interp->_weights[i_part] = NULL;

          cpl_spatial_interp->_n_computed_tgt[i_part] = 0;
          cpl_spatial_interp->_computed_tgt[i_part] = NULL;

          cpl_spatial_interp->_n_uncomputed_tgt[i_part] = 0;
          cpl_spatial_interp->_uncomputed_tgt[i_part] = NULL;

          cpl_spatial_interp->_n_involved_sources_tgt[i_part] = 0;
          cpl_spatial_interp->_involved_sources_tgt[i_part] = NULL;

          cpl_spatial_interp->_tgt_distance[i_part] = NULL;
          cpl_spatial_interp->_tgt_projected[i_part] = NULL;
          cpl_spatial_interp->_tgt_closest_elt_gnum[i_part] = NULL;
          cpl_spatial_interp->_elt_pts_inside_idx[i_part] = NULL;
          cpl_spatial_interp->_points_gnum[i_part] = NULL;
          cpl_spatial_interp->_points_coords[i_part] = NULL;
          cpl_spatial_interp->_points_uvw[i_part] = NULL;
          cpl_spatial_interp->_points_dist2[i_part] = NULL;
          cpl_spatial_interp->_points_projected_coords[i_part] = NULL;
          cpl_spatial_interp->_cell_vtx[i_part] = NULL;
          cpl_spatial_interp->_cell_vtx_idx[i_part] = NULL;

        }
      }
    }

    localization_init();

    localization_points_cloud_setting();
    
    localization_surface_setting();

    localization_compute();

    if (_ptsp != nullptr) {
      if (!_coupledCodeProperties->localCodeIs()) {
        PDM_part_to_part_free (_ptsp);
        _ptsp = nullptr;
      }
      else {

        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
          PDM_part_to_part_free (_ptsp);
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

    localization_get();

    localization_free();



    if (!_coupledCodeProperties->localCodeIs()) {
     // printf("_src_gnum %d :", _src_n_gnum[0]);
     // for (int i = 0; i < _src_n_gnum[0]; i++) {
     //   printf(" %ld", _src_gnum[0][i]);
     // }
     // printf("\n");

     // printf("_tgt_gnum %d :", _tgt_n_gnum[0]);
     // for (int i = 0; i < _tgt_n_gnum[0]; i++) {
     //   printf(" %ld", _tgt_gnum[0][i]);
     // }
     // printf("\n");

     // for (int i = 0; i < _src_n_gnum[0]; i++) {
     //   printf("_elt_pts_inside_idx %d %d: ", _elt_pts_inside_idx[0][i], _elt_pts_inside_idx[0][i+1]);
     //   for (int j = _elt_pts_inside_idx[0][i]; j < _elt_pts_inside_idx[0][i+1]; j++) {
     //     printf(" %ld", _points_gnum[0][j]);
     //   }
     //   printf("\n");
     // }
     // printf("\n");
     // fflush(stdout);
      if (_ptsp == NULL) {
        _ptsp = PDM_part_to_part_create ((const PDM_g_num_t **)_src_gnum,
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
        printf("ptp from mesh_location OK\n");
        fflush(stdout);
      }
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

        if (_ptsp == NULL) {
          if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
            _ptsp = PDM_part_to_part_create ((const PDM_g_num_t **)_src_gnum,
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
            _ptsp = PDM_part_to_part_create ((const PDM_g_num_t **)cpl_spatial_interp->_src_gnum,
                                             (const int *)cpl_spatial_interp->_src_n_gnum,
                                             _cplNPart,
                                             (const PDM_g_num_t **) _tgt_gnum,
                                             (const int *) _tgt_n_gnum,
                                             _nPart,
                                             (const int **)cpl_spatial_interp->_elt_pts_inside_idx,
                                             (const PDM_g_num_t **)cpl_spatial_interp->_points_gnum,
                                             _pdmCplComm);

          }
        }
        else {
          printf("ptp from mesh_location OK\n");
          fflush(stdout);
        }

        cpl_spatial_interp->_ptsp = _ptsp;
      } 
    }


  }

  void SpatialInterpLocation::interpolate (Field *referenceField, double **buffer) 
  {
    int nComponent = referenceField->nComponentGet();
    CWP_Dof_location_t referenceFieldType = referenceField->locationGet();
    CWP_Interpolation_t interpolationType = referenceField->interpolationTypeGet();
    const CWP_Field_storage_t storage = referenceField->storageTypeGet();

    // int  _nPart; !< Mesh partition number                                                    

    if (interpolationType == CWP_INTERPOLATION_USER) {

      CWP_Interp_from_location_t interpolationFunction = referenceField->interpolationFunctionGet();
      void *interpolationFunction_f = referenceField->fortranInterpolationFunctionGet();

      if (interpolationFunction != NULL) {

        for (int i_part = 0 ; i_part < _nPart ; i_part++) {
          double      *part_points_coords           = _points_coords[i_part];
          double      *part_points_uvw              = _points_uvw[i_part];
          double      *part_points_dist2            = _points_dist2[i_part];

          int         *part_weights_idx             = _weights_idx[i_part];
          double      *part_weights                 = _weights[i_part];

         (*interpolationFunction)(_cpl->entitiesDimGet(),
                                  _localCodeProperties->nameGet().c_str(),
                                  _mesh->nBlockGet(),
                                  _mesh->blocksTypeGet(),
                                  i_part,
                                  _mesh->getPartNVertex(i_part),
                                  _mesh->getVertexCoords(i_part),
                                  _mesh->getVertexGNum(i_part),
                                  _mesh->getPartNElts(i_part),
                                  _mesh->eltIdBlockGet(i_part),
                                  _mesh->eltInBlockGet(i_part),
                                  _cell_vtx_idx[i_part],
                                  _cell_vtx[i_part],
                                  _mesh->GNumEltsGet(i_part),
                                  _n_elt_weights[i_part],
                                  _elt_pts_inside_idx[i_part],
                                  part_points_coords,
                                  part_points_dist2,
                                  part_points_uvw,
                                  part_weights_idx,
                                  part_weights,
                                  nComponent,
                                  referenceFieldType,
                                  (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_SOURCE),
                                  buffer
                                  );
        }
    
      }

      else if (interpolationFunction_f != NULL) {


      }

    }

    else { 
      if (referenceFieldType == CWP_DOF_LOCATION_CELL_CENTER) {

        for (int i_part = 0; i_part < _nPart; i_part++) {
          int    *part_elt_pts_inside_idx = _elt_pts_inside_idx[i_part];
          double *referenceData           = (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_SOURCE);

          int    part_n_elt               = _mesh->getPartNElts(i_part);

          double *local_buffer = (double *) *buffer;

          int ival = 0;
          if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
            for (int i = 0; i < part_n_elt; i++) {
              for (int j = part_elt_pts_inside_idx[i]; j < part_elt_pts_inside_idx[i+1]; j++) {
                for (int k1 = 0; k1 < nComponent; k1++) {
                  ival = k1*part_elt_pts_inside_idx[part_n_elt] + j;
                  local_buffer[ival] = referenceData[part_n_elt*k1 + i];
                }
              }
            }
          }
          else { // if (storage == CWP_FIELD_STORAGE_INTERLACED) {
            for (int i = 0; i < part_n_elt; i++) {
              for (int j = part_elt_pts_inside_idx[i]; j < part_elt_pts_inside_idx[i+1]; j++) {
                for (int k1 = 0; k1 < nComponent; k1++) {
                  local_buffer[ival++] = referenceData[i*nComponent + k1];
                }
              }
            }
          }

        }
      }
      
      else if (referenceFieldType == CWP_DOF_LOCATION_NODE || referenceFieldType == CWP_DOF_LOCATION_USER) {

        for (int i_part = 0; i_part < _nPart; i_part++) {
          int         *part_elt_pts_inside_idx = _elt_pts_inside_idx[i_part];
          double      *referenceData           = (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_SOURCE);

          int         *part_weights_idx        = _weights_idx[i_part];
          double      *part_weights            = _weights[i_part];

          int          part_n_elt              = _mesh->getPartNElts(i_part);
          int          part_n_vtx              = _mesh->getPartNVertex(i_part);

          int         *connec_idx2             = _cell_vtx_idx[i_part];
          int         *connec2                 = _cell_vtx[i_part];

          int         *connec_idx = connec_idx2; 
          int         *connec = connec2;
          
          double *local_buffer = (double *) *buffer;

          int ival = 0;

          if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
            for (int i = 0; i < part_n_elt; i++) {
              for (int j = part_elt_pts_inside_idx[i]; j < part_elt_pts_inside_idx[i+1]; j++) {
                for (int k1 = 0; k1 < nComponent; k1++) {
                  ival = k1*part_elt_pts_inside_idx[part_n_elt] + j;
                  local_buffer[ival] = 0;
                  int k2 = connec_idx[i];
                  assert(connec_idx[i+1] - connec_idx[i] == part_weights_idx[j + 1] - part_weights_idx[j]);
                  for (int k = part_weights_idx[j]; k < part_weights_idx[j+1]; k++) {
                    int isom = connec[k2++] - 1;
                    local_buffer[ival] += part_weights[k] * referenceData[part_n_vtx*k1 + isom];
                  }
                }
              }
            }
          }
          else { // if (storage == CWP_FIELD_STORAGE_INTERLACED) {
            for (int i = 0; i < part_n_elt; i++) {
              for (int j = part_elt_pts_inside_idx[i]; j < part_elt_pts_inside_idx[i+1]; j++) {
                for (int k1 = 0; k1 < nComponent; k1++) {
                  local_buffer[ival] = 0;
                  int k2 = connec_idx[i];
                  assert(connec_idx[i+1] - connec_idx[i] == part_weights_idx[j + 1] - part_weights_idx[j]);
                  for (int k = part_weights_idx[j]; k < part_weights_idx[j+1]; k++) {
                    int isom = connec[k2++] - 1;
                    local_buffer[ival] += part_weights[k] * referenceData[isom*nComponent+k1];
                  }
                  ival++;
                }
              }
            }
          }
        }
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



}

/**
 * \endcond
 */
