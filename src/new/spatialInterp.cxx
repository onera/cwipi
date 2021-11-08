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
#include <map>
#include "spatialInterp.hxx"
#include "spatialInterp_i.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include <mpi.h>

#include "pdm_mesh_nodal.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_gnum.h"
#include "pdm_timer.h"
#include "pdm_gnum_location.h"
#include "pdm_binary_search.h"
#include "bftc_error.h"
#include "bftc_printf.h"
#include "cwp.h"
#include <limits>
#include <vector>
#include <algorithm>
#include <cmath>


/**
 * \cond
 */

namespace cwipi {

  SpatialInterp::SpatialInterp()
  {
  }


  SpatialInterp::~SpatialInterp()
  {
    if (_ptsp != nullptr) {
      PDM_part1_to_selected_part2_free (_ptsp);
      _ptsp = nullptr;
    }

    delete[] _n_elt_weights;
    delete[] _weights_idx;
    delete[] _weights;

    delete[] _n_computed_tgt;
    delete[] _computed_tgt;

    delete[] _n_uncomputed_tgt;
    delete[] _uncomputed_tgt;

    delete[] _src_n_gnum;
    delete[] _tgt_n_gnum;
    delete[] _src_gnum;
    delete[] _tgt_gnum;
  }


  void 
  SpatialInterp::init(
    Coupling                  *coupling, 
    CWP_Dof_location_t         localCodeDofLocation,
    CWP_Dof_location_t         cplCodeDofLocation,
    SpatialInterpExchDirection exchDirection 
  )
  {
    _cpl                    = coupling;
    _visu                   = coupling->visuGet();
    _mesh                   = coupling->meshGet();

    _localCodeDofLocation   = localCodeDofLocation;
    _coupledCodeDofLocation = cplCodeDofLocation;

    _localCodeProperties    = _cpl->localCodePropertiesGet();
    _coupledCodeProperties  = _cpl->coupledCodePropertiesGet();

    _cplComm                = _cpl->communicationGet()->cplCommGet();
    _pdmCplComm             = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_cplComm));

    _unionComm              = _cpl->communicationGet()->unionCommGet();
    _pdmUnionComm           = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_unionComm));

    _nPart                  = _mesh -> getNPart();

    _rootRankUnionComm      = _cpl->communicationGet()->unionCommLocCodeRootRanksGet();
    _cplRootRankUnionComm   = _cpl->communicationGet()->unionCommCplCodeRootRanksGet();

    _rootRankCplComm        = _cpl->communicationGet()->cplCommLocCodeRootRanksGet();
    _cplRootRankCplComm     = _cpl->communicationGet()->cplCommCplCodeRootRanksGet();

    _exchDirection          = exchDirection;

    //_connectableRanks_cpl = _cpl -> communicationGet() -> cplCommCplRanksGet();
    //_connectableRanks     = _cpl -> communicationGet() -> cplCommLocRanksGet();

    //_id     = _localCodeProperties   -> idGet();
    //_id_cpl = _coupledCodeProperties -> idGet();

    _send_buffer.resize(_cpl->fieldsGet()->size());  /*!< Send buffer (size = n_field) */
    _recv_buffer.resize(_cpl->fieldsGet()->size());  /*!< Recv buffer (size = n_field) */

    _send_request.resize(_cpl->fieldsGet()->size()); /*!< Send request (size = n_field) */
    _recv_request.resize(_cpl->fieldsGet()->size()); /*!< Recv request (size = n_field) */

    _send_adler.reserve(_cpl->fieldsGet()->size());
    _recv_adler.reserve(_cpl->fieldsGet()->size());

    for (int i = 0; i < _cpl->fieldsGet()->size(); i++) {
      _send_buffer[i]  = NULL;
      _recv_buffer[i]  = NULL;
      _send_request[i] = -1;
      _recv_request[i] = -1;
    }

    if (!_coupledCodeProperties->localCodeIs()) {
      _cpl->communicationGet()->iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                            1,
                                                                            (void *) &_nPart,
                                                                            -1,
                                                                            NULL,
                                                                            1,
                                                                            (void *) &_cplNPart,
                                                                            -1,
                                                                            NULL);
   
    }

    else {
      cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

      SpatialInterp *cpl_spatial_interp;

      if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet(); 

        cpl_spatial_interp = cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];

      }

      else {

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 

        cpl_spatial_interp = cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];

      }

      int cpl_cplNPart = cpl_cpl.meshGet()-> getNPart();
      _cpl->communicationGet()->iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                            1,
                                                                            (void *) &_nPart,
                                                                            1,
                                                                            (void *) &cpl_cplNPart,
                                                                            1,
                                                                            (void *) &_cplNPart,
                                                                            1,
                                                                            (void *) &(cpl_spatial_interp->_cplNPart));

    }

    _ptsp = NULL;

    if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
      _src_n_gnum = new int [_nPart];
      _src_gnum = new const PDM_g_num_t* [_nPart];

      _tgt_n_gnum = new int [_cplNPart];
      _tgt_gnum = new const PDM_g_num_t* [_cplNPart];

      for (int i_part = 0 ; i_part < _nPart ; i_part++) { 
        _src_n_gnum[i_part] = 0;
        _src_gnum[i_part] = nullptr;
      }

      for (int i_part = 0 ; i_part < _cplNPart ; i_part++) { 
        _tgt_n_gnum[i_part] = 0;
        _tgt_gnum[i_part] = nullptr;
      }
    }
    else {
      _src_n_gnum = new int [_cplNPart];
      _src_gnum = new const PDM_g_num_t* [_cplNPart];

      _tgt_n_gnum = new int [_nPart];
      _tgt_gnum = new const PDM_g_num_t* [_nPart];

      for (int i_part = 0 ; i_part < _cplNPart ; i_part++) { 
        _src_n_gnum[i_part] = 0;
        _src_gnum[i_part] = nullptr;
      }

      for (int i_part = 0 ; i_part < _nPart ; i_part++) { 
        _tgt_n_gnum[i_part] = 0;
        _tgt_gnum[i_part] = nullptr;
      }
    }

    printf("_nPart, _cplNPart : %d %d\n", _nPart, _cplNPart);
    fflush(stdout);

    _n_elt_weights = new int [_nPart];
    _weights_idx = new int* [_nPart];
    _weights = new double* [_nPart];

    _n_computed_tgt = new int [_nPart];
    _computed_tgt = new int* [_nPart];

    _n_uncomputed_tgt = new int [_nPart];
    _uncomputed_tgt = new int* [_nPart];

    for (int i = 0; i < _nPart; i++) {
      _weights_idx[i] = NULL;
      _weights[i] = NULL;
      _computed_tgt[i] = NULL;
      _uncomputed_tgt[i] = NULL;
      _n_uncomputed_tgt[i] = 0;
      _n_computed_tgt[i] = 0;
      _n_elt_weights[i] = 0;
    }

  }


/***************************************************************************/
/***************************************************************************/

  void SpatialInterp::issend(Field* referenceField) {

    if (!_coupledCodeProperties->localCodeIs()) {

      int           *n_elt1;
      int          **selected_part2_idx;
      PDM_g_num_t  **selected_part2;

      const int intId            = referenceField->fieldIDIntGet();
      const CWP_Type_t data_type = referenceField->dataTypeGet();
      const size_t s_data        = sizeof(double);
      const int stride           = referenceField->nComponentGet();
      const int tag              = intId; // Pas bon si les champs n'ont pas été définis dans le même ordre

      PDM_part1_to_selected_part2_selected_part2_get (_ptsp,
                                                      &n_elt1,
                                                      &selected_part2_idx,
                                                      &selected_part2);

      _send_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);
      _recv_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);

      for (int i = 0; i < _nPart; i++) {
        _send_buffer[intId][i] = (double *) malloc(sizeof(double) * stride * n_elt1[i]);
        _recv_buffer[intId][i] = nullptr;
      }

      if (_interpolation_time == CWP_SPATIAL_INTERP_AT_SEND) {
        interpolate (referenceField, _send_buffer[intId]);
      }

      uint32_t mpi_tag = _adler32 (referenceField->fieldIDGet().c_str(), referenceField->fieldIDGet().size()) % MPI_TAG_UB;

      int idx = -1;

      while (idx != -1) {
        if (_send_adler.size()) { 
          idx = PDM_binary_search_uint32t (mpi_tag,
                                       &(_send_adler[0]),
                                       _send_adler.size()) ;
          if (idx != 1) {
            mpi_tag = (mpi_tag + 1) % MPI_TAG_UB;
            if (mpi_tag == 0) {
              mpi_tag += 1;
            }
          }
        }
      }

      if (idx == -1) {
        std::vector<uint32_t>::iterator it  = _send_adler.begin();
        std::vector<uint32_t>::iterator it2 = _send_adler.end();
  
        while(it != _send_adler.end()) {
          if (*it > mpi_tag) {
            it2 = it;
            break;
          }
          it++;
        }
 
        _send_adler.insert (it2, mpi_tag);

      }  

      // Fake reveceive

      PDM_part1_to_selected_part2_irecv (_ptsp,
                                          s_data,
                                          stride,
                                (void **) _recv_buffer[intId],
                                          (int) mpi_tag,
                                         &(_recv_request[intId]));

      PDM_part1_to_selected_part2_issend (_ptsp,
                                          s_data,
                                          stride,
                                (void **) _send_buffer[intId],
                                          (int) mpi_tag,
                                         &(_send_request[intId]));


    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        int           *n_elt1;
        int          **selected_part2_idx;
        PDM_g_num_t  **selected_part2;

        const int intId            = referenceField->fieldIDIntGet();
        const CWP_Type_t data_type = referenceField->dataTypeGet();
        const size_t s_data        = sizeof(double);
        const int stride           = referenceField->nComponentGet();
        const int tag              = intId; // Pas bon si les champs n'ont pas été définis dans le même ordre

        PDM_part1_to_selected_part2_selected_part2_get (_ptsp,
                                                        &n_elt1,
                                                        &selected_part2_idx,
                                                        &selected_part2);

        _send_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);

        for (int i = 0; i < _nPart; i++) {
          _send_buffer[intId][i] = (double *) malloc(sizeof(double) * stride  * n_elt1[i]);
        }

        if (_interpolation_time == CWP_SPATIAL_INTERP_AT_SEND) {
          interpolate (referenceField, _send_buffer[intId]);
        }

        uint32_t mpi_tag = _adler32 (referenceField->fieldIDGet().c_str(), referenceField->fieldIDGet().size()) % MPI_TAG_UB;

        int idx = -1;

        while (idx != -1) {
          if (_send_adler.size()) { 
            idx = PDM_binary_search_uint32t (mpi_tag,
                                         &(_send_adler[0]),
                                         _send_adler.size()) ;
            if (idx != 1) {
              mpi_tag = (mpi_tag + 1) % MPI_TAG_UB;
              if (mpi_tag == 0) {
                mpi_tag += 1;
              }
            }
          }
        }

        if (idx == -1) {
          std::vector<uint32_t>::iterator it  = _send_adler.begin();
          std::vector<uint32_t>::iterator it2 = _send_adler.end();
    
          while(it != _send_adler.end()) {
            if (*it > mpi_tag) {
              it2 = it;
              break;
            }
            it++;
          }
   
          _send_adler.insert (it2, mpi_tag);

        }  

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 
        SpatialInterp *cpl_spatial_interp = cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];

        Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];


        const int cpl_intId = cpl_referenceField->fieldIDIntGet();
        cpl_spatial_interp->_send_buffer[cpl_intId] = (double **) malloc(sizeof(double *) * _cplNPart);
        cpl_spatial_interp->_recv_buffer[cpl_intId] = (double **) malloc(sizeof(double *) * _cplNPart);

        for (int i = 0; i < _cplNPart; i++) {
          cpl_spatial_interp->_recv_buffer[cpl_intId][i] = (double *) malloc(sizeof(double) * stride * cpl_spatial_interp->_n_computed_tgt[i]);
        }

        PDM_part1_to_selected_part2_irecv (_ptsp,
                                           s_data,
                                           stride,
                                 (void **) cpl_spatial_interp->_recv_buffer[cpl_intId],
                                     (int) mpi_tag,
                                          &(cpl_spatial_interp->_recv_request[cpl_intId]));

        PDM_part1_to_selected_part2_issend (_ptsp,
                                            s_data,
                                            stride,
                                  (void **) _send_buffer[intId],
                                      (int) mpi_tag,
                                           &(_send_request[intId]));


      }
    }
  }


  void SpatialInterp::waitIssend(Field* referenceField) {

    if (!_coupledCodeProperties->localCodeIs()) {

      const int intId            = referenceField->fieldIDIntGet();

      PDM_part1_to_selected_part2_irecv_wait (_ptsp, _recv_request[intId]);
      PDM_part1_to_selected_part2_issend_wait (_ptsp, _send_request[intId]);

    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        const int intId            = referenceField->fieldIDIntGet();

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 
        SpatialInterp *cpl_spatial_interp = cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];

        Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];

        const int cpl_intId = cpl_referenceField->fieldIDIntGet();

        PDM_part1_to_selected_part2_irecv_wait (_ptsp, cpl_spatial_interp->_recv_request[cpl_intId]);
        PDM_part1_to_selected_part2_issend_wait (_ptsp, _send_request[intId]);

      }
    }
  }


  void SpatialInterp::irecv(Field *recevingField) {

    // if (!_coupledCodeProperties->localCodeIs()) {

    //   int           *n_elt1;
    //   int          **selected_part2_idx;
    //   PDM_g_num_t  **selected_part2;

    //   const int intId            = referenceField->fieldIDIntGet();
    //   const CWP_Type_t data_type = referenceField->dataTypeGet();
    //   const size_t s_data        = sizeof(double);
    //   const int stride           = referenceField->nComponentGet();
    //   const int tag              = intId; // Pas bon si les champs n'ont pas été définis dans le même ordre

    //   PDM_part1_to_selected_part2_selected_part2_get (_ptsp,
    //                                                   &n_elt1,
    //                                                   &selected_part2_idx,
    //                                                   &selected_part2);

    //   _send_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);
    //   _recv_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);

    //   for (int i = 0; i < _nPart; i++) {
    //     _send_buffer[intId][i] = (double *) malloc(sizeof(double) * stride * n_elt1[i]);
    //     _recv_buffer[intId][i] = nullptr;
    //   }

    //   if (_interpolation_time == CWP_SPATIAL_INTERP_AT_SEND) {
    //     interpolate (referenceField, _send_buffer[intId]);
    //   }

    //   uint32_t mpi_tag = _adler32 (referenceField->fieldIDGet().c_str(), referenceField->fieldIDGet().size()) % MPI_TAG_UB;

    //   int idx = -1;

    //   while (idx != -1) {
    //     if (_send_adler.size()) { 
    //       idx = PDM_binary_search_uint32t (mpi_tag,
    //                                    &(_send_adler[0]),
    //                                    _send_adler.size()) ;
    //       if (idx != 1) {
    //         mpi_tag = (mpi_tag + 1) % MPI_TAG_UB;
    //         if (mpi_tag == 0) {
    //           mpi_tag += 1;
    //         }
    //       }
    //     }
    //   }

    //   if (idx == -1) {
    //     std::vector<uint32_t>::iterator it  = _send_adler.begin();
    //     std::vector<uint32_t>::iterator it2 = _send_adler.end();
  
    //     while(it != _send_adler.end()) {
    //       if (*it > mpi_tag) {
    //         it2 = it;
    //         break;
    //       }
    //       it++;
    //     }
 
    //     _send_adler.insert (it2, mpi_tag);

    //   }  

    //   // Fake reveceive

    //   PDM_part1_to_selected_part2_irecv (_ptsp,
    //                                       s_data,
    //                                       stride,
    //                             (void **) _recv_buffer[intId],
    //                                       (int) mpi_tag,
    //                                      &(_recv_request[intId]));

    //   PDM_part1_to_selected_part2_issend (_ptsp,
    //                                       s_data,
    //                                       stride,
    //                             (void **) _send_buffer[intId],
    //                                       (int) mpi_tag,
    //                                      &(_send_request[intId]));


    // }

    // else {
    //   if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
    //     cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

    //     int           *n_elt1;
    //     int          **selected_part2_idx;
    //     PDM_g_num_t  **selected_part2;

    //     const int intId            = referenceField->fieldIDIntGet();
    //     const CWP_Type_t data_type = referenceField->dataTypeGet();
    //     const size_t s_data        = sizeof(double);
    //     const int stride           = referenceField->nComponentGet();
    //     const int tag              = intId; // Pas bon si les champs n'ont pas été définis dans le même ordre

    //     PDM_part1_to_selected_part2_selected_part2_get (_ptsp,
    //                                                     &n_elt1,
    //                                                     &selected_part2_idx,
    //                                                     &selected_part2);

    //     _send_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);

    //     for (int i = 0; i < _nPart; i++) {
    //       _send_buffer[intId][i] = (double *) malloc(sizeof(double) * stride  * n_elt1[i]);
    //     }

    //     if (_interpolation_time == CWP_SPATIAL_INTERP_AT_SEND) {
    //       interpolate (referenceField, _send_buffer[intId]);
    //     }

    //     uint32_t mpi_tag = _adler32 (referenceField->fieldIDGet().c_str(), referenceField->fieldIDGet().size()) % MPI_TAG_UB;

    //     int idx = -1;

    //     while (idx != -1) {
    //       if (_send_adler.size()) { 
    //         idx = PDM_binary_search_uint32t (mpi_tag,
    //                                      &(_send_adler[0]),
    //                                      _send_adler.size()) ;
    //         if (idx != 1) {
    //           mpi_tag = (mpi_tag + 1) % MPI_TAG_UB;
    //           if (mpi_tag == 0) {
    //             mpi_tag += 1;
    //           }
    //         }
    //       }
    //     }

    //     if (idx == -1) {
    //       std::vector<uint32_t>::iterator it  = _send_adler.begin();
    //       std::vector<uint32_t>::iterator it2 = _send_adler.end();
    
    //       while(it != _send_adler.end()) {
    //         if (*it > mpi_tag) {
    //           it2 = it;
    //           break;
    //         }
    //         it++;
    //       }
   
    //       _send_adler.insert (it2, mpi_tag);

    //     }  

    //     std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 
    //     SpatialInterp *cpl_spatial_interp = cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];

    //     Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];


    //     const int cpl_intId = cpl_referenceField->fieldIDIntGet();
    //     cpl_spatial_interp->_send_buffer[cpl_intId] = (double **) malloc(sizeof(double *) * _cplNPart);
    //     cpl_spatial_interp->_recv_buffer[cpl_intId] = (double **) malloc(sizeof(double *) * _cplNPart);

    //     for (int i = 0; i < _cplNPart; i++) {
    //       cpl_spatial_interp->_recv_buffer[cpl_intId][i] = (double *) malloc(sizeof(double) * stride * cpl_spatial_interp->_n_computed_tgt[i]);
    //     }

    //     PDM_part1_to_selected_part2_irecv (_ptsp,
    //                                        s_data,
    //                                        stride,
    //                              (void **) cpl_spatial_interp->_recv_buffer[cpl_intId],
    //                                  (int) mpi_tag,
    //                                       &(cpl_spatial_interp->_recv_request[cpl_intId]));

    //     PDM_part1_to_selected_part2_issend (_ptsp,
    //                                         s_data,
    //                                         stride,
    //                               (void **) _send_buffer[intId],
    //                                   (int) mpi_tag,
    //                                        &(_send_request[intId]));


    //   }
    // }
  }


  // void SpatialInterp::waitIssend(Field* referenceField) {

  //   if (!_coupledCodeProperties->localCodeIs()) {

  //     const int intId            = referenceField->fieldIDIntGet();

  //     PDM_part1_to_selected_part2_irecv_wait (_ptsp, _recv_request[intId]);
  //     PDM_part1_to_selected_part2_issend_wait (_ptsp, _send_request[intId]);

  //   }

  //   else {
  //     if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
  //       cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

  //       const int intId            = referenceField->fieldIDIntGet();

  //       std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 
  //       SpatialInterp *cpl_spatial_interp = cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];

  //       Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];

  //       const int cpl_intId = cpl_referenceField->fieldIDIntGet();

  //       PDM_part1_to_selected_part2_irecv_wait (_ptsp, cpl_spatial_interp->_recv_request[cpl_intId]);
  //       PDM_part1_to_selected_part2_issend_wait (_ptsp, _send_request[intId]);

  //     }
  //   }
  // }


  void SpatialInterp::waitIrecv(Field* referenceField) {

  }

  
  uint32_t SpatialInterp::_adler32 
  (
    const void *buf,
    size_t buflength
  )
  {

    const uint8_t * buffer = (const uint8_t *)buf;

    uint32_t s1 = 1;
    uint32_t s2 = 0;

    for (size_t n = 0; n < buflength; n++) {
      s1 = (s1 + buffer[n]) % 65521;
      s2 = (s2 + s1) % 65521;
    }

    return (s2 << 16) | s1;
  }



/***************************************************************************/
/***************************************************************************/



} // end namespace cwipi

/**
 * \endcond
 */
