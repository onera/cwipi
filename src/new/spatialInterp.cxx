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

      int          **n_elt1;
      int          **selected_part2_idx;
      PDM_g_num_t  **selected_part2;

      PDM_part1_to_selected_part2_selected_part2_get (_ptsp,
                                                      &n_elt1,
                                                      &selected_part2_idx,
                                                      &selected_part2);

    //   // Allocation des buffer !!!

    //   // Remplissage buffer avec l'interpolation

    //   if (_interpolation_time == CWP_SPATIAL_INTERP_AT_SEND) {
    //     interpolate (referenceField);
    //   }

    //   const CWP_Type_t data_type = referenceField->dataTypeGet();
    //   const size_t s_data        = sizeof(double);
    //   const int stride           = referenceField->nComponentGet();
    //   const int fieldIntID       = referenceField->fieldIDIntGet();

    //   PDM_part1_to_selected_part2_issend (_ptsp,
    //                                       s_data,
    //                                       stride,
    //                             (void **) _send_buffer[fieldIntID],
    //                                       fieldIntID,
    //                                       &(_send_request[fieldIntID]));

    //   if (_interpolation_time == CWP_SPATIAL_INTERP_AT_RECV) {
    //     interpolate (referenceField);
    //   }

    // }

    // else {
    //   if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

    //   }
    // }

  }

  void SpatialInterp::waitIssend(Field* referenceField) {

  }


  void SpatialInterp::irecv(Field *recevingField) {



  // int  *ptp2_n_ref_gnum2;
  // int **ptp2_ref_gnum2;
  // PDM_part1_to_selected_part2_ref_gnum2_get (ptp2,
  //                                           &ptp2_n_ref_gnum2,
  //                                           &ptp2_ref_gnum2);

  // int  *ptp2_n_unref_gnum2;
  // int **ptp2_unref_gnum2;
  // PDM_part1_to_selected_part2_unref_gnum2_get (ptp2,
  //                                           &ptp2_n_unref_gnum2,
  //                                           &ptp2_unref_gnum2);


  // int         **ptp2_gnum1_come_from_idx;
  // PDM_g_num_t **ptp2_gnum1_come_from;
  // PDM_part1_to_selected_part2_gnum1_come_from_get (ptp2,
  //                                                 &ptp2_gnum1_come_from_idx,
  //                                                 &ptp2_gnum1_come_from);


      // for (int j = 0; j < n_ref_gnum2[i]; j++) {
      //   for (int k = gnum1_come_from_idx[i][j] ; k < gnum1_come_from_idx[i][j+1]; k++) {
      //     printf(" "PDM_FMT_G_NUM"", gnum1_come_from[i][k]);
      //   }
      //   printf ("\n");
      // }




    // _idx_target  .resize   (_nb_part + 1);
    // _idx_target[0] = 0;
    // for (int i_part = 0; i_part < _nb_part; i_part++) {
    //   _idx_target[i_part+1] = _idx_target[i_part] + _n_target[i_part];
    // }

    // int  dataTypeSize       = recevingField -> dataTypeSizeGet();

    // //Crée un buffer de réception et le stocke (alloue)
    // recevingField -> ReceptionBufferCreation(_n_tot_target);
    // //printf("n_tot_targer %i\n",_n_tot_target);
    // /* Loop on possibly intersecting distant ranks */
    // /*---------------------------------------------*/

    // //Réception des données sur toutes les partitions
    // void* data = recevingField -> recvBufferGet();

    // int nComponent = recevingField -> nComponentGet();

    // void* loc_v_ptr = data;

    // MPI_Request request;
    // int* displ_recv = (int*)malloc(sizeof(int)*cplComm_size);
    // int* count_recv = (int*)malloc(sizeof(int)*cplComm_size);


    // for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
    //   count_recv [i_proc]  =  dataTypeSize * nComponent  * ( _targets_localization_idx[ i_proc ][_nb_part_cpl] - _targets_localization_idx[i_proc][0]  );
    //   displ_recv [i_proc]  =  dataTypeSize * nComponent * _targets_localization_idx[i_proc][0];
    //   /* printf("displ_recv [%i] %i count_recv %i\n",i_proc,displ_recv [i_proc]/(dataTypeSize * nComponent),
    //           _targets_localization_idx[ i_proc ][_nb_part_cpl] - _targets_localization_idx[i_proc][0]); */
    // }

    // int* displ_send = (int*)malloc(sizeof(int)*cplComm_size);
    // int* count_send = (int*)malloc(sizeof(int)*cplComm_size);
    // for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
    //   count_send[i_proc]= 0;
    //   displ_send[i_proc]=  0;
    // }

    // void* send_buffer = NULL;
    // int tag =recevingField -> fieldIDIntGet();

    // MPI_Ialltoallv(send_buffer,count_send,displ_send,MPI_BYTE,
    //                loc_v_ptr,count_recv,displ_recv,MPI_BYTE,
    //                _cplComm,&request);

    // //printf("IRECV %s recevingField -> fieldIDIntGet() %i rank %i both %i request %i\n",
    // //recevingField ->fieldIDGet().c_str(),recevingField -> fieldIDIntGet(),cplComm_rank,_both_codes_are_local,request);

    // free(count_recv);
    // free(displ_recv);
    // free(count_send);
    // free(displ_send);
    // recevingField -> lastRequestAdd (tag,request);
  }

  void SpatialInterp::waitIrecv(Field* referenceField) {

  }




/***************************************************************************/
/***************************************************************************/



} // end namespace cwipi

/**
 * \endcond
 */
