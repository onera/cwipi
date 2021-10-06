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
#include <spatialInterp.hxx>
#include "coupling.hxx"
#include "coupling_i.hxx"
#include <mpi.h>

#include <pdm_mesh_nodal.h>
#include <pdm_dist_cloud_surf.h>
#include <pdm_gnum.h>
#include <pdm_timer.h>
#include <pdm_gnum_location.h>
#include <bftc_error.h>
#include <bftc_printf.h>
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
    _visu                   = coupling->visuGet();

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

    // _connectableRanks_cpl = _cpl -> communicationGet() -> cplCommCplRanksGet();
    // _connectableRanks     = _cpl -> communicationGet() -> cplCommLocRanksGet();

    //_id     = _localCodeProperties   -> idGet();
    //_id_cpl = _coupledCodeProperties -> idGet();


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

    printf("_nPart, _cplNPart : %d %d\n", _nPart, _cplNPart);
    fflush(stdout);

    // if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && slave == 0 ) ){

    //   if(cplComm_rank == _senderRank) {
    //     int tagsend = 2;
    //     int tagrecv = 2;
    //     MPI_Status status;
    //     MPI_Sendrecv(&_nb_part,1,MPI_INT,_senderRank_cpl,tagsend,&_nb_part_cpl,1,MPI_INT,_senderRank_cpl,tagrecv,_cplComm,&status);
    //   }
    // }

    // n_uncomputed_tgt.resize(_nb_part);

    // _coords_target =(double**)     malloc( sizeof(double*)     *_nb_part);
    // _coords_user_targets =(double**)     malloc( sizeof(double*)     *_nb_part);


    // for(int i_part=0;i_part<_nb_part;i_part++){
    //   _coords_user_targets[i_part] = NULL;
    //   _coords_target      [i_part] = NULL;
    // }

    // _gnum_user_targets   =(CWP_g_num_t**)malloc( sizeof(CWP_g_num_t*)*_nb_part);

    // _n_vtx          = (int*)malloc(sizeof(int)*_nb_part);
    // _n_elt          = (int*)malloc(sizeof(int)*_nb_part);
    // _n_target       = (int*)malloc(sizeof(int)*_nb_part);
    // _n_user_targets = (int*)malloc(sizeof(int)*_nb_part);
  }


/***************************************************************************/
/***************************************************************************/

  void SpatialInterp::issend(Field* referenceField) {

    // int  dataTypeSize       = referenceField -> dataTypeSizeGet();
    // int nComponent = referenceField -> nComponentGet();

    // int tag =referenceField -> fieldIDIntGet();
    // void* dist_v_ptr = NULL;

    // void* interpolatedFieldData = interpolate(referenceField);

    // dist_v_ptr = interpolatedFieldData;

    // MPI_Request request;

    // int* displ_send = (int*)malloc(sizeof(int)*cplComm_size);
    // int* count_send = (int*)malloc(sizeof(int)*cplComm_size);
    // for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
    //   count_send[i_proc]= dataTypeSize * nComponent * (_targets_localization_idx_cpl[i_proc][_nb_part]-_targets_localization_idx_cpl[i_proc][0]);
    //   displ_send[i_proc]=  dataTypeSize * nComponent * _targets_localization_idx_cpl[i_proc][0];
    // }

    // int* displ_recv = (int*)malloc(sizeof(int)*cplComm_size);
    // int* count_recv = (int*)malloc(sizeof(int)*cplComm_size);
    // for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
    //   count_recv[i_proc]  =  0 ;
    //   displ_recv [i_proc]  =  0 ;
    // }

    // void* recv_buffer = NULL;

    // MPI_Ialltoallv(dist_v_ptr ,count_send,displ_send,MPI_BYTE,
    //                recv_buffer,count_recv,displ_recv,MPI_BYTE,
    //                _cplComm,&request);


    // free(count_recv);
    // free(displ_recv);
    // free(count_send);
    // free(displ_send);
    // referenceField -> lastRequestAdd(tag,request);
  }

  void SpatialInterp::irecv(Field *recevingField) {
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





/***************************************************************************/
/***************************************************************************/



} // end namespace cwipi

/**
 * \endcond
 */
