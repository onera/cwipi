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

  void SpatialInterp::_transform_to_index(int** array,int l1, int l2) {

    for (int i_proc = 0; i_proc < l1; i_proc++) {
      int sav = array [i_proc][0];
      if(i_proc==0) array [i_proc][0]=0;
      else array [i_proc][0] = array [i_proc-1][l2];
      array [i_proc][l2]=0;
      for (int i_part = 1; i_part < l2+1; i_part++) {
        int sav2 = array [i_proc][i_part];
        array [i_proc][i_part] = sav + array [i_proc][i_part-1];
        sav = sav2;
      }
    }
  }

  void SpatialInterp::_transform_to_index(int* array,int l1) {

    int sav = array[0];
    array[0]=0;
    array [l1]=0;
    for (int i_part = 1; i_part < l1+1; i_part++) {
      int sav2 = array[i_part];
      array[i_part] = sav + array[i_part-1];
      sav = sav2;
    }
  }

  void SpatialInterp::user_targets_gnum_compute() {
    int coord_def = 1;
    for (int i=0; i<_nb_part; i++){
      if(_coords_user_targets[i] == NULL) {
        coord_def = 0;
        break;
      }
    }

    if(coord_def == 1) {
      int pdmGNum_handle_index  = PDM_gnum_create (3, _nb_part, PDM_FALSE, 1e-3, _pdm_localComm,
                                                   PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);
      for (int i_part=0; i_part<_nb_part; i_part++){
        PDM_gnum_set_from_coords (pdmGNum_handle_index, i_part, _n_user_targets[i_part], _coords_user_targets[i_part], NULL);
      }
      PDM_gnum_compute (pdmGNum_handle_index);
      for (int i_part=0; i_part<_nb_part; i_part++){
        _gnum_user_targets[i_part] = const_cast<CWP_g_num_t*>(PDM_gnum_get(pdmGNum_handle_index, i_part));
      }
    }
  }

/***************************************************************************/
/***************************************************************************/

  void SpatialInterp::issend(Field* referenceField) {

    int  dataTypeSize       = referenceField -> dataTypeSizeGet();
    int nComponent = referenceField -> nComponentGet();

    int tag =referenceField -> fieldIDIntGet();
    void* dist_v_ptr = NULL;

    void* interpolatedFieldData = interpolate(referenceField);

    dist_v_ptr = interpolatedFieldData;

    MPI_Request request;

    int* displ_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_send[i_proc]= dataTypeSize * nComponent * (_targets_localization_idx_cpl[i_proc][_nb_part]-_targets_localization_idx_cpl[i_proc][0]);
      displ_send[i_proc]=  dataTypeSize * nComponent * _targets_localization_idx_cpl[i_proc][0];
    }

    int* displ_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_recv[i_proc]  =  0 ;
      displ_recv [i_proc]  =  0 ;
    }

    void* recv_buffer = NULL;

    MPI_Ialltoallv(dist_v_ptr ,count_send,displ_send,MPI_BYTE,
                   recv_buffer,count_recv,displ_recv,MPI_BYTE,
                   _cplComm,&request);


    free(count_recv);
    free(displ_recv);
    free(count_send);
    free(displ_send);
    referenceField -> lastRequestAdd(tag,request);
  }

  void SpatialInterp::issend_p2p(Field* referenceField) {

    int  dataTypeSize       = referenceField -> dataTypeSizeGet();
    int nComponent = referenceField -> nComponentGet();

    int tag =referenceField -> fieldIDIntGet();
    void* dist_v_ptr = NULL;

    void* interpolatedFieldData = interpolate(referenceField);

    dist_v_ptr = interpolatedFieldData;

    int* displ_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_send[i_proc]= dataTypeSize * nComponent * (_targets_localization_idx_cpl[i_proc][_nb_part]-_targets_localization_idx_cpl[i_proc][0]);
      displ_send[i_proc]=  dataTypeSize * nComponent * _targets_localization_idx_cpl[i_proc][0];
    }

    std::vector<MPI_Request> sreq(_n_ranks_cpl,MPI_REQUEST_NULL);
    int tagcode = 157;

    int ind = 0;
    vector<int>::iterator it;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if(count_send[*it] > 0){
        MPI_Issend(  &( ( (char*)dist_v_ptr )[ displ_send[*it] ] ) ,  count_send[*it], MPI_BYTE,
                     *it, tagcode ,
                     _cplComm, &(sreq[ind]) );
      }
    }


    free(count_send);
    free(displ_send);
    referenceField -> lastRequestAdd_p2p(tag,sreq);
  }

  void SpatialInterp::waitIssend (Field* sendingField) {

    MPI_Status status;
    int tag;
    MPI_Request request;

    tag     = sendingField -> fieldIDIntGet();
    request = sendingField -> lastRequestGet(tag);

    // printf("%s sendingField -> fieldIDIntGet() %i rank %i both %i request %i\n",
    // sendingField ->fieldIDGet().c_str(),sendingField -> fieldIDIntGet(),_rank,_both_codes_are_local,request);
    MPI_Wait(&request, &status);

    // printf("After Wait issend %s sendingField -> fieldIDIntGet() %i rank %i both %i request %i\n",
    // sendingField ->fieldIDGet().c_str(),sendingField -> fieldIDIntGet(),_rank,_both_codes_are_local,request);

    if(_visu -> isCreated() && sendingField -> visuStatusGet() == CWP_STATUS_ON) {
      _visu -> WriterField(sendingField);
    }
  }

  void SpatialInterp::waitIssend_p2p (Field* sendingField) {

    std::vector<MPI_Status> sstatus(_n_ranks_cpl);
    int tag;
    tag     = sendingField -> fieldIDIntGet();
    std::vector<MPI_Request> sreq = sendingField -> lastRequestGet_p2p(tag);

    int ind = 0;
    std::vector<int>::iterator it;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if(sreq[ind]!=MPI_REQUEST_NULL)
        MPI_Wait( &(sreq[ind]), &(sstatus[ind]) );
    }

    if(_visu -> isCreated() && sendingField -> visuStatusGet() == CWP_STATUS_ON) {
      _visu -> WriterField(sendingField);
    }
  }

  void SpatialInterp::irecv(Field *recevingField) {
    _idx_target  .resize   (_nb_part + 1);
    _idx_target[0] = 0;
    for (int i_part = 0; i_part < _nb_part; i_part++) {
      _idx_target[i_part+1] = _idx_target[i_part] + _n_target[i_part];
    }

    int  dataTypeSize       = recevingField -> dataTypeSizeGet();

    //Crée un buffer de réception et le stocke (alloue)
    recevingField -> ReceptionBufferCreation(_n_tot_target);
    //printf("n_tot_targer %i\n",_n_tot_target);
    /* Loop on possibly intersecting distant ranks */
    /*---------------------------------------------*/

    //Réception des données sur toutes les partitions
    void* data = recevingField -> recvBufferGet();

    int nComponent = recevingField -> nComponentGet();

    void* loc_v_ptr = data;

    MPI_Request request;
    int* displ_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);


    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_recv [i_proc]  =  dataTypeSize * nComponent  * ( _targets_localization_idx[ i_proc ][_nb_part_cpl] - _targets_localization_idx[i_proc][0]  );
      displ_recv [i_proc]  =  dataTypeSize * nComponent * _targets_localization_idx[i_proc][0];
      /* printf("displ_recv [%i] %i count_recv %i\n",i_proc,displ_recv [i_proc]/(dataTypeSize * nComponent),
              _targets_localization_idx[ i_proc ][_nb_part_cpl] - _targets_localization_idx[i_proc][0]); */
    }

    int* displ_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_send[i_proc]= 0;
      displ_send[i_proc]=  0;
    }

    void* send_buffer = NULL;
    int tag =recevingField -> fieldIDIntGet();

    MPI_Ialltoallv(send_buffer,count_send,displ_send,MPI_BYTE,
                   loc_v_ptr,count_recv,displ_recv,MPI_BYTE,
                   _cplComm,&request);

    //printf("IRECV %s recevingField -> fieldIDIntGet() %i rank %i both %i request %i\n",
    //recevingField ->fieldIDGet().c_str(),recevingField -> fieldIDIntGet(),_rank,_both_codes_are_local,request);

    free(count_recv);
    free(displ_recv);
    free(count_send);
    free(displ_send);
    recevingField -> lastRequestAdd (tag,request);
  }

  void SpatialInterp::irecv_p2p(Field *recevingField) {
    _idx_target  .resize   (_nb_part + 1);
    _idx_target[0] = 0;
    for (int i_part = 0; i_part < _nb_part; i_part++) {
      _idx_target[i_part+1] = _idx_target[i_part] + _n_target[i_part];
    }

    int  dataTypeSize       = recevingField -> dataTypeSizeGet();

    //Crée un buffer de réception et le stocke (alloue)
    recevingField -> ReceptionBufferCreation(_n_tot_target);
    //printf("n_tot_targer %i\n",_n_tot_target);
    /* Loop on possibly intersecting distant ranks */
    /*---------------------------------------------*/

    //Réception des données sur toutes les partitions
    void* data = recevingField -> recvBufferGet();

    int nComponent = recevingField -> nComponentGet();

    void* loc_v_ptr = data;

    int* displ_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);


    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_recv [i_proc]  =  dataTypeSize * nComponent  * ( _targets_localization_idx[ i_proc ][_nb_part_cpl] - _targets_localization_idx[i_proc][0]  );
      displ_recv [i_proc]  =  dataTypeSize * nComponent * _targets_localization_idx[i_proc][0];
    }


    int tag =recevingField -> fieldIDIntGet();

    std::vector<MPI_Request> rreq(_n_ranks_cpl,MPI_REQUEST_NULL);
    int tagcode = 157;

    int ind = 0;
    vector<int>::iterator it;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if(count_recv[*it] > 0)
        MPI_Irecv(  &((char*) loc_v_ptr)[ displ_recv[*it] ] ,count_recv[*it], MPI_BYTE,
                    *it, tagcode ,
                    _cplComm, &(rreq[ind]) );
    }


    free(count_recv);
    free(displ_recv);
    recevingField -> lastRequestAdd_p2p (tag,rreq);
  }

  void SpatialInterp::waitIrecv (Field* recevingField) {
    MPI_Status status;


    int tag = recevingField -> fieldIDIntGet();
    MPI_Request request = recevingField -> lastRequestGet(tag);
    // printf("%s recevingField -> fieldIDIntGet() %i rank %i both %i request %i\n",
    // recevingField ->fieldIDGet().c_str(),recevingField -> fieldIDIntGet(),_rank,_both_codes_are_local,request);

    if(_both_codes_are_local == 0)
      MPI_Wait(&request, &status);

    //  printf("%s After recevingField -> fieldIDIntGet() %i rank %i both %i \n",
    //  recevingField ->fieldIDGet().c_str(),recevingField -> fieldIDIntGet(),_rank,_both_codes_are_local);
    //Récupère un pointeur vers le bloc de données reçues
    void*              recvData          = recevingField -> recvBufferGet  ();
    int                nComponent        = recevingField -> nComponentGet  ();
    int                dataTypeSize      = recevingField -> dataTypeSizeGet();

    //Reorganize by partition datas which are organized by sending processp
    std::vector<void*> userDataMem (_nb_part,NULL);
    for (int i_part=0;i_part<_nb_part;i_part++) {
      userDataMem [i_part] = recevingField -> dataGet(i_part);
      if(userDataMem[i_part] == NULL ) PDM_error(__FILE__, __LINE__, 0, "Reception memory has not been allocated.\n");
      n_uncomputed_tgt[i_part]=0;
    }

    for(int i_proc=0; i_proc<_n_ranks_cpl;i_proc++) {
      int distant_rank = (*_connectableRanks_cpl)[i_proc];
      /*     printf("itarget [%i] rank %i _nb_part_cpl %i %i %i\n",i_proc,distant_rank,_nb_part_cpl,
                 _targets_localization_idx[ distant_rank ][0],
                 _targets_localization_idx[ distant_rank ][_nb_part_cpl]);
       */    for (int itarget = _targets_localization_idx[ distant_rank ][0]; itarget < _targets_localization_idx[ distant_rank ][_nb_part_cpl]; itarget++) {
        //printf("itarget %i \n",itarget);
        // Index in the interpolated Data array
        int interpInd = itarget;
        int iel = _targets_localization_data[itarget].l_num_origin ;
        int lpart = _targets_localization_data[itarget].origin_part ;
        if(_targets_localization_data[itarget].distance != INFINITY) {
          //Index of the corresponding local reference Data.
          for (int k = 0; k < nComponent; k++) {
            memcpy((char*)userDataMem[lpart] + dataTypeSize * ( nComponent * iel + k ) ,
                   (char*)recvData + dataTypeSize * ( nComponent * interpInd + k ),
                   dataTypeSize);
          }//loop on k
        }
        else {
          n_uncomputed_tgt[lpart]++;
          for (int k = 0; k < nComponent; k++) {
            memcpy((char*)userDataMem[lpart] + dataTypeSize * ( nComponent * iel + k ) ,
                   (char*)recvData + dataTypeSize * ( nComponent * interpInd + k ),
                   dataTypeSize);
            *( (double*) ((char*)userDataMem[lpart] + dataTypeSize * ( nComponent * iel + k ) ) ) = -1.0;
          }//loop on k
        }
      }// loop on itarget
    }// loop on proc

    if(_cpl -> commTypeGet() == CWP_COMM_PAR_WITHOUT_PART){
      int index = 0;
      for (int i_part = 0; i_part < _nb_part; i_part++) {
        memcpy((char*)recvData + dataTypeSize * nComponent * index ,
               (char*)userDataMem[i_part]  ,
               nComponent * dataTypeSize * _n_target[i_part]);
        index += _n_target[i_part];
      }
      MPI_Bcast(recvData,dataTypeSize*nComponent*_n_tot_target,MPI_BYTE, _senderLocalRank, _localCodeProperties -> connectableCommGet());
    }

    if(_visu -> isCreated() && recevingField -> visuStatusGet() == CWP_STATUS_ON
       && (_cpl -> commTypeGet() == CWP_COMM_PAR_WITH_PART || (_cpl -> commTypeGet() == CWP_COMM_PAR_WITHOUT_PART && _rank == _senderRank) ) ) {
      _visu -> WriterField(recevingField);
    }
  }

  void SpatialInterp::waitIrecv_p2p (Field* recevingField) {
    std::vector<MPI_Status> rstatus(_n_ranks_cpl);

    int tag = recevingField -> fieldIDIntGet();
    std::vector<MPI_Request> rreq = recevingField -> lastRequestGet_p2p(tag);


    int ind = 0;
    vector<int>::iterator it;
    if(_both_codes_are_local == 0)
      for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
        if( rreq[ind] != MPI_REQUEST_NULL )
          MPI_Wait( &(rreq[ind]), &(rstatus[ind]) );
      }

    //Récupère un pointeur vers le bloc de données reçues
    void*              recvData          = recevingField -> recvBufferGet  ();
    int                nComponent        = recevingField -> nComponentGet  ();
    int                dataTypeSize      = recevingField -> dataTypeSizeGet();

    //Reorganize by partition datas which are organized by sending processp
    std::vector<void*> userDataMem (_nb_part,NULL);
    for (int i_part=0;i_part<_nb_part;i_part++) {
      userDataMem [i_part] = recevingField -> dataGet(i_part);
      if(userDataMem[i_part] == NULL ) PDM_error(__FILE__, __LINE__, 0, "Reception memory has not been allocated.\n");
      n_uncomputed_tgt[i_part]=0;
    }

    for(int i_proc=0; i_proc<_n_ranks_cpl;i_proc++) {
      int distant_rank = (*_connectableRanks_cpl)[i_proc];
      /*     printf("itarget [%i] rank %i _nb_part_cpl %i %i %i\n",i_proc,distant_rank,_nb_part_cpl,
                 _targets_localization_idx[ distant_rank ][0],
                 _targets_localization_idx[ distant_rank ][_nb_part_cpl]);
       */    for (int itarget = _targets_localization_idx[ distant_rank ][0]; itarget < _targets_localization_idx[ distant_rank ][_nb_part_cpl]; itarget++) {
        //printf("itarget %i \n",itarget);
        // Index in the interpolated Data array
        int interpInd = itarget;
        int iel = _targets_localization_data[itarget].l_num_origin ;
        int lpart = _targets_localization_data[itarget].origin_part ;
        if(_targets_localization_data[itarget].distance != INFINITY) {
          //Index of the corresponding local reference Data.
          for (int k = 0; k < nComponent; k++) {
            memcpy((char*)userDataMem[lpart] + dataTypeSize * ( nComponent * iel + k ) ,
                   (char*)recvData + dataTypeSize * ( nComponent * interpInd + k ),
                   dataTypeSize);
          }//loop on k
        }
        else {
          n_uncomputed_tgt[lpart]++;
          for (int k = 0; k < nComponent; k++) {
            memcpy((char*)userDataMem[lpart] + dataTypeSize * ( nComponent * iel + k ) ,
                   (char*)recvData + dataTypeSize * ( nComponent * interpInd + k ),
                   dataTypeSize);
            *( (double*) ((char*)userDataMem[lpart] + dataTypeSize * ( nComponent * iel + k ) ) ) = -1.0;
          }//loop on k
        }
      }// loop on itarget
    }// loop on proc

    if(_cpl -> commTypeGet() == CWP_COMM_PAR_WITHOUT_PART){
      int index = 0;
      for (int i_part = 0; i_part < _nb_part; i_part++) {
        memcpy((char*)recvData + dataTypeSize * nComponent * index ,
               (char*)userDataMem[i_part]  ,
               nComponent * dataTypeSize * _n_target[i_part]);
        index += _n_target[i_part];
      }
      MPI_Bcast(recvData,dataTypeSize*nComponent*_n_tot_target,MPI_BYTE, _senderLocalRank, _localCodeProperties -> connectableCommGet());
    }

    if(_visu -> isCreated() && recevingField -> visuStatusGet() == CWP_STATUS_ON
       && (_cpl -> commTypeGet() == CWP_COMM_PAR_WITH_PART || (_cpl -> commTypeGet() == CWP_COMM_PAR_WITHOUT_PART && _rank == _senderRank) ) ) {
      _visu -> WriterField(recevingField);
    }
  }

  void SpatialInterp::null_exchange_for_uncoupled_process() {

    int* displ_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_send[i_proc] = 0;
      displ_send[i_proc] = 0;
    }

    int* displ_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_recv[i_proc]  =  0;
      displ_recv [i_proc]  = 0;
    }

    MPI_Request request;

    void* send_buffer=NULL;
    void* recv_buffer=NULL;

    MPI_Ialltoallv(send_buffer, count_send, displ_send,MPI_BYTE,
                   recv_buffer, count_recv, displ_recv,MPI_BYTE,
                   _cplComm,&request);

    free(count_recv);
    free(displ_recv);
    free(count_send);
    free(displ_send);

    MPI_Status status;

    MPI_Wait(&request,&status);
  }

  void SpatialInterp::null_exchange_for_uncoupled_process_p2p() {

    void* send_buffer=NULL;
    void* recv_buffer=NULL;

    std::vector<MPI_Request> rreq(_n_ranks);
    std::vector<MPI_Request> sreq(_n_ranks_cpl);
    int tagcode = 157;

    int ind = 0;
    vector<int>::iterator it;
    /* for( it = _connectableRanks -> begin(); it!=_connectableRanks -> end(); it++,ind++){
       MPI_Irecv(  recv_buffer, 0, MPI_BYTE,
                   *it, tagcode ,
                   _cplComm, &(rreq[ind]) );
     }

     ind = 0;
     for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
       MPI_Issend( send_buffer ,0, MPI_BYTE,
                   *it, tagcode ,
                   _cplComm, &(sreq[ind]) );
     }*/

    std::vector<MPI_Status> sstatus(_n_ranks_cpl);
    std::vector<MPI_Status> rstatus(_n_ranks);

    ind = 0;
    /*    if(_both_codes_are_local == 0)
          for( it = _connectableRanks -> begin(); it!=_connectableRanks -> end(); it++,ind++){
            MPI_Wait( &(rreq[ind]), &(rstatus[ind]) );
          }

          for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
            MPI_Wait( &(sreq[ind]), &(sstatus[ind]) );
          }
     */
  }

  void SpatialInterp::both_codes_on_the_same_process_exchange(Field* referenceField,Field* recevingField) {

    /* Sending section */

    int  dataTypeSize       = referenceField -> dataTypeSizeGet();
    int nComponent = referenceField -> nComponentGet();

    int tag =referenceField -> fieldIDIntGet();
    void* dist_v_ptr = NULL;

    void* interpolatedFieldData = interpolate(referenceField);

    dist_v_ptr = interpolatedFieldData;

    int* displ_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_send[i_proc]= dataTypeSize * nComponent * (_targets_localization_idx_cpl[i_proc][_nb_part]-_targets_localization_idx_cpl[i_proc][0]);
      displ_send[i_proc]=  dataTypeSize * nComponent * _targets_localization_idx_cpl[i_proc][0];
    }

    /* Receving Section */

    _idx_target  .resize   (_spatial_interp_cpl -> _nb_part + 1);
    _idx_target[0] = 0;
    for (int i_part = 0; i_part < _spatial_interp_cpl -> _nb_part; i_part++) {
      _idx_target[i_part+1] = _idx_target[i_part] + _spatial_interp_cpl -> _n_target[i_part];
    }

    int  dataTypeSize_recv = recevingField -> dataTypeSizeGet();
    recevingField -> ReceptionBufferCreation(_spatial_interp_cpl  -> _n_tot_target);

    void* recv_ptr = recevingField -> recvBufferGet();

    int nComponent_recv = recevingField -> nComponentGet();

    int* displ_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_recv[i_proc]  =  dataTypeSize_recv * nComponent_recv
                             * ( _spatial_interp_cpl -> _targets_localization_idx[ i_proc ][_spatial_interp_cpl -> _nb_part_cpl] - _spatial_interp_cpl -> _targets_localization_idx[i_proc][0]  );
      displ_recv [i_proc]  =  dataTypeSize_recv * nComponent_recv * _spatial_interp_cpl -> _targets_localization_idx[i_proc][0];
    }

    MPI_Request request;

    MPI_Ialltoallv(dist_v_ptr ,count_send,displ_send,MPI_BYTE,
                   recv_ptr,count_recv,displ_recv,MPI_BYTE,
                   _cplComm,&request);

    //printf("BOTH %s referenceField -> fieldIDIntGet() %i rank %i request %i\n",
    //referenceField ->fieldIDGet().c_str(),referenceField -> fieldIDIntGet(),_rank,request);


    free(count_recv);
    free(displ_recv);
    free(count_send);
    free(displ_send);
    referenceField -> lastRequestAdd(tag,request);
  }

  void SpatialInterp::both_codes_on_the_same_process_exchange_p2p(Field* referenceField,Field* recevingField) {

    /* Sending section */

    int  dataTypeSize       = referenceField -> dataTypeSizeGet();
    int nComponent = referenceField -> nComponentGet();

    int tag =referenceField -> fieldIDIntGet();
    void* dist_v_ptr = NULL;

    void* interpolatedFieldData = interpolate(referenceField);

    dist_v_ptr = interpolatedFieldData;

    int* displ_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_send[i_proc]= dataTypeSize * nComponent * (_targets_localization_idx_cpl[i_proc][_nb_part]-_targets_localization_idx_cpl[i_proc][0]);
      displ_send[i_proc]=  dataTypeSize * nComponent * _targets_localization_idx_cpl[i_proc][0];
    }

    /* Receving Section */

    _idx_target  .resize   (_spatial_interp_cpl -> _nb_part + 1);
    _idx_target[0] = 0;
    for (int i_part = 0; i_part < _spatial_interp_cpl -> _nb_part; i_part++) {
      _idx_target[i_part+1] = _idx_target[i_part] + _spatial_interp_cpl -> _n_target[i_part];
    }

    int  dataTypeSize_recv = recevingField -> dataTypeSizeGet();
    recevingField -> ReceptionBufferCreation(_spatial_interp_cpl  -> _n_tot_target);

    void* recv_ptr = recevingField -> recvBufferGet();

    int nComponent_recv = recevingField -> nComponentGet();

    int* displ_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i_proc=0; i_proc < _n_ranks_g; i_proc++) {
      count_recv[i_proc]  =  dataTypeSize_recv * nComponent_recv
                             * ( _spatial_interp_cpl -> _targets_localization_idx[ i_proc ][_spatial_interp_cpl -> _nb_part_cpl] - _spatial_interp_cpl -> _targets_localization_idx[i_proc][0]  );
      displ_recv [i_proc]  =  dataTypeSize_recv * nComponent_recv * _spatial_interp_cpl -> _targets_localization_idx[i_proc][0];
    }

    std::vector<MPI_Request> rreq(_n_ranks,MPI_REQUEST_NULL);
    std::vector<MPI_Request> sreq(_n_ranks_cpl,MPI_REQUEST_NULL);
    int tagcode = 157;

    int ind = 0;
    vector<int>::iterator it;
    for( it = _connectableRanks -> begin(); it!=_connectableRanks -> end(); it++,ind++){
      if( count_recv[*it] > 0 )
        MPI_Irecv(  &((char*) recv_ptr)[ displ_recv[*it] ] ,count_recv[*it], MPI_BYTE,
                    *it, tagcode ,
                    _cplComm, &(rreq[ind]) );
    }

    ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if( count_send[*it] >0)
        MPI_Issend(  &((char*) dist_v_ptr)[ displ_send[*it] ] ,count_send[*it], MPI_BYTE,
                     *it, tagcode ,
                     _cplComm, &(sreq[ind]) );
    }

    free(count_recv);
    free(displ_recv);
    free(count_send);
    free(displ_send);
    referenceField -> lastRequestAdd_p2p(tag,sreq);
    recevingField -> lastRequestAdd_p2p(tag,rreq);
  }

  /***********************************************************
   ***********************************************************
   **            Data communication functions               **
   **                                                       **
   ***********************************************************
   ***********************************************************/

  void SpatialInterp::prepare_data_communication_send() {

    for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++)
      free(_process_and_partition_count[i_proc] );
    free(_process_and_partition_count);

    _targets_localization_data_count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    _targets_localization_data_disp_send  = (int*)malloc(sizeof(int)*_n_ranks_g);

    for (int i= 0; i < _n_ranks_g; i++) {
      _targets_localization_data_count_send[i] = sizeof(target_data)*(_targets_localization_idx[i][_nb_part_cpl] - _targets_localization_idx[i][0]);
      _targets_localization_data_disp_send [i] = sizeof(target_data)*_targets_localization_idx[i][0];
    }

  }

  void SpatialInterp::prepare_data_communication_recv() {

    _transform_to_index(_targets_localization_idx_cpl,_n_ranks_g,_nb_part);

    if(_targets_localization_data_cpl!=NULL) free(_targets_localization_data_cpl);
    _targets_localization_data_cpl = (target_data*)malloc(sizeof(target_data)*_targets_localization_idx_cpl[_n_ranks_g-1][_nb_part]);

    _targets_localization_data_count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    _targets_localization_data_disp_recv  = (int*)malloc(sizeof(int)*_n_ranks_g);

    for (int i= 0; i < _n_ranks_g; i++) {
      _targets_localization_data_count_recv[i] = sizeof(target_data)*(_targets_localization_idx_cpl[i][_nb_part] - _targets_localization_idx_cpl[i][0]);
      _targets_localization_data_disp_recv [i] = sizeof(target_data)*_targets_localization_idx_cpl[i][0];
    }

  }

  void SpatialInterp::data_communication_send() {

    int* count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* disp_recv  = (int*)malloc(sizeof(int)*_n_ranks_g);
    for (int i= 0; i < _n_ranks_g; i++) {
      count_recv[i] = 0;
      disp_recv [i] = 0;
    }
    void* recv_buffer_trash = NULL;

    MPI_Request sreq;
    MPI_Ialltoallv((void*)_targets_localization_data, _targets_localization_data_count_send, _targets_localization_data_disp_send, MPI_BYTE,
                   recv_buffer_trash, count_recv,     disp_recv, MPI_BYTE,
                   _cplComm,&sreq);


    MPI_Status stat;
    MPI_Wait(&sreq,&stat);

    free(count_recv);
    free(disp_recv );
  }

  void SpatialInterp::data_communication_send_p2p() {

    std::vector<MPI_Request> sreq(_n_ranks_cpl,MPI_REQUEST_NULL);
    std::vector<MPI_Status> sstatus(_n_ranks_cpl);

    int tagcode = 155;

    int ind = 0;
    vector<int>::iterator it;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      //printf("Send rank %i *it %i \n",_rank,*it);
      if(_targets_localization_data_count_send[*it]>0)
        MPI_Issend( (void*)&( ( (char*)_targets_localization_data ) [_targets_localization_data_disp_send[*it]] )    ,  _targets_localization_data_count_send[*it],  MPI_BYTE,
                    *it, tagcode ,
                    _cplComm, &(sreq[ind]) );
    }

    ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if( sreq[ind] != MPI_REQUEST_NULL )
        MPI_Wait( &(sreq[ind]), &(sstatus[ind]) );
    }

  }

  void SpatialInterp::data_communication_recv() {
    int* count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* disp_send  = (int*)malloc(sizeof(int)*_n_ranks_g);

    for (int i= 0; i < _n_ranks_g; i++) {
      count_send[i] = 0;
      disp_send [i] = 0;
    }

    void* sbuffer_trash = NULL;

    MPI_Request rreq;

    MPI_Ialltoallv(sbuffer_trash                       ,  count_send              , disp_send              , MPI_BYTE,
                   (void*)_targets_localization_data_cpl,  _targets_localization_data_count_recv, _targets_localization_data_disp_recv, MPI_BYTE,
                   _cplComm,&rreq);

    MPI_Status stat;
    MPI_Wait(&rreq,&stat);
    free(count_send);
    free(disp_send );

  }

  void SpatialInterp::data_communication_recv_p2p() {

    std::vector<MPI_Request> rreq(_n_ranks_cpl, MPI_REQUEST_NULL);
    std::vector<MPI_Status> rstatus(_n_ranks_cpl);

    int tagcode = 155;

    int ind = 0;
    vector<int>::iterator it;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if( _targets_localization_data_count_recv[*it] > 0 )
        MPI_Irecv( (void*)&( ((char*)_targets_localization_data_cpl)[ _targets_localization_data_disp_recv[*it] ] ),  _targets_localization_data_count_recv[*it],  MPI_BYTE,
                   *it, tagcode,
                   _cplComm, &(rreq[ind]) );
    }

    ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if( rreq[ind] != MPI_REQUEST_NULL )
        MPI_Wait( &(rreq[ind]), &(rstatus[ind]) );
    }

  }

  void SpatialInterp::both_data_communication() {

    MPI_Request req;
    MPI_Ialltoallv((void*) (_spatial_interp_cpl -> _targets_localization_data), _spatial_interp_cpl -> _targets_localization_data_count_send, _spatial_interp_cpl -> _targets_localization_data_disp_send, MPI_BYTE,
                   (void*)_targets_localization_data_cpl,  _targets_localization_data_count_recv, _targets_localization_data_disp_recv, MPI_BYTE,
                   _cplComm,&req);

    MPI_Status stat;
    MPI_Wait(&req,&stat);

  }

  void SpatialInterp::both_data_communication_p2p() {

    std::vector<MPI_Request> sreq(_n_ranks,MPI_REQUEST_NULL);
    std::vector<MPI_Status> sstatus(_n_ranks);
    std::vector<MPI_Request> rreq(_n_ranks_cpl, MPI_REQUEST_NULL);
    std::vector<MPI_Status> rstatus(_n_ranks_cpl);

    int tagcode = 155;


    vector<int>::iterator it;
    int ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if(  _targets_localization_data_count_recv[*it] > 0 )
        MPI_Irecv((void*)&( ((char*)_targets_localization_data_cpl)[ _targets_localization_data_disp_recv[*it] ] ),  _targets_localization_data_count_recv[*it],  MPI_BYTE,
                  *it, tagcode,
                  _cplComm, &(rreq[ind]) );
    }

    ind = 0;
    for( it = _connectableRanks -> begin(); it!=_connectableRanks -> end(); it++,ind++){
      if( _spatial_interp_cpl -> _targets_localization_data_count_send[*it] > 0 )
        MPI_Issend( (void*)&( ( (char*) _spatial_interp_cpl -> _targets_localization_data ) [_spatial_interp_cpl -> _targets_localization_data_disp_send[*it]] )    ,
                    _spatial_interp_cpl -> _targets_localization_data_count_send[*it],  MPI_BYTE,
                    *it, tagcode ,
                    _cplComm, &(sreq[ind]) );
    }


    ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if( rreq[ind] != MPI_REQUEST_NULL )
        MPI_Wait( &(rreq[ind]), &(rstatus[ind]) );
    }

    ind = 0;
    for( it = _connectableRanks -> begin(); it!=_connectableRanks -> end(); it++,ind++){
      if( sreq[ind] != MPI_REQUEST_NULL )
        MPI_Wait( &(sreq[ind]), &(sstatus[ind]) );
    }


  }

  void SpatialInterp::data_communication_null() {

    int* count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* disp_send  = (int*)malloc(sizeof(int)*_n_ranks_g);

    for (int i= 0; i < _n_ranks_g; i++) {
      count_send[i] = 0;
      disp_send [i] = 0;
    }
    void* sbuffer_trash = NULL;

    int* count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);
    int* disp_recv  = (int*)malloc(sizeof(int)*_n_ranks_g);

    for (int i= 0; i < _n_ranks_g; i++) {
      count_recv[i] = 0;
      disp_recv [i] = 0;
    }

    void* recv_buffer_trash = NULL;

    MPI_Request rreq;
    MPI_Ialltoallv(sbuffer_trash                       ,  count_send              , disp_send              , MPI_BYTE,
                   recv_buffer_trash, count_recv,     disp_recv, MPI_BYTE,
                   _cplComm,&rreq);

    MPI_Status stat;
    MPI_Wait(&rreq,&stat);

    free(count_send);
    free(disp_send );
    free(count_recv);
    free(disp_recv );

  }

  void SpatialInterp::data_communication_wait_send() {
    free(_targets_localization_data_count_send);
    free(_targets_localization_data_disp_send );
  }

  void SpatialInterp::data_communication_wait_recv() {
    free(_targets_localization_data_count_recv);
    free(_targets_localization_data_disp_recv );
    _n_tot_target_cpl    = _targets_localization_idx_cpl[_n_ranks_g-1][_nb_part];
  }

  void SpatialInterp::computeFree(){

    for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
      if(_targets_localization_idx_cpl != NULL) free(_targets_localization_idx_cpl[i_proc]);
      if(_targets_localization_idx     != NULL) free(_targets_localization_idx[i_proc]);
    }

    free(_targets_localization_data_cpl);
    if(_targets_localization_idx_cpl != NULL) free(_targets_localization_idx_cpl);
    if(_Texch_t == CWP_FIELD_EXCH_RECV) free(_targets_localization_data);
    if(_targets_localization_idx     != NULL) free(_targets_localization_idx);
  }

  /***********************************************************
   ***********************************************************
   **            Data index communication functions         **
   **                                                       **
   ***********************************************************
   ***********************************************************/

  void SpatialInterp::data_index_communication_send() {

    int* sbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part_cpl);
    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part<_nb_part_cpl;i_part++) {
        sbuffer[ i_proc * _nb_part_cpl + i_part ] = _process_and_partition_count[i_proc][i_part];
      }

    int* recvbuffer_trash = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part_cpl);

    MPI_Request sreq;

    MPI_Ialltoall(sbuffer, _nb_part_cpl, MPI_INT,
                  recvbuffer_trash, _nb_part_cpl, MPI_INT,
                  _cplComm,&sreq);


    MPI_Status stat;
    MPI_Wait(&sreq,&stat);

    free(sbuffer   );
    free(recvbuffer_trash   );
  }

  void SpatialInterp::data_index_communication_send_p2p() {

    int* sbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part_cpl);
    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part<_nb_part_cpl;i_part++) {
        sbuffer[ i_proc * _nb_part_cpl + i_part ] = _process_and_partition_count[i_proc][i_part];
      }


    std::vector<MPI_Request> sreq(_n_ranks_cpl);
    std::vector<MPI_Status> sstatus(_n_ranks_cpl);

    int tagcode = 152;

    int ind = 0;
    vector<int>::iterator it;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      MPI_Issend(&(sbuffer[*it * _nb_part_cpl]),_nb_part_cpl,MPI_INT,
                 *it, tagcode ,
                 _cplComm, &(sreq[ind]) );
      //printf("Send rank %i ind %i req %i tag %i it %i\n",_rank,ind,&(sreq[ind]),tagcode + *it,*it);
    }

    ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      MPI_Wait( &(sreq[ind]), &(sstatus[ind]) );
      //printf("Wait rank %i ind %i req %i tag %i\n",_rank,ind,&(sreq[ind]),tagcode + *it);
    }

    //printf("After Wait Send rank %i \n",_rank);

    free(sbuffer   );
  }

  void SpatialInterp::data_index_communication_recv() {

    std::vector<int> connectable = *_connectableRanks_cpl;

    int* recvbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part);
    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        recvbuffer[ i_proc * _nb_part + i_part ] = 0;
      }

    int *sendbuffer_trash = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part);
    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        sendbuffer_trash[ i_proc * _nb_part + i_part ] = 0;
      }

    MPI_Request rreq;

    MPI_Ialltoall(sendbuffer_trash, _nb_part, MPI_INT,
                  recvbuffer, _nb_part, MPI_INT,
                  _cplComm,&rreq);

    MPI_Status stat;
    MPI_Wait(&rreq,&stat);


    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        _targets_localization_idx_cpl[i_proc][i_part] = recvbuffer[ i_proc * _nb_part + i_part ];
      }

    if(sendbuffer_trash!=NULL) free(sendbuffer_trash   );
    if(recvbuffer!=NULL) free(recvbuffer   );
  }

  void SpatialInterp::data_index_communication_recv_p2p() {

    int* recvbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part);
    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        recvbuffer[ i_proc * _nb_part + i_part ] = 0;
      }


    std::vector<MPI_Request> rreq(_n_ranks_cpl);
    std::vector<MPI_Status> rstatus(_n_ranks_cpl);

    int tagcode = 152;

    int ind = 0;
    vector<int>::iterator it;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      MPI_Irecv(&(recvbuffer[*it * _nb_part]),_nb_part,MPI_INT,
                *it, tagcode,
                _cplComm, &(rreq[ind]) );
    }

    ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      MPI_Wait( &(rreq[ind]), &(rstatus[ind]) );
    }

    //printf("After Wait Recv rank %i \n",_rank);

    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        _targets_localization_idx_cpl[i_proc][i_part] = recvbuffer[ i_proc * _nb_part + i_part ];
      }

    free(recvbuffer   );
  }

  void SpatialInterp::both_index_communication_p2p() {

    int* sbuffer = (int*) malloc(sizeof(int)*_n_ranks_g * _nb_part );
    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part< _nb_part; i_part++) {
        sbuffer[ i_proc * _nb_part + i_part ] = (_spatial_interp_cpl -> _process_and_partition_count)[i_proc][i_part];
      }

    int* recvbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part);
    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        recvbuffer[ i_proc * _nb_part + i_part ] = 0;
      }

    std::vector<MPI_Request> sreq(_n_ranks);
    std::vector<MPI_Status> sstatus(_n_ranks);

    std::vector<MPI_Request> rreq(_n_ranks_cpl);
    std::vector<MPI_Status> rstatus(_n_ranks_cpl);

    int tagcode = 152;

    vector<int>::iterator it;
    int ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if(*it != _rank)
        MPI_Irecv(&(recvbuffer[*it * _nb_part]),_nb_part,MPI_INT,
                  *it, tagcode,
                  _cplComm, &(rreq[ind]) );
      else
        memcpy(&(recvbuffer[*it * _nb_part]), &(sbuffer[*it * _nb_part]), sizeof(int)*_nb_part );

      //       printf("After sendrecv rank %i ind %i it %i %i _n_ranks_cpl %i\n",_rank,ind,*it,rreq[ind],_n_ranks_cpl);

    }


    ind = 0;
    for( it = _connectableRanks -> begin(); it!=_connectableRanks -> end(); it++,ind++){
      if(*it != _rank)
        MPI_Issend(&(sbuffer[*it * _nb_part]),_nb_part,MPI_INT,
                   *it, tagcode ,
                   _cplComm, &(sreq[ind]) );
      else
        memcpy(&(recvbuffer[*it * _nb_part]), &(sbuffer[*it * _nb_part]), sizeof(int)*_nb_part );
    }




    ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      //printf("After rank %i MPI _wait ind %i it %i totot %i %i %i _n_ranks_cpl %i\n",_rank,ind,*it,toto,titi,rreq[ind],_n_ranks_cpl);
      if(*it != _rank)
        MPI_Wait( &(rreq[ind]), &(rstatus[ind]) );
    }

    ind = 0;
    for( it = _connectableRanks -> begin(); it!=_connectableRanks -> end(); it++,ind++){
      if(*it != _rank)
        MPI_Wait( &(sreq[ind]), &(sstatus[ind]) );
    }


    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        _targets_localization_idx_cpl[i_proc][i_part] = recvbuffer[ i_proc * _nb_part + i_part ];
      }

    free(sbuffer   );
    free(recvbuffer   );
  }

  void SpatialInterp::data_index_communication_null() {

    std::vector<int> connectable = *_connectableRanks_cpl;

    MPI_Request request;
    int* sbuffer = NULL;
    int* recvbuffer = NULL;
    if(_Texch_t == CWP_FIELD_EXCH_SEND){

      sbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part);
      for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
        for(int i_part=0;i_part<_nb_part;i_part++)
          sbuffer[ i_proc * _nb_part + i_part ] = 0;

      recvbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part);

      MPI_Ialltoall(sbuffer, _nb_part, MPI_INT,
                    recvbuffer, _nb_part, MPI_INT,
                    _cplComm,&request);

    }
    else if(_Texch_t == CWP_FIELD_EXCH_RECV) {

      // if(_targets_localization_idx==NULL) {
      _targets_localization_idx   =(int**)malloc(sizeof(int*)*_n_ranks_g);
      for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++)
        _targets_localization_idx [i_proc] = NULL;
      //  }
      _process_and_partition_count =(int**)malloc(sizeof(int*)*_n_ranks_g);

      for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
        // if(_targets_localization_idx [i_proc] == NULL)
        _targets_localization_idx [i_proc] = (int*) malloc( sizeof(int) * (1+_nb_part_cpl) );
        _process_and_partition_count [i_proc]=(int*)malloc(sizeof(int)*(1+_nb_part_cpl));

        for (int i_part = 0; i_part < _nb_part_cpl+1; i_part++) {
          _targets_localization_idx    [i_proc][i_part] = 0;
          _process_and_partition_count [i_proc][i_part] = 0;
        }
      }



      sbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part_cpl);
      for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
        for(int i_part=0;i_part<_nb_part_cpl;i_part++)
          sbuffer[ i_proc * _nb_part_cpl + i_part ] = 0;

      recvbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part_cpl);

      MPI_Ialltoall(sbuffer, _nb_part_cpl, MPI_INT,
                    recvbuffer, _nb_part_cpl, MPI_INT,
                    _cplComm,&request);

    }
    MPI_Status stat;
    MPI_Wait(&request,&stat);

    if(sbuffer!=NULL) free(sbuffer   );
    if(recvbuffer!=NULL) free(recvbuffer   );

  }

  void SpatialInterp::both_index_communication() {

    int* sbuffer = (int*) malloc(sizeof(int)*_n_ranks_g * _nb_part );
    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part< _nb_part; i_part++) {
        sbuffer[ i_proc * _nb_part + i_part ] = (_spatial_interp_cpl -> _process_and_partition_count)[i_proc][i_part];
      }

    int* recvbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part);

    MPI_Request sreq;

    MPI_Ialltoall(sbuffer, _nb_part, MPI_INT,
                  recvbuffer, _nb_part, MPI_INT,
                  _cplComm,&sreq);


    MPI_Status stat;
    MPI_Wait(&sreq,&stat);

    for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        _targets_localization_idx_cpl[i_proc][i_part] = recvbuffer[ i_proc * _nb_part + i_part ];
      }

    free(sbuffer   );
    free(recvbuffer   );
  }

  void SpatialInterp::_IAlltoallIndexSend(void* send_buffer,
                                     int* send_count,
                                     int* send_disp,
                                     MPI_Datatype type,
                                     MPI_Comm comm,
                                     std::vector<int> connectableRanks
                                    ){
      int comm_size;
      MPI_Comm_size(comm,&comm_size);
      int nranks = connectableRanks.size();

      _send_requests.resize(nranks,0);

      int tagsend = -1;
      if(localName == _codeVector[0]) {
        tagsend =0;
      }
      else {
       tagsend =1;
      }

      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = connectableRanks[i_rank];

        if(type == MPI_BYTE){
           int* sendptr = (int*)send_buffer;
          //  target_data* sendptr = (target_data*)send_buffer;
            MPI_Issend(&(  sendptr [ send_disp[distant_rank]/sizeof(int) /*recv_size[distant_rank]*/ ] ), send_count[distant_rank] * 1, type, distant_rank, tagsend,
                   comm,
                   &(_send_requests[i_rank]));
        }
      }//end for on i_rank
  }




 void SpatialInterp::_IAlltoallIndexRecv(void* recv_buffer,
                int* recv_count,
                int* recv_disp,
                MPI_Datatype type,
                MPI_Comm comm,
                std::vector<int> connectableRanks
                ){

      int comm_size;
      MPI_Comm_size(comm,&comm_size);
      int nranks = connectableRanks.size();

      _recv_requests.resize(nranks,0);

      int tagrecv = -1;
      if(localName == _codeVector[0]) {
        tagrecv = 1;
      }
      else {
        tagrecv = 0;
      }

      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = connectableRanks[i_rank];
        if(type == MPI_BYTE){
          int* recvptr = (int*)recv_buffer;
          //target_data* recvptr = (target_data*)recv_buffer;
          MPI_Irecv(&(  recvptr [ recv_disp[distant_rank]/sizeof(int)  /*recv_size[distant_rank]*/ ] ), recv_count[distant_rank] * 1,type, distant_rank, tagrecv,
                    comm,
                    &(_recv_requests[i_rank]));
        }
      }//end for on i_rank

  }

/***************************************************************************/
/***************************************************************************/



} // end namespace cwipi

/**
 * \endcond
 */
