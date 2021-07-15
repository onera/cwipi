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

  void SpatialInterp::init(Coupling *coupling, CWP_Dof_location_t pointsCloudLocation,bool slave) {
    _mesh   = coupling -> meshGet();
    _visu   = coupling -> visuGet();
    _pointsCloudLocation = pointsCloudLocation;
    _cpl = coupling;
    _localCodeProperties = _cpl -> localCodePropertiesGet();
    _coupledCodeProperties = _cpl -> coupledCodePropertiesGet();

    _slave = slave;
    _nb_part = _mesh -> getNPart();

    _cplComm = _cpl -> communicationGet() -> cplCommGet();
    _globalComm = _localCodeProperties -> globalCommGet();
    _localComm = _mesh -> getMPIComm();
    _pdm_cplComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_cplComm));

    localName   = _localCodeProperties -> nameGet();
    coupledName = _coupledCodeProperties -> nameGet();

    _senderRank     = _cpl -> communicationGet() -> unionCommLocCodeRootRanksGet();
    _senderRank_cpl = _cpl -> communicationGet() -> unionCommCplCodeRootRanksGet();

    _connectableRanks_cpl = _cpl -> communicationGet() -> cplCommCplRanksGet();
    _connectableRanks     = _cpl -> communicationGet() -> cplCommLocRanksGet();
    localComm_size = _connectableRanks->size();
    localComm_size_cpl = _connectableRanks_cpl->size();

    MPI_Group cplGroup, localGroup;
    MPI_Comm_group(_localComm, &localGroup);
    MPI_Comm_group(_cplComm, &cplGroup);
    MPI_Group_translate_ranks(cplGroup, 1, &_senderRank, localGroup, &_senderLocalRank);

    int globalComm_size, globalComm_rank;
    MPI_Comm_size(_globalComm, &globalComm_size);
    MPI_Comm_rank(_globalComm, &globalComm_rank);
    MPI_Comm_size(_cplComm, &cplComm_size);
    MPI_Comm_rank(_cplComm, &cplComm_rank);
    MPI_Comm_size(_localComm, &localComm_size);

    int comp = localName.compare(coupledName);
    _codeVector.resize(2);
    if(comp>0) {
      _codeVector[0] = localName  ;
      _codeVector[1] = coupledName;
    }
    else {
      _codeVector[1] = localName  ;
      _codeVector[0] = coupledName;
    }


    _both_codes_are_local=0;
    if(_localCodeProperties ->localCodeIs() && _coupledCodeProperties ->localCodeIs()) {
      _both_codes_are_local = 1;
      _nb_part_cpl = _nb_part;//fix?
    }


    _id     = _localCodeProperties   -> idGet();
    _id_cpl = _coupledCodeProperties -> idGet();

    if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && slave == 0 ) ){

      if(cplComm_rank == _senderRank) {
        int tagsend = 2;
        int tagrecv = 2;
        MPI_Status status;
        MPI_Sendrecv(&_nb_part,1,MPI_INT,_senderRank_cpl,tagsend,&_nb_part_cpl,1,MPI_INT,_senderRank_cpl,tagrecv,_cplComm,&status);
      }
    }

    _isActiveRank     = _localCodeProperties   -> isActiveRank();

    n_uncomputed_tgt.resize(_nb_part);

    _gnum_target   =(CWP_g_num_t**)malloc( sizeof(CWP_g_num_t*)*_nb_part);
    _coords_target =(double**)     malloc( sizeof(double*)     *_nb_part);
    _coords_user_targets =(double**)     malloc( sizeof(double*)     *_nb_part);


    for(int i_part=0;i_part<_nb_part;i_part++){
      _coords_user_targets[i_part] = NULL;
      _coords_target      [i_part] = NULL;
    }

    _gnum_user_targets   =(CWP_g_num_t**)malloc( sizeof(CWP_g_num_t*)*_nb_part);

    _n_vtx          = (int*)malloc(sizeof(int)*_nb_part);
    _n_elt          = (int*)malloc(sizeof(int)*_nb_part);
    _n_target       = (int*)malloc(sizeof(int)*_nb_part);
    _n_user_targets = (int*)malloc(sizeof(int)*_nb_part);
  }

  /***********************************************************
**           Mesh information functions                  **
***********************************************************/

  void SpatialInterp::mesh_info_get() {

    _n_tot_elt         =0;
    _n_tot_vtx         =0;
    _n_tot_user_targets=0;
    for(int i_part =0;i_part<_nb_part;i_part++) {

      if (_pointsCloudLocation == CWP_DOF_LOCATION_CELL_CENTER && _Texch_t == CWP_FIELD_EXCH_RECV ) {
        _n_target   [i_part]     = _mesh -> getPartNElts(i_part);
        _gnum_target[i_part]     = _mesh -> GNumEltsGet(i_part);
        _coords_target [i_part]  = _mesh -> eltCentersGet(i_part);
      }
      else if (_pointsCloudLocation == CWP_DOF_LOCATION_NODE && _Texch_t == CWP_FIELD_EXCH_RECV) {
        _n_target      [i_part]  = _mesh -> getPartNVertex (i_part);
        _gnum_target   [i_part]  = _mesh -> getVertexGNum  (i_part);
        _coords_target [i_part]  = _mesh -> getVertexCoords(i_part);
      }
      else if (_pointsCloudLocation == CWP_DOF_LOCATION_USER && _Texch_t == CWP_FIELD_EXCH_RECV ) {
        //   printf("info_mesh _n_target [i_part] %i _n_user_targets[i_part] %i\n",_n_target[i_part],_n_user_targets[i_part]);
        _n_target      [i_part]  = _n_user_targets     [i_part];
        //   printf("info_mesh _n_target [i_part] %i _n_user_targets[i_part] %i\n",_n_target[i_part],_n_user_targets[i_part]);
        _gnum_target   [i_part]  = _gnum_user_targets  [i_part];
        _coords_target [i_part]  = _coords_user_targets[i_part];
      }

      _n_elt[i_part]  = _mesh -> getPartNElts(i_part);
      _n_tot_elt+=_n_elt[i_part];

      _n_vtx[i_part]  = _mesh -> getPartNVertex(i_part);
      _n_tot_vtx+=_n_vtx[i_part];

      _n_tot_user_targets += _n_user_targets [i_part];

    } //end loop on i_part


    if (_pointsCloudLocation == CWP_DOF_LOCATION_CELL_CENTER && _Texch_t == CWP_FIELD_EXCH_RECV) {
      _n_tot_target = _n_tot_elt;
    }
    else if (_pointsCloudLocation == CWP_DOF_LOCATION_NODE && _Texch_t == CWP_FIELD_EXCH_RECV) {
      _n_tot_target = _n_tot_vtx;
    }
    else if (_pointsCloudLocation == CWP_DOF_LOCATION_USER && _Texch_t == CWP_FIELD_EXCH_RECV) {
      _n_tot_target = _n_tot_user_targets;
    }

    MPI_Barrier(_localComm);

    if(_cpl -> commTypeGet() == CWP_COMM_PAR_WITH_PART){
      /************* Elements ***********/
      CWP_g_num_t n_tot_elt_long = (CWP_g_num_t)_n_tot_elt;
      MPI_Reduce(&n_tot_elt_long,&_n_g_elt_over_part,1,MPI_LONG,MPI_SUM,0,_localComm);
      MPI_Bcast(&_n_g_elt_over_part,1,MPI_LONG,0,_localComm);

      /************* Vertices ***********/
      CWP_g_num_t n_tot_vtx_long = (CWP_g_num_t)_n_tot_vtx;
      MPI_Reduce(&n_tot_vtx_long,&_n_g_vtx_over_part,1,MPI_LONG,MPI_SUM,0,_localComm);
      MPI_Bcast(&_n_g_vtx_over_part,1,MPI_LONG,0,_localComm);
    }
    else {
      _n_g_vtx_over_part = (CWP_g_num_t)_n_tot_vtx;
      _n_g_elt_over_part = (CWP_g_num_t)_n_tot_elt;
    }

  }

  void SpatialInterp::mesh_cpl_info_get() {

    if(_slave==0){
      /*      Partition Number exchange           */

      if(cplComm_rank == _senderRank) {
        int tagsend = 2;
        int tagrecv = 2;
        MPI_Status status;
        MPI_Sendrecv(&_n_g_elt_over_part,1,MPI_LONG,_senderRank_cpl,tagsend,&_n_g_elt_cpl_over_part,1,MPI_LONG,_senderRank_cpl,tagrecv,_cplComm,&status);
        tagsend = 3;
        tagrecv = 3;
        MPI_Sendrecv(&_n_g_vtx_over_part,1,MPI_LONG,_senderRank_cpl,tagsend,&_n_g_vtx_cpl_over_part,1,MPI_LONG,_senderRank_cpl,tagrecv,_cplComm,&status);
      }

      CWP_g_num_t tmp1;

      if(_id < _id_cpl) {
        MPI_Bcast(&_n_g_elt_cpl_over_part,1,MPI_LONG,_senderRank,_cplComm);
        MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_cplComm);
      }
      else{
        MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_cplComm);
        MPI_Bcast(&_n_g_elt_cpl_over_part,1,MPI_LONG,_senderRank,_cplComm);
      }

      if(_id < _id_cpl) {
        MPI_Bcast(&_n_g_vtx_cpl_over_part,1,MPI_LONG,_senderRank,_cplComm);
        MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_cplComm);
      }
      else{
        MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_cplComm);
        MPI_Bcast(&_n_g_vtx_cpl_over_part,1,MPI_LONG,_senderRank,_cplComm);
      }
    }
    else{
      _n_g_elt_cpl_over_part = _spatial_interp_cpl->_n_g_elt_over_part;
      _n_g_vtx_cpl_over_part = _spatial_interp_cpl->_n_g_vtx_over_part;
    }
  }

  void SpatialInterp::info_mesh() {

    if(_both_codes_are_local == 0){
      if(_isActiveRank)  mesh_info_get();
      mesh_cpl_info_get();
    }
    else if(_Texch_t == CWP_FIELD_EXCH_SEND) {
      mesh_info_get();
      _spatial_interp_cpl -> mesh_info_get();

      mesh_cpl_info_get();
      _spatial_interp_cpl -> mesh_cpl_info_get();
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
      PDM_gen_gnum_t *pdmGNum_handle_index  = PDM_gnum_create (3, _nb_part, PDM_FALSE, 1e-3, _mesh->_pdm_localComm,
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

  void SpatialInterp::user_target_points_set(int i_part, int n_pts, double *coord) {
    if (_pointsCloudLocation != CWP_DOF_LOCATION_USER) PDM_error(__FILE__, __LINE__, 0, "You cannot use user_target_points_set for CWP_Dof_location_t different of CWP_DOF_LOCATION_USER.\n");
      else {
        _n_user_targets[i_part] = n_pts;
          _coords_user_targets[i_part] = coord;
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

    int* displ_send = (int*)malloc(sizeof(int)*cplComm_size);
    int* count_send = (int*)malloc(sizeof(int)*cplComm_size);
    for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
      count_send[i_proc]= dataTypeSize * nComponent * (_targets_localization_idx_cpl[i_proc][_nb_part]-_targets_localization_idx_cpl[i_proc][0]);
      displ_send[i_proc]=  dataTypeSize * nComponent * _targets_localization_idx_cpl[i_proc][0];
    }

    int* displ_recv = (int*)malloc(sizeof(int)*cplComm_size);
    int* count_recv = (int*)malloc(sizeof(int)*cplComm_size);
    for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
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

    int* displ_send = (int*)malloc(sizeof(int)*cplComm_size);
    int* count_send = (int*)malloc(sizeof(int)*cplComm_size);
    for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
      count_send[i_proc]= dataTypeSize * nComponent * (_targets_localization_idx_cpl[i_proc][_nb_part]-_targets_localization_idx_cpl[i_proc][0]);
      displ_send[i_proc]=  dataTypeSize * nComponent * _targets_localization_idx_cpl[i_proc][0];
    }

    std::vector<MPI_Request> sreq(localComm_size_cpl,MPI_REQUEST_NULL);
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

  void SpatialInterp::waitIssend_p2p (Field* sendingField) {

    std::vector<MPI_Status> sstatus(localComm_size_cpl);
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
    int* displ_recv = (int*)malloc(sizeof(int)*cplComm_size);
    int* count_recv = (int*)malloc(sizeof(int)*cplComm_size);


    for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
      count_recv [i_proc]  =  dataTypeSize * nComponent  * ( _targets_localization_idx[ i_proc ][_nb_part_cpl] - _targets_localization_idx[i_proc][0]  );
      displ_recv [i_proc]  =  dataTypeSize * nComponent * _targets_localization_idx[i_proc][0];
      /* printf("displ_recv [%i] %i count_recv %i\n",i_proc,displ_recv [i_proc]/(dataTypeSize * nComponent),
              _targets_localization_idx[ i_proc ][_nb_part_cpl] - _targets_localization_idx[i_proc][0]); */
    }

    int* displ_send = (int*)malloc(sizeof(int)*cplComm_size);
    int* count_send = (int*)malloc(sizeof(int)*cplComm_size);
    for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
      count_send[i_proc]= 0;
      displ_send[i_proc]=  0;
    }

    void* send_buffer = NULL;
    int tag =recevingField -> fieldIDIntGet();

    MPI_Ialltoallv(send_buffer,count_send,displ_send,MPI_BYTE,
                   loc_v_ptr,count_recv,displ_recv,MPI_BYTE,
                   _cplComm,&request);

    //printf("IRECV %s recevingField -> fieldIDIntGet() %i rank %i both %i request %i\n",
    //recevingField ->fieldIDGet().c_str(),recevingField -> fieldIDIntGet(),cplComm_rank,_both_codes_are_local,request);

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

    int* displ_recv = (int*)malloc(sizeof(int)*cplComm_size);
    int* count_recv = (int*)malloc(sizeof(int)*cplComm_size);


    for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
      count_recv [i_proc]  =  dataTypeSize * nComponent  * ( _targets_localization_idx[ i_proc ][_nb_part_cpl] - _targets_localization_idx[i_proc][0]  );
      displ_recv [i_proc]  =  dataTypeSize * nComponent * _targets_localization_idx[i_proc][0];
    }


    int tag =recevingField -> fieldIDIntGet();

    std::vector<MPI_Request> rreq(localComm_size_cpl,MPI_REQUEST_NULL);
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

  void SpatialInterp::waitIrecv_p2p (Field* recevingField) {
    std::vector<MPI_Status> rstatus(localComm_size_cpl);

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

    for(int i_proc=0; i_proc<localComm_size_cpl;i_proc++) {
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
       && (_cpl -> commTypeGet() == CWP_COMM_PAR_WITH_PART || (_cpl -> commTypeGet() == CWP_COMM_PAR_WITHOUT_PART && cplComm_rank == _senderRank) ) ) {
      _visu -> WriterField(recevingField);
    }
  }

  void SpatialInterp::null_exchange_for_uncoupled_process() {

    int* displ_send = (int*)malloc(sizeof(int)*cplComm_size);
    int* count_send = (int*)malloc(sizeof(int)*cplComm_size);
    for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
      count_send[i_proc] = 0;
      displ_send[i_proc] = 0;
    }

    int* displ_recv = (int*)malloc(sizeof(int)*cplComm_size);
    int* count_recv = (int*)malloc(sizeof(int)*cplComm_size);
    for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
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

  void SpatialInterp::both_codes_on_the_same_process_exchange_p2p(Field* referenceField,Field* recevingField) {

    /* Sending section */

    int  dataTypeSize       = referenceField -> dataTypeSizeGet();
    int nComponent = referenceField -> nComponentGet();

    int tag =referenceField -> fieldIDIntGet();
    void* dist_v_ptr = NULL;

    void* interpolatedFieldData = interpolate(referenceField);

    dist_v_ptr = interpolatedFieldData;

    int* displ_send = (int*)malloc(sizeof(int)*cplComm_size);
    int* count_send = (int*)malloc(sizeof(int)*cplComm_size);
    for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
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

    int* displ_recv = (int*)malloc(sizeof(int)*cplComm_size);
    int* count_recv = (int*)malloc(sizeof(int)*cplComm_size);
    for (int i_proc=0; i_proc < cplComm_size; i_proc++) {
      count_recv[i_proc]  =  dataTypeSize_recv * nComponent_recv
                             * ( _spatial_interp_cpl -> _targets_localization_idx[ i_proc ][_spatial_interp_cpl -> _nb_part_cpl] - _spatial_interp_cpl -> _targets_localization_idx[i_proc][0]  );
      displ_recv [i_proc]  =  dataTypeSize_recv * nComponent_recv * _spatial_interp_cpl -> _targets_localization_idx[i_proc][0];
    }

    std::vector<MPI_Request> rreq(localComm_size,MPI_REQUEST_NULL);
    std::vector<MPI_Request> sreq(localComm_size_cpl,MPI_REQUEST_NULL);
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

    for (int i_proc = 0; i_proc < cplComm_size; i_proc++)
      free(_process_and_partition_count[i_proc] );
    free(_process_and_partition_count);

    _targets_localization_data_count_send = (int*)malloc(sizeof(int)*cplComm_size);
    _targets_localization_data_disp_send  = (int*)malloc(sizeof(int)*cplComm_size);

    for (int i= 0; i < cplComm_size; i++) {
      _targets_localization_data_count_send[i] = sizeof(target_data)*(_targets_localization_idx[i][_nb_part_cpl] - _targets_localization_idx[i][0]);
      _targets_localization_data_disp_send [i] = sizeof(target_data)*_targets_localization_idx[i][0];
    }

  }

  void SpatialInterp::prepare_data_communication_recv() {

    _transform_to_index(_targets_localization_idx_cpl,cplComm_size,_nb_part);

    if(_targets_localization_data_cpl!=NULL) free(_targets_localization_data_cpl);
    _targets_localization_data_cpl = (target_data*)malloc(sizeof(target_data)*_targets_localization_idx_cpl[cplComm_size-1][_nb_part]);

    _targets_localization_data_count_recv = (int*)malloc(sizeof(int)*cplComm_size);
    _targets_localization_data_disp_recv  = (int*)malloc(sizeof(int)*cplComm_size);

    for (int i= 0; i < cplComm_size; i++) {
      _targets_localization_data_count_recv[i] = sizeof(target_data)*(_targets_localization_idx_cpl[i][_nb_part] - _targets_localization_idx_cpl[i][0]);
      _targets_localization_data_disp_recv [i] = sizeof(target_data)*_targets_localization_idx_cpl[i][0];
    }

  }

  void SpatialInterp::data_communication_send_p2p() {

    std::vector<MPI_Request> sreq(localComm_size_cpl,MPI_REQUEST_NULL);
    std::vector<MPI_Status> sstatus(localComm_size_cpl);

    int tagcode = 155;

    int ind = 0;
    vector<int>::iterator it;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      //printf("Send rank %i *it %i \n",cplComm_rank,*it);
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

  void SpatialInterp::data_communication_recv_p2p() {

    std::vector<MPI_Request> rreq(localComm_size_cpl, MPI_REQUEST_NULL);
    std::vector<MPI_Status> rstatus(localComm_size_cpl);

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

  void SpatialInterp::both_data_communication_p2p() {

    std::vector<MPI_Request> sreq(localComm_size,MPI_REQUEST_NULL);
    std::vector<MPI_Status> sstatus(localComm_size);
    std::vector<MPI_Request> rreq(localComm_size_cpl, MPI_REQUEST_NULL);
    std::vector<MPI_Status> rstatus(localComm_size_cpl);

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

    int* count_send = (int*)malloc(sizeof(int)*cplComm_size);
    int* disp_send  = (int*)malloc(sizeof(int)*cplComm_size);

    for (int i= 0; i < cplComm_size; i++) {
      count_send[i] = 0;
      disp_send [i] = 0;
    }
    void* sbuffer_trash = NULL;

    int* count_recv = (int*)malloc(sizeof(int)*cplComm_size);
    int* disp_recv  = (int*)malloc(sizeof(int)*cplComm_size);

    for (int i= 0; i < cplComm_size; i++) {
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
    _n_tot_target_cpl    = _targets_localization_idx_cpl[cplComm_size-1][_nb_part];
  }

  void SpatialInterp::computeFree(){

    for (int i_proc = 0; i_proc < cplComm_size; i_proc++) {
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

  void SpatialInterp::data_index_communication_send_p2p() {

    int* sbuffer = (int*) malloc(sizeof(int)*cplComm_size*_nb_part_cpl);
    for(int i_proc=0;i_proc<cplComm_size;i_proc++)
      for(int i_part=0;i_part<_nb_part_cpl;i_part++) {
        sbuffer[ i_proc * _nb_part_cpl + i_part ] = _process_and_partition_count[i_proc][i_part];
      }


    std::vector<MPI_Request> sreq(localComm_size_cpl);
    std::vector<MPI_Status> sstatus(localComm_size_cpl);

    int tagcode = 152;

    int ind = 0;
    vector<int>::iterator it;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      MPI_Issend(&(sbuffer[*it * _nb_part_cpl]),_nb_part_cpl,MPI_INT,
                 *it, tagcode ,
                 _cplComm, &(sreq[ind]) );
      //printf("Send rank %i ind %i req %i tag %i it %i\n",cplComm_rank,ind,&(sreq[ind]),tagcode + *it,*it);
    }

    ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      MPI_Wait( &(sreq[ind]), &(sstatus[ind]) );
      //printf("Wait rank %i ind %i req %i tag %i\n",cplComm_rank,ind,&(sreq[ind]),tagcode + *it);
    }

    //printf("After Wait Send rank %i \n",cplComm_rank);

    free(sbuffer   );
  }

  void SpatialInterp::data_index_communication_recv_p2p() {

    int* recvbuffer = (int*) malloc(sizeof(int)*cplComm_size*_nb_part);
    for(int i_proc=0;i_proc<cplComm_size;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        recvbuffer[ i_proc * _nb_part + i_part ] = 0;
      }


    std::vector<MPI_Request> rreq(localComm_size_cpl);
    std::vector<MPI_Status> rstatus(localComm_size_cpl);

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

    //printf("After Wait Recv rank %i \n",cplComm_rank);

    for(int i_proc=0;i_proc<cplComm_size;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        _targets_localization_idx_cpl[i_proc][i_part] = recvbuffer[ i_proc * _nb_part + i_part ];
      }

    free(recvbuffer   );
  }

  void SpatialInterp::both_index_communication_p2p() {

    int* sbuffer = (int*) malloc(sizeof(int)*cplComm_size * _nb_part );
    for(int i_proc=0;i_proc<cplComm_size;i_proc++)
      for(int i_part=0;i_part< _nb_part; i_part++) {
        sbuffer[ i_proc * _nb_part + i_part ] = (_spatial_interp_cpl -> _process_and_partition_count)[i_proc][i_part];
      }

    int* recvbuffer = (int*) malloc(sizeof(int)*cplComm_size*_nb_part);
    for(int i_proc=0;i_proc<cplComm_size;i_proc++)
      for(int i_part=0;i_part<_nb_part;i_part++) {
        recvbuffer[ i_proc * _nb_part + i_part ] = 0;
      }

    std::vector<MPI_Request> sreq(localComm_size);
    std::vector<MPI_Status> sstatus(localComm_size);

    std::vector<MPI_Request> rreq(localComm_size_cpl);
    std::vector<MPI_Status> rstatus(localComm_size_cpl);

    int tagcode = 152;

    vector<int>::iterator it;
    int ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      if(*it != cplComm_rank)
        MPI_Irecv(&(recvbuffer[*it * _nb_part]),_nb_part,MPI_INT,
                  *it, tagcode,
                  _cplComm, &(rreq[ind]) );
      else
        memcpy(&(recvbuffer[*it * _nb_part]), &(sbuffer[*it * _nb_part]), sizeof(int)*_nb_part );

      //       printf("After sendrecv rank %i ind %i it %i %i localComm_size_cpl %i\n",cplComm_rank,ind,*it,rreq[ind],localComm_size_cpl);

    }


    ind = 0;
    for( it = _connectableRanks -> begin(); it!=_connectableRanks -> end(); it++,ind++){
      if(*it != cplComm_rank)
        MPI_Issend(&(sbuffer[*it * _nb_part]),_nb_part,MPI_INT,
                   *it, tagcode ,
                   _cplComm, &(sreq[ind]) );
      else
        memcpy(&(recvbuffer[*it * _nb_part]), &(sbuffer[*it * _nb_part]), sizeof(int)*_nb_part );
    }




    ind = 0;
    for( it = _connectableRanks_cpl -> begin(); it!=_connectableRanks_cpl -> end(); it++,ind++){
      //printf("After rank %i MPI _wait ind %i it %i totot %i %i %i localComm_size_cpl %i\n",cplComm_rank,ind,*it,toto,titi,rreq[ind],localComm_size_cpl);
      if(*it != cplComm_rank)
        MPI_Wait( &(rreq[ind]), &(rstatus[ind]) );
    }

    ind = 0;
    for( it = _connectableRanks -> begin(); it!=_connectableRanks -> end(); it++,ind++){
      if(*it != cplComm_rank)
        MPI_Wait( &(sreq[ind]), &(sstatus[ind]) );
    }


    for(int i_proc=0;i_proc<cplComm_size;i_proc++)
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

      sbuffer = (int*) malloc(sizeof(int)*cplComm_size*_nb_part);
      for(int i_proc=0;i_proc<cplComm_size;i_proc++)
        for(int i_part=0;i_part<_nb_part;i_part++)
          sbuffer[ i_proc * _nb_part + i_part ] = 0;

      recvbuffer = (int*) malloc(sizeof(int)*cplComm_size*_nb_part);

      MPI_Ialltoall(sbuffer, _nb_part, MPI_INT,
                    recvbuffer, _nb_part, MPI_INT,
                    _cplComm,&request);

    }
    else if(_Texch_t == CWP_FIELD_EXCH_RECV) {

      // if(_targets_localization_idx==NULL) {
      _targets_localization_idx   =(int**)malloc(sizeof(int*)*cplComm_size);
      for (int i_proc = 0; i_proc < cplComm_size; i_proc++)
        _targets_localization_idx [i_proc] = NULL;
      //  }
      _process_and_partition_count =(int**)malloc(sizeof(int*)*cplComm_size);

      for (int i_proc = 0; i_proc < cplComm_size; i_proc++) {
        // if(_targets_localization_idx [i_proc] == NULL)
        _targets_localization_idx [i_proc] = (int*) malloc( sizeof(int) * (1+_nb_part_cpl) );
        _process_and_partition_count [i_proc]=(int*)malloc(sizeof(int)*(1+_nb_part_cpl));

        for (int i_part = 0; i_part < _nb_part_cpl+1; i_part++) {
          _targets_localization_idx    [i_proc][i_part] = 0;
          _process_and_partition_count [i_proc][i_part] = 0;
        }
      }



      sbuffer = (int*) malloc(sizeof(int)*cplComm_size*_nb_part_cpl);
      for(int i_proc=0;i_proc<cplComm_size;i_proc++)
        for(int i_part=0;i_part<_nb_part_cpl;i_part++)
          sbuffer[ i_proc * _nb_part_cpl + i_part ] = 0;

      recvbuffer = (int*) malloc(sizeof(int)*cplComm_size*_nb_part_cpl);

      MPI_Ialltoall(sbuffer, _nb_part_cpl, MPI_INT,
                    recvbuffer, _nb_part_cpl, MPI_INT,
                    _cplComm,&request);

    }
    MPI_Status stat;
    MPI_Wait(&request,&stat);

    if(sbuffer!=NULL) free(sbuffer   );
    if(recvbuffer!=NULL) free(recvbuffer   );

  }



/***************************************************************************/
/***************************************************************************/



} // end namespace cwipi

/**
 * \endcond
 */
