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
#include <geometry.hxx>
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

namespace cwipi {

  Geometry::Geometry()
    :_targets_localization_idx(NULL),
     _targets_localization_data(NULL),
     _targets_localization_data_cpl(NULL),
     _targets_localization_idx_cpl(NULL),
     _gnum_target(NULL),
     _coords_target(NULL),
     _n_vtx(NULL),
     _n_elt(NULL)
  {
  }


  Geometry::~Geometry()
  {
    free(_n_vtx);
    free(_n_elt);
    free(_gnum_target);
    free(_coords_target);
    computeFree();
  }
  
  
    
  void Geometry::init(Coupling *coupling, CWP_Field_value_t geometryLocation, int slave) {
    _mesh   = coupling -> meshGet();
    _visu   = coupling -> visuGet();
    _referenceFieldsDB = coupling -> fieldsGet();
    _geometryLocation = geometryLocation; 
    _cpl = coupling;
    _id_dist = -1;
    _localCodeProperties = _cpl -> localCodePropertiesGet();
    _coupledCodeProperties = _cpl -> coupledCodePropertiesGet();

    _slave = slave;
    _pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(_mesh -> getMPICommP()));
    _nb_part = _mesh -> getNPart();


    _unionComm = _cpl -> communicationGet() -> unionCommGet();
    _globalComm = _localCodeProperties -> globalCommGet();
    _localComm = _mesh -> getMPIComm();     
    _connectableComm =_localCodeProperties -> connectableCommGet();

    MPI_Comm_size(_unionComm,&_n_ranks_g);
       
    _pdm_globalComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_globalComm));
    _pdm_unionComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_unionComm));

    localName   = _localCodeProperties -> nameGet();
    coupledName = _coupledCodeProperties -> nameGet();
    
    _senderRank     = _cpl -> communicationGet() -> unionCommLocCodeRootRanksGet();
    _senderRank_cpl = _cpl -> communicationGet() -> unionCommCplCodeRootRanksGet();
    
    _connectableRanks_cpl = _cpl -> communicationGet() -> unionCommCplRanksGet();
    _connectableRanks     = _cpl -> communicationGet() -> unionCommLocRanksGet();
    _n_ranks_cpl = _connectableRanks_cpl->size();
    _n_ranks     = _connectableRanks->size();
   
   /* for(int i=0;i<_connectableRanks_cpl->size();i++)
      printf("_connectableRanks_cpl %i\n",(*_connectableRanks_cpl)[i]);
   */
   
    MPI_Comm_rank(_unionComm,&_rank);
    
    std::vector<int> connectable = (*_connectableRanks);
    std::vector<int>::iterator it = std::find(connectable.begin(), connectable.end(), _senderRank);
    _senderLocalRank = -1;
    if (it != connectable.end()) {
      _senderLocalRank = std::distance(connectable.begin(), it);
    }
    
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
    if(_localCodeProperties ->localCodeIs() && _coupledCodeProperties ->localCodeIs())
       _both_codes_are_local = 1; 
 
 
    _id     = _localCodeProperties   -> idGet();
    _id_cpl = _coupledCodeProperties -> idGet();

    if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && slave == 0 ) ){
                   
      if(_rank == _senderRank) {
        int tagsend = 2;
        int tagrecv = 2;
        MPI_Status status;
        MPI_Sendrecv(&_nb_part,1,MPI_INT,_senderRank_cpl,tagsend,&_nb_part_cpl,1,MPI_INT,_senderRank_cpl,tagrecv,_unionComm,&status); 
      }
    }

    CouplingDB* cplDB = _cpl -> couplingDBGet();
    string cplId = coupling -> IdGet();     
    if(_both_codes_are_local == 1 && slave == 0) {
      Coupling coupling_cpl = cplDB -> couplingGet(*_coupledCodeProperties,cplId);
      _geometry_cpl = coupling_cpl.geometryGet(_geometryLocation );
      _geometry_cpl -> _geometry_cpl = this;
    }
    
    int tmp1,tmp2;
   
    if(_both_codes_are_local == 0 ){
      if(_id < _id_cpl) {
        MPI_Bcast(&_nb_part_cpl,1,MPI_INT,_senderRank,_unionComm);
        MPI_Bcast(&tmp1,1,MPI_INT,_senderRank_cpl,_unionComm); 
      }
      else{
        MPI_Bcast(&tmp1,1,MPI_INT,_senderRank_cpl,_unionComm); 
        MPI_Bcast(&_nb_part_cpl,1,MPI_INT,_senderRank,_unionComm); 
      }
    }
    else if( slave == 0 ) {
      printf("_senderRank %i _senderRank_cpl %i\n",_senderRank,_senderRank_cpl);
      if(_id < _id_cpl) {
        MPI_Bcast(&_nb_part_cpl,1,MPI_INT,_senderRank,_unionComm);
        MPI_Bcast(&(_geometry_cpl -> _nb_part_cpl), 1,MPI_INT,_senderRank_cpl,_unionComm); 
      }
      else{
        MPI_Bcast(&(_geometry_cpl -> _nb_part_cpl), 1,MPI_INT,_senderRank_cpl,_unionComm); 
        MPI_Bcast(&_nb_part_cpl,1,MPI_INT,_senderRank,_unionComm);
      } 
    }

    _isCoupledRank     = _localCodeProperties   -> isCoupledRank();
    _isCoupledRank_cpl = _coupledCodeProperties -> isCoupledRank();
   
    n_uncomputed_tgt.resize(_nb_part);
 
    _gnum_target   =(CWP_g_num_t**)malloc( sizeof(CWP_g_num_t*)*_nb_part);
    _coords_target =(double**)     malloc( sizeof(double*)     *_nb_part); 
         
    _n_vtx    =(int*)malloc(sizeof(int)*_nb_part);  
    _n_elt    =(int*)malloc(sizeof(int)*_nb_part);     
    _n_target =(int*)malloc(sizeof(int)*_nb_part);
    
  }
  

/***************************************************************************/
/***************************************************************************/

  void Geometry::computeFree(){

   for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
     if(_targets_localization_idx_cpl != NULL) free(_targets_localization_idx_cpl[i_proc]);
     if(_targets_localization_idx     != NULL) free(_targets_localization_idx[i_proc]);
   }
   
   free(_targets_localization_data_cpl);
   if(_targets_localization_idx_cpl != NULL) free(_targets_localization_idx_cpl);
   if(_Texch_t == CWP_FIELD_EXCH_RECV) free(_targets_localization_data);
   if(_targets_localization_idx     != NULL) free(_targets_localization_idx);
  }


/***************************************************************************/
/***************************************************************************/

void Geometry::mesh_cpl_info_get() {
    if(_slave==0){
      /*      Partition Number exchange           */
    
      if(_rank == _senderRank) {
         int tagsend = 2;
         int tagrecv = 2;
         MPI_Status status;
         MPI_Sendrecv(&_n_g_elt_over_part,1,MPI_LONG,_senderRank_cpl,tagsend,&_n_g_elt_cpl_over_part,1,MPI_LONG,_senderRank_cpl,tagrecv,_unionComm,&status); 
         tagsend = 3;
         tagrecv = 3;
         MPI_Sendrecv(&_n_g_vtx_over_part,1,MPI_LONG,_senderRank_cpl,tagsend,&_n_g_vtx_cpl_over_part,1,MPI_LONG,_senderRank_cpl,tagrecv,_unionComm,&status); 
      }

      CWP_g_num_t tmp1,tmp2;
      
      if(_id < _id_cpl) {
         MPI_Bcast(&_n_g_elt_cpl_over_part,1,MPI_LONG,_senderRank,_unionComm);
         MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_unionComm); 
      }
      else{
        MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_unionComm); 
        MPI_Bcast(&_n_g_elt_cpl_over_part,1,MPI_LONG,_senderRank,_unionComm); 
      }
    
      if(_id < _id_cpl) {
        MPI_Bcast(&_n_g_vtx_cpl_over_part,1,MPI_LONG,_senderRank,_unionComm);
        MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_unionComm); 
      }
      else{
        MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_unionComm); 
        MPI_Bcast(&_n_g_vtx_cpl_over_part,1,MPI_LONG,_senderRank,_unionComm); 
      }
    }
    else{
      _n_g_elt_cpl_over_part = _geometry_cpl->_n_g_elt_over_part;
      _n_g_vtx_cpl_over_part = _geometry_cpl->_n_g_vtx_over_part;
    }
  }

/*************************************************/
  
  void Geometry::mesh_info_get() {

    int tag=0;
    
    _n_tot_elt=0;
    _n_tot_vtx=0;   
    for(int i_part =0;i_part<_nb_part;i_part++) {   
      if (_geometryLocation == CWP_FIELD_VALUE_CELL_POINT) {
        _n_target   [i_part]     = _mesh -> getPartNElts(i_part); 
        _gnum_target[i_part]     = _mesh -> GNumEltsGet(i_part);   
        _coords_target [i_part]  = _mesh -> eltCentersGet(i_part);             
      }

      if (_geometryLocation == CWP_FIELD_VALUE_NODE) {
        _n_target      [i_part]  = _mesh -> getPartNVertex (i_part);
        _gnum_target   [i_part]  = _mesh -> getVertexGNum  (i_part);
        _coords_target [i_part]  = _mesh -> getVertexCoords(i_part);        
      }      
    
      _n_elt[i_part]  = _mesh -> getPartNElts(i_part);
      _n_tot_elt+=_n_elt[i_part];
      
      _n_vtx[i_part]  = _mesh -> getPartNVertex(i_part);
      _n_tot_vtx+=_n_vtx[i_part];
    } //end loop on i_part 


    if (_geometryLocation == CWP_FIELD_VALUE_CELL_POINT) {
      _n_tot_target = _n_tot_elt;
    }
    else if (_geometryLocation == CWP_FIELD_VALUE_NODE) {
      _n_tot_target = _n_tot_vtx;
    }   

    /************* Elements ***********/
    CWP_g_num_t n_tot_elt_long = (CWP_g_num_t)_n_tot_elt;
    MPI_Reduce(&n_tot_elt_long,&_n_g_elt_over_part,1,MPI_LONG,MPI_SUM,0,_localComm);
    MPI_Bcast(&_n_g_elt_over_part,1,MPI_LONG,0,_localComm);

    /************* Vertices ***********/
    CWP_g_num_t n_tot_vtx_long = (CWP_g_num_t)_n_tot_vtx;
    MPI_Reduce(&n_tot_vtx_long,&_n_g_vtx_over_part,1,MPI_LONG,MPI_SUM,0,_connectableComm);
    MPI_Bcast(&_n_g_vtx_over_part,1,MPI_LONG,0,_connectableComm);
  }


/***************************************************************************/
/***************************************************************************/
  void Geometry::info_mesh(CWP_Field_exch_t _Texch_t) {

    if(_both_codes_are_local == 0){
      if(_isCoupledRank)  mesh_info_get();
      mesh_cpl_info_get();
    }
    else if(_Texch_t == CWP_FIELD_EXCH_SEND) {
      mesh_info_get();  
      _geometry_cpl -> mesh_info_get();

       mesh_cpl_info_get();
      _geometry_cpl -> mesh_cpl_info_get();
    }
  }

/***************************************************************************/
/***************************************************************************/

  void Geometry::compute(CWP_Field_exch_t Texch_t) {
    _Texch_t = Texch_t;

    /* Get informations about the local and the coupled meshes */
    info_mesh(_Texch_t);

    if(_both_codes_are_local == 0){
      if(_isCoupledRank) {   
      
        /*********************************/
        /*         Localization         **/
        /*********************************/      
      
        /*Surface and cloud points localization setting */
        if(_Texch_t == CWP_FIELD_EXCH_SEND ) localization_surface_setting(&_id_dist);        
        if(_Texch_t == CWP_FIELD_EXCH_RECV ) localization_points_cloud_setting(&_id_dist);

        /* Localization compute, get and free*/
        localization_compute        (_id_dist); 
        if(_Texch_t == CWP_FIELD_EXCH_RECV) localization_get(_id_dist)  ;
        PDM_dist_cloud_surf_free(_id_dist,1);

        /*********************************/
        /*  Communication tree building **/
        /*********************************/ 
        
        /* From a global number obtained the MPI rank and mesh partition of the element */
        /* Setting and request */                
        if(_Texch_t == CWP_FIELD_EXCH_RECV) broadcasting_request(&_id_gnum_location);
        if(_Texch_t == CWP_FIELD_EXCH_SEND) broadcasting_set    (&_id_gnum_location);
        /* Compute, get and free*/ 
        location_compute                   (_id_gnum_location);   
        if(_Texch_t == CWP_FIELD_EXCH_RECV) location_get(_id_gnum_location) ;
        PDM_gnum_location_free(_id_gnum_location,1);

        /* Initialization */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) filling_of_broadcasting_array();      
        if(_Texch_t == CWP_FIELD_EXCH_SEND) initialization_of_reception_array();
        
        /*  Communication  of the communication tree index */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) broadcasting_index_communication() ;
        if(_Texch_t == CWP_FIELD_EXCH_SEND) reception_index_communication() ;

        /*  Communication of the communication tree */
        /* Preparation */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) prepare_data_communication_send();
        if(_Texch_t == CWP_FIELD_EXCH_SEND) prepare_data_communication_recv() ; 
        /* MPI asynchronous communication */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) data_communication_send();
        if(_Texch_t == CWP_FIELD_EXCH_SEND) data_communication_recv();
        /* MPI Wait */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) data_communication_wait_send();
        if(_Texch_t == CWP_FIELD_EXCH_SEND) data_communication_wait_recv();

      }//end if isCoupledRank
      else {
        /*************************************************/
        /*  Localization for uncoupled ranks processes  **/
        /*************************************************/     
        localization_null_setting(&_id_dist);
        localization_compute     (_id_dist); 
        PDM_dist_cloud_surf_free(_id_dist,1);
        
        /***************************************************************/
        /*  Communication tree building for uncoupled ranks processes **/
        /***************************************************************/ 
        broadcasting_set_null(&_id_gnum_location);
        location_compute                   (_id_gnum_location);   
        PDM_gnum_location_free(_id_gnum_location,1);

        /*  Communication of the communication tree index for uncoupled ranks processes*/        
        broadcasting_index_null();
        /*  Communication of the communication tree for uncoupled ranks processes*/
        data_communication_null();
      }
    }
    else {
      if(_Texch_t == CWP_FIELD_EXCH_SEND) {
        _geometry_cpl -> _Texch_t = CWP_FIELD_EXCH_RECV;
        localization_surface_setting(&_id_dist);    
        localization_compute        (_id_dist);  

        localization_get_cpl        (_id_dist) ;
        PDM_dist_cloud_surf_free(_id_dist,1);
              
        broadcasting_set    (&_id_gnum_location);    
    
        location_compute                 (_id_gnum_location);     
        location_get_cpl (_id_gnum_location);

        PDM_gnum_location_free(_id_gnum_location,1);

        _geometry_cpl -> filling_of_broadcasting_array(); 
        initialization_of_reception_array();
        both_index_communication() ;

        _geometry_cpl ->prepare_data_communication_send();
        prepare_data_communication_recv() ;
        both_data_communication();
        _geometry_cpl -> data_communication_wait_send();
        data_communication_wait_recv();
      }//end if localName == _codeVector[0]
      
    }//end both_are_local

  }




  void Geometry::_IAlltoallIndexSend(void* send_buffer,
                                     int* send_count,
                                     int* send_disp,
                                     MPI_Datatype type, 
                                     MPI_Comm comm,
                                     std::vector<int> connectableRanks
                                    ){
      int rank=-1;
     
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




 void Geometry::_IAlltoallIndexRecv(void* recv_buffer,
                int* recv_count,
                int* recv_disp,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableRanks
                ){

      int rank=-1;
     
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

  void Geometry::irecv2(Field *recevingField) {
    _idx_target  .resize   (_nb_part + 1);
    _idx_target[0] = 0;
    for (int i_part = 0; i_part < _nb_part; i_part++) {
      _idx_target[i_part+1] = _idx_target[i_part] + _n_target[i_part];   
    }


    int  dataTypeSize       = recevingField -> dataTypeSizeGet(); 
    //Crée un buffer de réception et le stocke (alloue)
    recevingField -> ReceptionBufferCreation(_idx_target,_n_tot_target);
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
      count_recv[i_proc]  =  dataTypeSize * nComponent * ( _targets_localization_idx[ i_proc ][_nb_part_cpl] - _targets_localization_idx[i_proc][0]  );
      displ_recv [i_proc]  =  dataTypeSize * nComponent * _targets_localization_idx[i_proc][0];
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
                   _unionComm,&request);
    
    free(count_recv);
    free(displ_recv);
    free(count_send);
    free(displ_send);    
    recevingField -> lastRequestAdd (tag,request);
 }


/******************************************************/

  void Geometry::waitIrecv (Field* recevingField) {
    MPI_Status status;

    if(_both_codes_are_local == 0) {
      int tag = recevingField -> fieldIDIntGet();
      int request = recevingField -> lastRequestGet(tag);    
      MPI_Wait(&request, &status);
    }

    //Récupère un pointeur vers le bloc de données reçues
    void*              recvData          = recevingField -> recvBufferGet  ();
    int                nComponent        = recevingField -> nComponentGet  ();
    int                dataTypeSize      = recevingField -> dataTypeSizeGet(); 
    CWP_Field_value_t  recevingFieldType = recevingField -> typeGet        ();

    //Reorganize by partition datas which are organized by sending processp
    std::vector<void*> userDataMem (_nb_part,NULL);
    for (int i_part=0;i_part<_nb_part;i_part++) {
       userDataMem [i_part] = recevingField -> dataGet(i_part);
       if(userDataMem[i_part] == NULL ) PDM_error(__FILE__, __LINE__, 0, "Reception memory has not been allocated.\n");
       n_uncomputed_tgt[i_part]=0;
    }

    for(int i_proc=0; i_proc<_n_ranks_cpl;i_proc++) {
      int distant_rank = (*_connectableRanks_cpl)[i_proc];
      for (int itarget = _targets_localization_idx[ distant_rank ][0]; itarget < _targets_localization_idx[ distant_rank ][_nb_part_cpl]; itarget++) {  
        // Index in the interpolated Data array
        int interpInd = itarget;  
        int iel = _targets_localization_data[itarget].l_num_origin ;
        int lpart = _targets_localization_data[itarget].origin_part ;       
        if(_targets_localization_data[itarget].distance != INFINITY) {
          //Index of the corresponding local reference Data.
          for (int k = 0; k < nComponent; k++) {
            memcpy(userDataMem[lpart] + dataTypeSize * ( nComponent * iel + k ) ,
                    recvData + dataTypeSize * ( nComponent * interpInd + k ),
                    dataTypeSize);
          }//loop on k
        }
        else {
          n_uncomputed_tgt[lpart]++;
          for (int k = 0; k < nComponent; k++) {
            memcpy(userDataMem[lpart] + dataTypeSize * ( nComponent * iel + k ) ,
                   recvData + dataTypeSize * ( nComponent * interpInd + k ),
                   dataTypeSize);
           *( (double*) (userDataMem[lpart] + dataTypeSize * ( nComponent * iel + k ) ) ) = -1.0;
          }//loop on k
        }
      }// loop on itarget
    }// loop on proc

    if(_visu -> isCreated()) {
      _visu -> WriterField(recevingField);
    }
  }

/******************************************************/

  void Geometry::waitIssend (Field* sendingField) {

    MPI_Status status;
    int tag;
    MPI_Request request;
    
    tag     = sendingField -> fieldIDIntGet();
    request = sendingField -> lastRequestGet(tag);
    
    MPI_Wait(&request, &status);

    int nComponent = sendingField -> nComponentGet();
    if(_visu -> isCreated()) {
       _visu -> WriterField(sendingField);
    }
  }


} // end namespace cwipi



