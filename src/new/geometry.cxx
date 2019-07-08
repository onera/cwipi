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
     _n_vtx(NULL),
     _n_elt(NULL),
     _gnum_target(NULL),
     _coords_target(NULL),
     _n_g_elt(NULL),
     _n_g_vtx(NULL),
     _both_codes_are_local__array(NULL)
  {
  }


  Geometry::~Geometry()
  {
    free(_n_vtx);
    free(_n_elt);
    free(_gnum_target);
    free(_coords_target);
    free(_n_g_elt);
    free(_n_g_vtx);
    free(_both_codes_are_local__array);
    
    computeFree();
  }
  
  
    
  void Geometry::init(Coupling *coupling, CWP_Field_value_t geometryLocation) {
    _mesh   = coupling -> meshGet();
    _visu   = coupling -> visuGet();
    _referenceFieldsDB = coupling -> fieldsGet();
    _geometryLocation = geometryLocation; 
    _cpl = coupling;
    _id_dist1 = -1;
    _localCodeProperties = _cpl -> localCodePropertiesGet();
    _coupledCodeProperties = _cpl -> coupledCodePropertiesGet();


    _pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(_mesh -> getMPICommP()));
    _nb_part = _mesh -> getNPart();

    _globalComm = _localCodeProperties -> globalCommGet();
    _localComm = _mesh -> getMPIComm();     
    _connectableComm =_localCodeProperties -> connectableCommGet();

    MPI_Comm_size(_globalComm,&_n_ranks_g);
       
    _pdm_globalComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_globalComm));

    localName   = _localCodeProperties -> nameGet();
    coupledName = _coupledCodeProperties -> nameGet();
    
    _senderRank     = _localCodeProperties   -> rootRankGet();
    _senderRank_cpl = _coupledCodeProperties -> rootRankGet();


    _connectableRanks_cpl = _coupledCodeProperties -> connectableRanksGet();
    _connectableRanks     = _localCodeProperties   -> connectableRanksGet();
    _n_ranks_cpl = _connectableRanks_cpl->size();
    _n_ranks     = _connectableRanks->size();
   
   
    _intraRanks     =  _localCodeProperties -> intraRanksGet();
    _intraRanks_cpl =  _coupledCodeProperties -> intraRanksGet();
    _n_ranks_intra_cpl = _intraRanks_cpl->size();
    _n_ranks_intra     = _intraRanks->size();
    
  
    MPI_Comm_rank(_globalComm,&_rank);
    
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
 
 
    _both_codes_are_local__array = (int*)malloc(sizeof(int)*_n_ranks_g);
    
    int* _both_codes_are_local__array_tmp = (int*)malloc(sizeof(int)*_n_ranks);
    
    CouplingDB* cplDB = _cpl -> couplingDBGet();
    string cplId = _cpl -> IdGet();

    if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && cplDB -> couplingIs(*_coupledCodeProperties,cplId) ) ){
      MPI_Allgather(&_both_codes_are_local,1,MPI_INT,
                   _both_codes_are_local__array,1,MPI_INT,
                   _globalComm );
                   
      if(_rank == _senderRank) {
        int tagsend = 2;
        int tagrecv = 2;
        MPI_Status status;
        MPI_Sendrecv(&_nb_part,1,MPI_INT,_senderRank_cpl,tagsend,&_nb_part_cpl,1,MPI_INT,_senderRank_cpl,tagrecv,_globalComm,&status); 
      }
    }
    
    if(_both_codes_are_local == 1 && cplDB -> couplingIs(*_coupledCodeProperties,cplId)) {

      Coupling coupling_cpl = cplDB -> couplingGet(*_coupledCodeProperties,cplId);
      _geometry_cpl = coupling_cpl.geometryGet(_geometryLocation );
      _geometry_cpl -> _geometry_cpl = this;
    }
    
    int id     = _localCodeProperties -> idGet();
    int id_cpl = _coupledCodeProperties -> idGet();

    int tmp1,tmp2;
    if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && cplDB -> couplingIs(*_coupledCodeProperties,cplId) ) ){
      if(id < id_cpl) {
        MPI_Bcast(&_nb_part_cpl,1,MPI_INT,_senderRank,_globalComm);
        MPI_Bcast(&tmp1,1,MPI_INT,_senderRank_cpl,_globalComm); 
      }
      else{
        MPI_Bcast(&tmp1,1,MPI_INT,_senderRank_cpl,_globalComm); 
        MPI_Bcast(&_nb_part_cpl,1,MPI_INT,_senderRank,_globalComm); 
      }

    }
    else {
      _nb_part_cpl = _geometry_cpl -> _nb_part;
    }
    
    for(int i=0;i<_n_ranks_g;i++) {
   //    printf("rank %i _both_codes_are_local__array[%i] %i %i\n",_rank,i,_both_codes_are_local__array[i],_both_codes_are_local);
    }

    _isCoupledRank     = _localCodeProperties   -> isCoupledRank();
    _isCoupledRank_cpl = _coupledCodeProperties -> isCoupledRank();
   
    n_uncomputed_tgt.resize(_nb_part);
    
    _n_vtx    =(int*)malloc(sizeof(int)*_nb_part);  
    _n_elt    =(int*)malloc(sizeof(int)*_nb_part);     
    _n_target =(int*)malloc(sizeof(int)*_nb_part);
    _gnum_target   =(CWP_g_num_t**)malloc( sizeof(CWP_g_num_t*)*_nb_part);
    _coords_target =(double**)     malloc( sizeof(double*)     *_nb_part);
    
    _n_g_elt    =(CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*_nb_part);  
    _n_g_vtx    =(CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*_nb_part);  


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
      
    if (_geometryLocation == CWP_FIELD_VALUE_NODE) {
      _n_tot_target = _n_tot_vtx;
    }   

    int lsize;
    MPI_Comm_size(_localComm,&lsize);

   // if(_both_codes_are_local == 0 || (_both_codes_are_local == 0 && _Texch_t = CWP_FIELD_EXCH_SEND)){

      PDM_timer_t *t1 = PDM_timer_create();
      PDM_timer_init(t1);
      PDM_timer_resume(t1); 

      CWP_g_num_t n_tot_elt_long = (CWP_g_num_t)_n_tot_elt;
      MPI_Reduce(&n_tot_elt_long,&_n_g_elt_over_part,1,MPI_LONG,MPI_SUM,0,_localComm);
      MPI_Bcast(&_n_g_elt_over_part,1,MPI_LONG,0,_localComm);

      PDM_timer_hang_on(t1);
      double elapsRecv = 0.0;// PDM_timer_elapsed(t1);       

      /*************VTX ***********/
      CWP_g_num_t n_tot_vtx_long = (CWP_g_num_t)_n_tot_vtx;
      MPI_Reduce(&n_tot_vtx_long,&_n_g_vtx_over_part,1,MPI_LONG,MPI_SUM,0,_connectableComm);
      MPI_Bcast(&_n_g_vtx_over_part,1,MPI_LONG,0,_connectableComm);
   /* }
    else {
      _n_g_elt_over_part = 
    
    }*/ 
  }

  void Geometry::mesh_info_get2() {



 /*      PDM_timer_t *t1 = PDM_timer_create();
      PDM_timer_init(t1);
      PDM_timer_resume(t1);       
 */
   /*  _WaitSend(*_connectableRanks, &_send_requests);
    _WaitRecv(*_connectableRanks, &_recv_requests);
*/

    
   /* _WaitSend(senderRankV,&_send_requests);
     if(_rank == _senderRank) _WaitRecv(*_connectableRanks,&_recv_requests);
*/
/*   PDM_timer_hang_on(t1);
   double elapsWaitSend = PDM_timer_elapsed(t1);    
   printf("elapsWaitSend %12.5e\n",elapsWaitSend);  



    std::vector<int> senderRankV(1);
    senderRankV[0]= _senderRank;
 
       
    int tag = -1;
    if(_geometryLocation == CWP_FIELD_VALUE_CELL_POINT)
      tag = 1500;
    else if(_geometryLocation == CWP_FIELD_VALUE_NODE)
      tag = 1510;


     if(_geometryLocation == CWP_FIELD_VALUE_CELL_POINT)
      tag = 1480;
    else if(_geometryLocation == CWP_FIELD_VALUE_NODE)
      tag = 1490;*/
    
    
 /*   _WaitSend(senderRankV,&_send_requests2);
     if(_rank == _senderRank) _WaitRecv(*_connectableRanks,&_recv_requests2);
     
     
    int lsize;
    MPI_Comm_size(_localComm,&lsize);

    _n_g_elt_over_part = 0;
    _n_g_vtx_over_part = 0;
 
 
 
    int** array = (int**) malloc(sizeof(int*) * _n_ranks_g);
    for(int i=0;i<_n_ranks_g;i++){
      array[i] = (int*) malloc(sizeof(int) * 1);
    }

    int** array2 = (int**) malloc(sizeof(int*) * _n_ranks_g);
    for(int i=0;i<_n_ranks_g;i++){
      array2[i] = (int*) malloc(sizeof(int) * 1);
    }


    if(_rank == _senderRank){
    
      for(int i_proc =0;i_proc<lsize;i_proc++) { 
         for(int i_part =0;i_part<_nb_part;i_part++) {
           int disRank = (*_connectableRanks)[i_proc]; 
           _n_g_elt_over_part+=_n_g_elt_tmp[disRank][i_part];
           _n_g_vtx_over_part+=_n_g_vtx_tmp[disRank][i_part];
         }
       }

       for(int i=0;i<_n_ranks_g;i++) array [i][0] = _n_g_elt_over_part;
       for(int i=0;i<_n_ranks_g;i++) array2[i][0] = _n_g_vtx_over_part;
       
      _IAlltoall2Send(array,
                NULL,
                1,
                MPI_INT, 
                _globalComm,
                *_connectableRanks,
                &_send_requests,
                tag
                );
    
      _IAlltoall2Send(array2,
                NULL,
                1,
                MPI_INT, 
                _globalComm,
                *_connectableRanks,
                &_send_requests2,
                tag
                );
    }
    else {
    
    
          _IAlltoall2Recv(array,
                NULL,
                1,
                MPI_INT, 
                _globalComm,
                senderRankV,
                &_recv_requests,
                tag
                );
    
    
      _IAlltoall2Recv(array2,
                NULL,
                1,
                MPI_INT, 
                _globalComm,
                senderRankV,
                &_recv_requests2,
                tag
                );
    }




    if(_rank == _senderRank) _WaitSend(*_connectableRanks,&_send_requests);
    else _WaitRecv(senderRankV,&_recv_requests);

    if(_rank == _senderRank) _WaitSend(*_connectableRanks,&_send_requests2);
    else _WaitRecv(senderRankV,&_recv_requests2);
     */
 /*   free(_n_elt_tmp);
    free(_n_vtx_tmp);
    free(_n_g_elt_tmp);
    free(_n_g_vtx_tmp);*/
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
  void Geometry::info_mesh(CWP_Field_exch_t _Texch_t) {

    double elapsInfoGet,elapsInfoGet2,elapsInfoGetCpl;
    
    if(_both_codes_are_local == 0){
    
      PDM_timer_t *t1 = PDM_timer_create();
      PDM_timer_init(t1);
      PDM_timer_resume(t1); 
    
      if(_isCoupledRank)  mesh_info_get();
      PDM_timer_hang_on(t1);
      elapsInfoGet = PDM_timer_elapsed(t1);

      PDM_timer_init(t1); 
      PDM_timer_resume(t1); 
      
      if(_isCoupledRank) mesh_info_get2();
         PDM_timer_hang_on(t1);
      elapsInfoGet2 = PDM_timer_elapsed(t1);   

      PDM_timer_resume(t1);      
      mesh_cpl_info_get();
      PDM_timer_hang_on(t1);
      elapsInfoGetCpl = PDM_timer_elapsed(t1);  
      
    }
    else if(_Texch_t == CWP_FIELD_EXCH_SEND) {
      mesh_info_get();  
      _geometry_cpl -> mesh_info_get();

       mesh_info_get2();  
      _geometry_cpl -> mesh_info_get2();      

       mesh_cpl_info_get();
      _geometry_cpl -> mesh_cpl_info_get2();
    }

     double elapsInfoGet2Max,elapsInfoGetMax,elapsInfoGetCplMax;
     
    
    if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && _Texch_t == CWP_FIELD_EXCH_SEND) ){ 
     
     MPI_Allreduce(&elapsInfoGet, &elapsInfoGetMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     MPI_Allreduce(&elapsInfoGet2, &elapsInfoGet2Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);    
     MPI_Allreduce(&elapsInfoGetCpl, &elapsInfoGetCplMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     
     if (_rank == 0) {
        printf ("Info Time detail : Get %12.5e Get2 %12.5e GetCpl %12.5e _n_ranks %i \n",
        elapsInfoGetMax,elapsInfoGet2Max,elapsInfoGetCplMax,_n_ranks);
     }
     
     }

   //  while(1==1){}
   //  while(1==1){}
    
    
  /*
    if(_both_codes_are_local == 0){
      mesh_info_get();
      MPI_Barrier(_globalComm);
      mesh_cpl_info_get_send();
      mesh_cpl_info_get_recv();
      MPI_Barrier(_globalComm);  
    }
    else {
      if(_Texch_t == CWP_FIELD_EXCH_SEND) {
        mesh_info_get();
        _geometry_cpl -> mesh_info_get();
        MPI_Barrier(_globalComm);
      mesh_cpl_info_get_send();
      mesh_cpl_info_get_recv();
        _geometry_cpl -> mesh_cpl_info_get2();   
        MPI_Barrier(_globalComm);
      }
    }*/  
  }

/***************************************************************************/
/***************************************************************************/

  void Geometry::compute(CWP_Field_exch_t Texch_t) {
    _Texch_t = Texch_t;
    double elapsDataComm = 0.0;
    double elapsIndexComm =0.0;
    double elapsLocateSetting = 0.0;
    double elapsLocationTime = 0.0;
    double elapsInfo = 0.0;
 
    PDM_timer_t *t1 = PDM_timer_create();
    PDM_timer_init(t1);
    PDM_timer_resume(t1);  
 
    info_mesh(_Texch_t);
    
    PDM_timer_hang_on(t1);
    elapsInfo = PDM_timer_elapsed(t1);

    if(_both_codes_are_local == 0){
      if(_isCoupledRank) {   
      
        PDM_timer_t *t2 = PDM_timer_create();
        PDM_timer_init(t2);
        PDM_timer_resume(t2);  
      
        if(_Texch_t == CWP_FIELD_EXCH_SEND ) locate_setting_surface(&_id_dist1);        
        if(_Texch_t == CWP_FIELD_EXCH_RECV ) locate_setting_request(&_id_dist1);
        PDM_timer_hang_on(t2);
        elapsLocateSetting = PDM_timer_elapsed(t2);
        //printf("Before locate compute\n");
        locate_compute        (_id_dist1); 
        //printf("After locate compute\n");     
        PDM_dist_cloud_surf_dump_times(_id_dist1);

        if(_Texch_t == CWP_FIELD_EXCH_RECV) locate_get(_id_dist1)  ;
        PDM_dist_cloud_surf_free(_id_dist1,1);

        PDM_timer_t *t5 = PDM_timer_create();
        PDM_timer_init(t5);
        PDM_timer_resume(t5);  
         
        if(_Texch_t == CWP_FIELD_EXCH_RECV) broadcasting_request(&_id_gnum_location1);
        if(_Texch_t == CWP_FIELD_EXCH_SEND) broadcasting_set    (&_id_gnum_location1);

        location_compute                   (_id_gnum_location1);   
       
        if(_Texch_t == CWP_FIELD_EXCH_RECV) location_get(_id_gnum_location1) ;

        PDM_timer_hang_on(t5); 
        elapsLocationTime  = PDM_timer_elapsed(t5);

        PDM_gnum_location_free(_id_gnum_location1,1);



        if(_Texch_t == CWP_FIELD_EXCH_RECV) filling_of_broadcasting_array();      
        if(_Texch_t == CWP_FIELD_EXCH_SEND) initialization_of_reception_array();
        
        //MPI_Barrier(_globalComm);
        
        PDM_timer_t *t3 = PDM_timer_create();
        PDM_timer_init(t3);
        PDM_timer_resume(t3);  
        
        if(_Texch_t == CWP_FIELD_EXCH_RECV) broadcasting_index_communication() ;
        if(_Texch_t == CWP_FIELD_EXCH_SEND) reception_index_communication() ;
  
        PDM_timer_hang_on(t3); 
        elapsIndexComm = PDM_timer_elapsed(t3);
        
        PDM_timer_t *t4 = PDM_timer_create();
        PDM_timer_init(t4);
        PDM_timer_resume(t4);  


   
        if(_Texch_t == CWP_FIELD_EXCH_RECV) prepare_data_communication_send();
        if(_Texch_t == CWP_FIELD_EXCH_SEND) prepare_data_communication_recv() ; 

        if(_Texch_t == CWP_FIELD_EXCH_RECV) data_communication_send();
        if(_Texch_t == CWP_FIELD_EXCH_SEND) data_communication_recv();

        if(_Texch_t == CWP_FIELD_EXCH_RECV) data_communication_wait_send();
        if(_Texch_t == CWP_FIELD_EXCH_SEND) data_communication_wait_recv();

        PDM_timer_hang_on(t4);
        elapsDataComm = PDM_timer_elapsed(t4);
        
      }//end if isCoupledRank
      else {
        locate_setting_null(&_id_dist1);
        
        locate_compute     (_id_dist1); 
        PDM_dist_cloud_surf_dump_times(_id_dist1);
        PDM_dist_cloud_surf_free(_id_dist1,1);
        broadcasting_set_null(&_id_gnum_location1);
        location_compute                   (_id_gnum_location1);   
        PDM_gnum_location_free(_id_gnum_location1,1);
      //  MPI_Barrier(_globalComm);
        broadcasting_index_null();
        
        data_communication_null();
        
      }
    }
    else {
      if(_Texch_t == CWP_FIELD_EXCH_SEND) {
        _geometry_cpl -> _Texch_t = CWP_FIELD_EXCH_RECV;
        locate_setting_surface(&_id_dist1);    
        //printf("Before locate compute"); 
        locate_compute        (_id_dist1);  
        //printf("After locate compute");     
	PDM_dist_cloud_surf_dump_times(_id_dist1);

        
        locate_get_cpl        (_id_dist1) ;
        PDM_dist_cloud_surf_free(_id_dist1,1);
              
        broadcasting_set    (&_id_gnum_location1);    
    
        location_compute                 (_id_gnum_location1);     
        location_get_cpl (_id_gnum_location1);

        PDM_gnum_location_free(_id_gnum_location1,1);

        _geometry_cpl -> filling_of_broadcasting_array(); 
        initialization_of_reception_array();

        both_index_communication() ;
        _geometry_cpl ->prepare_data_communication_send();
        prepare_data_communication_recv() ;
 
        _geometry_cpl -> data_communication_send();
        data_communication_recv();
       
        _geometry_cpl -> data_communication_wait_send();
        data_communication_wait_recv();
      }//end if localName == _codeVector[0]
      
    }//end both_are_local

     double elapsDataCommMax,elapsIndexCommMax,elapsLocateSettingMax,elapsLocationTimeMax,elapsInfoMax;

    if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && _Texch_t == CWP_FIELD_EXCH_SEND) ){      
     MPI_Allreduce(&elapsLocateSetting, &elapsLocateSettingMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     MPI_Allreduce(&elapsIndexComm, &elapsIndexCommMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);    
     MPI_Allreduce(&elapsDataComm, &elapsDataCommMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     MPI_Allreduce(&elapsLocationTime, &elapsLocationTimeMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     MPI_Allreduce(&elapsInfo, &elapsInfoMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     
     if (_rank == 0) {
        printf ("Geom Compute Time : LocateSetting %12.5e IndexComm %12.5e DataComm %12.5e Location Time %12.5e InfoMesh %12.5e \n",
        elapsLocateSettingMax,elapsIndexCommMax,elapsDataCommMax,elapsLocationTimeMax,elapsInfoMax);
     }
       
    }    
    
        // while(1==1){}
  }




void Geometry::mesh_cpl_info_get() {
       
    /*      Partition Number exchange           */

     int tsize;
     MPI_Comm_size(_globalComm,&tsize);
    // while(1==1){}
 
 
    PDM_timer_t *t1 = PDM_timer_create();
   PDM_timer_init(t1);
   PDM_timer_resume(t1);  
    
    if(_rank == _senderRank) {
        int tagsend = 2;
        int tagrecv = 2;
        MPI_Status status;
        MPI_Sendrecv(&_n_g_elt_over_part,1,MPI_LONG,_senderRank_cpl,tagsend,&_n_g_elt_cpl_over_part,1,MPI_LONG,_senderRank_cpl,tagrecv,_globalComm,&status); 
        tagsend = 3;
        tagrecv = 3;
        MPI_Sendrecv(&_n_g_vtx_over_part,1,MPI_LONG,_senderRank_cpl,tagsend,&_n_g_vtx_cpl_over_part,1,MPI_LONG,_senderRank_cpl,tagrecv,_globalComm,&status); 
    }
    
    
    PDM_timer_hang_on(t1);
    double elapsAlltoAll = PDM_timer_elapsed(t1);
    double elapsAlltoAllMax = 0.0;   
    MPI_Reduce(&elapsAlltoAll, &elapsAlltoAllMax, 1, MPI_DOUBLE, MPI_MAX,0,_globalComm);
    
    
    int id     = _localCodeProperties -> idGet();
    int id_cpl = _coupledCodeProperties -> idGet();

    CWP_g_num_t tmp1,tmp2;
    
   PDM_timer_t *t2 = PDM_timer_create();
   PDM_timer_init(t2);
   PDM_timer_resume(t2);  
    
    if(id < id_cpl) {
        MPI_Bcast(&_n_g_elt_cpl_over_part,1,MPI_LONG,_senderRank,_globalComm);
        MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_globalComm); 
    }
    else{
      MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_globalComm); 
      MPI_Bcast(&_n_g_elt_cpl_over_part,1,MPI_LONG,_senderRank,_globalComm); 
    }
    
    if(id < id_cpl) {
      MPI_Bcast(&_n_g_vtx_cpl_over_part,1,MPI_LONG,_senderRank,_globalComm);
      MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_globalComm); 
    }
    else{
      MPI_Bcast(&tmp1,1,MPI_LONG,_senderRank_cpl,_globalComm); 
      MPI_Bcast(&_n_g_vtx_cpl_over_part,1,MPI_LONG,_senderRank,_globalComm); 
    }
      

 
    PDM_timer_hang_on(t2);
    double elapsbcast = PDM_timer_elapsed(t2);
    double elapsbcastMax = 0.0;   
    MPI_Reduce(&elapsbcast, &elapsbcastMax, 1, MPI_DOUBLE, MPI_MAX,0,_globalComm);
    
    if (_rank==0) {
      printf ("get_cpl details Time : %12.5e %12.5e\n",elapsAlltoAllMax,elapsbcastMax);
    }  
      
      
  }



  void Geometry::mesh_cpl_info_get2() {
       
    /*      Partition Number exchange           */

    if(_both_codes_are_local == 1) {
      _n_g_elt_cpl_over_part = _geometry_cpl->_n_g_elt_over_part;
    }
        
    if(_both_codes_are_local == 1) {
      _n_g_vtx_cpl_over_part = _geometry_cpl->_n_g_vtx_over_part;
    }
  }
  

  void Geometry::mesh_cpl_info_get_send() {
       
    /*      Partition Number exchange           */

    _nb_part_cpl=-1;
       
    MPI_Status status;
    
    int sender = (*_connectableRanks)[_n_ranks-1];
    int sender_cpl = (*_connectableRanks)[_n_ranks-1];
      
     int tsize;

     MPI_Request srequest[_n_ranks_g];
     MPI_Request request,requests;
     
     MPI_Comm_size(_globalComm,&tsize);
     
     int tag = 0;    
     switch(_geometryLocation){
       case CWP_FIELD_VALUE_CELL_POINT:
         tag = 1;
       case CWP_FIELD_VALUE_NODE:
         tag = 2;        
     }
     
     if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && localName == _codeVector[0])){
     
       for(int i=0;i<_n_ranks_g;i++) {
         int distant_rank = i;//(*_connectableRanks_cpl)[i];
         int tag2=tag*10000+distant_rank;
         if(distant_rank != _rank)
           MPI_Issend(&_both_codes_are_local, 1, MPI_INT,
                      distant_rank, tag2,
                      _globalComm,&srequest[i]);      
         
          tag2=tag*10000+_rank;
       }
     }

     if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && localName == _codeVector[0])){
       for(int i=0;i<_n_ranks_cpl;i++) {
         int distant_rank = i;//(*_connectableRanks_cpl)[i];
         if(distant_rank != _rank)
           MPI_Wait(&srequest[i],&status);      
       }
     }
     else
       _both_codes_are_local__array = _geometry_cpl -> _both_codes_are_local__array;
       
     _both_codes_are_local__array[_rank]=_both_codes_are_local;

     MPI_Barrier(_globalComm);
       
     tag+=100;
     if(_rank == _senderRank ){
       for(int i=0;i<_n_ranks_cpl;i++) {
         int distant_rank = (*_connectableRanks_cpl)[i];
         int tag2=tag*10000+distant_rank;
         if( _both_codes_are_local__array[ distant_rank ] == 0 ) {
           MPI_Issend(&_nb_part, 1, MPI_INT,
                      distant_rank, tag2,
                      _globalComm,&srequest[i]);      
         }
       }
     }

    if(_both_codes_are_local == 0) 
      MPI_Wait(&request,&status); 

    if(_rank == _senderRank ){
      for(int i=0;i<_n_ranks_cpl;i++) {        
        if( _both_codes_are_local__array[ (*_connectableRanks_cpl)[i] ] == 0 ) 
          MPI_Wait(&srequest[i],&status);            
      }
    }
        
    MPI_Barrier(_globalComm);
        
    if(_both_codes_are_local == 1) {
      _nb_part_cpl = _geometry_cpl->_nb_part;
    }

    tag+=100;
    /*   Number of elements over all processes and partitions exchange                  */
     if(_rank == _senderRank ){
       for(int i=0;i<_n_ranks_cpl;i++) {
         int distant_rank = (*_connectableRanks_cpl)[i];
         int tag2=tag*10000+distant_rank;
         if( _both_codes_are_local__array[distant_rank] == 0 ) {
           MPI_Issend(&_n_g_elt_over_part, 1, MPI_INT,
                      distant_rank, tag2,
                      _globalComm,&srequest[i]);      
         }
       }
     }

    if(_rank == _senderRank ){
      for(int i=0;i<_n_ranks_cpl;i++) {        
        if( _both_codes_are_local__array[ (*_connectableRanks_cpl)[i] ] == 0 ) 
          MPI_Wait(&srequest[i],&status);            
      }
    }
        MPI_Barrier(_globalComm);
    if(_both_codes_are_local == 1) {
      _n_g_elt_cpl_over_part = _geometry_cpl->_n_g_elt_over_part;
    }
 
      /*   Number of elements over all processes and partitions exchange                  */
     tag+=100;
      if(_rank == _senderRank ){
       for(int i=0;i<_n_ranks_cpl;i++) {
         int distant_rank = (*_connectableRanks_cpl)[i];
         int tag2=tag*10000+distant_rank;
         if( _both_codes_are_local__array[ distant_rank ] == 0 ) {
           MPI_Issend(&_n_g_vtx_over_part, 1, MPI_INT,
                      distant_rank, tag2,
                      _globalComm,&srequest[i]);      
         }
       }
     }

    if(_rank == _senderRank ){
      for(int i=0;i<_n_ranks_cpl;i++) {        
        if( _both_codes_are_local__array[ (*_connectableRanks_cpl)[i] ] == 0 ) 
          MPI_Wait(&srequest[i],&status);            
      }
    }
    
    if(_both_codes_are_local == 1) {
      _n_g_vtx_cpl_over_part = _geometry_cpl->_n_g_vtx_over_part;
    }
  }



  void Geometry::mesh_cpl_info_get_recv() {
       
    /*      Partition Number exchange   */

    _nb_part_cpl=-1;
       
    MPI_Status status;
    
    int sender = (*_connectableRanks)[_n_ranks-1];
    int sender_cpl = (*_connectableRanks)[_n_ranks-1];
      
     int tsize;

     MPI_Request rrequest[_n_ranks_g];
     MPI_Request request,requests,requestr;
     
     MPI_Comm_size(_globalComm,&tsize);
     
     int tag = 0;    
     switch(_geometryLocation){
       case CWP_FIELD_VALUE_CELL_POINT:
         tag = 1;
       case CWP_FIELD_VALUE_NODE:
         tag = 2;        
     }
     
     if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && localName == _codeVector[0])){
     
       for(int i=0;i<_n_ranks_g;i++) {
         int distant_rank = i;//(*_connectableRanks_cpl)[i];
         int tag2=tag*10000+_rank;
         if(distant_rank != _rank)
            MPI_Irecv(&(_both_codes_are_local__array[i]), 1, MPI_INT,
                     distant_rank, tag2,
                     _globalComm,&rrequest[i]);   
       }
     }

     if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && localName == _codeVector[0])){
       for(int i=0;i<_n_ranks_cpl;i++) {
         int distant_rank = i;//(*_connectableRanks_cpl)[i];
         if(distant_rank != _rank)
           MPI_Wait(&rrequest[i],&status);  
       }
     }
     else
       _both_codes_are_local__array = _geometry_cpl -> _both_codes_are_local__array;
       
     _both_codes_are_local__array[_rank]=_both_codes_are_local;

     MPI_Barrier(_globalComm);
       
     tag+=100;

    if(_both_codes_are_local == 0) {
     int tag2=tag*10000+_rank;
     MPI_Irecv(&_nb_part_cpl, 1, MPI_INT,
               _senderRank_cpl, tag2,
               _globalComm,&request);   
   }

    if(_both_codes_are_local == 0) 
      MPI_Wait(&request,&status); 

    MPI_Barrier(_globalComm);
        
    if(_both_codes_are_local == 1) {
      _nb_part_cpl = _geometry_cpl->_nb_part;
    }

    tag+=100;
    /*   Number of elements over all processes and partitions exchange                  */

    if(_both_codes_are_local == 0) {
      int tag2 = tag*10000+_rank;
      MPI_Irecv(&_n_g_elt_cpl_over_part, 1, MPI_INT,
               _senderRank_cpl, tag2,
               _globalComm,&request);   
    }

    if(_both_codes_are_local == 0) 
      MPI_Wait(&request,&status); 

        MPI_Barrier(_globalComm);
    if(_both_codes_are_local == 1) {
      _n_g_elt_cpl_over_part = _geometry_cpl->_n_g_elt_over_part;
    }
 
      /*   Number of elements over all processes and partitions exchange                  */
     tag+=100;

    if(_both_codes_are_local == 0) {
     int tag2=tag*10000+_rank;
     MPI_Irecv(&_n_g_vtx_cpl_over_part, 1, MPI_INT,
               _senderRank_cpl, tag2,
               _globalComm,&request);   
   }

    if(_both_codes_are_local == 0) 
      MPI_Wait(&request,&status); 

    if(_both_codes_are_local == 1) {
      _n_g_vtx_cpl_over_part = _geometry_cpl->_n_g_vtx_over_part;
    }
  }



 void Geometry::_IAlltoall2Send(int** send_buffer,
                int* send_size,
                int send_stride,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableRanks,
                std::vector<int>* send_requests,
                int tag
                ){

      int nranks = connectableRanks.size();
      _send_requests.resize(nranks,0);

      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = connectableRanks[i_rank];

        if(type == MPI_INT){
          int** ptr_send = ( int**)send_buffer;
            MPI_Issend(&(  ptr_send [distant_rank] [0] ), send_stride * 1, type, distant_rank, tag,
                   comm,
                   &(_send_requests[i_rank]));        
        }
  //      printf("rank %i send_requests[%i] %i\n",_rank,i_rank,_send_requests[i_rank]);
      }//end for on i_rank
  }


  void Geometry::_IAlltoall2Recv(int** recv_buffer,
                int* recv_size,
                int recv_stride,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableRanks,
                std::vector<int>* recv_requests,
                int tag
                ){

      int nranks = connectableRanks.size();

      _recv_requests.resize(nranks,0);
        
      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = connectableRanks[i_rank];
        if(type == MPI_INT){
          int** ptr = (int**)recv_buffer;
          MPI_Irecv(&(  ptr [distant_rank] [0] ), recv_stride * 1,type, distant_rank, tag,
                  comm,
                  &(_recv_requests[i_rank])); 
        }   
 //       printf("rank %i recv_requests[%i] %i\n",_rank,i_rank,_recv_requests[i_rank]);
      }//end for on i_rank
  }





  void Geometry::_IAlltoall2(int** send_buffer,
                int* send_size,
                int send_stride,
                int** recv_buffer,
                int* recv_size,
                int recv_stride,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableRanks
                ){

      int rank=-1;
      MPI_Comm_rank(comm,&rank);
     
      int comm_size;
      MPI_Comm_size(comm,&comm_size);
      int nranks = connectableRanks.size();
      if(recv_size==NULL) {
        recv_size = (int*) malloc(sizeof(int)*comm_size);
        for(int i=0;i<comm_size;i++){
          recv_size[i]=1;
        } 
      }

      if(send_size==NULL) {
        send_size = (int*) malloc(sizeof(int)*comm_size);
        for(int i=0;i<comm_size;i++){
          send_size[i]=1;
        } 
      }

      _send_requests.resize(nranks,0);
      _recv_requests.resize(nranks,0);

      int tagsend = -1;
      int tagrecv = -1;
      if(localName == _codeVector[0]) {
        tagsend =0; tagrecv = 1;
      }
      else {
       tagsend =1; tagrecv = 0;
      }

      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = connectableRanks[i_rank];
        
       // printf("rank %i distant_rank %i\n",rank,distant_rank);

        if(type == MPI_INT){
          int** ptr_send = ( int**)send_buffer;
            MPI_Issend(&(  ptr_send [distant_rank] [0] ), send_stride * 1/*send_size[distant_rank]*/, type, distant_rank, tagsend,
                   comm,
                   &(_send_requests[i_rank]));        
   
           int** ptr = (int**)recv_buffer;
 
           MPI_Irecv(&(  ptr [distant_rank] [0] ), recv_stride * 1/*recv_size[distant_rank]*/,type, distant_rank, tagrecv,
                  comm,
                  &(_recv_requests[i_rank])); 
        }   

      }//end for on i_rank
    free(send_size);
    free(recv_size);
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



 void Geometry::_IAlltoallIndex(void* send_buffer,
                int* send_count,
                int* send_disp,
                void* recv_buffer,
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

      _send_requests.resize(nranks,0);
      _recv_requests.resize(nranks,0);

      int tagsend = -1;
      int tagrecv = -1;
      if(localName == _codeVector[0]) {
        tagsend =0; tagrecv = 1;
      }
      else {
       tagsend =1; tagrecv = 0;
      }

      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = connectableRanks[i_rank];
        
       // printf("rank %i distant_rank %i\n",rank,distant_rank);

        if(type == MPI_BYTE){
        
            target_data* sendptr = (target_data*)send_buffer;
            MPI_Issend(&(  sendptr [ send_disp[distant_rank] / sizeof(int/*target_data*/) ] ), send_count[distant_rank] * 1, type, distant_rank, tagsend,
                   comm,
                   &(_send_requests[i_rank]));        
 
           target_data* recvptr = (target_data*)recv_buffer;
           
           MPI_Irecv(&(  recvptr [ recv_disp[distant_rank] / sizeof(/*target_data*/int) ] ), recv_count[distant_rank] * 1,type, distant_rank, tagrecv,
                  comm,
                  &(_recv_requests[i_rank])); 
        }   

      }//end for on i_rank

  }


  void Geometry::_Wait(std::vector<int> ranks,
                           std::vector<int>* send_requests,
                           std::vector<int>* recv_requests){
                
      MPI_Status status;  
      int nranks = ranks.size();
      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = ranks[i_rank];
     
        if(_rank != distant_rank ) {  
          MPI_Wait(&((*send_requests)[i_rank]), &status); 
          MPI_Wait(&((*recv_requests)[i_rank]), &status); 
        }
      }//end for on i_rank      
      send_requests -> resize(0,0);
      recv_requests -> resize(0,0);
  }


  void Geometry::_WaitSend(std::vector<int> ranks,
                           std::vector<int>* send_requests){
                
      MPI_Status status;  
      int nranks = ranks.size();
      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = ranks[i_rank];
      //  printf("rank %i Wait send_requests[%i] %i\n",_rank,i_rank,_send_requests[i_rank]);     
        if(_rank != distant_rank ) {  
          MPI_Wait(&(_send_requests[i_rank]), &status); 
        }

      }//end for on i_rank      
      _send_requests.resize(0,0);
  }

  void Geometry::_WaitRecv(std::vector<int> ranks,
                           std::vector<int>* recv_requests){
                
      MPI_Status status;  
      int nranks = ranks.size();
       
      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = ranks[i_rank];
     //   printf("rank %i Wait recv_requests[%i] %i\n",_rank,i_rank,_recv_requests[i_rank]);     
        if(_rank != distant_rank ) {  
          MPI_Wait(&(_recv_requests[i_rank]), &status);
        }

      }//end for on i_rank      
      _recv_requests.resize(0,0);
  }

  void Geometry::_Wait(){
                
      MPI_Status status;  
      int nranks = (*_connectableRanks_cpl).size();
       
      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = (*_connectableRanks_cpl)[i_rank];
     
        if(_rank != distant_rank ) {  
          MPI_Wait(&(_send_requests[i_rank]), &status); 
          MPI_Wait(&(_recv_requests[i_rank]), &status);
        }
      }//end for on i_rank      
      _send_requests.resize(0,0);
      _recv_requests.resize(0,0);
  }

  void Geometry::_IBcast(void* send_buffer,
                     int send_size,
                     int send_stride,
                     void* recv_buffer,
                int recv_size,
                int recv_stride,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableranks,
                int nranks,
                int rootRank){
      MPI_Status status;
 
      int rank;
      MPI_Comm_rank(comm,&rank);
      int tag = 0;    
      switch(_geometryLocation){
        case CWP_FIELD_VALUE_CELL_POINT:
          tag = 1;
        case CWP_FIELD_VALUE_NODE:
          tag = 2;        
      }

      int ind_proc = 0;
      
      std::vector<int> send_requests(nranks,0);
      std::vector<int> recv_requests(nranks,0);
            
      for(int i_rank=0;i_rank<nranks;i_rank++) {
        if( rank == rootRank){ 
          int distant_rank = connectableranks[i_rank];
        //  if(distant_rank != rootRank) {
            MPI_Issend(send_buffer, send_stride * send_size, type, distant_rank, tag,
                   comm,
                   &(send_requests[i_rank]));
          }
        //}
      }//end for on i_rank

      int distant_rank = rootRank;
      
      if(type == MPI_DOUBLE){
        double* ptr = (double*)recv_buffer;
      //  if(rootRank!=rank) 
          MPI_Irecv(&(  ptr[ind_proc] ), recv_stride * recv_size,type, rootRank, tag,
                  comm,
                  &(recv_requests[0]));  
      }   

      if(type == MPI_LONG_LONG_INT){
        CWP_g_num_t* ptr = (CWP_g_num_t*)recv_buffer;
      //  if(rootRank!=rank) 
          MPI_Irecv(&(  ptr[ind_proc] ), recv_stride * recv_size,type, rootRank, tag,
                  comm,
                 &(recv_requests[0]));  
      }   

      if(type == MPI_INT){
        int* ptr = (int*)recv_buffer;
     //   if(rank != rootRank  ) {
          
          MPI_Irecv(&(  ptr[ind_proc] ), recv_stride * recv_size,type, rootRank, tag,
                  comm,
                  &(recv_requests[0]));  
      //   }
      }   
       

      if( rank == rootRank ){ 
        for(int i_rank=0;i_rank<nranks;i_rank++) {
        
          int distant_rank = connectableranks[i_rank];
      
       //  if(distant_rank != rootRank) 
            MPI_Wait(&(send_requests[i_rank]), &status);
        }
      }

   //  if(rank != rootRank) 
        MPI_Wait(&(recv_requests[0]), &status);
  }



/***************************************************************************/
/***************************************************************************/

  void Geometry::irecv(Field *recevingField) {
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

    for (int i_proc = 0; i_proc < _n_ranks_cpl; i_proc++) {
      int tag =recevingField -> fieldIDIntGet();
      int distant_rank = (*_connectableRanks_cpl)[i_proc];
      void* loc_v_ptr = data + dataTypeSize * nComponent * _targets_localization_idx[distant_rank][0] ;

      MPI_Request request;

      int longueur =  dataTypeSize * nComponent * ( _targets_localization_idx[ distant_rank ][_nb_part_cpl] - _targets_localization_idx[distant_rank][0]  );
      //printf("Recv from %i to %i start %i longueur %i\n",_rank,i_proc,nComponent*_targets_localization_idx[distant_rank][0],longueur);

      MPI_Irecv(loc_v_ptr, longueur, MPI_BYTE, distant_rank, tag,
                _globalComm,
                &request);

      recevingField -> lastRequestAdd (i_proc,request);
   } //end for
 }

/******************************************************/

  void Geometry::waitIrecv (Field* recevingField) {

    MPI_Status status;

    for (int i_proc=0; i_proc < _n_ranks_cpl; i_proc++) {
      int request = recevingField -> lastRequestGet(i_proc);
      MPI_Wait(&request, &status);
    } //i_proc loop

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

    //TODO: A travailler pour optimiser le temps de communication
    // Qui communique avec qui ? Dans quel ordre ?

    for (int i_proc=0; i_proc < _n_ranks_cpl; i_proc++) {
      int request = sendingField -> lastRequestGet(i_proc);
      MPI_Wait(&request, &status);
    } //i_proc loop

    int nComponent = sendingField -> nComponentGet();
    if(_visu -> isCreated()) {
       _visu -> WriterField(sendingField);
    }
  }


} // end namespace cwipi



