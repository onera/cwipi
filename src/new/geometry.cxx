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
#include <pdm_mesh_dist.h>
#include <pdm_gnum.h>
#include <pdm_gnum_location.h>
#include <bftc_error.h>
#include <bftc_printf.h>
#include "cwp.h"



namespace cwipi {

  Geometry::Geometry()
  {
  }

  Geometry::~Geometry()
  {
  }
  
  
    
  void Geometry::init(Coupling *coupling, CWP_Field_value_t geometryLocation) {
    _mesh   = coupling -> meshGet();
    _visu   = coupling -> visuGet();
    _referenceFieldsDB = coupling -> fieldsGet();
    _geometryLocation = geometryLocation; 
    _cpl = coupling;

    _localCodeProperties = _cpl -> localCodePropertiesGet();
    _coupledCodeProperties = _cpl -> coupledCodePropertiesGet();


    _pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(_mesh -> getMPICommP()));
    _nb_part = _mesh -> getNPart();

    _globalComm = _localCodeProperties -> globalCommGet();
    _localComm = _mesh -> getMPIComm();     

    MPI_Comm_rank(_globalComm,&_rank);
    MPI_Comm_rank(_mesh -> getMPIComm(),&_localRank);

    MPI_Comm_size(_globalComm,&_n_ranks_g);
       
    _pdm_globalComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_globalComm));



    localName   = _localCodeProperties -> nameGet();
    coupledName = _coupledCodeProperties -> nameGet();

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
    
    printf("_both_codes_are_local %i rank %i\n",_both_codes_are_local,_rank);
    CouplingDB* cplDB = _cpl -> couplingDBGet();
    string cplId = _cpl -> IdGet();

    if(_both_codes_are_local == 1 && cplDB -> couplingIs(*_coupledCodeProperties,cplId)) {

      Coupling coupling_cpl = cplDB -> couplingGet(*_coupledCodeProperties,cplId);
      printf("coupling -> IdGet() %s\n",_cpl -> IdGet().c_str());
      _geometry_cpl_cell_point = coupling_cpl.geometryGet(CWP_FIELD_VALUE_CELL_POINT);
      _geometry_cpl_cell_point -> _geometry_cpl_cell_point = this;
    }
    
    _connectableRanks_cpl = _coupledCodeProperties -> connectableRanksGet();
    _connectableRanks     = _localCodeProperties   -> connectableRanksGet();
    _n_ranks_cpl = _connectableRanks_cpl->size();
    _n_ranks     = _connectableRanks->size();
         
   
  }
  
  
  void Geometry::mesh_info_get() {

    int tag=0;
    
    _n_vtx    =(int*)malloc(sizeof(int)*_nb_part);  
    _n_elt    =(int*)malloc(sizeof(int)*_nb_part);     

    _n_tot_elt=0;
    _n_tot_vtx=0;   
    for(int i_part =0;i_part<_nb_part;i_part++) {   
      _n_elt[i_part]  = _mesh -> getPartNElts(i_part);
      _n_tot_elt+=_n_elt[i_part];
      _n_vtx[i_part]  = _mesh -> getPartNVertex(i_part);
      _n_tot_vtx+=_n_vtx[i_part];
          printf("HHH5 TEST %i %i\n",i_part,_mesh -> getPartNElts(i_part));
    } 

    _n_g_elt    =(int*)malloc(sizeof(int)*_nb_part);  

    printf("Before MPI_Allreduce rank %i\n",_rank);
    MPI_Allreduce(_n_elt, _n_g_elt, _nb_part, MPI_INT,MPI_SUM,_localComm);
    printf("After MPI_Allreduce rank %i\n",_rank);


   _n_g_elt_over_part = 0;
   for(int i_part =0;i_part<_nb_part;i_part++) { 
     _n_g_elt_over_part+=_n_g_elt[i_part];
   }

    _n_g_vtx    =(int*)malloc(sizeof(int)*_nb_part);  
   
    printf("Before MPI_AllreduceVtx\n");
    MPI_Allreduce(_n_vtx, _n_g_vtx, _nb_part, MPI_INT,MPI_SUM,_localComm);
    printf("After MPI_AllreduceVtx\n");

  _n_g_vtx_over_part = 0;
  for(int i_part =0;i_part<_nb_part;i_part++) { 
    _n_g_vtx_over_part+=_n_g_vtx[i_part];
  }
  
  /* if(_both_codes_are_local==0 || (_both_codes_are_local==1 )){  
  _both_codes_are_local__array = (int*)malloc(sizeof(int)*_n_ranks_g);
   MPI_Alltoall(&_both_codes_are_local,1,MPI_INT,
               _both_codes_are_local__array,1,MPI_INT,
               _globalComm);
               
    /*}
*/
 /*   if(_both_codes_are_local == 1 && localName == _codeVector[1]) {
      _geometry_cpl_cell_point->_both_codes_are_local__array = _both_codes_are_local__array;
    }
*/
 }



/***************************************************************************/
/***************************************************************************/

  void Geometry::compute(int *n_uncomputed_tgt) {
  
   // if() {

    
    if (_geometryLocation == CWP_FIELD_VALUE_CELL_POINT) {

             
      if(_both_codes_are_local == 0){
         mesh_info_get();
         printf("After mesh_info_get() %i\n",_rank);  
         mesh_cpl_info_get();
         printf("After mesh_cpl_info_get() %i\n",_rank);  
        // 
         if(localName == _codeVector[0]) {

           locate_cell_point_setting_surface(_id_dist1);
         }
         if(localName == _codeVector[1]) {
           locate_cell_point_setting_request(_id_dist1);
         }

         locate_cell_point_compute          (_id_dist1)    ; 
         
          MPI_Barrier(_globalComm);
         if (localName == _codeVector[1]) locate_cell_point_get(_id_dist1)  ;
         
         PDM_mesh_dist_free(_id_dist1,1);
         
         printf("After  locate_cell_point_get_cpl 1\n");
         if(localName == _codeVector[1]) {
           redistribution_cell_point_request(_id_gnum_location1);
         }
         
         if(localName == _codeVector[0]) {
           redistribution_cell_point_set    (_id_gnum_location1);
         }

         location_compute                   (_id_gnum_location1);   
         
          if(localName == _codeVector[1]) location_get(_id_gnum_location1) ;

          PDM_gnum_location_free(_id_gnum_location1,1);
           
         if(localName == _codeVector[0]) {
           locate_cell_point_setting_request(_id_dist2);
         }
 
          if(localName == _codeVector[1]) {
           locate_cell_point_setting_surface(_id_dist2);
         }      
         
         locate_cell_point_compute          (_id_dist2)  ;
                     
         if (localName == _codeVector[0]) locate_cell_point_get(_id_dist2)  ;    
         PDM_mesh_dist_free(_id_dist2,1);
                     printf("After  locate_cell_point_get_cpl 2\n");
         if(localName == _codeVector[0]) {
           redistribution_cell_point_request(_id_gnum_location2);
         } 
         if(localName == _codeVector[1]) {
           redistribution_cell_point_set    (_id_gnum_location2);
         }    
       
         location_compute                   (_id_gnum_location2);   
         printf("After location_compute 2 %i\n",_rank);
         MPI_Barrier(_globalComm);
         if(localName == _codeVector[0]) location_get(_id_gnum_location2)      ;
         printf("After location_get 2 %i\n",_rank);
         PDM_gnum_location_free(_id_gnum_location2,1);
         printf("After PDM_gnum_location_free 2 %i\n",_rank);
         MPI_Barrier(_globalComm);
         printf("Before redistribution_cell_point_filling_of_redistribution_array %i\n",_rank);
         redistribution_cell_point_filling_of_redistribution_array(); 
         printf("After redistribution_cell_point_filling_of_redistribution_array %i\n",_rank);
         MPI_Barrier(_globalComm);
         printf("Before redistribution_index_communication rank %i\n",_rank);
         redistribution_index_communication(0) ;
         _Wait()    ;     
             
         MPI_Barrier(_globalComm); 
         printf("After redistribution_index_communication rank %i\n",_rank);
         redistribution_communication() ;

         MPI_Barrier(_globalComm); 
         
         redistribution_communication2() ;//while(1==1){}
         printf("After redistribution_communication rank %i\n",_rank);//while(1==1){}
         MPI_Barrier(_globalComm); 
         redistribution_wait_and_targets_array_filling();
      }
      else {
      
        if(localName == _codeVector[0]) {
            printf("Before mesh_info_get()\n");  
            mesh_info_get();
            _geometry_cpl_cell_point -> mesh_info_get();
            printf("After mesh_info_get() %i\n",_rank);  
            mesh_cpl_info_get();
            _geometry_cpl_cell_point -> mesh_cpl_info_get2();   
            printf("After mesh_cpl_info_get() %i\n",_rank);

            locate_cell_point_setting_surface(_id_dist1);
            locate_cell_point_compute        (_id_dist1);
            MPI_Barrier(_globalComm);
            locate_cell_point_get_cpl        (_id_dist1) ;
            printf("After  locate_cell_point_get_cpl 1\n");
            PDM_mesh_dist_free(_id_dist1,1);
               
            redistribution_cell_point_set    (_id_gnum_location1);
            location_compute                 (_id_gnum_location1);
            printf("After  location_compute 1\n");
            location_get_cpl                 (_id_gnum_location1);
            printf("After  location_get_cpl 1\n");
            PDM_gnum_location_free(_id_gnum_location1,1);

            locate_cell_point_setting_request(_id_dist2);
            printf("After  ocate_cell_point_setting_request2\n");  
            locate_cell_point_compute        (_id_dist2);
            printf("After  locate_cell_point_compute\n");  
       
            locate_cell_point_get            (_id_dist2) ;

            PDM_mesh_dist_free(_id_dist2,1);
            
            printf("After locate_cell_point_get 2\n");
            redistribution_cell_point_request(_id_gnum_location2);
            location_compute                 (_id_gnum_location2)  ;
            MPI_Barrier(_globalComm);
            location_get                     (_id_gnum_location2);
            printf("After location_get 2 %i\n",_rank);
            PDM_gnum_location_free(_id_gnum_location2,1);
            printf("After PDM_gnum_location_free 2 %i\n",_rank);
            MPI_Barrier(_globalComm);
            
            printf("Before redistribution_cell_point_filling_of_redistribution_arra %i\n",_rank); 
            redistribution_cell_point_filling_of_redistribution_array();  
            _geometry_cpl_cell_point -> redistribution_cell_point_filling_of_redistribution_array(); 
            printf("After redistribution_cell_point_filling_of_redistribution_arra %i\n",_rank); 
          
            MPI_Barrier(_globalComm);
            
          printf("Before redistribution_index_communication %i\n",_rank);  
          redistribution_index_communication(0)    ;

          printf("After redistribution_index_communication %i\n",_rank);  
          _geometry_cpl_cell_point -> redistribution_index_communication(1)    ;
          printf("After redistribution_index_communication CPL %i\n",_rank);  

          _Wait()    ;
          printf("After _Wait()    ; %i\n",_rank);  
          _geometry_cpl_cell_point -> _Wait()    ;
         printf("After _geometry_cpl_cell_point -> _Wait()    ; %i\n",_rank);  
          
          MPI_Barrier(_globalComm);
          redistribution_communication()          ;
          printf("After redistribution_communication %i\n",_rank); 

          _geometry_cpl_cell_point -> redistribution_communication()          ;
          printf("After redistribution_communication CPL %i\n",_rank); //while(1==1){}
          MPI_Barrier(_globalComm); 
          redistribution_communication2();
          _geometry_cpl_cell_point -> redistribution_communication2();
          MPI_Barrier(_globalComm); 
          redistribution_wait_and_targets_array_filling();
          _geometry_cpl_cell_point -> redistribution_wait_and_targets_array_filling();

        }

      }
    }
    if (_geometryLocation == CWP_FIELD_VALUE_NODE) {
    //  locate_node();
   //   redistribution_node();
    }
  }



  void Geometry::mesh_cpl_info_get2() {
       
    /*      Partition Number exchange           */

    _nb_part_cpl=-1;
       
    MPI_Status status;
      
    _both_codes_are_local__array = (int*)malloc(sizeof(int)*_n_ranks_g);
 
     printf("Before_both_codes_are_local MPI_INT MPI_Allgather %i %i\n",_rank,_both_codes_are_local);     
  
     _both_codes_are_local__array = _geometry_cpl_cell_point -> _both_codes_are_local__array;
     printf("After MPI_INT MPI_Allgather %i \n",_rank);   

        
    if(_both_codes_are_local == 1) {
      _nb_part_cpl = _geometry_cpl_cell_point->_nb_part;
    }
 
        
    if(_both_codes_are_local == 1) {
      _n_g_elt_cpl_over_part = _geometry_cpl_cell_point->_n_g_elt_over_part;
    }
        
    if(_both_codes_are_local == 1) {
      _n_g_vtx_cpl_over_part = _geometry_cpl_cell_point->_n_g_vtx_over_part;
    }
 
    printf("After info_mesh %i\n",_rank);               
    
  }
  


  void Geometry::mesh_cpl_info_get() {
       
    /*      Partition Number exchange           */

    _nb_part_cpl=-1;
       
    MPI_Status status;
    
    int sender = (*_connectableRanks)[_n_ranks-1];
    int sender_cpl = (*_connectableRanks)[_n_ranks-1];
      
    _both_codes_are_local__array = (int*)malloc(sizeof(int)*_n_ranks_g);
 
     printf("Before MPI_INT MPI_Allgather %i \n",_rank);     
  
     if(_both_codes_are_local == 0 || (_both_codes_are_local == 1 && localName == _codeVector[0]))
       MPI_Allgather(&_both_codes_are_local,1,MPI_INT,
                     _both_codes_are_local__array,1,MPI_INT,
                     _globalComm);
    
     else
       _both_codes_are_local__array = _geometry_cpl_cell_point -> _both_codes_are_local__array;
     printf("After MPI_INT MPI_Allgather %i \n",_rank);   
//while(1==1){}
     int senderRank=0;
     while( _both_codes_are_local__array[ (*_connectableRanks)[senderRank] ] == 1 && senderRank < _n_ranks) {
       senderRank++;
     } 
     if (senderRank == _n_ranks) senderRank--;
     senderRank = (*_connectableRanks)[senderRank ];

     int senderRank_cpl=0;
     while( _both_codes_are_local__array[ (*_connectableRanks_cpl)[senderRank_cpl] ] == 1 && senderRank_cpl < _n_ranks_cpl ) {
       senderRank_cpl++;
     } 
     if (senderRank_cpl == _n_ranks_cpl) senderRank_cpl--;
     senderRank_cpl = (*_connectableRanks_cpl)[senderRank_cpl] ;
            
     printf("AfterTTTMPI_INT _IBcast(&_nb_part %i %i\n",senderRank,senderRank_cpl);     
      
     MPI_Request srequest[_n_ranks_cpl];
     MPI_Request request,requests,requestr;
       
     int tag = 1523;
      
     if(_rank == senderRank ){
       for(int i=0;i<_n_ranks_cpl;i++) {
         if( _both_codes_are_local__array[ (*_connectableRanks_cpl)[i] ] == 0 ) {
           MPI_Issend(&_nb_part, 1, MPI_INT,
                      (*_connectableRanks_cpl)[i], tag,
                      _globalComm,&srequest[i]);      
         }
       }
     }

    if(_both_codes_are_local == 0) {
     MPI_Irecv(&_nb_part_cpl, 1, MPI_INT,
               senderRank_cpl, tag,
               _globalComm,&request);   
   }

    if(_both_codes_are_local == 0) 
      MPI_Wait(&request,&status); 

    if(_rank == senderRank ){
      for(int i=0;i<_n_ranks_cpl;i++) {        
        if( _both_codes_are_local__array[ (*_connectableRanks_cpl)[i] ] == 0 ) 
          MPI_Wait(&srequest[i],&status);            
      }
    }
        
    if(_both_codes_are_local == 1) {
      _nb_part_cpl = _geometry_cpl_cell_point->_nb_part;
    }
 

    tag++;
    /*   Number of elements over all processes and partitions exchange                  */
     if(_rank == senderRank ){
       for(int i=0;i<_n_ranks_cpl;i++) {
         if( _both_codes_are_local__array[ (*_connectableRanks_cpl)[i] ] == 0 ) {
           MPI_Issend(&_n_g_elt_over_part, 1, MPI_INT,
                      (*_connectableRanks_cpl)[i], tag,
                      _globalComm,&srequest[i]);      
         }
       }
     }

    if(_both_codes_are_local == 0) {
     MPI_Irecv(&_n_g_elt_cpl_over_part, 1, MPI_INT,
               senderRank_cpl, tag,
               _globalComm,&request);   
   }

    if(_both_codes_are_local == 0) 
      MPI_Wait(&request,&status); 

    if(_rank == senderRank ){
      for(int i=0;i<_n_ranks_cpl;i++) {        
        if( _both_codes_are_local__array[ (*_connectableRanks_cpl)[i] ] == 0 ) 
          MPI_Wait(&srequest[i],&status);            
      }
    }
        
    if(_both_codes_are_local == 1) {
      _n_g_elt_cpl_over_part = _geometry_cpl_cell_point->_n_g_elt_over_part;
    }
 
      /*   Number of elements over all processes and partitions exchange                  */
    tag++;
      if(_rank == senderRank ){
       for(int i=0;i<_n_ranks_cpl;i++) {
         if( _both_codes_are_local__array[ (*_connectableRanks_cpl)[i] ] == 0 ) {
           MPI_Issend(&_n_g_vtx_over_part, 1, MPI_INT,
                      (*_connectableRanks_cpl)[i], tag,
                      _globalComm,&srequest[i]);      
         }
       }
     }

    if(_both_codes_are_local == 0) {
     MPI_Irecv(&_n_g_vtx_cpl_over_part, 1, MPI_INT,
               senderRank_cpl, tag,
               _globalComm,&request);   
   }

    if(_both_codes_are_local == 0) 
      MPI_Wait(&request,&status); 

    if(_rank == senderRank ){
      for(int i=0;i<_n_ranks_cpl;i++) {        
        if( _both_codes_are_local__array[ (*_connectableRanks_cpl)[i] ] == 0 ) 
          MPI_Wait(&srequest[i],&status);            
      }
    }
        
    if(_both_codes_are_local == 1) {
      _n_g_vtx_cpl_over_part = _geometry_cpl_cell_point->_n_g_vtx_over_part;
    }
 
    printf("After info_mesh %i\n",_rank);           
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
        
            target_data_elt* sendptr = (target_data_elt*)send_buffer;
            MPI_Issend(&(  sendptr [ send_disp[distant_rank] / sizeof(target_data_elt) ] ), send_count[distant_rank] * 1, type, distant_rank, tagsend,
                   comm,
                   &(_send_requests[i_rank]));        
 
           target_data_elt* recvptr = (target_data_elt*)recv_buffer;
           
           MPI_Irecv(&(  recvptr [ recv_disp[distant_rank] / sizeof(target_data_elt) ] ), recv_count[distant_rank] * 1,type, distant_rank, tagrecv,
                  comm,
                  &(_recv_requests[i_rank])); 
        }   

      }//end for on i_rank

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
      int ind_proc = 0;
 
      
      std::vector<int> send_requests(nranks,0);
      std::vector<int> recv_requests(nranks,0);
            
      for(int i_rank=0;i_rank<nranks;i_rank++) {
        if( rank == rootRank){ 
          int distant_rank = connectableranks[i_rank];
          printf("rrrr %i\n",distant_rank);
          if(distant_rank != rootRank) {
            MPI_Issend(send_buffer, send_stride * send_size, type, distant_rank, tag,
                   comm,
                   &(send_requests[i_rank]));
                   
            printf("testibcast rank %i send i_rank %i distant_rank %i\n",rank,i_rank,distant_rank);     
          }
        }
      }//end for on i_rank

      int distant_rank = rootRank;
      
      if(type == MPI_DOUBLE){
        double* ptr = (double*)recv_buffer;
        if(rootRank!=rank) 
          MPI_Irecv(&(  ptr[ind_proc] ), recv_stride * recv_size,type, rootRank, tag,
                  comm,
                  &(recv_requests[0]));  
      }   

      if(type == MPI_LONG_LONG_INT){
        CWP_g_num_t* ptr = (CWP_g_num_t*)recv_buffer;
        if(rootRank!=rank) 
          MPI_Irecv(&(  ptr[ind_proc] ), recv_stride * recv_size,type, rootRank, tag,
                  comm,
                 &(recv_requests[0]));  
      }   

      if(type == MPI_INT){
        int* ptr = (int*)recv_buffer;
        if(rank != rootRank  ) {
          
          MPI_Irecv(&(  ptr[ind_proc] ), recv_stride * recv_size,type, rootRank, tag,
                  comm,
                  &(recv_requests[0]));  
         printf("testibcast rank %i recv distant_rank %i\n",rank,distant_rank);         
         }
      }   
       

      if( rank == rootRank ){ 
        for(int i_rank=0;i_rank<nranks;i_rank++) {
        
          int distant_rank = connectableranks[i_rank];
          printf("testibcast waitsend rootRank %i rank %i i_rank %i distant_rank %i\n",rootRank,rank,i_rank,distant_rank);
      
         if(distant_rank != rootRank) 
            MPI_Wait(&(send_requests[i_rank]), &status);
        }
      }

     if(rank != rootRank) 
        MPI_Wait(&(recv_requests[0]), &status);
  }


  void Geometry::_IAlltoall(void* send_buffer,
                int* send_size,
                int send_stride,
                void* recv_buffer,
                int* recv_size,
                int recv_stride,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableRanks_cpl,
                int nranks){
      MPI_Status status;
      int tag = 0;
      int ind_proc = 0;
      int ind_proc_send = 0;
            
      std::vector<int> send_requests(nranks,0);
      std::vector<int> recv_requests(nranks,0);
      
      if(recv_size==NULL) {
        recv_size = (int*) malloc(sizeof(int)*nranks);
        for(int i=0;i<nranks;i++){
          recv_size[i]=1;
        } 
      }
      
      int inc_send = 1;
      if(send_size==NULL) {
        inc_send=0;
        send_size = (int*) malloc(sizeof(int)*nranks);
        for(int i=0;i<nranks;i++){
          send_size[i]=1;
        } 
      }
      
      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = connectableRanks_cpl[i_rank];

        if(type == MPI_DOUBLE){
        
          double* ptr_send = (double*)send_buffer;
          MPI_Issend(&( ptr_send[ind_proc_send] ), send_stride * send_size[i_rank], type, distant_rank, tag,
                   comm,
                   &(send_requests[i_rank]));
        
          double* ptr = (double*)recv_buffer;
        
          MPI_Irecv(&(  ptr[ind_proc] ), recv_stride * recv_size[i_rank],type, distant_rank, tag,
                  comm,
                  &(recv_requests[i_rank]));  
        }   

        if(type == MPI_LONG_LONG_INT){
        
          CWP_g_num_t* ptr_send = (CWP_g_num_t*)send_buffer;
          MPI_Issend(&( ptr_send[ind_proc_send] ), send_stride * send_size[i_rank], type, distant_rank, tag,
                   comm,
                   &(send_requests[i_rank]));
        
        
          CWP_g_num_t* ptr = (CWP_g_num_t*)recv_buffer;
        
          MPI_Irecv(&(  ptr[ind_proc] ), recv_stride * recv_size[i_rank],type, distant_rank, tag,
                  comm,
                  &(recv_requests[i_rank]));  
        }   

        if(type == MPI_INT){
        
          int* ptr_send = (int*)send_buffer;
          MPI_Issend( &( ptr_send[ind_proc_send] ) , send_stride * send_size[i_rank], type, distant_rank, tag,
                   comm,
                   &(send_requests[i_rank]));
        
          int* ptr = (int*)recv_buffer;
        
          MPI_Irecv(&(  ptr[ind_proc] ), recv_stride * recv_size[i_rank],type, distant_rank, tag,
                  comm,
                  &(recv_requests[i_rank]));  
        }   


       if(inc_send!=0) {      
         ind_proc_send += send_stride * send_size[i_rank];       
       }
       ind_proc += recv_stride * recv_size[i_rank];
      }//end for on i_rank

      for(int i_rank=0;i_rank<nranks;i_rank++) {
      
        int distant_rank = connectableRanks_cpl[i_rank];
        if(i_rank !=distant_rank) {
          MPI_Wait(&(send_requests[i_rank]), &status);
          MPI_Wait(&(recv_requests[i_rank]), &status);
        }
      }//end for on i_rank 
  }



/***************************************************************************/
/***************************************************************************/

  void Geometry::irecv(Field<double> *recevingField) {

    _idx_elt  .resize   (_nb_part + 1);
    _idx_elt[0] = 0;
    for (int i_part = 0; i_part < _nb_part; i_part++) {
      _idx_elt[i_part+1] = _idx_elt[i_part] + _n_elt[i_part];   
    }

      //Crée un buffer de réception et le stocke (alloue)
      recevingField -> ReceptionBufferCreation(_idx_elt,_n_tot_elt);
      /* Loop on possibly intersecting distant ranks */
      /*---------------------------------------------*/

      //Réception des données sur toutes les partitions
      double* data = recevingField -> recvBufferGet();

      CWP_Type_t dataType = CWP_DOUBLE;

      int size = sizeof(dataType);
      int nComponent = recevingField -> nComponentGet();

      for (int i_proc = 0; i_proc < _n_ranks_cpl; i_proc++) {

        int tag =0;


        int distant_rank = (*_connectableRanks_cpl)[i_proc];
        void* loc_v_ptr = &(data[ nComponent*_targets_cpl_idx_cpl[distant_rank][0] ]);

        MPI_Request request;
        
        

        int longueur = nComponent * ( _targets_cpl_idx_cpl[ distant_rank ][_nb_part] - _targets_cpl_idx_cpl[distant_rank][0]  );
        printf("Recv from %i to %i start longueur %i\n",_rank,i_proc,longueur);

        MPI_Irecv(loc_v_ptr, longueur, MPI_DOUBLE, distant_rank, tag,
                  _globalComm,
                  &request);

        recevingField -> lastRequestAdd (i_proc,request);
       // printf("recv request %i\n",request);
          
     } //end for

  }

/******************************************************/

  void Geometry::waitIrecv (Field<double>* recevingField) {

    MPI_Status status;

    //TODO: A travailler pour optimiser le temps de communication
    // Qui communique avec qui ? Dans quel ordre ?
    for (int i_proc=0; i_proc < _n_ranks_cpl; i_proc++) {
      int request = recevingField -> lastRequestGet(i_proc);

     printf("request %i %i\n",request,i_proc);

      MPI_Wait(&request, &status);
      
    } //i_proc loop

    //Récupère un pointeur vers le bloc de données reçues
    double* recvData = recevingField -> recvBufferGet();
    int nComponent = recevingField -> nComponentGet();
    CWP_Field_value_t   recevingFieldType = recevingField -> typeGet        ();

    //Reorganize by partition datas which are organized by sending processp
    
    double** userDataMem = (double**)malloc(sizeof(double*)*_nb_part);
    for (int i_part=0;i_part<_nb_part;i_part++) {
       userDataMem[i_part] = recevingField -> dataGet(i_part);
       if(userDataMem[i_part] == NULL ) PDM_error(__FILE__, __LINE__, 0, "Reception memory has not been allocated.\n");
   }

   for(int i_proc=0; i_proc<_n_ranks_cpl;i_proc++) {
     int distant_rank = (*_connectableRanks_cpl)[i_proc];
     if (recevingFieldType == CWP_FIELD_VALUE_CELL_POINT) {
       for (int itarget = _targets_cpl_idx_cpl[ distant_rank ][0]; itarget < _targets_cpl_idx_cpl[ distant_rank ][_nb_part_cpl]; itarget++) {  
         // Index in the interpolated Data array
         int interpInd = itarget;  
             
         //Index of the corresponding local reference Data.
         int iel = _targets_cpl_cpl[itarget].l_num_origin ;
         int lpart = _targets_cpl_cpl[itarget].inc ;
         
         for (int k = 0; k < nComponent; k++) {
              userDataMem[lpart][ nComponent * iel + k ] = recvData[ nComponent * interpInd + k  ];
         }//loop on k
       }// loop on itarget
     } // if referenceFieldType == CWP_FIELD_VALUE_CELL_POINT   
          
         if (recevingFieldType == CWP_FIELD_VALUE_NODE) {

           
         } // if referenceFieldType == CWP_FIELD_VALUE_NODE           
  }// loop on proc


    if(_visu -> isCreated()) {
        printf("_visu -> WriterField(recevingField);\n");
       _visu -> WriterField(recevingField);
    }
     //  
  }


/******************************************************/

  void Geometry::waitIssend (Field<double>* sendingField) {

    MPI_Status status;

    //TODO: A travailler pour optimiser le temps de communication
    // Qui communique avec qui ? Dans quel ordre ?

    for (int i_proc=0; i_proc < _n_ranks_cpl; i_proc++) {
      
      int request = sendingField -> lastRequestGet(i_proc);

      //On attend la fin de l'échange
      MPI_Wait(&request, &status);

    } //i_proc loop
//while(1==1){}
    int nComponent = sendingField -> nComponentGet();
//while(1==1){}
    if(_visu -> isCreated()) {
       printf("_visu -> WriterField(sendingField);\n");    
       _visu -> WriterField(sendingField);
    }

  }


} // end namespace cwipi



