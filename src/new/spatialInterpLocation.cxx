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
#include <spatialInterpLocation.hxx>
#include <spatialInterp.hxx>
#include <mpi.h>
#include <pdm_mpi.h>
#include <pdm_mesh_nodal.h>
#include <pdm_dist_cloud_surf.h>
#include <pdm_mesh_location.h>
#include <pdm_gnum.h>
#include <pdm_gnum_location.h>
#include <pdm_geom_elem.h>
#include <bftc_error.h>
#include <bftc_printf.h>
#include "cwp.h"
#include <limits>
#include <algorithm>
#include <cmath>
#include <pdm_timer.h>

/**
 * \cond
 */

namespace cwipi {

  void _transform_to_index(int** array,int l1, int l2) {

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

  void _transform_to_index(int* array,int l1) {

    int sav = array[0];
    array[0]=0;
    array [l1]=0;
    for (int i_part = 1; i_part < l1+1; i_part++) {
      int sav2 = array[i_part];
      array[i_part] = sav + array[i_part-1];
      sav = sav2;
    }
  }



  SpatialInterpLocation::SpatialInterpLocation()
  :SpatialInterp::SpatialInterp(),
   _targets_localization_idx(NULL),
   _targets_localization_data(NULL),
   _targets_localization_idx_cpl(NULL),
   _targets_localization_data_cpl(NULL),
   _gnum_target(NULL),
   _coords_target(NULL),
   _n_vtx(NULL),
   _n_elt(NULL),
   _pdmGNum_handle_index(-1),
   _weights_src_idx(NULL),
   _weights_src(NULL)
  {
  }




  SpatialInterpLocation::~SpatialInterpLocation()
  {
    free(_n_vtx);
    free(_n_elt);
    free(_gnum_target);
    free(_coords_target);

    free(_n_user_targets);
    free(_gnum_user_targets);
    free(_coords_user_targets);

    if (_weights_src_idx != NULL) {
      for (int i = 0; i < _nb_part; i++) {
        free (_weights_src_idx[i]);
      }
      free (_weights_src_idx);
    }

    if (_weights_src != NULL) {
      for (int i = 0; i < _nb_part; i++) {
        free (_weights_src[i]);
      }
      free (_weights_src);
    }
    computeFree();
  }



  typedef enum {
    CWP_LOCATION_DIST_CLOUD_SURF,
    CWP_LOCATION_MESH_LOCATION_OCTREE,
    CWP_LOCATION_MESH_LOCATION_DBBTREE
  } CWP_Location_method_t;

  const double CWP_MESH_LOCATION_BBOX_TOLERANCE = 1.e-3;

static CWP_Location_method_t _get_location_method
  (
   )
  {
    CWP_Location_method_t location_method = CWP_LOCATION_DIST_CLOUD_SURF; // default location method
    char *env_location_method;
    env_location_method = getenv ("CWP_LOCATION_METHOD");
    if (env_location_method != NULL) {

      int i_location_method = atoi(env_location_method);
      if (i_location_method == 0) {
        location_method = CWP_LOCATION_DIST_CLOUD_SURF;
      }
      else if (i_location_method == 1) {
        location_method = CWP_LOCATION_MESH_LOCATION_OCTREE;
      }
      else if (i_location_method == 2) {
        location_method = CWP_LOCATION_MESH_LOCATION_DBBTREE;
      }

      if (0) {
        printf("CWP_LOCATION_METHOD = %d\n", location_method);
      }
    }

    return location_method;
  }

  static double _get_location_tolerance
  (
   )
  {
    double tolerance = CWP_MESH_LOCATION_BBOX_TOLERANCE;
    char *env_location_tolerance;
    env_location_tolerance = getenv ("CWP_LOCATION_TOLERANCE");

    if (env_location_tolerance != NULL) {
      tolerance = atof(env_location_tolerance);
    }

    return tolerance;
  }


  void SpatialInterpLocation::spatialInterpWeightsCompute(CWP_Field_exch_t Texch_t) {
    _Texch_t = Texch_t;
    /*
      In case of withOutPart the user provided not null data only on the root rank (senderRank).
    */

    if(_both_codes_are_local == 0 ){
      if(_Texch_t == CWP_FIELD_EXCH_RECV && _pointsCloudLocation == CWP_DOF_LOCATION_USER)
        user_targets_gnum_compute();
    }
    else {
      if(_Texch_t == CWP_FIELD_EXCH_SEND) {
        _spatial_interp_cpl -> _Texch_t =  CWP_FIELD_EXCH_RECV;
        if(_isCoupledRank && _pointsCloudLocation == CWP_DOF_LOCATION_USER)
          _spatial_interp_cpl -> user_targets_gnum_compute();
      }
    }
    /* Get informations about the local and the coupled meshes */
    info_mesh();
    if( (_both_codes_are_local == 0 && _cpl -> commTypeGet() == CWP_COMM_PAR_WITH_PART    && _isCoupledRank)
     || (_both_codes_are_local == 0 && _cpl -> commTypeGet() == CWP_COMM_PAR_WITHOUT_PART && _isCoupledRank && _rank == _senderRank) ){
        /*********************************/
        /*         Localization         **/
        /*********************************/

        /*Surface and cloud points localization setting */
        if(_Texch_t == CWP_FIELD_EXCH_SEND ) localization_surface_setting     (&_id_dist);
        if(_Texch_t == CWP_FIELD_EXCH_RECV ) localization_points_cloud_setting(&_id_dist);

      /* Localization compute, get and free*/
      localization_compute        (_id_dist);
      if(_Texch_t == CWP_FIELD_EXCH_RECV) localization_get(_id_dist)  ;
      //PDM_dist_cloud_surf_free(_id_dist,1);
      CWP_Location_method_t location_method = _get_location_method();
        if (location_method == CWP_LOCATION_DIST_CLOUD_SURF) {
          PDM_dist_cloud_surf_free (_id_dist, 1);
        }
        else if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE ||
                   location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {
          PDM_mesh_location_free (_id_dist, 1);
        }
        else {
          PDM_error (__FILE__, __LINE__, 0, "Unknown location method.\n");
        }


        /*********************************/
        /*  Communication tree building **/
        /*********************************/

        /* From a global number obtained the MPI rank and mesh partition of the element */
        /* Setting and request */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) triplet_location_request(&_id_gnum_location);
        if(_Texch_t == CWP_FIELD_EXCH_SEND) triplet_location_set    (&_id_gnum_location);
        /* Compute, get and free*/
        triplet_location_compute                   (_id_gnum_location);
        if(_Texch_t == CWP_FIELD_EXCH_RECV) triplet_location_get(_id_gnum_location) ;
        PDM_gnum_location_free(_id_gnum_location,1);

        /* Initialization */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) filling_of_sending_communication_tree_array();
        /* targets_localization_idx_cpl allocation and init */
        if(_Texch_t == CWP_FIELD_EXCH_SEND) initialization_of_receving_communication_tree_array();

        /*  Communication  of the communication tree index */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) data_index_communication_send_p2p() ;
        if(_Texch_t == CWP_FIELD_EXCH_SEND) data_index_communication_recv_p2p() ;

        /*  Communication of the communication tree */
        /* Preparation */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) prepare_data_communication_send();
        if(_Texch_t == CWP_FIELD_EXCH_SEND) prepare_data_communication_recv() ;

        /* MPI asynchronous communication */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) data_communication_send_p2p();
        if(_Texch_t == CWP_FIELD_EXCH_SEND) data_communication_recv_p2p();

        /* MPI Wait */
        if(_Texch_t == CWP_FIELD_EXCH_RECV) data_communication_wait_send();
        if(_Texch_t == CWP_FIELD_EXCH_SEND) data_communication_wait_recv();

      }
      else if( _both_codes_are_local == 1 && _cpl -> commTypeGet() == CWP_COMM_PAR_WITH_PART && _isCoupledRank ) {
        if(_Texch_t == CWP_FIELD_EXCH_SEND) {
           _spatial_interp_cpl -> _Texch_t = CWP_FIELD_EXCH_RECV;
          localization_surface_setting(&_id_dist);
          localization_compute        (_id_dist);

        localization_get_cpl        (_id_dist) ;
        //PDM_dist_cloud_surf_free(_id_dist,1);
        CWP_Location_method_t location_method = _get_location_method();
        if (location_method == CWP_LOCATION_DIST_CLOUD_SURF) {
          PDM_dist_cloud_surf_free (_id_dist, 1);
        }
        else if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE ||
                   location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {
          PDM_mesh_location_free (_id_dist, 1);
        }
        else {
          PDM_error (__FILE__, __LINE__, 0, "Unknown location method.\n");
        }

          triplet_location_set    (&_id_gnum_location);

          triplet_location_compute  (_id_gnum_location);
          triplet_location_get_cpl (_id_gnum_location);

          PDM_gnum_location_free(_id_gnum_location,1);

          _spatial_interp_cpl -> filling_of_sending_communication_tree_array();
          initialization_of_receving_communication_tree_array();

          both_index_communication_p2p() ;

          _spatial_interp_cpl ->prepare_data_communication_send();
          prepare_data_communication_recv() ;


          both_data_communication_p2p();

          _spatial_interp_cpl -> data_communication_wait_send();
          data_communication_wait_recv();

        }
     }
     else if(_both_codes_are_local == 0 && _cpl -> commTypeGet() == CWP_COMM_PAR_WITH_PART && !_isCoupledRank ) {
        /*************************************************/
        /*  Localization for uncoupled ranks processes  **/
        /*************************************************/
        if(_Texch_t == CWP_FIELD_EXCH_SEND) localization_null_setting_send(&_id_dist);
        if(_Texch_t == CWP_FIELD_EXCH_RECV) localization_null_setting_recv(&_id_dist);

        localization_compute     (_id_dist);

      //PDM_dist_cloud_surf_free(_id_dist,1);
      CWP_Location_method_t location_method = _get_location_method();
      if (location_method == CWP_LOCATION_DIST_CLOUD_SURF) {
        PDM_dist_cloud_surf_free (_id_dist, 1);
      }
      else if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE ||
               location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {
        PDM_mesh_location_free (_id_dist, 1);
      }
      else {
        PDM_error (__FILE__, __LINE__, 0, "Unknown location method.\n");
      }

        /***************************************************************/
        /*  Communication tree building for uncoupled ranks processes **/
        /***************************************************************/

        if(_Texch_t == CWP_FIELD_EXCH_SEND) initialization_of_receving_communication_tree_array();

        if(_Texch_t == CWP_FIELD_EXCH_RECV) triplet_location_null_recv(&_id_gnum_location);
        if(_Texch_t == CWP_FIELD_EXCH_SEND) triplet_location_null_send(&_id_gnum_location);

        triplet_location_compute                   (_id_gnum_location);
        PDM_gnum_location_free(_id_gnum_location,1);

        /*  Communication of the communication tree index for uncoupled ranks processes*/
        data_index_communication_null();
        /*  Communication of the communication tree for uncoupled ranks processes*/
        data_communication_null();
      }
  }







void SpatialInterpLocation::issend(Field* referenceField) {

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


void SpatialInterpLocation::issend_p2p(Field* referenceField) {

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



  void SpatialInterpLocation::user_target_points_set(int i_part, int n_pts, double* coord) {
    if( !(_pointsCloudLocation == CWP_DOF_LOCATION_USER ) )
      PDM_error(__FILE__, __LINE__, 0, "You cannot use user_target_points_set for CWP_Dof_location_t different of CWP_DOF_LOCATION_USER.\n");
    else {
      _n_user_targets     [i_part] = n_pts;
      _coords_user_targets[i_part] = coord;

     /* printf("_n_user_targets [%i] %i _coords_user_targets[i_part][0] %f\n",i_part,_n_user_targets[i_part],_coords_user_targets[i_part][0]);
        for(int i=0;i<3*_n_user_targets     [i_part];i++)
        printf("coords_user_targets[%i][%i] %f\n",i_part,i,_coords_user_targets[i_part][i]);
     */
    }
  }







  void SpatialInterpLocation::waitIssend (Field* sendingField) {

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

  void SpatialInterpLocation::waitIssend_p2p (Field* sendingField) {

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




  void SpatialInterpLocation::irecv(Field *recevingField) {
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


 void SpatialInterpLocation::irecv_p2p(Field *recevingField) {
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







  void SpatialInterpLocation::waitIrecv (Field* recevingField) {
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




  void SpatialInterpLocation::waitIrecv_p2p (Field* recevingField) {
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










void SpatialInterpLocation::null_exchange_for_uncoupled_process() {

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

  void SpatialInterpLocation::null_exchange_for_uncoupled_process_p2p() {

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








  void SpatialInterpLocation::both_codes_on_the_same_process_exchange(Field* referenceField,Field* recevingField) {

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

  void SpatialInterpLocation::both_codes_on_the_same_process_exchange_p2p(Field* referenceField,Field* recevingField) {

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





  void SpatialInterpLocation::init(Coupling *coupling, CWP_Dof_location_t pointsCloudLocation, int slave) {
    _mesh   = coupling -> meshGet();
    _visu   = coupling -> visuGet();
    _referenceFieldsDB = coupling -> fieldsGet();
    _pointsCloudLocation = pointsCloudLocation;
    _cpl = coupling;
    _id_dist = -1;
    _localCodeProperties = _cpl -> localCodePropertiesGet();
    _coupledCodeProperties = _cpl -> coupledCodePropertiesGet();

    _slave = slave;
    _pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(_mesh -> getMPICommP()));
    _nb_part = _mesh -> getNPart();


    _cplComm = _cpl -> communicationGet() -> cplCommGet();
    _globalComm = _localCodeProperties -> globalCommGet();
    _localComm = _mesh -> getMPIComm();
    _connectableComm =_localCodeProperties -> connectableCommGet();


    printf("_cplComm : %d\n", _cplComm);
    MPI_Comm_size(_cplComm,&_n_ranks_g);

    _pdm_globalComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_globalComm));
    _pdm_cplComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_cplComm));
    _pdm_connectableComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_connectableComm));

    localName   = _localCodeProperties -> nameGet();
    coupledName = _coupledCodeProperties -> nameGet();

    _senderRank     = _cpl -> communicationGet() -> unionCommLocCodeRootRanksGet();
    _senderRank_cpl = _cpl -> communicationGet() -> unionCommCplCodeRootRanksGet();

    _connectableRanks_cpl = _cpl -> communicationGet() -> cplCommCplRanksGet();
    _connectableRanks     = _cpl -> communicationGet() -> cplCommLocRanksGet();
    _n_ranks_cpl = _connectableRanks_cpl->size();
    _n_ranks     = _connectableRanks->size();

   /* for(int i=0;i<_connectableRanks_cpl->size();i++)
      printf("_connectableRanks_cpl %i\n",(*_connectableRanks_cpl)[i]);
   */

    MPI_Comm_rank(_cplComm,&_rank);


    MPI_Group globalGroup,unionGroup;
    MPI_Comm_group(_globalComm, &globalGroup);

    MPI_Group intraGroup = _localCodeProperties -> connectableGroupGet();
    MPI_Comm_group(_cplComm, &unionGroup);


    MPI_Group_translate_ranks(unionGroup, 1, &_senderRank,
                              intraGroup ,    &_senderLocalRank);

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

      if(_rank == _senderRank) {
        int tagsend = 2;
        int tagrecv = 2;
        MPI_Status status;
        MPI_Sendrecv(&_nb_part,1,MPI_INT,_senderRank_cpl,tagsend,&_nb_part_cpl,1,MPI_INT,_senderRank_cpl,tagrecv,_cplComm,&status);
      }
    }

    CouplingDB* cplDB = _cpl -> couplingDBGet();
    string cplId = coupling -> IdGet();
    if(_both_codes_are_local == 1 && slave == 0) {
      Coupling coupling_cpl = cplDB -> couplingGet(*_coupledCodeProperties,cplId);
      _spatial_interp_cpl = dynamic_cast<SpatialInterpLocation*>( coupling_cpl.spatialInterpGet(_pointsCloudLocation ) );
      _spatial_interp_cpl -> _spatial_interp_cpl = this;
    }

    int tmp1;
    //printf("rank %i _senderRank %i _senderRank_cpl %i\n",_rank,_senderRank, _senderRank_cpl);

    if(_both_codes_are_local == 0 ){
      if(_id < _id_cpl) {
        MPI_Bcast(&_nb_part_cpl,1,MPI_INT,_senderRank,_cplComm);
        MPI_Bcast(&tmp1,1,MPI_INT,_senderRank_cpl,_cplComm);
      }
      else{
        MPI_Bcast(&tmp1,1,MPI_INT,_senderRank_cpl,_cplComm);
        MPI_Bcast(&_nb_part_cpl,1,MPI_INT,_senderRank,_cplComm);
      }
    }
    else if( slave == 0 ) {
      //printf("_senderRank %i _senderRank_cpl %i\n",_senderRank,_senderRank_cpl);
      if(_id < _id_cpl) {
        MPI_Bcast(&_nb_part_cpl,1,MPI_INT,_senderRank,_cplComm);
        MPI_Bcast(&(_spatial_interp_cpl -> _nb_part_cpl), 1,MPI_INT,_senderRank_cpl,_cplComm);
      }
      else{
        MPI_Bcast(&(_spatial_interp_cpl -> _nb_part_cpl), 1,MPI_INT,_senderRank_cpl,_cplComm);
        MPI_Bcast(&_nb_part_cpl,1,MPI_INT,_senderRank,_cplComm);
      }
    }




    _isCoupledRank     = _localCodeProperties   -> isCoupledRank();
    _isCoupledRank_cpl = _coupledCodeProperties -> isCoupledRank();

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




  void SpatialInterpLocation::mesh_info_get() {

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

    MPI_Barrier(_connectableComm);

    if(_cpl -> commTypeGet() == CWP_COMM_PAR_WITH_PART){
      /************* Elements ***********/
      CWP_g_num_t n_tot_elt_long = (CWP_g_num_t)_n_tot_elt;
      MPI_Reduce(&n_tot_elt_long,&_n_g_elt_over_part,1,MPI_LONG,MPI_SUM,0,_connectableComm);
      MPI_Bcast(&_n_g_elt_over_part,1,MPI_LONG,0,_connectableComm);

      /************* Vertices ***********/
      CWP_g_num_t n_tot_vtx_long = (CWP_g_num_t)_n_tot_vtx;
      MPI_Reduce(&n_tot_vtx_long,&_n_g_vtx_over_part,1,MPI_LONG,MPI_SUM,0,_connectableComm);
      MPI_Bcast(&_n_g_vtx_over_part,1,MPI_LONG,0,_connectableComm);
    }
    else {
      _n_g_vtx_over_part = (CWP_g_num_t)_n_tot_vtx;
      _n_g_elt_over_part = (CWP_g_num_t)_n_tot_elt;
    }

  }







void SpatialInterpLocation::mesh_cpl_info_get() {

    if(_slave==0){
      /*      Partition Number exchange           */

      if(_rank == _senderRank) {
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






  void SpatialInterpLocation::info_mesh() {

    if(_both_codes_are_local == 0){
      if(_isCoupledRank)  mesh_info_get();
      mesh_cpl_info_get();
    }
    else if(_Texch_t == CWP_FIELD_EXCH_SEND) {
       mesh_info_get();
      _spatial_interp_cpl -> mesh_info_get();

       mesh_cpl_info_get();
      _spatial_interp_cpl -> mesh_cpl_info_get();
    }
  }



      /***********************************************************
       ***********************************************************
       **                                                       **
       **            Localization object functions              **
       **                                                       **
       ***********************************************************
       ***********************************************************/


  void SpatialInterpLocation::localization_points_cloud_setting(int* id_dist) {

    CWP_Location_method_t location_method = _get_location_method();

    if (location_method == CWP_LOCATION_DIST_CLOUD_SURF) {
      /* Paradigm mesh localisation _distance creation */
      *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_cplComm );

      PDM_dist_cloud_surf_n_part_cloud_set(*id_dist,   0, _nb_part);

      PDM_dist_cloud_surf_surf_mesh_global_data_set (*id_dist,
                                                     _n_g_elt_cpl_over_part,
                                                     _n_g_vtx_cpl_over_part,
                                                     _nb_part_cpl);

      for (int i_part = 0; i_part < _nb_part; i_part++) {

        CWP_g_num_t* gnum_target          = _gnum_target  [i_part];
        double*      coords_target        = _coords_target[i_part];

        PDM_dist_cloud_surf_cloud_set (*id_dist,
                                       0,
                                       i_part,
                                       _n_target[i_part],
                                       coords_target,
                                       gnum_target
                                       );
      }

      if(_both_codes_are_local == 0) {
        for(int i_part =0; i_part<_nb_part_cpl; i_part++) {
          int n_elt_null = 0;
          int*         connecIdx = (int*)malloc(sizeof(int)*(1+n_elt_null));
          int*         connec    = (int*)malloc(sizeof(int)*n_elt_null);

          double*      coords    = (double*)malloc(3*sizeof(double)*n_elt_null);
          CWP_g_num_t* gnum_vtx  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*n_elt_null);
          CWP_g_num_t* gnum_elt  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*n_elt_null);

          connecIdx[0]=0;

          PDM_dist_cloud_surf_surf_mesh_part_set (*id_dist,
                                                  i_part,
                                                  0,
                                                  connecIdx,
                                                  connec,
                                                  gnum_elt,
                                                  0,
                                                  coords,
                                                  gnum_vtx);
        }
      }
      else {
        for(int i_part =0; i_part<_nb_part_cpl; i_part++) {
          Mesh* mesh_cpl = _spatial_interp_cpl -> _mesh;
          int*         connecIdx_cpl = mesh_cpl -> connecIdxGet(i_part);
          int*         connec_cpl    = mesh_cpl -> connecGet(i_part);

          int          n_vtx_cpl     = mesh_cpl -> getPartNVertex(i_part);
          int          n_elts_cpl    = mesh_cpl -> getPartNElts(i_part);
          double*      coords_cpl    = mesh_cpl -> getVertexCoords(i_part);
          CWP_g_num_t* gnum_vtx_cpl  = mesh_cpl -> getVertexGNum(i_part);
          CWP_g_num_t* gnum_elt_cpl  = mesh_cpl -> GNumEltsGet(i_part);

          PDM_dist_cloud_surf_surf_mesh_part_set (*id_dist,
                                                  i_part,
                                                  n_elts_cpl,
                                                  connecIdx_cpl,
                                                  connec_cpl,
                                                  gnum_elt_cpl,
                                                  n_vtx_cpl,
                                                  coords_cpl,
                                                  gnum_vtx_cpl);

        }//loop on part
      }//end if
    } // End if location_method == CWP_LOCATION_DIST_CLOUD_SURF

    else if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE ||
             location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {

      *id_dist = PDM_mesh_location_create (PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_cplComm);

      PDM_mesh_location_method_t mesh_location_method;
      if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE) {
        mesh_location_method = PDM_MESH_LOCATION_OCTREE;
      }
      else if (location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {
        mesh_location_method = PDM_MESH_LOCATION_DBBTREE;
      }

      PDM_mesh_location_method_set (*id_dist, mesh_location_method);
      PDM_mesh_location_tolerance_set (*id_dist, _get_location_tolerance());


      PDM_mesh_location_n_part_cloud_set (*id_dist, 0, _nb_part);

      PDM_mesh_location_mesh_global_data_set (*id_dist,
                                              _nb_part_cpl);

      for (int i_part = 0; i_part < _nb_part; i_part++) {

        CWP_g_num_t *gnum_target   = _gnum_target[i_part];
        double      *coords_target = _coords_target[i_part];

        PDM_mesh_location_cloud_set (*id_dist,
                                     0,
                                     i_part,
                                     _n_target[i_part],
                                     coords_target,
                                     gnum_target);
      }

      if (_both_codes_are_local == 0) {
        for (int i_part = 0; i_part < _nb_part_cpl; i_part++) {
          int n_elt_null = 0;
          int         *face_edge_idx = (int*) malloc (sizeof(int) * (n_elt_null + 1));
          int         *face_edge     = (int*) malloc (sizeof(int) * n_elt_null);
          int         *edge_vtx_idx  = (int*) malloc (sizeof(int) * (n_elt_null + 1));
          int         *edge_vtx      = (int*) malloc (sizeof(int) * n_elt_null);

          double      *coords        = (double*) malloc (sizeof(double) * n_elt_null * 3);
          CWP_g_num_t *vtx_gnum      = (CWP_g_num_t*) malloc (sizeof(CWP_g_num_t) * n_elt_null);
          CWP_g_num_t *edge_gnum     = (CWP_g_num_t*) malloc (sizeof(CWP_g_num_t) * n_elt_null);
          CWP_g_num_t *face_gnum     = (CWP_g_num_t*) malloc (sizeof(CWP_g_num_t) * n_elt_null);

          face_edge_idx[0] = 0;
          edge_vtx_idx[0] = 0;

          PDM_mesh_location_part_set_2d (*id_dist,
                                         i_part,
                                         0,
                                         face_edge_idx,
                                         face_edge,
                                         face_gnum,
                                         0,
                                         edge_vtx_idx,
                                         edge_vtx,
                                         edge_gnum,
                                         0,
                                         coords,
                                         vtx_gnum);
        }
      }

      else {
        Mesh *mesh_cpl = _spatial_interp_cpl -> _mesh;

        /*int mesh_nodal_id = mesh_cpl->getPdmNodalIndex();

        PDM_mesh_location_shared_nodal_mesh_set (*id_dist,
        mesh_nodal_id);*/
        for (int i_part = 0; i_part < _nb_part_cpl; i_part++) {
          int n_vtx  = mesh_cpl -> getPartNVertex(i_part);
          int n_face = mesh_cpl -> getNFace(i_part);
          int n_edge = mesh_cpl -> getNEdge(i_part);

          CWP_g_num_t *vtx_gnum  = mesh_cpl -> getVertexGNum(i_part);
          CWP_g_num_t *face_gnum = mesh_cpl -> GNumEltsGet(i_part);
          CWP_g_num_t *edge_gnum = NULL;//unused

          double *coords        = mesh_cpl -> getVertexCoords(i_part);
          int    *face_edge_idx = mesh_cpl -> getFaceEdgeIndex(i_part);
          int    *face_edge     = mesh_cpl -> getFaceEdge(i_part);
          int    *edge_vtx_idx  = mesh_cpl -> getEdgeVtxIndex(i_part);
          int    *edge_vtx      = mesh_cpl -> getEdgeVtx(i_part);


          PDM_mesh_location_part_set_2d (*id_dist,
                                         i_part,
                                         n_face,
                                         face_edge_idx,
                                         face_edge,
                                         face_gnum,
                                         n_edge,
                                         edge_vtx_idx,
                                         edge_vtx,
                                         edge_gnum,
                                         n_vtx,
                                         coords,
                                         vtx_gnum);
        }
      }

    } // End if location_method == CWP_LOCATION_MESH_LOCATION_*

    else {
      PDM_error (__FILE__, __LINE__, 0, "Unknown location method.\n");
    }

  }







  void SpatialInterpLocation::localization_null_setting_send(int* id_dist) {

    CWP_Location_method_t location_method = _get_location_method();

    if (location_method == CWP_LOCATION_DIST_CLOUD_SURF) {
      /* Paradigm mesh localisation _distance creation */
      *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_cplComm );

      PDM_dist_cloud_surf_n_part_cloud_set(*id_dist,   0, _nb_part_cpl);

      PDM_dist_cloud_surf_surf_mesh_global_data_set (*id_dist          ,
                                                     _n_g_elt_over_part,
                                                     _n_g_vtx_over_part,
                                                     _nb_part          );

      // printf("ENULL send %I64d %I64d _nb_part %i _nb_part_cpl %i\n",_n_g_elt_cpl_over_part,_n_g_vtx_cpl_over_part,_nb_part,_nb_part_cpl);

      for(int i_part =0;i_part<_nb_part_cpl;i_part++) {
        int n_elt_null = 0;
        double*      coords_null    = (double*)malloc(3*sizeof(double)*n_elt_null);
        CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*n_elt_null);

        PDM_dist_cloud_surf_cloud_set (*id_dist,
                                       0,
                                       i_part,
                                       0,
                                       NULL,//coords_null,
                                       NULL//gnum_elt_null
                                       );
      }

      for(int i_part =0; i_part<_nb_part; i_part++) {
        int n_elt_null = 0;
        int*         connecIdx_null = (int*)malloc(sizeof(int)*(1+n_elt_null));
        int*         connec_null    = (int*)malloc(sizeof(int)*n_elt_null);

        double*      coords_null    = (double*)malloc(3*sizeof(double)*n_elt_null);
        CWP_g_num_t* gnum_vtx_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*n_elt_null);
        CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*n_elt_null);

       // connecIdx_null[0]=0;

        PDM_dist_cloud_surf_surf_mesh_part_set (*id_dist,
                                                i_part        ,
                                                0             ,
                                                connecIdx_null,
                                                connec_null   ,
                                                gnum_elt_null ,
                                                0             ,
                                                coords_null   ,
                                                gnum_vtx_null );

      }
    } // End if location_method == CWP_LOCATION_DIST_CLOUD_SURF

    else if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE ||
             location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {

      *id_dist = PDM_mesh_location_create (PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_cplComm);

      PDM_mesh_location_method_t mesh_location_method;
      if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE) {
        mesh_location_method = PDM_MESH_LOCATION_OCTREE;
      }
      else if (location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {
        mesh_location_method = PDM_MESH_LOCATION_DBBTREE;
      }

      PDM_mesh_location_method_set (*id_dist, mesh_location_method);
      PDM_mesh_location_tolerance_set (*id_dist, _get_location_tolerance());


      PDM_mesh_location_n_part_cloud_set (*id_dist, 0, _nb_part_cpl);

      PDM_mesh_location_mesh_global_data_set (*id_dist,
                                              _nb_part);

      for (int i_part = 0; i_part < _nb_part_cpl; i_part++) {
        PDM_mesh_location_cloud_set (*id_dist,
                                     0,
                                     i_part,
                                     0,
                                     NULL,
                                     NULL);
      }


      for (int i_part = 0; i_part < _nb_part; i_part++) {
        int n_elt_null = 0;
        int         *face_edge_idx = (int*) malloc (sizeof(int) * (n_elt_null + 1));
        int         *face_edge     = (int*) malloc (sizeof(int) * n_elt_null);
        int         *edge_vtx_idx  = (int*) malloc (sizeof(int) * (n_elt_null + 1));
        int         *edge_vtx      = (int*) malloc (sizeof(int) * n_elt_null);

        double      *coords        = (double*) malloc (sizeof(double) * n_elt_null * 3);
        CWP_g_num_t *vtx_gnum      = (CWP_g_num_t*) malloc (sizeof(CWP_g_num_t) * n_elt_null);
        CWP_g_num_t *edge_gnum     = (CWP_g_num_t*) malloc (sizeof(CWP_g_num_t) * n_elt_null);
        CWP_g_num_t *face_gnum     = (CWP_g_num_t*) malloc (sizeof(CWP_g_num_t) * n_elt_null);

        face_edge_idx[0] = 0;
        edge_vtx_idx[0] = 0;

        PDM_mesh_location_part_set_2d (*id_dist,
                                       i_part,
                                       0,
                                       face_edge_idx,
                                       face_edge,
                                       face_gnum,
                                       0,
                                       edge_vtx_idx,
                                       edge_vtx,
                                       edge_gnum,
                                       0,
                                       coords,
                                       vtx_gnum);
      }
    } // End if location_method == CWP_LOCATION_MESH_LOCATION_*

    else {
      PDM_error (__FILE__, __LINE__, 0, "Unknown location method.\n");
    }
  }







  void SpatialInterpLocation::localization_null_setting_recv(int* id_dist) {

    CWP_Location_method_t location_method = _get_location_method();

    if (location_method == CWP_LOCATION_DIST_CLOUD_SURF) {
      /* Paradigm mesh localisation _distance creation */
      *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_cplComm );

      PDM_dist_cloud_surf_n_part_cloud_set(*id_dist,   0, _nb_part);

      PDM_dist_cloud_surf_surf_mesh_global_data_set (*id_dist,
                                                     _n_g_elt_cpl_over_part,
                                                     _n_g_vtx_cpl_over_part,
                                                     _nb_part_cpl);

      //printf("ENULL %I64d %I64d _nb_part %i _nb_part_cpl %i\n",_n_g_elt_cpl_over_part,_n_g_vtx_cpl_over_part,_nb_part,_nb_part_cpl);

      for(int i_part =0;i_part<_nb_part;i_part++) {
        int n_elt_null = 1;
        double*      coords_null    = (double*)malloc(3*sizeof(double)*n_elt_null);
        CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*n_elt_null);

        PDM_dist_cloud_surf_cloud_set (*id_dist,
                                       0,
                                       i_part,
                                       0,
                                       NULL,//coords_null,
                                       NULL//gnum_elt_null
                                       );
      }

      for(int i_part =0; i_part<_nb_part_cpl; i_part++) {
        int n_elt_null = 0;
        int*         connecIdx_null = (int*)malloc(sizeof(int)*(1+n_elt_null));
        int*         connec_null    = (int*)malloc(sizeof(int)*n_elt_null);

        double*      coords_null    = (double*)malloc(3*sizeof(double)*n_elt_null);
        CWP_g_num_t* gnum_vtx_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*n_elt_null);
        CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*n_elt_null);

        connecIdx_null[0]=0;

        PDM_dist_cloud_surf_surf_mesh_part_set (*id_dist,
                                                i_part,
                                                0,
                                                connecIdx_null,
                                                connec_null,
                                                gnum_elt_null,
                                                0,
                                                coords_null,
                                                gnum_vtx_null);
      }
    } // End if location_method == CWP_LOCATION_DIST_CLOUD_SURF

    else if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE ||
             location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {

      *id_dist = PDM_mesh_location_create (PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_cplComm);

      PDM_mesh_location_method_t mesh_location_method;
      if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE) {
        mesh_location_method = PDM_MESH_LOCATION_OCTREE;
      }
      else if (location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {
        mesh_location_method = PDM_MESH_LOCATION_DBBTREE;
      }

      PDM_mesh_location_method_set (*id_dist, mesh_location_method);
      PDM_mesh_location_tolerance_set (*id_dist, _get_location_tolerance());


      PDM_mesh_location_n_part_cloud_set (*id_dist, 0, _nb_part);

      PDM_mesh_location_mesh_global_data_set (*id_dist,
                                              _nb_part_cpl);


      for (int i_part = 0; i_part < _nb_part; i_part++) {
        PDM_dist_cloud_surf_cloud_set (*id_dist,
                                       0,
                                       i_part,
                                       0,
                                       NULL,
                                       NULL);
      }

      for (int i_part = 0; i_part < _nb_part_cpl; i_part++) {
        int n_elt_null = 0;
        int         *face_edge_idx = (int*) malloc (sizeof(int) * (n_elt_null + 1));
        int         *face_edge     = (int*) malloc (sizeof(int) * n_elt_null);
        int         *edge_vtx_idx  = (int*) malloc (sizeof(int) * (n_elt_null + 1));
        int         *edge_vtx      = (int*) malloc (sizeof(int) * n_elt_null);

        double      *coords        = (double*) malloc (sizeof(double) * n_elt_null * 3);
        CWP_g_num_t *vtx_gnum      = (CWP_g_num_t*) malloc (sizeof(CWP_g_num_t) * n_elt_null);
        CWP_g_num_t *edge_gnum     = (CWP_g_num_t*) malloc (sizeof(CWP_g_num_t) * n_elt_null);
        CWP_g_num_t *face_gnum     = (CWP_g_num_t*) malloc (sizeof(CWP_g_num_t) * n_elt_null);

        face_edge_idx[0] = 0;
        edge_vtx_idx[0] = 0;

        PDM_mesh_location_part_set_2d (*id_dist,
                                       i_part,
                                       0,
                                       face_edge_idx,
                                       face_edge,
                                       face_gnum,
                                       0,
                                       edge_vtx_idx,
                                       edge_vtx,
                                       edge_gnum,
                                       0,
                                       coords,
                                       vtx_gnum);
      }
    } // End if location_method == CWP_LOCATION_MESH_LOCATION_*

    else {
      PDM_error (__FILE__, __LINE__, 0, "Unknown location method.\n");
    }
  }






  void SpatialInterpLocation::localization_surface_setting(int* id_dist) {

    CWP_Location_method_t location_method = _get_location_method();

    if (location_method == CWP_LOCATION_DIST_CLOUD_SURF) {
      /* Paradigm mesh localisation _distance creation */
      *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_cplComm );

      //printf("_nb_part_cpl = %d\n", _nb_part_cpl);
      PDM_dist_cloud_surf_n_part_cloud_set(*id_dist, 0, _nb_part_cpl);
      //printf("_n_g_elt_over_part %i _n_g_vtx_over_part %i\n",_n_g_elt_over_part,_n_g_vtx_over_part);
      PDM_dist_cloud_surf_surf_mesh_global_data_set (*id_dist,
                                                     _n_g_elt_over_part,
                                                     _n_g_vtx_over_part,
                                                     _nb_part);

      for(int i_part =0;i_part<_nb_part;i_part++) {
        int* connecIdx = _mesh -> connecIdxGet(i_part);
        int* connec = _mesh -> connecGet(i_part);

        int n_vtx  = _mesh -> getPartNVertex(i_part);
        int n_elts  = _mesh -> getPartNElts(i_part);
        double* coords = _mesh -> getVertexCoords(i_part);
        CWP_g_num_t* gnum_vtx = _mesh -> getVertexGNum(i_part);
        CWP_g_num_t* gnum_elt = _mesh -> GNumEltsGet(i_part);

        PDM_dist_cloud_surf_surf_mesh_part_set (*id_dist,
                                                i_part,
                                                n_elts,
                                                connecIdx,
                                                connec,
                                                gnum_elt,
                                                n_vtx,
                                                coords,
                                                gnum_vtx);
      }

      if(_both_codes_are_local == 0) {
        for(int i_part =0; i_part<_nb_part_cpl; i_part++) {

          int n_elt_null = 0;
          double*      coords    = (double*)malloc(3*sizeof(double)*n_elt_null);
          CWP_g_num_t* gnum_elt  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*n_elt_null);

          PDM_dist_cloud_surf_cloud_set (*id_dist,
                                         0,
                                         i_part,
                                         n_elt_null,
                                         coords ,
                                         gnum_elt
                                         );
        }//loop on part
      }
      else {
        for(int i_part =0; i_part<_nb_part_cpl; i_part++) {
          int          n_target_cpl      = _spatial_interp_cpl -> _n_target     [i_part];
          CWP_g_num_t* gnum_target_cpl   = _spatial_interp_cpl -> _gnum_target  [i_part];
          double*      coords_target_cpl = _spatial_interp_cpl -> _coords_target[i_part];
          PDM_dist_cloud_surf_cloud_set (*id_dist,
                                         0,
                                         i_part,
                                         n_target_cpl,
                                         coords_target_cpl,
                                         gnum_target_cpl
                                         );
        }//loop on part
      } //end of if

    } // End if location_method == CWP_LOCATION_DIST_CLOUD_SURF

    else if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE ||
             location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {

      *id_dist = PDM_mesh_location_create (PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_cplComm);

      PDM_mesh_location_method_t mesh_location_method;
      if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE) {
        mesh_location_method = PDM_MESH_LOCATION_OCTREE;
      }
      else if (location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {
        mesh_location_method = PDM_MESH_LOCATION_DBBTREE;
      }

      PDM_mesh_location_method_set (*id_dist, mesh_location_method);
      PDM_mesh_location_tolerance_set (*id_dist, _get_location_tolerance());


      PDM_mesh_location_n_part_cloud_set (*id_dist, 0, _nb_part_cpl);

      PDM_mesh_location_mesh_global_data_set (*id_dist,
                                              _nb_part);
      /*int mesh_nodal_id = _mesh->getPdmNodalIndex();

      PDM_mesh_location_shared_nodal_mesh_set (*id_dist,
      mesh_nodal_id);*/
      for (int i_part = 0; i_part < _nb_part; i_part++) {
          int n_vtx  = _mesh -> getPartNVertex(i_part);
          int n_face = _mesh -> getNFace(i_part);
          int n_edge = _mesh -> getNEdge(i_part);

          CWP_g_num_t *vtx_gnum  = _mesh -> getVertexGNum(i_part);
          CWP_g_num_t *face_gnum = _mesh -> GNumEltsGet(i_part);
          CWP_g_num_t *edge_gnum = NULL;//unused

          double *coords        = _mesh -> getVertexCoords(i_part);
          int    *face_edge_idx = _mesh -> getFaceEdgeIndex(i_part);
          int    *face_edge     = _mesh -> getFaceEdge(i_part);
          int    *edge_vtx_idx  = _mesh -> getEdgeVtxIndex(i_part);
          int    *edge_vtx      = _mesh -> getEdgeVtx(i_part);

          PDM_mesh_location_part_set_2d (*id_dist,
                                         i_part,
                                         n_face,
                                         face_edge_idx,
                                         face_edge,
                                         face_gnum,
                                         n_edge,
                                         edge_vtx_idx,
                                         edge_vtx,
                                         edge_gnum,
                                         n_vtx,
                                         coords,
                                         vtx_gnum);
        }

      if (_both_codes_are_local == 0) {
        for (int i_part = 0; i_part < _nb_part_cpl; i_part++) {

          int n_elt_null = 0;
          double      *coords   = (double*)      malloc (sizeof(double)      * n_elt_null * 3);
          CWP_g_num_t *gnum_elt = (CWP_g_num_t*) malloc (sizeof(CWP_g_num_t) * n_elt_null);

          PDM_mesh_location_cloud_set (*id_dist,
                                       0,
                                       i_part,
                                       n_elt_null,
                                       coords ,
                                       gnum_elt);
        }
      }
      else {
        for (int i_part = 0; i_part < _nb_part_cpl; i_part++) {
          int          n_target_cpl      = _spatial_interp_cpl -> _n_target     [i_part];
          CWP_g_num_t *gnum_target_cpl   = _spatial_interp_cpl -> _gnum_target  [i_part];
          double      *coords_target_cpl = _spatial_interp_cpl -> _coords_target[i_part];
          PDM_mesh_location_cloud_set (*id_dist,
                                       0,
                                       i_part,
                                       n_target_cpl,
                                       coords_target_cpl,
                                       gnum_target_cpl);
        }
      }

    } // End if location_method == CWP_LOCATION_MESH_LOCATION_*

    else {
      PDM_error (__FILE__, __LINE__, 0, "Unknown location method.\n");
    }
  }







  void SpatialInterpLocation::localization_compute(int id_dist) {

    CWP_Location_method_t location_method = _get_location_method();

    if (location_method == CWP_LOCATION_DIST_CLOUD_SURF) {
      PDM_dist_cloud_surf_compute (id_dist);
      PDM_dist_cloud_surf_dump_times (id_dist);
    }

    else if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE ||
             location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {
      PDM_mesh_location_compute (id_dist);
      PDM_mesh_location_dump_times (id_dist);
    }

    else {
      PDM_error (__FILE__, __LINE__, 0, "Unknown location method.\n");
    }
  }






  void SpatialInterpLocation::localization_get_cpl(int id_dist) {

    CWP_Location_method_t location_method = _get_location_method();

    if (location_method == CWP_LOCATION_DIST_CLOUD_SURF) {
      _spatial_interp_cpl -> _distance           = (double**)malloc(sizeof(double*) * _nb_part_cpl);
      _spatial_interp_cpl -> _projected          = (double**)malloc(sizeof(double*) * _nb_part_cpl);
      _spatial_interp_cpl -> _closest_elt_gnum   = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*) * _nb_part_cpl);

      for(int i_part =0;i_part<_nb_part_cpl;i_part++) {
        int          n_target_cpl    = _spatial_interp_cpl -> _n_target[i_part];
        PDM_dist_cloud_surf_get (id_dist,
                                 0,
                                 i_part,
                                 &(_spatial_interp_cpl -> _distance[i_part]),
                                 &(_spatial_interp_cpl -> _projected[i_part]),
                                 &(_spatial_interp_cpl -> _closest_elt_gnum[i_part]));

        for(int i=0;i<n_target_cpl;i++){
          if(_spatial_interp_cpl -> _closest_elt_gnum[i_part][i]>CWP_g_num_t(_n_g_elt_over_part)
             || _spatial_interp_cpl -> _closest_elt_gnum[i_part][i]<CWP_g_num_t(1)
             || _spatial_interp_cpl -> _distance[i_part][i]>0.01){
            _spatial_interp_cpl -> _closest_elt_gnum[i_part][i]=CWP_g_num_t(1);
            _spatial_interp_cpl -> _distance[i_part][i]=INFINITY;
          }
        }
      } //end loop on i_part

    } // End if location_method == CWP_LOCATION_DIST_CLOUD_SURF

    else if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE ||
             location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {
      _spatial_interp_cpl -> _distance          = (double **)      malloc (sizeof(double *)      * _nb_part_cpl);
      _spatial_interp_cpl -> _projected         = (double **)      malloc (sizeof(double *)      * _nb_part_cpl);
      _spatial_interp_cpl -> _closest_elt_gnum  = (CWP_g_num_t **) malloc (sizeof(CWP_g_num_t *) * _nb_part_cpl);

      int          useless_n_points;
      double      *useless_coord;
      PDM_g_num_t *useless_g_num;

      int    *weights_idx; //?
      double *weights;     //?

      for (int i_part = 0; i_part < _nb_part_cpl; i_part++) {
        int n_target_cpl = _spatial_interp_cpl -> _n_target[i_part];

        PDM_mesh_location_get (id_dist,
                               0,
                               i_part,
                               &useless_n_points,
                               &useless_coord,
                               &useless_g_num,
                               &(_spatial_interp_cpl -> _closest_elt_gnum[i_part]),
                               &weights_idx,
                               &weights,
                               &(_spatial_interp_cpl -> _projected[i_part]));

        _spatial_interp_cpl -> _distance[i_part] = (double *) malloc (sizeof(double) * n_target_cpl);
        //double *coords_target = _spatial_interp_cpl -> _coords_target[i_part];

        for (int i = 0; i < n_target_cpl; i++) {
          _spatial_interp_cpl -> _distance[i_part][i] = 0.;
          /*for (int j = 0; j < 3; j++) {
            double d = _projected[i_part][3*i + j] - coords_target[3*i + j];
            _spatial_interp_cpl -> _distance[i_part][i] += d*d;
            }*/

          if (_spatial_interp_cpl -> _closest_elt_gnum[i_part][i] > CWP_g_num_t(_n_g_elt_cpl_over_part) ||
              _spatial_interp_cpl -> _closest_elt_gnum[i_part][i] < CWP_g_num_t(1) ||
              _spatial_interp_cpl -> _distance [i_part][i] > 0.1) {
            _spatial_interp_cpl -> _closest_elt_gnum[i_part][i] = CWP_g_num_t(1);
            _spatial_interp_cpl -> _distance[i_part][i] = INFINITY;
          }
        }
      }

    } // End if location_method == CWP_LOCATION_MESH_LOCATION_*

    else {
      PDM_error (__FILE__, __LINE__, 0, "Unknown location method.\n");
    }
  }







  void SpatialInterpLocation::localization_get(int id_dist) {

    CWP_Location_method_t location_method = _get_location_method();

    if (location_method == CWP_LOCATION_DIST_CLOUD_SURF) {
      _distance           = (double**)malloc(sizeof(double*) * _nb_part);
      _projected          = (double**)malloc(sizeof(double*) * _nb_part);
      _closest_elt_gnum   = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*) * _nb_part);

      for (int i_part = 0; i_part < _nb_part; i_part++) {

        PDM_dist_cloud_surf_get (id_dist,
                                 0,
                                 i_part,
                                 &(_distance [i_part]),
                                 &(_projected[i_part]),
                                 &(_closest_elt_gnum[i_part]));


        for (int i = 0; i < _n_target[i_part]; i++){
          if (_closest_elt_gnum[i_part][i] > CWP_g_num_t(_n_g_elt_cpl_over_part) ||
              _closest_elt_gnum[i_part][i] < CWP_g_num_t(1) ||
              _distance [i_part][i] > 0.1) {
            _closest_elt_gnum[i_part][i] = CWP_g_num_t(1);
            _distance[i_part][i] = INFINITY;
          }
        }
      }

    } // End if location_method == CWP_LOCATION_DIST_CLOUD_SURF

    else if (location_method == CWP_LOCATION_MESH_LOCATION_OCTREE ||
             location_method == CWP_LOCATION_MESH_LOCATION_DBBTREE) {
      _distance          = (double **)      malloc (sizeof(double *)      * _nb_part);
      _projected         = (double **)      malloc (sizeof(double *)      * _nb_part);
      _closest_elt_gnum  = (CWP_g_num_t **) malloc (sizeof(CWP_g_num_t *) * _nb_part);

      int          useless_n_points;
      double      *useless_coord;
      PDM_g_num_t *useless_g_num;

      int    *weights_idx; //?
      double *weights;     //?

      for (int i_part = 0; i_part < _nb_part; i_part++) {

        PDM_mesh_location_get (id_dist,
                               0,
                               i_part,
                               &useless_n_points,
                               &useless_coord,
                               &useless_g_num,
                               &(_closest_elt_gnum[i_part]),
                               &weights_idx,
                               &weights,
                               &(_projected[i_part]));

        _distance[i_part] = (double *) malloc (sizeof(double) * _n_target[i_part]);
        double *coords_target = _coords_target[i_part];

        for (int i = 0; i < _n_target[i_part]; i++) {
          _distance[i_part][i] = 0.;
          for (int j = 0; j < 3; j++) {
            double d = _projected[i_part][3*i + j] - coords_target[3*i + j];
            _distance[i_part][i] += d*d;
          }

          if (_closest_elt_gnum[i_part][i] > CWP_g_num_t(_n_g_elt_cpl_over_part) ||
              _closest_elt_gnum[i_part][i] < CWP_g_num_t(1) ||
              _distance [i_part][i] > 0.1) {
            _closest_elt_gnum[i_part][i] = CWP_g_num_t(1);
            _distance[i_part][i] = INFINITY;
          }
        }
      }

    } // End if location_method == CWP_LOCATION_MESH_LOCATION_*

    else {
      PDM_error (__FILE__, __LINE__, 0, "Unknown location method.\n");
    }
 }// End locate_cell_point





      /***********************************************************
       ***********************************************************
       **                                                       **
       **   Process, partition, num triplet location from       **
       **           global numbering functions                  **
       **                                                       **
       ***********************************************************
       ***********************************************************/




 void SpatialInterpLocation::triplet_location_request(int* id_gnum_location) {

  *id_gnum_location = PDM_gnum_location_create(_nb_part_cpl,_nb_part, _pdm_cplComm);

  for(int i_part =0;i_part<_nb_part;i_part++) {
    for(int i=0; i<_n_target[i_part]; i++) {
 //    if(_distance[i_part][i] == INFINITY ) {
  /*    printf("_closest_elt_gnum[%i][%i] rank %i %I64d coords %f %f %f _distance %f N %i\n",
      i_part,i,_rank,_closest_elt_gnum[i_part][i],
      _projected[i_part][3*i],_projected[i_part][3*i+1],_projected[i_part][3*i+2],
      _distance[i_part][i],
      _n_target[i_part]);
  */
    // }
    }
    PDM_gnum_location_requested_elements_set(*id_gnum_location,i_part, _n_target[i_part],&(_closest_elt_gnum[i_part][0]));
  }

  if(_both_codes_are_local == 0) {
    for(int i_part =0; i_part<_nb_part_cpl; i_part++) {
      CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*0);
      PDM_gnum_location_elements_set(*id_gnum_location,i_part,0, gnum_elt_null);
    }
  }
  else {
    for(int i_part =0; i_part<_nb_part_cpl; i_part++) {
      Mesh* mesh_cpl = _spatial_interp_cpl -> _mesh;
      CWP_g_num_t* gnum_elt_cpl = mesh_cpl -> GNumEltsGet(i_part);
      int          n_elt_cpl    = mesh_cpl -> getPartNElts(i_part);

      PDM_gnum_location_elements_set(*id_gnum_location,i_part, n_elt_cpl,gnum_elt_cpl);
    }//loop on part
  }//end if
 }






 void SpatialInterpLocation::triplet_location_set(int* id_gnum_location) {

  *id_gnum_location = PDM_gnum_location_create(_nb_part,_nb_part_cpl, _pdm_cplComm);

  for(int i_part =0;i_part<_nb_part;i_part++) {

    CWP_g_num_t* gnum_elt = _mesh -> GNumEltsGet(i_part);
    /*printf("rank %i _n_elt[%i]  _nb_part %i _nb_part_cpl %i _both_codes_are_local %i %i\n",
    _rank,i_part,_nb_part,_nb_part_cpl,_both_codes_are_local,_n_g_vtx_cpl_over_part);
    printf("rank %i _n_elt[%i] %i _nb_part %i _nb_part_cpl %i\n",_rank,i_part,_n_elt[i_part],_nb_part,_nb_part_cpl);
    */
    PDM_gnum_location_elements_set(*id_gnum_location,i_part, _n_elt[i_part],gnum_elt);
  }

  if(_both_codes_are_local == 0) {
    for(int i_part =0; i_part<_nb_part_cpl; i_part++) {
      CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*0);
      PDM_gnum_location_requested_elements_set(*id_gnum_location,i_part, 0,gnum_elt_null);
    }
  }
  else {
    for(int i_part =0; i_part<_nb_part_cpl; i_part++) {
      int n_target_cpl = _spatial_interp_cpl -> _n_target[i_part];
      PDM_gnum_location_requested_elements_set(*id_gnum_location,i_part, n_target_cpl, _spatial_interp_cpl -> _closest_elt_gnum[i_part]);
   }//loop on part
  }//end if
 }





void SpatialInterpLocation::triplet_location_null_send(int* id_gnum_location) {

  *id_gnum_location = PDM_gnum_location_create(_nb_part,_nb_part_cpl, _pdm_cplComm);

  for(int i_part =0;i_part<_nb_part;i_part++) {
    CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*0);
    PDM_gnum_location_elements_set(*id_gnum_location,i_part, 0,gnum_elt_null);
  }

  for(int i_part =0; i_part<_nb_part_cpl; i_part++) {
    CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*0);
    PDM_gnum_location_requested_elements_set(*id_gnum_location,i_part, 0,gnum_elt_null);
  }//loop on part
 }





void SpatialInterpLocation::triplet_location_null_recv(int* id_gnum_location) {

  *id_gnum_location = PDM_gnum_location_create(_nb_part_cpl,_nb_part, _pdm_cplComm);

  for(int i_part =0;i_part<_nb_part_cpl;i_part++) {
    CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*0);
    PDM_gnum_location_elements_set(*id_gnum_location,i_part, 0,gnum_elt_null);
  }

  for(int i_part =0; i_part<_nb_part; i_part++) {
    CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*0);
    PDM_gnum_location_requested_elements_set(*id_gnum_location,i_part, 0,gnum_elt_null);
  }//loop on part
 }






 void SpatialInterpLocation::triplet_location_compute(int id_gnum_location) {
    PDM_gnum_location_compute(id_gnum_location);
 }





 void SpatialInterpLocation::triplet_location_get(int id_gnum_location) {

  _target_proc_part_num_idx =(int**)malloc(sizeof(int*)*_nb_part);
  _target_proc_part_num     =(int**)malloc(sizeof(int*)*_nb_part);

  for (int i_part = 0; i_part < _nb_part; i_part++) {

    PDM_gnum_location_get(id_gnum_location,
                            i_part,
                            &(_target_proc_part_num_idx[i_part]),
                            &(_target_proc_part_num    [i_part])
                            );
   }
 }






  void SpatialInterpLocation::triplet_location_get_cpl(int id_gnum_location) {
    _spatial_interp_cpl ->_target_proc_part_num_idx =(int**)malloc(sizeof(int*)*_nb_part_cpl);
    _spatial_interp_cpl ->_target_proc_part_num     =(int**)malloc(sizeof(int*)*_nb_part_cpl);

    for(int i_part =0;i_part<_nb_part_cpl;i_part++) {
      PDM_gnum_location_get(id_gnum_location,
                            i_part,
                            &(_spatial_interp_cpl ->_target_proc_part_num_idx[i_part]),
                            &(_spatial_interp_cpl ->_target_proc_part_num[i_part])
                            );
    }
  }// End locate_cell_point





      /***********************************************************
       ***********************************************************
       **            Communication tree array functions         **
       **                                                       **
       ***********************************************************
       ***********************************************************/




 void SpatialInterpLocation:: initialization_of_receving_communication_tree_array() {
   if(_targets_localization_idx_cpl == NULL){
     _targets_localization_idx_cpl = (int**)malloc(  sizeof(int*) *_n_ranks_g);
     for(int i =0; i<_n_ranks_g;i++)
       _targets_localization_idx_cpl[i] = NULL;
   }

   for(int i =0; i<_n_ranks_g;i++) {
     if(_targets_localization_idx_cpl[i] == NULL)
       _targets_localization_idx_cpl[i] = (int*)malloc(sizeof(int)*(1+_nb_part));
     for (int i_part = 0; i_part < _nb_part+1; i_part++) {
       _targets_localization_idx_cpl[i][i_part] = 0;
     }
   }

 }






 void SpatialInterpLocation::filling_of_sending_communication_tree_array() {

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

   for (int i_part = 0; i_part < _nb_part; i_part++) {
     for(int k=0;k<_n_target[i_part];k++){
     //  printf("_target_proc_part_num_idx[i_part][%i] rank %i %i _n_target[i_part] %i\n",k,_rank,_target_proc_part_num_idx[i_part][k],_n_target[i_part]);
     //  printf("_target_proc_part_num[i_part][%i] rank %i %i\n",k,_rank,_target_proc_part_num[i_part][ _target_proc_part_num_idx[i_part][k] ]);
       int elt_proc = _target_proc_part_num[i_part][ _target_proc_part_num_idx[i_part][k] ];
       int elt_part = _target_proc_part_num[i_part][ _target_proc_part_num_idx[i_part][k] + 1];

       _targets_localization_idx    [elt_proc][elt_part]++;
       _process_and_partition_count [elt_proc][elt_part]++;
    }
  }//end i_part

 /*  for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
     for (int i_part = 0; i_part < _nb_part_cpl+1; i_part++) {
       printf("VERIF _targets_localization_idx[%i][%i] _rank %i %i _Texch_t %i\n",i_proc,i_part,_rank,_targets_localization_idx [i_proc][i_part],_Texch_t);
     }
   }

*/
  _transform_to_index(_targets_localization_idx,_n_ranks_g,_nb_part_cpl);

/*
   for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
     for (int i_part = 0; i_part < _nb_part_cpl+1; i_part++) {
       printf("VERIFINDEX _targets_localization_idx[%i][%i] _rank %i %i _Texch_t %i\n",i_proc,i_part,_rank,_targets_localization_idx [i_proc][i_part],_Texch_t);
     }
   }
*/
  int** idx_proc = (int**)malloc(sizeof(int*)* _n_ranks_g);
  for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
   idx_proc[i_proc] = (int*)malloc(sizeof(int)* _nb_part_cpl);
   for (int i_part = 0; i_part < _nb_part_cpl; i_part++)
      idx_proc[i_proc][i_part]=0;
   }


  if(_targets_localization_data!=NULL) free(_targets_localization_data);
  _targets_localization_data = (target_data*)malloc(sizeof(target_data)*_targets_localization_idx[_n_ranks_g-1][_nb_part_cpl]);

  for (int i_part = 0; i_part < _nb_part; i_part++) {
    for(int k=0;k<_n_target[i_part];k++){
       int elt_proc = _target_proc_part_num[i_part][ _target_proc_part_num_idx[i_part][k]    ];
       int elt_part = _target_proc_part_num[i_part][ _target_proc_part_num_idx[i_part][k] + 1];
       int num = _target_proc_part_num[i_part][ _target_proc_part_num_idx[i_part][k] + 2] - 1;

       int idx = _targets_localization_idx[elt_proc][elt_part] + idx_proc[elt_proc][elt_part];
       //Local Numbering
       _targets_localization_data [idx].lnum              = num;
       //Coupled numbering
       _targets_localization_data [idx ].origin_part      = i_part;
       _targets_localization_data [idx ].projectedX       = _projected [i_part][ 3*k    ];
       _targets_localization_data [idx ].projectedY       = _projected [i_part][ 3*k +1 ];
       _targets_localization_data [idx ].projectedZ       = _projected [i_part][ 3*k +2 ];
       _targets_localization_data [idx ].distance         = _distance  [i_part][   k    ];
       //Closest local element numbering
       _targets_localization_data [idx ].closest_elt_gnum  = _closest_elt_gnum[i_part][ k ];
       //Coupled process origin
       _targets_localization_data [idx].origin_proc       = _rank ;
       //Coupled origin partition
       _targets_localization_data [idx ].closest_elt_part  = elt_part;
       _targets_localization_data [idx ].l_num_origin      = k;
       idx_proc[elt_proc][elt_part]++;

 /*     printf("_targets_localization_data rank %i [%i] lnum %i origin_part %i closest_elt_gnum %i elt_proc %i elt_part %i distance %f proj %f\n",
       _rank,k,
       _targets_localization_data [idx ].lnum,
       _targets_localization_data [idx ].origin_part,
       int(_targets_localization_data [idx ].closest_elt_gnum),
       elt_proc,elt_part,
       _targets_localization_data [idx ].distance,
       _targets_localization_data [idx ].projectedX);*/
      }
    }//end i_part

    for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
      free(idx_proc[i_proc]);
    }

    for (int i_part = 0; i_part < _nb_part; i_part++) {
      free(_distance [i_part]);
      free(_projected[i_part]);
      free(_closest_elt_gnum[i_part]);
      free(_target_proc_part_num_idx[i_part]);
      free(_target_proc_part_num[i_part]);
    }

    free(_target_proc_part_num_idx      );
    free(_target_proc_part_num          );
    free(idx_proc);
    free(_distance );
    free(_projected);
    free(_closest_elt_gnum);

 }





      /***********************************************************
       ***********************************************************
       **            Data index communication functions         **
       **                                                       **
       ***********************************************************
       ***********************************************************/




 void SpatialInterpLocation::data_index_communication_send() {

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






 void SpatialInterpLocation::data_index_communication_send_p2p() {

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





 void SpatialInterpLocation::data_index_communication_recv() {

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




 void SpatialInterpLocation::data_index_communication_recv_p2p() {

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




 void SpatialInterpLocation::both_index_communication_p2p() {

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



 void SpatialInterpLocation::data_index_communication_null() {

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





 void SpatialInterpLocation::both_index_communication() {

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








      /***********************************************************
       ***********************************************************
       **            Data communication functions               **
       **                                                       **
       ***********************************************************
       ***********************************************************/





 void SpatialInterpLocation::prepare_data_communication_send() {

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





 void SpatialInterpLocation::prepare_data_communication_recv() {

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





 void SpatialInterpLocation::data_communication_send() {

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



 void SpatialInterpLocation::data_communication_send_p2p() {

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




 void SpatialInterpLocation::data_communication_recv() {
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



 void SpatialInterpLocation::data_communication_recv_p2p() {

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






 void SpatialInterpLocation::both_data_communication() {

   MPI_Request req;
   MPI_Ialltoallv((void*) (_spatial_interp_cpl -> _targets_localization_data), _spatial_interp_cpl -> _targets_localization_data_count_send, _spatial_interp_cpl -> _targets_localization_data_disp_send, MPI_BYTE,
                  (void*)_targets_localization_data_cpl,  _targets_localization_data_count_recv, _targets_localization_data_disp_recv, MPI_BYTE,
                 _cplComm,&req);

   MPI_Status stat;
   MPI_Wait(&req,&stat);

  }


 void SpatialInterpLocation::both_data_communication_p2p() {

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






 void SpatialInterpLocation::data_communication_null() {

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




 void SpatialInterpLocation::data_communication_wait_send() {
   free(_targets_localization_data_count_send);
   free(_targets_localization_data_disp_send );
 }





 void SpatialInterpLocation::data_communication_wait_recv() {
   free(_targets_localization_data_count_recv);
   free(_targets_localization_data_disp_recv );
   _n_tot_target_cpl    = _targets_localization_idx_cpl[_n_ranks_g-1][_nb_part];
 }







  void SpatialInterpLocation::user_targets_gnum_compute() {
    int coord_def = 1;
    for (int i=0; i<_nb_part; i++){
      if(_coords_user_targets[i] == NULL) {
        coord_def = 0;
        break;
      }
    }

    if(coord_def == 1) {
      _pdmGNum_handle_index  = PDM_gnum_create (3, _nb_part, PDM_FALSE, 1e-3, _pdm_connectableComm);
      for (int i_part=0; i_part<_nb_part; i_part++){
        PDM_gnum_set_from_coords (_pdmGNum_handle_index, i_part, _n_user_targets[i_part], _coords_user_targets[i_part], NULL);
      }
      PDM_gnum_compute (_pdmGNum_handle_index);
      for (int i_part=0; i_part<_nb_part; i_part++){
        _gnum_user_targets[i_part] = const_cast<CWP_g_num_t*>(PDM_gnum_get(_pdmGNum_handle_index, i_part));
      }
    }
  }





  void SpatialInterpLocation::computeFree(){

   for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
     if(_targets_localization_idx_cpl != NULL) free(_targets_localization_idx_cpl[i_proc]);
     if(_targets_localization_idx     != NULL) free(_targets_localization_idx[i_proc]);
   }

   free(_targets_localization_data_cpl);
   if(_targets_localization_idx_cpl != NULL) free(_targets_localization_idx_cpl);
   if(_Texch_t == CWP_FIELD_EXCH_RECV) free(_targets_localization_data);
   if(_targets_localization_idx     != NULL) free(_targets_localization_idx);
  }





  void* SpatialInterpLocation::interpolate(Field* referenceField) {

    int                 nComponent         = referenceField -> nComponentGet  ();
    CWP_Dof_location_t   referenceFieldType = referenceField -> typeGet        ();
    int                 dataTypeSize       = referenceField -> dataTypeSizeGet();
    CWP_Interpolation_t interpolationType  = referenceField -> interpolationTypeGet();
    void               *interpolatedData   = referenceField -> sendBufferGet  ();

    if (interpolatedData != NULL) free(interpolatedData);
    interpolatedData = (void*) malloc( dataTypeSize * nComponent*_n_tot_target_cpl);

    if (_weights_src_idx == NULL) {
      _weights_src_idx  = (int **) malloc(sizeof(int*)* _nb_part);
      for (int i = 0; i < _nb_part; i++) {
        _weights_src_idx[i] = NULL;
      }
      _weights_src  = (double **) malloc(sizeof(double*)* _nb_part);
      for (int i = 0; i < _nb_part; i++) {
        _weights_src[i] = NULL;
      }
    }


    for(int i_part=0;i_part<_nb_part;i_part++){

      int weights_src_empty = 0;
      int l_weights_src = -1;
      int s_weights_src = -1;

      if (_weights_src_idx[i_part] == NULL) {
        weights_src_empty = 1;
        _weights_src_idx[i_part] = (int *) malloc (sizeof(int) *
                 (_targets_localization_idx_cpl[_n_ranks_g-1][i_part+1] + 1));
        _weights_src_idx[i_part][0] = 0;

        s_weights_src =
          4 * _targets_localization_idx_cpl[_n_ranks_g-1][i_part+1];
        _weights_src[i_part] =
          (double *) malloc (sizeof(double) * s_weights_src);
        l_weights_src = 0;
      }

      void* referenceData = referenceField -> dataGet(i_part);
      // For a cell center field : give the value of the located cell

      for(int i_proc=0;i_proc<_n_ranks_g;i_proc++){
        if(interpolationType == CWP_INTERPOLATION_USER) {

            CWP_Interp_from_location_t interpolationFunction  = referenceField -> interpolationFunctionGet();

            int n_tgt = _targets_localization_idx_cpl[i_proc][i_part+1] - _targets_localization_idx_cpl[i_proc][i_part];

            int* connecIdx = _mesh -> connecIdxGet(i_part);
            int* connec    = _mesh -> connecGet(i_part);
            double* coords = _mesh -> getVertexCoords(i_part);

            int    *tgt_pts_location = (int*)    malloc(sizeof(int)   *n_tgt);
            int    *tgt_pts_location_p1 = (int*)    malloc(sizeof(int)   *n_tgt);
            double *tgt_pts_dist     = (double*) malloc(sizeof(double)*n_tgt);
            double *tgt_pts_projected_coords = (double*) malloc(3 * sizeof(double)*n_tgt);
            int    *tgt_pts_bary_coords_idx = NULL;
            double *tgt_pts_bary_coords     = NULL;

            int i=0;
            for (int itarget = _targets_localization_idx_cpl[i_proc][i_part]; itarget < _targets_localization_idx_cpl[i_proc][i_part+1]; itarget++) {
              tgt_pts_location[i] = _targets_localization_data_cpl[itarget].lnum     ;
              tgt_pts_location_p1[i] = _targets_localization_data_cpl[itarget].lnum + 1 ;
              tgt_pts_dist    [i] = _targets_localization_data_cpl[itarget].distance ;
              double x_target = _targets_localization_data_cpl[itarget].projectedX;
              double y_target = _targets_localization_data_cpl[itarget].projectedY;
              double z_target = _targets_localization_data_cpl[itarget].projectedZ;

              tgt_pts_projected_coords[3*i   ] = x_target;
              tgt_pts_projected_coords[3*i +1] = y_target;
              tgt_pts_projected_coords[3*i +2] = z_target;
              i++;
            }

            PDM_geom_elem_compute_polygon_barycentric_coordinates(n_tgt,
                                                                  tgt_pts_location_p1,
                                                                  tgt_pts_projected_coords,
                                                                  connecIdx,
                                                                  connec,
                                                                  coords,
                                                                  &tgt_pts_bary_coords_idx,
                                                                  &tgt_pts_bary_coords
                                                                  );

            void* tmpData = (char*) interpolatedData + dataTypeSize * nComponent * _targets_localization_idx_cpl[i_proc][i_part];

            (*interpolationFunction)( CWP_INTERFACE_SURFACE   ,
                                   _n_vtx[i_part]          ,
                                   _n_elt[i_part]          ,
                                   n_tgt                   ,
                                   coords                  ,
                                   connecIdx               ,
                                   connec                  ,
                                   tgt_pts_projected_coords,
                                   tgt_pts_location        ,
                                   tgt_pts_dist            ,
                                   tgt_pts_bary_coords_idx ,
                                   tgt_pts_bary_coords     ,
                                   nComponent              ,
                                   referenceFieldType      ,
                                   referenceData           ,
                                   referenceFieldType      ,
                                   tmpData
                                 );

             if(tgt_pts_location    != NULL ) free(tgt_pts_location             );
             if(tgt_pts_location_p1 != NULL ) free(tgt_pts_location_p1          );
             if(tgt_pts_dist        != NULL ) free(tgt_pts_dist                 );
             if(tgt_pts_projected_coords != NULL ) free(tgt_pts_projected_coords);
             if(tgt_pts_bary_coords_idx  != NULL ) free(tgt_pts_bary_coords_idx );
             if(tgt_pts_bary_coords      != NULL ) free(tgt_pts_bary_coords     );


        }
        else{
         if (referenceFieldType == CWP_DOF_LOCATION_CELL_CENTER) {

          for (int itarget = _targets_localization_idx_cpl[i_proc][i_part]; itarget < _targets_localization_idx_cpl[i_proc][i_part+1]; itarget++) {
            //Index of the corresponding local reference Data.
            int iel = _targets_localization_data_cpl[itarget].lnum ;
            // Index in the interpolated Data array
            int interpInd = itarget;
          //  printf("iel %i itarget %i refData %f _n_tot_target_cpl %i\n",iel,itarget,referenceData[iel],_n_tot_target_cpl);
            for (int k = 0; k < nComponent; k++) {
              memcpy( (char*)interpolatedData + dataTypeSize * ( nComponent*interpInd + k ) ,
                      (char*)referenceData + dataTypeSize * ( nComponent*iel + k )          ,
                      dataTypeSize);
              /*printf("interpolatedData[ nComponent*interpInd + k  ] %f i_part %i i_proc %i nComponent*interpInd + k %i
                        nComponent*iel + k %i referenceData[nComponent*iel + k ] %f\n",
                        interpolatedData[ nComponent*interpInd + k  ],i_part,i_proc,
                        nComponent*interpInd + k,nComponent*iel + k),referenceData[nComponent*iel + k ]);
               */
            }
          } // loop on itarget

        } // if referenceFieldType == CWP_DOF_LOCATION_CELL_CENTER
        else if (referenceFieldType == CWP_DOF_LOCATION_NODE) {

          int* connecIdx = _mesh -> connecIdxGet(i_part);
          int* connec    = _mesh -> connecGet(i_part);

          double* coords = _mesh -> getVertexCoords(i_part);
          //printf("_targets_localization_idx_cpl[%i][%i] rank %i %i %i\n",
         // i_proc,i_part,_rank,_targets_localization_idx_cpl[i_proc][i_part],_targets_localization_idx_cpl[i_proc][i_part+1]);

          for (int itarget = _targets_localization_idx_cpl[i_proc][i_part]; itarget < _targets_localization_idx_cpl[i_proc][i_part+1]; itarget++) {
            // Index in the interpolated Data array
            int interpInd = itarget;
            int iel = _targets_localization_data_cpl[itarget].lnum ;
            int ielP1 = iel+1;

            double value = 0.0;
            if(_targets_localization_data_cpl[itarget].distance != INFINITY ) {
              double x_target = _targets_localization_data_cpl[itarget].projectedX;
              double y_target = _targets_localization_data_cpl[itarget].projectedY;
              double z_target = _targets_localization_data_cpl[itarget].projectedZ;

              double tgtCoords[3]    = {x_target,y_target,z_target};
              double *barCoords      = NULL;
           /*   printf("iel %i itarget %i target %f %f %f _n_tot_target_cpl %i conneIDX \n",
              iel,itarget,x_target,y_target,z_target,_n_tot_target_cpl);*/

              if (weights_src_empty) {
                int    *_barCoordsIndex = NULL;
                double *_barCoords      = NULL;
                PDM_geom_elem_compute_polygon_barycentric_coordinates(1,
                                                                      &ielP1,
                                                                      tgtCoords,
                                                                      connecIdx,
                                                                      connec,
                                                                      coords,
                                                                      &_barCoordsIndex,
                                                                      &_barCoords
                                                                      );
                int n_elt = connecIdx[iel+1] - connecIdx[iel];
                _weights_src_idx[i_part][itarget+1] =
                  _weights_src_idx[i_part][itarget] + n_elt;

                if (s_weights_src <= _weights_src_idx[i_part][itarget+1]) {
                  s_weights_src *= 2;
                  _weights_src[i_part] = (double *)realloc((void *)(_weights_src[i_part]),
                    sizeof(double) * s_weights_src);
                }

                for (int i = 0; i < n_elt; i++) {
                  _weights_src[i_part][_weights_src_idx[i_part][itarget] +i] =
                    _barCoords[i];
                }

                free (_barCoordsIndex);
                free (_barCoords);
              }

              barCoords = &(_weights_src[i_part][_weights_src_idx[i_part][itarget]]);

              for (int k = 0; k < nComponent; k++) {
                value = 0.0;
                int k1=0;
                for (int i_vtx = connecIdx[iel]; i_vtx < connecIdx[iel+1]; i_vtx++) {
                   value +=  barCoords[k1++] * (*(double*)( (char*)referenceData + dataTypeSize * (nComponent * (connec[i_vtx]-1) + k) ) );
                }
                memcpy((char*)interpolatedData + dataTypeSize * ( nComponent * interpInd + k), &value, dataTypeSize);
              }//end k component loop

            }
            else {
              for (int k = 0; k < nComponent; k++) {
                value = 1000.0;
                memcpy((char*)interpolatedData + dataTypeSize * ( nComponent * interpInd + k), &value, dataTypeSize);
              }
            }
          } // loop on itarget
        } // if referenceFieldType == CWP_DOF_LOCATION_NODE || referenceFieldType == CWP_DOF_LOCATION_USER
        }  // else if(interpolationType == CWP_INTERPOLATION_TYPE_USER)
      } //Loop on i_proc
    } // loop on i_part

    referenceField -> sendBufferSet(interpolatedData);
    return interpolatedData;
  }





}//end_cwipi

/**
 * \endcond
 */
