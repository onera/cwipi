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

  SpatialInterpLocation::SpatialInterpLocation()
  :SpatialInterp::SpatialInterp()
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
          PDM_dist_cloud_surf_free (_id_dist);
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
          PDM_dist_cloud_surf_free (_id_dist);
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
        PDM_dist_cloud_surf_free (_id_dist);
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

  void SpatialInterpLocation::init(Coupling *coupling, CWP_Dof_location_t pointsCloudLocation, int slave) {
    SpatialInterp::init(coupling, pointsCloudLocation, slave);

    _id_dist = -1;

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
      *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_MESH_SETTED,
                                               1,
                                               _pdm_cplComm,
                                               PDM_OWNERSHIP_UNGET_RESULT_IS_FREE );

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

      *id_dist = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm);

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
      *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm,
                                               PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

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

      *id_dist = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm);

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
      *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm,
                                               PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

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

      *id_dist = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm);

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
      *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm,
                                               PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

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

      *id_dist = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED, 1, _pdm_cplComm);

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
