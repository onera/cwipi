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
#include <geomLocation.hxx>
#include <geometry.hxx>
#include <mpi.h>
#include <pdm_mpi.h>
#include <pdm_mesh_nodal.h>
#include <pdm_dist_cloud_surf.h>
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

namespace cwipi {

  GeomLocation::GeomLocation()
  :Geometry::Geometry()
  {
  }
  
  GeomLocation::~GeomLocation()
  {
  }

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

  void* GeomLocation::interpolate(Field* referenceField) {

    int                 nComponent         = referenceField -> nComponentGet  ();
    CWP_Field_value_t   referenceFieldType = referenceField -> typeGet        ();
    CWP_Field_storage_t storage            = referenceField -> storageTypeGet ();
    void               *interpolatedData   = referenceField -> sendBufferGet  ();
    int                 dataTypeSize       = referenceField -> dataTypeSizeGet(); 
    
    if (interpolatedData != NULL) free(interpolatedData);
    interpolatedData = (void*) malloc( dataTypeSize * nComponent*_n_tot_target_cpl);


/*  |               proc 1             ||             proc 2              ||
    | part 1  | part 2 |  ... | part N || part 1  | part 2 |  ... | part N||
    | 1 2 3 4   6 7 8 9   ...             123  124  125 126   ...

*/
   
    for(int i_part=0;i_part<_nb_part;i_part++){
    
      void* referenceData = referenceField -> dataGet(i_part);
      // For a cell center field : give the value of the located cell
      
      for(int i_proc=0;i_proc<_n_ranks_g;i_proc++){

        if (referenceFieldType == CWP_FIELD_VALUE_CELL_POINT) {
          for (int itarget = _targets_localization_idx_cpl[i_proc][i_part]; itarget < _targets_localization_idx_cpl[i_proc][i_part+1]; itarget++) {
            //Index of the corresponding local reference Data.
            int iel = _targets_localization_data_cpl[itarget].lnum ;
            // Index in the interpolated Data array
            int interpInd = itarget;
          //  printf("iel %i itarget %i refData %f _n_tot_target_cpl %i\n",iel,itarget,referenceData[iel],_n_tot_target_cpl);
            for (int k = 0; k < nComponent; k++) {
              memcpy( interpolatedData + dataTypeSize * ( nComponent*interpInd + k ) ,
                      referenceData + dataTypeSize * ( nComponent*iel + k )          ,
                      dataTypeSize);
              //printf("interpolatedData[ nComponent*interpInd + k  ] %f i_part %i i_proc %i nComponent*interpInd + k %i nComponent*iel + k %i referenceData[nComponent*iel + k ] %f\n",
              //interpolatedData[ nComponent*interpInd + k  ],i_part,i_proc,nComponent*interpInd + k,nComponent*iel + k),referenceData[nComponent*iel + k ];
            }    
          } // loop on itarget
        } // if referenceFieldType == CWP_FIELD_VALUE_CELL_POINT

        if (referenceFieldType == CWP_FIELD_VALUE_NODE) {
          int* connecIdx = _mesh -> connecIdxGet(i_part);
          int* connec    = _mesh -> connecGet(i_part);
          
          int n_vtx      = _mesh -> getPartNVertex(i_part);
          int n_elts     = _mesh -> getPartNElts(i_part);
          double* coords = _mesh -> getVertexCoords(i_part);
        /*  printf("_targets_localization_idx_cpl[%i][%i] rank %i %i %i\n",
          i_proc,i_part,_rank,_targets_localization_idx_cpl[i_proc][i_part],_targets_localization_idx_cpl[i_proc][i_part+1]);
          */
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
               
              double tgtCoords[3] = {x_target,y_target,z_target};
              int *barCoordsIndex =NULL;
              double *barCoords = NULL;
           /*   printf("iel %i itarget %i target %f %f %f _n_tot_target_cpl %i conneIDX \n",
              iel,itarget,x_target,y_target,z_target,_n_tot_target_cpl);*/
              int conv_bar_compute = PDM_geom_elem_compute_polygon_barycentric_coordinates(1,
                                          	    		                         &ielP1,
                                                                			  tgtCoords,
                                                         		                  connecIdx,
                                                                   			  connec,
                                                           				  coords,
                                                  			                 &barCoordsIndex,
                                                                  			 &barCoords
                                                         				 );               
              for (int k = 0; k < nComponent; k++) {
                value = 0.0;
                for (int i_vtx = connecIdx[iel]; i_vtx < connecIdx[iel+1]; i_vtx++) {
                   value +=  barCoords[i_vtx - connecIdx[iel] ] * (*(double*)( referenceData + dataTypeSize * (nComponent * (connec[i_vtx]-1) + k) ) );
                }
                memcpy(interpolatedData + dataTypeSize * ( nComponent * interpInd + k), &value, dataTypeSize);
              }//end k component loop
              
              if(barCoordsIndex != NULL) free(barCoordsIndex);
              if(barCoords      != NULL) free(barCoords     );
            }
            else {
              for (int k = 0; k < nComponent; k++) {
                value = 1000.0;
                memcpy(interpolatedData + dataTypeSize * ( nComponent * interpInd + k), &value, dataTypeSize);
              }
            }
          } // loop on itarget
        } // if referenceFieldType == CWP_FIELD_VALUE_NODE
      } //Loop on i_proc
    } // loop on i_part
                
    referenceField -> sendBufferSet(interpolatedData);
    return interpolatedData;
  }
     
     
/*****************************************************************/

void GeomLocation::issend(Field* referenceField) {
      
      int  dataTypeSize       = referenceField -> dataTypeSizeGet(); 
      int nComponent = referenceField -> nComponentGet();

      int tag =referenceField -> fieldIDIntGet();
      void* dist_v_ptr = NULL;
      int dist_rank = -1;
      //printf("Avant interpolate |%s|\n",referenceFieldID.c_str());
      void* interpolatedFieldData = interpolate(referenceField);  
      //printf("Après interpolate |%s|\n",referenceFieldID.c_str());

        /* Loop on possibly intersecting distant ranks */
        /*---------------------------------------------*/
      for (int i_proc = 0; i_proc < _n_ranks_cpl; i_proc++) {    
        int n_points_dist = 0;
        
        int distant_rank = (*_connectableRanks_cpl)[i_proc];
        
        dist_v_ptr = interpolatedFieldData + dataTypeSize * nComponent * _targets_localization_idx_cpl[distant_rank][0];

        int request=-1;

        //printf("start _distantPartProcTargetIdx[0][%i] %i\n",i_proc,_targets_localization_idx_cpl[distant_rank][0]);

        int longueur = dataTypeSize * nComponent * (_targets_localization_idx_cpl[distant_rank][_nb_part]-_targets_localization_idx_cpl[distant_rank][0]);//_n_targets_dist_proc[i_proc];
        
        //printf("Send from %i to %i start %i longueur %i\n",_rank,distant_rank,nComponent*_targets_localization_idx_cpl[distant_rank][0],longueur);

        MPI_Issend(dist_v_ptr, longueur, MPI_BYTE, distant_rank, tag,
                   _globalComm,
                   &request);

        referenceField -> lastRequestAdd(i_proc,request);
      } /* End of loop on possibly intersecting ranks */
  }  
  
  
 
  void GeomLocation::locate_setting_surface(int* id_dist) {

    /* Paradigm mesh localisation _distance creation */
    *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_globalComm );
    
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
        Mesh* mesh_cpl = _geometry_cpl -> meshGet();
        int          n_target_cpl      = _geometry_cpl -> nTargetGet(i_part);
        CWP_g_num_t* gnum_target_cpl   = _geometry_cpl -> gnumTargetGet(i_part);
        double*      coords_target_cpl = _geometry_cpl -> coordsTargetGet(i_part);

        PDM_dist_cloud_surf_cloud_set (*id_dist,
                              0,
                              i_part,
                              n_target_cpl,
                              coords_target_cpl,
                              gnum_target_cpl
                             );   
      }//loop on part   
    } //end of if
  }



void GeomLocation::locate_setting_null(int* id_dist) {

    /*
     TODO: Intéressant pour la suite
     Inclure dans mesh_dist et autre la possibilité de définir plusieurs surfaces ...

   */

    /* Paradigm mesh localisation _distance creation */
    *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_globalComm );

    PDM_dist_cloud_surf_n_part_cloud_set(*id_dist,   0, _nb_part);  

    PDM_dist_cloud_surf_surf_mesh_global_data_set (*id_dist,
                                             _n_g_elt_cpl_over_part,
                                             _n_g_vtx_cpl_over_part,
                                             _nb_part_cpl);  
                                             
    //printf("ENULL %I64d %I64d \n",_n_g_elt_cpl_over_part,_n_g_vtx_cpl_over_part);

    for(int i_part =0;i_part<_nb_part;i_part++) {   
      int n_elt_null = 0;    
      double*      coords    = (double*)malloc(3*sizeof(double)*n_elt_null);
      CWP_g_num_t* gnum_elt  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*n_elt_null);

      PDM_dist_cloud_surf_cloud_set (*id_dist,
                              0,
                              i_part,
                              0,
                              coords,
                              gnum_elt
                             );
    }
 
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









  void GeomLocation::locate_setting_request(int* id_dist) {

    /*
     TODO: Intéressant pour la suite
     Inclure dans mesh_dist et autre la possibilité de définir plusieurs surfaces ...

   */

    /* Paradigm mesh localisation _distance creation */
    *id_dist   = PDM_dist_cloud_surf_create( PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_globalComm );

    PDM_dist_cloud_surf_n_part_cloud_set(*id_dist,   0, _nb_part);  

    PDM_dist_cloud_surf_surf_mesh_global_data_set (*id_dist,
                                             _n_g_elt_cpl_over_part,
                                             _n_g_vtx_cpl_over_part,
                                             _nb_part_cpl);  

    for(int i_part =0;i_part<_nb_part;i_part++) {   

      CWP_g_num_t* gnum_target          = gnumTargetGet  (i_part);
      double*      coords_target        = coordsTargetGet(i_part);

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
        Mesh* mesh_cpl = _geometry_cpl -> meshGet();
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
 }

  void GeomLocation::locate_compute(int id_dist) {
    PDM_dist_cloud_surf_compute(id_dist);
  }

 void GeomLocation::broadcasting_request(int* id_gnum_location) {
 
  *id_gnum_location = PDM_gnum_location_create(_nb_part_cpl,_nb_part, _pdm_globalComm);

  for(int i_part =0;i_part<_nb_part;i_part++) {    
    for(int i=0; i<_n_target[i_part]; i++) {
 //    if(_distance[i_part][i] == INFINITY ) {
  /*    printf("_closest_elt_gnum[%i][%i] rank %i %I64d coords %f %f %f _distance %f N %i\n",
      i_part,i,_rank,_closest_elt_gnum[i_part][i],
      _coords_target[i_part][3*i],_coords_target[i_part][3*i+1],_coords_target[i_part][3*i+2],
      _distance[i_part][i],
      _n_target[i_part]);
 */   // }
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
      CWP_g_num_t* gnum_target_cpl = _geometry_cpl -> gnumTargetGet(i_part);
      int          n_target_cpl    = _geometry_cpl -> nTargetGet   (i_part);

      Mesh* mesh_cpl = _geometry_cpl -> meshGet();
      CWP_g_num_t* gnum_elt_cpl = mesh_cpl -> GNumEltsGet(i_part);     
      int          n_elt_cpl    = mesh_cpl -> getPartNElts(i_part);
       
      PDM_gnum_location_elements_set(*id_gnum_location,i_part, n_elt_cpl,gnum_elt_cpl);    
    }//loop on part
  }//end if
 }


 void GeomLocation::broadcasting_set(int* id_gnum_location) {
 
  *id_gnum_location = PDM_gnum_location_create(_nb_part,_nb_part_cpl, _pdm_globalComm);

  for(int i_part =0;i_part<_nb_part;i_part++) {    

    CWP_g_num_t* gnum_target = gnumTargetGet(i_part);
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
      int n_target_cpl = _geometry_cpl -> nTargetGet(i_part);
      double*      coords_target        =  _geometry_cpl ->  coordsTargetGet(i_part);
      CWP_g_num_t* closest_elt          =  _geometry_cpl -> closestEltGnumGet  (i_part);
      double* dist_target               =  _geometry_cpl -> distanceTargetGet(i_part);
   //   CWP_g_num_t* toto= ((*_geometry_cpl)._closest_elt_gnum)[i_part];
     // double* tata= ((*_geometry_cpl)._distance)[i_part];
    //   for(int i=0; i< n_target_cpl; i++) {
        // if(dist_target[i] >= 1000.0 || gnum_target[i]>10000 || gnum_target[i]<0) {
           
   /*        printf("CPL_closest_elt_gnum[i_part][%i] rank %i %I64d %I64d coords %f %f %f _distance %f N %i\n",
           i,_rank,
           toto[i],
           closest_elt[i],
            coords_target[3*i],
            coords_target[3*i+1],
            coords_target[3*i+2],
            tata[i],
           n_target_cpl);*/
       //  }
        
       PDM_gnum_location_requested_elements_set(*id_gnum_location,i_part, n_target_cpl, _geometry_cpl -> closestEltGnumGet  (i_part));
   }//loop on part
  }//end if
 }
 /*************************************************************************************/
 /*************************************************************************************/
  /*************************************************************************************/


void GeomLocation::broadcasting_set_null(int* id_gnum_location) {
 
  *id_gnum_location = PDM_gnum_location_create(_nb_part,_nb_part_cpl, _pdm_globalComm);

  for(int i_part =0;i_part<_nb_part;i_part++) {    
    CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*0);
    PDM_gnum_location_elements_set(*id_gnum_location,i_part, 0,gnum_elt_null);      
  }

  for(int i_part =0; i_part<_nb_part_cpl; i_part++) {     
    CWP_g_num_t* gnum_elt_null  = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)*0);
    PDM_gnum_location_requested_elements_set(*id_gnum_location,i_part, 0,gnum_elt_null);
  }//loop on part
 }


  void GeomLocation::locate_get_cpl(int id_dist) {
    _geometry_cpl -> _distance           = (double**)malloc(sizeof(double*) * _nb_part_cpl);
    _geometry_cpl -> _projected          = (double**)malloc(sizeof(double*) * _nb_part_cpl);
    _geometry_cpl -> _closest_elt_gnum   = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*) * _nb_part_cpl);
    
    for(int i_part =0;i_part<_nb_part_cpl;i_part++) {     
      int          n_target_cpl    = _geometry_cpl -> nTargetGet(i_part);
      PDM_dist_cloud_surf_get (id_dist,
                         0,
                         i_part,
                         &(_geometry_cpl -> _distance[i_part]),
                         &(_geometry_cpl -> _projected[i_part]),
                         &(_geometry_cpl -> _closest_elt_gnum[i_part]));
      
      for(int i=0;i<n_target_cpl;i++){
        if(_geometry_cpl -> _closest_elt_gnum[i_part][i]>CWP_g_num_t(_n_g_elt_over_part) 
           || _geometry_cpl -> _closest_elt_gnum[i_part][i]<CWP_g_num_t(1) 
           || _geometry_cpl -> _distance[i_part][i]>0.01){
          _geometry_cpl -> _closest_elt_gnum[i_part][i]=CWP_g_num_t(1);
          _geometry_cpl -> _distance[i_part][i]=INFINITY;
        }       
      }         
    } //end loop on i_part
  }// End locate_cell_point




  void GeomLocation::locate_get(int id_dist) {
  
    _distance           = (double**)malloc(sizeof(double*) * _nb_part);
    _projected          = (double**)malloc(sizeof(double*) * _nb_part);
    _closest_elt_gnum   = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*) * _nb_part);

    for(int i_part =0;i_part<_nb_part;i_part++) {     
    
        PDM_dist_cloud_surf_get (id_dist,
                         0,
                         i_part,
                         &(_distance [i_part]),
                         &(_projected[i_part]),
                         &(_closest_elt_gnum[i_part]));
                         

       for(int i=0;i<_n_target[i_part];i++){
         if(_closest_elt_gnum[i_part][i]>CWP_g_num_t(_n_g_elt_cpl_over_part) || _closest_elt_gnum[i_part][i]<CWP_g_num_t(1) || _distance [i_part][i]>0.1){
           _closest_elt_gnum[i_part][i]=CWP_g_num_t(1);
           _distance [i_part][i]=INFINITY;
           
         }       
       }         
    }
 }// End locate_cell_point


 /*************************************************************************************/
 /*************************************************************************************/
  /*************************************************************************************/
 
 void GeomLocation::location_compute(int id_gnum_location) {
    PDM_gnum_location_compute(id_gnum_location);  
 }


 void GeomLocation::location_get(int id_gnum_location) {

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



  void GeomLocation::location_get_cpl(int id_gnum_location) {
    _geometry_cpl ->_target_proc_part_num_idx =(int**)malloc(sizeof(int*)*_nb_part_cpl);
    _geometry_cpl ->_target_proc_part_num     =(int**)malloc(sizeof(int*)*_nb_part_cpl);

    for(int i_part =0;i_part<_nb_part_cpl;i_part++) {     
      PDM_gnum_location_get(id_gnum_location,
                            i_part,
                            &(_geometry_cpl ->_target_proc_part_num_idx[i_part]),
                            &(_geometry_cpl ->_target_proc_part_num[i_part])
                            );  
    }
  }// End locate_cell_point

 void GeomLocation::filling_of_broadcasting_array() {


   PDM_timer_t *t1 = PDM_timer_create();
   PDM_timer_init(t1);
   PDM_timer_resume(t1);  

   if(_targets_localization_idx==NULL) {
    _targets_localization_idx   =(int**)malloc(sizeof(int*)*_n_ranks_g);
    for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) 
      _targets_localization_idx [i_proc] = NULL;
   }
   _localization_count_comm_proc =(int**)malloc(sizeof(int*)*_n_ranks_g);

   for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
     if(_targets_localization_idx [i_proc] == NULL) 
       _targets_localization_idx [i_proc]=(int*)malloc(sizeof(int)*(1+_nb_part_cpl));
       _localization_count_comm_proc [i_proc]=(int*)malloc(sizeof(int)*(1+_nb_part_cpl));
     for (int i_part = 0; i_part < _nb_part_cpl+1; i_part++) {
       _targets_localization_idx [i_proc][i_part]=0;
       _localization_count_comm_proc [i_proc][i_part]=0;
     }
   }

   for (int i_part = 0; i_part < _nb_part; i_part++) {
     for(int k=0;k<_n_target[i_part];k++){
  //     printf("_target_proc_part_num_idx[i_part][%i] rank %i %i\n",k,_rank,_target_proc_part_num_idx[i_part][k]);
  //     printf("_target_proc_part_num[i_part][%i] rank %i %i\n",k,_rank,_target_proc_part_num[i_part][ _target_proc_part_num_idx[i_part][k] ]);
       int elt_proc = _target_proc_part_num[i_part][ _target_proc_part_num_idx[i_part][k] ];
       int elt_part = _target_proc_part_num[i_part][ _target_proc_part_num_idx[i_part][k] + 1];
      
       _targets_localization_idx [elt_proc][elt_part]++;     
       _localization_count_comm_proc [elt_proc][elt_part]++;   
    }
  }//end i_part

  _transform_to_index(_targets_localization_idx,_n_ranks_g,_nb_part_cpl);

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
 
 /*        printf("_targets_localization_data rank %i [%i] lnum %i origin_part %i closest_elt_gnum %i elt_proc %i elt_part %i distance %f proj %f\n",
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
    
       PDM_timer_hang_on(t1);
   double elapsAlltoAll = PDM_timer_elapsed(t1);
   double elapsAlltoAllMax = 0.0;   
   MPI_Reduce(&elapsAlltoAll, &elapsAlltoAllMax, 1, MPI_DOUBLE, MPI_MAX,0,_connectableComm);
   int rankconn;
   MPI_Comm_rank(_connectableComm,&rankconn);
   if (rankconn==0) {
        printf ("Filling Time : %12.5e \n",elapsAlltoAllMax);
   }
    
    
 }

 void GeomLocation:: initialization_of_reception_array() {
   if(_targets_localization_idx_cpl == NULL){
     _targets_localization_idx_cpl = (int**)malloc(  sizeof(int*) *_n_ranks_g);      
     for(int i =0; i<_n_ranks_g;i++) 
       _targets_localization_idx_cpl[i] = NULL;
   }
  
   for(int i =0; i<_n_ranks_g;i++) {
     if(_targets_localization_idx_cpl[i] == NULL) 
       _targets_localization_idx_cpl[i] = (int*)malloc(sizeof(int)*(1+_nb_part)); 
     for (int i_part = 0; i_part < _nb_part+1; i_part++) {
       _targets_localization_idx_cpl[i][i_part]=0;
     }
   }
   
 } 


 void GeomLocation::broadcasting_index_null() {
 
    int tag = -1;
    if(_geometryLocation == CWP_FIELD_VALUE_CELL_POINT)
      tag = 1520;
    else if(_geometryLocation == CWP_FIELD_VALUE_NODE)
      tag = 1530;

   std::vector<int> connectable = *_connectableRanks_cpl;
   int id     = _localCodeProperties -> idGet();
   int id_cpl = _coupledCodeProperties -> idGet();

   int* sbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part_cpl);
   for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
     for(int i_part=0;i_part<_nb_part_cpl;i_part++)
       sbuffer[ i_proc * _nb_part_cpl + i_part ] = 0;

   MPI_Barrier(_globalComm);

   int* recvbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part_cpl);
   MPI_Request request;
   
   MPI_Ialltoall(sbuffer, _nb_part_cpl, MPI_INT, 
                 recvbuffer, _nb_part_cpl, MPI_INT,
                 _globalComm,&request);
   
   MPI_Status stat;
   MPI_Wait(&request,&stat);
   
   free(sbuffer   );
   free(recvbuffer   );
   
 } 

 void GeomLocation::broadcasting_index_communication_async() {
 } 

 void GeomLocation::both_index_communication() {

    int tag = -1;
    if(_geometryLocation == CWP_FIELD_VALUE_CELL_POINT)
      tag = 1520;
    else if(_geometryLocation == CWP_FIELD_VALUE_NODE)
      tag = 1530;

   int* sbuffer = (int*) malloc(sizeof(int)*_n_ranks_g * _nb_part );
   for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
     for(int i_part=0;i_part< _nb_part; i_part++) {
       sbuffer[ i_proc * _nb_part_cpl + i_part ] = (_geometry_cpl -> _localization_count_comm_proc)[i_proc][i_part];
     }

   int* recvbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part);
   MPI_Barrier(_globalComm);
   
   PDM_timer_t *t1 = PDM_timer_create();
   PDM_timer_init(t1);
   PDM_timer_resume(t1);  
   
   MPI_Request sreq;
   
   MPI_Ialltoall(sbuffer, _nb_part, MPI_INT, 
                 recvbuffer, _nb_part, MPI_INT,
                 _globalComm,&sreq);     

   
   MPI_Status stat;
   MPI_Wait(&sreq,&stat);

   PDM_timer_hang_on(t1);
   double elapsAlltoAll = PDM_timer_elapsed(t1);
   double elapsAlltoAllMax = 0.0;   
   MPI_Reduce(&elapsAlltoAll, &elapsAlltoAllMax, 1, MPI_DOUBLE, MPI_MAX,0,_connectableComm);
   int rankconn;
   MPI_Comm_rank(_connectableComm,&rankconn);
   
   
   for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
     for(int i_part=0;i_part<_nb_part;i_part++) 
       _targets_localization_idx_cpl[i_proc][i_part] = recvbuffer[ i_proc * _nb_part + i_part ];   
   
   if (rankconn==0) {
        printf ("AlltoAllSend Time : %12.5e %i\n",elapsAlltoAllMax,_n_ranks_g);
   }
   
   free(sbuffer   );
   free(recvbuffer   );
 } 



 void GeomLocation::broadcasting_index_communication() {
 
 
    int tag = -1;
    if(_geometryLocation == CWP_FIELD_VALUE_CELL_POINT)
      tag = 1520;
    else if(_geometryLocation == CWP_FIELD_VALUE_NODE)
      tag = 1530;

   int id     = _localCodeProperties -> idGet();
   int id_cpl = _coupledCodeProperties -> idGet();

   int* sbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part_cpl);
   for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
     for(int i_part=0;i_part<_nb_part_cpl;i_part++) {
       sbuffer[ i_proc * _nb_part_cpl + i_part ] = _localization_count_comm_proc[i_proc][i_part];
     }

   int* recvbuffer_trash = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part);
   MPI_Barrier(_globalComm);
   
   PDM_timer_t *t1 = PDM_timer_create();
   PDM_timer_init(t1);
   PDM_timer_resume(t1);  
   
   MPI_Request sreq;
   
   MPI_Ialltoall(sbuffer, _nb_part_cpl, MPI_INT, 
                 recvbuffer_trash, _nb_part, MPI_INT,
                 _globalComm,&sreq);    
                   
   
   MPI_Status stat;
   MPI_Wait(&sreq,&stat);

   PDM_timer_hang_on(t1);
   double elapsAlltoAll = PDM_timer_elapsed(t1);
   double elapsAlltoAllMax = 0.0;   
   MPI_Reduce(&elapsAlltoAll, &elapsAlltoAllMax, 1, MPI_DOUBLE, MPI_MAX,0,_connectableComm);
   int rankconn;
   MPI_Comm_rank(_connectableComm,&rankconn);
   
   
   if (rankconn==0) {
        printf ("AlltoAllSend Time : %12.5e %i\n",elapsAlltoAllMax,_n_ranks_g);
   }
   
   free(sbuffer   );
   free(recvbuffer_trash   );
 } 


 void GeomLocation::reception_index_communication() {
    int tag = -1;
    if(_geometryLocation == CWP_FIELD_VALUE_CELL_POINT)
      tag = 1520;
    else if(_geometryLocation == CWP_FIELD_VALUE_NODE)
      tag = 1530;

   std::vector<int> connectable = *_connectableRanks_cpl;
   int id     = _localCodeProperties -> idGet();
   int id_cpl = _coupledCodeProperties -> idGet();

   int* recvbuffer = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part);
   for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
     for(int i_part=0;i_part<_nb_part;i_part++) {
       recvbuffer[ i_proc * _nb_part + i_part ] = 0;   
     }
   
   int *sendbuffer_trash = (int*) malloc(sizeof(int)*_n_ranks_g*_nb_part_cpl);
   for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
     for(int i_part=0;i_part<_nb_part_cpl;i_part++) {
       sendbuffer_trash[ i_proc * _nb_part_cpl + i_part ] = 0;   
     }
     
   MPI_Barrier(_globalComm);
   
   PDM_timer_t *t1 = PDM_timer_create();
   PDM_timer_init(t1);
   PDM_timer_resume(t1);  
   
   MPI_Request rreq;
   
   MPI_Ialltoall(sendbuffer_trash, _nb_part_cpl, MPI_INT, 
                 recvbuffer, _nb_part, MPI_INT,
                 _globalComm,&rreq);   

   MPI_Status stat;
   MPI_Wait(&rreq,&stat);

   PDM_timer_hang_on(t1);
   double elapsAlltoAll = PDM_timer_elapsed(t1);
   double elapsAlltoAllMax = 0.0;   
   MPI_Reduce(&elapsAlltoAll, &elapsAlltoAllMax, 1, MPI_DOUBLE, MPI_MAX,0,_connectableComm);
   
   int rankconn;
   MPI_Comm_rank(_connectableComm,&rankconn);
   
   if (rankconn==0) {
        printf ("AlltoAllRecv Time : %12.5e \n",elapsAlltoAllMax);
   }
   
   for(int i_proc=0;i_proc<_n_ranks_g;i_proc++)
     for(int i_part=0;i_part<_nb_part;i_part++) {
       _targets_localization_idx_cpl[i_proc][i_part] = recvbuffer[ i_proc * _nb_part + i_part ];   
     }

   if(sendbuffer_trash!=NULL) free(sendbuffer_trash   );
   if(recvbuffer!=NULL) free(recvbuffer   );
 } 



 void GeomLocation::reception_index_communication_async() {
 }
 
 
 void GeomLocation::prepare_data_communication_send() {
  
   for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) 
     free(_localization_count_comm_proc[i_proc] );
   free(_localization_count_comm_proc);
  
   _localization_count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
   _localization_disp_send  = (int*)malloc(sizeof(int)*_n_ranks_g); 
 
   for (int i= 0; i < _n_ranks_g; i++) { 
     _localization_count_send[i] = sizeof(target_data)*(_targets_localization_idx[i][_nb_part_cpl] - _targets_localization_idx[i][0]);
     _localization_disp_send [i] = sizeof(target_data)*_targets_localization_idx[i][0];
   }

 }


 void GeomLocation::prepare_data_communication_recv() {
  
   _transform_to_index(_targets_localization_idx_cpl,_n_ranks_g,_nb_part);
  
   if(_targets_localization_data_cpl!=NULL) free(_targets_localization_data_cpl);
   _targets_localization_data_cpl = (target_data*)malloc(sizeof(target_data)*_targets_localization_idx_cpl[_n_ranks_g-1][_nb_part]);        

   _localization_count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);         
   _localization_disp_recv  = (int*)malloc(sizeof(int)*_n_ranks_g);         
 
   for (int i= 0; i < _n_ranks_g; i++) { 
     _localization_count_recv[i] = sizeof(target_data)*(_targets_localization_idx_cpl[i][_nb_part] - _targets_localization_idx_cpl[i][0]); 
     _localization_disp_recv [i] = sizeof(target_data)*_targets_localization_idx_cpl[i][0];
   }

 }
 
 void GeomLocation::data_communication_send() {
 
 
 
   int* count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);         
   int* disp_recv  = (int*)malloc(sizeof(int)*_n_ranks_g);         
 
   for (int i= 0; i < _n_ranks_g; i++) { 
     count_recv[i] = 0; 
     disp_recv [i] = 0;
   }
   void* recv_buffer_trash = NULL;
  
   MPI_Request sreq;   
   MPI_Ialltoallv((void*)_targets_localization_data, _localization_count_send, _localization_disp_send, MPI_BYTE, 
                  recv_buffer_trash, count_recv,     disp_recv, MPI_BYTE,
                  _globalComm,&sreq);     

   
   MPI_Status stat;
   MPI_Wait(&sreq,&stat);
   
   free(count_recv);
   free(disp_recv );
  /*
  _IAlltoallIndexSend((void*)_targets_localization_data, _localization_count_send, _localization_disp_send,
                  MPI_BYTE,
                  _globalComm, *_connectableRanks_cpl);*/
  }

 void GeomLocation::data_communication_recv() {
   int* count_send = (int*)malloc(sizeof(int)*_n_ranks_g);         
   int* disp_send  = (int*)malloc(sizeof(int)*_n_ranks_g);         
 
   for (int i= 0; i < _n_ranks_g; i++) { 
     count_send[i] = 0; 
     disp_send [i] = 0;
   }
  
     void* sbuffer_trash = NULL;
  
   MPI_Request rreq;   
   
   MPI_Ialltoallv(sbuffer_trash                       ,  count_send              , disp_send              , MPI_BYTE, 
                 (void*)_targets_localization_data_cpl,  _localization_count_recv, _localization_disp_recv, MPI_BYTE,
                 _globalComm,&rreq);     
   
   MPI_Status stat;
   MPI_Wait(&rreq,&stat); 
   free(count_send);
   free(disp_send );
 
 
 /* 
  _IAlltoallIndexRecv((void*)_targets_localization_data_cpl, _localization_count_recv, _localization_disp_recv, 
                  MPI_BYTE,
                  _globalComm, *_connectableRanks_cpl);
                  
  */
                  
  }


 void GeomLocation::both_data_communication() {
   
  
   MPI_Request req;   
   MPI_Ialltoallv((void*) _geometry_cpl -> _targets_localization_data, _geometry_cpl -> _localization_count_send, _geometry_cpl -> _localization_disp_send, MPI_BYTE,  
                  (void*)_targets_localization_data_cpl,  _localization_count_recv, _localization_disp_recv, MPI_BYTE,
                 _globalComm,&req);     
   
   MPI_Status stat;
   MPI_Wait(&req,&stat); 
   
  }


 void GeomLocation::data_communication_null() {
   
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
                 _globalComm,&rreq);     
   
   MPI_Status stat;
   MPI_Wait(&rreq,&stat); 
   
   free(count_send);
   free(disp_send );
   free(count_recv);
   free(disp_recv );
 
   
 
 /* 
  _IAlltoallIndexRecv((void*)_targets_localization_data_cpl, _localization_count_recv, _localization_disp_recv, 
                  MPI_BYTE,
                  _globalComm, *_connectableRanks_cpl);
                  
  */
                  
  }
  
 void GeomLocation::data_communication_wait_send() {  
  // _WaitSend(*_connectableRanks_cpl,&_send_requests);
  
   free(_localization_count_send);
   free(_localization_disp_send );
 }


 void GeomLocation::data_communication_wait_recv() {  
 //  _WaitRecv(*_connectableRanks_cpl,&_recv_requests);
  
   free(_localization_count_recv);
   free(_localization_disp_recv );
  
   _n_tot_target_cpl    = _targets_localization_idx_cpl[_n_ranks_g-1][_nb_part];
 }
 
/*******************************************************************************************************************/
/*******************************************************************************************************************/
/*******************************************************************************************************************/


  
}//end_cwipi
  

