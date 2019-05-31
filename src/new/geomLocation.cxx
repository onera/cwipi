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
#include <pdm_mesh_dist.h>
#include <pdm_gnum.h>
#include <pdm_gnum_location.h>
#include <bftc_error.h>
#include <bftc_printf.h>
#include "cwp.h"
#include <limits>

#include <cmath>

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

  double* GeomLocation::interpolate(Field <double>* referenceField) {

    int                 nComponent         = referenceField -> nComponentGet  ();
    CWP_Field_value_t   referenceFieldType = referenceField -> typeGet        ();
    CWP_Field_storage_t storage            = referenceField -> storageTypeGet ();
    double             *interpolatedData   = referenceField -> sendBufferGet  ();

    
    if (interpolatedData == NULL) interpolatedData = (double*) malloc(sizeof(double)/**nComponent*/*_n_tot_target_cpl);


/*  |               proc 1             ||             proc 2              ||
    | part 1  | part 2 |  ... | part N || part 1  | part 2 |  ... | part N||
    | 1 2 3 4   6 7 8 9   ...             123  124  125 126   ...

*/
   
    for(int i_part=0;i_part<_nb_part;i_part++){
      double* referenceData = referenceField -> dataGet(i_part);
      // For a cell center field : give the value of the located cell
      
      for(int i_proc=0;i_proc<_n_ranks_g;i_proc++){

        if (referenceFieldType == CWP_FIELD_VALUE_CELL_POINT) {
          
          for (int itarget = _targets_cpl_idx[i_proc][i_part]; itarget < _targets_cpl_idx[i_proc][i_part+1]; itarget++) {

            //Index of the corresponding local reference Data.
            
            int iel = _targets_cpl[itarget].lnum ;
            // Index in the interpolated Data array
            int interpInd = itarget;
            printf("iel %i itarget %i refData %f _n_tot_target_cpl %i\n",iel,itarget,referenceData[iel],_n_tot_target_cpl);
            for (int k = 0; k < nComponent; k++)
              {interpolatedData[ nComponent*interpInd + k  ] = referenceData[nComponent*iel + k ];
               //printf("interpolatedData[ nComponent*interpInd + k  ] %f i_part %i i_proc %i nComponent*interpInd + k %i nComponent*iel + k %i referenceData[nComponent*iel + k ] %f\n",
               //interpolatedData[ nComponent*interpInd + k  ],i_part,i_proc,nComponent*interpInd + k,nComponent*iel + k),referenceData[nComponent*iel + k ];
              }    
            } // loop on itarget
         } // if referenceFieldType == CWP_FIELD_VALUE_CELL_POINT

        if (referenceFieldType == CWP_FIELD_VALUE_NODE) {
        
          int* connecIdx = _mesh -> connecIdxGet(i_part);
          int* connec = _mesh -> connecGet(i_part);
          
          int n_vtx  = _mesh -> getPartNVertex(i_part);
          int n_elts  = _mesh -> getPartNElts(i_part);
          double* coords = _mesh -> getVertexCoords(i_part);

           printf("_targets_cpl_idx[%i][%i] rank %i %i %i\n",i_proc,i_part,_rank,_targets_cpl_idx[i_proc][i_part],_targets_cpl_idx[i_proc][i_part+1]);
           for (int itarget = _targets_cpl_idx[i_proc][i_part]; itarget < _targets_cpl_idx[i_proc][i_part+1]; itarget++) {
             //Index of the corresponding local reference Data.
            
                 
         
             // Index in the interpolated Data array
             int interpInd = itarget;
             int iel = _targets_cpl[itarget].lnum ;
             
             if(_targets_cpl[itarget].distance != INFINITY ) {

               double x_target = _targets_cpl[itarget].projectedX;
               double y_target = _targets_cpl[itarget].projectedY;
               double z_target = _targets_cpl[itarget].projectedZ;    
                              
              double sum_vtx = 0.0;
             // printf("iel %i itarget %i refData %f _n_tot_target_cpl %i\n",iel,itarget,referenceData[iel],_n_tot_target_cpl);

              interpolatedData[ interpInd  ] = 0.0;
              
              for (int i_vtx = connecIdx[iel]; i_vtx < connecIdx[iel+1]; i_vtx++) {

               double x_vtx = coords[3*(connec[i_vtx]-1)  ];
               double y_vtx = coords[3*(connec[i_vtx]-1)+1];
               double z_vtx = coords[3*(connec[i_vtx]-1)+2];
               
               double bx = abs(x_target - x_vtx);
               double by = abs(y_target - y_vtx);
               double bz = abs(z_target - z_vtx);
               
               
               double dvtx = sqrt(bx*bx +by*by + bz*bz);
               sum_vtx += dvtx;
                
               //Interpolation linéaire
               
               interpolatedData[ interpInd  ] += dvtx * referenceData[connec[i_vtx]-1] ;
               }
               if(sum_vtx!=0) interpolatedData[ interpInd  ] /= sum_vtx;
             //  interpolatedData[ interpInd  ] = coords[3*(connec[connecIdx[iel]]-1)  ];
              }
              else {
                interpolatedData[ interpInd  ] = 1000.0;
              }
            } // loop on itarget
         } // if referenceFieldType == CWP_FIELD_VALUE_NODE
         
       } //Loop on i_proc

    } // loop on i_part


    referenceField -> sendBufferSet(interpolatedData);
    return interpolatedData;
  }
     
     
/*****************************************************************/

void GeomLocation::issend(Field <double>* referenceField) {
      
      int size = sizeof(double);
      int nComponent = referenceField -> nComponentGet();

      int tag =atoi((referenceField -> fieldIDGet()).c_str());
      double* dist_v_ptr = NULL;
      int dist_rank = -1;
      //On va supposer pour le moment que les données contenues dans interpolatedFieldData
      // sont contigues en mémoire 
      //printf("Avant interpolate |%s|\n",referenceFieldID.c_str());
      double* interpolatedFieldData;
     
      interpolatedFieldData = interpolate(referenceField);  
      
      //printf("Après interpolate |%s|\n",referenceFieldID.c_str());

        /* Loop on possibly intersecting distant ranks */
        /*---------------------------------------------*/
      for (int i_proc = 0; i_proc < _n_ranks_cpl; i_proc++) {    
        int n_points_dist = 0;
        std::vector<double*> v_interpolatedFieldData(_nb_part,NULL);
        
        int distant_rank = (*_connectableRanks_cpl)[i_proc];
        
        dist_v_ptr = &(interpolatedFieldData[nComponent*_targets_cpl_idx[distant_rank][0]]);

        int request=-1;

        printf("start _distantPartProcTargetIdx[0][%i] %i\n",i_proc,_targets_cpl_idx[distant_rank][0]);

        int longueur = nComponent * (_targets_cpl_idx[distant_rank][_nb_part]-_targets_cpl_idx[distant_rank][0]);//_n_targets_dist_proc[i_proc];
        
        printf("Send from %i to %i start %i longueur %i\n",_rank,distant_rank,nComponent*_targets_cpl_idx[distant_rank][0],longueur);

        MPI_Issend(dist_v_ptr, longueur, MPI_DOUBLE, distant_rank, tag,
                   _globalComm,
                   &request);

        referenceField -> lastRequestAdd(i_proc,request);

        //On met dans v_interpolatedFieldData le début de chaque bout de partition
        //associé au proc i_proc ça correspond à l'envoie request
        //referenceField -> localInterpFieldStorage(request,v_interpolatedFieldData);

        } /* End of loop on possibly intersecting ranks */

  }  
  
  
 
  void GeomLocation::locate_setting_surface(int id_dist) {

    /* Paradigm mesh localisation _distance creation */

    id_dist   = PDM_mesh_dist_create( PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_globalComm );

    PDM_mesh_dist_n_part_cloud_set(id_dist,   0, _nb_part_cpl);  

    PDM_mesh_dist_surf_mesh_global_data_set (id_dist,
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

      PDM_mesh_dist_surf_mesh_part_set (id_dist,
                                      i_part,
                                      n_elts,
                                      connecIdx,
                                      connec,
                                      gnum_elt,
                                      n_vtx,
                                      coords,
                                      gnum_vtx);     
    }
 
 
    for(int i_part =0; i_part<_nb_part_cpl; i_part++) {     

    
      if(_both_codes_are_local == 0) {
      
        int* tmp1,*tmp2,*tmp5;  
        CWP_g_num_t* tmp3,*tmp4;
        double* tmp6; 
 
 
 
        int          n_target      =  nTargetGet(i_part);
        CWP_g_num_t* gnum_target   =  gnumTargetGet(i_part);
        double*      coords_target =  coordsTargetGet(i_part); 
        
        PDM_mesh_dist_cloud_set (id_dist,
                              0,
                              i_part,
                              0,
                              NULL ,
                              NULL
                             );                              
      }
      else {
        Mesh* mesh_cpl = _geometry_cpl -> meshGet();
        int          n_target_cpl      = _geometry_cpl -> nTargetGet(i_part);
        CWP_g_num_t* gnum_target_cpl   = _geometry_cpl -> gnumTargetGet(i_part);
        double*      coords_target_cpl = _geometry_cpl -> coordsTargetGet(i_part);

        PDM_mesh_dist_cloud_set (id_dist,
                              0,
                              i_part,
                              n_target_cpl,
                              coords_target_cpl,
                              gnum_target_cpl
                             );   
     }   
   }
  
  }

  void GeomLocation::locate_setting_request(int id_dist) {

    /*
     TODO: Intéressant pour la suite
     Inclure dans mesh_dist et autre la possibilité de définir plusieurs surfaces ...

   */

    /* Paradigm mesh localisation _distance creation */
    id_dist   = PDM_mesh_dist_create( PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_globalComm );

    PDM_mesh_dist_n_part_cloud_set(id_dist,   0, _nb_part);  

    PDM_mesh_dist_surf_mesh_global_data_set (id_dist,
                                           _n_g_elt_cpl_over_part,
                                           _n_g_vtx_cpl_over_part,
                                           _nb_part_cpl);  

    for(int i_part =0;i_part<_nb_part;i_part++) {   

      CWP_g_num_t* gnum_target          = gnumTargetGet  (i_part);
      double*      coords_target        = coordsTargetGet(i_part);

      PDM_mesh_dist_cloud_set (id_dist,
                              0,
                              i_part,
                              _n_target[i_part],
                              coords_target,
                              gnum_target
                             );
    }
 
 
    for(int i_part =0; i_part<_nb_part_cpl; i_part++) {     
    
      if(_both_codes_are_local == 0) {
        int*         connecIdx = _mesh -> connecIdxGet(i_part);
        int*         connec    = _mesh -> connecGet(i_part);
        
        double*      coords    = _mesh -> getVertexCoords(i_part);
        CWP_g_num_t* gnum_vtx  = _mesh -> getVertexGNum(i_part);
        CWP_g_num_t* gnum_elt  = _mesh -> GNumEltsGet(i_part);     

  /* TODO: Plante quand on met les pointeurs à NULL dans surf_mesh_part_set ou quand on utilise les tmp.
  */
        
        PDM_mesh_dist_surf_mesh_part_set (id_dist,
                                          i_part,
                                          0,
                                          connecIdx,
                                          connec,
                                          gnum_elt,
                                          0,
                                          coords,
                                          gnum_vtx);                                 
      }
      else {
        Mesh* mesh_cpl = _geometry_cpl -> meshGet();
        int*         connecIdx_cpl = mesh_cpl -> connecIdxGet(i_part);
        int*         connec_cpl    = mesh_cpl -> connecGet(i_part);

        int          n_vtx_cpl     = mesh_cpl -> getPartNVertex(i_part);
        int          n_elts_cpl    = mesh_cpl -> getPartNElts(i_part);
        double*      coords_cpl    = mesh_cpl -> getVertexCoords(i_part);
        CWP_g_num_t* gnum_vtx_cpl  = mesh_cpl -> getVertexGNum(i_part);
        CWP_g_num_t* gnum_elt_cpl  = mesh_cpl -> GNumEltsGet(i_part);     
  
        PDM_mesh_dist_surf_mesh_part_set (id_dist,
                                          i_part,
                                          n_elts_cpl,
                                          connecIdx_cpl,
                                          connec_cpl,
                                          gnum_elt_cpl,
                                          n_vtx_cpl,
                                          coords_cpl,
                                          gnum_vtx_cpl);      
      
       }   
     }
 }

  void GeomLocation::locate_compute(int id_dist) {
    PDM_mesh_dist_compute(id_dist);
  }

 void GeomLocation::broadcasting_request(int id_gnum_location) {
 
  id_gnum_location = PDM_gnum_location_create(_nb_part,_nb_part_cpl, _pdm_globalComm);

  for(int i_part =0;i_part<_nb_part;i_part++) {    
    for(int i=0; i<_n_target[i_part]; i++) {
     if(_distance[i_part][i] == INFINITY ) {
 /*     printf("_closest_elt_gnum[i_part][%i] rank %i %I64d coords %f %f %f _distance %f N %i\n",
      i,_rank,_closest_elt_gnum[i_part][i],
      _coords_target[i_part][3*i],_coords_target[i_part][3*i+1],_coords_target[i_part][3*i+2],
      _distance[i_part][i],
      _n_target[i_part]);*/
     }
    } 
    PDM_gnum_location_requested_elements_set(id_gnum_location,i_part, _n_target[i_part],_closest_elt_gnum[i_part]);   
  }
 
  for(int i_part =0; i_part<_nb_part_cpl; i_part++) {     
    
      if(_both_codes_are_local == 0) {
        CWP_g_num_t* gnum_elt = _mesh -> GNumEltsGet(i_part);     
        PDM_gnum_location_elements_set(id_gnum_location,i_part, 0,gnum_elt);    
      }     
     else {  
       CWP_g_num_t* gnum_target_cpl = _geometry_cpl -> gnumTargetGet(i_part);
       int          n_target_cpl    = _geometry_cpl -> nTargetGet   (i_part);

       Mesh* mesh_cpl = _geometry_cpl -> meshGet();
       CWP_g_num_t* gnum_elt_cpl = mesh_cpl -> GNumEltsGet(i_part);     
       int          n_elt_cpl    = mesh_cpl -> getPartNElts(i_part);
       
       PDM_gnum_location_elements_set(id_gnum_location,i_part, n_elt_cpl,gnum_elt_cpl);    
     }
  }
 }


 void GeomLocation::broadcasting_set(int id_gnum_location) {
 
  id_gnum_location = PDM_gnum_location_create(_nb_part,_nb_part_cpl, _pdm_globalComm);

  for(int i_part =0;i_part<_nb_part_cpl;i_part++) {    

    CWP_g_num_t* gnum_target = gnumTargetGet(i_part);
    CWP_g_num_t* gnum_elt = _mesh -> GNumEltsGet(i_part);     
    
    
    PDM_gnum_location_elements_set(id_gnum_location,i_part, _n_elt[i_part],gnum_elt);      
  }

  for(int i_part =0; i_part<_nb_part_cpl; i_part++) {     
    
      if(_both_codes_are_local == 0) {
      
        CWP_g_num_t* gnum_target = gnumTargetGet(i_part);
        CWP_g_num_t* gnum_elt = _mesh -> GNumEltsGet(i_part);     
        PDM_gnum_location_requested_elements_set(id_gnum_location,i_part, 0,gnum_elt);
      }     
     else {
       int n_target_cpl = _geometry_cpl -> nTargetGet(i_part);
       double*      coords_target        =  _geometry_cpl ->  coordsTargetGet(i_part);
       CWP_g_num_t* closest_elt          =  _geometry_cpl -> closestEltGnumGet  (i_part);
       double* dist_target               =  _geometry_cpl -> distanceTargetGet(i_part);
       CWP_g_num_t* toto= ((*_geometry_cpl)._closest_elt_gnum)[i_part];
       double* tata= ((*_geometry_cpl)._distance)[i_part];
       for(int i=0; i< n_target_cpl; i++) {
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
       } 
       PDM_gnum_location_requested_elements_set(id_gnum_location,i_part, n_target_cpl, _geometry_cpl -> closestEltGnumGet  (i_part));
     }
  }
 }
 /*************************************************************************************/
 /*************************************************************************************/
  /*************************************************************************************/

  void GeomLocation::locate_get_cpl(int id_dist) {
    _geometry_cpl -> _distance           = (double**)malloc(sizeof(double*) * _nb_part);
    _geometry_cpl -> _projected          = (double**)malloc(sizeof(double*) * _nb_part);
    _geometry_cpl -> _closest_elt_gnum   = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*) * _nb_part);
    
    for(int i_part =0;i_part<_nb_part;i_part++) {     
        Mesh* mesh_cpl = _geometry_cpl -> meshGet();     
        int          n_target_cpl    = _geometry_cpl -> nTargetGet(i_part);
        CWP_g_num_t* tmp         ;//    = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t) * n_target_cpl);
        
        double     * tmp2        ;//    =(double*)     malloc( sizeof(double) * n_target_cpl);
        double     * tmp3        ;//    =(double*)     malloc(3* sizeof(double) * n_target_cpl);
        
        PDM_mesh_dist_get (id_dist,
                         0,
                         i_part,
                         &tmp2,
                         &tmp3,
                         &tmp);

       for(int i=0;i<n_target_cpl;i++){
         if(tmp[i]>CWP_g_num_t(_n_g_elt_over_part) || tmp[i]<CWP_g_num_t(1) || tmp2[i]>0.01){
         //  printf("CPL_distance[%i] %f gnum %I64d\n",i,tmp2[i],tmp[i]);
           tmp[i]=CWP_g_num_t(1);
           tmp2[i]=INFINITY;
         }       
       }         
       _geometry_cpl -> _distance        [i_part] = (double*)     malloc( sizeof(double) * n_target_cpl);
       _geometry_cpl -> _projected       [i_part] = (double*)     malloc(3 * sizeof(double) * n_target_cpl);
       _geometry_cpl -> _closest_elt_gnum[i_part] = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t) * n_target_cpl);
       
       memcpy((_geometry_cpl -> _distance        )[i_part],tmp2, sizeof(double)      * n_target_cpl );
       memcpy((_geometry_cpl -> _projected       )[i_part],tmp3, 3 * sizeof(double)      * n_target_cpl );
       memcpy((_geometry_cpl -> _closest_elt_gnum)[i_part],tmp,      sizeof(CWP_g_num_t) * n_target_cpl );
    }

 }// End locate_cell_point




  void GeomLocation::locate_get(int id_dist) {
  
    _distance           = (double**)malloc(sizeof(double*) * _nb_part);
    _projected          = (double**)malloc(sizeof(double*) * _nb_part);
    _closest_elt_gnum   = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*) * _nb_part);

    for(int i_part =0;i_part<_nb_part;i_part++) {     
    
        CWP_g_num_t* tmp ;
        double     * tmp2;
        double     * tmp3;
        PDM_mesh_dist_get (id_dist,
                         0,
                         i_part,
                         &tmp2,
                         &tmp3,
                         &tmp);
                         
       _distance [i_part] = (double*)malloc( sizeof(double) * _n_target[i_part]);
       _projected[i_part] = (double*)malloc(3 * sizeof(double) * _n_target[i_part]);
       _closest_elt_gnum[i_part] = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t) * _n_target[i_part]);


       for(int i=0;i<_n_target[i_part];i++){
         if(tmp[i]>CWP_g_num_t(_n_g_elt_cpl_over_part) || tmp[i]<CWP_g_num_t(1) || tmp2[i]>0.1){
           printf("_distance[%i] %f gnum %I64d\n",i,tmp2[i],tmp[i]);
           tmp[i]=CWP_g_num_t(1);
           tmp2[i]=INFINITY;
           
         }       
       }         
       
       memcpy(_distance        [i_part],tmp2,  sizeof(double)      * _n_target[i_part] );
       memcpy(_projected       [i_part],tmp3, 3 * sizeof(double)      * _n_target[i_part] );
       memcpy(_closest_elt_gnum[i_part],tmp,      sizeof(CWP_g_num_t) * _n_target[i_part] );
    }

 }// End locate_cell_point


 /*************************************************************************************/
 /*************************************************************************************/
  /*************************************************************************************/
 
   void GeomLocation::location_compute(int id_gnum_location) {
    PDM_gnum_location_compute(id_gnum_location);  
  }


 void GeomLocation::location_get(int id_gnum_location) {
   
  for (int i_part = 0; i_part < _nb_part; i_part++) {

    int* tmp;
    int* tmp2;

    PDM_gnum_location_get(id_gnum_location,
                            i_part,
                            &tmp,
                            &tmp2
                            );  


    for(int u =0;u<1+_n_target[i_part];u++) {     
      ( _location_idx)[i_part][u]=tmp[u];
    }
    
    for(int u =0;u<3*_n_target[i_part];u++) { 
    ( _location)[i_part][u]=tmp2[u];
    }

   // memcpy(_location_idx[i_part],tmp ,(1+_n_target[i_part]) * sizeof(int));
   // memcpy(_location    [i_part],tmp2,3* _n_target[i_part] * sizeof(int));
  
   } 
 } 



  void GeomLocation::location_get_cpl(int id_gnum_location) {
    for(int i_part =0;i_part<_nb_part;i_part++) {     
      int* tmp;
      int* tmp2;
      int          n_target_cpl    =  _geometry_cpl -> nTargetGet(i_part);
    
      PDM_gnum_location_get(id_gnum_location,
                            i_part,
                            &tmp,
                            &tmp2
                            );  
    
    
      for(int u =0;u<1+n_target_cpl;u++) {     
        (_geometry_cpl -> _location_idx)[i_part][u]=tmp[u];
      }
      //memcpy(_geometry_cpl -> _location_idx[i_part], tmp,   (1+n_target_cpl)*sizeof(int) );
      for(int u =0;u<3*n_target_cpl;u++) { 
        (_geometry_cpl -> _location)[i_part][u]=tmp2[u];
      }
      //memcpy(_geometry_cpl -> _location    [i_part], tmp2, 3*n_target_cpl   *sizeof(int) );
    }// end i_part loop

 }// End locate_cell_point

 void GeomLocation::broadcasting_filling_of_broadcasting_array() {


  _location_idx_comm_proc   =(int**)malloc(sizeof(int*)*_n_ranks_g);
  _location_count_comm_proc =(int**)malloc(sizeof(int*)*_n_ranks_g);

  for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
    _location_idx_comm_proc [i_proc]=(int*)malloc(sizeof(int)*(1+_nb_part));
    _location_count_comm_proc [i_proc]=(int*)malloc(sizeof(int)*(1+_nb_part));
    for (int i_part = 0; i_part < _nb_part+1; i_part++) {
      _location_idx_comm_proc [i_proc][i_part]=0;
      _location_count_comm_proc [i_proc][i_part]=0;
    }
  }


  for (int i_part = 0; i_part < _nb_part_cpl; i_part++) {
    for(int k=0;k<_n_target[i_part];k++){
       
  //     printf("_location_idx[i_part][%i] rank %i %i\n",k,_rank,_location_idx[i_part][k]);
  //     printf("_location[i_part][%i] rank %i %i\n",k,_rank,_location[i_part][ _location_idx[i_part][k] ]);
       int elt_proc = _location[i_part][ _location_idx[i_part][k] ];
       int elt_part = _location[i_part][ _location_idx[i_part][k] + 1];
      
       _location_idx_comm_proc [elt_proc][elt_part]++;     
       _location_count_comm_proc [elt_proc][elt_part]++;   
    }
  }//end i_part

  _transform_to_index(_location_idx_comm_proc,_n_ranks_g,_nb_part);

  int** idx_proc = (int**)malloc(sizeof(int*)* _n_ranks_g);
  for (int i_proc = 0; i_proc < _n_ranks_g; i_proc++) {
   idx_proc[i_proc] = (int*)malloc(sizeof(int)* _nb_part_cpl);
   for (int i_part = 0; i_part < _nb_part; i_part++)
      idx_proc[i_proc][i_part]=0;
   }


  _location_comm_proc = (target_data*)malloc(sizeof(target_data)*_location_idx_comm_proc[_n_ranks_g-1][_nb_part]);  

  for (int i_part = 0; i_part < _nb_part; i_part++) {
    for(int k=0;k<_n_target[i_part];k++){   
       int elt_proc = _location[i_part][ _location_idx[i_part][k]    ];
       int elt_part = _location[i_part][ _location_idx[i_part][k] + 1];
       int num = _location[i_part][ _location_idx[elt_part][k] + 2] - 1;
      
       int idx = _location_idx_comm_proc[elt_proc][elt_part] + idx_proc[elt_proc][elt_part];
       //Local Numbering
       _location_comm_proc [idx].lnum              = num;
       //Coupled numbering
       _location_comm_proc [idx ].origin_part      = i_part;
       _location_comm_proc [idx ].projectedX       = _projected [i_part][ 3*k    ];
       _location_comm_proc [idx ].projectedY       = _projected [i_part][ 3*k +1 ];
       _location_comm_proc [idx ].projectedZ       = _projected [i_part][ 3*k +2 ];       
       _location_comm_proc [idx ].distance         = _distance  [i_part][   k    ];    
       //Closest local element numbering
       _location_comm_proc [idx ].closest_elt_gnum  = _closest_elt_gnum[i_part][ k ];
       //Coupled process origin
       _location_comm_proc [idx].origin_proc       = _rank ;  
       //Coupled origin partition
       _location_comm_proc [idx ].closest_elt_part  = elt_part;  
       _location_comm_proc [idx ].l_num_origin      = k;  
       idx_proc[elt_proc][elt_part]++; 
 
 /*        printf("_location_comm_proc rank %i [%i] lnum %i origin_part %i closest_elt_gnum %i elt_proc %i elt_part %i distance %f proj %f\n",
       _rank,k,
       _location_comm_proc [idx ].lnum,
       _location_comm_proc [idx ].origin_part,
       int(_location_comm_proc [idx ].closest_elt_gnum),
       elt_proc,elt_part,
       _location_comm_proc [idx ].distance,
       _location_comm_proc [idx ].projectedX);*/
    }
  }//end i_part

  _location_idx_proc_recv = (int**)malloc(  sizeof(int*) *_n_ranks_g);      
  
  for(int i =0; i<_n_ranks_g;i++) {
    _location_idx_proc_recv[i] = (int*)malloc(sizeof(int)*(1+_nb_part)); 
    for (int i_part = 0; i_part < _nb_part+1; i_part++) {
      _location_idx_proc_recv[i][i_part]=0;
    }
  }
 } 


 void GeomLocation::broadcasting_index_communication() {
    _IAlltoall2(
      _location_count_comm_proc,
      NULL, 
      _nb_part_cpl,
      _location_idx_proc_recv,
      NULL,
      _nb_part,
      MPI_INT,
      _globalComm,
     *_connectableRanks_cpl
      ); 
 } 


 void GeomLocation::broadcasting_communication() {

  _transform_to_index(_location_idx_proc_recv,_n_ranks_g,_nb_part);
  
  _location_recv = (target_data*)malloc(sizeof(target_data)*_location_idx_proc_recv[_n_ranks_g-1][_nb_part_cpl]);        

  _location_count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);         
  _location_count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
  _location_disp_recv  = (int*)malloc(sizeof(int)*_n_ranks_g);         
  _location_disp_send  = (int*)malloc(sizeof(int)*_n_ranks_g); 
 
  for (int i= 0; i < _n_ranks_g; i++) { 
    _location_count_send[i] = sizeof(target_data)*(_location_idx_comm_proc[i][_nb_part_cpl] - _location_idx_comm_proc[i][0]);
    _location_count_recv[i] = sizeof(target_data)*(_location_idx_proc_recv[i][_nb_part_cpl] - _location_idx_proc_recv[i][0]);

    _location_disp_send [i] = sizeof(target_data)*_location_idx_comm_proc[i][0];
    _location_disp_recv [i] = sizeof(target_data)*_location_idx_proc_recv[i][0];
  }

 }
 
 void GeomLocation::broadcasting_communication2() {
  
  _IAlltoallIndex((void*)_location_comm_proc, _location_count_send, _location_disp_send,
                  (void*)_location_recv, _location_count_recv, _location_disp_recv, 
                  MPI_BYTE,
                  _globalComm, *_connectableRanks_cpl);
  }
  

 void GeomLocation::broadcasting_wait_and_targets_array_filling() {  
    _Wait();
  
   _targets_cpl       = _location_recv;
   
   //Continent l'index et donc le nombre de cible pour chaque proc et chaque partition local
   _targets_cpl_idx     = _location_idx_proc_recv;
      
   _targets_cpl_cpl     = _location_comm_proc;
   
   _targets_cpl_idx_cpl = _location_idx_comm_proc;

   _n_tot_target_cpl    = _location_idx_proc_recv[_n_ranks_g-1][_nb_part];

 }
 
/*******************************************************************************************************************/
/*******************************************************************************************************************/
/*******************************************************************************************************************/


  
}//end_cwipi
  

