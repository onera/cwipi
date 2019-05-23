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

  void _IAlltoallIndex2vtx(target_data_vtx* send_buffer,
                int** send_idx,
                int send_stride,
                target_data_vtx* recv_buffer,
                int** recv_idx,
                int recv_stride,
                MPI_Comm comm,
                std::vector<int> connectableRanks_cpl,
                int _n_ranks_cpl,
                int _nb_part_cpl){
      MPI_Status status;
      int tag = 0;
            
      std::vector<int> send_requests(_n_ranks_cpl,0);
      std::vector<int> recv_requests(_n_ranks_cpl,0);

      for(int i_rank=0;i_rank<_n_ranks_cpl;i_rank++) {
        int distant_rank = connectableRanks_cpl[i_rank];
        
        printf("distant_rank %i\n",distant_rank);

        int ind_proc      = recv_idx[i_rank][0];
        int ind_proc_send = send_idx[i_rank][0];
        
          target_data_vtx* ptr_send = ( target_data_vtx*)send_buffer;
          int send_size = (send_idx[i_rank][_nb_part_cpl] - send_idx[i_rank][0])*sizeof(target_data_vtx);
          MPI_Issend(&(  ptr_send [ind_proc_send] ), send_stride * send_size, MPI_BYTE, distant_rank, tag,
                   comm,
                   &(send_requests[i_rank]));        
   
          int recv_size = (recv_idx[i_rank][_nb_part_cpl] - recv_idx[i_rank][0])*sizeof(target_data_vtx);
          target_data_vtx* ptr = (target_data_vtx*)recv_buffer;
         printf("IND_PROC %i %i\n",send_size,recv_size);
          MPI_Irecv(&(  ptr [ind_proc] ), recv_stride * recv_size,MPI_BYTE, distant_rank, tag,
                  comm,
                  &(recv_requests[i_rank]));  

       
       
       
      }//end for on i_rank

      for(int i_rank=0;i_rank<_n_ranks_cpl;i_rank++) {
        MPI_Wait(&(send_requests[i_rank]), &status);
        MPI_Wait(&(recv_requests[i_rank]), &status);
      }//end for on i_rank 
  }



  void _IAlltoallIndex2(target_data* send_buffer,
                int** send_idx,
                int send_stride,
                target_data* recv_buffer,
                int** recv_idx,
                int recv_stride,
                MPI_Comm comm,
                std::vector<int> connectableRanks_cpl,
                int _n_ranks_cpl,
                int _nb_part_cpl){
      MPI_Status status;
      int tag = 0;
            
      std::vector<int> send_requests(_n_ranks_cpl,0);
      std::vector<int> recv_requests(_n_ranks_cpl,0);

      for(int i_rank=0;i_rank<_n_ranks_cpl;i_rank++) {
        int distant_rank = connectableRanks_cpl[i_rank];
        
        printf("distant_rank %i\n",distant_rank);

        int ind_proc      = recv_idx[i_rank][0];
        int ind_proc_send = send_idx[i_rank][0];
        
          target_data* ptr_send = ( target_data*)send_buffer;
          int send_size = (send_idx[i_rank][_nb_part_cpl] - send_idx[i_rank][0])*sizeof(target_data);
          MPI_Issend(&(  ptr_send [ind_proc_send] ), send_stride * send_size, MPI_BYTE, distant_rank, tag,
                   comm,
                   &(send_requests[i_rank]));        
   
          int recv_size = (recv_idx[i_rank][_nb_part_cpl] - recv_idx[i_rank][0])*sizeof(target_data);
          target_data* ptr = (target_data*)recv_buffer;
         printf("IND_PROC %i %i\n",send_size,recv_size);
          MPI_Irecv(&(  ptr [ind_proc] ), recv_stride * recv_size,MPI_BYTE, distant_rank, tag,
                  comm,
                  &(recv_requests[i_rank]));  

       
       
       
      }//end for on i_rank

      for(int i_rank=0;i_rank<_n_ranks_cpl;i_rank++) {
        MPI_Wait(&(send_requests[i_rank]), &status);
        MPI_Wait(&(recv_requests[i_rank]), &status);
      }//end for on i_rank 
  }


 void _IAlltoallIndex1vtx(target_data_vtx** send_buffer,
                int** send_idx,
                int send_stride,
                target_data_vtx* recv_buffer,
                int** recv_idx,
                int recv_stride,
                MPI_Comm comm,
                std::vector<int> connectableRanks_cpl,
                int _n_ranks_cpl,
                int _nb_part){
      MPI_Status status;
      int tag = 0;
      int ind_proc = 0;
      int ind_proc_send = 0;
            
      std::vector<int> send_requests(_n_ranks_cpl,0);
      std::vector<int> recv_requests(_n_ranks_cpl,0);

      for(int i_rank=0;i_rank<_n_ranks_cpl;i_rank++) {
        int distant_rank = connectableRanks_cpl[i_rank];
        
        printf("distant_rank %i\n",distant_rank);


        

          target_data_vtx** ptr_send = ( target_data_vtx**)send_buffer;
          int send_size = (send_idx[i_rank][_nb_part] - send_idx[i_rank][0])*sizeof(target_data_vtx);
          MPI_Issend(&(  ptr_send [i_rank] [0] ), send_stride * send_size, MPI_BYTE, distant_rank, tag,
                   comm,
                   &(send_requests[i_rank]));        
   
          int recv_size = (recv_idx[i_rank][_nb_part] - recv_idx[i_rank][0])*sizeof(target_data_vtx);
          ind_proc = recv_idx[i_rank][0];
          
          target_data_vtx* ptr = (target_data_vtx*)recv_buffer;
          printf("IND_PROC %i %i\n",send_size,recv_size);
          MPI_Irecv(&(  ptr [ind_proc] ), recv_stride * recv_size,MPI_BYTE, distant_rank, tag,
                  comm,
                  &(recv_requests[i_rank]));  

       ind_proc_send = send_stride * send_idx[i_rank][0];
      }//end for on i_rank

      for(int i_rank=0;i_rank<_n_ranks_cpl;i_rank++) {
        MPI_Wait(&(send_requests[i_rank]), &status);
        MPI_Wait(&(recv_requests[i_rank]), &status);
      }//end for on i_rank 
  }



  void _IAlltoallIndex1(target_data** send_buffer,
                int** send_idx,
                int send_stride,
                target_data* recv_buffer,
                int** recv_idx,
                int recv_stride,
                MPI_Comm comm,
                std::vector<int> connectableRanks_cpl,
                int _n_ranks_cpl,
                int _nb_part){
      MPI_Status status;
      int tag = 0;
      int ind_proc = 0;
      int ind_proc_send = 0;
            
      std::vector<int> send_requests(_n_ranks_cpl,0);
      std::vector<int> recv_requests(_n_ranks_cpl,0);

      for(int i_rank=0;i_rank<_n_ranks_cpl;i_rank++) {
        int distant_rank = connectableRanks_cpl[i_rank];
        
        printf("distant_rank %i\n",distant_rank);


        

          target_data** ptr_send = ( target_data**)send_buffer;
          int send_size = (send_idx[distant_rank][_nb_part] - send_idx[distant_rank][0])*sizeof(target_data);
          MPI_Issend(&(  ptr_send [distant_rank] [0] ), send_stride * send_size, MPI_BYTE, distant_rank, tag,
                   comm,
                   &(send_requests[i_rank]));        
   
          int recv_size = (recv_idx[distant_rank][_nb_part] - recv_idx[distant_rank][0])*sizeof(target_data);
          ind_proc = recv_idx[distant_rank][0];
          
          target_data* ptr = (target_data*)recv_buffer;
          printf("IND_PROC %i %i\n",send_size,recv_size);
          MPI_Irecv(&(  ptr [ind_proc] ), recv_stride * recv_size,MPI_BYTE, distant_rank, tag,
                  comm,
                  &(recv_requests[i_rank]));  

       ind_proc_send = send_stride * send_idx[distant_rank][0];
      }//end for on i_rank

      for(int i_rank=0;i_rank<_n_ranks_cpl;i_rank++) {
        MPI_Wait(&(send_requests[i_rank]), &status);
        MPI_Wait(&(recv_requests[i_rank]), &status);
      }//end for on i_rank 
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

          
          for (int itarget = _targets_cpl_idx[i_proc][i_part]; itarget < _targets_cpl_idx[i_proc][i_part+1]; itarget++) {

            //Index of the corresponding local reference Data.
            int iel = _targets_vtx_cpl[itarget].lnum ;
            
            // Index in the interpolated Data array
            int interpInd = itarget;
            
            double x_target = _targets_vtx_cpl[itarget].projectedX;
            double y_target = _targets_vtx_cpl[itarget].projectedY;
            double z_target = _targets_vtx_cpl[itarget].projectedZ;    
                              
            double sum_vtx = 0.0;
             printf("iel %i itarget %i refData %f _n_tot_target_cpl %i\n",iel,itarget,referenceData[iel],_n_tot_target_cpl);

              interpolatedData[ interpInd  ] = 0.0;
              double dvtxmin = 100000;
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
               
             /*    printf("dvtx %f %i\n",dvtx,iel);
                 if(dvtx < dvtxmin) {
                    interpolatedData[ interpInd  ] += referenceData[connec[i_vtx]-1] ;
                    dvtxmin=dvtx;
                 }*/
               }
               interpolatedData[ interpInd  ] /= sum_vtx;
            
            
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

      int tag=0;
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
  
  
 
  void GeomLocation::locate_cell_point_setting_surface(int id_dist) {

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
        
        PDM_mesh_dist_cloud_set (id_dist,
                              0,
                              i_part,
                              0,
                              tmp6,
                              tmp3
                             );                              
      }
      else {
        Mesh* mesh_cpl = _geometry_cpl_cell_point -> meshGet();
        int          n_target_cpl      = _geometry_cpl_cell_point -> nTargetGet(i_part);
        CWP_g_num_t* gnum_elt_cpl      = mesh_cpl -> GNumEltsGet(i_part);     
        double*      centers_cpl       = mesh_cpl -> eltCentersGet(i_part);

        PDM_mesh_dist_cloud_set (id_dist,
                              0,
                              i_part,
                              n_target_cpl,
                              centers_cpl,
                              gnum_elt_cpl
                             );   
     }   
   }
  
  }

  void GeomLocation::locate_cell_point_setting_request(int id_dist) {

    /*
     TODO: Intéressant pour la suite
     Inclure dans mesh_dist et autre la possibilité de définir plusieurs surfaces ...

   */

    /* Paradigm mesh localisation _distance creation */
    // Envoyeur _codeVector[1] Surface _codeVector[0]
    id_dist   = PDM_mesh_dist_create( PDM_MESH_NATURE_SURFACE_MESH, 1, _pdm_globalComm );

    PDM_mesh_dist_n_part_cloud_set(id_dist,   0, _nb_part);  

    PDM_mesh_dist_surf_mesh_global_data_set (id_dist,
                                           _n_g_elt_cpl_over_part,
                                           _n_g_vtx_cpl_over_part,
                                           _nb_part_cpl);  

    for(int i_part =0;i_part<_nb_part;i_part++) {   

      CWP_g_num_t* gnum_elt = _mesh -> GNumEltsGet(i_part);
      double* centers       = _mesh -> eltCentersGet(i_part);

      PDM_mesh_dist_cloud_set (id_dist,
                              0,
                              i_part,
                              _n_target[i_part],
                              centers,
                              gnum_elt
                             );
    }
 
 
    for(int i_part =0; i_part<_nb_part_cpl; i_part++) {     
    
      if(_both_codes_are_local == 0) {
 /*       int* tmp1 = (int*) malloc(sizeof(int)*1000);
        int* tmp2 = (int*) malloc(sizeof(int)*1000);
         
        CWP_g_num_t* tmp3 = (CWP_g_num_t*) malloc(sizeof(CWP_g_num_t)*1000);
        CWP_g_num_t* tmp4 = (CWP_g_num_t*) malloc(sizeof(CWP_g_num_t)*1000);
        double* tmp6 = (double*) malloc(sizeof(double)*1000);
*/

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
        Mesh* mesh_cpl = _geometry_cpl_cell_point -> meshGet();
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

  void GeomLocation::locate_cell_point_compute(int id_dist) {
    PDM_mesh_dist_compute(id_dist);
  }

 void GeomLocation::redistribution_cell_point_request(int id_gnum_location) {
 
  id_gnum_location = PDM_gnum_location_create(_nb_part,_nb_part_cpl, _pdm_globalComm);

  for(int i_part =0;i_part<_nb_part;i_part++) {    
    CWP_g_num_t* gnum_elt = _mesh -> GNumEltsGet(i_part);
    PDM_gnum_location_requested_elements_set(id_gnum_location,i_part, _n_target[i_part],_closest_elt_gnum[i_part]);
  }

 
  for(int i_part =0; i_part<_nb_part_cpl; i_part++) {     
    
      if(_both_codes_are_local == 0) {
        CWP_g_num_t* gnum_elt = _mesh -> GNumEltsGet(i_part);
        PDM_gnum_location_elements_set(id_gnum_location,i_part, 0,gnum_elt);    
      }     
     else {
       Mesh* mesh_cpl = _geometry_cpl_cell_point -> meshGet();     
       CWP_g_num_t* gnum_elt_cpl = mesh_cpl -> GNumEltsGet(i_part);
       int          n_elts_cpl    = mesh_cpl -> getPartNElts(i_part);
       
       PDM_gnum_location_elements_set(id_gnum_location,i_part, n_elts_cpl,gnum_elt_cpl);    
     }
  }
 }


 void GeomLocation::redistribution_cell_point_set(int id_gnum_location) {
 
  id_gnum_location = PDM_gnum_location_create(_nb_part,_nb_part_cpl, _pdm_globalComm);

  for(int i_part =0;i_part<_nb_part_cpl;i_part++) {    

    CWP_g_num_t* gnum_elt = _mesh -> GNumEltsGet(i_part);
    
    PDM_gnum_location_elements_set(id_gnum_location,i_part, _n_target[i_part],gnum_elt);      
  }

  for(int i_part =0; i_part<_nb_part_cpl; i_part++) {     
    
      if(_both_codes_are_local == 0) {
      
        CWP_g_num_t* gnum_elt = _mesh -> GNumEltsGet(i_part);
        PDM_gnum_location_requested_elements_set(id_gnum_location,i_part, 0,gnum_elt);
      }     
     else {
       Mesh* mesh_cpl = _geometry_cpl_cell_point -> meshGet();     
       
       int          n_elts_cpl    = mesh_cpl -> getPartNElts(i_part);

       PDM_gnum_location_requested_elements_set(id_gnum_location,i_part, n_elts_cpl, (_geometry_cpl_cell_point -> _closest_elt_gnum)[i_part]);
     }
  }
 }
 /*************************************************************************************/
 /*************************************************************************************/
  /*************************************************************************************/

  void GeomLocation::locate_cell_point_get_cpl(int id_dist) {
    _geometry_cpl_cell_point -> _distance           = (double**)malloc(sizeof(double*) * _nb_part);
    _geometry_cpl_cell_point -> _projected          = (double**)malloc(sizeof(double*) * _nb_part);
    _geometry_cpl_cell_point -> _closest_elt_gnum   = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*) * _nb_part);
    
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
        
       Mesh* mesh_cpl = _geometry_cpl_cell_point -> meshGet();     
       int          n_target_cpl    = _geometry_cpl_cell_point -> nTargetGet(i_part);
                         
       _geometry_cpl_cell_point -> _distance [i_part] = (double*)malloc(3 * sizeof(double) * n_target_cpl);
       _geometry_cpl_cell_point -> _projected[i_part] = (double*)malloc(3 * sizeof(double) * n_target_cpl);
       _geometry_cpl_cell_point -> _closest_elt_gnum[i_part] = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t) * n_target_cpl);
       
       memcpy(_geometry_cpl_cell_point -> _distance        [i_part],tmp2, 3 * sizeof(double)      * n_target_cpl );
       memcpy(_geometry_cpl_cell_point -> _projected       [i_part],tmp3, 3 * sizeof(double)      * n_target_cpl );
       memcpy(_geometry_cpl_cell_point -> _closest_elt_gnum[i_part],tmp,      sizeof(CWP_g_num_t) * n_target_cpl );
    }

 }// End locate_cell_point




  void GeomLocation::locate_cell_point_get(int id_dist) {
  
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
                         
       _distance [i_part] = (double*)malloc(3 * sizeof(double) * _n_target[i_part]);
       _projected[i_part] = (double*)malloc(3 * sizeof(double) * _n_target[i_part]);
       _closest_elt_gnum[i_part] = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t) * _n_target[i_part]);
       
       memcpy(_distance        [i_part],tmp2, 3 * sizeof(double)      * _n_target[i_part] );
       memcpy(_projected       [i_part],tmp3, 3 * sizeof(double)      * _n_target[i_part] );
       memcpy(_closest_elt_gnum[i_part],tmp,      sizeof(CWP_g_num_t) * _n_target[i_part] );
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
 
  _location_idx =(int**)malloc(sizeof(int*)*_nb_part);
  _location     =(int**)malloc(sizeof(int*)*_nb_part);

  for (int i_part = 0; i_part < _nb_part; i_part++) {
    _location_idx[i_part]=NULL;
    _location    [i_part]=NULL;


    int* tmp;
    int* tmp2;
    int   n_elts    = _mesh -> getPartNElts(i_part);
    
    PDM_gnum_location_get(id_gnum_location,
                            i_part,
                            &tmp,
                            &tmp2
                            );  
 

    _location_idx[i_part] = (int*) malloc(sizeof(int) * (1+n_elts));
    _location[i_part]     = (int*) malloc(3* sizeof(int) * (n_elts));
    
    memcpy(_location_idx[i_part],tmp,(1+n_elts)*sizeof(int));
    memcpy(_location[i_part],tmp2,3*(n_elts)*sizeof(int));
  
   } 
 } 



  void GeomLocation::location_get_cpl(int id_gnum_location) {
   
    _geometry_cpl_cell_point -> _location_idx =(int**)malloc(sizeof(int*)*_nb_part);
    _geometry_cpl_cell_point -> _location     =(int**)malloc(sizeof(int*)*_nb_part);
    Mesh* mesh_cpl = _geometry_cpl_cell_point -> meshGet();
    
    for(int i_part =0;i_part<_nb_part;i_part++) {     
       int* tmp;
       int* tmp2;
            
       int          n_elts_cpl    = mesh_cpl -> getPartNElts(i_part);
    
    
      _geometry_cpl_cell_point -> _location_idx[i_part]=NULL;
      _geometry_cpl_cell_point -> _location    [i_part]=NULL;
      PDM_gnum_location_get(id_gnum_location,
                            i_part,
                            &tmp,
                            &tmp2
                            );  
 
    _geometry_cpl_cell_point ->_location_idx[i_part] = (int*) malloc(sizeof(int) * (1+n_elts_cpl));
    _geometry_cpl_cell_point -> _location[i_part] = (int*) malloc(3* sizeof(int) * n_elts_cpl);
    
    memcpy(_geometry_cpl_cell_point -> _location_idx[i_part], tmp,   (1+n_elts_cpl)*sizeof(int) );
    memcpy(_geometry_cpl_cell_point -> _location    [i_part], tmp2,3*n_elts_cpl*sizeof(int) );
      
    }

 }// End locate_cell_point

 void GeomLocation::redistribution_cell_point_filling_of_redistribution_array() {


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
    int   n_elts    = _mesh -> getPartNElts(i_part);
    for(int k=0;k<n_elts;k++){
       
      // printf("_location_idx[i_part][%i] rank %i %i\n",k,_rank,_location_idx[i_part][k]);
      // printf("_location[i_part][%i] rank %i %i\n",k,_rank,_location[i_part][ _location_idx[i_part][k] ]);
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
       //Closest local element numbering
       _location_comm_proc [idx ].closest_elt_gnum  = _closest_elt_gnum[i_part][ k ];
       //Coupled process origin
       _location_comm_proc [idx].origin_proc       = _rank ;  
       //Coupled origin partition
       _location_comm_proc [idx ].closest_elt_part  = elt_part;  
       _location_comm_proc [idx ].l_num_origin      = k;  
       idx_proc[elt_proc][elt_part]++; 
      // printf("_location_comm_proc [%i].closest_elt_gnum %i\n",int(_location_comm_proc [idx ].closest_elt_gnum));
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


 void GeomLocation::redistribution_index_communication(int tag) {
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


 void GeomLocation::redistribution_communication() {

  _transform_to_index(_location_idx_proc_recv,_n_ranks_g,_nb_part);

  // printf("_location_idx_proc_recv[_n_ranks_g-1][_nb_part_cpl] rank %i %i\n",_rank,_location_idx_proc_recv[_n_ranks_g-1][_nb_part_cpl]);

  _location_recv = (target_data*)malloc(sizeof(target_data)*_location_idx_proc_recv[_n_ranks_g-1][_nb_part_cpl]);        

  _location_count_recv = (int*)malloc(sizeof(int)*_n_ranks_g);         
  _location_count_send = (int*)malloc(sizeof(int)*_n_ranks_g);
  _location_disp_recv  = (int*)malloc(sizeof(int)*_n_ranks_g);         
  _location_disp_send  = (int*)malloc(sizeof(int)*_n_ranks_g); 
 
  for (int i= 0; i < _n_ranks_g; i++) { 
    _location_count_send[i] = sizeof(target_data)*(_location_idx_comm_proc[i][_nb_part_cpl] - _location_idx_comm_proc[i][0]);
    _location_count_recv[i] = sizeof(target_data)*(_location_idx_proc_recv[i][_nb_part_cpl] - _location_idx_proc_recv[i][0]);
    
    //printf("location_count rank %i i %i send %i recv %i size %i\n",
    //_rank,i,_location_count_send[i]/sizeof(target_data),_location_count_recv[i]/sizeof(target_data),_location_idx_proc_recv[_n_ranks_g-1][_nb_part_cpl]);
    _location_disp_send [i] = sizeof(target_data)*_location_idx_comm_proc[i][0];
    _location_disp_recv [i] = sizeof(target_data)*_location_idx_proc_recv[i][0];
  }

 }
 
 void GeomLocation::redistribution_communication2() {
  
  _IAlltoallIndex((void*)_location_comm_proc, _location_count_send, _location_disp_send,
                  (void*)_location_recv, _location_count_recv, _location_disp_recv, 
                  MPI_BYTE,
                  _globalComm, *_connectableRanks_cpl);
  }
  

 void GeomLocation::redistribution_wait_and_targets_array_filling() {  
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
  

