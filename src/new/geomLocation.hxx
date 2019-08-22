#ifndef __GEOMLOCATION_H__
#define __GEOMLOCATION_H__
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

#include "mesh.hxx"
#include "geometry.hxx"
#include "field.hxx"

namespace cwipi {

  class GeomLocation: public Geometry 
    {
  
    public:
      GeomLocation();
      
      virtual ~GeomLocation();

      private:     
      void localization_points_cloud_setting(int* id_dist);
      void localization_null_setting   (int* id_dist);          
      void localization_surface_setting(int* id_dist);     
      void locate_setting_code_1_local() ;
      
      void localization_compute(int id_dist)             ;
      void localization_get(int id_dist)      ;
      void localization_get_cpl(int id_dist)      ;
      
      void broadcasting_request (int* id_gnum_location) ;
      void broadcasting_set     (int* id_gnum_location) ;
      void broadcasting_set_null(int* id_gnum_location) ;      
      void location_compute                 (int id_gnum_location) ;
                  
      void location_get(int id_gnum_location)      ;
      void location_get_cpl(int id_gnum_location)  ;

        void filling_of_broadcasting_array ();
        void init(Coupling *coupling, CWP_Field_value_t geometryLocation,int slave) ;     
        void mesh_info_get();
        void mesh_cpl_info_get();  
        void info_mesh(CWP_Field_exch_t _Texch_t) ;  
      void initialization_of_reception_array ();
      void broadcasting_index_communication()    ;
      void broadcasting_index_communication_async()    ;      
      void reception_index_communication()    ;
      void both_index_communication()    ;
      void reception_index_communication_async()    ;
      
      void broadcasting_index_null();

      void prepare_data_communication_send()  ;
      void prepare_data_communication_recv()  ;

      void data_communication_send()          ;
      void data_communication_recv()          ;
      void data_communication_null()          ;
      void both_data_communication()          ;
      
      void data_communication_wait_send();
      void data_communication_wait_recv();    
           void computeFree();
           
      public:


      void compute(CWP_Field_exch_t Texch_t) ;
                 
      void issend(Field* sendingField);

    void 
    waitIssend
    (Field* sendingField);


    void 
    irecv
    (Field* recevingField
    );


    void 
    waitIrecv
    (Field* recevingField);    
   

      void exchange_null() ;
      
      void both_exchange(Field* sendingField,
                                 Field* recevingField
                                 ) ;
      void* interpolate (Field* referenceField);  

    private:
       GeomLocation              *_geometry_cpl          ;    
   

     
  protected:

     //Pointer to other objects
    CWP_Field_value_t                    _geometryLocation      ;

    /* Localization data */
    double      **_distance         ; 
    double      **_projected        ; 
    CWP_g_num_t **_closest_elt_gnum ;  

    double      **_distance_cpl         ; 
    double      **_projected_cpl        ; 
    CWP_g_num_t **_closest_elt_gnum_cpl ;  

    int         **_targets_localization_idx     ;
    int         **_targets_localization_idx_cpl ;
    target_data  *_targets_localization_data    ;
    target_data  *_targets_localization_data_cpl;  
    
    //TODO: To delete and replace by using other members    
    std::vector<int>   _idx_target        ;    /*!< Index of the number of target by partition */

    /* Displacement and count for all_to_all MPI communication of targets_localization_data */
    int* _targets_localization_data_count_recv;         
    int* _targets_localization_data_count_send ;
    int* _targets_localization_data_disp_recv ;         
    int* _targets_localization_data_disp_send ; 

    /* Triplet global numbering, MPI process, partition results */
    int** _process_and_partition_count ; /*!< Element count by MPI process rank and partition */

    int** _target_proc_part_num_idx ; /*!< Index array of triplet process, partition, numbering for each target */
    int** _target_proc_part_num     ; /*!< Array of triplet process, partition, numbering for each target */

   /* Mesh informations */

   CWP_g_num_t _n_g_elt_over_part    ; /*!< Number of element of the process (over all the partitions)              */
   CWP_g_num_t _n_g_vtx_over_part    ; /*!< Number of vertices of the process (over all the partitions)             */
   CWP_g_num_t _n_g_elt_cpl_over_part; /*!< Number of coupled code element of the process (over all the partitions) */
   CWP_g_num_t _n_g_vtx_cpl_over_part; /*!< Number of coupled code vertices of the process (over all the partitions)*/

   int  _n_tot_target    ; /*!< Target total number on the process                                      */
   int  _n_tot_target_cpl; /*!< Number of coupled code target received by the process for interpolation */    
   int *_n_target        ; /*!< Target total number on the process by partition                         */
        
   int *_n_vtx           ; /*!< Vertice total number on the process by partition                        */
   int  _n_tot_vtx       ; /*!< Vertice total number on the process                                     */
      
   int *_n_elt           ; /*!< Element total number on the process by partition                        */
   int  _n_tot_elt       ; /*!< Element total number on the process                                     */
 
   int  _nb_part_cpl     ; /*!< Coupled code mesh partition number                                      */
   int  _nb_part         ; /*!< Mesh partition number                                                   */
 
   CWP_g_num_t  **_gnum_target  ;  /*<! Target global numbering by partition */
   double       **_coords_target;  /*<! Target coordinates by partition */
  
    /* Paradigm structure identifier */
    int _id_dist         ; /*!< Identifier for the localization object of paradigm */
    int _id_gnum_location; /*!< Identifier for the global numbering to (process,partition,numbering) triplet object of paradigm */

  }; //end GeomLocation
  

  
}
#endif //__GEOMLOCATION_H__


