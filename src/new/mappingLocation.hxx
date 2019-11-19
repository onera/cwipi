#ifndef __MAPPINGLOCATION_H__
#define __MAPPINGLOCATION_H__
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
#include "mapping.hxx"
#include "field.hxx"

namespace cwipi {

  class MappingLocation: public Mapping 
    {
  
    public:


    /**
      *
      * \brief Mapping location constructor.
      *
      */

      MappingLocation();


    /**
      *
      * \brief Mapping location destructor.
      *
      */
      
      virtual ~MappingLocation();


    /**
      *
      * \brief Compute of the coupling mapping. Localization and communication
      *        tree building.
      *
      * \param [in] Texch_t    Type of exchange for the mapping (sending or reception).
      *
      */

      void compute(CWP_Field_exch_t Texch_t) ;


    /**
      *
      * \brief Non-blocking sending communication including interpolation.
      *
      * \param [in] sendingField    Pointer to the referenceField used for interpolation.
      *
      */

      void issend(Field* referenceField);
      void issend_p2p(Field* referenceField);
    /**
      *
      * \brief Wait for non-blocking sending communication.
      *
      * \param [in] sendingField    Pointer to the referenceField used for interpolation.
      *
      */

      void waitIssend(Field* referenceField);
      void waitIssend_p2p(Field* referenceField);


    /**
      *
      * \brief Non-blocking receving communication.
      *
      * \param [in] recevingField    Pointer to the resulting interpolated field.
      *
      */

      void irecv(Field* recevingField);
      void irecv_p2p(Field* recevingField);
    /**
      *
      * \brief Wait for non-blocking receving communication.
      *
      * \param [in] recevingField     Pointer to the resulting interpolated field.
      *
      */

      void waitIrecv (Field* recevingField);    
      void waitIrecv_p2p (Field* recevingField);  

    /**
      *
      * \brief Null exchange communication for uncoupled MPI process.
      *
      */ 
   
      void null_exchange_for_uncoupled_process() ;


    /**
      *
      * \brief Exchange communication in case where local and coupled codes are
      *        on the same MPI process.
      *
      * \param [in] sendingField      Pointer to the referenceField used for interpolation.
      * \param [in] recevingField     Pointer to the resulting interpolated field.
      *
      */

      void both_codes_on_the_same_process_exchange (Field* sendingField,
                                                    Field* recevingField
                                                    ) ;

    private:     


    /**
      *
      * \brief Initialization of the Mapping object.
      *
      * \param [in] coupling            Pointer the coupling object.
      * \param [in] pointsCloudLocation Location of the cloud of points.
      * \param [in] coupling            Pointer the coupling object.
      *
      */
    
      void init (Coupling *coupling, CWP_Field_value_t pointsCloudLocation,int slave) ;   
  

      /***********************************************************
       **           Mesh information functions                  **
       ***********************************************************/

    /**
      *
      * \brief Get informations from the code mesh to use 
      * in Mapping object.
      *
      */
        
      void mesh_info_get();

    /**
      *
      * \brief Get informations from the coupled code mesh to use
      * in Mapping object.
      *
      */

      void mesh_cpl_info_get(); 

    /**
      *
      * \brief Get informations from local and coupled code mesh 
      *  to use in Mapping object.
      *
      */
 
      void info_mesh() ;  


      /***********************************************************
       ***********************************************************
       **                                                       **
       **            Localization object functions              **
       **                                                       **
       ***********************************************************
       ***********************************************************/
 

   /**
      *
      * \brief Setting of the points cloud for localization.
      *
      * \param [out] id_dist   Localization object identifier.
      *
      */
     
      void localization_points_cloud_setting (int* id_dist) ;


    /**
      *
      * \brief Setting of the surface mesh and cloud points at 
      * null for the localization object in a case of sending 
      * code i.e. code which interpolate reference field. 
      *
      * \param [out] id_dist   Localization object identifier.
      *
      */

      void localization_null_setting_send (int* id_dist) ;          


    /**
      *
      * \brief Setting of the surface mesh and cloud points at 
      * null for the localization object in a case of receving 
      * code i.e. code which provides cloud points for interpolation.
      *
      * \param [out] id_dist   Localization object identifier.
      *
      */
  
      void localization_null_setting_recv (int* id_dist) ;     


    /**
      *
      * \brief Setting of the surface mesh and cloud points at 
      * null for the localization object in a case of receving 
      * code i.e. code which provides cloud points for interpolation.
      *
      * \param [out] id_dist   Localization object identifier.
      *
      */     
      
      void localization_surface_setting (int* id_dist) ;     
      

    /**
      *
      * \brief Compute of localization of a points cloud on a surface
      *  mesh through the localization object.
      *
      * \param [int] id_dist   Localization object identifier.
      *
      */   

      void localization_compute (int id_dist) ;


    /**
      *
      * \brief Get localization results from localization object.
      *
      * \param [int] id_dist   Localization object identifier.
      *
      */ 

      void localization_get (int id_dist) ;


    /**
      *
      * \brief Get localization results from localization object 
      * from the coupled code in the case where the both codes are on 
      * the same process. 
      *
      * \param [int] id_dist   Localization object identifier.
      *
      */ 

      void localization_get_cpl (int id_dist) ;


      /***********************************************************
       ***********************************************************
       **                                                       **
       **   Process, partition, num triplet location from       **
       **           global numbering functions                  **
       **                                                       **
       ***********************************************************
       ***********************************************************/
 

   /**
      *
      * \brief Setting of requested global numbering for process, partition,
      *        num triplet location from global numbering object.
      *        
      * \param [in] id_gnum_location  process, partition, num triplet location 
      *              from global numbering identifier.
      *
      */
      
      void triplet_location_request (int* id_gnum_location) ;


   /**
      *
      * \brief Setting of researched global numbering for process, partition,
      *        num triplet location from global numbering object.
      *        
      * \param [in] id_gnum_location  rocess, partition, num triplet location 
      *              from global numbering identifier.
      *
      */

      void triplet_location_set (int* id_gnum_location) ;


    /**
      *
      * \brief Setting of researched global numbering for process, partition,
      *  num triplet location from global numbering object in a case of sending 
      * code i.e. code which interpolates provided cloud points.
      *
      * \param [in] id_gnum_location  rocess, partition, num triplet location 
      *              from global numbering identifier.
      *
      */

      void triplet_location_null_send (int* id_gnum_location) ;   

    /**
      *
      * \brief Setting of researched global numbering for process, partition,
      *  num triplet location from global numbering object in a case of receving 
      * code i.e. code which provides cloud points for interpolation.
      *
      * \param [in] id_gnum_location  rocess, partition, num triplet location 
      *              from global numbering identifier.
      *
      */
   
      void triplet_location_null_recv (int* id_gnum_location) ;  


    /**
      *
      * \brief Compute of process, partition, num triplet location 
      *        from global numbering object.
      *
      * \param [in] id_gnum_location    Process, partition, num triplet location 
      *                                 from global numbering identifier.
      *
      */
           
      void triplet_location_compute  (int id_gnum_location) ;



    /**
      *
      * \brief Get process, partition, num triplet location
      *        the case where the both codes are on 
      *        the same process. 
      *
      * \param [in] id_gnum_location     Process, partition, num triplet location 
      *                                  from global numbering identifier.
      *
      */ 
                  
      void triplet_location_get(int id_gnum_location)      ;


    /**
      *
      * \brief Get process, partition, num triplet location
      * from the coupled code in the case where the both codes are on 
      * the same process. 
      *
      * \param [in] id_gnum_location     Process, partition, num triplet location 
      *                                  from global numbering identifier.
      */ 

      void triplet_location_get_cpl(int id_gnum_location)  ;


      /***********************************************************
       ***********************************************************
       **            Communication tree array functions         **
       **                                                       **
       ***********************************************************
       ***********************************************************/

    /**
      *
      * \brief Initialization of the communication tree array
      *        containing localization informations of the coupled
      *        mesh point cloud.
      *
      */ 

      void initialization_of_receving_communication_tree_array ();


    /**
      *
      * \brief Filling of the communication tree array
      *        containing localization informations of the 
      *        mesh point cloud.
      *
      */ 

      void filling_of_sending_communication_tree_array ();



      /***********************************************************
       ***********************************************************
       **            Data index communication functions         **
       **                                                       **
       ***********************************************************
       ***********************************************************/

    /**
      *
      * \brief Send of the communication tree array index
      *        containing localization informations of the 
      *        mesh point cloud.
      *
      */ 

      void data_index_communication_send()    ;


      void data_index_communication_send_p2p()    ;



    /**
      *
      * \brief Reception of the communication tree array index
      *        containing localization informations of the 
      *        coupled mesh point cloud.
      *
      */ 

      void data_index_communication_recv()    ;

      void data_index_communication_recv_p2p()    ;

    /**
      *
      * \brief Send and reception of the communication tree 
      *        array index containing localization informations 
      *        in a case where the both are on the same MPI process. 
      *
      */ 

      void both_index_communication()    ;
      void both_index_communication_p2p()    ;


    /**
      *
      * \brief Null communication the communication tree 
      *        array index for uncoupled MPI process. 
      *
      */ 

      void data_index_communication_null();

      /***********************************************************
       ***********************************************************
       **            Data communication functions               **
       **                                                       **
       ***********************************************************
       ***********************************************************/

      void prepare_data_communication_send()  ;
      void prepare_data_communication_recv()  ;

      void data_communication_send()          ;
      void data_communication_recv()          ;
      void data_communication_null()          ;
      void both_data_communication()          ;
 
      void data_communication_send_p2p()      ;
      void data_communication_recv_p2p()      ;
      void both_data_communication_p2p()      ;
      
      void data_communication_wait_send()     ;
      void data_communication_wait_recv()     ;    
        
      void computeFree();

      /***********************************************************
       **         User definde cloud points functions           **
       ***********************************************************/

      void user_target_points_set(int i_part, int n_pts, double* coord);
      void user_targets_gnum_compute();

    /**
      *
      * \brief Interpolation of a point cloud on a reference field.
      *
      * \param [in]   referenceField   Reference field pointer
      *
      */
            
      void* interpolate (Field* referenceField); 

        
      MappingLocation    *_mapping_cpl            ;  /*!< Coupled code mapping object (for both codes are local case) */
   
      CWP_Field_value_t    _pointsCloudLocation   ;  /*!< Type of points cloud treated by this mapping instance (cell centers, vertices or user defined) */
 
       /* Localization data */
 
      double      **_distance                     ;  /*!< Distance to the closest element surface by partition */
      double      **_projected                    ;  /*!< Projected point coordinates (on the closest element surface) */
      CWP_g_num_t **_closest_elt_gnum             ;  /*!< Closest element global numbering */

      int         **_targets_localization_idx     ;  /*!< Data index (by process and by partition) of target localization*/
      target_data  *_targets_localization_data    ;  /*!< Data of target localization */      
      int         **_targets_localization_idx_cpl ;  /*!< Data index (by process and by partition) of the received target localization*/
      target_data  *_targets_localization_data_cpl;  /*!< Data of the received target localization */
    
      //TODO: To delete and replace by using other members    
      std::vector<int>   _idx_target              ;  /*!< Index of the number of target by partition */

      /* Displacement and count for all_to_all MPI communication of targets_localization_data */

      int* _targets_localization_data_count_recv  ;  /* Counts for all_to_all MPI communication of targets_localization_data (reception) */       
      int* _targets_localization_data_count_send  ;  /* Counts for all_to_all MPI communication of targets_localization_data (sending) */   
      int* _targets_localization_data_disp_recv   ;  /* Displacements for all_to_all MPI communication of targets_localization_data (reception) */         
      int* _targets_localization_data_disp_send   ;  /* Displacements for all_to_all MPI communication of targets_localization_data (sending) */   

      /* Triplet global numbering, MPI process, partition results */

      int** _process_and_partition_count          ;  /*!< Element count by MPI process rank and partition */
      int** _target_proc_part_num_idx             ;  /*!< Index array of triplet process, partition, numbering for each target */
      int** _target_proc_part_num                 ;  /*!< Array of triplet process, partition, numbering for each target */

      /* Mesh informations */

      CWP_g_num_t  **_gnum_target                 ;  /*<! Target global numbering by partition */
      double       **_coords_target               ;  /*<! Target coordinates by partition */

      CWP_g_num_t _n_g_elt_over_part              ;  /*!< Number of element of the process (over all the partitions)              */
      CWP_g_num_t _n_g_vtx_over_part              ;  /*!< Number of vertices of the process (over all the partitions)             */
      CWP_g_num_t _n_g_elt_cpl_over_part          ;  /*!< Number of coupled code element of the process (over all the partitions) */
      CWP_g_num_t _n_g_vtx_cpl_over_part          ;  /*!< Number of coupled code vertices of the process (over all the partitions)*/

      int  _n_tot_target                          ;  /*!< Target total number on the process                                       */
      int  _n_tot_target_cpl                      ;  /*!< Number of coupled code target received by the process for interpolation  */    
      int *_n_target                              ;  /*!< Target total number on the process by partition                          */
         
      int *_n_vtx                                 ;  /*!< Vertice total number on the process by partition                         */
      int  _n_tot_vtx                             ;  /*!< Vertice total number on the process                                      */
      
      int *_n_elt                                 ;  /*!< Element total number on the process by partition                         */
      int  _n_tot_elt                             ;  /*!< Element total number on the process                                      */
 
      int  _nb_part_cpl                           ;  /*!< Coupled code mesh partition number                                       */
      int  _nb_part                               ;  /*!< Mesh partition number                                                    */
 
      /* Paradigm structure identifier */
    
      int _id_dist                                ;  /*!< Identifier for the localization object of paradigm */
      int _id_gnum_location                       ;  /*!< Identifier for the global numbering to (process,partition,numbering) triplet object of paradigm */

      /* user targets definition for CWP_FIELD_VALUE_USER field type */
     
      int*          _n_user_targets               ;  /*!< Number of targets defined by the user for CWP_FIELD_VALUE_USER field type        */
      int           _n_tot_user_targets           ;  /*!< Total number of targets defined by the user for CWP_FIELD_VALUE_USER field type  */
      double**      _coords_user_targets          ;  /*!< Target coordinates defined by the user for CWP_FIELD_VALUE_USER field type       */
      CWP_g_num_t** _gnum_user_targets            ;  /*!< Target global numbering defined by the user for CWP_FIELD_VALUE_USER field type  */
    
      int _pdmGNum_handle_index;

  }; //end MappingLocation
  

  
}
#endif //__MAPPINGLOCATION_H__


