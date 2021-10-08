#ifndef __SPATIAL_INTERP_H__
#define __SPATIAL_INTERP_H__
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

#include <cmath>

#include "mesh.hxx"
#include "field.hxx"
#include "codeProperties.hxx"
#include "coupling.hxx"
#include "pdm_part1_to_selected_part2.h"

/**
 * \cond
 */
namespace cwipi {

// A conserver et à renommer !!!

  typedef enum {
    CWP_SPATIAL_INTERP_AT_SEND,
    CWP_SPATIAL_INTERP_AT_RECV
  } CWP_SpatialInterp_time_t;

// A supprimer !!!

  struct target_data {
    int          lnum    ;
    int          origin_part     ;
    CWP_g_num_t  closest_elt_gnum;
    int          origin_proc     ;
    int          closest_elt_part;
    int          l_num_origin    ;
    double       projectedX       ;
    double       projectedY       ;
    double       projectedZ       ;
    double       distance       ;
  };


  typedef enum {
    SPATIAL_INTERP_EXCH_SEND,
    SPATIAL_INTERP_EXCH_RECV
  } SpatialInterpExchDirection;

// A supprimer !!!

  static const char *CWP_Field_exch_t_str [] = {"CWP_FIELD_EXCH_SEND","CWP_FIELD_EXCH_RECV","CWP_FIELD_EXCH_SENDRECV"};

  class Mesh;
  class Field;
  class Visu;

  /**
   *
   * \class SpatialInterp spatialInterp.hxx "spatialInterp.hxx"
   * \brief SpatialInterp algorithm
   *
   *  This class computes the spatial interpolation weights of points cloud into a mesh and
   *  builds a communication graph to transmit interpolated fieldsDouble on
   *  points cloud from fieldsDouble defined on the mesh.
   *
   */

  class SpatialInterp {

  public:

    /**
     * \brief Constructor
     *
     */

    SpatialInterp();

    /**
     * \brief Destructor
     *
     */

    virtual ~SpatialInterp();

    virtual void 
    init (
      Coupling           *coupling, 
      CWP_Dof_location_t localCodeDofLOcation,
      CWP_Dof_location_t cplCodeDofLOcation,
      SpatialInterpExchDirection exchDirection );

    virtual void weightsCompute()  = 0; // A changer de nom
    
    virtual void* interpolate (Field* referenceField) = 0;

    void issend(Field* referenceField);

    void irecv(Field* recevingField);


    /**
     *
     * \brief Return the number of uncomputed targets
     *
     * \return                Number of uncomputed targets
     *
     */

    int
    nUncomputedTargetsGet(int i_part) const;

    /**
     *
     * \brief Return uncomputed targets
     *
     * \return                Uncomputed targets
     *
     */

    const int *
    uncomputedTargetsGet(int i_part) const;

    /**
     *
     * \brief Return the number of computed targets
     *
     * \return                Number of computed targets
     */

    int
    nComputedTargetsGet(int i_part) const;

    /**
     *
     * \brief Return computed targets
     *
     *
     * \return                Computed targets
     *
     */

    const int *
    computedTargetsGet(int i_part) const;


  protected:

    Coupling                   *_cpl;
    Mesh                       *_mesh;                  /*!< Interface Mesh */

    Visu                       *_visu;                  /*!< Visualization object */
    CodeProperties             *_localCodeProperties;   
    CodeProperties             *_coupledCodeProperties; 

    CWP_Dof_location_t        _localCodeDofLocation;     /*!< Type of points cloud treated by this mapping instance (cell centers, vertices or user defined) */
    CWP_Dof_location_t        _coupledCodeDofLocation;   /*!< Type of points cloud treated by this mapping instance (cell centers, vertices or user defined) */

    SpatialInterpExchDirection  _exchDirection;  /*!< Spatial interpolation (for both codes are local case) */
   
    PDM_part1_to_selected_part2_t *_ptsp;
    
    CWP_SpatialInterp_time_t     _interpolation_time      ;

    int  _nPart; /*!< Mesh partition number                                                    */
    int  _cplNPart;  /*!< Coupled mesh partition number                                                    */

    int _rootRankUnionComm   ;
    int _cplRootRankUnionComm;

    int _rootRankCplComm   ;
    int _cplRootRankCplComm;

    MPI_Comm _cplComm;
    PDM_MPI_Comm _pdmCplComm;
    MPI_Comm _unionComm;
    PDM_MPI_Comm _pdmUnionComm;

    MPI_Comm _localComm;          // Processus involved in the coupling for the local code




  // A conserver ou supprimer 
  protected:
    /* code Properties */
    int _id;
    int _id_cpl;
    string coupledName;
    string localName;


    int  _nPart_cpl                           ;  /*!< Coupled code mesh partition number                                       */


    /* informations about MPI process (rank) */
    int cplComm_rank;       // Rank in cplComm
    int cplComm_size;       // Size of cplComm
    int localComm_size;     // Size of localComm
    int localComm_size_cpl; // Size of localComm of the coupled code


    int *_n_elt_weights;
    int **_weights_idx;
    double **_weights;

    int *_n_computed_tgt;
    int **_computed_tgt;

    int *_n_uncomputed_tgt;
    int **_uncomputed_tgt;


  // A supprimer
  protected:
    bool _both_codes_are_local;
    bool _slave;


    int _senderRank;
    int _senderRank_cpl;
    int _senderLocalRank;


    vector<string> _codeVector;

    std::vector<int>* _connectableRanks_cpl;
    std::vector<int>* _connectableRanks    ;


    std::vector<int> n_uncomputed_tgt;

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

    int  _n_tot_target                          ;  /*!< Target total number on the process                                       */
    int  _n_tot_target_cpl                      ;  /*!< Number of coupled code target received by the process for interpolation  */

    /* user targets definition for CWP_DOF_LOCATION_USER field type */

    int*          _n_user_targets               ;  /*!< Number of targets defined by the user for CWP_DOF_LOCATION_USER field type        */
    int           _n_tot_user_targets           ;  /*!< Total number of targets defined by the user for CWP_DOF_LOCATION_USER field type  */
    CWP_g_num_t** _gnum_user_targets            ;  /*!< Target global numbering defined by the user for CWP_DOF_LOCATION_USER field type  */

    CWP_g_num_t _n_g_elt_over_part              ;  /*!< Number of element of the process (over all the partitions)              */
    CWP_g_num_t _n_g_vtx_over_part              ;  /*!< Number of vertices of the process (over all the partitions)             */
    CWP_g_num_t _n_g_elt_cpl_over_part          ;  /*!< Number of coupled code element of the process (over all the partitions) */
    CWP_g_num_t _n_g_vtx_cpl_over_part          ;  /*!< Number of coupled code vertices of the process (over all the partitions)*/


    int *_n_vtx                                 ;  /*!< Vertice total number on the process by partition                         */
    int  _n_tot_vtx                             ;  /*!< Vertice total number on the process                                      */

    int *_n_elt                                 ;  /*!< Element total number on the process by partition                         */
    int  _n_tot_elt                             ;  /*!< Element total number on the process                                      */


  // A supprimer ou Rendre privé 
  public:

    static void _transform_to_index(int** array,int l1, int l2);
    void mesh_info_get();
    void mesh_cpl_info_get();
    void issend_p2p(Field* referenceField);
    void info_mesh() ;
    void waitIssend_p2p(Field* referenceField);
    void irecv_p2p(Field* recevingField);
    void waitIrecv_p2p (Field* recevingField);
    void null_exchange_for_uncoupled_process() ;
    void both_codes_on_the_same_process_exchange_p2p (Field* sendingField,
                                                      Field* recevingField);
    void data_index_communication_recv_p2p()    ;
    void both_index_communication_p2p()    ;
    void data_index_communication_null();

    void prepare_data_communication_send()  ;
    void prepare_data_communication_recv()  ;

    void data_communication_null()          ;

    void data_communication_send_p2p()      ;
    void data_communication_recv_p2p()      ;
    void both_data_communication_p2p()      ;

    void data_communication_wait_send()     ;
    void data_communication_wait_recv()     ;

    void computeFree()                      ;

    void data_index_communication_send_p2p()    ;

  };

}

#endif //__SPATIAL_INTERP_H__
