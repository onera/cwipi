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

/**
 * \cond
 */
namespace cwipi {
    typedef enum {
        CWP_INTERP_AT_SEND,
        CWP_INTERP_AT_RECV
    } CWP_INTERP_TIME;

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

    static void _transform_to_index(int** array,int l1, int l2);

    virtual void init(Coupling *coupling, CWP_Dof_location_t pointsCloudLocation, bool slave);

    /***********************************************************
     **           Mesh information functions                  **
     ***********************************************************/

    /**
      *
      * \brief Get informations from the code mesh to use
      * in SpatialInterp object.
      *
      */

    void mesh_info_get();

    /**
      *
      * \brief Get informations from the coupled code mesh to use
      * in SpatialInterp object.
      *
      */

    void mesh_cpl_info_get();

    /**
      *
      * \brief Get informations from local and coupled code mesh
      *  to use in SpatialInterp object.
      *
      */

    void info_mesh() ;

    virtual void spatialInterpWeightsCompute(CWP_Field_exch_t Texch_t) =0;

    virtual void* interpolate (Field* referenceField) = 0;

    void user_target_points_set(int i_part, int n_pts, double* coord);

    void user_targets_gnum_compute();


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

    void both_codes_on_the_same_process_exchange_p2p (Field* sendingField,
                                                      Field* recevingField
    ) ;

    /**
  *
  * \brief Reception of the communication tree array index
  *        containing localization informations of the
  *        coupled mesh point cloud.
  *
  */

    void data_index_communication_recv_p2p()    ;

    /**
      *
      * \brief Send and reception of the communication tree
      *        array index containing localization informations
      *        in a case where the both are on the same MPI process.
      *
      */

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

    void data_communication_null()          ;

    void data_communication_send_p2p()      ;
    void data_communication_recv_p2p()      ;
    void both_data_communication_p2p()      ;

    void data_communication_wait_send()     ;
    void data_communication_wait_recv()     ;

    void computeFree()                      ;

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

    void data_index_communication_send_p2p()    ;

    /**
     *
     * \brief Return the number of uncomputed targets
     *
     * \return                Number of uncomputed targets
     *
     */

    int
    nUncomputedTargetsGet(int i_part) {
      return n_uncomputed_tgt[i_part];
    }

    /**
     *
     * \brief Return uncomputed targets
     *
     * \return                Uncomputed targets
     *
     */

    inline const int *
    uncomputedTargetsGet() const;

    /**
     *
     * \brief Return the number of computed targets
     *
     * \return                Number of computed targets
     */

    inline int
    nComputedTargetsGet() const;

    /**
     *
     * \brief Return computed targets
     *
     *
     * \return                Computed targets
     *
     */

    inline const int *
    computedTargetsGet() const;


    Coupling                            *_cpl                   ;
    Mesh                                *_mesh                  ;  /*!< Interface Mesh */

    bool _both_codes_are_local{};
    bool _slave{};
    CWP_Field_exch_t _Texch_t{};

    int *_n_target                                              ;  /*!< Target total number on the process by partition */

    /* Mesh informations */
    CWP_g_num_t  **_gnum_target                                 ;  /*!< Target global numbering by partition */
    double       **_coords_target                               ;  /*!< Target coordinates by partition */

protected:
    //Pointer to other objects
    Visu                                *_visu{}                  ;    /*!< Visualization object */
    CodeProperties                      *_localCodeProperties{}   ;
    CodeProperties                      *_coupledCodeProperties{} ;

    CWP_INTERP_TIME                     interpolation_time      ;

    SpatialInterp                         *_spatial_interp_cpl{}  ;  /*!< Spatial interpolation (for both codes are local case) */

    CWP_Dof_location_t                  _pointsCloudLocation{}    ;  /*!< Type of points cloud treated by this mapping instance (cell centers, vertices or user defined) */

    /* code Properties */
    int _id{};
    int _id_cpl{};
    string coupledName{};
    string localName{};

    int _senderRank{};
    int _senderRank_cpl{};
    int _senderLocalRank{};

    int  _nb_part_cpl{}                           ;  /*!< Coupled code mesh partition number                                       */
    int  _nb_part{}                               ;  /*!< Mesh partition number                                                    */

    /** MPI processes informations **/

    /* MPI Communicators */
    MPI_Comm _globalComm{};         // Gathers every processus
    MPI_Comm _cplComm{};            // Processus involved in the coupling in either code
    MPI_Comm _localComm{};          // Processus involved in the coupling for the local code
    PDM_MPI_Comm  _pdm_cplComm{};   // _cplComm for Paradigm

    int **_weights_src_idx{};
    double **_weights_src{};

    vector<string> _codeVector;

    std::vector<int>* _connectableRanks_cpl{};
    std::vector<int>* _connectableRanks{}    ;

    /* informations about MPI process (rank) */
    bool _isActiveRank{};
    int cplComm_rank{};       // Rank in cplComm
    int cplComm_size{};       // Size of cplComm
    int localComm_size{};     // Size of localComm
    int localComm_size_cpl{}; // Size of localComm of the coupled code

    std::vector<int> n_uncomputed_tgt;

    int         **_targets_localization_idx{}     ;  /*!< Data index (by process and by partition) of target localization*/
    target_data  *_targets_localization_data{}    ;  /*!< Data of target localization */
    int         **_targets_localization_idx_cpl{} ;  /*!< Data index (by process and by partition) of the received target localization*/
    target_data  *_targets_localization_data_cpl{};  /*!< Data of the received target localization */

    //TODO: To delete and replace by using other members
    std::vector<int>   _idx_target              ;  /*!< Index of the number of target by partition */

    /* Displacement and count for all_to_all MPI communication of targets_localization_data */

    int* _targets_localization_data_count_recv{}  ;  /* Counts for all_to_all MPI communication of targets_localization_data (reception) */
    int* _targets_localization_data_count_send{}  ;  /* Counts for all_to_all MPI communication of targets_localization_data (sending) */
    int* _targets_localization_data_disp_recv{}   ;  /* Displacements for all_to_all MPI communication of targets_localization_data (reception) */
    int* _targets_localization_data_disp_send{}   ;  /* Displacements for all_to_all MPI communication of targets_localization_data (sending) */

    /* Triplet global numbering, MPI process, partition results */

    int** _process_and_partition_count{}          ;  /*!< Element count by MPI process rank and partition */

    int  _n_tot_target{}                          ;  /*!< Target total number on the process                                       */
    int  _n_tot_target_cpl{}                      ;  /*!< Number of coupled code target received by the process for interpolation  */

    /* user targets definition for CWP_DOF_LOCATION_USER field type */

    int*          _n_user_targets{}               ;  /*!< Number of targets defined by the user for CWP_DOF_LOCATION_USER field type        */
    int           _n_tot_user_targets{}           ;  /*!< Total number of targets defined by the user for CWP_DOF_LOCATION_USER field type  */
    CWP_g_num_t** _gnum_user_targets{}            ;  /*!< Target global numbering defined by the user for CWP_DOF_LOCATION_USER field type  */

    CWP_g_num_t _n_g_elt_over_part{}              ;  /*!< Number of element of the process (over all the partitions)              */
    CWP_g_num_t _n_g_vtx_over_part{}              ;  /*!< Number of vertices of the process (over all the partitions)             */
    CWP_g_num_t _n_g_elt_cpl_over_part{}          ;  /*!< Number of coupled code element of the process (over all the partitions) */
    CWP_g_num_t _n_g_vtx_cpl_over_part{}          ;  /*!< Number of coupled code vertices of the process (over all the partitions)*/

    /* Mesh informations */
    double**      _coords_user_targets{}          ;  /*!< Target coordinates defined by the user for CWP_DOF_LOCATION_USER field type       */

    int *_n_vtx{}                                 ;  /*!< Vertice total number on the process by partition                         */
    int  _n_tot_vtx{}                             ;  /*!< Vertice total number on the process                                      */

    int *_n_elt{}                                 ;  /*!< Element total number on the process by partition                         */
    int  _n_tot_elt{}                             ;  /*!< Element total number on the process                                      */
  };

    /**
     * \endcond
     */

}

#endif //__SPATIAL_INTERP_H__
