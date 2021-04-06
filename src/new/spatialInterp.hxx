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

#include "mesh.hxx"
#include "field.hxx"
#include "codeProperties.hxx"
#include "coupling.hxx"

/**
 * \cond
 */
namespace cwipi {


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
  friend class SpatialInterpLocation;
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

    static void _transform_to_index(int* array,int l1);

    virtual void init(Coupling *coupling, CWP_Dof_location_t pointsCloudLocation,int slave) =0;

    virtual void spatialInterpWeightsCompute(CWP_Field_exch_t Texch_t) =0;

    virtual void* interpolate (Field* referenceField) = 0;

    virtual void user_target_points_set(int i_part, int n_pts, double* coord) =0;

    void user_targets_gnum_compute();

    /**
     *
     * \brief Exchange data field with the coupled application with blocking
     *        communications.
     *
     * This function exchanges interpolated fieldsDouble between coupled codes.
     *
     * \warning  The size of tgt_field_id size is n_computed_tgt.
     *           If \f$ n\_uncomputed\_tgt \ne n\_tgt\_pts \f$,
     *           user himself must set values for uncomputed target points.
     *
     * \param [in]  src                       Source field (NULL -> no sending)
     * \param [in]  tgt                       Target field (NULL -> no receiving)
     * \param [in]  ptFortranInterpolationFct Fortran user interpolation (or NULL)
     * \param [out] n_uncomputed_tgt          Number of uncomputed target
     *
     * \return                                Exchange status
     *
     */

    CWP_Err_t
    sendRecv
    (Field *src,
     Field *tgt,
     void      *ptFortranInterpolationFct,
     int       *n_uncomputed_tgt);

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
    void null_exchange_for_uncoupled_process_p2p() ;

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

    void data_index_communication_send()    ;


    void data_index_communication_send_p2p()    ;

    /**
     * \brief Setting user target points
     *
     * This function must be called if the nature of receiving fieldsDouble
     * is \ref CWP_DOF_LOCATION_USER
     *
     * \param [in]  n_pts   Number of points
     * \param [in]  coords   Coordinates (size = 3 * n_pts)
     *
     */

    void
    userTgtPtsSet
    (const int            n_pts,
     double               coords[]);

    /**
     *
     * \brief Setting of an user interpolation from location.
     *
     * This function takes into account an user interpolation function written with
     * void (*\ref CWP_Interp_from_location_t) interface.
     *
     * \param [in] fct        Function
     *
     */

    void
    InterpUser
    (CWP_Interp_from_location_t fct);

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

    void _IAlltoallIndexSend(void* send_buffer,
                             int* send_count,
                             int* send_disp,
                             MPI_Datatype type,
                             MPI_Comm comm,
                             std::vector<int> connectableRanks
                            );

    void _IAlltoallIndexRecv(void* recv_buffer,
                             int* recv_count,
                             int* recv_disp,
                             MPI_Datatype type,
                             MPI_Comm comm,
                             std::vector<int> connectableRanks
                             );

    //TODO: Acess function
    int _both_codes_are_local;

    int slaveGet() {
      return  _slave;
    }

   int bothLocalGet() {
      return  _both_codes_are_local;
    }

  protected:
     //Pointer to other objects
    Mesh                                *_mesh                  ;    /*!< Interface Mesh       */
    Visu                                *_visu                  ;    /*!< Visualization object */
    CodeProperties                      *_localCodeProperties   ;
    CodeProperties                      *_coupledCodeProperties ;
    Coupling                            *_cpl                   ;

    SpatialInterp                       *_spatial_interp_cpl    ;  /*!< Spatial interpolation (for both codes are local case) */

    CWP_Field_exch_t                     _Texch_t               ;

    /* code Properties */
    int _id;
    int _id_cpl;
    string coupledName;
    string localName;

    int _senderRank;
    int _senderRank_cpl;
    int _senderLocalRank;

    int  _nb_part_cpl                           ;  /*!< Coupled code mesh partition number                                       */
    int  _nb_part                               ;  /*!< Mesh partition number                                                    */

    /** MPI processes informations **/

   /* MPI Communicators */
   MPI_Comm _globalComm;         // Gathers every processus
   MPI_Comm _cplComm;            // Processus involved in the coupling in either code
   MPI_Comm _localComm;          // Processus involved in the coupling for the local code
   PDM_MPI_Comm  _pdm_cplComm;   // _cplComm for Paradigm
   PDM_MPI_Comm  _pdm_localComm; // _localComm for Paradigm

   vector<string> _codeVector;

   int  _slave;

   std::vector<int>* _connectableRanks_cpl;
   std::vector<int>* _connectableRanks    ;

   /* informations about MPI process (rank) */
   bool _isCoupledRank;

   int _rank;
   int  _n_ranks    ;
   int  _n_ranks_cpl;
   int  _n_ranks_g    ;

   /* MPI Request */
   std::vector<MPI_Request> _send_requests;
   std::vector<MPI_Request> _recv_requests;

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
    int *_n_target                              ;  /*!< Target total number on the process by partition                          */

    /* user targets definition for CWP_DOF_LOCATION_USER field type */

    int*          _n_user_targets               ;  /*!< Number of targets defined by the user for CWP_DOF_LOCATION_USER field type        */
    int           _n_tot_user_targets           ;  /*!< Total number of targets defined by the user for CWP_DOF_LOCATION_USER field type  */
    double**      _coords_user_targets          ;  /*!< Target coordinates defined by the user for CWP_DOF_LOCATION_USER field type       */
    CWP_g_num_t** _gnum_user_targets            ;  /*!< Target global numbering defined by the user for CWP_DOF_LOCATION_USER field type  */

  };

    /**
     * \endcond
     */

}

#endif //__SPATIAL_INTERP_H__
