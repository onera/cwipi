#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__
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

  class Mesh;
  class Field;  
  class Visu;
  /** 
   * \class Geometry geometry.hxx "geometry.hxx"
   * \brief Geometry algorithm
   *
   *  This class computes the geometry algotrithm of points cloud into a mesh and
   *  builds a communication graph to transmit interpolated fieldsDouble on 
   *  points cloud from fieldsDouble defined on the mesh.
   * 
   */

  class Geometry {
  friend class GeomLocation;
  public:

    /**
     * \brief Constructor
     *
     */

    Geometry();

    /**
     * \brief Destructor
     *
     */

    virtual ~Geometry();

    virtual void init(Coupling *coupling, CWP_Field_value_t geometryLocation,int slave) =0;

    virtual void compute(CWP_Field_exch_t Texch_t) =0;

    virtual void user_target_points_set(int i_part, int n_pts, double* coord) =0;
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
     * \brief Sending of data field to the coupled application with nonblocking 
     *        communications.
     *
     * This function sends interpolated field to the coupled code. 
     * 
     * \param [in]  sendingField                      Sending field    
     *
     *
     */

    virtual void issend (Field* sendingField) = 0;


    virtual void null_exchange_for_uncoupled_process () = 0;

    virtual void both_codes_on_the_same_process_exchange (Field* sendingField,Field* recevingField) = 0;

    /**
     *
     * \brief Waiting of the end of exchange related to request.
     *
     * This function waits the end of exchange related to request
     * from \ref CWP_Issend
     * 
     */

    virtual void waitIssend(Field* sendingField) = 0;

    /**
     *
     * \brief Receiving of Data field from the coupled application with nonblocking 
     *        communications.
     *
     * This function receives interpolated field from the coupled code 
     * 
     * \param [in]  recevingField       Receving field   
     *
     *
     */

    virtual void irecv (Field* recevingField) = 0;

    /**
     *
     * \brief Waiting of the end of exchange related to request.
     *
     * This function waits the end of exchange related to request 
     * from \ref CWP_Irecv
     * 
     * \param [in] request    Request to wait the end of exchange
     *
     */

    virtual void waitIrecv(Field* recevingField) = 0;    
   

    /**
     * \brief Setting user target points
     *
     * This function must be called if the nature of receiving fieldsDouble 
     * is \ref CWP_FIELD_VALUE_USER_TO_NODE
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
     * void (*\ref CWP_Interp_from_target_proc_part_num_t) interface.
     * 
     * \param [in] fct        Function
     *
     */

    void 
    InterpUser
    (CWP_Interp_from_target_proc_part_num_t fct);

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

     
  protected:
    
    Geometry &operator=(const Geometry &other);  /*!< Assigment operator not available */
    Geometry (const Geometry& other);            /*!< Copy constructor not available */

     //Pointer to other objects
    Mesh                                *_mesh                  ;    /*!< Interface Mesh       */
    Visu                                *_visu                  ;    /*!< Visualization object */
    CodeProperties                      *_localCodeProperties   ;
    CodeProperties                      *_coupledCodeProperties ;
    Coupling                            *_cpl                   ;
    
    std::map <std::string,Field*>       *_referenceFieldsDB     ;
    CWP_Field_exch_t                     _Texch_t               ;

    /* code Properties */
    int _id;
    int _id_cpl;
    string coupledName;
    string localName; 

    int _senderRank;  
    int _senderRank_cpl;

    int _senderLocalRank;

   /** MPI processes informations **/
   
   /* MPI Communicators */
   MPI_Comm _globalComm ;
   MPI_Comm _unionComm ;   
   MPI_Comm _localComm  ;
   MPI_Comm _connectableComm  ;   
   PDM_MPI_Comm  _pdm_connectableComm ;   
   PDM_MPI_Comm  _pdm_localComm ;
   PDM_MPI_Comm  _pdm_globalComm ;
   PDM_MPI_Comm  _pdm_unionComm ;
   
   vector<string> _codeVector;

   int  _slave;   

   const std::vector<int>* _connectableRanks_cpl;
   const std::vector<int>* _connectableRanks    ;  

   /* informations about MPI process (rank) */
   bool _isCoupledRank;
   bool _isCoupledRank_cpl;

   int _rank;    
   int  _n_ranks    ;
   int  _n_ranks_cpl;
   int  _n_ranks_g    ;
   
   /* MPI Request */
   std::vector<MPI_Request> _send_requests;
   std::vector<MPI_Request> _recv_requests;  
   
   std::vector<int> n_uncomputed_tgt;
   
  };



}

#endif //__GEOMETRY_H__
