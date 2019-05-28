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
#include "field.hpp"
#include "codeProperties.hxx"
#include "coupling.hxx"

namespace cwipi {


  typedef struct target_data {
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

  typedef struct target_data_vtx {
    int          lnum    ;
    int          inc     ;
    CWP_g_num_t  closest_elt_gnum;
    int          origin_proc     ;
    int          closest_elt_part;
    int          l_num_origin    ;
    double       projectedX       ;     
    double       projectedY       ;
    double       projectedZ       ;      
  };


  class Mesh;

  /** 
   * \class Geometry geometry.hxx "geometry.hxx"
   * \brief Geometry algorithm
   *
   *  This class computes the geometry algotrithm of points cloud into a mesh and
   *  builds a communication graph to transmit interpolated fields on 
   *  points cloud from fields defined on the mesh.
   * 
   */

  class Geometry {
    
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


     //Depends on the Geometry type
     // Ca sera dans les classes concrètes
     virtual void locate_setting_request(int id_dist) =0;
     virtual void locate_setting_surface(int id_dist) =0;    
     virtual void locate_compute(int id_dist)             =0;
     virtual void locate_get(int id_dist)      =0;
     virtual void locate_get_cpl(int id_dist)      =0;     
     
     
      virtual void broadcasting_request(int id_gnum_location)  =0;
      virtual void broadcasting_set    (int id_gnum_location)  =0;
      virtual void location_compute                 (int id_gnum_location)  =0;
     
     virtual void location_get(int id_gnum_location)      =0;
     virtual void location_get_cpl(int id_gnum_location)      =0;
     virtual void broadcasting_filling_of_broadcasting_array() =0;
     virtual void broadcasting_index_communication() =0;
     
     virtual void broadcasting_communication() =0;
     virtual void broadcasting_communication2() =0;
     virtual void broadcasting_wait_and_targets_array_filling() =0;
     
     //Depends on the Geometry type
     // Ca sera dans les classes concrètes
     virtual double* interpolate(Field <double>* referenceField) =0;  

     void init(Coupling *coupling, CWP_Field_value_t geometryLocation);
     void mesh_info_get();
     void mesh_cpl_info_get();
     void mesh_cpl_info_get2();
     void compute(int *n_uncomputed_tgt);
     inline Geometry* getCoupledGeometry();


     int nTargetGet(int i_part) {
       return _n_target[i_part];
     }

    /**
     *
     * \brief Exchange data field with the coupled application with blocking 
     *        communications.
     *
     * This function exchanges interpolated fields between coupled codes. 
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
    (Field <double> *src,
     Field <double> *tgt,
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

    virtual void  
    issend
    (Field <double>* sendingField
    ) = 0;

    /**
     *
     * \brief Waiting of the end of exchange related to request.
     *
     * This function waits the end of exchange related to request
     * from \ref CWP_Issend
     * 
     */

    void 
    waitIssend
    (Field <double>* sendingField);

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


    void 
    irecv
    (Field <double>* recevingField
    );

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

    void 
    waitIrecv
    (Field <double>* recevingField);    
   

    /**
     * \brief Setting user target points
     *
     * This function must be called if the nature of receiving fields 
     * is \ref CWP_FIELD_VALUE_USER
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

    inline int 
    nUncomputedTargetsGet() const;

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

    Mesh*
    meshGet() {
      return _mesh;
    }

    CWP_g_num_t*
    closestEltGnumGet(int i_part) {
      return _closest_elt_gnum[i_part];
    }    
    
    void
    closestEltGnumSet(CWP_g_num_t* closest_elt_gnum, int i_part) {
      _closest_elt_gnum[i_part]=closest_elt_gnum;
    }    

    void
    closestEltGnumSet(CWP_g_num_t closest_elt_gnum, int i_part,int i_el) {
      _closest_elt_gnum[i_part][i_el]=closest_elt_gnum;
    }    


    CWP_g_num_t*
    gnumTargetGet(int i_part) {
      return _gnum_target[i_part];
    }  


    double*
    distanceTargetGet(int i_part) {
      return _distance[i_part];
    }  
        
    
    double*
    coordsTargetGet(int i_part) {
      return _coords_target[i_part];
    }  
          

    

 void _IAlltoallIndex(void* send_buffer,
                int* send_count,
                int* send_disp,
                void* recv_buffer,
                int* recv_count,
                int* recv_disp,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableRanks
                );


  void _IAlltoall2(int** send_buffer,
                int* send_size,
                int send_stride,
                int** recv_buffer,
                int* recv_size,
                int recv_stride,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableRanks
                );


    void _Wait();      

    double      **_distance         ; 
    double      **_projected        ; 
    CWP_g_num_t **_closest_elt_gnum ;  


  int** _location_count_comm_proc ;
  int** _location_idx_proc_recv   ; 

  int** _location_idx ;
  int** _location     ;


  int* _location_count_recv;         
  int* _location_count_send ;
  int* _location_disp_recv ;         
  int* _location_disp_send ; 


  protected:
    
    Geometry &operator=(const Geometry &other);  /*!< Assigment operator not available */
    Geometry (const Geometry& other);            /*!< Copy constructor not available */

  protected:
  
    // Informations about locations (from the local proc send to a distant proc)
    // Local Targets

    std::vector<int>                      _idx_target        ;    


void _IBcast(void* send_buffer,
                     int send_size,
                     int send_stride,
                     void* recv_buffer,
                int recv_size,
                int recv_stride,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableRanks_cpl,
                int _n_ranks_cpl,
                int rootRank);
 
  void _IAlltoall(void* send_buffer,
                int* send_size,
                int send_stride,
                void* recv_buffer,
                int* recv_size,
                int recv_stride,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableRanks_cpl,
                int _n_ranks_cpl);

    
    //Distant Targets (from a distant proc send to the local proc)



    Mesh                                *_mesh;      /*!< Interface Mesh */
    Visu                                *_visu;
    CodeProperties                      *_localCodeProperties;
    CodeProperties                      *_coupledCodeProperties;
    std::map <std::string,Field<double>*>    *_referenceFieldsDB;
    std::map <int,Field<double>*>        _requestFieldsDB;
    CWP_Field_value_t                    _geometryLocation;
    Coupling                            *_cpl;
    
    int** _targets_cpl_idx;
    target_data* _targets_cpl;
    target_data_vtx* _targets_vtx_cpl;
    int** _targets_cpl_idx_cpl;
    target_data* _targets_cpl_cpl;
    target_data_vtx* _targets_vtx_cpl_cpl;       
    
    int  _option;
    int  _n_tot_target;
    int  _n_tot_target_cpl;

    int* _n_targets_dist_proc; 


    const std::vector<int>* _connectableRanks_cpl;
    const std::vector<int>* _connectableRanks    ;  
    int  _n_ranks    ;
    int  _n_ranks_cpl;
    int  _n_ranks_g    ;

    
    int _nb_part;
    int _nb_part_cpl;


    int* _n_vtx;   
    int* _n_vtx_cpl;   
    int* _n_elt;
    int* _n_elt_cpl;
    
    int* _n_target;
    CWP_g_num_t** _gnum_target;
    double** _coords_target; 

    
    int _n_tot_elt;
    int _n_tot_vtx;
    
    int _n_g_elt_over_part;
    int _n_g_vtx_over_part;
    int _n_g_elt_cpl_over_part;
    int _n_g_vtx_cpl_over_part;  
    
    int**                          _n_targets_dist_proc_dist_part;
    int**                          _n_targets_recv_dist_proc_loc_part    ;
    int**                          _idx_targets_recv_dist_proc_loc_part  ;
    
    int* _n_recv_targets_loc_from_dist_proc;
 
    int* _n_elt_exch_recv_by_proc_by_part;
    int* _n_vtx_exch_recv_by_proc_by_part;
        
    int* _n_elt_exch_recv_tot_per_part;
    int* _n_vtx_exch_recv_tot_per_part; 

   int _n_tot_elt_exch_cpl;
   int _n_tot_vtx_exch_cpl;


   int* _n_g_elt;
   int* _n_g_vtx;
 
   double* _centers_conc ;
   double* _coords_conc ;
   CWP_g_num_t* _gnum_elt_conc ;
   CWP_g_num_t* _gnum_vtx_conc ;
   int* _lnum_elt_conc ;
   int* _lnum_vtx_conc ;


  double* _centers_conc_cpl       ;
  double* _coords_conc_cpl        ;
  CWP_g_num_t* _gnum_elt_conc_cpl ;
  int* _lnum_elt_conc_cpl         ; 
  CWP_g_num_t* _gnum_vtx_conc_cpl ;
  int* _lnum_vtx_conc_cpl         ; 

 
    

    double      **_distance_cpl          ; 
    double      **_projected_cpl         ; 
    CWP_g_num_t **_closest_elt_gnum_cpl ;  

  int _id_dist1 ;
  int _id_dist2 ;

  int _id_gnum_location1;
  int _id_gnum_location2;


  int** _location_idx_comm_proc   ;
  target_data* _location_comm_proc;
  target_data* _location_recv;



   MPI_Comm _globalComm ;
   MPI_Comm _localComm  ;
   PDM_MPI_Comm  _pdm_localComm ;
   PDM_MPI_Comm  _pdm_globalComm ;

   Geometry* _geometry_cpl;

   int _both_codes_are_local; 
   int* _both_codes_are_local__array;
   
   int _rank;
   int _localRank;   
   vector<string> _codeVector;
   string coupledName;
   string localName; 
   
   std::vector<int> _send_requests;
   std::vector<int> _recv_requests;   
   
   
  };

  Geometry* Geometry::getCoupledGeometry() {
  
  }



}

#endif //__GEOMETRY_H__
