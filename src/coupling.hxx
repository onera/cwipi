#ifndef __COUPLING_H__
#define __COUPLING_H__
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

#include <string>
#include <map>
#include <vector>

#include <fvmc_nodal.h>

#include "cwp.h"
#include "bftc_printf.h"
#include "communication.hxx"
#include "couplingDB.hxx"
#include "couplingDB_i.hxx"

//#include "geometry.hxx"
//#include "visualization.hxx"
//#include "support.hxx"
//#include "field.hpp"

using namespace std;

namespace cwipi {

  class CodeProperties;

  class Mesh;

  class LocationToDistantMesh;

  class LocationToLocalMesh;

  /** 
   * \class Coupling coupling.hxx "coupling.hxx"
   * \brief Coupling between two codes.
   *
   *  This class defines a coupling object and the associated communcations
   *  between two codes  
   * 
   */

  class Coupling {

  public:

    /**
     * \brief Constructor.
     *
     * This function creates a coupling object and defines its properties.
     *
     * \param [in]  cplId                       Coupling identifier
     * \param [in]  comm_type                    Communication type
     * \param [in]  localCodeProperties          Local code properties
     * \param [in]  coupledCodeProperties        Coupled code properties
     * \param [in]  geom_algo                    Geometric algorithm
     * \param [in]  support_type                 Support type
     * \param [in]  n_part                       Number of interface partitions 
     * \param [in]  moving_status                Support moving status
     * \param [in]  recv_freq_type               Type of receiving frequency
     * \param [in]  cplDB                        Coupling data base where it coupling is stored
     *
     */

    Coupling
    (
     const string                &cplId,
     const CWP_Comm_t           commType,
     const CodeProperties        &localCodeProperties,
     const CodeProperties        &coupledCodeProperties,
     const CWP_Geom_t           geomAlgo,
     const CWP_Support_t        supportType,
     const int                    nPart,
     const CWP_Displacement_t  movingStatus,
     const CWP_Freq_t           recvFreqType,
     CouplingDB                 &cplDB 
    );

    /**
     * \brief Destructor.
     *
     */

    virtual ~Coupling();

    /**
     * \brief data exchange <b>(Not implemented yet)</b> 
     *
     * Exchange depending on exchange frequency
     * 
     */

    void
    exchange
    (
    )
    {
      bftc_printf("\n cwipi error : CWIPI_exchange not implemented yet\n");
      exit(1);
    }

    /**
     *
     * \brief Return the number of uncomputed targets
     * 
     * \return                Number of uncomputed targets
     */

    inline int 
    nUncomputedTargetsGet
    (
    ) const;

    /**
     *
     * \brief Return uncomputed targets
     *
     * \return                Uncomputed targets
     */

    inline const int *
    uncomputedTargetsGet
    (
    ) const;

    /**
     *
     * \brief Return the number of computed targets
     * 
     * \return                Number of computed targets
     */

    inline int 
    nComputedTargetsGet
    (
    ) const;

    /**
     *
     * \brief Return computed targets
     * 
     * \return                Computed targets
     */

    inline const int *
    computedTargetsGet
    (
    ) const;

    /*----------------------------------------------------------------------------*
     * Methods about exchange frequency                                           *
     *----------------------------------------------------------------------------*/

    /**
     * \brief Setting receiving frequency.
     *
     * This function set receiving frequency. It must be used when
     * the type of receiving frequency is \ref FREQ_RELATED_N_TIME_STEP
     *
     * \param [in]  n_step     Frequency in steps number
     *
     */

    inline void
    recvFreqSet
    (
     int n_step
    );

    /**
     * \brief Setting the next receiving time.
     *
     * This function set the next receiving time. It must be used when
     * the type of receiving frequency is \ref FREQ_ASYNCHRONOUS
     *
     * \param [in]  next_time     Next receiving time
     *
     */

    inline void
    recvNextTimeSet
    (
     double next_time
    );

    /*----------------------------------------------------------------------------*
     * Methods about geometry                                                     *
     *----------------------------------------------------------------------------*/

    /**
     * \brief Computation geometry                                  
     *
     * This function compute geometry 
     *
     * \param [out] n_uncomputed_tgt    Number of uncomputed target
     *
     */

    void 
    geomCompute
    (
     int *n_uncomputed_tgt
    );

    /**
     * \brief Setting geometry properties
     *
     * This function set the geometric algorithm properties.
     *
     * \param [in]       fmt       Format with the syntax : "prop1, prop2, ..."
     * \param [in,out]   pa        List of properties values
     *
     */

    void 
    geomPropertiesSet
    (
     const char *fmt,
     va_list    *pa
    );

    /*----------------------------------------------------------------------------*
     * Methods about visualization                                                *
     *----------------------------------------------------------------------------*/

    /**
     * \brief Enable visualization output
     *
     * This function enable visualization output.
     *
     * \param [in]  freq             Output frequency
     * \param [in]  format           Output format to visualize exchanged fields
     *                               on the coupled mesh. Choice between :
     *                               - "EnSight Gold"
     *                               - "MED_ficher"
     *                               - "CGNS"
     *                               .
     * \param [in]  format_option   Output options "opt1, opt2, ..." :
     *                         - text               output text files
     *                         - binary             output binary files (default)
     *                         - big_endian         force binary files
     *                                              to big-endian
     *                         - discard_polygons   do not output polygons
     *                                              or related values
     *                         - discard_polyhedra  do not output polyhedra
     *                                              or related values
     *                         - divide_polygons    tesselate polygons
     *                                              with triangles
     *                         - divide_polyhedra   tesselate polyhedra
     *                                              with tetrahedra and pyramids
     *                                              (adding a vertex near
     *                                               each polyhedron's center)
     *                         .
     *
     */

    inline void 
    visuSet
    (
     const int   freq,
     const char *format,
     const char *format_option
    );

    /*----------------------------------------------------------------------------*
     * Methods about User target points                                           *
     *----------------------------------------------------------------------------*/

    /**
     * \brief Setting user target points
     *
     * This function must be called if the nature of receiving fields 
     * is \ref CWP_FIELD_LOCATION_USER
     *
     * \param [in]  n_pts   Number of points
     * \param [in]  coord   Coordinates (size = 3 * \ref n_pts)          
     *
     */

    void 
    userTgtPtsSet
    (
     const int n_pts,
     double    coord[]
    );

    /*----------------------------------------------------------------------------*
     * Methods  about Support                                                     *
     *----------------------------------------------------------------------------*/

    /**
     * \brief Setting vertices
     *
     * This method set partition vertices
     *
     * \param [in]  i_part      Current partition
     * \param [in]  n_pts       Number of points
     * \param [in]  coord       Coordinates (size = 3 * \ref n_pts)          
     * \param [in]  parent_num  Pointer to parent element number (or NULL)
     *
     */

    inline void 
    supportVtcsSet
    (
     const int          i_part,
     const int          n_pts,
     const double       coord[],
     const CWP_g_num_t parent_num[]
    );

    /**
     * \brief End setting support
     *
     * This function finalizes the support building
     *
     */

    inline void 
    supportEndSet
    (
    );

    /**
     * \brief Adding a connectivity block to the geometric support
     *
     * This function adds a connectivity block to the geometric support for
     * \ref CWP_SUPPORT_MESH support type. Definition of element connectivity is :
     *
     *  - edge (\ref CWP_BLOCK_EDGE2) :
     *
     *   \code
     *       1 x-------x 2
     *   \endcode
     *
     *  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
     *
     *   \code 
     *       1 x-------x 3
     *          \     /
     *           \   /
     *            \ /
     *             x 2
     *   \endcode
     *
     *  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
     *
     *   \code
     *          4 x-------x 3
     *           /       /
     *          /       /
     *       1 x-------x2
     *   \endcode
     *
     *   - tetrahedron (\ref CWP_BLOCK_CELL_TETRA4) : 
     *
     *   \code
     *             x 4
     *            /|\
     *           / | \
     *          /  |  \
     *       1 x- -|- -x 3
     *          \  |  /
     *           \ | /
     *            \|/
     *             x 2
     *   \endcode
     *
     *   - pyramid (\ref CWP_BLOCK_CELL_PYRAM5) :
     *
     *   \code
     *              5 x
     *               /|\
     *              //| \
     *             // |  \
     *          4 x/--|---x 3
     *           //   |  /
     *          //    | /
     *       1 x-------x 2
     *   \endcode
     *
     *  - prism (\ref CWP_BLOCK_CELL_PRISM6) :
     *
     *   \code
     *       4 x-------x 6
     *         |\     /|
     *         | \   / |
     *       1 x- \-/ -x 3
     *          \ 5x  /
     *           \ | /
     *            \|/
     *             x 2
     *   \endcode
     *
     *  -  hexaedron (\ref CWP_BLOCK_CELL_HEXA8) :
     *
     *   \code
     *          8 x-------x 7
     *           /|      /|
     *          / |     / |
     *       5 x-------x6 |
     *         | 4x----|--x 3
     *         | /     | /
     *         |/      |/
     *       1 x-------x 2
     *   \endcode
     *
     * \param [in]  i_part      Current partition
     * \param [in]  block_type  Block type
     * \param [in]  n_elts      Number of elements
     * \param [in]  connec      Connectivity (size = n_vertex_elt * \ref n_elts)          
     * \param [in]  parent_num  Pointer to parent element number (or NULL)
     *
     */

    inline void 
    supportBlockAdd
    (
     const int           i_part,
     const CWP_Block_t block_type,
     const int           n_elts,
     const int           connec[],
     const CWP_g_num_t  parent_num[]
    );

    /**
     * \brief Adding a polygon connectivity block to the geometric support
     *
     * This function adds a polygon connectivity block to the geometric support for
     * \ref CWP_SUPPORT_MESH support type.
     *
     * \param [in]  i_part      Current partition
     * \param [in]  n_elts      Number of elements
     * \param [in]  connec_idx  Connectivity index (connec_id[0] = 0 and 
     *                          size = \ref n_elts + 1)          
     * \param [in]  connec      Connectivity (size = connec_id[n_elts] * \ref n_elts)          
     * \param [in]  parent_num  Pointer to parent element number (or NULL)
     *
     */

    inline void 
    supportFPolyBlockAdd
    (
     const int            i_part,
     const CWP_Block_t  block_type,
     const int            n_elts,
     const int            connec[],
     const CWP_g_num_t   parent_num[]
    );

    /**
     * \brief Adding a polyhedron connectivity block to the geometric support
     *
     * This function add a connectivity block to the geometric support if support
     * type is only \ref CWP_SUPPORT_MESH. Definition of element connectivity is :
     *
     * \param [in]  i_part            Current partition
     * \param [in]  n_elts            Number of elements
     * \param [in]  cell_face_idx     Polyhedron to face index 
     *                                (src_poly_cell_face_idx[0] = 0 and
     *                                 size = n_elts + 1)
     * \param [in]  cell_face_connec  Polyhedron to face connectivity 
     *                                (size = cell_face_idx[n_elts])
     * \param [in]  n_faces           Number of faces      
     * \param [in]  face_vtx_idx      Polyhedron face to vertex index 
     *                                (face_vertex_idx[0] = 0 and
     *                                 size_idx = max(cell_face_connec) + 1)
     * \param [in]  face_vtx_connec   Polyhedron face to vertex connectivity
     *                                (size = face_vertex_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)
     *
     */

    inline void 
    supportCPolyBlockAdd
    (
     const int           i_part,
     const int           n_elts,
     const int           cell_face_idx[],
     const int           cell_face[],
     const int           n_faces,
     const int           face_vtx_idx[],
     const int           face_vtx[],
     const CWP_g_num_t  parent_num[]
    );

    /**
     * \brief Geometric support removal                                  
     *
     * This function delete the geometric support  
     *
     */

    inline void 
    supportDel
    (
    );

    /**
     * \brief Map a fvm nodal as support mesh                                  
     *
     * This function  map a fvm nodal as support mesh
     *
     * \param [in]  i_part            Current partition
     * \param [in]  fvmc_nodal        fvm nodal mes     
     *
     */

    inline void 
    fvmcNodalShared
    (
     const int           i_part,
           fvmc_nodal_t *fvmc_nodal
    );

    /*----------------------------------------------------------------------------*
     * Methods about field                                                        *
     *----------------------------------------------------------------------------*/

    /**
     *
     * \brief Create a new field
     * 
     * \param [in]  field_id       Field id
     * \param [in]  data_type      Data type          
     * \param [in]  storage        Storage type          
     * \param [in]  n_component    Number of componenent
     * \param [in]  nature         Nature
     * \param [in]  visu_status    Visualization status
     * 
     */

    inline void
    fieldCreate
    (
     const string                &field_id,
     const CWP_Type_t     data_type,
     const CWP_Field_storage_t  storage,
     const int                    n_component,
     const CWP_Field_value_t   nature,
     const CWP_Status_t         visu_status
    );

    /**
     *
     * \brief Set data mapping
     * 
     * \param [in]  field_id       Field identifier
     * \param [in]  data           Storage array (Mapping)
     * 
     */

    inline void
    fieldMappingset
    (
     const string &field_id,
           double data[]   
    );

    /**
     *
     * \brief Get nunmber of field components
     * 
     * \param [in]   field_id       Field identifier
     *
     * \return                      number of field components
     * 
     */

    inline int
    fieldNComponentGet
    (
     const string &field_id
    );

    /**
     *
     * \brief Get field nature
     * 
     * \param [in]   field_id       Field identifier
     *
     * \return                      Field data type
     * 
     */

    inline CWP_Field_value_t
    fieldNatureGet
    (
     const string &field_id
    );
    
    /**
     *
     * \brief Get field data type
     * 
     * \param [in]   field_id       Field identifier
     *
     * \return                      Field data type
     * 
     */

    inline CWP_Type_t
    fieldTypeGet
    (
     const string &field_id
    );
    
    /**
     *
     * \brief Get field storage type
     * 
     * \param [in]   field_id       Field identifier
     * 
     */

    inline CWP_Field_storage_t
    fieldStorageGet
    (
     const string &field_id
    );

    /**
     *
     * \brief Removing a field
     * 
     * \param [in]  field_id       Field identifier
     * 
     */

    inline void
    fieldDel
    (
     const string &field_id
    );

    /*----------------------------------------------------------------------------*
     * Methods about exchange                                                     *
     *----------------------------------------------------------------------------*/

   /**
     * \brief Exchange data field with the coupled code with blocking 
     *        communications.
     *
     * This function exchanges interpolated fields between coupled codes. 
     * 
     * \warning  The size of \ref tgt_field_id size is n_computed_tgt. 
     *           If \f$ n\_uncomputed\_tgt \ne n\_tgt\_pts \f$,
     *           user himself must set values for uncomputed target points.
     *
     * \param [in]  src_id                    Source field (NULL -> no sending)
     * \param [in]  tgt_id                    Target field (NULL -> no receiving)
     * \param [in]  ptFortranInterpolationFct Fortran user interpolation (or NULL)
     * \param [out] n_uncomputed_tgt          Number of uncomputed target
     *
     * \return                                Exchange status
     *
     */

    inline CWP_Err_t 
    sendRecv
    (
     string &src_field_id,
     string &tgt_field_id,
     void   *ptFortranInterpolationFct,
     int    *n_uncomputed_tgt
    );

    /**
     *
     * \brief Sending of data field to the coupled code with nonblocking 
     *        communications.
     *
     * This function sends interpolated field to the coupled code. 
     * 
     * \param [in]  src_id                    Source field    
     * \param [in]  ptFortranInterpolationFct Fortran user interpolation (or NULL)
     *
     * \param [out] request                   Request to call by 
     *                                        \ref CWP_Wait_issend to wait 
     *                                        the end of exchange
     *
     */

    inline void 
    issend
    (
     string &src_id,
     void   *ptFortranInterpolationFct,
     int    *request
    );

    /**
     *
     * \brief Waiting of the end of exchange related to \ref request.
     *
     * This function waits the end of exchange related to \ref request
     * from \ref CWP_Issend
     * 
     * \param [in] request    Request to wait the end of exchange
     *
     */

    void 
    waitIssend
    (
     int request
    );

    /**
     *
     * \brief Receiving of Data field from the coupled code with nonblocking 
     *        communications.
     *
     * This function receives interpolated field from the coupled code 
     * 
     * \param [in]  tgt       Target field   
     *
     * \param [out] request   Request to call by \ref CWP_Wait_irecv  
     *                        to wait the end of exchange
     *
     */

    void 
    irecv
    (
     string &tgt,
     int     *request
    );

    /**
     *
     * \brief Waiting of the end of exchange related to \ref request.
     *
     * This function waits the end of exchange related to \ref request 
     * from \ref CWP_Irecv
     * 
     * \param [in] request    Request to wait the end of exchange
     *
     */

    void 
    waitIrecv
    (
     int request
    );    

    /*----------------------------------------------------------------------------*
     * methods about user interpolation                                           *
     *----------------------------------------------------------------------------*/

    /**
     *
     * \brief Setting of an user interpolation from location.
     *
     * This function takes into account an user interpolation function written with
     * \ref void (*CWP_Interp_from_location_t) interface.
     * 
     * \param [in] fct        Function
     *
     */

    void 
    interpFromLocSet
    (
     CWP_Interp_from_location_t fct
    );

    /**
     *
     * \brief Setting of a FORTRAN user interpolation from location.
     *
     * This function takes into account an user interpolation function written
     * in FORTRAN .
     * 
     * \param [in] fct        Function
     *
     */

    void 
    interpFromLocSetF
    (
     void       *fct
    );

    /**
     *
     * \brief Setting of an user interpolation from intersection.
     *
     * This function takes into account an user interpolation function written with
     * \ref void (*CWP_Interp_from_intersec_t) interface.
     * 
     * \param [in] fct        Function
     *
     */

    void 
    interpFromInterSet   
    (
     CWP_Interp_from_intersec_t fct
    );

    /**
     *
     * \brief Setting of a FORTRAN user interpolation from intersection.
     *
     * This function takes into account an user interpolation function written
     * in FORTRAN .
     * 
     * \param [in] fct        Function
     *
     */

    void 
    interpFromInterSetF      
    (
     void       *fct
    );

    /**
     *
     * \brief Setting of an user interpolation from closest points
     *
     * This function takes into account an user interpolation function written with
     * \ref void (*CWP_Interp_from_closest_pts_t) interface.
     * 
     * \param [in] fct        Function
     *
     */

    void 
    interpFromClosestSet    
    (
     CWP_Interp_from_closest_pts_t fct
    );

    /**
     *
     * \brief Setting of a FORTRAN user interpolation from closest points
     *
     * This function takes into account an user interpolation function written
     * in FORTRAN .
     * 
     * \param [in] fct        Function
     *
     */

    void 
    interpFromClosestSetF      
    (
     void *fct
    );


    /**
     *
     * \brief Return communication type
     *
     * 
     * \param [in] fct        Function
     *
     */

    inline CWP_Comm_t 
    commTypeGet      
    (
    );
    
    
    
    
  private:

    Coupling();

  private:
    const string                     _cplId;                 /*!< Coupling identifier */
    CWP_Comm_t                       _commType;              /*!< Communication type */ 
          Communication              &_communication;         /*!< Communication */ 
    const CodeProperties             &_localCodeProperties;   /*!< Local code properties */
    const CodeProperties             &_coupledCodeProperties; /*!< Coupled code properties */
    //    Geometry                   &_geometry;              /*!< Geometric algorithm */
    //    Support                    &_support;               /*!< Geometric support */
    const CWP_Freq_t                _recvFreqType;          /*!< Receiving frequency type */
    //        Visualization              *_visu;                  /*!< Visualization */
          double                      _recvFreq;              /*!< Receiving frequency */
          double                      _recvNextTime;          /*!< Next receiving time */
    //map < string, Field<double> * >  &_fields;                /*!< Fields storage */
          CouplingDB                 &_cplDB;                  /*!< Coupling Data base */
  }; 

  
  

  /**
   *
   * \brief Return communication type
   *
   * 
   * \param [in] fct        Function
   *
   */

  CWP_Comm_t 
  Coupling::commTypeGet      
  (
  )
  {
    return _commType;
  }

}

#endif //__COUPLING_PROPERTIES_H__
