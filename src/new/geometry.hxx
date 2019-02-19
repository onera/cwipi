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

namespace cwipi {

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

    /**
     * \brief Set meshInterf
     *
     * This function defines the geometry meshInterf geometry.
     */

    void 
    meshInterfSet
    (const int        i_part,
     const int        n_pts,
     double           coords[],
     CWP_g_num_t      parent_num[]);

   
    /**
     * \brief Finalize the mesh setting after blocks and coordinates additions
     *
     */    

    void 
    meshInterfEndSet
    (
    );
     
     

    /**
     * \brief Adding a connectivity block to the geometric meshInterf
     *
     * This function adds a connectivity block to the geometric meshInterf.
     *   Definition of element connectivity is :
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
     * \param [in]  connec      Connectivity (size = n_vertex_elt *n_elts)          
     * \param [in]  parent_num  Pointer to parent element number (or NULL)
     *
     */

    void 
    meshInterfBlockAdd
    (const int            i_part,
     const CWP_Block_t    block_type,
     const int            n_elts,
     int                  connec[],
     CWP_g_num_t          parent_num[]);

    /**
     * \brief Adding a polygon connectivity block to the geometric meshInterf
     *
     * This function adds a polygon connectivity block to the geometric meshInterf.
     *
     * \param [in]  i_part      Current partition
     * \param [in]  block_type  Type of the block elements
     * \param [in]  n_elts      Number of elements
     * \param [in]  connec_idx  Connectivity index (connec_id[0] = 0 and 
     *                          size =n_elts + 1)          
     * \param [in]  connec      Connectivity (size = connec_id[n_elts] *n_elts)        
     * \param [in]  parent_num  Pointer to parent element number (or NULL)
     *
     */

    void 
    meshInterfFPolyBlockAdd
    (const int          i_part,
     const CWP_Block_t  block_type,
     const int          n_elts,
     int                connec_idx[],
     int                connec[],
     CWP_g_num_t        parent_num[]);

    /**
     *
     * \brief Adding a polyhedron connectivity block to the geometric meshInterf
     *
     * This function add a connectivity block to the geometric meshInterf.
     *
     * Definition of element connectivity is :
     *
     * \param [in]  i_part            Current partition
     * \param [in]  n_elts            Number of elements
     * \param [in]  cell_face_idx     Polyhedron to face index 
     *                                (src_poly_cell_face_idx[0] = 0 and
     *                                 size = n_elts + 1)
     * \param [in]  cell_face          Polyhedron to face connectivity 
     *                                (size = cell_face_idx[n_elts])
     * \param [in]  n_faces           Number of faces      
     * \param [in]  face_vtx_idx      Polyhedron face to vertex index 
     *                                (face_vertex_idx[0] = 0 and
     *                                 size_idx = max(cell_face_connec) + 1)
     * \param [in]  face_vtx           Polyhedron face to vertex connectivity
     *                                (size = face_vertex_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)
     *
     */

    void 
    meshInterfCPolyBlockAdd
    (const int     i_part,
     const int     n_elts,
     int           cell_face_idx[],
     int           cell_face[],
     const int     n_faces,
     int           face_vtx_idx[],
     int           face_vtx[],
     CWP_g_num_t   parent_num[]);

    /**
     *
     * \brief Map a fvm nodal as meshInterf mesh                                  
     *
     * This function  map a fvm nodal as meshInterf mesh
     *
     * \param [in]  i_part            Current partition
     * \param [in]  fvmc_nodal        fvm nodal mes     
     *
     */

    void 
    fvmcNodalShared
    (const int      i_part,
     fvmc_nodal_t  *fvmc_nodal);

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
     * \param [in]  src                       Source field    
     * \param [in]  ptFortranInterpolationFct Fortran user interpolation (or NULL)
     *
     * \param [out] request                   Request to call by 
     *                                        \ref CWP_Wait_issend to wait 
     *                                        the end of exchange
     *
     */

    void 
    issend
    (Field <double> &src,
     void           *ptFortranInterpolationFct,
     int            *request);

    /**
     *
     * \brief Waiting of the end of exchange related to request.
     *
     * This function waits the end of exchange related to request
     * from \ref CWP_Issend
     * 
     * \param [in] request    Request to wait the end of exchange
     *
     */

    void 
    waitIssend
    (int request);

    /**
     *
     * \brief Receiving of Data field from the coupled application with nonblocking 
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
    (Field <double>  &tgt,
     int             *request);

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
    (int request);    

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

  private:
    
    Geometry &operator=(const Geometry &other);  /*!< Assigment operator not available */
    Geometry (const Geometry& other);            /*!< Copy constructor not available */

  private:
    std::vector<double>         *_tmpVertexField;  /*!< Evite une allocation a chaque extrapolation */
    std::vector<double>         *_tmpDistantField; /*!< TODO: Fusionner _tmpDistantField utiliser pour exchange
                                                        et les comm asynchrones */
    Mesh                        *_meshInterf;      /*!< Interface Mesh */
   // std::map<int, std::vector<double> * > &_tmpDistantFieldsIssend; //TODO: temporaire : A revoir lors
                                                                    // de la restructuration
 //   std::map<int, const double * >        &_tmpLocalFieldsIrecv;
 //   std::map<int, const char * >          &_tmpExchangeNameIrecv;
 //   std::map<int, int >                   &_tmpStrideIrecv;
//    std::map<int, int >                   &_tmpTimeStepIrecv;
 //   std::map<int, double >                &_tmpTimeValueIrecv;
 //   std::map<int, const char * >          &_tmpFieldNameIrecv;

  };

}

#endif //__GEOMETRY_H__
