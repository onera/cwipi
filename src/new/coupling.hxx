#ifndef __COUPLING_H__
#define __COUPLING_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2012-2017  ONERA

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

#include "cwp.h"
#include "pdm_printf.h"
#include "communication.hxx"
#include "couplingDB.hxx"
#include "couplingDB_i.hxx"
#include "mesh.hxx"

#include "spatialInterp.hxx"
#include "visualization.hxx"
#include "field.hxx"


/**
 * \cond
 */
using namespace std;

namespace cwipi {

  class CodeProperties;
  class SpatialInterp;
  class Mesh;
  class Field;
  class Visu;
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
     * \param [in]  cplId                        Coupling identifier
     * \param [in]  commType                     Communication type
     * \param [in]  localCodeProperties          Local code properties
     * \param [in]  coupledCodeProperties        Coupled code properties
     * \param [in]  spatialInterpAlgo                     SpatialInterp algorithm
     * \param [in]  nPart                        Number of interface partitions
     * \param [in]  movingStatus                 Mesh moving status
     * \param [in]  recvFreqType                 Type of receiving frequency
     * \param [in]  cplDB                        Coupling data base where it coupling is stored
     *
     */

    Coupling (
     const string                &cplId,
     const CWP_Comm_t             commType,
     const CodeProperties        &localCodeProperties,
     const CodeProperties        &coupledCodeProperties,
     const CWP_Interface_t        entities_dim,
     const CWP_Spatial_interp_t   spatialInterpAlgo,
     const int                    nPart,
     const CWP_Dynamic_mesh_t     movingStatus,
     const CWP_Time_exch_t        recvFreqType,
     CouplingDB                  &cplDB
    );

    /**
     * \brief Destructor.
     *
     */

    virtual ~Coupling();

    /*----------------------------------------------------------------------------*
     * Methods about exchange frequency                                           *
     *----------------------------------------------------------------------------*/

    /**
     * \brief Setting receiving frequency.
     *
     * This function set receiving frequency. It must be used when
     * the type of receiving frequency is \ref CWP_TIME_EXCH_N_TIME_STEP
     *
     * \param [in]  n_step     Frequency in steps number
     *
     */

    inline void
    recvFreqSet (
      int n_step
    );

    /**
     * \brief Setting the next receiving time.
     *
     * This function set the next receiving time. It must be used when
     * the type of receiving frequency is \ref CWP_TIME_EXCH_ASYNCHRONOUS
     *
     * \param [in]  next_time     Next receiving time
     *
     */

    inline void
    recvNextTimeSet (
      double next_time
    );

    /*----------------------------------------------------------------------------*
     * Methods about spatial interpolation                                        *
     *----------------------------------------------------------------------------*/

    /**
     * \brief Computation spatial interpolation weights
     *
     * This function compute spatial interpolation weights
     *
     * \param [out] n_uncomputed_tgt    Number of uncomputed target
     *
     */

    void
    spatialInterpWeightsCompute ();

    /**
     * \brief Set the spatial interpolation properties.
     *
     *
     * \param [in]       fmt       Format with the syntax : "prop1, prop2, ..."
     * \param [in,out]   pa        List of properties values
     *
     */

    inline void
    spatialInterpPropertiesSet (
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
     * \param [in]  format           Output format to visualize exchanged fieldsDouble
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

    void
    visuSet (
      const int               freq,
      const CWP_Visu_format_t format,
      const char             *format_option
    );


    /**
     *
     * \brief Return the Visu object associated to this coupling
     *
     * \return Visu object pointer
     *
     */

    inline Visu* 
    visuGet ();

    /*----------------------------------------------------------------------------*
     * Methods  about mesh                                                     *
     *----------------------------------------------------------------------------*/

    /**
     * \brief Setting vertices
     *
     * This method set partition vertices
     *
     * \param [in]  i_part      Current partition
     * \param [in]  n_pts       Number of points
     * \param [in]  coord       Coordinates (size = 3 * n_pts)
     * \param [in]  global_num  Pointer to global element number (or NULL)
     *
     */

    inline void
    meshVtcsSet (
      const int          i_part,
      const int          n_pts,
      double             coord[],
      CWP_g_num_t        global_num[]
    );


   /**
    * \brief Add a block to the interface mesh.
    *
    *
    * \param [in]  block_type       Block type
    *
    * \return block identifier
    */

    inline int
    meshBlockAdd (
      const CWP_Block_t     block_type
    );


    /**
     * \brief Set a standard block to the interface mesh
     *
     * This function adds a connectivity block to the geometric support.
     *
     *  Definition of element connectivity is :
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
     * \param [in]  block_id    Block identifier
     * \param [in]  n_elts      Number of elements
     * \param [in]  connec      Connectivity (size = n_vertex_elt * n_elts)
     * \param [in]  global_num  Pointer to global element numbering (or NULL)
     *
     */

    inline void
    meshStdBlockSet(
      const int           i_part,
      const int           block_id,
      const int           n_elts,
      int                 connec[],
      CWP_g_num_t       global_num[]
    );


    /**
     * \brief Set a generic high order block to the interface mesh
     *
     *
     * \param [in]  i_part      Partition identifier
     * \param [in]  block_id    Block identifier
     * \param [in]  n_elts      Number of elements
     * \param [in]  order       SpatialInterp order
     * \param [in]  connec      Connectivity (size = n_vertex_elt * n_elts)
     * \param [in]  global_num  Pointer to global element number (or NULL)
     *
     */

    inline void
    meshHighOrderBlockSet (
      const int           i_part,
      const int           block_id,
      const int           n_elts,
      const int           order,
      int                 connec[],
      CWP_g_num_t         global_num[]
    );

    /**
     * \brief Set the connectivity of a polygon block in a mesh interface partition.
     *
     *
     * \param [in]  i_part      Current partition
     * \param [in]  block_id    Block identifier
     * \param [in]  n_elts      Number of elements
     * \param [in]  connec_idx  Connectivity index (connec_id[0] = 0 and
     *                          size = n_elts + 1)
     * \param [in]  connec      Connectivity (size = connec_id[n_elts] * n_elts)
     * \param [in]  global_num  Pointer to global element number (or NULL)
     *
     */

    inline void
    meshFPolyBlockSet (
      const int            i_part,
      const int            block_id,
      const int            n_elts,
      int                  connec_idx[],
      int                  connec[],
      CWP_g_num_t          global_num[]
    );

    /**
     * \brief Set the connectivity of a polyhedron block in a mesh interface partition.
     *
     * Definition of element connectivity is :
     *
     * \param [in]  i_part            Current partition
     * \param [in]  block_id          Block identifier
     * \param [in]  n_elts            Number of elements
     * \param [in]  connec_cells_idx  Polyhedron to face index
     *                                (src_poly_cell_face_idx[0] = 0 and
     *                                 size = n_elts + 1)
     * \param [in]  connec_cells      Polyhedron to face connectivity
     *                                (size = cell_face_idx[n_elts])
     * \param [in]  n_faces           Number of faces
     * \param [in]  connec_faces_idx  Polyhedron face to vertex index
     *                                (face_vertex_idx[0] = 0 and
     *                                size_idx = max(cell_face_connec) + 1)
     * \param [in]  connec_faces      Polyhedron face to vertex connectivity
     *                                (size = face_vertex_idx[size_idx - 1])
     * \param [in]  global_num        Pointer to global element number (or NULL)
     *
     */

    inline void
    meshCPolyBlockSet (
      const int           i_part,
      const int           block_id,
      const int           n_elts,
      const int           n_faces,
      int                 connec_faces_idx[],
      int                 connec_faces[],
      int                 connec_cells_idx[],
      int                 connec_cells[],
      CWP_g_num_t         global_num[]
    );


    /**
     * \brief Adding a polyhedron block to the mesh from
     * a face-to-cell connectivity and a vertices-to-faces connectivity.
     *
     * This function adds a polyhedron 3D block to the mesh from
     * a face-to-cell connectivity and a vertices-to-faces connectivity.
     *
     * \param [in]  i_part            Current partition
     * \param [in]  n_cells           Number of elements
     * \param [in]  cell_face_idx     Polyhedron to face index
     *                                (src_poly_cell_face_idx[0] = 0 and
     *                                 size = n_elts + 1)
     * \param [in]  cell_face         Polyhedron to face connectivity
     *                                (size = cell_face_idx[n_elts])
     * \param [in]  n_faces           Number of faces
     * \param [in]  face_vtx_idx      Polyhedron vertices to faces index
     *                                (face_vtx_idx[0] = 0 and
     *                                 size_idx = max(face_vtx) + 1)
     * \param [in]  face_vtx          Polyhedron vertices to faces connectivity
     *                                (size = face_vtx_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)
     *
     */

    inline void
    meshFromCellFaceSet(
      const int   i_part,
      const int   n_cells,
      int         cell_face_idx[],
      int         cell_face[],
      int         n_faces,
      int         face_vtx_idx[],
      int         face_vtx[],
      CWP_g_num_t parent_num[]
    );


    /**
     * \brief Adding a polygon 2D block to the mesh from
     * a vertices-to-faces connectivity and a edge-to-face connectivity.
     *
     * This function add a polygon 2D block to the mesh from
     * a vertices-to-faces connectivity and a edge-to-face connectivity.
     *
     * \param [in]  i_part            Current partition
     * \param [in]  n_faces           Number of faces
     * \param [in]  face_edge_idx     Polygon vertices to faces index
     *                                (face_edge_idx[0] = 0 and
     *                                size_idx = max(face_edge) + 1)
     * \param [in]  face_edge         Polyhegon vertices to face connectivity
     *                                (size = face_edge_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)
     * \param [in]  n_edges           Number of edges
     * \param [in]  edge_vtx_idx      Vertices to edges connectivity index
     *                                (edge_vtx_idx[0] = 0 and
     *                                size_idx = max(edge_vtx) + 1)
     * \param [in]  edge_vtx          Polygon vertices to edges connectivity
     *                                (size = edge_vtx_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)
     *
     */

    inline void
    meshFromFacesEdgeSet(
      const int   i_part,
      const int   n_faces,
      int         face_edge_idx[],
      int         face_edge[],
      const int   n_edges,
      int         edge_vtx_idx[],
      int         edge_vtx[],
      CWP_g_num_t parent_num[]
    );


    /**
     * \brief SpatialInterp mesh removal
     *
     * This function delete the  mesh
     *
     */

    inline void
    meshDel();


    /**
     *
     * \brief Finalize mesh description (Computation of entities global numbers if not given by the user)
     *
     */

    inline void 
    meshFinalize();

    /**
     *
     * \brief Return the mesh object associated to this coupling
     *
     * \return Mesh object pointer
     *
     */

    inline Mesh* 
    meshGet();


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
     * \param [in]  exch_type      Exchange type
     * \param [in]  visu_status    Visualization status
     *
     */

    void 
    fieldCreate (
     const string               &field_id,
     const CWP_Type_t           data_type,
     const CWP_Field_storage_t  storage,
     const int                  n_component,
     const CWP_Dof_location_t    nature,
     const CWP_Field_exch_t     exch_type,
     const CWP_Status_t         visu_status
    );


    /**
     * \brief Return if a field identifier exists
     *
     * \param [in]  field_id         Field identifier
     *
     * \return status
     */

    bool
    fieldIs (
     const string &field_id
    );


   /**
    *
    * \brief Set Field data
    *
    * \param [in]  field_id       Field identifier
    * \param [in]  data           Storage array (mapping)
    *
    */

    void fieldDataSet (
      const std::string &field_id,
      int i_part,
      const CWP_Field_map_t   map_type,
      void *data
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
    fieldNComponentGet (
      const string &field_id
    );


    /**
     *
     * \brief Get location of the degrees of freedom
     *
     * \param [in]   field_id       Field identifier
     *
     * \return                      Field data type
     *
     */

    inline CWP_Dof_location_t
    fieldDofLOcationGet (
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
    fieldStorageGet (
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
    fieldDel (
     const string &field_id
    );

    /*----------------------------------------------------------------------------*
     * Methods about exchange                                                     *
     *----------------------------------------------------------------------------*/

    /**
     * \brief data exchange <b>(Not implemented yet)</b>
     *
     * Exchange depending on exchange frequency
     *
     */

    void
    exchange ();


   /**
     * \brief Exchange data field with the coupled code with blocking
     *        communications.
     *
     * This function exchanges interpolated fieldsDouble between coupled codes.
     *
     * \warning  The size of tgt_field_id size is n_computed_tgt.
     *           If \f$ n\_uncomputed\_tgt \ne n\_tgt\_pts \f$,
     *           user himself must set values for uncomputed target points.
     *
     * \param [in]  src_field_id              Source field (NULL -> no sending)
     * \param [in]  tgt_field_id              Target field (NULL -> no receiving)
     * \param [in]  ptFortranInterpolationFct Fortran user interpolation (or NULL)
     * \param [out] n_uncomputed_tgt          Number of uncomputed target
     *
     * \return                                Exchange status
     *
     */

    void
    sendrecv (
      const string &field_id
    );

    /**
     *
     * \brief Sending of data field to the coupled code with nonblocking
     *        communications.
     *
     * This function sends interpolated field to the coupled code.
     *
     * \param [in]  src_id                    Source field
     *
     */

    void
    issend (
      const string &src_field_id
    );


    /**
     *
     * \brief Waiting of the end of exchange related to request.
     *
     * This function waits the end of exchange related to request
     * from \ref CWP_Issend
     *
     * \param [in] src_id                    Source field
     *
     */

    void
    waitIssend (
      const string &src_field_id
    );

    /**
     *
     * \brief Receiving of Data field from the coupled code with nonblocking
     *        communications.
     *
     * This function receives interpolated field from the coupled code
     *
     * \param [in]  receving_field_id       Target field ID
     *
     *
     */

    void
    irecv (
      const string &receving_field_id
    );

    /**
     *
     * \brief Waiting of the end of exchange related to request.
     *
     * This function waits the end of exchange related to request
     * from \ref CWP_Irecv
     *
     * \param [in]  receving_field_id       Target field ID
     *
     */

    void
    waitIrecv (
      const string &receving_field_id
    );


    /**
     *
     * \brief Return the number of uncomputed targets
     *
     * \return                Number of uncomputed targets
     */

    inline int
    nUncomputedTargetsGet(
      const string &field_id,
      const int  i_part
    ) const;

    /**
     *
     * \brief Return uncomputed targets
     *
     * \return                Uncomputed targets
     */

    inline const int *
    uncomputedTargetsGet (
      const string &field_id,
      const int  i_part
    ) const;

    /**
     *
     * \brief Return the number of computed targets
     *
     * \return                Number of computed targets
     */

    inline int
    nComputedTargetsGet (
      const string &field_id,
      const int  i_part
    ) const;

    /**
     *
     * \brief Return computed targets
     *
     * \return                Computed targets
     */

    inline const int *
    computedTargetsGet (
      const string &field_id,
      const int  i_part
    ) const;

    /*----------------------------------------------------------------------------*
     * methods about user interpolation                                           *
     *----------------------------------------------------------------------------*/

    /**
     *
     * \brief Setting of an user interpolation from location.
     *
     * This function takes into account an user interpolation function written with
     *  void (* \ref CWP_Interp_from_location_t) interface.
     *
     * \param [in] fct        Function
     *
     */

    inline void
    interpFromLocSet ( 
      const string field_id,
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

    inline void
    interpFromLocSetF (
      void       *fct
    );

    /**
     *
     * \brief Setting of an user interpolation from intersection.
     *
     * This function takes into account an user interpolation function written with
     *  void (* \ref CWP_Interp_from_intersect_t) interface.
     *
     * \param [in] fct        Function
     *
     */

    inline void
    interpFromInterSet (
      CWP_Interp_from_intersect_t fct
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

    inline void
    interpFromInterSetF (
      void       *fct
    );

    /**
     *
     * \brief Setting of an user interpolation from closest points
     *
     * This function takes into account an user interpolation function written with
     *  void (* \ref CWP_Interp_from_closest_pts_t) interface.
     *
     * \param [in] fct        Function
     *
     */

    inline void
    interpFromClosestSet (
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

    inline void
    interpFromClosestSetF (
      void *fct
    );


    /*----------------------------------------------------------------------------*
     * methods about attributes                                                   *
     *----------------------------------------------------------------------------*/

    /**
     *
     * \brief Return communication type
     *
     * \return CWP_Comm_t Communication Type
     *
     */

    inline CWP_Comm_t
    commTypeGet();


    /**
     *
     * \brief Return the dimesnion of geometric entities of this coupling
     *
     * \return Entities dimension 
     *
     */

    inline CWP_Interface_t 
    entitiesDimGet();


    /**
     *
     * \brief Return the local code fields defined for this coupling
     *
     * \return Fields 
     *
     */

    inline std::map < string, Field * >* 
    fieldsGet();


    /**
     *
     * \brief Return the spatial interpolation objects according to direction of the exchange
     * 
     * \param [in] exchDirection     Direction of the exchange
     *
     * \return Local  spatial interpolation objects
     *
     */

    inline std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > ,SpatialInterp*>* 
    spatialInterpGet(
      CWP_Field_exch_t exchDirection
    );


    /**
     *
     * \brief Return the spatial interpolation according to parameters
     * 
     * \param [in] localLocation     Local location of degrees of freedom
     * \param [in] cplLocation       Coupled location of degrees of freedom
     * \param [in] exchDirection     Direction of the exchange
     *
     * \return Local code properties
     *
     */

    inline SpatialInterp* 
    spatialInterpGet (
      CWP_Dof_location_t localLocation, 
      CWP_Dof_location_t cplLocation, 
      CWP_Field_exch_t   exchDirection
    );

    /**
     *
     * \brief Return the local code properties
     *
     * \return Local code properties
     *
     */
    
    inline CodeProperties* 
    localCodePropertiesGet();
    

    /**
     *
     * \brief Return the coupled code properties
     *
     * \return Local code properties
     *
     */
    
    inline CodeProperties* 
    coupledCodePropertiesGet();

    /**
     *
     * \brief Return the communication way 
     *
     * \return Local code properties
     *
     */

    inline Communication* 
    communicationGet();

    /**
     *
     * \brief Return the couplings data base
     *
     * \return Coupling data base
     *  
     */

    inline CouplingDB*  
    couplingDBGet();
    

    /**
     *
     * \brief Return the coupling id
     *
     * \return id
     *  
     */
 
    inline string       
    IdGet();

    /**
     *
     * \brief Return sendSpatial map
     *
     * \return id
     *  
     */
 
    inline std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &    
    sendSpatialInterpGet();


    /**
     *
     * \brief Return recvSpatial map
     *
     * \return id
     *  
     */
 
    inline std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &    
    recvSpatialInterpGet();

    /*----------------------------------------------------------------------------*
     * methods about user target                                                  *
     *----------------------------------------------------------------------------*/


    /**
     * \brief Define user target points 
     * 
     * \param [in]  iPart     Current partition
     * \param [in]  nPts      Number of points
     * \param [in]  coords    Coordinates
     * \param [in]  gNum      Global Number (or NULL)
     *
     */

    inline void
    userTargetSet (
      const int         iPart,
      const int         nPts,
      const double      coords[],
      const CWP_g_num_t gNum[]
    );


    /**
     * \brief Return global number of user targets 
     * 
     * \param [in]  iPart     Current partition
     * 
     * 
     * \return Global Number 
     *
     */

    inline const CWP_g_num_t *
    userTargetGNumGet (
      const int         iPart
    );


    /**
     * \brief Return coords of user targets 
     * 
     * \param [in]  iPart     Current partition
     * 
     * 
     * \return Coordinates 
     *
     */

    inline const double *
    userTargetCoordsGet (
      const int         iPart
    ) const ;


    /**
     * \brief Return number of user targets 
     * 
     * \param [in]  iPart     Current partition
     * 
     * 
     * \return Number of user targets
     *
     */

    inline int
    userTargetNGet (
      const int         iPart
    );


  public:

    // A supprimer 

    CWP_g_num_t*
    globalNumGet(int id_block,int i_part);

  private:


    /**
     *
     * \brief Compute user target global number (if not given by user)
     *
     */

    void 
    userTargetGnumCompute();


    /**
     *
     * \brief unused default constructor
     *
     */

    Coupling();

  private:
    const string                            _cplId;                 /*!< Coupling identifier */
          CWP_Comm_t                        _commType;              /*!< Communication type */
          Communication                    &_communication;         /*!< Communication */
    const CodeProperties                   &_localCodeProperties;   /*!< Local code properties */
    const CodeProperties                   &_coupledCodeProperties; /*!< Coupled code properties */
    const CWP_Interface_t                  _entities_dim;           /*!< Mesh entities dimension */
          Mesh                             &_mesh;                  /*!< SpatialInterp mesh */
    const CWP_Time_exch_t                   _recvFreqType;          /*!< Receiving frequency type */
          Visu                             &_visu;                  /*!< Visualization */
          double                            _recvFreq;              /*!< Receiving frequency */
          double                            _recvNextTime;          /*!< Next receiving time */
          std::map < string, Field * >     &_fields;                /*!< Fields Data Base */
          CouplingDB                       &_cplDB;                 /*!< Coupling Data base */
          int*                              _iteration;             /*!< Interation ???? */
          CWP_Dynamic_mesh_t                _displacement;          /*!< Type of mesh displacement */ 
    const CWP_Spatial_interp_t              _spatialInterpAlgo;     /*!< Spatial intepolation algorithm */
    const int                               _nPart;                 /*!< Number of partitions */  

          int                              *_userTargetN;           /*!< Number of user targets on by partition (size number partitions of the mesh) */
    const CWP_g_num_t                     **_userTargetGnum;        /*!< Target global numbering by partition (size number partitions of the mesh) */
          CWP_g_num_t                     **_localUserTargetGnum;   /*!< Target global numbering by partition (used if _gnum_user_target is not setted by user) */
    const double                          **_userTargetCoord;       /*!< Target coordinates by partition (size number partitions of the mesh) */

    std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &_spatial_interp_send; /*!< local sent Spatial interpolation objects 
                                                                                                                  to associate with receive distant spatial interpolatiol */
    std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &_spatial_interp_recv; /*!< local receive Spatial interpolation objects 
                                                                                                                  to associate with sent distant spatial interpolatiol */
  };
}

/**
 * \endcond
 */

#endif //__COUPLING_H__
