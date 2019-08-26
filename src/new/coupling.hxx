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

#include <fvmc_nodal.h>

#include "cwp.h"
#include "bftc_printf.h"
#include "communication.hxx"
#include "couplingDB.hxx"
#include "couplingDB_i.hxx"
#include "mesh.hxx"

#include "geometry.hxx"
#include "visualization.hxx"
#include "field.hxx"

using namespace std;

namespace cwipi {

  class CodeProperties;
  class Geometry;
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
     * \param [in]  geomAlgo                     Geometric algorithm
     * \param [in]  nPart                        Number of interface partitions 
     * \param [in]  movingStatus                 Mesh moving status
     * \param [in]  recvFreqType                 Type of receiving frequency
     * \param [in]  cplDB                        Coupling data base where it coupling is stored
     *
     */

    Coupling
    (
     const string                &cplId,
     const CWP_Comm_t             commType,
     const CodeProperties        &localCodeProperties,
     const CodeProperties        &coupledCodeProperties,
     const CWP_Geom_t            geomAlgo,
     const int                    nPart,
     const CWP_Displacement_t     movingStatus,
     const CWP_Freq_t             recvFreqType,
     CouplingDB                  &cplDB 
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

    int 
    nUncomputedTargetsGet
    (
      const CWP_Field_value_t geometryLocation,
      const int  i_part
    ) ;

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
     * the type of receiving frequency is \ref CWP_FREQ_RELATED_N_TIME_STEP
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
     * the type of receiving frequency is \ref CWP_FREQ_ASYNCHRONOUS
     *
     * \param [in]  next_time     Next receiving time
     *
     */

    void
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
     CWP_Field_value_t geometryLocation,
     CWP_Field_exch_t exchange_type
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
    visuSet
    (
     const int               freq,
     const CWP_Visu_format_t format,
     const char             *format_option
    );

    /*----------------------------------------------------------------------------*
     * Methods about User target points                                           *
     *----------------------------------------------------------------------------*/

    /**
     * \brief Setting user target points
     *
     * This function must be called if the nature of receiving fieldsDouble 
     * is \ref CWP_FIELD_VALUE_USER_TO_NODE
     *
     * \param [in]  n_pts   Number of points
     * \param [in]  coord   Coordinates (size = 3 * n_pts)          
     *
     */

    void 
    userTgtPtsSet
    (
     const int i_part,
     const int n_pts,
     double    coord[]
    );

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

    void 
    meshVtcsSet
    (
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
    
    int 
    meshBlockAdd
    (
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

    void 
    meshStdBlockSet
    (
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
     * \param [in]  order       Geometric order
     * \param [in]  connec      Connectivity (size = n_vertex_elt * n_elts)          
     * \param [in]  global_num  Pointer to global element number (or NULL)
     *
     */
    
   /* void 
    meshHighOrderBlockSet
    (
     const int           i_part,
     const int           block_id,
     const int           n_elts,
     const int           order,
     int                 connec[],
     CWP_g_num_t         global_num[]);

*/
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
     
    void 
    meshFPolyBlockSet
    (
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

    void 
    meshCPolyBlockSet
    (
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
     * \brief Adding a polyhedron block to the geometric mesh from
     * a face-to-cell connectivity and a vertices-to-faces connectivity.
     *
     * This function add a polyhedron 3D block to the geometric mesh from
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

    void
    meshFromCellFaceSet(const int   i_part,
                        const int   n_cells,
                        int         cell_face_idx[],
                        int         cell_face[],
                        int         n_faces,
                        int         face_vtx_idx[],
                        int         face_vtx[],
                        CWP_g_num_t parent_num[]); 




    /**
     * \brief Adding a polygon 2D block to the geometric mesh from
     * a vertices-to-faces connectivity and a edge-to-face connectivity.
     *
     * This function add a polygon 2D block to the geometric mesh from
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
     
    void
    meshFromFacesEdgeSet(const int   i_part,
                         const int   n_faces,
                         int         face_edge_idx[],
                         int         face_edge[],
                         const int   n_edges,
                         int         edge_vtx_idx[],
                         int         edge_vtx[],
                         CWP_g_num_t parent_num[]); 
     
     
    /**
     * \brief Geometric mesh removal                                  
     *
     * This function delete the geometric mesh  
     *
     */
     
    void 
    meshDel
    (
    );

    /**
     * \brief Map a fvm nodal as mesh mesh                                  
     *
     * This function  map a fvm nodal as mesh mesh
     *
     * \param [in]  i_part            Current partition
     * \param [in]  fvmc_nodal        fvm nodal mesh 
     *
     */

   /* void 
    fvmcNodalShared
    (
     const int           i_part,
     fvmc_nodal_t       *fvmc_nodal
    );*/

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
    fieldCreate
    (
     const string               &field_id,
     const CWP_Type_t           data_type,
     const CWP_Field_storage_t  storage,
     const int                  n_component,
     const CWP_Field_value_t    nature,
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
    fieldIs
    (
     const string &field_id
    );
    
    

    /**
     *
     * \brief Set data mapping
     * 
     * \param [in]  field_id       Field identifier
     * \param [in]  data           Storage array (Mapping)
     * 
     */

  void fieldDataSet
  (
    const std::string &field_id,
    int i_part,
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

    int
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

    CWP_Field_value_t
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
  
    CWP_Field_value_t
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

    CWP_Field_storage_t
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

    void
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
     *
     */

    void 
    issend
    (
     string &src_field_id
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
    waitIssend
    (
     string &src_field_id
    );



    CWP_g_num_t* 
    globalNumGet(int id_block,int i_part);

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
    irecv
    (
     string &receving_field_id
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
    waitIrecv
    (
     string &receving_field_id
    );    

    /*----------------------------------------------------------------------------*
     * methods about user interpolation                                           *
     *----------------------------------------------------------------------------*/

    /**
     *
     * \brief Setting of an user interpolation from location.
     *
     * This function takes into account an user interpolation function written with
     *  void (* \ref CWP_Interp_from_target_proc_part_num_t) interface.
     * 
     * \param [in] fct        Function
     *
     */

    void 
    interpFromLocSet
    (
     CWP_Interp_from_target_proc_part_num_t fct
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
     *  void (* \ref CWP_Interp_from_intersec_t) interface.
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
     *  void (* \ref CWP_Interp_from_closest_pts_t) interface.
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
     * \return CWP_Comm_t Communication Type
     *
     */

    inline CWP_Comm_t 
    commTypeGet      
    (
    );
    
    /**
     *
     * \brief Return the Visu object pointer handling visualization
     *
     * \return Visu object pointer 
     *
     */
     
    inline Visu* visuGet ();    
    
    inline Mesh* meshGet();
    
    inline std::map < string, Field * >* fieldsGet();
    inline std::map <CWP_Field_value_t,Geometry*>* geometryGet();
    inline CodeProperties* localCodePropertiesGet();
    inline CodeProperties* coupledCodePropertiesGet();

    inline Communication* communicationGet();

    inline Geometry*    geometryGet(CWP_Field_value_t field_value_t) ;
    inline CouplingDB*  couplingDBGet();
    inline string       IdGet();
    
    void meshFinalize();
        
  private:

    Coupling();

  private:
    const string                            _cplId;                 /*!< Coupling identifier */
          CWP_Comm_t                        _commType;              /*!< Communication type */ 
          Communication                    &_communication;         /*!< Communication */ 
    const CodeProperties                   &_localCodeProperties;   /*!< Local code properties */
    const CodeProperties                   &_coupledCodeProperties; /*!< Coupled code properties */
    std::map <CWP_Field_value_t,Geometry*> &_geometry;              /*!< Geometric algorithm */
          Mesh                             &_mesh;                  /*!< Geometric mesh */
    const CWP_Freq_t                        _recvFreqType  ;        /*!< Receiving frequency type */
          Visu                             &_visu;                  /*!< Visualization */
          double                            _recvFreq;              /*!< Receiving frequency */
          double                            _recvNextTime;          /*!< Next receiving time */
    std::map < string, Field * >           &_fields;          /*!< Fields Data Base */
          CouplingDB                       &_cplDB;                 /*!< Coupling Data base */
          int*                              _iteration;
          CWP_Displacement_t                _displacement;
  }; 



  string Coupling::IdGet(){
     return _cplId;
  }
  
  CodeProperties* Coupling::localCodePropertiesGet() {
    return const_cast<CodeProperties*>(&_localCodeProperties);
  }

  Communication* Coupling::communicationGet() {
    return const_cast<Communication*>(&_communication);
  }

  Geometry* Coupling::geometryGet(CWP_Field_value_t field_value_t) {
  
    std::map <CWP_Field_value_t,Geometry*> ::iterator p;
    p = _geometry.find(field_value_t);
    if (p == _geometry.end()) 
      PDM_error(__FILE__, __LINE__, 0, "Geometry not found.\n");
    return p->second;
  }


  
  CouplingDB* Coupling::couplingDBGet() {
    return &_cplDB;
  }


  CodeProperties* Coupling::coupledCodePropertiesGet() {
    return const_cast<CodeProperties*>(&_coupledCodeProperties);
  }

  Visu* Coupling::visuGet() {
     return &_visu;
  }    
    
  Mesh* Coupling::meshGet() {
     return &_mesh;
  }    

  std::map <CWP_Field_value_t,Geometry*>* Coupling::geometryGet() {
     return &_geometry;
  }    

  std::map < string, Field * >* Coupling::fieldsGet() {
     return &_fields;
  }  



   



    /**
     *
     * \brief Return communication type
     *
     * \return CWP_Comm_t Communication Type
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
