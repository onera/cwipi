#ifndef __MESH_H__
#define __MESH_H__
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

#include <vector>
#include <map>

#include <mpi.h>

#include <pdm_mesh_nodal.h>
#include <bftc_error.h>

#include "cwp.h"




namespace cwipi {



/**
 * \struct _block
 * \brief  Block element adresses              
 *
 * Structure containing block element adresses                 
 */

typedef struct {
    int              _id_part;               /*!< Partition identifier */
    int              _id_block;              /*!< Block identifier */
    CWP_Block_t      _blockType;             /*!< Block Type */
    int              _nElts;                 /*!< Number of elements */
    int*             _connec_idx;            /*!< Connectivity index */
    int*             _connec;                /*!< Connectivity */
    int              _nFaces;                /*!< Number of faces */             
    int*             _connec_faces_idx;      /*!< Faces connectivity indexes (or NULL) */
    int*             _connec_faces;          /*!< Faces connectivity or NULL */
    CWP_g_num_t*     _global_num_block;      /*!< Block global Numbering */
    int*             _parent_num;            /*!< Block parent numbering */
    bool             _isBlockFinalized;      /*!< True if the block is Finalized (set) false otherwise. */
} _block;



  /** 
   * \class Mesh mesh.hxx "mesh.hxx"
   * \brief Geometry mesh
   *
   *  This class computes defines th geometry mesh
   * 
   */

  class Mesh {
    
  public:

    /**
     * \brief Mesh constructor
     * 
     * Construct the CWIPI mesh by using paradigm nodal methods.
     * 
     * \param [in] localComm Coupling Communicator.
     * \param [in] npart     Number of mesh partitions.
     *
     */
 
    Mesh(const MPI_Comm &localComm,
         int npart);


    /**
     * \brief Mesh destructor
     *
     * Destroy the Mesh object.
     *
     */

     virtual ~Mesh();
     
     /**
     * \brief Set the mesh coordinates
     *
     * \param [in] i_part       Index of the mesh partition
     * \param [in] n_pts        Number of partition vertices
     * \param [in] coords       Array containing the verticies coordinates
     * \param [in] global_num   Global numbering of the vertices
     *
     */
   
     void nodal_coord_set(const int   i_part,
                          const int   n_pts,
                          double      coords[],
                          CWP_g_num_t global_num[]); 
                          
     /**
     * \brief Finalize the mesh setting after blocks and coordinates additions
     *
     */
       
     void endSet();
   
    /**
     * \brief Mesh deletion and free memory
     *
     */
   
     void meshDel();
     
     /**
     * \brief Addition of a block (set of cells) to the mesh partition
     *
     * This function add a block to the geometric mesh.
     *
     * \param [in] i_part       Index of the mesh partition
     * \param [in] block_type  Type of the block addition     
     * \param [in] n_elts      Number of block elements
     * \param [in] connec_idx  Vertices to elements connectivity index
     * \param [in] connec      Vertices to elements connectivity 
     * \param [in] n_faces     Number of faces elements (or NULL)
     * \param [in] face_vtx_idx Vertices to faces connectivity index (or NULL)
     * \param [in] face_vtx      Vertices to faces connectivity (or NULL)
     * \param [in] global_num  Global numbering of the vertices in the block (or NULL)
     * \param [in] parent_num  Parent numbering in the block (or NULL)
     *
     */

     int blockAdd(const int                  i_part,
                       const CWP_Block_t      block_type,
                       const int              n_elts,
                       int                    connec_idx[],
                       int                    connec[],
                       const int              n_faces,
                       int                    face_vtx_idx[],
                       int                    face_vtx[],    
                       CWP_g_num_t            global_num[],
                       int                    parent_num[]);

     /**
     * \brief Get the Mesh Block Type
     *
     * This function gets the mesh block type.
     *
     * \param [in] id_block    Block identifier
     * \return CWP_Block_t Mesh  Block Type   
     *
     */

     CWP_Block_t Mesh_nodal_block_type_get(const int   id_block
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
     
     void fromCellFaceSet(const int   i_part,
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
                       
     void fromFacesEdgeSet(const int   i_part,
                           const int   n_faces,
                           int         face_edge_idx[],
                           int         face_edge[],
                           const int   n_edges,
                           int         edge_vtx_idx[],
                           int         edge_vtx[],
                           CWP_g_num_t parent_num[]);
                               

    /**
    * \brief Update the block database of the id_part partition.
    *
    * This function updates the block database of the id_part partition.
    *
    * \param [in] id_part   Partition identifier. 
    *
    */

    void updateBlockDB(int id_part);

    /**
    * \brief Get the vertices number of the partition i_part
    *
    * This function gets the vertices number of the partition i_part.
    *
    * \param [in] i_part   Partition identifier
    * \return Vertex number of i_part partition.   
    *
    */

    inline int& getPartNVertex(int i_part) const;

    /**
    * \brief Get the vertices coordinates of the i_part partition
    *
    * This function gets the vertices number of the i_part partition.
    *
    * \param [in] i_part   Partition identifier
    * \return Vertices coordinates of the i_part partition.   
    *
    */

    inline const std::vector<double*> getVertexCoords(int i_part) const;


    /**
    * \brief Get the paradigm nodal object associated with the mesh.
    *
    * This function gets the paradigm nodal object associated with the mesh.
    *
    * \return Paradigm nodal object associated with the mesh.  
    *
    */

    inline PDM_Mesh_nodal_t& getPdmNodal() const;

    /**
    * \brief Get the number of elements of the id_block block.
    *
    * This function gets the number of elements of the id_block block.
    *
    * \param [in] id_block   Block identifier
    * \return                Number of elements of the id_block block
    *
    */

    inline int getBlockNElts(int id_block);

    /**
    * \brief Get the number of elements of the id_part partition
    *
    * This function gets the number of elements of the id_part partition.
    *
    * \param [in] id_part    Partition identifier
    * \return                Number of elements of the id_part partition
    *
    */

    inline const int getPartNElts(int id_part) const;

    /**
    * \brief Get the number of polyhedra in a partition
    *
    * This function gets the number of polyhedra in a partition.
    *
    * \param [in] i_part   Partition identifier
    * \return              Number of polydra of the i_part partition
    *
    */

    inline const int getPartNPolyhedra(int i_part) const;

    /**
    * \brief Get a block element connectivity index
    *
    * This function gets the element connectivity index of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Element connectivity index of the id_block block
    *
    */

    inline int* getEltConnectivityIndex(int id_block);

    /**
    * \brief Get a block element connectivity
    *
    * This function gets the element connectivity of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Element connectivity of the id_block block
    *
    */
    
    inline int* getEltConnectivity     (int id_block);
 
     /**
    * \brief Get a block face connectivity index
    *
    * This function gets the face connectivity index of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Face connectivity index of the id_block block
    *
    */
    
    inline int* getFaceConnectivityIndex(int id_block);

    /**
    * \brief Get a block face connectivity
    *
    * This function gets the face connectivity of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Face connectivity of the id_block block
    *
    */

    inline int* getFaceConnectivity     (int id_block);
    

   // inline const std::vector<double>& getVolume();

   // inline const std::vector<double>& getCellCenterCoords();


  private:
    
    const MPI_Comm                         &_localComm;       /*!< Communicator */
    int                                     _nDim;            /*!< Entities dimensions */  
    int                                     _nBlocks;         /*!< Number of blocks of the mesh */
    int*                                    _blocks_id;       /*!< List of block identifiers */
    std::vector<int*>                       _blocks_id_part;  /*!< Communicator */
    int                                     _order;           /*!< Mesh order */
    std::vector<int>                        _nVertex;         /*!< Number of vertices for each partition  */
    std::vector<int>                        _nElts;           /*!< Number of elements for each partition  */
    std::vector<double*>                    _coords;          /*!< Vertices coordinate for each partition  */
    std::map< int, CWP_g_num_t*>            _global_num;      /*!< Global numbering for each partition  */
    std::map< int, CWP_g_num_t*>            _parent_num;      /*!< Parent numbering for each partition  */
    int                                     _npart;           /*!< Number of partition  */
    int                                     _pdmNodal_handle_index;  /*!< Mesh (nodal) index for paradigm handler */
    int                                     _pdmGNum_handle_index;   /*!< Global number index for paradigm handler   */
    bool                                    _isNodalFinalized;       /*!< True if all the block of the nodal are finalized   */
    PDM_Mesh_nodal_t                       *_pdmNodal;           /*!< Pointer to the paradigm mesh (nodal) object   */
    std::map<int,_block>                    _blocks;        /*!< Blocks database  */

    
    
    
  //   Mesh &operator=(const Mesh &other);  /*!< Assigment operator not available */
  //   Mesh (const Mesh& other);            /*!< Copy constructor not available */

  };

  const std::vector<int>& Mesh::getPartNVertex(int i_part) const
  {
    return _nVertex[i_part];
  }

  const std::vector<double*> Mesh::getVertexCoords(int i_part)  const
  {
    return _coords[i_part];
  }

  PDM_Mesh_nodal_t& Mesh::getPdmNodal() const
  {
    return *_pdmNodal;
  }


  int Mesh::getBlockNElts(int id_block)
  {
    return _blocks[id_block]._nElts;
  }

  const int Mesh::getPartNElts(int id_part) const
  {
    return _nElts[id_part];
  }

  const int Mesh::getPartNPolyhedra(int i_part) const
  {
   // return _nPolyhedra[i_part];
  }

  inline int* Mesh::getEltConnectivityIndex(int id_block)
  {
    return _blocks[id_block]._connec_idx;
  }

  inline int* Mesh::getEltConnectivity(int id_block)
  { 
    return _blocks[id_block]._connec;
  }
  
  
  inline int* Mesh::getFaceConnectivity(int id_block)
  { 
    return _blocks[id_block]._connec_faces;
  }

  inline int* Mesh::getFaceConnectivityIndex(int id_block)
  { 
    return _blocks[id_block]._connec_faces_idx;
  }





}

#endif //__Mesh_H__
