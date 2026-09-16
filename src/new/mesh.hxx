#ifndef __Mesh_H__
#define __Mesh_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2026  ONERA

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

#include "block.hxx"
#include "blockFP.hxx"
#include "blockCP.hxx"
#include "blockStd.hxx"
#include "blockHO.hxx"

#include "cwp.h"

#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_printf.h"

#include <mpi.h>

#include <vector>

/**
 * \cond
 */

namespace cwipi {

  class Coupling;

  /**
   * \class Mesh mesh.hxx "mesh.hxx"
   * \brief Interface mesh.
   *
   * This class defines the interface mesh objects
   *
   */
  class Mesh {

  public:

    /**
     * \brief Mesh constructor.
     *
     * Construct the CWIPI mesh by using paradigm nodal methods
     *
     * \param [in]  local_comm   Coupling communicator
     * \param [in]  n_part       Number of mesh partitions
     * \param [in]  displacement Type of mesh displacement
     * \param [in]  cpl          Pointer to the Coupling object
     *
     */
    Mesh
    (
      const MPI_Comm     &local_comm,
      int                 n_part,
      CWP_Dynamic_mesh_t  displacement,
      Coupling           *cpl
    );


    /**
     * \brief Mesh destructor.
     *
     * Destroy the Mesh object
     *
     */
    virtual ~Mesh();


    /**
     * \brief Set the vertex coordinates.
     *
     * \param [in]  i_part     Index of the mesh partition
     * \param [in]  n_vtx      Number of partition vertices
     * \param [in]  coords     Array containing the vertices coordinates
     * \param [in]  global_num Global numbering of the vertices or NULL
     *
     */
    void
    coordSet
    (
      const int         i_part,
      const int         n_vtx,
            double      coords[],
            CWP_g_num_t global_num[]
    );


    /**
     * \brief Add a polyhedron block to the geometric mesh from
     *        a face-to-cell connectivity and a vertices-to-faces connectivity.
     *
     * \param [in]  i_part        Current partition
     * \param [in]  n_cells       Number of elements
     * \param [in]  cell_face_idx Polyhedron to face index
     *                            (cell_face_idx[0] = 0 and
     *                             size = n_cells + 1)
     * \param [in]  cell_face     Polyhedron to face connectivity
     *                            (size = cell_face_idx[n_cells])
     * \param [in]  n_faces       Number of faces
     * \param [in]  face_vtx_idx  Face to vertex index
     *                            (face_vtx_idx[0] = 0 and
     *                             size = max(face_vtx) + 1)
     * \param [in]  face_vtx      Face to vertex connectivity
     *                            (size = face_vtx_idx[size - 1])
     * \param [in]  global_num    Cell global numbering (or NULL)
     *
     */
    void
    fromCellFaceSet
    (
      const int         i_part,
      const int         n_cells,
            int         cell_face_idx[],
            int         cell_face[],
            int         n_faces,
            int         face_vtx_idx[],
            int         face_vtx[],
            CWP_g_num_t global_num[]
    );


    /**
     * \brief Add a polygon 2D block to the geometric mesh from
     *        a face-to-edge connectivity and an edge-to-vertex connectivity.
     *
     * \param [in]  i_part        Current partition
     * \param [in]  n_faces       Number of faces
     * \param [in]  face_edge_idx Face to edge index
     *                            (face_edge_idx[0] = 0 and
     *                             size = n_faces + 1)
     * \param [in]  face_edge     Face to edge connectivity
     *                            (size = face_edge_idx[n_faces])
     * \param [in]  n_edges       Number of edges
     * \param [in]  edge_vtx      Edge to vertex connectivity
     *                            (size = edge_vtx[n_edges])
     * \param [in]  global_num    Face global numbering (or NULL)
     *
     */
    void
    fromFaceEdgeSet
    (
      const int         i_part,
      const int         n_faces,
            int         face_edge_idx[],
            int         face_edge[],
      const int         n_edges,
            int         edge_vtx[],
            CWP_g_num_t global_num[]
    );


    /**
     * \brief Add a polygon 2D block to the geometric mesh from
     *        a face-to-vertex connectivity.
     *
     * \param [in]  i_part       Current partition
     * \param [in]  n_faces      Number of faces
     * \param [in]  face_vtx_idx Face to vertex index
     * \param [in]  face_vtx     Face to vertex connectivity
     * \param [in]  global_num   Face global numbering (or NULL)
     *
     */
    void
    fromFaceVtxSet
    (
      const int         i_part,
      const int         n_faces,
            int         face_vtx_idx[],
            int         face_vtx[],
            CWP_g_num_t global_num[]
    );


    /**
     * \brief Add a polyhedron 3D block to the geometric mesh from
     *        a cell-to-vertex connectivity.
     *
     * \param [in]  i_part       Current partition
     * \param [in]  n_cells      Number of cells
     * \param [in]  cell_vtx_idx Cell to vertex index
     * \param [in]  cell_vtx     Cell to vertex connectivity
     * \param [in]  global_num   Cell global numbering (or NULL)
     *
     */
    void
    fromCellVtxSet
    (
      const int         i_part,
      const int         n_cells,
            int         cell_vtx_idx[],
            int         cell_vtx[],
            CWP_g_num_t global_num[]
    );


    /**
     * \brief Addition of a block (set of cells) to the mesh partition.
     *
     * \param [in]  block_type Type of block
     *
     * \return Block identifier
     */
    int
    blockAdd
    (
      const CWP_Block_t block_type
    );


    /**
     * \brief Set a standard block to the interface mesh.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  block_id   Block identifier
     * \param [in]  n_elts     Number of elements
     * \param [in]  connec     Elements to vertices connectivity
     * \param [in]  global_num Element global numbering (or NULL)
     *
     */
    void
    stdBlockSet
    (
      const int         i_part,
      const int         block_id,
      const int         n_elts,
            int         connec[],
            CWP_g_num_t global_num[]
    );


    /**
     * \brief Set a high-order block to the interface mesh.
     *
     * \param [in]  i_part      Partition identifier
     * \param [in]  block_id    Block identifier
     * \param [in]  n_elts      Number of elements
     * \param [in]  connec      Elements to vertices connectivity
     * \param [in]  global_num  Element global numbering (or NULL)
     * \param [in]  order       Element order
     * \param [in]  ho_ordering HO ordering name
     *
     */
    void
    hoBlockSet
    (
      const int          i_part,
      const int          block_id,
      const int          n_elts,
            int          connec[],
            CWP_g_num_t  global_num[],
      const int          order,
      const char        *ho_ordering
    );


    /**
     * \brief Set a face polygon block to the interface mesh.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  block_id   Block identifier
     * \param [in]  n_elts     Number of elements
     * \param [in]  connec_idx Elements to vertices connectivity index
     * \param [in]  connec     Elements to vertices connectivity
     * \param [in]  global_num Element global numbering (or NULL)
     *
     */
    void
    poly2DBlockSet
    (
      const int         i_part,
      const int         block_id,
      const int         n_elts,
            int         connec_idx[],
            int         connec[],
            CWP_g_num_t global_num[]
    );


    /**
     * \brief Set a face polyhedron block to the interface mesh.
     *
     * \param [in]  i_part           Partition identifier
     * \param [in]  block_id         Block identifier
     * \param [in]  n_elts           Number of elements
     * \param [in]  n_faces          Number of faces
     * \param [in]  connec_faces_idx Faces to vertices connectivity index
     * \param [in]  connec_faces     Faces to vertices connectivity
     * \param [in]  connec_cells_idx Cells to faces connectivity index
     * \param [in]  connec_cells     Cells to faces connectivity
     * \param [in]  global_num       Element global numbering (or NULL)
     *
     */
    void
    poly3DBlockSet
    (
      const int         i_part,
      const int         block_id,
      const int         n_elts,
      const int         n_faces,
            int         connec_faces_idx[],
            int         connec_faces[],
            int         connec_cells_idx[],
            int         connec_cells[],
            CWP_g_num_t global_num[]
    );


    /**
     * \brief Get a standard block from the interface mesh.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  block_id   Block identifier
     * \param [out] n_elts     Number of elements
     * \param [out] connec     Elements to vertices connectivity
     * \param [out] global_num Element global numbering (or NULL)
     *
     */
    void
    stdBlockGet
    (
      const int           i_part,
      const int           block_id,
            int          *n_elts,
            int         **connec,
            CWP_g_num_t **global_num
    );


    /**
     * \brief Get a high-order block from the interface mesh.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  block_id   Block identifier
     * \param [out] n_elts     Number of elements
     * \param [out] order      Element order
     * \param [out] connec     Elements to vertices connectivity
     * \param [out] global_num Element global numbering (or NULL)
     *
     */
    void
    hoBlockGet
    (
      const int           i_part,
      const int           block_id,
            int          *n_elts,
            int          *order,
            int         **connec,
            CWP_g_num_t **global_num
    );


    /**
     * \brief Get a face polygon block from the interface mesh.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  block_id   Block identifier
     * \param [out] n_elts     Number of elements
     * \param [out] connec_idx Elements to vertices connectivity index
     * \param [out] connec     Elements to vertices connectivity
     * \param [out] global_num Element global numbering (or NULL)
     *
     */
    void
    poly2DBlockGet
    (
      const int           i_part,
      const int           block_id,
            int          *n_elts,
            int         **connec_idx,
            int         **connec,
            CWP_g_num_t **global_num
    );


    /**
     * \brief Get a face polyhedron block.
     *
     * \param [in]  i_part           Partition identifier
     * \param [in]  block_id         Block identifier
     * \param [out] n_elts           Number of elements
     * \param [out] n_faces          Number of faces
     * \param [out] connec_faces_idx Faces to vertices connectivity index
     * \param [out] connec_faces     Faces to vertices connectivity
     * \param [out] connec_cells_idx Cells to faces connectivity index
     * \param [out] connec_cells     Cells to faces connectivity
     * \param [out] global_num       Element global numbering (or NULL)
     *
     */
    void
    poly3DBlockGet
    (
      const int           i_part,
      const int           block_id,
            int          *n_elts,
            int          *n_faces,
            int         **connec_faces_idx,
            int         **connec_faces,
            int         **connec_cells_idx,
            int         **connec_cells,
            CWP_g_num_t **global_num
    );


    /**
     * \brief Get the block type.
     *
     * \param [in]  id_block Block identifier
     * \return Block type
     *
     */
    inline CWP_Block_t
    blockTypeGet
    (
      int id_block
    ) const;


    /**
     * \brief Get the number of partitions.
     *
     * \return Number of partitions
     *
     */
    inline int
    nPartGet() const;


    /**
     * \brief Get the number of vertices of the partition.
     *
     * \param [in]  i_part Partition identifier
     * \return Number of vertices of the partition
     *
     */
    inline int
    nVtxGet
    (
      int i_part
    ) const;


    /**
     * \brief Get the number of elements of the partition.
     *
     * \param [in]  i_part Partition identifier
     * \return Number of elements of the partition
     *
     */
    int
    nEltsGet
    (
      int i_part
    ) const;


    /**
     * \brief Get the vertices coordinates of the partition.
     *
     * \param [in]  i_part Partition identifier
     * \return Vertices coordinates of the partition
     *
     */
    inline double*
    coordGet
    (
      int i_part
    ) const;


    /**
     * \brief Get the vertices global numbering of the partition.
     *
     * \param [in]  i_part Partition identifier
     * \return Vertices global numbering of the partition
     *
     */
    inline CWP_g_num_t*
    vertexGlobalNumGet
    (
      int i_part
    ) const;


    /**
     * \brief Get the global numbering of a block.
     *
     * \param [in]  id_block Block identifier
     * \param [in]  i_part   Partition identifier
     * \return Global numbering of the block
     *
     */
    inline CWP_g_num_t*
    blockGlobalNumGet
    (
      int id_block,
      int i_part
    ) const;


    /**
     * \brief Get the global numbering of elements for a partition.
     *
     * \param [in]  i_part Partition identifier
     * \return Global numbering of elements
     *
     */
    CWP_g_num_t*
    elementGlobalNumGet
    (
      int i_part
    );


    /**
     * \brief Get the element centers for a partition.
     *
     * \param [in]  i_part Partition identifier
     * \return Element centers
     *
     */
    double*
    elementCentersGet
    (
      int i_part
    );


    /**
     * \brief Compute the element centers for a partition.
     *
     * \param [in]  i_part Partition identifier
     *
     */
    void
    elementCentersCompute
    (
      int i_part
    );


    /**
     * \brief Return pdmNodal handler identifier.
     *
     * \return Pdm nodal handler identifier
     *
     */
    inline PDM_part_mesh_nodal_t*
    pdmNodalIndexGet() const;


    /**
     * \brief Return MPI communicator pointer.
     *
     * \return MPI communicator pointer
     *
     */
    inline const MPI_Comm*
    mpiCommGet();


    /**
     * \brief Finalize the geometry computation.
     *
     */
    void
    geomFinalize();


    /**
     * \brief Update the block database.
     *
     */
    void
    blockDBUpdate();


    /**
     * \brief Mesh deletion and free memory.
     *
     */
    void
    meshDel();


    /**
     * \brief Get the nodal standard type block size.
     *
     * \param [in]  type Block type
     * \return Size of the nodal standard type block
     *
     */
    int
    blockStdTypeSizeGet
    (
      CWP_Block_t type
    );

  private:

    const MPI_Comm                &_localComm;            /*!< Communicator */
    int                           _npart;                 /*!< Number of partitions */
    int                           _nBlocks;               /*!< Number of blocks of the mesh */
    int*                          _blocks_id;             /*!< List of block identifiers */
    CWP_Block_t*                  _blocksType;            /*!< Blocks type */
    std::vector<cwipi::Block*>    _blockDB;               /*!< Blocks database */

    std::vector<int>              _nVertex;               /*!< Number of vertices for each partition */
    std::vector<double*>          _coords;                /*!< Vertices coordinate for each partition */
    std::vector<CWP_g_num_t*>     _global_num_vtx;        /*!< Global vertices numbering for each partition */

    std::vector<int>              _nElts;                 /*!< Number of elements for each partition */
    std::vector<double*>          _elt_centers;           /*!< Element centers for each partition */
    std::vector<int*>             _elt_id_block;          /*!< Element to block mapping for each partition */
    std::vector<int*>             _elt_in_block;          /*!< Element index in block for each partition */

    CWP_Dynamic_mesh_t            _displacement;          /*!< Type of mesh displacement */
    Coupling                      *_cpl;                  /*!< Pointer to the Coupling object */

    std::vector<int>              _nCell;                 /*!< Number of cells for each partition */
    std::vector<int*>             _cellFaceIdx;           /*!< Cell to face index for each partition */
    std::vector<int*>             _cellFace;              /*!< Cell to face connectivity for each partition */
    std::vector<int*>             _cellVtxIdx;            /*!< Cell to vertex index for each partition */
    std::vector<int*>             _cellVtx;               /*!< Cell to vertex connectivity for each partition */
    std::vector<CWP_g_num_t*>     _cellLNToGN;            /*!< Cell local to global numbering for each partition */

    std::vector<int>              _nFace;                 /*!< Number of faces for each partition */
    std::vector<int*>             _faceEdgeIdx;           /*!< Face to edge index for each partition */
    std::vector<int*>             _faceEdge;              /*!< Face to edge connectivity for each partition */
    std::vector<int*>             _faceVtxIdx;            /*!< Face to vertex index for each partition */
    std::vector<int*>             _faceVtx;               /*!< Face to vertex connectivity for each partition */
    std::vector<CWP_g_num_t*>     _faceLNToGN;            /*!< Face local to global numbering for each partition */

    std::vector<int>              _nEdge;                 /*!< Number of edges for each partition */
    std::vector<int*>             _edgeVtx;               /*!< Edge to vertex connectivity for each partition */

    PDM_MPI_Comm                  _pdm_localComm;         /*!< Paradigm local MPI communicator */

    int                           _faceEdgeMethod;        /*!< Face to edge method identifier */
    int                           _faceVtxMethod;         /*!< Face to vertex method identifier */
    int                           _cellFaceMethod;        /*!< Cell to face method identifier */
    int                           _cellVtxMethod;         /*!< Cell to vertex method identifier */

    PDM_part_mesh_nodal_t        *_pdmNodal_handle_index; /*!< Mesh (nodal) index for paradigm handler */
    PDM_geometry_kind_t           _geom_kind;             /*!< Geometry kind */

    bool                          _isVtxGnumComputed;     /*!< Flag: vertex global numbering computed */
    bool                          _isEltGnumComputed;     /*!< Flag: element global numbering computed */
    bool                          _isFaceGnumComputed;    /*!< Flag: face global numbering computed */

  };

  inline int
  Mesh::nPartGet() const
  {
    return _npart;
  }


  inline PDM_part_mesh_nodal_t*
  Mesh::pdmNodalIndexGet() const
  {
    return _pdmNodal_handle_index;
  }


  inline const MPI_Comm*
  Mesh::mpiCommGet()
  {
    return &_localComm;
  }


  inline int
  Mesh::nVtxGet
  (
    int i_part
  ) const
  {
    return _nVertex[i_part];
  }


  inline double*
  Mesh::coordGet
  (
    int i_part
  ) const
  {
    return _coords[i_part];
  }


  inline CWP_g_num_t*
  Mesh::vertexGlobalNumGet
  (
    int i_part
  ) const
  {
    return _global_num_vtx[i_part];
  }


  inline CWP_Block_t
  Mesh::blockTypeGet
  (
    int id_block
  ) const
  {
    return _blockDB[id_block]->blockTypeGet();
  }


  inline CWP_g_num_t*
  Mesh::blockGlobalNumGet
  (
    int id_block,
    int i_part
  ) const
  {
    return _blockDB[id_block]->globalNumGet(i_part);
  }

}

/**
 * \endcond
 */

#endif // __Mesh_H__