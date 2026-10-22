#ifndef __MESH_H__
#define __MESH_H__
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

#include "cwp.h"

#include "block.hxx"
#include "blockCP.hxx"
#include "blockFP.hxx"
#include "blockHO.hxx"
#include "blockStd.hxx"

#include "pdm_part_mesh_nodal.h"

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
     * \param [in]  global_num Global numbering of the vertices (or NULL)
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
     * \param [in]  connec     Element to vertex connectivity
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
     * \brief Get a standard block from the interface mesh.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  block_id   Block identifier
     * \param [out] n_elts     Number of elements
     * \param [out] connec     Element to vertex connectivity
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
     * \brief Set a high-order block to the interface mesh.
     *
     * \param [in]  i_part      Partition identifier
     * \param [in]  block_id    Block identifier
     * \param [in]  n_elts      Number of elements
     * \param [in]  connec      Element to vertex connectivity
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
     * \brief Get a high-order block from the interface mesh.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  block_id   Block identifier
     * \param [out] n_elts     Number of elements
     * \param [out] order      Element order
     * \param [out] connec     Element to vertex connectivity
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
     * \brief Set a face polygon block to the interface mesh.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  block_id   Block identifier
     * \param [in]  n_elts     Number of elements
     * \param [in]  connec_idx Element to vertex connectivity index
     * \param [in]  connec     Element to vertex connectivity
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
     * \brief Get a face polygon block from the interface mesh.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  block_id   Block identifier
     * \param [out] n_elts     Number of elements
     * \param [out] connec_idx Element to vertex connectivity index
     * \param [out] connec     Element to vertex connectivity
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
     * \brief Set a face polyhedron block to the interface mesh.
     *
     * \param [in]  i_part           Partition identifier
     * \param [in]  block_id         Block identifier
     * \param [in]  n_elts           Number of elements
     * \param [in]  n_faces          Number of faces
     * \param [in]  connec_faces_idx Face to vertex connectivity index
     * \param [in]  connec_faces     Face to vertex connectivity
     * \param [in]  connec_cells_idx Cell to face connectivity index
     * \param [in]  connec_cells     Cell to face connectivity
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
     * \brief Get a face polyhedron block from the interface mesh.
     *
     * \param [in]  i_part           Partition identifier
     * \param [in]  block_id         Block identifier
     * \param [out] n_elts           Number of elements
     * \param [out] n_faces          Number of faces
     * \param [out] connec_faces_idx Face to vertex connectivity index
     * \param [out] connec_faces     Face to vertex connectivity
     * \param [out] connec_cells_idx Cell to face connectivity index
     * \param [out] connec_cells     Cell to face connectivity
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
    ) const;


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
     * \brief Return part mesh nodal identifier.
     *
     * \return Part mesh nodal identifier
     *
     */
    inline PDM_part_mesh_nodal_t*
    pdmPartMeshNodalGet() const;


    /**
     * \brief Finalize the geometry computation.
     *
     */
    void
    geomFinalize();


    /**
     * \brief Mesh deletion and free memory.
     *
     */
    void
    meshDel();

  private:

    /**
     * \brief Update the block database.
     *
     */
    void
    _blockDBUpdate();


    /**
     * \brief Return the mesh dimension from the interface entity dimension.
     *
     * \return Mesh dimension (0, 1, 2 or 3)
     */
    int
    _meshDimensionGet() const;


    int                           _n_part;                /*!< Number of partitions */
    int                           _n_blocks;              /*!< Number of blocks of the mesh */
    std::vector<cwipi::Block*>    _block_DB;              /*!< Blocks database */

    std::vector<int>              _n_vtx;                 /*!< Number of vertices for each partition */
    std::vector<double*>          _coords;                /*!< Vertices coordinate for each partition */
    std::vector<CWP_g_num_t*>     _global_num_vtx;        /*!< Vertices global numbering for each partition */

    std::vector<int>              _n_elts;                /*!< Number of elements for each partition */
    std::vector<double*>          _elt_centers;           /*!< Element centers for each partition */

    CWP_Dynamic_mesh_t            _displacement;          /*!< Type of mesh displacement */
    Coupling                     *_cpl;                   /*!< Pointer to the Coupling object */

    std::vector<int>              _n_cells;               /*!< Number of cells for each partition */
    std::vector<int*>             _cell_face_idx;         /*!< Cell to face index for each partition */
    std::vector<int*>             _cell_face;             /*!< Cell to face connectivity for each partition */
    std::vector<int*>             _cell_vtx_idx;          /*!< Cell to vertex index for each partition */
    std::vector<int*>             _cell_vtx;              /*!< Cell to vertex connectivity for each partition */
    std::vector<CWP_g_num_t*>     _global_num_cells;      /*!< Cell global numbering for each partition */

    std::vector<int>              _n_faces;               /*!< Number of faces for each partition */
    std::vector<int*>             _face_edge_idx;         /*!< Face to edge index for each partition */
    std::vector<int*>             _face_edge;             /*!< Face to edge connectivity for each partition */
    std::vector<int*>             _face_vtx_idx;          /*!< Face to vertex index for each partition */
    std::vector<int*>             _face_vtx;              /*!< Face to vertex connectivity for each partition */
    std::vector<CWP_g_num_t*>     _global_num_faces;      /*!< Face global numbering for each partition */

    std::vector<int>              _n_edges;               /*!< Number of edges for each partition */
    std::vector<int*>             _edge_vtx;              /*!< Edge to vertex connectivity for each partition */

    PDM_MPI_Comm                  _local_comm_pdm;        /*!< Paradigm local MPI communicator */

    PDM_part_mesh_nodal_t        *_pmn;                   /*!< Paradigm part mesh nodal identifier */
    PDM_geometry_kind_t           _geom_kind;             /*!< Geometry kind */

    bool                          _is_face_edge_method;   /*!< Flag: face to edge   method */
    bool                          _is_face_vtx_method;    /*!< Flag: face to vertex method */
    bool                          _is_cell_face_method;   /*!< Flag: cell to face   method */
    bool                          _is_cell_vtx_method;    /*!< Flag: cell to vertex method */

    bool                          _is_vtx_gnum_computed;  /*!< Flag: vertex  global numbering computed */
    bool                          _is_elt_gnum_computed;  /*!< Flag: element global numbering computed */
    bool                          _is_face_gnum_computed; /*!< Flag: face    global numbering computed */

  };

  inline CWP_Block_t
  Mesh::blockTypeGet
  (
    int id_block
  ) const
  {
    return _block_DB[id_block]->blockTypeGet();
  }


  inline int
  Mesh::nPartGet() const
  {
    return _n_part;
  }


  inline int
  Mesh::nVtxGet
  (
    int i_part
  ) const
  {
    return _n_vtx[i_part];
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


  inline CWP_g_num_t*
  Mesh::blockGlobalNumGet
  (
    int id_block,
    int i_part
  ) const
  {
    return _block_DB[id_block]->globalNumGet(i_part);
  }


  inline PDM_part_mesh_nodal_t*
  Mesh::pdmPartMeshNodalGet() const
  {
    return _pmn;
  }

}

/**
 * \endcond
 */

#endif // __MESH_H__