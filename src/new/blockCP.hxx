#ifndef __BLOCKCP_H__
#define __BLOCKCP_H__
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

#include "cwp.h"

#include <vector>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class BlockCP blockCP.hxx "blockCP.hxx"
   * \brief Polyhedron 3D Block Mesh definition.
   *
   *  This class defines the polyhedron 3D Block objects.
   *  It inherits from the Block class through a Factory.
   *
   */
  class BlockCP :
    public Block {

  public:

    /**
     * \brief Default constructor.
     *
     */
    BlockCP();


    /**
     * \brief Destructor.
     *
     */
    virtual ~BlockCP();


    /**
     * \brief Add a block to the mesh.
     *
     * \param [in]  block_type Type of the block
     * \param [in]  mesh       Mesh object owning the block
     *
     */
    virtual void
    blockAdd
    (
      const CWP_Block_t  block_type,
            Mesh        *mesh
    ) override;


    /**
     * \brief Set a CWIPI block in a partition.
     *
     * \param [in]  i_part           Partition identifier
     * \param [in]  n_elts           Number of elements of the block
     * \param [in]  n_faces          Number of faces of the block
     * \param [in]  connec_faces_idx Face to vertex connectivity index
     * \param [in]  connec_faces     Face to vertex connectivity
     * \param [in]  connec_cells_idx Cell to face connectivity index
     * \param [in]  connec_cells     Cell to face connectivity
     * \param [in]  global_num       Global numbering of the block
     *
     */
    void
    blockSet
    (
      const int          i_part,
      const int          n_elts,
      const int          n_faces,
            int         *connec_faces_idx,
            int         *connec_faces,
            int         *connec_cells_idx,
            int         *connec_cells,
            CWP_g_num_t *global_num
    );


    /**
     * \brief Get a CWIPI block in a partition.
     *
     * \param [in]  i_part           Partition identifier
     * \param [out] n_elts           Number of elements of the block
     * \param [out] n_faces          Number of faces of the block
     * \param [out] connec_faces_idx Face to vertex connectivity index
     * \param [out] connec_faces     Face to vertex connectivity
     * \param [out] connec_cells_idx Cell to face connectivity index
     * \param [out] connec_cells     Cell to face connectivity
     * \param [out] global_num       Global numbering of the block
     *
     */
    void
    blockGet
    (
      const int           i_part,
            int          *n_elts,
            int          *n_faces,
            int         **connec_faces_idx,
            int         **connec_faces,
            int         **connec_cells_idx,
            int         **connec_cells,
            CWP_g_num_t **global_num
    );


    /**
     * \brief Return the cell to face connectivity for a partition.
     *
     * \param [in]  i_part Partition identifier
     *
     * \return Connectivity array for the partition
     *
     */
    inline int*
    cellConnecGet
    (
      const int i_part
    ) const;


    /**
     * \brief Return the cell to face connectivity index for a partition.
     *
     * \param [in]  i_part Partition identifier
     *
     * \return Connectivity index array for the partition
     *
     */
    inline int*
    cellConnecIdxGet
    (
      const int i_part
    ) const;


    /**
     * \brief Return the number of faces for a partition (Cells_POLY_3D).
     *
     * \param [in]  i_part Partition identifier
     *
     * \return Number of faces for the partition
     *
     */
    inline int
    nFacesGet
    (
      const int i_part
    ) const;


    /**
     * \brief Return the number of faces for all partitions.
     *
     * \return Vector of face counts for each partition
     *
     */
    inline const std::vector<int>&
    nFacesGet() const;


    /**
     * \brief Return the face to vertex connectivity for a partition.
     *
     * \param [in]  i_part Partition identifier
     *
     * \return Faces connectivity array for the partition
     *
     */
    inline int*
    faceConnecGet
    (
      const int i_part
    ) const;


    /**
     * \brief Return the face to vertex connectivity for all partitions.
     *
     * \return Vector of connectivity arrays for each partition
     *
     */
    inline const std::vector<int*>&
    faceConnecGet() const;


    /**
     * \brief Return the face to vertex connectivity index for a partition.
     *
     * \param [in]  i_part Partition identifier
     *
     * \return Faces connectivity index array for the partition
     *
     */
    inline int*
    faceConnecIdxGet
    (
      const int i_part
    ) const;


    /**
     * \brief Return the face to vertex connectivity index for all partitions.
     *
     * \return Vector of connectivity index arrays for each partition
     *
     */
    inline const std::vector<int*>&
    faceConnecIdxGet() const;


    /**
     * \brief Set the faces global numbering for a partition.
     *
     * \param [in]  i_part           Partition identifier
     * \param [in]  global_num_faces Faces global numbering array
     *
     */
    inline void
    faceGlobalNumSet
    (
      const int          i_part,
            CWP_g_num_t *global_num_faces
    );


    /**
     * \brief Set the faces global numbering for all partitions.
     *
     * \param [in]  global_num_faces Vector of global numbering arrays for each partition
     *
     */
    inline void
    faceGlobalNumSet
    (
      const std::vector<CWP_g_num_t*>& global_num_faces
    );


    /**
     * \brief Return the faces global numbering for a partition.
     *
     * \param [in]  i_part Partition identifier
     *
     * \return Faces global numbering array for the partition
     *
     */
    inline CWP_g_num_t*
    faceGlobalNumGet
    (
      const int i_part
    ) const;


    /**
     * \brief Free the faces global numbering.
     *
     * \param [in]  i_part Partition identifier
     *
     */
    void
    faceGlobalNumFree
    (
      const int i_part
    );


    /**
     * \brief Finalize geometry processing.
     *
     */
    void
    geomFinalize() override;

  private:

    std::vector<int>          _n_faces;          /*!< Number of faces for each partition */
    std::vector<int*>         _connec_faces_idx; /*!< Faces connectivity index for each partition */
    std::vector<int*>         _connec_faces;     /*!< Faces connectivity for each partition */
    std::vector<int*>         _connec_cells_idx; /*!< Cells connectivity index for each partition */
    std::vector<int*>         _connec_cells;     /*!< Cells connectivity for each partition */
    std::vector<CWP_g_num_t*> _global_num_faces; /*!< Faces global numbering */

  };

  inline int*
  BlockCP::cellConnecGet
  (
    const int i_part
  ) const
  {
    return _connec_cells[i_part];
  }


  inline int*
  BlockCP::cellConnecIdxGet
  (
    const int i_part
  ) const
  {
    return _connec_cells_idx[i_part];
  }


  inline int
  BlockCP::nFacesGet
  (
    const int i_part
  ) const
  {
    return _n_faces[i_part];
  }


  inline const std::vector<int>&
  BlockCP::nFacesGet() const
  {
    return _n_faces;
  }


  inline int*
  BlockCP::faceConnecGet
  (
    const int i_part
  ) const
  {
    return _connec_faces[i_part];
  }


  inline const std::vector<int*>&
  BlockCP::faceConnecGet() const
  {
    return _connec_faces;
  }


  inline int*
  BlockCP::faceConnecIdxGet
  (
    const int i_part
  ) const
  {
    return _connec_faces_idx[i_part];
  }


  inline const std::vector<int*>&
  BlockCP::faceConnecIdxGet() const
  {
    return _connec_faces_idx;
  }


  inline void
  BlockCP::faceGlobalNumSet
  (
    const int          i_part,
          CWP_g_num_t *global_num_faces
  )
  {
    _global_num_faces[i_part] = global_num_faces;
  }


  inline void
  BlockCP::faceGlobalNumSet
  (
    const std::vector<CWP_g_num_t*>& global_num_faces
  )
  {
    _global_num_faces = global_num_faces;
  }


  inline CWP_g_num_t*
  BlockCP::faceGlobalNumGet
  (
    const int i_part
  ) const
  {
    return _global_num_faces[i_part];
  }

}

/**
 * \endcond
 */

#endif // __BLOCKCP_H__