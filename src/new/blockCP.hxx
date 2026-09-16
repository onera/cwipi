#ifndef __BlockCP_H__
#define __BlockCP_H__
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
     * \param [in]  blockType Type of the block
     * \param [in]  mesh      Mesh object owning the block
     */
    virtual void
    BlockAdd
    (
      CWP_Block_t  blockType,
      Mesh        *mesh
    );


    /**
     * \brief Set a CWIPI block in a partition.
     *
     * \param [in]  i_part           Partition identifier
     * \param [in]  n_elts           Number of elements of the block in the partition
     * \param [in]  n_faces          Number of faces of the block in the partition
     * \param [in]  connec_faces_idx Vertices to faces connectivity index
     * \param [in]  connec_faces     Vertices to faces connectivity
     * \param [in]  connec_cells_idx Faces to cells connectivity index
     * \param [in]  connec_cells     Faces to cells connectivity
     * \param [in]  global_num       Global element numbering
     */
    void
    blockSet
    (
      int          i_part,
      int          n_elts,
      int          n_faces,
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
     * \param [out] n_elts           Number of elements of the block in the partition
     * \param [out] n_faces          Number of faces of the block in the partition
     * \param [out] connec_faces_idx Vertices to faces connectivity index
     * \param [out] connec_faces     Vertices to faces connectivity
     * \param [out] connec_cells_idx Faces to cells connectivity index
     * \param [out] connec_cells     Faces to cells connectivity
     * \param [out] global_num       Global element numbering
     */
    void
    blockGet
    (
      int           i_part,
      int          *n_elts,
      int          *n_faces,
      int         **connec_faces_idx,
      int         **connec_faces,
      int         **connec_cells_idx,
      int         **connec_cells,
      CWP_g_num_t **global_num
    );


    /**
     * \brief Return the element connectivity (Standard or Face_Poly_2D)
     *        or cells-faces connectivity (Cells_POLY_3D) for each partition.
     *
     * \return Vector of connectivity arrays for each partition
     */
    inline std::vector<int*>
    ConnecGet();


    /**
     * \brief Return the element connectivity index (Face_Poly_2D)
     *        or cells-faces connectivity index (Cells_POLY_3D) for each partition.
     *
     * \return Vector of connectivity index arrays for each partition
     */
    inline std::vector<int*>
    ConnecIDXGet();


    /**
     * \brief Return the element faces connectivity (Cells_POLY_3D)
     *        for each partition.
     *
     * \return Vector of faces connectivity arrays for each partition
     */
    inline std::vector<int*>
    ConnecFacesGet();


    /**
     * \brief Return the element faces connectivity index (Cells_POLY_3D)
     *        for each partition.
     *
     * \return Vector of faces connectivity index arrays for each partition
     */
    inline std::vector<int*>
    ConnecFacesIDXGet();


    /**
     * \brief Return the number of faces for each partition (Cells_POLY_3D).
     *
     * \return Vector of face counts for each partition
     */
    inline std::vector<int>
    NFacesGet();


    /**
     * \brief Return the faces global numbering for each partition.
     *
     * \return Vector of faces global numbering arrays
     */
    inline std::vector<CWP_g_num_t*>
    FacesGNumGet();


    /**
     * \brief Set the faces global numbering.
     *
     * \param [in] face_ln_to_gn Faces global numbering vector
     */
    inline void
    FacesGNumSet
    (
      std::vector<CWP_g_num_t*> face_ln_to_gn
    );


    /**
     * \brief Finalize geometry processing.
     *
     */
    void
    geomFinalize();


    /**
     * \brief Free the faces global numbering.
     *
     * \param [in] i_part Partition identifier
     */
    void
    FacesGNumFree
    (
      int i_part
    );

  private:

    std::vector<int>          _n_faces;          /*!< Number of faces for each partition */
    std::vector<int*>         _connec_faces_idx; /*!< Faces connectivity index for each partition */
    std::vector<int*>         _connec_faces;     /*!< Faces connectivity for each partition */
    std::vector<int*>         _connec_cells_idx; /*!< Cells connectivity index for each partition */
    std::vector<int*>         _connec_cells;     /*!< Cells connectivity for each partition */
    std::vector<CWP_g_num_t*> _face_ln_to_gn;    /*!< Face global numbering */

  };

  inline std::vector<int*>
  BlockCP::ConnecGet()
  {
    return _connec_cells;
  }


  inline std::vector<int*>
  BlockCP::ConnecIDXGet()
  {
    return _connec_cells_idx;
  }


  inline std::vector<int*>
  BlockCP::ConnecFacesGet()
  {
    return _connec_faces;
  }


  inline std::vector<int*>
  BlockCP::ConnecFacesIDXGet()
  {
    return _connec_faces_idx;
  }


  inline std::vector<int>
  BlockCP::NFacesGet()
  {
    return _n_faces;
  }


  inline std::vector<CWP_g_num_t*>
  BlockCP::FacesGNumGet()
  {
    return _face_ln_to_gn;
  }


  inline void
  BlockCP::FacesGNumSet
  (
    std::vector<CWP_g_num_t*> face_ln_to_gn
  )
  {
    _face_ln_to_gn = face_ln_to_gn;
  }

}

/**
 * \endcond
 */

#endif // __BlockCP_H__