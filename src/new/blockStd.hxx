#ifndef __BlockStd_H__
#define __BlockStd_H__
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

#include <cstdlib>
#include <vector>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class BlockStd blockStd.hxx "blockStd.hxx"
   * \brief Standard Block Mesh definition.
   *
   *  This class defines the standard Block objects.
   *  It inherits from the Block class through a Factory.
   *
   */
  class BlockStd :
    public Block {

  public:

    /**
     * \brief Default constructor.
     *
     */
    BlockStd();


    /**
     * \brief Destructor.
     *
     */
    virtual ~BlockStd();


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
     * \param [in]  i_part     Partition identifier
     * \param [in]  n_elts     Number of elements of the block in the partition
     * \param [in]  connec     Elements connectivity
     * \param [in]  global_num Mesh global numbering of the block
     */
    void
    blockSet
    (
      int          i_part,
      int          n_elts,
      int         *connec,
      CWP_g_num_t *global_num
    );


    /**
     * \brief Get a CWIPI block.
     *
     * \param [in]  i_part     Partition identifier
     * \param [out] n_elts     Number of elements of the block in the partition
     * \param [out] connec     Elements connectivity
     * \param [out] global_num Mesh global numbering of the block
     */
    void
    blockGet
    (
      int           i_part,
      int          *n_elts,
      int         **connec,
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
     * \brief Return the element connectivity for a specific partition.
     *
     * \param [in] i_part Partition identifier
     * \return Connectivity array for the specified partition
     */
    inline int*
    ConnecGet
    (
      int i_part
    );


    /**
     * \brief Finalize geometry processing.
     *
     *  This method should be called after all geometric data has been set.
     */
    void
    geomFinalize();

  private:

    std::vector<int*> _connec; /*!< Connectivity for each partition */

  };

  inline std::vector<int*>
  BlockStd::ConnecGet()
  {
    return _connec;
  }


  inline int*
  BlockStd::ConnecGet
  (
    int i_part
  )
  {
    return _connec[i_part];
  }

}

/**
 * \endcond
 */

#endif // __BlockStd_H__