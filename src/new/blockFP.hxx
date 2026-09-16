#ifndef __BlockFP_H__
#define __BlockFP_H__
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
   * \class BlockFP blockFP.hxx "blockFP.hxx"
   * \brief Polygon 2D Block Mesh definition.
   *
   *  This class defines the polygon 2D Block objects.
   *  It inherits from the Block class through a Factory.
   *
   */
  class BlockFP :
    public Block {

  public:

    /**
     * \brief Default constructor.
     *
     */
    BlockFP();


    /**
     * \brief Destructor.
     *
     */
    virtual ~BlockFP();


    /**
     * \brief Add a block to the mesh.
     *
     * \param [in]  block_type Type of the block
     * \param [in]  mesh       Mesh object owning the block
     */
    virtual void
    blockAdd
    (
      CWP_Block_t  block_type,
      Mesh        *mesh
    );


    /**
     * \brief Set a CWIPI block in a partition.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  n_elts     Number of elements of the block in the partition
     * \param [in]  connec_idx Elements connectivity index
     * \param [in]  connec     Elements connectivity
     * \param [in]  global_num Global element numbering
     */
    void
    blockSet
    (
      int          i_part,
      int          n_elts,
      int         *connec_idx,
      int         *connec,
      CWP_g_num_t *global_num
    );


    /**
     * \brief Get a CWIPI block in a partition.
     *
     * \param [in]  i_part     Partition identifier
     * \param [out] n_elts     Number of elements of the block in the partition
     * \param [out] connec_idx Elements connectivity index
     * \param [out] connec     Elements connectivity
     * \param [out] global_num Global element numbering
     */
    void
    blockGet
    (
      int           i_part,
      int          *n_elts,
      int         **connec_idx,
      int         **connec,
      CWP_g_num_t **global_num
    );


    /**
     * \brief Return the connectivity for a partition.
     *
     * \param [in] i_part Partition identifier
     * \return Connectivity array
     */
    inline int*
    connecGet
    (
      int i_part
    ) const;


    /**
     * \brief Return the connectivity index for a partition.
     *
     * \param [in] i_part Partition identifier
     * \return Connectivity index array
     */
    inline int*
    connecIdxGet
    (
      int i_part
    ) const;


    /**
     * \brief Finalize geometry processing.
     *
     */
    void
    geomFinalize();

  private:

    std::vector<int*> _connec_idx; /*!< Connectivity index for each partition */
    std::vector<int*> _connec;     /*!< Connectivity for each partition */

  };

  inline int*
  BlockFP::connecGet
  (
    int i_part
  ) const
  {
    return _connec[i_part];
  }


  inline int*
  BlockFP::connecIdxGet
  (
    int i_part
  ) const
  {
    return _connec_idx[i_part];
  }

}

/**
 * \endcond
 */

#endif // __BlockFP_H__