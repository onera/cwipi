#ifndef __BLOCKHO_H__
#define __BLOCKHO_H__
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
#include "blockStd.hxx"

#include "cwp.h"

#include <vector>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class BlockHO blockHO.hxx "blockHO.hxx"
   * \brief Higher-Order Block Mesh definition.
   *
   *  This class defines the higher-order Block objects.
   *  It inherits from the BlockStd class.
   *
   */
  class BlockHO :
    public BlockStd {

  public:

    /**
     * \brief Default constructor.
     *
     */
    BlockHO();


    /**
     * \brief Destructor.
     *
     */
    virtual ~BlockHO();


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
    ) override;


    /**
     * \brief Set a CWIPI block in a partition.
     *
     * \param [in]  i_part      Partition identifier
     * \param [in]  n_elts      Number of elements of the block
     * \param [in]  connec      Element to vertex connectivity
     * \param [in]  global_num  Global numbering of the block
     * \param [in]  order       Element order
     * \param [in]  ho_ordering HO ordering name
     */
    void
    blockSet
    (
            int          i_part,
            int          n_elts,
            int         *connec,
            CWP_g_num_t *global_num,
      const int          order,
      const char        *ho_ordering
    );


    /**
     * \brief Get a CWIPI block.
     *
     * \param [in]  i_part      Partition identifier
     * \param [out] n_elts      Number of elements of the block
     * \param [out] connec      Element to vertex connectivity
     * \param [out] global_num  Global numbering of the block
     * \param [out] order       Element order
     * \param [out] ho_ordering HO ordering name
     */
    void
    blockGet
    (
      int           i_part,
      int          *n_elts,
      int         **connec,
      CWP_g_num_t **global_num,
      int          *order,
      char        **ho_ordering
    );


    /**
     * \brief Return the connectivity for a partition.
     *
     * \param [in]  i_part Partition identifier
     * \return Connectivity array
     */
    inline int*
    connecGet
    (
      int i_part
    ) const;


    /**
     * \brief Return the IJK connectivity for a partition.
     *
     * \param [in]  i_part Partition identifier
     * \return IJK connectivity array
     */
    inline int*
    connecIJKGet
    (
      int i_part
    ) const;


    /**
     * \brief Return the element order.
     *
     * \return Element order
     */
    inline int
    orderGet() const;


    /**
     * \brief Return the higher-order ordering scheme.
     *
     * \return High-order ordering name
     */
    inline char*
    hoOrderingGet() const;


    /**
     * \brief Finalize geometry processing.
     *
     */
    void
    geomFinalize() override;

  private:

    /**
     * \brief Reorder the connectivity.
     *
     */
    void
    _reorder();


    int                _order;       /*!< Element order */
    char              *_ho_ordering; /*!< HO ordering name */
    std::vector<int*>  _connec;      /*!< Connectivity in user order for each partition */
    std::vector<int*>  _connec_ijk;  /*!< Connectivity in Lagrange (IJK) order for each partition */

  };

  inline int*
  BlockHO::connecGet
  (
    int i_part
  ) const
  {
    return _connec[i_part];
  }


  inline int*
  BlockHO::connecIJKGet
  (
    int i_part
  ) const
  {
    return _connec_ijk[i_part];
  }


  inline int
  BlockHO::orderGet() const
  {
    return _order;
  }


  inline char*
  BlockHO::hoOrderingGet() const
  {
    return _ho_ordering;
  }

}

/**
 * \endcond
 */

#endif // __BLOCKHO_H__