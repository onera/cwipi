#ifndef __BlockHO_H__
#define __BlockHO_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

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

#include <cstdlib>
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
    ~BlockHO();


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
     * \param [in]  i_part      Partition identifier
     * \param [in]  n_elts      Number of elements of the block in the partition
     * \param [in]  connec      Elements connectivity
     * \param [in]  global_num  Global element numbering
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
     * \param [in]  i_part        Partition identifier
     * \param [out] n_elts        Number of elements of the block in the partition
     * \param [out] connec        Elements connectivity
     * \param [out] global_num    Global element numbering
     * \param [out] order         Element order
     * \param [out] ho_ordering   HO ordering name
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
     * \brief Return the element connectivity for each partition.
     *
     * \return Vector of connectivity arrays for each partition
     */
    inline std::vector<int*>
    ConnecGet();


    /**
     * \brief Return the IJK connectivity for each partition.
     *
     * \return Vector of IJK connectivity arrays for each partition
     */
    inline std::vector<int*>
    ConnecIJKGet();


    /**
     * \brief Return the connectivity for a partition.
     *
     * \param [in] i_part Partition identifier
     * \return Connectivity array
     */
    inline int*
    ConnecGet
    (
      int i_part
    );


    /**
     * \brief Return the IJK connectivity for a partition.
     *
     * \param [in] i_part Partition identifier
     * \return IJK connectivity array
     */
    inline int*
    ConnecIJKGet
    (
      int i_part
    );


    /**
     * \brief Return the element order.
     *
     * \return Element order
     */
    inline int
    OrderGet();


    /**
     * \brief Return the higher-order ordering scheme.
     *
     * \return Higher-order ordering string
     */
    inline char*
    HOOrderingGet();


    /**
     * \brief Finalize geometry processing.
     *
     */
    void
    geomFinalize();


    /**
     * \brief Reorder the connectivity.
     *
     */
    void
    reorder();

  private:

    int                _order;            /*!< Element order */
    char              *_ho_ordering;      /*!< HO ordering name */
    std::vector<int*>  _connec;           /*!< Connectivity in user order for each partition */
    std::vector<int*>  _connec_ijk;       /*!< Connectivity in Lagrange (IJK) order for each partition */

  };

  inline int*
  BlockHO::ConnecGet
  (
    int i_part
  )
  {
    return _connec[i_part];
  }


  inline std::vector<int*>
  BlockHO::ConnecGet()
  {
    return _connec;
  }


  inline int*
  BlockHO::ConnecIJKGet
  (
    int i_part
  )
  {
    return _connec_ijk[i_part];
  }


  inline std::vector<int*>
  BlockHO::ConnecIJKGet()
  {
    return _connec_ijk;
  }


  inline int
  BlockHO::OrderGet()
  {
    return _order;
  }


  inline char*
  BlockHO::HOOrderingGet()
  {
    return _ho_ordering;
  }

}

/**
 * \endcond
 */

#endif // __BlockHO_H__