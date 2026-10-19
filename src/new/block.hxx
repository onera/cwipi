#ifndef __BLOCK_H__
#define __BLOCK_H__
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

#include "cwp_priv.h"

CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wcast-function-type")
#include "cwp.h"
CWP_GCC_SUPPRESS_WARNING_POP

#include "pdm_gnum.h"

#include <vector>

/**
 * \cond
 */

namespace cwipi {

  class Mesh;

  /**
   * \class Block block.hxx "block.hxx"
   * \brief Mesh mother Block definition.
   *
   *  This class defines the Block objects.
   *  Blocks contain elements, are a part of
   *  a Mesh and are defined on Mesh partitions.
   *
   */
  class Block {

  public:

    /**
     * \brief Default constructor.
     *
     */
    Block();


    /**
     * \brief Destructor.
     *
     */
    virtual ~Block();


    /**
     * \brief Add a block to the mesh.
     *
     * \param [in]  block_type Type of the block
     * \param [in]  mesh       Mesh object owning the block
     */
    virtual void
    blockAdd
    (
      const CWP_Block_t  block_type,
            Mesh        *mesh
    );


    /**
     * \brief Return the CWIPI block type.
     *
     * \return Block type
     */
    inline CWP_Block_t
    blockTypeGet() const;


    /**
     * \brief Return true if the block is standard.
     *
     * \return True if the block is standard
     */
    inline bool
    blockTypeStdIs() const;


    /**
     * \brief Set the CWIPI block identifier.
     *
     * \param [in]  block_id Block identifier
     */
    inline void
    blockIdCwipiSet
    (
      const int block_id
    );


    /**
     * \brief Return the CWIPI block identifier.
     *
     * \return Block identifier
     */
    inline int
    blockIdCwipiGet() const;


    /**
     * \brief Set the PDM block identifier.
     *
     * \param [in]  block_id Block identifier
     */
    inline void
    blockIdPdmSet
    (
      const int block_id
    );


    /**
     * \brief Return the PDM block identifier.
     *
     * \return Block identifier
     */
    inline int
    blockIdPdmGet() const;


    /**
     * \brief Return the number of elements for a partition.
     *
     * \param [in]  i_part Partition identifier
     * \return Number of elements
     */
    inline int
    nEltsGet
    (
      const int i_part
    ) const;


    /**
     * \brief Set the element global numbering.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  global_num Global numbering array
     * \param [in]  owner      Ownership flag
     *
     * \note If owner is PDM_OWNERSHIP_KEEP, the Block takes ownership
     *       and will free the memory when the block is destroyed.
     */
    void
    globalNumSet
    (
      const int              i_part,
            CWP_g_num_t     *global_num,
      const PDM_ownership_t  owner
    );


    /**
     * \brief Return the element global numbering for a partition.
     *
     * \param [in]  i_part Partition identifier
     * \return Global numbering array
     */
    inline CWP_g_num_t*
    globalNumGet
    (
      const int i_part
    ) const;


    /**
     * \brief Free the element global numbering.
     *
     * \param [in]  i_part Partition identifier
     */
    void
    globalNumFree
    (
      const int i_part
    );


    /**
     * \brief Set the parent numbering.
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  parent_num Parent numbering array
     * \param [in]  owner      Ownership flag
     *
     * \note If owner is PDM_OWNERSHIP_KEEP, the Block takes ownership
     *       and will free the memory when the block is destroyed.
     */
    void
    parentNumSet
    (
      const int              i_part,
            int             *parent_num,
      const PDM_ownership_t  owner
    );


    /**
     * \brief Return the parent numbering.
     *
     * \param [in]  i_part Partition identifier
     * \return Parent numbering array
     */
    inline int*
    parentNumGet
    (
      const int i_part
    ) const;


    /**
     * \brief Free the parent numbering.
     *
     * \param [in]  i_part Partition identifier
     */
    void
    parentNumFree
    (
      const int i_part
    );


    /**
     * \brief Return the element centers.
     *
     * \param [in]  i_part Partition identifier
     * \return Element centers array
     */
    double*
    elementCentersGet
    (
      const int i_part
    );


    /**
     * \brief Finalize geometry processing.
     *
     *  Pure virtual method to be implemented by derived classes.
     */
    virtual void
    geomFinalize() = 0;

  protected:

    CWP_Block_t                _block_type;       /*!< Block type */
    int                        _block_id_pdm;     /*!< Block identifier in PDM */
    int                        _block_id_cwipi;   /*!< Block identifier in CWIPI */
    Mesh                      *_mesh;             /*!< Pointer to the mesh object owning the block */
    int                        _n_part;           /*!< Number of partitions */
    std::vector<int>           _n_elts;           /*!< Number of elements for each partition */
    std::vector<double*>       _elt_centers;      /*!< Element centers */
    std::vector<CWP_g_num_t*>  _global_num;       /*!< Global numbering in the mesh */
    std::vector<int*>          _parent_num;       /*!< Parent numbering */
    PDM_ownership_t            _global_num_owner; /*!< Owner of global numbers */
    PDM_ownership_t            _parent_num_owner; /*!< Owner of parent numbers */

  };

  inline CWP_Block_t
  Block::blockTypeGet() const
  {
    return _block_type;
  }


  inline bool
  Block::blockTypeStdIs() const
  {
    return (_block_type == CWP_BLOCK_NODE  ||
            _block_type == CWP_BLOCK_EDGE2       ||
            _block_type == CWP_BLOCK_FACE_TRIA3  ||
            _block_type == CWP_BLOCK_FACE_QUAD4  ||
            _block_type == CWP_BLOCK_CELL_TETRA4 ||
            _block_type == CWP_BLOCK_CELL_HEXA8  ||
            _block_type == CWP_BLOCK_CELL_PRISM6 ||
            _block_type == CWP_BLOCK_CELL_PYRAM5);
  }


  inline void
  Block::blockIdCwipiSet
  (
    const int block_id
  )
  {
    _block_id_cwipi = block_id;
  }


  inline int
  Block::blockIdCwipiGet() const
  {
    return _block_id_cwipi;
  }


  inline void
  Block::blockIdPdmSet
  (
    const int block_id
  )
  {
    _block_id_pdm = block_id;
  }


  inline int
  Block::blockIdPdmGet() const
  {
    return _block_id_pdm;
  }


  inline int
  Block::nEltsGet
  (
    const int i_part
  ) const
  {
    return _n_elts[i_part];
  }


  inline CWP_g_num_t*
  Block::globalNumGet
  (
    const int i_part
  ) const
  {
    return _global_num[i_part];
  }


  inline int*
  Block::parentNumGet
  (
    const int i_part
  ) const
  {
    return _parent_num[i_part];
  }

}

/**
 * \endcond
 */

#endif // __BLOCK_H__