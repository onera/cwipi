#ifndef __Block_H__
#define __Block_H__
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

#include "cwp_priv.h"
CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wcast-function-type")
#include "cwp.h"
CWP_GCC_SUPPRESS_WARNING_POP
#include "pdm_error.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include <map>
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
     * \brief Return the CWIPI block type.
     *
     * \return The CWIPI block type
     */
    inline CWP_Block_t
    blockTypeGet();


    /**
     * \brief Return the block identifier.
     *
     * \return Block identifier
     */
    inline int
    blockIDCWIPIGet();


    /**
     * \brief Set the CWIPI block identifier.
     *
     * \param [in] block_id Block identifier
     */
    inline void
    blockIDCWIPISet
    (
      int block_id
    );


    /**
     * \brief Set the PDM block identifier.
     *
     * \param [in] block_id Block identifier
     */
    inline void
    blockIDPDMSet
    (
      int block_id
    );


    /**
     * \brief Return the PDM block identifier.
     *
     * \return Block identifier
     */
    inline int
    blockIDPDMGet();


    /**
     * \brief Return the number of elements for a partition.
     *
     * \param [in] i_part Partition identifier
     * \return Number of elements
     */
    inline int
    NEltsGet
    (
      int i_part
    );


    /**
     * \brief Return the element global numbering for a partition.
     *
     * \param [in] i_part Partition identifier
     * \return Global numbering array
     */
    inline CWP_g_num_t*
    GNumMeshGet
    (
      int i_part
    );


    /**
     * \brief Set the element global numbering.
     *
     * \param [in] i_part  Partition identifier
     * \param [in] gnum    Global numbering array
     * \param [in] owner   Ownership flag
     */
    void
    GNumMeshSet
    (
      int              i_part,
      CWP_g_num_t     *gnum,
      PDM_ownership_t  owner
    );


    /**
     * \brief Free the element global numbering.
     *
     * \param [in] i_part Partition identifier
     */
    void
    GNumMeshFree
    (
      int i_part
    );


    /**
     * \brief Return the element centers.
     *
     * \param [in] i_part Partition identifier
     * \return Element centers array
     */
    const double*
    eltCentersGet
    (
      int i_part
    );


    /**
     * \brief Set the parent numbering.
     *
     * \param [in] i_part     Partition identifier
     * \param [in] parent_num Parent numbering array
     * \param [in] owner      Ownership flag
     */
    void
    blockSetParentNum
    (
      int             i_part,
      int            *parent_num,
      PDM_ownership_t owner
    );


    /**
     * \brief Return the parent numbering.
     *
     * \param [in] i_part Partition identifier
     * \return Parent numbering array
     */
    inline int*
    ParentNumGet
    (
      int i_part
    );


    /**
     * \brief Free the parent numbering.
     *
     * \param [in] i_part Partition identifier
     */
    void
    ParentNumFree
    (
      int i_part
    );


    /**
     * \brief Finalize geometry processing.
     *
     *  Pure virtual method to be implemented by derived classes.
     */
    virtual void
    geomFinalize() = 0;

  private:

    /**
     * \brief Assignment operator.
     *
     */
    Block& operator=
    (
      const Block &other
    );


    /**
     * \brief Copy constructor.
     *
     */
    Block
    (
      const Block& other
    );


  protected:

    CWP_Block_t                 _blockType;          /*!< Block type */
    MPI_Comm                   *_localComm;          /*!< Communicator */
    int                         _n_part;             /*!< Number of partitions */
    std::vector<int>            _n_elt;              /*!< Number of elements for each partition */
    std::vector<double*>        _cells_center;       /*!< Cell centers */
    int                         _block_id_pdm;       /*!< Block identifier in PDM */
    int                         _block_id_cwipi;     /*!< Block identifier in CWIPI */
    std::vector<CWP_g_num_t*>   _global_num;         /*!< Global numbering in the mesh */
    std::vector<int*>           _parent_num;         /*!< Parent numbering */
    Mesh                       *_mesh;               /*!< Pointer to the mesh object owning the block */
    PDM_ownership_t             _owner_gnum;         /*!< Owner of global numbers */
    PDM_ownership_t             _owner_parent_num;   /*!< Owner of parent numbers */

  };

  inline CWP_Block_t
  Block::blockTypeGet()
  {
    return _blockType;
  }


  inline int
  Block::blockIDCWIPIGet()
  {
    return _block_id_cwipi;
  }


  inline void
  Block::blockIDCWIPISet
  (
    int block_id
  )
  {
    _block_id_cwipi = block_id;
  }


  inline void
  Block::blockIDPDMSet
  (
    int block_id
  )
  {
    _block_id_pdm = block_id;
  }


  inline int
  Block::blockIDPDMGet()
  {
    return _block_id_pdm;
  }


  inline int
  Block::NEltsGet
  (
    int i_part
  )
  {
    return _n_elt[i_part];
  }


  inline CWP_g_num_t*
  Block::GNumMeshGet
  (
    int i_part
  )
  {
    return _global_num[i_part];
  }


  inline int*
  Block::ParentNumGet
  (
    int i_part
  )
  {
    return _parent_num[i_part];
  }

}

/**
 * \endcond
 */

#endif // __Block_H__