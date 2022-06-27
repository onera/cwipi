#ifndef __BLOCK_H__
#define __BLOCK_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2013-2017  ONERA

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
#include "cwp_priv.h"

#include "pdm_mesh_nodal.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include <map>
#include <vector>
#include "pdm_error.h"


/**
 * \cond
 */

namespace cwipi {

  /**
   * \class Block block.hxx "block.hxx"
   * \brief Mesh mother Block definition
   *
   *  This class defines the Block objects.
   *  Blocks contain elements, are a part of
   *  a Mesh and are defined on Mesh partitions.
   *
   */

  class Block {

  public :

    /**
     *
     * \brief Constructor
     *
     */

    Block();

    /**
     * \brief Destructor.
     *
     */

    virtual ~Block();


    /**
     * \brief Block Type conversion from Paradigm to CWIPI
     *
     * Converts a Paradigm block type to a CWIPI block type.
     *
     * \param [in] PDM_block_type A paradigm block type.
     *
     * \return A CWIPI block type.
     *
     */

    CWP_Block_t          CwpBlockTypeFromPdmBlockType (PDM_Mesh_nodal_elt_t PDM_block_type);

    /**
     * \brief Block Type conversion from CWIPI to Paradigm
     *
     * Converts a CWIPI block type to a Paradigm block type.
     *
     * \param [in] CWP_block_type        A CWIPI block type.
     * \return A Paradigm block type.
     *
     */

    PDM_Mesh_nodal_elt_t PdmBlockTypeFromCwpBlockType (CWP_Block_t CWP_block_type);


    /**
     * \brief Add and Set the CWIPI block from a Paradigm block
     *
     * \param [in] pdm_id_block A block identifier from Paradigm
     * \param [in] mesh         Pointer to the Mesh object owning the block.
     *
     *
     */

    virtual void FromPDMBlock(int pdm_id_block, void* mesh)
    {
      CWP_UNUSED(pdm_id_block  );
      CWP_UNUSED(mesh );
      PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
    }

    /**
     * \brief Set a CWIPI block in a partition
     *
     * \param [in]  i_part          Partition identifier
     * \param [in]  n_elts          Number of elements of the block in the partition.
     * \param [in]  connec          Elements connectivity.
     * \param [in]  mesh_global_num Mesh global numbering of the block
     *
     */

    virtual void blockSet(int i_part,int n_elts,int* connec,CWP_g_num_t* mesh_global_num)
    {
      CWP_UNUSED(i_part          );
      CWP_UNUSED(n_elts          );
      CWP_UNUSED(connec          );
      CWP_UNUSED(mesh_global_num      );
      PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
    }

    /**
     * \brief Set a CWIPI block in a partition
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  n_elts     Number of elements of the block in the partition.
     * \param [in]  connec_idx Elements connectivity index
     * \param [in]  connec     Elements connectivity
     * \param [in]  global_num Mesh global numbering of the block
     *
     */

    virtual void blockSet(int i_part,int n_elts,
                          int* connec_idx,
                          int* connec,
                          CWP_g_num_t* global_num)
    {
      CWP_UNUSED(i_part          );
      CWP_UNUSED(n_elts          );
      CWP_UNUSED(connec_idx      );
      CWP_UNUSED(connec          );
      CWP_UNUSED(global_num      );
      PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
    }

    /**
     * \brief Set a CWIPI block in a partition
     *
     * \param [in]  i_part     Partition identifier
     * \param [in]  n_elts     Number of elements of the block in the partition
     * \param [in]  n_faces    Number of faces of the block in the partition
     * \param [in]  connec_faces_idx Vertices to faces connectivity index
     * \param [in]  connec_faces     Vertices to faces connectivity
     * \param [in]  connec_cells_idx Faces to cells connectivity index
     * \param [in]  connec_cells     Faces to cells connectivity
     * \param [in]  global_num Mesh  Global numbering of the block
     *
     */

    virtual void blockSet(int i_part,int n_elts,int n_faces,
                          int* connec_faces_idx,
                          int* connec_faces,
                          int* connec_cells_idx,
                          int* connec_cells,
                          CWP_g_num_t* global_num)
    {
      CWP_UNUSED(i_part          );
      CWP_UNUSED(n_elts          );
      CWP_UNUSED(n_faces         );
      CWP_UNUSED(connec_faces_idx);
      CWP_UNUSED(connec_faces    );
      CWP_UNUSED(connec_cells_idx);
      CWP_UNUSED(connec_cells    );
      CWP_UNUSED(global_num      );
      PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
    }

    /**
     *
     * \brief Block addition
     *
     * Add a block to the mesh.
     *
     * \param [in] blockType              Type of the block
     * \param [in] mesh                   The Mesh object owning the block
     */

    void BlockAdd(CWP_Block_t blockType, void* mesh);

    /**
     *
     * \brief Return the CWIPI block type.
     *
     * \return The CWIPI block type.
     */

    inline CWP_Block_t blockTypeGet();


    /**
     *
     * \brief Return the block identifier.
     *
     */

    inline int         blockIDGet();
    inline int         blockIDCWIPIGet();

    inline void blockIDCWIPISet(int block_id);

    inline void blockIDPDMSet(int block_id);


    inline int blockIDPDMGet();
    /**
     *
     * \brief return the number of element for each partition
     *
     *
     */

     inline std::vector <int>  NEltsGet();

     inline int  NEltsGet(int i_part);

    /**
     *
     * \brief return the element global numbering for each partition
     *
     *
     */

    CWP_g_num_t* GNumMeshGet(int i_part);

    void GNumMeshSet(int i_part, CWP_g_num_t* gnum);

    CWP_g_num_t* GNumBlockGet(int i_part);

    const double* eltCentersGet(int i_part);

    /**
     *
     * \brief return the number of partition
     *
     *
     */

     inline int NPartGet();


    /**
     *
     * \brief return true is the global numbering (inside the mesh) is defined
     *
     *
     */

    inline bool globNumDefined();


    /**
     *
     * \brief return true is the Block is set on all the partition (ready for global numbering computation).
     *
     * It returns true if the block is set on all the partition (ready for global numbering computation).
     *
     */

    inline bool isSet();

    /**
     *
     * \brief Return true if the Block is already in Paradigm block database and false otherwise.
     *
     */

    inline bool inPDMDB();


    /**
     *
     * \brief Set the _inPDMDB variable to true, and so indicate that the block is already in the PDM database.
     *
     */

    inline void SetinPDMDB();

    /**
     *
     * \brief return the element connectivity (Standard or Face_Poly_2D CWP_Block_t) or cells-faces connectivity (Cells_POLY_3D)
     * for each partition.
     *
     *
     */

     virtual std::map<int,int*>  ConnecGet()
     {
       PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
       std::map<int,int*> null;
       return null;
     }

     virtual int*  ConnecGet(int i_part)
     {
       CWP_UNUSED (i_part);
       PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
       return NULL;
     }


     virtual void geomFinalize(int already_in_pdm)
     {
       CWP_UNUSED (already_in_pdm);      
       PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
     }


     virtual void GNumBlockSet(int i_part, CWP_g_num_t* global_num)
     {
       CWP_UNUSED (i_part);
       CWP_UNUSED (global_num);
       PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
     }

    /**
     *
     * \brief return the element connectivity index (Face_Poly_2D CWP_Block_t) or cells-faces connectivity index (Cells_POLY_3D)
     * for each partition.
     *
     *
     */

     virtual std::map<int,int*>  ConnecIDXGet()
     {
       PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
       std::map<int,int*> null;
       return null;
     }

     virtual int*  ConnecIDXGet(int i_part)
     {
       CWP_UNUSED (i_part);
       PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
       return NULL;
     }


    /**
     *
     * \brief return the element faces connectivity (Cells_POLY_3D)
     * for each partition.
     *
     *
     */

     virtual std::map<int,int*>  ConnecFacesGet()
     {
       PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
       std::map<int,int*> null;
       return null;
     }

    /**
     *
     * \brief return the element faces connectivity index (Cells_POLY_3D)
     * for each partition.
     *
     *
     */

     virtual std::map<int,int*>  ConnecFacesIDXGet()
     {
       PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
       std::map<int,int*> null;
       return null;
     }



    /**
     *
     * \brief return the number of faces for each partition (Cells_POLY_3D)
     *
     *
     */

     virtual std::map <int,int>  NFacesGet()
     {
       PDM_error(__FILE__, __LINE__, 0, "This function is not available for this type of Block.\n");
       std::map<int,int> null;
       return null;
     }


  private :

    /**
     *
     * \brief Assigment operator
     *
     */

    Block
    &operator=
    (const Block &other);

    /**
     *
     * \brief Copy constructor
     *
     */

    Block (const Block& other);

  protected:

    CWP_Block_t                _blockType;              /*!< Block Type */
    PDM_Mesh_nodal_t          *_pdmNodal_handle_index;  /*!< PDM Nodal Index */
    MPI_Comm                  *_localComm;              /*!< Communicator */
    int                        _n_part;                 /*!< Number of partitions */
    int                        _n_part_def;             /*!< Number of partitions where the block is defined */
    std::vector <int>          _n_elt;                  /*!< Number of elements for each partition */
    std::vector <double*>      _cells_center;     /*!< Cell centers */
    std::vector <int>          _part_id;                /*!< Partition where the block is defined  */
    int                        _block_id_pdm;           /*!< Block identifier */
    int                        _block_id_cwipi;         /*!< Block identifier */
    std::vector <CWP_g_num_t*> _global_num;             /*!< Global numbering in the Mesh  */
    std::vector <CWP_g_num_t*> _global_num_computed;    /*!< Global numbering computed in the Mesh  */
    std::vector <CWP_g_num_t*> _global_num_block;       /*!< Global numbering in the Block */
    void                      *_mesh;                   /*!< Pointer to the mesh object owning the block */
    std::vector<bool>          _isSet;                  /*!< Set or not for each partition */
    bool                       _isGNumSet;              /*!< Global Numbering set or not for each partition */
    bool                       _inPDMDB;                /*!< Indicate the Block is already in the Paradigm database */

  };

  bool
  Block::globNumDefined() {
    return _isGNumSet;
  }

  std::vector <int>
  Block::NEltsGet() {
    return _n_elt;
  }


  int
  Block::NEltsGet(int i_part) {
    return _n_elt[i_part];
  }

  int
  Block::NPartGet() {
    return _n_part;
  }



  bool Block::inPDMDB() {
    return _inPDMDB;
  }

  void Block::SetinPDMDB() {
    _inPDMDB=true;
  }


  bool
  Block:: isSet() {
    for(int i=0;i<_n_part;i++) {
      if(_isSet[i] == false)
        return false;
    }
    return true;
  }

  CWP_Block_t
  Block::blockTypeGet(){
    return _blockType;
  }


  int
  Block::blockIDGet(){
    return _block_id_pdm;
  }

  int
  Block::blockIDCWIPIGet(){
    return _block_id_cwipi;
  }

  void
  Block::blockIDCWIPISet(int block_id){
    _block_id_cwipi = block_id;
  }

  void
  Block::blockIDPDMSet(int block_id){
    _block_id_pdm = block_id;
  }

  int 
  Block::blockIDPDMGet(){
    return _block_id_pdm;
  }

}
#endif //__BLOCK_H__

/**
 * \endcond
 */
