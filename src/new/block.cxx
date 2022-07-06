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

#include "pdm_mesh_nodal.h"
#include "pdm_gnum.h"
#include "block.hxx"
#include <map>
#include <vector>
#include "mesh.hxx"

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  Block::Block()
         :_isGNumSet(false),
          _inPDMDB(false)
  {

  }


  void Block::BlockAdd(CWP_Block_t blockType, Mesh* mesh)
  {
     _mesh     = mesh;
     _localComm            = const_cast<MPI_Comm*>(static_cast<Mesh*>(mesh)->getMPICommP());

     _blockType             = blockType;
     _n_part                = static_cast<Mesh*>(mesh)->getNPart();

     _global_num         .resize(_n_part,NULL);
     _n_elt              .resize(_n_part);
     _cells_center       .resize(_n_part, NULL);
     _isSet              .resize(_n_part);

  }


  Block::~Block(){
   for (int i = 0; i < _n_part; i++) {
      if (_cells_center[i] != NULL) {
        free (_cells_center[i]);
      }
    }

  }


  const double* Block::eltCentersGet(int i_part) {

    if (_cells_center[i_part] == NULL) {
      PDM_Mesh_nodal_cell_centers_compute (_mesh->getPdmNodalIndex(), _block_id_pdm, i_part, PDM_OWNERSHIP_KEEP);
    }
    return PDM_Mesh_cell_centers_get (_mesh->getPdmNodalIndex(), _block_id_pdm, i_part);
  }

  CWP_g_num_t*
  Block::GNumMeshGet(int i_part) {
    return _global_num[i_part];
  }

  void
  Block::GNumMeshSet(int i_part,CWP_g_num_t* gnum) {
    _global_num[i_part] = gnum;
  }


  PDM_Mesh_nodal_elt_t Block::PdmBlockTypeFromCwpBlockType(
                                                           CWP_Block_t CWP_block_type
                                                          )
   {
     PDM_Mesh_nodal_elt_t block_type;
     switch (CWP_block_type) {

       case CWP_BLOCK_NODE: 
       block_type = PDM_MESH_NODAL_POINT;
       break;

       case CWP_BLOCK_EDGE2: block_type = PDM_MESH_NODAL_BAR2;
       break;

       case CWP_BLOCK_FACE_TRIA3: block_type =  PDM_MESH_NODAL_TRIA3;
       break;

       case CWP_BLOCK_FACE_QUAD4: block_type = PDM_MESH_NODAL_QUAD4;
       break;

       case CWP_BLOCK_CELL_TETRA4: block_type = PDM_MESH_NODAL_TETRA4;
       break;

       case CWP_BLOCK_FACE_POLY: block_type = PDM_MESH_NODAL_POLY_2D;
       break;

       case CWP_BLOCK_CELL_HEXA8: block_type = PDM_MESH_NODAL_HEXA8;
       break;

       case CWP_BLOCK_CELL_PYRAM5: block_type = PDM_MESH_NODAL_PYRAMID5;
       break;

       case CWP_BLOCK_CELL_PRISM6: block_type = PDM_MESH_NODAL_PRISM6;
       break;

       case CWP_BLOCK_CELL_POLY: block_type = PDM_MESH_NODAL_POLY_3D;
       break;

       default:block_type = PDM_MESH_NODAL_POINT;
               PDM_error(__FILE__, __LINE__, 0, "No referenced CWP_Block_t.\n");
      }
    return block_type;

  }

  CWP_Block_t Block::CwpBlockTypeFromPdmBlockType(PDM_Mesh_nodal_elt_t PDM_block_type)
   {
     CWP_Block_t block_type;
     switch (PDM_block_type) {

       case PDM_MESH_NODAL_POINT: block_type = CWP_BLOCK_NODE;
       break;

       case PDM_MESH_NODAL_BAR2: block_type = CWP_BLOCK_EDGE2;
       break;

       case PDM_MESH_NODAL_TRIA3: block_type = CWP_BLOCK_FACE_TRIA3;
       break;

       case PDM_MESH_NODAL_QUAD4: block_type = CWP_BLOCK_FACE_QUAD4;
       break;

       case PDM_MESH_NODAL_TETRA4: block_type = CWP_BLOCK_CELL_TETRA4;
       break;

       case PDM_MESH_NODAL_POLY_2D: block_type = CWP_BLOCK_FACE_POLY;
       break;

       case PDM_MESH_NODAL_HEXA8: block_type = CWP_BLOCK_CELL_HEXA8;
       break;

       case PDM_MESH_NODAL_PYRAMID5: block_type = CWP_BLOCK_CELL_PYRAM5;
       break;

       case PDM_MESH_NODAL_PRISM6: block_type = CWP_BLOCK_CELL_PRISM6;
       break;

       case PDM_MESH_NODAL_POLY_3D: block_type = CWP_BLOCK_CELL_POLY;
       break;

       default: block_type = CWP_BLOCK_NODE;
                PDM_error(__FILE__, __LINE__, 0, "No referenced PDM_Mesh_nodal_elt_t.\n");
      }
    return block_type;

   }




}

/**
 * \endcond
 */
