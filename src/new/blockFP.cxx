

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

#include "blockFP.hxx"
#include "cwp.h"
#include <pdm_mesh_nodal.h>
#include <pdm_gnum.h>
#include <bftc_error.h>
#include <bftc_printf.h>
#include <map>
#include <vector>
#include <mesh.hxx>

/**
 * \cond
 */

namespace cwipi {
  BlockFP::BlockFP()
     :Block::Block()
  {

  }


  BlockFP::~BlockFP()
  {

  }

  void BlockFP::FromPDMBlock(int pdm_id_block, void* mesh){

     _mesh     = mesh;
     _pdmNodal_handle_index = static_cast<Mesh*>(mesh) -> getPdmNodalIndex();
     _localComm            = const_cast<MPI_Comm*>(static_cast<Mesh*>(mesh) -> getMPICommP());

     PDM_Mesh_nodal_elt_t PDM_block_type = PDM_Mesh_nodal_block_type_get(_pdmNodal_handle_index,pdm_id_block);
     _blockType = CwpBlockTypeFromPdmBlockType (PDM_block_type);
     BlockAdd(_blockType, mesh);

     for(int id_part=0;id_part < _n_part;id_part++){
        int nElts = PDM_Mesh_nodal_block_n_elt_get(_pdmNodal_handle_index,
                                                   pdm_id_block,
                                                   id_part );
        _block_id_pdm = pdm_id_block;
        int* connec     = NULL;
        int* connec_idx = NULL;
        PDM_Mesh_nodal_block_poly2d_get  (_pdmNodal_handle_index,
                                          _block_id_pdm,
                                          id_part,
                                          &connec_idx,
                                          &connec);

         blockSet(id_part,nElts,connec_idx,connec,NULL);
       }

  }


  void BlockFP::blockSet(int i_part,
                         int n_elt,
                         int* connec_idx,
                         int* connec,
                         CWP_g_num_t* mesh_global_num){

     double* _cells_center_part = (double*)malloc (sizeof(double) * 3 * n_elt);

     _global_num [i_part] = mesh_global_num;

     _isSet[i_part] = true;
     _n_elt[i_part] = n_elt;
     _part_id.push_back(i_part);
     _n_part_def++;
     _cells_center[i_part] = _cells_center_part;
     _connec     .insert    ( std::pair < int, int* > (i_part,connec));
     _connec_idx.insert     ( std::pair < int, int* > (i_part,connec_idx));

  }

  void BlockFP::geomFinalize(int already_in_pdm){

     _pdmNodal_handle_index = static_cast<Mesh*>(_mesh) -> getPdmNodalIndex();

     if(already_in_pdm ==0)
      _block_id_pdm = PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                                PDM_FALSE,
                                                PdmBlockTypeFromCwpBlockType(_blockType));

    for(int i_part = 0; i_part<_n_part; i_part++){


       if(already_in_pdm ==0)
         PDM_Mesh_nodal_block_poly2d_set (_pdmNodal_handle_index,
                                          _block_id_pdm,
                                          i_part,
                                          _n_elt     [i_part],
                                          _connec_idx[i_part],
                                          _connec    [i_part],
                                          _global_num[i_part],
                                          NULL);

      Visu* visu = ((Mesh*)_mesh) -> getVisu();
      if(visu -> isCreated() && ((Mesh*)_mesh) -> getDisplacement() == CWP_DYNAMIC_MESH_STATIC) {
        visu -> GeomBlockPoly2D ( ((Mesh*)_mesh) -> getIdVisu( _block_id_cwipi ),
                                  i_part,
                                  _n_elt[i_part] ,
                                  _connec_idx [i_part],
                                  _connec     [i_part],
                                  _global_num [i_part]);


     }
    } //end i_part


  }
}


/**
 * \endcond
 */
