

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

#include "blockStd.hxx"
#include "cwp.h"
#include <pdm_mesh_nodal.h>
#include <pdm_gnum.h>
#include <bftc_error.h>
#include <bftc_printf.h>
#include <map>
#include <vector>

#include <mesh.hxx>

namespace cwipi {
  BlockStd::BlockStd()
     :Block::Block()
  {
  
  }


  BlockStd::~BlockStd()
  {
  
  }

  void BlockStd::FromPDMBlock(int pdm_id_block, void* mesh){
  
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
        _block_id=pdm_id_block;         
        int* connec = NULL;      
                                
        PDM_Mesh_nodal_block_std_get (_pdmNodal_handle_index,
                                      _block_id,
                                      id_part,
                                      &connec); 
                                      
         blockSet(id_part,nElts,connec,NULL);
                
       }
  }


  void BlockStd::blockSet(int i_part,int n_elt,int* connec,CWP_g_num_t* mesh_global_num){

     double* _cells_center_part = (double*)malloc (sizeof(double) * 3 * n_elt);

     if (mesh_global_num != NULL) {
       _global_num [i_part] = mesh_global_num;
     }
     else {
       _global_num [i_part] = _global_num_computed[i_part];
     }

     _isSet[i_part] = true;
     _n_elt[i_part] = n_elt;
     _part_id.push_back(i_part);
     _n_part_def=_n_part_def+1;
     _cells_center[i_part] = _cells_center_part;
     _connec.insert    ( std::pair < int, int* > (i_part,connec));

     int* _blocks_id = PDM_Mesh_nodal_blocks_id_get(_pdmNodal_handle_index);

     if (not inPDMDB() ) PDM_Mesh_nodal_block_std_set(_pdmNodal_handle_index,
                                                           _block_id,
                                                           i_part,    
                                                           n_elt,    
                                                           connec,   
                                                           _global_num [i_part],
                                                           NULL);   
                                                           
   // printf("_block_id %i PDM_Mesh_nodal_n_cell_get %i\n",_block_id,PDM_Mesh_nodal_n_cell_get(_pdmNodal_handle_index,i_part));

  }

}









