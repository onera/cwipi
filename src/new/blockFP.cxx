

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
        _block_id = pdm_id_block;         
        int* connec     = NULL;       
        int* connec_idx = NULL;                                         
        PDM_Mesh_nodal_block_poly2d_get  (_pdmNodal_handle_index,
                                          _block_id,
                                          id_part,
                                          &connec_idx,
                                          &connec);
                                       
         PDM_Mesh_nodal_elt_t PDM_block_type = PDM_Mesh_nodal_block_type_get(_pdmNodal_handle_index,
                                                                              pdm_id_block);

         blockSet(id_part,nElts,connec_idx,connec,NULL);         
       }

  }


  void BlockFP::blockSet(int i_part,int n_elt,int* connec_idx,int* connec,CWP_g_num_t* global_num){

     double* _cells_center_part = (double*)malloc (sizeof(double) * 3 * n_elt);
     double* coords = static_cast<Mesh*>(_mesh) -> getCoordinates(i_part);

     int n_vtx=-1;

     for (int i_dim = 0; i_dim < 3; i_dim++) {
       int ind_vtx=0;
       for (int i_elt =0; i_elt < n_elt; i_elt++) {
         _cells_center_part[3*i_elt+i_dim] = 0.0;
         n_vtx = connec_idx[i_elt+1] - connec_idx[i_elt];
         for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
           _cells_center_part[3*i_elt+i_dim] += coords[3*connec[ ind_vtx ] + i_dim];
           ind_vtx++; 
         } //i_vtx
         _cells_center_part[3*i_elt+i_dim] /= double(n_vtx);
       } //i_elts
     } //i_dim

     PDM_gnum_set_from_coords (_pdmGNum_handle_index, i_part, n_elt, _cells_center_part, NULL);
     
     _isSet[i_part] = true;
     _n_elt[i_part] = n_elt;
     _part_id.push_back(i_part);
     _n_part_def=_n_part_def+1;
     _cells_center[i_part] = _cells_center_part;
     _connec_idx.insert( std::pair < int, int* >         (i_part,connec_idx) );  
     _connec.insert    ( std::pair < int, int* >         (i_part,connec)     );
                               

                              
     if( isSet() ) {
       PDM_gnum_compute (_pdmGNum_handle_index);
       for (int i = 0;i<_n_part;i++) {

         _global_num[i] = const_cast<CWP_g_num_t*> (PDM_gnum_get (_pdmGNum_handle_index, i));
         if (not inPDMDB() ) PDM_Mesh_nodal_block_poly2d_set (_pdmNodal_handle_index,
                                          _block_id,
                                          i,    
                                          _n_elt[i],  
                                          _connec_idx[i], 
                                          _connec[i],   
                                          _global_num[i],
                                          NULL);  
       } //i
       PDM_Mesh_nodal_g_num_in_block_compute(_pdmNodal_handle_index,_block_id);
       
       for (int i = 0;i<_n_part;i++) {
          _global_num_block[i] = PDM_Mesh_nodal_block_g_num_get (_pdmNodal_handle_index,
                                                                 _block_id,
                                                                 i );
       } // i
       
       _isGNumSet = true;
       
     }//isSet()
  }
}

