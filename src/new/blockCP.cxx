

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

#include "blockCP.hxx"
#include "cwp.h"
#include <pdm_mesh_nodal.h>
#include <pdm_gnum.h>
#include <bftc_error.h>
#include <bftc_printf.h>
#include <map>
#include <vector>
#include <mesh.hxx>
#include <pdm_priv.h>

namespace cwipi {
  BlockCP::BlockCP()
     :Block::Block()
  {
  
  }


  BlockCP::~BlockCP()
  {
  
  }

  void BlockCP::FromPDMBlock(int pdm_id_block, void* mesh){
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
        int* connec_faces     = NULL;       
        int* connec_faces_idx = NULL;        
        int* connec_cells     = NULL;       
        int* connec_cells_idx = NULL;     
        int  n_faces = 0;      
                     
        PDM_Mesh_nodal_block_poly3d_get  (_pdmNodal_handle_index,
                                          _block_id,
                                          id_part,
                                          &n_faces,
                                          &connec_faces_idx,
                                          &connec_faces,
                                          &connec_cells_idx,
                                          &connec_cells);

         blockSet(id_part,nElts,n_faces, connec_faces_idx, 
                  connec_faces,connec_cells_idx,connec_cells,
                  NULL);
           
       }
  }


  void BlockCP::blockSet(int i_part,int n_elts,int n_faces,
                        int* connec_faces_idx, 
                        int* connec_faces,
                        int* connec_cells_idx,
                        int* connec_cells,
                        CWP_g_num_t* global_num) {

     UNUSED(global_num );
     double* _cells_center_part = (double*)malloc (sizeof(double) * 3 * n_elts);
     double* faces_center_part = (double*)malloc (sizeof(double) * 3 * n_faces);
     double* coords = static_cast<Mesh*>(_mesh) -> getCoordinates(i_part);

     int n_vtx=-1;

     
     for (int i_dim = 0; i_dim < 3; i_dim++) {
       int ind_vtx=0; 
       for (int i_faces =0; i_faces < n_faces; i_faces++) {
         faces_center_part[3*i_faces+i_dim] = 0.0;
         n_vtx = connec_faces_idx[i_faces+1] - connec_faces_idx[i_faces];
         for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
           faces_center_part[3*i_faces+i_dim] += coords[3*connec_faces[ ind_vtx ] + i_dim];
           ind_vtx++;  
         } //i_vtx
         faces_center_part[3*i_faces+i_dim] /= double(n_vtx);
       } //i_faces
     } //i_dim
     int nelt_faces=-1;

     for (int i_dim = 0; i_dim < 3; i_dim++) {    
       int ind_faces=0;     
       for (int i_elt =0; i_elt < n_elts; i_elt++) {
         _cells_center_part[3*i_elt+i_dim] = 0.0;
         nelt_faces = connec_cells_idx[i_elt+1] - connec_cells_idx[i_elt];
         for (int i_faces = 0; i_faces < nelt_faces; i_faces++) {
           connec_cells[ ind_faces ] = PDM_ABS(connec_cells[ ind_faces ]);
           _cells_center_part[3*i_elt+i_dim] += faces_center_part[3*connec_cells[ ind_faces ] + i_dim]; 
           ind_faces++; 
         } //i_faces
         _cells_center_part[3*i_elt+i_dim] /= double(nelt_faces);
       } //i_elt
     } //i_dim
          
     PDM_gnum_set_from_coords (_pdmGNum_handle_index, i_part, n_elts, _cells_center_part, NULL);
     
     _isSet[i_part] = true;
     _n_elt[i_part] = n_elts;
     _part_id.push_back(i_part);
     _n_part_def=_n_part_def+1;
     _cells_center[i_part] = _cells_center_part;
     _connec_cells_idx.insert( std::pair < int, int* >         (i_part,connec_cells_idx) );  
     _connec_cells.insert    ( std::pair < int, int* >         (i_part,connec_cells)     );   
     _connec_faces_idx.insert( std::pair < int, int* >         (i_part,connec_faces_idx) );  
     _connec_faces.insert    ( std::pair < int, int* >         (i_part,connec_faces)     );

    
     if( isSet() ) {
       PDM_gnum_compute (_pdmGNum_handle_index);
       for (int i = 0;i<_n_part;i++) {

         _global_num[i] = const_cast<CWP_g_num_t*> (PDM_gnum_get (_pdmGNum_handle_index, i));
         if (not inPDMDB() ) PDM_Mesh_nodal_block_poly3d_set (_pdmNodal_handle_index,
                                      _block_id,
                                      i,    
                                      _n_elt[i],
                                      _n_faces[i],
                                      _connec_faces_idx[i],
                                      _connec_faces[i],
                                      _connec_cells_idx[i], 
                                      _connec_cells[i],   
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

