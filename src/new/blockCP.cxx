

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


  void BlockCP::blockSet(int i_part,int n_elts,int n_faces,
                        int* connec_faces_idx, 
                        int* connec_faces,
                        int* connec_cells_idx,
                        int* connec_cells,
                        CWP_g_num_t* global_num) {

     CWP_UNUSED(global_num );
     double* _cells_center_part = (double*)malloc (sizeof(double) * 3 * n_elts);
     double* faces_center_part = (double*)malloc (sizeof(double) * 3 * n_faces);
     double* coords = static_cast<Mesh*>(_mesh)->getCoordinates(i_part);

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
          
     
     _isSet[i_part] = true;
     _n_elt[i_part] = n_elts;
     _part_id.push_back(i_part);
     _n_part_def=_n_part_def+1;
     _cells_center[i_part] = _cells_center_part;
     _connec_cells_idx.insert( std::pair < int, int* >         (i_part,connec_cells_idx) );  
     _connec_cells.insert    ( std::pair < int, int* >         (i_part,connec_cells)     );   
     _connec_faces_idx.insert( std::pair < int, int* >         (i_part,connec_faces_idx) );  
     _connec_faces.insert    ( std::pair < int, int* >         (i_part,connec_faces)     );

  
  }



    void BlockCP::geomFinalize(){

      for(int i_part = 0; i_part<_n_part; i_part++){

        Visu* visu = ((Mesh*)_mesh)->getVisu();
        if(visu->isCreated() && ((Mesh*)_mesh)->getDisplacement() == CWP_DYNAMIC_MESH_STATIC) {
          visu->GeomBlockPoly3D( ((Mesh*)_mesh)->getIdVisu( _block_id_cwipi ),
                                    i_part,
                                    _n_elt[i_part] ,
                                    _n_faces[i_part],
                                    _connec_faces_idx [i_part],
                                    _connec_faces     [i_part],
                                    _connec_cells_idx[i_part],
                                    _connec_cells    [i_part],
                                    _global_num [i_part]);


        }
      } //end i_part


    }
}

