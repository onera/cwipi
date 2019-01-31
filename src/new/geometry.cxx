/*
  This file is part of the CWIPI library. 

  Copyright (C) 2012  ONERA

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

#include <vector>
#include <map>
#include <geometry.hxx>
#include <mpi.h>

#include <pdm_mesh_nodal.h>
#include <pdm_gnum.h>
#include <bftc_error.h>
#include <bftc_printf.h>
#include "cwp.h"



namespace cwipi {

  Geometry::Geometry()
  :
  _tmpVertexField (NULL),
  _tmpDistantField(NULL)
  {
  }
  
  Geometry::~Geometry()
  {
  }
  
  void Geometry::meshInterfSet (const char            *cpl_id,
                                const int              i_part,
                                const int              n_pts,
                                double           coords[],
                                CWP_g_num_t      parent_num[])
  {
     _meshInterf->nodal_coord_set(i_part,
                                  n_pts,
                                  coords,
                                  parent_num);
  }
  
  
  void Geometry::meshInterfBlockAdd (const int             i_part,
                                     const CWP_Block_t     block_type,
                                     const int             n_elts,
                                     int             connec[],
                                     CWP_g_num_t     parent_num[])
  {
    int id_block = _meshInterf->blockAdd(i_part,
                                         block_type,
                                         n_elts,
                                         NULL,
                                         connec,
                                         -1,
                                         NULL,
                                         NULL,    
                                         NULL,
                                         parent_num);
  }
  



  void Geometry::meshInterfFPolyBlockAdd( const int          i_part,
                                          const CWP_Block_t  block_type,
                                          const int          n_elts,
                                          int                connec_idx[],
                                          int                connec[],
                                          CWP_g_num_t        parent_num[])
  {
    int id_block = _meshInterf-> blockAdd(i_part,
                                          CWP_BLOCK_FACE_POLY,
                                           n_elts,
                                           connec_idx,
                                           connec,
                                           -1,
                                           NULL,
                                           NULL,                            
                                           NULL,
                                           parent_num);
  }
  

  void Geometry::meshInterfCPolyBlockAdd( const int     i_part,
                                          const int     n_elts,
                                          int           cell_face_idx[],
                                          int           cell_face[],
                                          const int     n_faces,
                                          int           face_vtx_idx[],
                                          int           face_vtx[],
                                          CWP_g_num_t   parent_num[])
  {
    int id_block = _meshInterf-> blockAdd(i_part,
                                          CWP_BLOCK_FACE_POLY,
                                           n_elts,
                                           cell_face_idx,
                                           cell_face,
                                           n_faces,
                                           face_vtx_idx,
                                           face_vtx,                        
                                           NULL,
                                           parent_num);
  }
  
  
  
  void Geometry::fvmcNodalShared( const int      i_part,
                                  fvmc_nodal_t  *fvmc_nodal)
  {
  
  
  
  
  }
  
  
  

}
