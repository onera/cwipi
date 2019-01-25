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
#include <mesh.hxx>

#include <mpi.h>

#include <pdm_mesh_nodal.h>
#include <bftc_error.h>
#include <bftc_printf.h>
#include "cwp.h"





namespace cwipi {



  Mesh::Mesh(const MPI_Comm &localComm,
         const int npart)
    : _localComm(localComm),
      _nDim(-1), _order(-1), _nVertex(NULL), _nElts(NULL), _coords(NULL),
      //_hoOrdering (NULL),
      _isNodalFinalized(false),
      _pdmNodal(NULL)
  {
    // pdm_nodal building
    _pdmNodal_handle_index = PDM_Mesh_nodal_create(npart,localComm);
    _nVertex.resize(npart);
    _coords.resize(npart);
    _nElts.resize(npart);
    _npart=PDM_Mesh_nodal_n_part_get(_pdmNodal_handle_index);
  }


  Mesh::~Mesh() 
  {
#if defined(DEBUG) && 0
    cout << "destroying mesh of partition  : TODO" << endl;
#endif

  }
  
  
  
   void Mesh::nodal_coord_set(const int   i_part,
                        const int   n_vtx,
                        double      coords[],
                        CWP_g_num_t global_num[])
   {
      PDM_Mesh_nodal_coord_set(_pdmNodal_handle_index,
                               i_part,    
                               n_vtx,    
                               coords,   
                               global_num);
                               
      _coords[i_part]   = const_cast<double*>(  PDM_Mesh_nodal_vertices_get(_pdmNodal_handle_index,
                                                        i_part) );
      _nVertex[i_part]  = PDM_Mesh_nodal_n_vertices_get(_pdmNodal_handle_index,
                                                        i_part);
      
      _global_num.insert( std::pair < int, CWP_g_num_t* > (i_part,global_num) );

      if(global_num==NULL) bftc_printf("Partition %i global numbering must be computed.",i_part);

   }
  
  
   void Mesh::endSet()
   { 
     if(_coords.size()==0) bftc_error(__FILE__, __LINE__, 0, "Set coordinates vertices before finalizing.\n");
     for(int i_part=1;i_part<=_npart;i_part++)
       {if(_coords[i_part]==NULL) bftc_error(__FILE__, __LINE__, 0, "Set all the paritions coordinates vertices before finalizing.\n");
       }
     
     std::map< id_part_block, int* >::iterator connec_it = _connect.begin();
     while (connec_it != _connect.end())
        {
           id_part_block id = connec_it->first;
          
           int id_part=id.get_id_part();
           int id_block=id.get_id_block();
           
           const CWP_g_num_t *g_num;
           
           if(_global_num_block[id]==NULL)
             {  
               PDM_Mesh_nodal_g_num_in_block_compute(_pdmNodal_handle_index,
                                                 id_block);  
           
               g_num = PDM_Mesh_nodal_block_g_num_get(_pdmNodal_handle_index,
                                                                     id_block,
                                                                     id_part);
                                                   
             }
             
           PDM_Mesh_nodal_block_std_set(_pdmNodal_handle_index,
                                        id_block,
                                        id_part,    
                                        _nElts[id_part],    
                                        connec_it->second,   
                                        const_cast <CWP_g_num_t *> (g_num),
                                        NULL);
           
           
           if(_global_num.find(id_part)==_global_num.end()) 
              _global_num.insert( std::pair < int, CWP_g_num_t* > (id_part,NULL) ); 
           
           connec_it++;
        }
           
        


     std::map< int, CWP_g_num_t* >::iterator gnum_it = _global_num.begin();
     while (gnum_it != _global_num.end())
        {  int id_part=gnum_it->first; 
           const CWP_g_num_t *g_num = PDM_Mesh_nodal_vertices_g_num_get(_pdmNodal_handle_index,
                                                                        id_part);               
           _global_num[id_part]=const_cast <CWP_g_num_t *>(g_num);
        }

     _isNodalFinalized=true;
   }
   
   
   void Mesh::meshDel()
   {
     PDM_Mesh_nodal_free(_pdmNodal_handle_index);
   }
   
  
   void Mesh::blockAdd(const int              i_part,
                       const Block_Addition_t add_type,
                       const CWP_Block_t      block_type,
                       const int              n_elts,
                       int                    connec_idx[],
                       int                    connec[],
                       CWP_g_num_t            global_num[],
                       CWP_g_num_t            parent_num[]
                       )
   {
      int id_block=-1;
      switch (block_type) {

      case CWP_BLOCK_NODE :
        id_block=PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                            PDM_TRUE,
                                            PDM_MESH_NODAL_POINT);
        break;
      case CWP_BLOCK_EDGE2 :
        id_block=PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                            PDM_TRUE,
                                            PDM_MESH_NODAL_BAR2);                          
        break;
      case CWP_BLOCK_FACE_TRIA3 :
        id_block=PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                            PDM_TRUE,
                                            PDM_MESH_NODAL_TRIA3);                          
        break;
      case CWP_BLOCK_FACE_QUAD4 :
        id_block=PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                            PDM_TRUE,
                                            PDM_MESH_NODAL_QUAD4);                          
        break;
      case CWP_BLOCK_CELL_TETRA4 :
        id_block=PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                            PDM_TRUE,
                                            PDM_MESH_NODAL_TETRA4);                          
        break;
      case CWP_BLOCK_FACE_POLY :
        id_block=PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                            PDM_TRUE,
                                            PDM_MESH_NODAL_POLY_2D);                          
        break;  
      case CWP_BLOCK_CELL_HEXA8 :
        id_block=PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                            PDM_TRUE,
                                            PDM_MESH_NODAL_HEXA8);                          
        break;
      case CWP_BLOCK_CELL_PYRAM5 :
        id_block=PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                            PDM_TRUE,
                                            PDM_MESH_NODAL_PYRAMID5);                          
        break;  
      case CWP_BLOCK_CELL_PRISM6 :
        id_block=PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                            PDM_TRUE,
                                            PDM_MESH_NODAL_PRISM6);                          
        break;
      case CWP_BLOCK_CELL_POLY :
        id_block=PDM_Mesh_nodal_block_add(_pdmNodal_handle_index,
                                            PDM_TRUE,
                                            PDM_MESH_NODAL_POLY_3D);                          
        break;      
      default :
        bftc_error(__FILE__, __LINE__, 0, "%s is not a valid CWIPI Block Type\n"
                   , block_type);
        break;
      }                                        
          
      id_part_block id(i_part,id_block);
      _connect.insert          ( std::pair < id_part_block, int* > (id,connec) );
      _connect_idx.insert      ( std::pair < id_part_block, int* > (id,connec_idx) );
      _add_type.insert         ( std::pair < id_part_block, Block_Addition_t > (id,add_type) );
      _global_num_block.insert ( std::pair < id_part_block, CWP_g_num_t* > (id,global_num) );
      _parent_num_block.insert ( std::pair < id_part_block, CWP_g_num_t* > (id,parent_num) );
      _nElts[i_part]=PDM_Mesh_nodal_n_cell_get(_pdmNodal_handle_index,i_part);     
   
                   
   }




}
