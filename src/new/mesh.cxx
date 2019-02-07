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
#include <pdm_gnum.h>
#include <bftc_error.h>
#include <bftc_printf.h>
#include "cwp.h"

namespace cwipi {



  /**
    * \brief Mesh constructor
    * 
    * Construct the CWIPI mesh by using paradigm nodal methods.
    * 
    * \param [in] npart Number of mesh partitions.
    *
    */
     
  Mesh::Mesh(const MPI_Comm &localComm,
         const int npart)
    : _localComm(localComm),
      _nDim(-1), _order(-1),_nBlocks(0),_nVertex(NULL), _coords(NULL),
      //_hoOrdering (NULL),
      _isNodalFinalized(false),
      _pdmNodal(NULL)
  { PDM_MPI_Comm pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&localComm));
    // pdm_nodal building
    _pdmNodal_handle_index = PDM_Mesh_nodal_create     (npart,pdm_localComm);
    _pdmGNum_handle_index  = PDM_gnum_create           (3, npart, PDM_FALSE, 1e-3, pdm_localComm);
    _npart                 = PDM_Mesh_nodal_n_part_get (_pdmNodal_handle_index); 
    _nVertex.resize(npart);
    _nElts.resize(npart);
    _coords .resize(npart);
    
  }


  Mesh::~Mesh() 
  {
#if defined(DEBUG) && 0
    cout << "destroying mesh of partition  : TODO" << endl;
#endif

  }
  
  
  void Mesh::updateBlockDB(int id_part)
  {
     _nBlocks   = PDM_Mesh_nodal_n_blocks_get (_pdmNodal_handle_index);
     _blocks_id = PDM_Mesh_nodal_blocks_id_get(_pdmNodal_handle_index);

     for(int i=0;i<_nBlocks;i++)
       { int id_block = _blocks_id[i];
       
         std::map<int,_block>::iterator It = _blocks.find(id_block);
         if (It != _blocks.end())
            if(It->second._id_part != id_part) continue;
         bftc_printf("WP_Block_t   block_type = Mesh_nodal_bl %i\n",i);  
         CWP_Block_t   block_type = Mesh_nodal_block_type_get     (id_block ); 
         int           nElts      = PDM_Mesh_nodal_block_n_elt_get(_pdmNodal_handle_index,
                                                                   id_block,
                                                                   id_part );

         bftc_printf("int*  connec_idx        = NULL; %i\n",i);  
         int*  connec_idx        = NULL;
         int*  connec            = NULL;
         int    nFaces           = 0   ;
         int*   n_faces          = NULL;
         int*   connec_faces_idx = NULL;  
         int*   connec_faces     = NULL;
         bool   isBlockFinalized = false;
         CWP_g_num_t*  global_num_block = NULL;
         int*  parent_num = NULL;

         switch(block_type) {
             case PDM_MESH_NODAL_POLY_2D: 
                       PDM_Mesh_nodal_block_poly2d_get(_pdmNodal_handle_index,
                                                                        id_block,
                                                                        id_part,
                                                                        &connec_idx,
                                                                        &connec); 
             break;
             case PDM_MESH_NODAL_POLY_3D: 
                       PDM_Mesh_nodal_block_poly3d_get(_pdmNodal_handle_index,
                                                                        id_block,
                                                                        id_part,
                                                                        &nFaces,
                                                                        &connec_faces_idx,
                                                                        &connec_faces,
                                                                        &connec_idx,
                                                                        &connec); 
             break;    
             default: PDM_Mesh_nodal_block_std_get (_pdmNodal_handle_index,
                                                    id_block,
                                                    id_part,
                                                    &connec); 
             }
             
         if(connec!=NULL) isBlockFinalized = true;
         bftc_printf("if (isBlockFinalized) %i %i\n",i,isBlockFinalized);  
         if (isBlockFinalized)
         { 
           bftc_printf("CWP_g_num_t*  parent_num\n");                                                          
           parent_num
                    = PDM_Mesh_nodal_block_parent_num_get (_pdmNodal_handle_index,
                                                             id_block,
                                                             id_part );
                                                           
          if(parent_num!=NULL)
             {bftc_printf("PDM_Mesh_nodal_g_num_in_block_compute %i block_type %i\n",i,block_type);  
              PDM_Mesh_nodal_g_num_in_block_compute(_pdmNodal_handle_index,
                                                    id_block);
                                                     
              bftc_printf("CWP_g_num_t*  global_num_block\n");       
              global_num_block
                      = PDM_Mesh_nodal_block_g_num_get (_pdmNodal_handle_index,
                                                        id_block,
                                                        id_part );
                                                     
             }
           else
             {bftc_printf("parent_num==NULL %i %i\n",i,_nBlocks);
             }                                                                           
         }         
           bftc_printf("block myBloc %i %i\n",i,_nBlocks); 
         _block myBlock = {id_part,
                           id_block,
                           block_type,
                           nElts,
                           connec_idx,
                           connec,
                           nFaces,
                           n_faces,
                           connec_faces_idx,
                           connec_faces,
                           global_num_block,
                           parent_num,
                           isBlockFinalized};
             
         It = _blocks.find(id_block);
         if (It != _blocks.end())
             It->second = myBlock;
         else
             _blocks.insert( std::pair < int, _block > (id_block,myBlock) );  
     }
           bftc_printf("_nElts[i]=PDM_Mesh_nodal_n_\n"); 
     for(int i=0;i<_npart;i++)
       {
         _nElts[i]=PDM_Mesh_nodal_n_cell_get(_pdmNodal_handle_index,i);
       }
    
   }
  
  

   void Mesh::nodal_coord_set(const int   i_part,
                              const int   n_vtx,
                              double      coords[],
                              CWP_g_num_t global_num[])
   {
      if(global_num==NULL)
        {
         PDM_gnum_set_from_coords (_pdmGNum_handle_index, i_part, n_vtx, coords, NULL);
         PDM_gnum_compute (_pdmGNum_handle_index);
         global_num =const_cast<CWP_g_num_t*>(PDM_gnum_get (_pdmGNum_handle_index, i_part));
         _global_num.insert( std::pair < int, CWP_g_num_t* > (i_part,global_num) );
        }
        
      PDM_Mesh_nodal_coord_set(_pdmNodal_handle_index,
                               i_part,    
                               n_vtx,    
                               coords,   
                               global_num);
                               
      _coords[i_part]   = const_cast<double*>(  PDM_Mesh_nodal_vertices_get(_pdmNodal_handle_index,
                                                        i_part) );
      _nVertex[i_part]  = PDM_Mesh_nodal_n_vertices_get(_pdmNodal_handle_index,
                                                        i_part);
   }
  
  
   CWP_Block_t Mesh::Mesh_nodal_block_type_get(const int id_block
                                        )
   {
     PDM_Mesh_nodal_elt_t PDM_block_type = PDM_Mesh_nodal_block_type_get(_pdmNodal_handle_index,
                                                                   id_block);
                                                                   
     switch (PDM_block_type) {

       case PDM_MESH_NODAL_POINT: return CWP_BLOCK_NODE;
       break;
       
       case PDM_MESH_NODAL_BAR2: return CWP_BLOCK_EDGE2;
       break;
                        
       case PDM_MESH_NODAL_TRIA3: return CWP_BLOCK_FACE_TRIA3;
       break;

       case PDM_MESH_NODAL_QUAD4: return CWP_BLOCK_FACE_QUAD4;
       break;
                       
       case PDM_MESH_NODAL_TETRA4: return CWP_BLOCK_CELL_TETRA4;
       break;

       case PDM_MESH_NODAL_POLY_2D: return CWP_BLOCK_FACE_POLY;
       break;
       
       case PDM_MESH_NODAL_HEXA8: return CWP_BLOCK_CELL_HEXA8;
       break;

       case PDM_MESH_NODAL_PYRAMID5: return CWP_BLOCK_CELL_PYRAM5;
       break;
       
       case PDM_MESH_NODAL_PRISM6: return CWP_BLOCK_CELL_PRISM6;
       break;
       
       case PDM_MESH_NODAL_POLY_3D: return CWP_BLOCK_CELL_POLY;
       break;
       
      }
  
   }
  
  
   void Mesh::endSet()
   { if(_coords.size()==0) bftc_error(__FILE__, __LINE__, 0, 
              "Set coordinates vertices before finalizing.\n");
              
     for(int i_part=0;i_part<_npart;i_part++)
        if(_coords[i_part]==NULL) bftc_error(__FILE__, __LINE__, 0, 
              "Set all the partitions coordinates vertices before finalizing.\n");

     
     std::map< int,_block >::iterator block_it = _blocks.begin();
     while (block_it != _blocks.end())
        {
           _block myblock = block_it->second;
          
           int id_part  = myblock._id_part;
           int id_block = myblock._id_block;
           
           PDM_Mesh_nodal_elt_t pdm_block_type
              =PDM_Mesh_nodal_block_type_get(_pdmNodal_handle_index,
                                             id_block);
                                          
           int* connec_idx       = myblock._connec_idx;
           int* connec           = myblock._connec;
           int* connec_faces     = myblock._connec_faces;
           int* connec_faces_idx = myblock._connec_faces_idx;
           CWP_g_num_t* global_num = _global_num[id_part];
           int* parent_num = myblock._parent_num;
           CWP_g_num_t* global_num_block = myblock._global_num_block;
           
           const CWP_g_num_t *g_num;
           
           switch (pdm_block_type) {
                                           
           case PDM_MESH_NODAL_POLY_2D :
              bftc_printf("block_poly2d_set\n");
              PDM_Mesh_nodal_block_poly2d_set (_pdmNodal_handle_index,
                                               id_block,
                                               id_part,    
                                               myblock._nElts,
                                               connec_idx, 
                                               connec,   
                                               global_num,
                                               parent_num);
              bftc_printf("After block_poly2d_set\n");
              break;
                                               
           case PDM_MESH_NODAL_POLY_3D :
              bftc_printf("block_poly3d_set\n");
              PDM_Mesh_nodal_block_poly3d_set (_pdmNodal_handle_index,
                                               id_block,
                                               id_part,    
                                               myblock._nElts,
                                               myblock._nFaces,
                                               connec_faces_idx, 
                                               connec_faces,   
                                               connec_idx, 
                                               connec,   
                                               global_num,
                                               parent_num);
              break;
                                               
           default :
              PDM_Mesh_nodal_block_std_set(_pdmNodal_handle_index,
                                           id_block,
                                           id_part,    
                                           myblock._nElts,    
                                           connec,   
                                           global_num,
                                           parent_num);
              break;
           } 
           bftc_printf("updateBlockDB(id_part);\n");
           updateBlockDB(id_part);
           block_it++;
        }
        
     
   }
   
   
   void Mesh::meshDel()
   {
     PDM_Mesh_nodal_free(_pdmNodal_handle_index);
   }
   
  
   int Mesh::blockAdd(const int              i_part,
                       const CWP_Block_t      block_type,
                       const int              n_elts,
                       int                    connec_idx[],
                       int                    connec[],
                       const int              n_faces,
                       int                    face_vtx_idx[],
                       int                    face_vtx[],    
                       CWP_g_num_t            global_num[],
                       int                    parent_num[])
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
      
      _block myBlock = {i_part,
                        id_block,
                        block_type,
                        n_elts,
                        connec_idx,
                        connec,
                        n_faces,
                        NULL,
                        face_vtx_idx,
                        face_vtx,
                        global_num,
                        parent_num,
                        false};   
      _blocks.insert( std::pair < int, _block > (id_block,myBlock));                                                         
      return id_block;   
                  
   }
   

     void Mesh::fromCellFaceSet(const int   i_part,
                        const int   n_cells,
                        int         cell_face_idx[],
                        int         cell_face[],
                        int         n_faces,
                        int         face_vtx_idx[],
                        int         face_vtx[],
                        CWP_g_num_t global_num[])
     {
       int face_vtx_nb[n_faces];
       for (int i=0;i<n_faces;i++)
         {
           face_vtx_nb[i] = face_vtx_idx[i+1] - face_vtx_idx[i] ;
         }
      
       int cell_face_nb[n_cells];
       for (int i=0;i<n_cells;i++)
         {
           cell_face_nb[i] = cell_face_idx[i+1] - cell_face_idx[i] ;
         }
      
       PDM_Mesh_nodal_cell3d_cellface_add(_pdmNodal_handle_index,
                                          i_part,
                                          n_cells,
                                          n_faces,
                                          face_vtx_idx,
                                          face_vtx_nb, 
                                          face_vtx,
                                          cell_face_idx,
                                          cell_face_nb, 
                                          cell_face,
                                          _global_num[i_part]);    
                                          
        updateBlockDB(i_part);                                       
                               
     }
   

     void Mesh::fromFacesEdgeSet(const int   i_part,
                                 const int   n_faces,
                                 int         face_edge_idx[],
                                 int         face_edge[],
                                 const int   n_edges,
                                 int         edge_vtx_idx[],
                                 int         edge_vtx[],
                                 CWP_g_num_t parent_num[])
     {
       int edge_vtx_nb[n_edges];
       for (int i=0;i<n_edges;i++)
         {
           edge_vtx_nb[i] = edge_vtx_idx[i+1] - edge_vtx_idx[i] ;
         }
      
       int face_edge_nb[n_faces];
       for (int i=0;i<n_faces;i++)
         {
           face_edge_nb[i] = face_edge_idx[i+1] - face_edge_idx[i] ;
         }
       
       PDM_Mesh_nodal_cell2d_celledge_add(_pdmNodal_handle_index,
                                          i_part,
                                          n_faces,
                                          n_edges,
                                          edge_vtx_idx,
                                          edge_vtx_nb, //Number of vertices for each edge
                                          edge_vtx,
                                          face_edge_idx,
                                          face_edge_nb, //Number of edges for each faces
                                          face_edge,
                                          parent_num);              
       updateBlockDB(i_part);

     }


}
