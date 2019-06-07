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
#include "factory.hpp"
#include "block.hxx"



    int _Mesh_nodal_block_std_type_size_get(CWP_Block_t type) {
   
      
      switch (type) {

       case CWP_BLOCK_NODE: return 1;
       break;
       
       case CWP_BLOCK_EDGE2: return 2;
       break;
                        
       case CWP_BLOCK_FACE_TRIA3: return 3;
       break;

       case CWP_BLOCK_FACE_QUAD4: return 4;
       break;
                       
       case CWP_BLOCK_CELL_TETRA4: return 4;
       break;

       case CWP_BLOCK_CELL_HEXA8: return 8;
       break;

       case CWP_BLOCK_CELL_PYRAM5: return 5;
       break;
 
       case CWP_BLOCK_CELL_PRISM6: return 6;
       break;
       
      }
   }


namespace cwipi {

  /**
   * \typedef FB
   * 
   * \brief Block Factory
   *
   *  A block \ref Factory wich makes \ref Block 
   *  class objects.
   *  The type of block objects build depends on the
   *  block type \ref CWP_Block_t .
   *
   */
   
  typedef Factory<Block, CWP_Block_t> FB;

  /**
    * \brief Mesh constructor
    * 
    * Construct the CWIPI mesh by using paradigm nodal methods.
    * 
    * \param [in] npart Number of mesh partitions.
    *
    */
     


  Mesh::Mesh(const MPI_Comm &localComm,
             Visu* visu,
             const int npart,
             const CWP_Displacement_t displacement) 
             : _localComm(localComm),
               _nDim(-1), _order(-1),_nBlocks(0),
               _pdmGNum_handle_index(-1),
                //_hoOrdering (NULL),
               _pdmNodal(NULL),_visu(visu),
               _displacement(displacement) { 

    PDM_MPI_Comm _pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&localComm));
   // _pdm_localComm=&pdm_localComm;
    // pdm_nodal building
    _pdmNodal_handle_index = PDM_Mesh_nodal_create     (npart,_pdm_localComm);
    _pdmGNum_handle_index  = PDM_gnum_create           (3, npart, PDM_FALSE, 1e-3, _pdm_localComm);
    _npart                 = PDM_Mesh_nodal_n_part_get (_pdmNodal_handle_index); 
    _nVertex   .resize(npart);
    _nElts     .resize(npart);
    _connec_idx.resize(npart,NULL);
    _connec    .resize(npart,NULL);
    _gnum_elt  .resize(npart,NULL);
    _elt_centers  .resize(npart,NULL);
    
  //  _coords .resize(npart,NULL);
  }

  Mesh::~Mesh() 
  {
#if defined(DEBUG) && 0
    cout << "destroying mesh of partition  : TODO" << endl;
#endif


  }




  void Mesh::eltCentersCompute(int i_part){

      int nb_part = getNPart();
      int n_elt_part = getPartNElts(i_part);

      if(_elt_centers[i_part] == NULL)
        _elt_centers[i_part] = (double*)malloc(sizeof(double)*3* n_elt_part);


      int ind_idx=0;
      for(int i=0;i<_nBlocks;i++){
        int id_block = _blocks_id[i];
        CWP_Block_t  block_type = Mesh_nodal_block_type_get  (id_block );
        int n_elt = _blockDB[id_block] -> NEltsGet()[i_part];

          const double* elt_centers_block = _blockDB[id_block] -> eltCentersGet(i_part);

          for(int j=0;j<n_elt;j++){
            for(int k=0;k<3;k++){
              _elt_centers[i_part][3*ind_idx+k]=elt_centers_block[3*j+k];
            }
            ind_idx++;
          }//end loop j
          
      }//end loop on block

  }
  
  void Mesh::connecCompute(int i_part){

        int nb_part = getNPart();
        int n_elt_part = getPartNElts(i_part);

        int connec_size=0;
        for(int i=0;i<_nBlocks;i++){
          int id_block = _blocks_id[i];
          CWP_Block_t  block_type = Mesh_nodal_block_type_get  (id_block );
          
          if(block_type!=CWP_BLOCK_CELL_POLY and block_type!=CWP_BLOCK_FACE_POLY){
            int block_type_size = _Mesh_nodal_block_std_type_size_get(block_type);
            int n_elt = _blockDB[id_block] -> NEltsGet()[i_part];
            connec_size+=n_elt*block_type_size;   
          }
          if(block_type == CWP_BLOCK_FACE_POLY){
            int n_elt = _blockDB[id_block] -> NEltsGet()[i_part];
            int* connec_idx = _blockDB[id_block] -> ConnecIDXGet()[i_part];
            connec_size+=connec_idx[n_elt];   
          }
          
        }//end loop on block
        

      _connec_idx[i_part] = (int*)malloc(sizeof(int)* n_elt_part+1);              
      _connec[i_part] = (int*)malloc(sizeof(int)* connec_size);
      _gnum_elt[i_part] = (CWP_g_num_t*)malloc(sizeof(CWP_g_num_t)* n_elt_part);

      int ind_idx=0;
      _connec_idx[i_part][0]=0;
      for(int i=0;i<_nBlocks;i++){
        int id_block = _blocks_id[i];
        CWP_Block_t  block_type = Mesh_nodal_block_type_get  (id_block );
        int n_elt = _blockDB[id_block] -> NEltsGet()[i_part];
        
        
        if(block_type!=CWP_BLOCK_CELL_POLY and block_type!=CWP_BLOCK_FACE_POLY){
          int block_type_size = _Mesh_nodal_block_std_type_size_get(block_type);
                   
          std::map<int,int*> connect = _blockDB[id_block] -> ConnecGet();
          int* partconnect = connect[i_part];
          CWP_g_num_t* gnum_block = _blockDB[id_block] -> GNumMeshGet(i_part);

          for(int j=0;j<n_elt;j++){
            _connec_idx[i_part][ind_idx+1]=_connec_idx[i_part][ind_idx]+block_type_size;
            _gnum_elt[i_part][ind_idx]=gnum_block[j];
            
            for(int k=0;k<block_type_size;k++){
              _connec[i_part][ _connec_idx[i_part][ind_idx] +k] = partconnect[ j*block_type_size + k ];
            }
            ind_idx++;
          }//end loop j
        }
        
        if( block_type==CWP_BLOCK_FACE_POLY){

          std::map<int,int*> connectIdx = _blockDB[id_block] -> ConnecIDXGet();        
          std::map<int,int*> connect    = _blockDB[id_block] -> ConnecGet();
          int* partconnect    = connect   [i_part];
          int* partconnectIdx = connectIdx[i_part];
          CWP_g_num_t* gnum_block = _blockDB[id_block] -> GNumMeshGet(i_part);

          for(int j=0;j<n_elt;j++){
            _connec_idx[i_part][ind_idx+1]=_connec_idx[i_part][ind_idx]+ (partconnectIdx[j+1]-partconnectIdx[j]);
            _gnum_elt[i_part][ind_idx]=gnum_block[j];
            
            for(int k = _connec_idx[i_part][ind_idx] ;k<_connec_idx[i_part][ind_idx+1]; k++){
              _connec[i_part][k] = partconnect[ partconnectIdx[j] + k - _connec_idx[i_part][ind_idx] ];
            }
            ind_idx++;
          }//end loop j
        }        
         
      }//end loop on block
     
  }


 CWP_g_num_t* Mesh::GNumEltsGet(int i_part){
    if(_connec_idx[i_part]==NULL) {
      connecCompute(i_part);
    }//end if NULL
  
    return _gnum_elt[i_part];    
  }


 double* Mesh::eltCentersGet(int i_part){
   if(_elt_centers[i_part]==NULL || _displacement != CWP_DISPLACEMENT_STATIC) {
     eltCentersCompute(i_part);
   }//end if NULL
  
   return _elt_centers[i_part];
 }




 int Mesh::GNVerticeGet(int i_part){
     int n_vert = getPartNVertex(i_part);
     
     //PDM_MPI_Comm pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm(_localComm);
 
     int n_vert_global=0;

     PDM_MPI_Allreduce(&n_vert,&n_vert_global,1,PDM_MPI_INT, PDM_MPI_SUM,_pdm_localComm);
     
     return n_vert_global;
 }
 
    
   // PDM_MPI_Allreduce(elts, som_elts, 5, PDM_MPI_INT, PDM_MPI_SUM, mesh->pdm_mpi_comm);
    
    
 int Mesh::GNEltGet(int i_part){
     int n_elt = getPartNElts(i_part);
     
     //PDM_MPI_Comm pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm((void*)_localComm));
 
     int n_elt_global=0;

     PDM_MPI_Allreduce(&n_elt,  &n_elt_global, 1,  PDM_MPI_INT, PDM_MPI_SUM,_pdm_localComm);
     
     return n_elt_global;
 }




  int* Mesh::connecIdxGet(int i_part) {
    if(_connec_idx[i_part]==NULL) {
      connecCompute(i_part);
    }//end if NULL
  
    return _connec_idx[i_part];    
  }
  
  int* Mesh::connecGet(int i_part) {
    if(_connec[i_part]==NULL) {
      connecCompute(i_part);
    }//end if NULL
  
    return _connec[i_part];    
  }
  
  
  void Mesh::updateBlockDB()
  {
     _nBlocks     = PDM_Mesh_nodal_n_blocks_get (_pdmNodal_handle_index);
     _blocks_id   = PDM_Mesh_nodal_blocks_id_get(_pdmNodal_handle_index);
     std::map<int,Block*>::iterator it;
     
     for(int i=0;i<_nBlocks;i++){
        int id_block = _blocks_id[i];
        it = _blockDB.find(id_block);
        if(it == _blockDB.end()) {
           CWP_Block_t  block_type = Mesh_nodal_block_type_get  (id_block ); 
           Block *newBlock = FB::getInstance().CreateObject(block_type);
           newBlock -> SetinPDMDB();
           newBlock -> FromPDMBlock(id_block,this);
           _blockDB.insert    ( std::pair < int, Block* >    (id_block,newBlock) );  
        }//endif
     }//endfor
     for(int i=0;i<_npart;i++){
       _nElts[i]  = PDM_Mesh_nodal_n_cell_get(_pdmNodal_handle_index,i);
     }
     
   }


   void Mesh::nodal_coord_set(const int   i_part,
                              const int   n_vtx,
                              double      coords[],
                              CWP_g_num_t global_num[])
   {  
   
      _coords .push_back (coords);
      _nVertex[i_part]  = n_vtx;
      PDM_gnum_set_from_coords (_pdmGNum_handle_index, i_part, n_vtx, coords, NULL);
      
      if(coordsDefined() and global_num==NULL)
        {
         PDM_gnum_compute (_pdmGNum_handle_index);

         CWP_g_num_t* global_num_part = NULL;
         
         for(int i_part2=0;i_part2<_npart;i_part2++) {
           global_num_part =const_cast<CWP_g_num_t*>(PDM_gnum_get (_pdmGNum_handle_index, i_part2));

           PDM_Mesh_nodal_coord_set(_pdmNodal_handle_index,
                                    i_part2,    
                                    _nVertex[i_part2],    
                                    _coords[i_part2],   
                                    global_num_part);
         
           _global_num.insert( std::pair < int, CWP_g_num_t* > (i_part2,global_num_part) );

           if(_visu -> isCreated() && _displacement == CWP_DISPLACEMENT_STATIC) {
             _visu -> GeomCoordSet(i_part2,
                                   _nVertex[i_part2],
                                   _coords[i_part2],
                                   global_num_part);  
           }
         
         }//loop i_part2
       }//endif coordsDefined() and global_num==NULL
   }
  






  
   CWP_Block_t Mesh::Mesh_nodal_block_type_get(const int id_block) {
   
      PDM_Mesh_nodal_elt_t pdm_id_block = PDM_Mesh_nodal_block_type_get(_pdmNodal_handle_index,id_block);
      switch (pdm_id_block) {

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
  
  
  
  void Mesh::stdBlockSet( const int              i_part,
                          const int              block_id,
                          const int              n_elts,
                          int                    connec[],  
                          CWP_g_num_t            global_num[]
                         ) {

     _blockDB [block_id] -> blockSet(i_part,n_elts,connec,global_num);

     _nElts[i_part]  = PDM_Mesh_nodal_n_cell_get(_pdmNodal_handle_index,
                                                   i_part);
                                                   
     _blocks_id = PDM_Mesh_nodal_blocks_id_get(_pdmNodal_handle_index);
     _nBlocks   = PDM_Mesh_nodal_n_blocks_get (_pdmNodal_handle_index);
     
     if(_visu -> isCreated() && _displacement == CWP_DISPLACEMENT_STATIC) {

        _visu -> GeomBlockStdSet (_id_visu[block_id],
                                  i_part,
                                  n_elts,
                                  connec,
                                  global_num);
      }  
 }


   void Mesh::poly2DBlockSet( const int              i_part,
                              const int              block_id,
                              const int              n_elts,
                              int                    connec_idx[],
                              int                    connec[], 
                              CWP_g_num_t            global_num[]
                            )
   {
     if(_coords[i_part]==NULL) bftc_error(__FILE__, __LINE__, 0, 
            "Set the partition coordinates vertices before finalizing.\n");
     
     _blockDB [block_id] -> blockSet(i_part,
                                     n_elts,
                                     connec_idx,
                                     connec,
                                     global_num);
     
     _nElts[i_part]  = PDM_Mesh_nodal_n_cell_get(_pdmNodal_handle_index,
                                                   i_part);
                                                   
     _blocks_id = PDM_Mesh_nodal_blocks_id_get(_pdmNodal_handle_index);
     _nBlocks   = PDM_Mesh_nodal_n_blocks_get (_pdmNodal_handle_index);     
     

     if(_visu -> isCreated() && _displacement == CWP_DISPLACEMENT_STATIC) {
     
        _visu -> GeomBlockPoly2D (_id_visu[block_id],
                                  i_part,
                                  n_elts,
                                  connec_idx,
                                  connec, 
                                  global_num);                                                          
     }
   }



/**********************************************************************/
  void Mesh::geomFinalize() {
    int g_num_computation_required = 0;
    std::map<int,cwipi::Block*>::iterator it = _blockDB.begin();
    while(it != _blockDB.end()) {
      for(int i_part =0;i_part<_npart;i_part++) {    
         int block_id = it -> second -> blockIDGet();
         CWP_g_num_t* global_num = it -> second -> GNumMeshGet(i_part);
         if(global_num == NULL) g_num_computation_required = 1;    
         if(g_num_computation_required == 1) break;
       }
     if(g_num_computation_required == 1) break;
     it++;
    }

    if(g_num_computation_required == 1) {
      PDM_Mesh_nodal_g_num_in_mesh_compute(_pdmNodal_handle_index);
      std::map<int,cwipi::Block*>::iterator it = _blockDB.begin();
      
      
        while(it != _blockDB.end()) {
          for(int i_part =0;i_part<_npart;i_part++) {
            int block_id = it -> second -> blockIDGet();
            CWP_g_num_t* global_num = it -> second -> GNumMeshGet(i_part);
            if(_visu -> isCreated() && _displacement == CWP_DISPLACEMENT_STATIC) _visu -> GeomBlockGNumMeshSet (_id_visu[block_id],
                                           i_part,
                                           global_num);
          } //Loop on i_part
        it++;
        }//Loop on blockDB
      
    } // end if g_num_computation_required
         
    if(_visu -> isCreated() && _displacement == CWP_DISPLACEMENT_STATIC ) _visu -> GeomWrite();   
                            
  } 
           

/**********************************************************************/


   void Mesh::poly3DBlockSet( const int              i_part,
                              const int              block_id,
                              const int              n_elts,
                              const int              n_faces,
                              int                    connec_faces_idx[],
                              int                    connec_faces[],
                              int                    connec_cells_idx[],
                              int                    connec_cells[], 
                              CWP_g_num_t            global_num[]
                            ) {
     if(_coords[i_part]==NULL) bftc_error(__FILE__, __LINE__, 0, 
            "Set the partition coordinates vertices before finalizing.\n");
     
     _blockDB [block_id] -> blockSet(i_part,n_elts,n_faces,
                                     connec_faces_idx,connec_faces,
                                     connec_cells_idx,connec_cells,
                                     NULL);   
                                     
     if(_visu -> isCreated() && _displacement == CWP_DISPLACEMENT_STATIC) {

        _visu -> GeomBlockPoly3D   (_id_visu[block_id],
                                    i_part,
                                    n_elts,
                                    n_faces,
                                    connec_faces_idx,
                                    connec_faces,
                                    connec_cells_idx,
                                    connec_cells, 
                                    _global_num[i_part]);                        
      }
                                    
                                     
     for(int i=0;i<_npart;i++){
       _nElts[i]  = PDM_Mesh_nodal_n_cell_get(_pdmNodal_handle_index,
                                                   i_part);
     }                          
   }
   
  
   void Mesh::meshDel()
   {
     PDM_Mesh_nodal_free(_pdmNodal_handle_index);
   }
   
  
   int Mesh::blockAdd(const CWP_Block_t      block_type
                     )
   {      
      Block *myBlock = FB::getInstance().CreateObject(block_type);
      myBlock -> BlockAdd(block_type,this);
      int block_id = myBlock -> blockIDGet();
      _blockDB.insert ( std::pair <int,Block*> (block_id,myBlock) );                        
      
     _blocks_id = PDM_Mesh_nodal_blocks_id_get(_pdmNodal_handle_index);
     _nBlocks   = PDM_Mesh_nodal_n_blocks_get (_pdmNodal_handle_index); 
      
      if(_visu -> isCreated() && _displacement == CWP_DISPLACEMENT_STATIC) {
        int id_visu = _visu -> GeomBlockAdd(block_type);
        _id_visu.insert(std::pair <int,int> (block_id,id_visu));
      }

      return block_id;   
                  
   }
  
     void Mesh::fromCellFaceSet(const int   i_part,
                        const int   n_cells,
                        int         cell_face_idx[],
                        int         cell_face[],
                        int         n_faces,
                        int         face_vtx_idx[],
                        int         face_vtx[],
                        CWP_g_num_t global_num[]) {
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
       updateBlockDB();                                                          
     }
   

     void Mesh::fromFacesEdgeSet(const int   i_part,
                                 const int   n_faces,
                                 int         face_edge_idx[],
                                 int         face_edge[],
                                 const int   n_edges,
                                 int         edge_vtx_idx[],
                                 int         edge_vtx[],
                                 CWP_g_num_t global_num[]) {
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
                                          NULL);              
       updateBlockDB();

     }

}
