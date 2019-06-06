#ifndef __MESH_H__
#define __MESH_H__
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

#include <mpi.h>

#include <pdm_mesh_nodal.h>
#include <bftc_error.h>

#include "block.hxx"
#include "cwp.h"
#include "visualization.hxx"




namespace cwipi {


  /** 
   * \class Mesh mesh.hxx "mesh.hxx"
   * \brief Interface mesh
   *
   *  This class defines the interface mesh objects.
   * 
   */
  class Coupling;
  class Visu;
  class Mesh {
    
  public:



    /**
     * \brief Mesh constructor
     * 
     * Construct the CWIPI mesh by using paradigm nodal methods.
     * 
     * \param [in] localComm Coupling Communicator.
     * \param [in] npart     Number of mesh partitions.
     * \param [in] visu      Pointer to the Visu object
     *
     */
 
    Mesh(const MPI_Comm &localComm,
              Visu* visu,
              int npart,
              CWP_Displacement_t displacement);


    /**
     * \brief Mesh destructor
     *
     * Destroy the Mesh object.
     *
     */

     virtual ~Mesh();
     
     /**
     * \brief Set the mesh coordinates
     *
     * \param [in] i_part       Index of the mesh partition
     * \param [in] n_pts        Number of partition vertices
     * \param [in] coords       Array containing the verticies coordinates
     * \param [in] global_num   Global numbering of the vertices
     *
     */
   
     void nodal_coord_set(const int   i_part,
                          const int   n_pts,
                          double      coords[],
                          CWP_g_num_t global_num[]); 
                          
    /**
     * \brief Mesh deletion and free memory
     *
     */
   
     void meshDel();
     
     /**
     * \brief Addition of a block (set of cells) to the mesh partition
     *
     * This function add a block to the geometric mesh.
     *
     * \param [in] block_type  Type of the block addition     
     *
     * \return block_id  Block Identifier
     */

     int blockAdd(const CWP_Block_t  block_type);
                 
                 
     /**
     * \brief Set a standard block to the interface mesh
     *
     * \param [in] i_part      Partition identifier
     * \param [in] block_id    Block identifier  
     * \param [in] n_elts      Number of block elements
     * \param [in] connec      Vertices to elements connectivity 
     * \param [in] global_num  Global numbering of the vertices in the block (or NULL)
     *
     */
                 
     void stdBlockSet( const int              i_part,
                       const int              block_id,
                       const int              n_elts,
                       int                    connec[], 
                       CWP_g_num_t            global_num[]
                     );
                     
     /**
     * \brief Set a face polygon block to the interface mesh
     *
     * \param [in] i_part      Partition identifier
     * \param [in] block_id    Block identifier
     * \param [in] n_elts      Number of block elements
     * \param [in] connec_idx  Vertices to elements connectivity index
     * \param [in] connec      Vertices to elements connectivity 
     * \param [in] global_num  Global numbering of the vertices in the block (or NULL)
     *
     */
                 
     void poly2DBlockSet( const int              i_part,
                          const int              block_id,
                          const int              n_elts,
                          int                    connec_idx[],
                          int                    connec[], 
                          CWP_g_num_t            global_num[]
                        );


     /**
     * \brief Set a face polhedron block to the interface mesh
     *
     * \param [in] i_part            Partition identifier
     * \param [in] block_id          Block identifier
     * \param [in] n_elts            Number of block elements
     * \param [in] n_faces           Number of faces elements (or NULL)
     * \param [in] connec_faces_idx  Vertices to faces connectivity index
     * \param [in] connec_faces      Vertices to faces connectivity 
     * \param [in] connec_cells_idx  Faces to cells connectivity index
     * \param [in] connec_cells      Faces to cells connectivity 
     * \param [in] global_num        Global numbering of the vertices in the block (or NULL)
     *
     */
                 
     void poly3DBlockSet( const int              i_part,
                          const int              block_id,
                          const int              n_elts,
                          const int              n_faces,
                          int                    connec_faces_idx[],
                          int                    connec_faces[],
                          int                    connec_cells_idx[],
                          int                    connec_cells[], 
                          CWP_g_num_t            global_num[]
                        );               

     /**
     * \brief Get the Mesh Block Type
     *
     * This function gets the mesh block type.
     *
     * \param [in] id_block    Block identifier
     * \return CWP_Block_t Mesh  Block Type   
     *
     */

     CWP_Block_t Mesh_nodal_block_type_get(const int   id_block
                                          );
                                          
    /**
     * \brief Adding a polyhedron block to the geometric mesh from
     * a face-to-cell connectivity and a vertices-to-faces connectivity.
     *
     * This function add a polyhedron 3D block to the geometric mesh from
     * a face-to-cell connectivity and a vertices-to-faces connectivity.
     *
     * \param [in]  i_part            Current partition
     * \param [in]  n_cells           Number of elements
     * \param [in]  cell_face_idx     Polyhedron to face index 
     *                                (src_poly_cell_face_idx[0] = 0 and
     *                                 size = n_elts + 1)
     * \param [in]  cell_face         Polyhedron to face connectivity 
     *                                (size = cell_face_idx[n_elts])
     * \param [in]  n_faces           Number of faces      
     * \param [in]  face_vtx_idx      Polyhedron vertices to faces index 
     *                                (face_vtx_idx[0] = 0 and
     *                                 size_idx = max(face_vtx) + 1)
     * \param [in]  face_vtx          Polyhedron vertices to faces connectivity
     *                                (size = face_vtx_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)
     *
     */
     
     void fromCellFaceSet(const int   i_part,
                          const int   n_cells,
                          int         cell_face_idx[],
                          int         cell_face[],
                          int         n_faces,
                          int         face_vtx_idx[],
                          int         face_vtx[],
                          CWP_g_num_t parent_num[]);  

    /**
     * \brief Adding a polygon 2D block to the geometric mesh from
     * a vertices-to-faces connectivity and a edge-to-face connectivity.
     *
     * This function add a polygon 2D block to the geometric mesh from
     * a vertices-to-faces connectivity and a edge-to-face connectivity.
     *
     * \param [in]  i_part            Current partition
     * \param [in]  n_faces           Number of faces      
     * \param [in]  face_edge_idx     Polygon vertices to faces index 
     *                                (face_edge_idx[0] = 0 and
     *                                size_idx = max(face_edge) + 1)
     * \param [in]  face_edge         Polyhegon vertices to face connectivity
     *                                (size = face_edge_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)
     * \param [in]  n_edges           Number of edges      
     * \param [in]  edge_vtx_idx      Vertices to edges connectivity index 
     *                                (edge_vtx_idx[0] = 0 and
     *                                size_idx = max(edge_vtx) + 1)
     * \param [in]  edge_vtx          Polygon vertices to edges connectivity
     *                                (size = edge_vtx_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)     
     *
     */
                       
     void fromFacesEdgeSet(const int   i_part,
                           const int   n_faces,
                           int         face_edge_idx[],
                           int         face_edge[],
                           const int   n_edges,
                           int         edge_vtx_idx[],
                           int         edge_vtx[],
                           CWP_g_num_t parent_num[]);
                               

    /**
    * \brief Update the block database.
    *
    * This function updates the block database .
    *
    *
    */

    void updateBlockDB();

    /**
    * \brief Get the vertices number of the partition i_part
    *
    * This function gets the vertices number of the partition i_part.
    *
    * \param [in] i_part   Partition identifier
    * \return Vertex number of i_part partition.   
    *
    */

    inline int getPartNVertex(int i_part) const;

    /**
    * \brief Get the vertices coordinates of the i_part partition
    *
    * This function gets the vertices number of the i_part partition.
    *
    * \param [in] i_part   Partition identifier
    * \return Vertices coordinates of the i_part partition.   
    *
    */

    inline double* getVertexCoords(int i_part);


    inline CWP_g_num_t* getVertexGNum(int i_part);


    /**
    * \brief Get the paradigm nodal object associated with the mesh.
    *
    * This function gets the paradigm nodal object associated with the mesh.
    *
    * \return Paradigm nodal object associated with the mesh.  
    *
    */

    inline PDM_Mesh_nodal_t& getPdmNodal() const;

    /**
    * \brief Get the number of elements of the id_block block.
    *
    * This function gets the number of elements of the id_block block.
    *
    * \param [in] id_block   Block identifier
    * \return                Number of elements of the id_block block
    *
    */

    inline int getBlockNElts(int id_block, int i_part);

    /**
    * \brief Get the number of elements of the id_part partition
    *
    * This function gets the number of elements of the id_part partition.
    *
    * \param [in] id_part    Partition identifier
    * \return                Number of elements of the id_part partition
    *
    */

    inline int getPartNElts(int id_part) const;

    /**
    * \brief Get the number of polyhedra in a partition
    *
    * This function gets the number of polyhedra in a partition.
    *
    * \param [in] i_part   Partition identifier
    * \return              Number of polydra of the i_part partition
    *
    */

    inline int getPartNPolyhedra(int i_part) const;

    /**
    * \brief Get a block element connectivity index
    *
    * This function gets the element connectivity index of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Element connectivity index of the id_block block
    *
    */

    inline int* getEltConnectivityIndex(int id_block,int i_part);

    /**
    * \brief Get a block element connectivity
    *
    * This function gets the element connectivity of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Element connectivity of the id_block block
    *
    */
    
    inline int* getEltConnectivity     (int id_block,int i_part);
 
     /**
    * \brief Get a block face connectivity index
    *
    * This function gets the face connectivity index of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Face connectivity index of the id_block block
    *
    */
    
    inline int* getFaceConnectivityIndex(int id_block);

    /**
    * \brief Get a block face connectivity
    *
    * This function gets the face connectivity of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Face connectivity of the id_block block
    *
    */

    inline int* getFaceConnectivity     (int id_block);
    
    /**
    * \brief True if coordinates are defined on all partitions False otherwise.
    *
    *
    */
    
    inline bool coordsDefined ();

   // inline const std::vector<double>& getVolume();

   // inline const std::vector<double>& getCellCenterCoords();

   /**
    * \brief Return pdmNodal Handler identifier
    *
    *
    */
   inline int getPdmNodalIndex();
   
   /**
    * \brief Return Coordinates of the i_part partition.
    *
    * \param [in] i_part partition identifier.
    * \return Coordinates of the i_part partition.
    *
    */
    
   inline double* getCoordinates(int i_part);
   

   /**
    * \brief Return MPI communicator
    *
    *
    */
    
   inline const MPI_Comm getMPIComm();
   
   /**
    * \brief Return MPI communicator pointer
    *
    *
    */
    
   inline const MPI_Comm* getMPICommP();
   
   /**
    * \brief Return number of partitions.
    *
    *
    */
    
   inline int getNPart();

   void geomFinalize();



   /**
    * \brief Set the Visu pointer object
    *
    * \param [in] visu Pointer to the Visu object
    *
    */
    
   inline void setVisu(Visu* visu);

   CWP_g_num_t* globalNumGet(int id_block,int i_part) {
      std::map<int,cwipi::Block*>::iterator it;
      it = _blockDB.find(id_block);
   /*   if(it != _blockDB.end()) 
        {Block* toto = it->second;
         CWP_g_num_t* tata = toto -> GNumMeshGet(i_part,0);
         
         return tata;
        }
       else
         printf("Pas trouvé\n");*/
   }
   
   void connecCompute(int i_part);   
   int* connecIdxGet(int i_part);
   int* connecGet(int i_part);
   
   int GNVerticeGet(int i_part);
   int GNEltGet(int i_part);
   CWP_g_num_t* GNumEltsGet(int i_part);
   double* eltCentersGet(int i_part);
  void eltCentersCompute(int i_part);
   
   int* blockDBGet() {
     return _blocks_id;
   }

   int nBlockGet() {
     return _nBlocks;
   } 


   CWP_Block_t blockTypeGet(int id_block) {
     return _blockDB[id_block] -> blockTypeGet(); 
   } 

   CWP_g_num_t* gnumMeshBlockGet(int id_block,int i_part) {
     return _blockDB[id_block] -> GNumMeshGet(i_part);
   } 


   CWP_g_num_t* gnumInsideBlockGet(int id_block,int i_part) {
     return _blockDB[id_block] -> GNumBlockGet(i_part);
   } 
   
  private:
    
    const MPI_Comm                          &_localComm;              /*!< Communicator */
    PDM_MPI_Comm                             _pdm_localComm;
    int                                     _nDim;                   /*!< Entities dimensions */  
    int                                     _nBlocks;                /*!< Number of blocks of the mesh */
    int*                                    _blocks_id;              /*!< List of block identifiers */
    int                                     _order;                  /*!< Mesh order */
    std::vector<int>                        _nVertex;                /*!< Number of vertices for each partition  */
    std::vector<int>                        _nElts;                  /*!< Number of elements for each partition  */
    std::vector<double*>                    _coords;                 /*!< Vertices coordinate for each partition  */
    std::vector<int*>                       _connec_idx;
    std::vector<int*>                       _connec;
    std::vector<CWP_g_num_t*>               _gnum_elt;
    std::vector<double*>                    _elt_centers;
    
    std::map< int, CWP_g_num_t*>            _global_num;             /*!< Global numbering for each partition  */
    int                                     _npart;                  /*!< Number of partition  */
    int                                     _pdmNodal_handle_index;  /*!< Mesh (nodal) index for paradigm handler */
    int                                     _pdmGNum_handle_index;   /*!< Global number index for paradigm handler   */
    PDM_Mesh_nodal_t                       *_pdmNodal;               /*!< Pointer to the paradigm mesh (nodal) object   */
    std::map<int,cwipi::Block*>             _blockDB;                /*!< Blocks database  */
    Visu                                   *_visu;                   /*!< Pointer to the Visu object */
    std::map<int,int>                       _id_visu;                /*!< Map of the PDM_Writer block identifier */  
    CWP_Displacement_t                      _displacement;          /*!< Type of mesh displacement */  
    
  //   Mesh &operator=(const Mesh &other);  /*!< Assigment operator not available */
  //   Mesh (const Mesh& other);            /*!< Copy constructor not available */

  };

 


  void Mesh::setVisu(Visu* visu) {
    _visu = visu;
  }

  int Mesh::getNPart() {
    return _npart;
  }

  double* Mesh::getCoordinates(int i_part) {
    return _coords[i_part];
  }

  int Mesh::getPdmNodalIndex() {
    return _pdmNodal_handle_index;
  }
 
  const MPI_Comm Mesh::getMPIComm() {
    return _localComm;
  }

  const MPI_Comm* Mesh::getMPICommP() {
    return &_localComm;
  }


  bool Mesh::coordsDefined () {
  
    bool ans = true;
    /*for (int i=0;i<_npart;i++) {
      if(_coords[i]==NULL) {
        printf("coorde de %i is NULL\n",i);
        ans = false;
      }
    }*/
    return (_coords.size()==_npart);
  }


  int Mesh::getPartNVertex(int i_part) const
  {
    return _nVertex[i_part];
  }

  double* Mesh::getVertexCoords(int i_part)
  {
    return _coords[i_part];
  }

  CWP_g_num_t* Mesh::getVertexGNum(int i_part)
  {
    return _global_num[i_part];
  }


  PDM_Mesh_nodal_t& Mesh::getPdmNodal() const
  {
    return *_pdmNodal;
  }


  int Mesh::getBlockNElts(int id_block,int i_part)
  {
    return _blockDB[id_block] -> NEltsGet()[i_part];
  }

  int Mesh::getPartNElts(int id_part) const
  {
    return _nElts[id_part];
  }

  int Mesh::getPartNPolyhedra(int i_part) const
  {
   // return _nPolyhedra[i_part];
  }

  int* Mesh::getEltConnectivityIndex(int id_block,int i_part)
  {
    return _blockDB[id_block] -> ConnecIDXGet()[i_part];
  }

  int* Mesh::getEltConnectivity(int id_block,int i_part)
  { 
    return _blockDB[id_block] -> ConnecGet()[i_part];
  }
  
  
  inline int* Mesh::getFaceConnectivity(int id_block)
  { 
 //   return _blocks[id_block]._connec_faces;
  }

  inline int* Mesh::getFaceConnectivityIndex(int id_block)
  { 
  //  return _blocks[id_block]._connec_faces_idx;
  }





}

#endif //__Mesh_H__
