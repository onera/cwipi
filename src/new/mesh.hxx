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

#include "cwp.h"




namespace cwipi {


typedef struct {
    int              _id_part;
    int              _id_block;
    CWP_Block_t      _blockType; 
    int              _nElts;
    int*             _connec_idx;
    int*             _connec;
    int              _nFaces;
    int*             _n_faces;
    int*             _connec_faces_idx;  
    int*             _connec_faces;
    CWP_g_num_t*     _global_num_block;
    int*             _parent_num;
    bool             _isBlockFinalized;
} _block;



  /** 
   * \class Mesh 
   *        Mesh.hxx 
   *        "Mesh.hxx"
   *
   * \brief Geometry Mesh
   *
   *  This class computes defines th geometry Mesh (mesh)
   * 
   */

  class Mesh {
    
  public:

    /**
     * \brief Mesh constructor
     * 
     * Construct the CWIPI mesh by using paradigm nodal methods.
     * 
     * \param [in] localComm Coupling Communicator.
     * \param [in] npart     Number of mesh partitions.
     *
     */
 
    Mesh(const MPI_Comm &localComm,
         int npart);


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
     * \param [in] ipart       Index of the mesh partition
     * \param [in] n_vtx       Number of partition vertices
     * \param [in] coords      Array containing the verticies coordinates
     * \param [in] global_num  Global numbering of the vertices
     *
     */
   
     void nodal_coord_set(const int   i_part,
                          const int   n_pts,
                          double      coords[],
                          CWP_g_num_t global_num[]); 
                          
     /**
     * \brief Finalize the mesh setting after blocks and coordinates additions
     *
     */
       
     void endSet();
   
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
     * \param [in] i_part       Index of the mesh partition
     * \param [in] block_type  Type of the block addition     
     * \param [in] n_elts      Number of block elements
     * \param [in] connec_idx  Vertices to elements connectivity index
     * \param [in] connec      Vertices to elements connectivity 
     * \param [in] n_faces     Number of faces elements (or NULL)
     * \param [in] face_vtx_idx Vertices to faces connectivity index (or NULL)
     * \param [in] face_vtx      Vertices to faces connectivity (or NULL)
     * \param [in] global_num  Global numbering of the vertices in the block (or NULL)
     * \param [in] parent_num  Parent numbering in the block (or NULL)
     *
     */

     int blockAdd(const int                  i_part,
                       const CWP_Block_t      block_type,
                       const int              n_elts,
                       int                    connec_idx[],
                       int                    connec[],
                       const int              n_faces,
                       int                    face_vtx_idx[],
                       int                    face_vtx[],    
                       CWP_g_num_t            global_num[],
                       int                    parent_num[]);

     CWP_Block_t Mesh_nodal_block_type_get(const int   id_block
                                          );
                                          
                                          
     void fromFacesEdgeSet(const int   i_part,
                               const int   n_faces,
                               int         face_edge_idx[],
                               int         face_edge[],
                               const int   n_edges,
                               int         edge_vtx_idx[],
                               int         edge_vtx[],
                               CWP_g_num_t parent_num[]);
                               
     void fromCellFaceSet(const int   i_part,
                          const int   n_cells,
                          int         cell_face_idx[],
                          int         cell_face[],
                          int         n_faces,
                          int         face_vtx_idx[],
                          int         face_vtx[],
                          CWP_g_num_t parent_num[]); 


    void updateBlockDB(int id_part);

    inline const std::vector<int>& getNVertex(int i_part) const;

    inline const std::vector<double*> getVertexCoords() const;

    inline PDM_Mesh_nodal_t& getPdmNodal() const;

    inline int getBlockNElts(int id_block);

    inline const int getPartNElts(int id_part) const;

    inline const std::vector<int>& getNPolyhedra(int i_part) const;

    inline int* getEltConnectivityIndex(int id_block);

    inline int* getEltConnectivity     (int i_block);
    
    inline int* getFaceConnectivityIndex(int id_block);

    inline int* getFaceConnectivity     (int i_block);
    

   // inline const std::vector<double>& getVolume();

   // inline const std::vector<double>& getCellCenterCoords();


  private:
    
    // TODO: renommer _nDim par entitesDim
    const MPI_Comm                         &_localComm;
    int                                     _nDim;
    int                                     _nBlocks;
    int*                                    _blocks_id;
    std::vector<int*>                       _blocks_id_part;
    int                                     _order;
    std::vector<int>                        _nVertex;
    std::vector<int>                        _nElts;    
    std::vector<double*>                    _coords;
    std::map< int, CWP_g_num_t*>            _global_num;
    std::map< int, CWP_g_num_t*>            _parent_num;
    int                                     _npart;
    int                                     _pdmNodal_handle_index;
    int                                     _pdmGNum_handle_index;
    bool                                    _isNodalFinalized;
    PDM_Mesh_nodal_t                       *_pdmNodal;
    std::map<int,_block>                    _blocks;

    
    
    
  //   Mesh &operator=(const Mesh &other);  /*!< Assigment operator not available */
  //   Mesh (const Mesh& other);            /*!< Copy constructor not available */

  // private:
  //   std::vector<double>         *_tmpVertexField;  // Evite une allocation a chaque extrapolation
  //   std::vector<double>         *_tmpDistantField; //TODO: Fusionner _tmpDistantField utiliser pour exchange
  //                                          // et les comm asynchrones
  //   std::map<int, std::vector<double> * > &_tmpDistantFieldsIssend; //TODO: temporaire : A revoir lors
  //                                                                   // de la restructuration
  //   std::map<int, const double * >        &_tmpLocalFieldsIrecv;
  //   std::map<int, const char * >          &_tmpExchangeNameIrecv;
  //   std::map<int, int >                   &_tmpStrideIrecv;
  //   std::map<int, int >                   &_tmpTimeStepIrecv;
  //   std::map<int, double >                &_tmpTimeValueIrecv;
  //   std::map<int, const char * >          &_tmpFieldNameIrecv;
  };

  const std::vector<int>& Mesh::getNVertex(int i_part) const
  {
    return _nVertex;
  }

  const std::vector<double*> Mesh::getVertexCoords()  const
  {
    return _coords;
  }

  PDM_Mesh_nodal_t& Mesh::getPdmNodal() const
  {
    return *_pdmNodal;
  }


  int Mesh::getBlockNElts(int id_block)
  {
    return _blocks[id_block]._nElts;
  }

  const int Mesh::getPartNElts(int id_part) const
  {
    return _nElts[id_part];
  }

  const std::vector<int>& Mesh::getNPolyhedra(int i_part) const
  {
   // return _nPolyhedra[i_part];
  }

  inline int* Mesh::getEltConnectivityIndex(int id_block)
  {
    return _blocks[id_block]._connec_idx;
  }

  inline int* Mesh::getEltConnectivity(int id_block)
  { 
    return _blocks[id_block]._connec;
  }
  
  
  inline int* Mesh::getFaceConnectivity(int id_block)
  { 
    return _blocks[id_block]._connec_faces;
  }

  inline int* Mesh::getFaceConnectivityIndex(int id_block)
  { 
    return _blocks[id_block]._connec_faces_idx;
  }





}

#endif //__Mesh_H__
