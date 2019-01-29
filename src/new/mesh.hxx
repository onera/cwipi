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


class id_part_block{

  public:
    inline id_part_block(int id_part,int id_block):_id_part(id_part),_id_block(id_block)
      {
      }
      
    inline ~id_part_block()
      {
      }
      
    inline int get_id_part()
      {return _id_part;
      }
      
    inline int get_id_block()
      {return _id_block;
      }

      
    inline bool operator < (id_part_block const &a) const
      {if(_id_part<a._id_part) return true;
       else
         {if(_id_part>a._id_part) return false;
          else
            {if(_id_block>a._id_block) return false;
             else return true;
            }
         }
      }
      
    inline bool operator > (id_part_block const &a) const
      {if(_id_part>a._id_part) return true;
       else
         {if(_id_part<a._id_part) return false;
          else
            {if(_id_block<a._id_block) return false;
             else return true;
            }
         }
      }
      
    inline bool operator == (id_part_block const &a)
      {if(_id_part==a._id_part and _id_block==a._id_block) return true;
       else return false;
      }
        
    int _id_part;
    int _id_block;
  private:
} ;





namespace cwipi {

  /** 
   * \class Mesh Mesh.hxx "Mesh.hxx"
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
     * \param [in] npart Number of mesh partitions.
     *
     */
 
    Mesh(const MPI_Comm &localComm,
         int npart);


    /**
     * \brief Mesh destructor
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
     * \param [in] ipart       Index of the mesh partition
     * \param [in] block_type  Type of the block addition     
     * \param [in] block_type  Block type i.e. Type of the block elements
     * \param [in] n_elts      Number of block elements
     * \param [in] connec_idx  Vertices to elements connecivity index
     * \param [in] connec      Vertices to elements connecivity 
     * \param [in] global_num  Global numbering of the vertices in the block
     * \param [in] parent_num  Parent numbering in the block
     */
   
     void blockAdd(const int                  i_part,
                       const CWP_Block_t      block_type,
                       const int              n_elts,
                       int                    connec_idx[],
                       int                    connec[],
                       const int              n_faces,
                       int                    face_vtx_idx[],
                       int                    face_vtx[],    
                       CWP_g_num_t            global_num[],
                       CWP_g_num_t            parent_num[]);

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

                               
             
  private:
    
    // TODO: renommer _nDim par entitesDim
    const MPI_Comm                         &_localComm;
    int                                     _nDim;
    int                                     _order;
    std::vector<int>                        _nVertex;
    std::vector<int>                        _nFaces;
    std::vector<int>                        _nElts;
    std::vector<double*>                    _coords;
    std::map< int, CWP_g_num_t*>            _global_num;
    std::map< id_part_block, CWP_g_num_t*>  _global_num_block;
    std::map< int, CWP_g_num_t*>            _parent_num;
    std::map< id_part_block, CWP_g_num_t*>  _parent_num_block;
    int                                     _npart;
    int                                     _pdmNodal_handle_index;
    int                                     _pdmGNum_handle_index;
    bool                                    _isNodalFinalized;
    PDM_Mesh_nodal_t                       *_pdmNodal;
    std::map< id_part_block, int* >         _connec;
    std::map< id_part_block, int* >         _connec_idx;
    std::map< id_part_block, int* >         _connec_faces;
    std::map< id_part_block, int* >         _connec_faces_idx;
    
    
    
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
  //   Mesh                                *Mesh;
  };

}

#endif //__Mesh_H__
