#ifndef __blockCP_H__
#define __blockCP_H__
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

#include "block.hxx"
#include "cwp.h"
#include <map>
#include <vector>

namespace cwipi {

  /** 
   * \class BlockCP blockCP.hxx "blockCP.hxx"
   * \brief Polyhedron 3D Block Mesh definition
   *
   *  This class defines the polyhedron 3D Block objects.
   *  It inherits from the Block class through a Factory.
   * 
   */

  class BlockCP :
    public Block {
    
    public:
        
     /**
      * \brief Destructor.
      *
      */
    
      BlockCP();


     /**
      * \brief Destructor.
      *
      */
    
      virtual ~BlockCP();
    
       /**
       * \brief Set a CWIPI block in a partition
       * 
       * \param [in]  i_part     Partition identifier
       * \param [in]  n_elts     Number of elements of the block in the partition
       * \param [in]  n_faces    Number of faces of the block in the partition
       * \param [in]  connec_faces_idx Vertices to faces connectivity index
       * \param [in]  connec_faces     Vertices to faces connectivity
       * \param [in]  connec_cells_idx Faces to cells connectivity index
       * \param [in]  connec_cells     Faces to cells connectivity
       * \param [in]  global_num Mesh  Global numbering of the block
       *
       */ 
             
       virtual void blockSet(int i_part,int n_elts,int n_faces,
                             int* connec_faces_idx, 
                             int* connec_faces,
                             int* connec_cells_idx,
                             int* connec_cells,
                             CWP_g_num_t* global_num);


       /**
        *
        * \brief return the element connectivity (Standard or Face_Poly_2D CWP_Block_t) or cells-faces connectiviy (Cells_POLY_3D)
        * for each partition.
        * 
        *
        */
     
        inline virtual std::map<int,int*>  ConnecGet();

       /**
        *
        * \brief return the element connectivity index (Face_Poly_2D CWP_Block_t) or cells-faces connectivity index (Cells_POLY_3D)
        * for each partition.
        * 
        *
        */
     
        inline virtual std::map<int,int*>  ConnecIDXGet();


       /**
        *
        * \brief return the element faces connectivity (Cells_POLY_3D)
        * for each partition.
        * 
        *
        */
     
        inline virtual std::map<int,int*>  ConnecFacesGet();


       /**
        *
        * \brief return the element faces connectivity index (Cells_POLY_3D)
        * for each partition.
        * 
        *
        */
     
        inline virtual std::map<int,int*>  ConnecFacesIDXGet();

       /**
        *
        * \brief return the number of faces for each partition (Cells_POLY_3D)
        * 
        *
        */
     
        inline virtual std::map<int,int >  NFacesGet();

        void geomFinalize();

    private:
      std::map<int,int >          _n_faces;             /*!< Number of faces for each partition */
      std::map<int,int*>          _connec_faces_idx;    /*!< Faces connectivity Index for each partition */
      std::map<int,int*>          _connec_faces;        /*!< Faces connectivity for each partition */
      std::map<int,int*>          _connec_cells_idx;    /*!< Cells onnectivity Index for each partition */
      std::map<int,int*>          _connec_cells;        /*!< Cells connectivity for each partition */    
      
  }; //BlockCP Class


  std::map<int,int*>  BlockCP::ConnecGet() {
    return _connec_cells;
  }

  std::map<int,int*>  BlockCP::ConnecIDXGet() {
  
    return _connec_cells_idx;
  }

  std::map<int,int*>  BlockCP::ConnecFacesGet() {
    return _connec_faces;
  }

  std::map<int,int*>  BlockCP::ConnecFacesIDXGet() {
    return _connec_faces_idx;
  }

  std::map<int,int >  BlockCP::NFacesGet() {
    return _n_faces;
  }



} //namespace cwipi

#endif //blockCP

