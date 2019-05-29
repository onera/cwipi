#ifndef __BlockFP_H__
#define __BlockFP_H__
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
   * \class BlockFP blockFP.hxx "blockFP.hxx"
   * \brief Polygon 2D Block Mesh definition
   *
   *  This class defines the polygon 2D Block objects.
   *  It inherits from the Block class through a Factory.
   * 
   */

  class BlockFP :
    public Block {
    
    public:
    
    
     /**
      * \brief Destructor.
      *
      */
      
      BlockFP();


     /**
      * \brief Destructor.
      *
      */
    
      virtual ~BlockFP();
     
      /**
       * \brief Set a CWIPI block in a partition 
       * 
       * \param [in]  i_part     Partition identifier
       * \param [in]  n_elts     Number of elements of the block in the partition.
       * \param [in]  connec_idx Elements connectivity index
       * \param [in]  connec     Elements connectivity
       * \param [in]  global_num Mesh global numbering of the block
       *
       */    
 
       virtual void blockSet(int i_part,int n_elts,
                             int* connec_idx,
                             int* connec,
                             CWP_g_num_t* mesh_global_num);
                             
      /**
       * \brief Add and Set the CWIPI block from a Paradigm block 
       * 
       * \param [in] pdm_id_block A block identifier from Paradigm
       * \param [in] mesh         Pointer to the Mesh object owning the block.
       *
       *
       */    
              
       virtual void FromPDMBlock(int pdm_id_block, void* mesh);

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
        
           
    private:
      std::map<int,int*>          _connec_idx;          /*!< Connectivity Index for each partition */
      std::map<int,int*>          _connec;              /*!< Connectivity for each partition */
     
  };


  std::map<int,int*>  BlockFP::ConnecGet() {
  
    return _connec;
  }

  std::map<int,int*>  BlockFP::ConnecIDXGet() {
  
    return _connec_idx;
  }

}


#endif //BlockFP

