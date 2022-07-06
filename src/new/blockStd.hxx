#ifndef __BlockStd_H__
#define __BlockStd_H__
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
#include <cstdlib>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class BlockStd blockStd.hxx "blockStd.hxx"
   * \brief Standard Block Mesh definition
   *
   *  This class defines the standard Block objects.
   *  It inherits from the Block class through a Factory.
   *
   */

  class BlockStd :
    public Block {

    public:


     /**
      * \brief Destructor.
      *
      */

      BlockStd();

     /**
      * \brief Destructor.
      *
      */

      ~BlockStd();

      /**
       * \brief Set a CWIPI block in a partition
       *
       * \param [in]  i_part     Partition identifier
       * \param [in]  n_elts     Number of elements of the block in the partition.
       * \param [in]  connec     Elements connectivity.
       * \param [in]  global_num Mesh global numbering of the block
       *
       */

       void blockSet(int i_part,int n_elts,int* connec,CWP_g_num_t* global_num);

       /**
        *
        * \brief return the element connectivity (Standard or Face_Poly_2D CWP_Block_t) or cells-faces connectiviy (Cells_POLY_3D)
        * for each partition.
        *
        *
        */

        inline std::map<int,int*>  ConnecGet();

        inline int*  ConnecGet(int i_part);

        inline int*  ConnecIDXGet(int i_part);

        bool  gnumRequired(){
           for(int i_part = 0; i_part<_n_part; i_part++){
             if(_global_num [i_part] == NULL )
               return true;
           }

           return false;
        }

        void geomFinalize();

    private:

      std::map<int,int*>          _connec;              /*!< Connectivity for each partition */
      std::map<int,int*>          _connec_idx;           /*!< Connectivity index for each partition */

  };



  int*  BlockStd::ConnecGet(int i_part) {

    return _connec[i_part];
  }


  int*  BlockStd::ConnecIDXGet(int i_part) {
     int n_vtx_per_cell = 0;
     if (_blockType == CWP_BLOCK_FACE_TRIA3)
       n_vtx_per_cell = 3;
     else if (_blockType == CWP_BLOCK_FACE_QUAD4)
       n_vtx_per_cell = 4;
     else if (_blockType == CWP_BLOCK_CELL_TETRA4)
       n_vtx_per_cell = 4;
     else if (_blockType == CWP_BLOCK_CELL_PYRAM5)
       n_vtx_per_cell = 5;
     else if (_blockType == CWP_BLOCK_CELL_PRISM6)
       n_vtx_per_cell = 6;
     else if (_blockType == CWP_BLOCK_CELL_HEXA8)
       n_vtx_per_cell = 8;
     else{
       PDM_error(__FILE__, __LINE__, 0, "Unknown block type in BlockStd::ConnecIDXGet.\n");
     }

     _connec_idx[i_part] = (int*)malloc(sizeof(int)*(1 + _n_elt[i_part]));
     for (int i=0; i<_n_elt[i_part]+1; i++)
        _connec_idx[i_part][i] = n_vtx_per_cell * i;

     return _connec_idx[i_part];
  }


  std::map<int,int*>  BlockStd::ConnecGet() {

    return _connec;
  }

}

/**
 * \endcond
 */


#endif //BlockStd
