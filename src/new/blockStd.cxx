

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

#include "blockStd.hxx"
#include "cwp.h"
#include <pdm_mesh_nodal.h>
#include <pdm_gnum.h>
#include <pdm_error.h>
#include <pdm_printf.h>
#include <map>
#include <vector>

#include <mesh.hxx>

/**
 * \cond
 */

namespace cwipi {
  BlockStd::BlockStd()
     :Block::Block()
  {

  }


  BlockStd::~BlockStd()
  {
    _connec.clear();
    _connec_idx.clear();
  }

  void BlockStd::blockSet(int i_part,int n_elt,int* connec,CWP_g_num_t* mesh_global_num){

     double* _cells_center_part = (double*)malloc (sizeof(double) * 3 * n_elt);

     _global_num [i_part] = mesh_global_num;

     _isSet[i_part] = true;
     _n_elt[i_part] = n_elt;
     _part_id.push_back(i_part);
     _n_part_def=_n_part_def+1;
     _cells_center[i_part] = _cells_center_part;
     _connec.insert    ( std::pair < int, int* > (i_part,connec));
  }


  void BlockStd::geomFinalize(){

    for(int i_part = 0; i_part<_n_part; i_part++){

      Visu* visu = ((Mesh*)_mesh)->getVisu();
      if(visu->isCreated() && ((Mesh*)_mesh)->getDisplacement() == CWP_DYNAMIC_MESH_STATIC) {
        visu->GeomBlockStdSet ( ((Mesh*)_mesh)->getIdVisu( _block_id_cwipi ),
                                  i_part,
                                  _n_elt[i_part] ,
                                  _connec[i_part],
                                  _global_num [i_part]);
      }
    }
  }


}


/**
 * \endcond
 */
