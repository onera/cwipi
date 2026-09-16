/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2026  ONERA

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

#include "blockFP.hxx"

#include "cwp.h"
#include "mesh.hxx"

#include "pdm_gnum.h"

#include <vector>

/**
 * \cond
 */

namespace cwipi {

  BlockFP::BlockFP()
    : Block::Block()
  {
  }


  BlockFP::~BlockFP()
  {
  }


  void
  BlockFP::BlockAdd
  (
    CWP_Block_t  blockType,
    Mesh        *mesh
  )
  {
    Block::BlockAdd(blockType, mesh);
    _connec_idx.resize(_n_part, nullptr);
    _connec    .resize(_n_part, nullptr);
  }


  void
  BlockFP::blockSet
  (
    int          i_part,
    int          n_elt,
    int         *connec_idx,
    int         *connec,
    CWP_g_num_t *global_num
  )
  {
    _global_num[i_part] = global_num;
    _n_elt     [i_part] = n_elt;
    _connec    [i_part] = connec;
    _connec_idx[i_part] = connec_idx;
  }


  void
  BlockFP::blockGet
  (
    int           i_part,
    int          *n_elts,
    int         **connec_idx,
    int         **connec,
    CWP_g_num_t **global_num
  )
  {
    *global_num = _global_num[i_part];
    *n_elts     = _n_elt     [i_part];
    *connec     = _connec    [i_part];
    *connec_idx = _connec_idx[i_part];
  }


  void
  BlockFP::geomFinalize()
  {
  }

}

/**
 * \endcond
 */