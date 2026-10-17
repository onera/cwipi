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

#include "blockHO.hxx"

#include "cwp.h"

#include "mesh.hxx"

#include "pdm_error.h"
#include "pdm_ho_ordering.h"
#include "pdm_mesh_nodal.h"

#include <cstdlib>
#include <cstring>
#include <vector>

/**
 * \cond
 */

namespace cwipi {

  BlockHO::BlockHO()
  : BlockStd::BlockStd()
  {
    _order       = 0;
    _ho_ordering = nullptr;
  }


  BlockHO::~BlockHO()
  {
    _connec.clear();
    for (size_t i = 0; i < _connec_ijk.size(); i++) {
      if (_connec_ijk[i] != nullptr) {
        free(_connec_ijk[i]);
      }
    }
    _connec_ijk.clear();
    if (_ho_ordering != nullptr) {
      free(_ho_ordering);
    }
  }


  void
  BlockHO::blockAdd
  (
    CWP_Block_t  block_type,
    Mesh        *mesh
  )
  {
    BlockStd::blockAdd(block_type, mesh);
    _connec    .resize(_n_part, nullptr);
    _connec_ijk.resize(_n_part, nullptr);
  }


  void
  BlockHO::blockSet
  (
          int          i_part,
          int          n_elts,
          int         *connec,
          CWP_g_num_t *global_num,
    const int          order,
    const char        *ho_ordering
  )
  {
    _n_elts    [i_part] = n_elts;
    _connec    [i_part] = connec;
    _global_num[i_part] = global_num;

    PDM_Mesh_nodal_elt_t elt_type  = static_cast<PDM_Mesh_nodal_elt_t>(_block_type);
    int                  elt_n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);
    _connec_ijk[i_part] = static_cast<int *>(malloc(sizeof(int) * n_elts * elt_n_vtx));
    memcpy(_connec_ijk[i_part], _connec[i_part], sizeof(int) * n_elts * elt_n_vtx);

    _order = order;
    if (_ho_ordering != nullptr) {
      free(_ho_ordering);
      _ho_ordering = nullptr;
    }
    if (ho_ordering != nullptr) {
      _ho_ordering = static_cast<char*>(malloc(sizeof(char) * (strlen(ho_ordering) + 1)));
      strcpy(_ho_ordering, ho_ordering);
    }
  }


  void
  BlockHO::blockGet
  (
    int           i_part,
    int          *n_elts,
    int         **connec,
    CWP_g_num_t **global_num,
    int          *order,
    char        **ho_ordering
  )
  {
    *n_elts      = _n_elts    [i_part];
    *connec      = _connec    [i_part];
    *global_num  = _global_num[i_part];
    *order       = _order;
    *ho_ordering = _ho_ordering;
  }


  void
  BlockHO::geomFinalize()
  {
    _reorder();
  }


  void
  BlockHO::_reorder()
  {
    if (PDM_ho_ordering_id_get(_ho_ordering) < 0) {
      PDM_error(__FILE__, __LINE__, 0,
                "CWP_Mesh_interf_ho_ordering_from_IJK_set must be called to provide CWIPI with the user's high-order numbering\n");
    }

    PDM_Mesh_nodal_elt_t elt_type = static_cast<PDM_Mesh_nodal_elt_t>(_block_type);

    for (int ipart = 0; ipart < _n_part; ipart++) {
      PDM_Mesh_nodal_reorder_elt_vtx(elt_type,
                                     _order,
                                     _ho_ordering,
                                     nullptr,
                                     _n_elts    [ipart],
                                     _connec    [ipart],
                                     _connec_ijk[ipart]);
    }
  }

}

/**
 * \endcond
 */