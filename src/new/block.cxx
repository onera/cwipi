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

#include "block.hxx"

#include "cwp.h"

#include "mesh.hxx"

#include "pdm_part_mesh_nodal.h"

#include <vector>

/**
 * \cond
 */

namespace cwipi {

  Block::Block()
    : _global_num_owner(PDM_OWNERSHIP_USER)
  {
  }


  Block::~Block()
  {
    for (int i = 0; i < _n_part; i++) {
      if (_global_num_owner == PDM_OWNERSHIP_KEEP && _global_num[i] != nullptr) {
        free(_global_num[i]);
      }
    }
  }


  void
  Block::blockAdd
  (
    CWP_Block_t  block_type,
    Mesh        *mesh
  )
  {
    _mesh = mesh;

    _block_type = block_type;
    _n_part     = static_cast<Mesh*>(mesh)->nPartGet();

    _n_elts     .resize(_n_part);
    _elt_centers.resize(_n_part, nullptr);
    _global_num .resize(_n_part, nullptr);
    _parent_num .resize(_n_part, nullptr);
  }


  double*
  Block::elementCentersGet
  (
    int i_part
  )
  {
    PDM_part_mesh_nodal_section_elt_center_reset(_mesh->pdmPartMeshNodalGet(),
                                                 _block_id_pdm,
                                                 i_part);

    PDM_part_mesh_nodal_section_elt_center_compute(_mesh->pdmPartMeshNodalGet(),
                                                   _block_id_pdm,
                                                   i_part,
                                                   PDM_OWNERSHIP_KEEP);

    _elt_centers[i_part] = const_cast<double*>(PDM_part_mesh_nodal_section_elt_center_get(_mesh->pdmPartMeshNodalGet(),
                                                                                          _block_id_pdm,
                                                                                          i_part,
                                                                                          PDM_OWNERSHIP_KEEP));

    return _elt_centers[i_part];
  }


  void
  Block::globalNumSet
  (
    int              i_part,
    CWP_g_num_t     *global_num,
    PDM_ownership_t  owner
  )
  {
    _global_num[i_part] = global_num;
    _global_num_owner   = owner;
  }


  void
  Block::globalNumFree
  (
    int i_part
  )
  {
    if (_global_num_owner == PDM_OWNERSHIP_KEEP) {
      if (_global_num[i_part] != nullptr) {
        free(_global_num[i_part]);
        _global_num[i_part] = nullptr;
      }
    }
  }


  void
  Block::parentNumSet
  (
    int              i_part,
    int             *parent_num,
    PDM_ownership_t  owner
  )
  {
    _parent_num[i_part] = parent_num;
    _parent_num_owner   = owner;
  }


  void
  Block::parentNumFree
  (
    int i_part
  )
  {
    if (_parent_num_owner == PDM_OWNERSHIP_KEEP) {
      if (_parent_num[i_part] != nullptr) {
        free(_parent_num[i_part]);
        _parent_num[i_part] = nullptr;
      }
    }
  }


}

/**
 * \endcond
 */