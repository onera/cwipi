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

#include "partData.hxx"

#include "cwp.h"

#include "pdm_error.h"

#include <cassert>
#include <map>

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  PartData::PartData
  (
    string                part_data_id,
    CWP_PartData_exch_t   exch_type,
    CWP_g_num_t         **global_num,
    int                  *n_elts,
    int                   n_part
  )
  : _part_data_id(part_data_id),
    _exch_type(exch_type),
    _ptp(nullptr),
    _global_num(global_num),
    _n_elts(n_elts),
    _n_part(n_part)
  {
  }


  PartData::~PartData()
  {
  }


  void
  PartData::dataSet
  (
    int      exch_id,
    size_t   s_unit,
    void   **data
  )
  {
    map<int, void **>::iterator it = _data.find(exch_id);

    if (it != _data.end()) {
      PDM_error("PartData '%s': Data pointer with exch_type %d already set for exch_id %d",
                _part_data_id.c_str(), _exch_type, exch_id);
    }

    pair<int, void **> new_pair_data(exch_id, data);
    _data.insert(new_pair_data);
    pair<int, size_t> new_pair_s_unit(exch_id, s_unit);
    _s_unit.insert(new_pair_s_unit);
  }


  void
  PartData::requestSet
  (
    int exch_id,
    int request
  )
  {
    map<int, int>::iterator it = _request.find(exch_id);

    if (it != _request.end()) {
      PDM_error("PartData '%s': Request pointer with exch_type %d already set for exch_id %d",
                _part_data_id.c_str(), _exch_type, exch_id);
    }

    pair<int, int> new_pair(exch_id, request);
    _request.insert(new_pair);
  }


  int
  PartData::requestGet
  (
    int exch_id
  )
  {
    map<int, int>::iterator it = _request.find(exch_id);

    if (it == _request.end()) {
      PDM_error("PartData '%s': Request for exch_type %d undefined for exch_id %d",
                _part_data_id.c_str(), _exch_type, exch_id);
    }

    return it->second;
  }


  void
  PartData::recvDataFilter
  (
    int exch_id
  )
  {
    assert(_exch_type == CWP_PARTDATA_RECV);

    map<int, void **>::iterator it = _data.find(exch_id);

    if (it == _data.end()) {
      PDM_error("PartData '%s': Recv data pointer not set for exch_id %d",
                _part_data_id.c_str(), exch_id);
    }

    map<int, size_t>::iterator it_s_unit = _s_unit.find(exch_id);
    size_t s_unit = it_s_unit->second;

    int         **come_from_idx = nullptr;
    PDM_g_num_t **come_from     = nullptr;
    PDM_part_to_part_gnum1_come_from_get(_ptp,
                                         &come_from_idx,
                                         &come_from);

    int n_part1 = 0;
    int n_part2 = 0;
    PDM_part_to_part_n_part_get(_ptp,
                                &n_part1,
                                &n_part2);

    unsigned char **data = reinterpret_cast<unsigned char **>(it->second);

    for (int i_part = 0; i_part < n_part2; i_part++) {
      int  n_ref = 0;
      int *ref   = nullptr;
      PDM_part_to_part_ref_lnum2_single_part_get(_ptp,
                                                 i_part,
                                                 &n_ref,
                                                 &ref);
      for (int i = 0; i < n_ref; i++) {
        int idx = come_from_idx[i_part][i];
        for (size_t j = 0; j < s_unit; j++) {
          data[i_part][s_unit*i + j] = data[i_part][s_unit*idx + j];
        }
      }
    }
  }


  void
  PartData::exchClear
  (
    int exch_id
  )
  {
    _data   .erase(exch_id);
    _s_unit .erase(exch_id);
    _request.erase(exch_id);
  }

}

/**
 * \endcond
 */