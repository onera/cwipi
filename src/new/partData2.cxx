/*
  This file is part of the CWIPI library.

  Copyright (C) 2023  ONERA

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

#include <map>
#include <vector>
#include <sstream>

#include <partData2.hxx>
#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_part_to_part.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_error.h"

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  /**
    * \brief Constructor
    *
    */
  PartData2::PartData2(std::string           part_data_id,
                       CWP_PartData_exch_t   exch_type,
                       CWP_g_num_t         **gnum_elt,
                       int                  *n_elt,
                       int                   n_part):
  _part_data_id(part_data_id),
  _ptp(NULL),
  _s_unit(*(new map < int, int >()))
  // _part1_to_part2_idx(NULL)
  {
    log_trace("PartData2 %s, exch_type %d, n_part %d\n", part_data_id.c_str(), exch_type, n_part);
    _gnum_elt[exch_type] = gnum_elt;
    _n_elt   [exch_type] = n_elt;
    _n_part  [exch_type] = n_part;

    // _n_exch  [exch_type] = 0;

    // if (exch_type == CWP_PARTDATA_SEND) {
    //   _part1_to_part2_idx = malloc(sizeof(int *) * n_part);
    //   for (int i = 0; i < n_part; i++) {
    //     _part1_to_part2_idx[i] = PDM_array_new_idx_from_const_stride_int(1, n_elt[i]);
    //   }
    // }
    // else {
    //   _part1_to_part2_idx = NULL;
    // }

    // _data   [exch_type] = new std::map<int, void **>();
    // _request[exch_type] = new std::vector<int>();
  }


  /**
    * \brief Destructor
    *
    */
  PartData2::~PartData2()
  {
    // for (int i = 0; i < 2; i++) {
    //   delete &_data[i];
    //   delete &_request[i];
    // }
    // log_trace("~PartData2 %s\n", _part_data_id.c_str());

    // delete &_s_unit;
  }


  void
  PartData2::data_set(int                   request,
                      CWP_PartData_exch_t   exch_type,
                      void                **data)
  {
    map<int, void **>::iterator it = _data[exch_type].find(request);

    if (it != _data[exch_type].end()) {
      PDM_error(__FILE__, __LINE__, 0, "Data pointer with exch_type %d already set for request %d\n", exch_type, request);
    }

    pair<int, void **> new_pair(request, data);
    _data[exch_type].insert(new_pair);
  };

  bool
  PartData2::data_get(int                    request,
                      CWP_PartData_exch_t    exch_type,
                      void                ***data)
  {
    map<int, void **>::iterator it = _data[exch_type].find(request);

    if (it == _data[exch_type].end()) {
      return false;
    }
    else {
      *data = it->second;
      return true;
    }
  }


  void
  PartData2::recv_data_filter(int request)
  {
    map<int, void **>::iterator it = _data[CWP_PARTDATA_RECV].find(request);

    if (it == _data[CWP_PARTDATA_RECV].end()) {
      PDM_error(__FILE__, __LINE__, 0, "Recv data pointer with was not set for request %d\n", request);
    }

    map<int, int>::iterator it_s_unit = _s_unit.find(request);
    int s_unit = it_s_unit->second;

    int         **come_from_idx = NULL;
    PDM_g_num_t **come_from     = NULL;
    PDM_part_to_part_gnum1_come_from_get(_ptp,
                                         &come_from_idx,
                                         &come_from);

    int n_part1 = 0;
    int n_part2 = 0;
    PDM_part_to_part_n_part_get(_ptp,
                                &n_part1,
                                &n_part2);


    unsigned char **data = (unsigned char **) it->second;

    for (int ipart = 0; ipart < n_part2; ipart++) {
      int  n_ref = 0;
      int *ref   = NULL;
      PDM_part_to_part_ref_lnum2_single_part_get(_ptp,
                                                 ipart,
                                                 &n_ref,
                                                 &ref);
      for (int i = 0; i < n_ref; i++) {
        int idx = come_from_idx[ipart][i];
        for (int j = 0; j < s_unit; j++) {
          data[ipart][s_unit*i + j] = data[ipart][s_unit*idx + j];
        }
      }
    }
  };


  void
  PartData2::request_clear(int                 request,
                           CWP_PartData_exch_t exch_type)
  {
    _data[exch_type].erase(request);

    map<int, int>::iterator it = _s_unit.find(request);
    if (it != _s_unit.end()) {
      _s_unit.erase(request);
    }
  }


}

/**
 * \endcond
 */
