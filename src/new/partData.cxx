/*
  This file is part of the CWIPI library.

  Copyright (C) 2022-2023  ONERA

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
#include <sstream>

#include <partData.hxx>
#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_part_to_part.h"
#include "pdm_logging.h"

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  PartData::PartData(std::string           part_data_id,
                     CWP_PartData_exch_t   exch_type,
                     CWP_g_num_t         **gnum_elt,
                     int                  *n_elt,
                     int                   n_part):
  _part_data_id(part_data_id),
  _part1_to_part2_data(NULL),
  _part1_to_part2_idx(NULL),
  _send_request(-2),
  _n_send_calls(0),
  _recv_request(-2),
  _n_recv_calls(0),
  _recv_buffer(NULL),
  _filtered_gnum1_come_from(NULL)
  {
    if (exch_type == CWP_PARTDATA_SEND) {
      _gnum_elt1 = gnum_elt;
      _n_elt1    = n_elt;
      _n_part1   = n_part;
    }

    else if (exch_type == CWP_PARTDATA_RECV) {
      _gnum_elt2 = gnum_elt;
      _n_elt2    = n_elt;
      _n_part2   = n_part;
    }
  }

  PartData::~PartData()
  {
    log_trace("~PartData %s\n", _part_data_id.c_str());
  }

  uint32_t PartData::_adler32
  (
   const void *buf,
   size_t buflength
   )
  {

    const uint8_t * buffer = (const uint8_t *)buf;

    uint32_t s1 = 1;
    uint32_t s2 = 0;

    for (size_t n = 0; n < buflength; n++) {
      s1 = (s1 + buffer[n]) % 65521;
      s2 = (s2 + s1) % 65521;
    }

    return (s2 << 16) | s1;
  }

}

/**
 * \endcond
 */
