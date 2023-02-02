/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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

#include <globalData.hxx>
#include "cwp.h"
#include "cwp_priv.h"

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  GlobalData::GlobalData(std::string     global_data_id,
                         size_t          s_send_entity,
                         int             send_stride,
                         int             n_send_entity,
                         void           *send_data):
  _global_data_id(global_data_id),
  _s_send_entity(s_send_entity),
  _send_stride(send_stride),
  _n_send_entity(n_send_entity),
  _send_data(send_data)
  {
  }

  GlobalData::GlobalData(std::string     global_data_id,
                         size_t         *s_recv_entity,
                         int            *recv_stride,
                         int            *n_recv_entity,
                         void          **recv_data):
  _global_data_id(global_data_id),
  _s_recv_entity(s_recv_entity),
  _recv_stride(recv_stride),
  _n_recv_entity(n_recv_entity),
  _recv_data(recv_data)
  {
  }

  GlobalData::~GlobalData()
  {
    if (_send_data != NULL) free(_send_data);
    if (*_recv_data != NULL) free(*_recv_data);
  }

}

/**
 * \endcond
 */
