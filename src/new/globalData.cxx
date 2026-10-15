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

#include "globalData.hxx"

#include "cwp.h"

#include <mpi.h>

/**
 * \cond
 */

namespace cwipi {

  GlobalData::GlobalData
  (
    size_t s_entity,
    int    stride,
    int    n_entity
  )
  : _s_entity(s_entity),
    _stride(stride),
    _n_entity(n_entity),
    _data_send_request(MPI_REQUEST_NULL),
    _data_recv_request(MPI_REQUEST_NULL),
    _send_data(nullptr),
    _recv_data(nullptr)
  {
  }


  GlobalData::~GlobalData()
  {
  }

}

/**
 * \endcond
 */