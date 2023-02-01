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

#include <sstream>
#include <mesh.hxx>
#include <map>

#include <coupling.hxx>
#include "cwp.h"
#include "cwp_priv.h"
#include "pdm_writer.h"

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  GlobalData::GlobalData(size_t          s_send_entity,
                         int             send_stride,
                         int             n_send_entity,
                         void           *send_data):
  _s_send_entity(s_send_entity),
  _send_stride(send_stride),
  _n_send_entity(n_send_entity),
  _send_data(send_data):
  {
  }

  GlobalData::~GlobalData()
  {
    free(send_data);
  }

}

/**
 * \endcond
 */
