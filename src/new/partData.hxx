#ifndef __PARTDATA_H__
#define __PARTDATA_H__
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
#include <iostream>

#include "cwp.h"
#include "cwp_priv.h"

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class PartData partData.hxx "partData.hxx"
   * \brief Abstract partionned data
   *
   *  This class is partionned data abstract interface
   *
   */


  class PartData {

  public:

    /**
     * \brief Constructor
     *
     */
    PartData() {}

    PartData(std::string     global_data_id,
               size_t          s_send_entity,
               int             send_stride,
               int             n_send_entity,
               void           *send_data);

    /**
     * \brief Destructor
     *
     */
    ~PartData();

    // TO DO: setters, getters etc

  private:

    std::string     _part_data_id;

    // TO DO: variables

  };

}

/**
 * \endcond
 */

#endif //__PARTDATA_H__
