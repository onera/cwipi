#ifndef __GLOBALDATA_H__
#define __GLOBALDATA_H__
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

#include "cwp.h"
#include "cwp_priv.h"

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class GlobalData globalData.hxx "globalData.hxx"
   * \brief Abstract global data
   *
   *  This class is global data abstract interface
   *
   */


  class GlobalData {

  public:

    /**
     * \brief Constructor
     *
     */
    GlobalData() {}

    // send

    GlobalData(std::string     global_data_id,
               size_t          s_entity,
               int             stride,
               int             n_entity,
               void           *data);

    // recv

    GlobalData(std::string     global_data_id,
               size_t         *s_recv_entity,
               int            *recv_stride,
               int            *n_recv_entity,
               void          **recv_data);

    /**
     * \brief Destructor
     *
     */
    ~GlobalData();

    // Getters

    inline MPI_Request *
    global_request_get()
    {
      return _global_request;
    }

    inline MPI_Request *
    data_request_get()
    {
      return _data_request;
    }

    // send

    inline size_t
    s_send_entity_get()
    {
      return _s_send_entity;
    }

    inline int
    send_stride_get()
    {
      return _send_stride;
    }

    inline int
    n_send_entity_get()
    {
      return _n_send_entity;
    }

    inline void *
    send_data_get()
    {
      return _send_data;
    }

    // recv

    inline size_t *
    s_recv_entity_get()
    {
      return _s_recv_entity;
    }

    inline int *
    recv_stride_get()
    {
      return _recv_stride;
    }

    inline int *
    n_recv_entity_get()
    {
      return _n_recv_entity;
    }

    inline void **
    recv_data_get()
    {
      return _recv_data;
    }

  private:

    std::string     _global_data_id;
    MPI_Request    *_global_request;
    MPI_Request    *_data_request;

    // send
    size_t          _s_send_entity;
    int             _send_stride;
    int             _n_send_entity;
    void           *_send_data;

    // recv
    size_t         *_s_recv_entity;
    int            *_recv_stride;
    int            *_n_recv_entity;
    void          **_recv_data;

  };

}

/**
 * \endcond
 */

#endif //__GLOBALDATA_H__
