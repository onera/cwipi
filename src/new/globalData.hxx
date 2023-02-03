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
#include <iostream>

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
               size_t          s_send_entity,
               int             send_stride,
               int             n_send_entity,
               void           *send_data);

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

    // send

    inline size_t&
    s_send_entity_get()
    {
      return _s_send_entity;
    }

    inline int&
    send_stride_get()
    {
      return _send_stride;
    }

    inline int&
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

    // request

    inline MPI_Request
    s_entity_request_get()
    {
      return _s_entity_request;
    }

    inline MPI_Request
    stride_request_get()
    {
      return _stride_request;
    }

    inline MPI_Request
    n_entity_request_get()
    {
      return _n_entity_request;
    }

    inline MPI_Request
    data_request_get()
    {
      return _data_request;
    }

    inline void
    s_entity_request_set
    (
     MPI_Request req
    )
    {
      _s_entity_request = req;
    }

    inline void
    stride_request_set
    (
     MPI_Request req
    )
    {
      _stride_request = req;
    }

    inline void
    n_entity_request_set
    (
     MPI_Request req
    )
    {
      _n_entity_request = req;
    }

    inline void
    data_request_set
    (
     MPI_Request req
    )
    {
      _data_request = req;
    }

  private:

    std::string     _global_data_id;
    MPI_Request     _s_entity_request;
    MPI_Request     _stride_request;
    MPI_Request     _n_entity_request;
    MPI_Request     _data_request;

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
