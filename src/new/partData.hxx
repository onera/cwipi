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

#include "pdm_part_to_part.h"

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

    PartData(std::string           part_data_id,
             CWP_PartData_exch_t   exch_type,
             CWP_g_num_t         **gnum_elt,
             int                  *n_elt,
             int                   n_part);

    /**
     * \brief Destructor
     *
     */
    ~PartData();

    // Getters

    inline PDM_part_to_part_t *
    get_ptp()
    {
      return _ptp;
    }

    inline size_t
    get_s_data()
    {
      return _s_data;
    }

    inline int
    get_n_components()
    {
      return _n_components;
    }

    inline CWP_g_num_t  **
    get_gnum_elt1()
    {
      return _gnum_elt1;
    }

    inline int *
    get_n_elt1()
    {
      return _n_elt1;
    }

    inline int
    get_n_part1()
    {
      return _n_part1;
    }

    inline void **
    get_part1_to_part2_data()
    {
      return _part1_to_part2_data;
    }

    inline int **
    get_part1_to_part2_idx()
    {
      return _part1_to_part2_idx;
    }

    inline int *
    get_request1()
    {
      return _request1;
    }

    inline CWP_g_num_t  **
    get_gnum_elt2()
    {
      return _gnum_elt2;
    }

    inline int *
    get_n_elt2()
    {
      return _n_elt2;
    }

    inline int
    get_n_part2()
    {
      return _n_part2;
    }

    inline void ***
    get_part2_data()
    {
      return _part2_data;
    }

    inline int *
    get_request2()
    {
      return _request2;
    }

    inline CWP_g_num_t **
    get_filtered_gnum1_come_from()
    {
      return _filtered_gnum1_come_from;
    }

    inline void **
    get_recv_buffer()
    {
      return _recv_buffer;
    }

    // Setters

    inline void
    set_ptp
    (
     PDM_part_to_part_t *ptp
    )
    {
      _ptp = ptp;
    }

    inline void
    set_s_data
    (
     size_t s_data
    )
    {
      _s_data = s_data;
    }

    inline void
    set_n_components
    (
     int n_components
    )
    {
      _n_components = n_components;
    }

    inline void
    set_part1_to_part2_data
    (
     void ** part1_to_part2_data
    )
    {
      _part1_to_part2_data = part1_to_part2_data;
    }

    inline void
    set_part1_to_part2_idx
    (
     int ** part1_to_part2_idx
    )
    {
      _part1_to_part2_idx = part1_to_part2_idx;
    }

    inline void
    set_request1
    (
     int * request1
    )
    {
      _request1 = request1;
    }

    inline void
    set_part2_data
    (
     void *** part2_data
    )
    {
      _part2_data = part2_data;
    }

    inline void
    set_request2
    (
     int * request2
    )
    {
      _request2 = request2;
    }

    inline void
    set_filtered_gnum1_come_from
    (
     CWP_g_num_t ** filtered_gnum1_come_from
    )
    {
      _filtered_gnum1_come_from = filtered_gnum1_come_from;
    }

    inline void
    set_recv_buffer
    (
     void ** recv_buffer
    )
    {
      _recv_buffer = recv_buffer;
    }

    int
    get_tag
    (
     const std::string    global_data_id,
     MPI_Comm        comm
    );

  private:

    uint32_t _adler32
    (
     const void *buf,
     size_t buflength
    );

    // global
    std::string         _part_data_id;
    PDM_part_to_part_t *_ptp;
    size_t              _s_data;
    int                 _n_components;

    // send
    CWP_g_num_t  **_gnum_elt1;
    int           *_n_elt1;
    int            _n_part1;
    void         **_part1_to_part2_data;
    int          **_part1_to_part2_idx;
    int           *_request1;

    // recv
    CWP_g_num_t  **_gnum_elt2;
    int           *_n_elt2;
    int            _n_part2;
    void        ***_part2_data;
    int           *_request2;

    // intern
    void         **_recv_buffer;
    CWP_g_num_t  **_filtered_gnum1_come_from;

  };

}

/**
 * \endcond
 */

#endif //__PARTDATA_H__
