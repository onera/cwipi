#ifndef __PARTDATA2_H__
#define __PARTDATA2_H__

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
#include <iostream>

#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_part_to_part.h"

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class PartData2 partData2.hxx "partData2.hxx"
   * \brief Abstract partionned data
   *
   *  This class is partionned data abstract interface
   *
   */


  class PartData2 {

  public:

    /**
     * \brief Constructors
     *
     */
    // PartData2() {}

    PartData2(std::string           part_data_id,
              CWP_PartData_exch_t   exch_type,
              CWP_g_num_t         **gnum_elt,
              int                  *n_elt,
              int                   n_part);

    /**
     * \brief Destructor
     *
     */
    ~PartData2();

    void
    data_set(int                   request,
             CWP_PartData_exch_t   exch_type,
             void                **data);

    bool
    data_get(int                    request,
             CWP_PartData_exch_t    exch_type,
             void                ***data);

    void
    recv_data_filter(int request);

    void
    request_clear(int                 request,
                  CWP_PartData_exch_t exch_type);

    /* awkward... (why does the destructor get called early??) */
    inline void s_unit_delete()
    {
      delete &_s_unit;
    }


    inline PDM_part_to_part_t *
    ptp_get()
    {
      return _ptp;
    }

    inline void
    ptp_set(PDM_part_to_part_t *ptp)
    {
      _ptp = ptp;
    }

    inline CWP_g_num_t **
    gnum_elt_get(CWP_PartData_exch_t exch_type)
    {
      return _gnum_elt[exch_type];
    }

    inline int *
    n_elt_get(CWP_PartData_exch_t exch_type)
    {
      return _n_elt[exch_type];
    }

    inline int
    n_part_get(CWP_PartData_exch_t exch_type)
    {
      return _n_part[exch_type];
    }

    inline int **
    part1_to_part2_idx_get()
    {
      int          *n_elt1             = NULL;
      int         **part1_to_part2_idx = NULL;
      PDM_g_num_t **part1_to_part2     = NULL;
      PDM_part_to_part_part1_to_part2_get(_ptp,
                                          &n_elt1,
                                          &part1_to_part2_idx,
                                          &part1_to_part2);

      return part1_to_part2_idx;
    }

    // inline int
    // n_exch_get(CWP_PartData_exch_t exch_type) {
    //   return _n_exch[exch_type];
    // }

    // inline void
    // n_exch_increment(CWP_PartData_exch_t exch_type) {
    //   _n_exch[exch_type]++;
    // }

    // inline void
    // part1_to_part2_idx_set(part1_to_part2_idx)
    // {
    //   _part1_to_part2_idx = part1_to_part2_idx;
    // }


  private:

    std::string              _part_data_id;
    PDM_part_to_part_t      *_ptp;

    CWP_g_num_t            **_gnum_elt[2]; // redundant with ptp -> to remove
    int                     *_n_elt[2];    // redundant with ptp -> to remove
    int                      _n_part[2];   // redundant with ptp -> to remove
    std::map<int, void **>   _data[2];
    // std::vector<int>         _request[2];

    /* Internal */
    std::map<int, int>      &_s_unit;
    // int _n_exch[2];
    // int **_part1_to_part2_idx;
  };

}

/**
 * \endcond
 */

#endif //__PARTDATA2_H__
