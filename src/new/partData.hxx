#ifndef __PARTDATA_H__
#define __PARTDATA_H__
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

#include "cwp.h"

#include "pdm_part_to_part.h"

#include <map>
#include <string>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class PartData partData.hxx "partData.hxx"
   * \brief Abstract partitioned data.
   *
   * This class is the partitioned data abstract interface
   *
   */
  class PartData {

  public:

    /**
     * \brief Constructor.
     *
     * \param [in]  part_data_id Part data identifier
     * \param [in]  exch_type    Exchange type
     * \param [in]  global_num   Global numbering
     * \param [in]  n_elts       Number of elements
     * \param [in]  n_part       Number of partitions
     *
     * The global_num and n_elts pointers are borrowed (not owned) by the PartData
     */
    PartData
    (
      std::string           part_data_id,
      CWP_PartData_exch_t   exch_type,
      CWP_g_num_t         **global_num,
      int                  *n_elts,
      int                   n_part
    );


    /**
     * \brief Destructor.
     *
     */
    ~PartData();


    /**
     * \brief Set the data.
     *
     * \param [in]  exch_id Exchange identifier
     * \param [in]  s_unit  Datatype size in bytes
     * \param [in]  data    Data pointer
     */
    void
    dataSet
    (
      int      exch_id,
      size_t   s_unit,
      void   **data
    );


    /**
     * \brief Set the request.
     *
     * \param [in]  exch_id Exchange identifier
     * \param [in]  request Request
     */
    void
    requestSet
    (
      int exch_id,
      int request
    );


    /**
     * \brief Get the request.
     *
     * \param [in]  exch_id Exchange identifier
     * \return Request
     */
    int
    requestGet
    (
      int exch_id
    );


    /**
     * \brief Filter the received data in-place.
     *
     * \param [in]  exch_id Exchange identifier
     */
    void
    recvDataFilter
    (
      int exch_id
    );


    /**
     * \brief Clear an exchange.
     *
     * \param [in]  exch_id Exchange identifier
     */
    void
    exchClear
    (
      int exch_id
    );


    /**
     * \brief Get the part-to-part object.
     *
     * \return Part-to-part object
     */
    inline PDM_part_to_part_t*
    pdmPartToPartGet() const;


    /**
     * \brief Set the part-to-part object.
     *
     * \param [in]  ptp Part-to-part object
     */
    inline void
    pdmPartToPartSet
    (
      PDM_part_to_part_t *ptp
    );


    /**
     * \brief Get the global numbering.
     *
     * \return Global numbering
     */
    inline CWP_g_num_t**
    globalNumGet() const;


    /**
     * \brief Get the number of elements.
     *
     * \return Number of elements
     */
    inline int*
    nEltsGet() const;


    /**
     * \brief Get the number of partitions.
     *
     * \return Number of partitions
     */
    inline int
    nPartGet() const;


    /**
     * \brief Get the exchange type.
     *
     * \return Exchange type
     */
    inline CWP_PartData_exch_t
    exchTypeGet() const;

  private:

    std::string              _part_data_id; /*!< Part data identifier */
    CWP_PartData_exch_t      _exch_type;    /*!< Exchange type */
    PDM_part_to_part_t      *_ptp;          /*!< Part-to-part object */

    CWP_g_num_t            **_global_num;   /*!< Global numbering */
    int                     *_n_elts;       /*!< Number of elements */
    int                      _n_part;       /*!< Number of partitions */
    std::map<int, void **>   _data;         /*!< Data per exchange */

    std::map<int, size_t>    _s_unit;       /*!< Datatype size in bytes per exchange */
    std::map<int, int>       _request;      /*!< Request per exchange */
  };

  inline PDM_part_to_part_t*
  PartData::pdmPartToPartGet() const
  {
    return _ptp;
  }


  inline void
  PartData::pdmPartToPartSet
  (
    PDM_part_to_part_t *ptp
  )
  {
    _ptp = ptp;
  }


  inline CWP_g_num_t**
  PartData::globalNumGet() const
  {
    return _global_num;
  }


  inline int*
  PartData::nEltsGet() const
  {
    return _n_elts;
  }


  inline int
  PartData::nPartGet() const
  {
    return _n_part;
  }


  inline CWP_PartData_exch_t
  PartData::exchTypeGet() const
  {
    return _exch_type;
  }

}

/**
 * \endcond
 */

#endif // __PARTDATA_H__