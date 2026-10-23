#ifndef __GLOBALDATA_H__
#define __GLOBALDATA_H__
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

#include <mpi.h>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class GlobalData globalData.hxx "globalData.hxx"
   * \brief Abstract global data.
   *
   * This class is global data abstract interface
   *
   */
  class GlobalData {

  public:

    /**
     * \brief Constructor.
     *
     * \param [in]  s_entity Size of datatype in bytes
     * \param [in]  stride   Stride of data
     * \param [in]  n_entity Number of entities
     *
     */
    GlobalData
    (
      const size_t s_entity,
      const int    stride,
      const int    n_entity
    );


    /**
     * \brief Destructor.
     *
     */
    virtual ~GlobalData();


    /**
     * \brief Get size of datatype in bytes.
     *
     * \return Reference to size of datatype in bytes
     */
    inline const size_t&
    sEntityGet() const;


    /**
     * \brief Get stride.
     *
     * \return Reference to stride
     */
    inline const int&
    strideGet() const;


    /**
     * \brief Get number of entities.
     *
     * \return Reference to number of entities
     */
    inline const int&
    nEntityGet() const;


    /**
     * \brief Set send data pointer.
     *
     * \param [in]  send_data Pointer to send data
     */
    inline void
    sendDataSet
    (
      void *send_data
    );


    /**
     * \brief Get send data pointer.
     *
     * \return Pointer to send data
     */
    inline void*
    sendDataGet() const;


    /**
     * \brief Set receive data pointer.
     *
     * \param [in]  recv_data Pointer to receive data
     */
    inline void
    recvDataSet
    (
      void *recv_data
    );


    /**
     * \brief Get receive data pointer.
     *
     * \return Pointer to receive data
     */
    inline void*
    recvDataGet() const;


    /**
     * \brief Set send MPI request.
     *
     * \param [in]  req MPI request for send
     */
    inline void
    dataSendRequestSet
    (
      MPI_Request req
    );


    /**
     * \brief Get send MPI request.
     *
     * \return MPI request for send
     */
    inline MPI_Request
    dataSendRequestGet() const;


    /**
     * \brief Set receive MPI request.
     *
     * \param [in]  req MPI request for receive
     */
    inline void
    dataRecvRequestSet
    (
      const MPI_Request req
    );


    /**
     * \brief Get receive MPI request.
     *
     * \return MPI request for receive
     */
    inline MPI_Request
    dataRecvRequestGet() const;

  private:

    size_t       _s_entity;          /*!< Size of datatype in bytes */
    int          _stride;            /*!< Stride of data */
    int          _n_entity;          /*!< Number of entities */
    void        *_send_data;         /*!< Send data buffer */
    void        *_recv_data;         /*!< Receive data buffer */
    MPI_Request  _data_send_request; /*!< Send MPI request */
    MPI_Request  _data_recv_request; /*!< Receive MPI request */

  };

  inline const size_t&
  GlobalData::sEntityGet() const
  {
    return _s_entity;
  }


  inline const int&
  GlobalData::strideGet() const
  {
    return _stride;
  }


  inline const int&
  GlobalData::nEntityGet() const
  {
    return _n_entity;
  }


  inline void
  GlobalData::sendDataSet
  (
    void *send_data
  )
  {
    _send_data = send_data;
  }


  inline void*
  GlobalData::sendDataGet() const
  {
    return _send_data;
  }


  inline void
  GlobalData::recvDataSet
  (
    void *recv_data
  )
  {
    _recv_data = recv_data;
  }


  inline void*
  GlobalData::recvDataGet() const
  {
    return _recv_data;
  }


  inline void
  GlobalData::dataSendRequestSet
  (
    const MPI_Request req
  )
  {
    _data_send_request = req;
  }


  inline MPI_Request
  GlobalData::dataSendRequestGet() const
  {
    return _data_send_request;
  }


  inline void
  GlobalData::dataRecvRequestSet
  (
    MPI_Request req
  )
  {
    _data_recv_request = req;
  }


  inline MPI_Request
  GlobalData::dataRecvRequestGet() const
  {
    return _data_recv_request;
  }

}

/**
 * \endcond
 */

#endif // __GLOBALDATA_H__