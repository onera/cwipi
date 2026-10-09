#ifndef __COMMUNICATION_H__
#define __COMMUNICATION_H__
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

#include "codeProperties.hxx"
#include "couplingDB.hxx"

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <string>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class Communication communication.hxx "communication.hxx"
   * \brief Communication abstract interface.
   *
   * This class is the abstract interface of communication instances
   * and manages MPI communication between coupled codes
   *
   */
  class Communication {

  public:

    /**
     * \brief Constructor.
     *
     */
    Communication();


    /**
     * \brief Destructor.
     *
     */
    virtual ~Communication();


    /**
     * \brief Initialize coupling communicators.
     *
     * \param [in]  localCodeProperties Local code properties
     * \param [in]  cplCodeProperties   Coupled code properties
     * \param [in]  cplId               Coupling identifier
     * \param [in]  cplDB               Coupling database
     */
    void
    init
    (
      const CodeProperties &localCodeProperties,
      const CodeProperties &cplCodeProperties,
      const std::string    &cplId,
            CouplingDB     &cplDB
    );


    /**
     * \brief Initialize coupling communicators.
     *
     * \param [in]  cplCodeComm Coupled code communication
     */
    void
    init
    (
      Communication &cplCodeComm
    );


    /**
     * \brief Get a tag from an Adler code.
     *
     * \param [in]  name   Data identifier
     * \param [in]  offset Offset
     * \return Tag
     */
    uint32_t
    tagGet
    (
      const std::string name,
            int         offset
    );


    /**
     * \brief Exchange data between two coupled codes through the union communicator.
     *        All code ranks receive data (sendrecv between root ranks then bcast).
     *
     * \param [in]  s_data          Size of data values to exchange
     * \param [in]  n_send_data     Number of data to send
     * \param [in]  send_data       Array of data to send
     * \param [in]  n_send_data_cpl Number of data to send (if current rank is shared between codes)
     * \param [in]  send_data_cpl   Array of data to send (if current rank is shared between codes)
     * \param [in]  n_recv_data     Number of data to receive
     * \param [out] recv_data       Array of data to receive
     * \param [in]  n_recv_data_cpl Number of data to receive (if current rank is shared between codes)
     * \param [out] recv_data_cpl   Array of data to receive (if current rank is shared between codes)
     */
    void
    globalDataBetweenCodesThroughUnionCommExch
    (
      size_t  s_data,
      int     n_send_data,
      void   *send_data,
      int     n_send_data_cpl,
      void   *send_data_cpl,
      int     n_recv_data,
      void   *recv_data,
      int     n_recv_data_cpl,
      void   *recv_data_cpl
    );


    /**
     * \brief Non blocking send of global data array.
     *
     * \param [in]  global_data_id    Global data identifier
     * \param [out] data_send_request Send request
     * \param [in]  s_entity          Size of data values
     * \param [in]  stride            Stride
     * \param [in]  n_entity          Number of entities
     * \param [in]  send_data         Array of data to send
     */
    void
    globalDataBetweenCodesThroughUnionCommIssend
    (
      const std::string  global_data_id,
            MPI_Request *data_send_request,
            size_t       s_entity,
            int          stride,
            int          n_entity,
            void        *send_data
    );


    /**
     * \brief Non blocking receive of global data array.
     *
     * \param [in]  global_data_id    Global data identifier
     * \param [out] data_recv_request Receive request
     * \param [in]  s_entity          Size of data values
     * \param [in]  stride            Stride
     * \param [in]  n_entity          Number of entities
     * \param [out] recv_data         Array of data to receive
     */
    void
    globalDataBetweenCodesThroughUnionCommIrecv
    (
      const std::string  global_data_id,
            MPI_Request *data_recv_request,
            size_t       s_entity,
            int          stride,
            int          n_entity,
            void        *recv_data
    );


    /**
     * \brief Non blocking send wait of global data array.
     *
     * \param [in]  data_send_request Send request
     * \param [in]  s_entity          Size of data values
     * \param [in]  stride            Stride
     * \param [in]  n_entity          Number of entities
     * \param [in]  send_data         Array of data to send
     * \param [out] recv_data         Array of data to receive
     */
    void
    globalDataBetweenCodesThroughUnionCommWaitIssend
    (
      MPI_Request  data_send_request,
      size_t       s_entity,
      int          stride,
      int          n_entity,
      void        *send_data,
      void        *recv_data
    );


    /**
     * \brief Non blocking receive wait of global data array.
     *
     * \param [in]  data_recv_request Receive request
     * \param [in]  s_entity          Size of data values
     * \param [in]  stride            Stride
     * \param [in]  n_entity          Number of entities
     * \param [in]  send_data         Array of data to send
     * \param [out] recv_data         Array of data to receive
     */
    void
    globalDataBetweenCodesThroughUnionCommWaitIrecv
    (
      MPI_Request  data_recv_request,
      size_t       s_entity,
      int          stride,
      int          n_entity,
      void        *send_data,
      void        *recv_data
    );


    /**
     * \brief Get the union communicator.
     *
     * \return Union communicator
     */
    inline MPI_Comm
    unionCommGet() const;


    /**
     * \brief Get the coupling communicator.
     *
     * \return Coupling communicator
     */
    inline MPI_Comm
    cplCommGet() const;


    /**
     * \brief Get the root rank of the coupled code into the union communicator.
     *
     * \return Root rank of the coupled code into the union communicator
     */
    inline int
    unionCommCplCodeRootRanksGet() const;


    /**
     * \brief Get the root rank of the local code into the union communicator.
     *
     * \return Root rank of the local code into the union communicator
     */
    inline int
    unionCommLocCodeRootRanksGet() const;


    /**
     * \brief Return the communicator type.
     *
     * \return Communicator type
     */
    virtual CWP_Comm_t
    commTypeGet() = 0;

  private:

    /**
     * \brief Compute an Adler code.
     *
     * \param [in]  buf       Buffer
     * \param [in]  buflength Buffer length
     * \return Adler code
     */
    uint32_t
    _adler32Compute
    (
      const void   *buf,
            size_t  buflength
    );

  protected:

    /**
     * \brief Build the coupling communicator from \ref CWP_Comm_t.
     *
     * \param [in]  cplCodeCommType Coupling communicator type
     */
    virtual void
    cplCommCreate
    (
      CWP_Comm_t cplCodeCommType
    ) = 0;

    const CodeProperties *_local_code_properties;          /*!< Pointer to the local code properties */
    const CodeProperties *_cpl_code_properties;            /*!< Pointer to the coupled code properties */

    int                    _tag;                           /*!< Tag for MPI */
    MPI_Comm               _union_comm;                    /*!< Union communicator between coupled codes */
    MPI_Comm               _cpl_comm;                      /*!< Coupling communicator */
    int                    _loc_code_root_rank_union_comm; /*!< Root rank of the local code into the union communicator */
    int                    _cpl_code_root_rank_union_comm; /*!< Root rank of the coupled code into the union communicator */

  };

  inline MPI_Comm
  Communication::unionCommGet() const
  {
    return _union_comm;
  }


  inline MPI_Comm
  Communication::cplCommGet() const
  {
    return _cpl_comm;
  }


  inline int
  Communication::unionCommCplCodeRootRanksGet() const
  {
    return _cpl_code_root_rank_union_comm;
  }


  inline int
  Communication::unionCommLocCodeRootRanksGet() const
  {
    return _loc_code_root_rank_union_comm;
  }

}

/**
 * \endcond
 */

#endif // __COMMUNICATION_H__