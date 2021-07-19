#ifndef __COMMUNICATION_H__
#define __COMMUNICATION_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2013-2017  ONERA

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

#include "codeProperties.hxx"
#include "couplingDB.hxx"


/**
 * \cond
 */

using namespace std;

namespace cwipi {

  /**
   * \class Communication communication.hxx "communication.hxx"
   * \brief Communication abstract interface
   *
   *  This class is the abstract interface of communication instances.
   *  A communication instance manage MPI communication between coupled codes
   */

  class Communication {

  public :

    /**
     *
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
     *
     * \brief Initialize coupling communicators.
     *
     * \param [in]  localCodeProperties   Local code properties
     * \param [in]  cplCodeProperties     Coupled code properties
     * \param [in]  cplId                 Coupling identifier
     * \param [in]  cplDB                 Coupling DataBase
     *
     */

    void
    init
    (
     const CodeProperties &localCodeProperties,
     const CodeProperties &cplCodeProperties,
     const string         &cplId,
     CouplingDB           &cplDB
     );

    /**
     *
     * \brief Initialize coupling communicators.
     *
     * \param [in]  cplCodeComm           Coupled code communication
     *
     */

    void
    init
    (
     Communication &cplCodeComm
    );


    /**
     *
     * \brief Exchange Data between two coupled code through Union Comminucator
     *        All code ranks receive data (sendRecv between root rank then Bcast)
     *
     * \param [in]    s_data           Size of a data tot exchange
     * \param [in]    n_send_data      Number of data to send
     * \param [in]    send_data        Array of data to send
     * \param [in]    n_recv_data      Number of data to receive
     * \param [inout] recv_data        Array of data to receive
     * \param [inout] request          MPI Request
     *
     */

    void
    iexchGlobalDataBetweenCodesThroughUnionCom
    (
     size_t       s_data,
     int          n_send_data,
     void        *send_data,
     int          n_recv_data,
     void        *recv_data,
     MPI_Request *request
    );



    MPI_Comm
    unionCommGet
    (
    );

    std::vector<int>*
    unionCommCplRanksGet
    (
    );

    std::vector<int>*
    unionCommLocRanksGet
    (
    );

    MPI_Comm
    cplCommGet
    (
    );

    std::vector<int>*
    cplCommCplRanksGet
    (
    );

    std::vector<int>*
    cplCommLocRanksGet
    (
    );


    int
    unionCommCplCodeRootRanksGet
    (
    );

    int
    unionCommLocCodeRootRanksGet
    (
    );

    int
    cplCommCplCodeRootRanksGet
    (
    );

    int
    cplCommLocCodeRootRanksGet
    (
    );

    /**
     *
     * \brief Return the communicator type
     *
     */

    virtual CWP_Comm_t
    commTypeGet() = 0;

    /**
     *
     * \brief Synchronise
     *
     */
/*
    virtual void
    sync
    (
     void *tab,
     MPI_Datatype mpiType,
     int tabSize
    ) = 0;
*/
  protected :

    /**
     *
     * \brief Building coupling communicator from \ref CWP_Comm_t.
     *
     */

    virtual void _cplCommCreate(CWP_Comm_t cplCodeCommType) = 0;


  private :

    /**
     *
     * \brief Assigment operator
     *
     */

    Communication
    &operator=
    (const Communication &other);

    /**
     *
     * \brief Copy constructor
     *
     */

    Communication (const Communication& other);

  protected:

    const CodeProperties *_localCodeProperties; /*!< Pointer to the local code properties */
    const CodeProperties *_cplCodeProperties;   /*!< Pointer to the coupled code properties */

    int       _tag;                     /*!< Tag for MPI */
    MPI_Group _unionGroup;              /*!< Union grou between coupled codes */
    MPI_Comm  _unionComm;               /*!< Union communicator between coupled codes */
    vector <int>  *_unionCommCplRanks;
    vector <int>  *_unionCommLocRanks;

    vector <int>  *_cplCommCplRanks;
    vector <int>  *_cplCommLocRanks;

    MPI_Group _cplGroup;                  /*!< Coupling group (part of merger inter communicator */
    MPI_Comm _cplComm;                    /*!< Coupling communicator (part of merger inter communicator */
    int      _locCodeRootRankUnionComm;   /*!< Root rank associated to the local code into the union communicator */
    int      _cplCodeRootRankUnionComm;   /*!< Root rank associated to the coupled code into the union communicator */
    int      _locCodeRootRankCplComm;   /*!< Root rank associated to the local code into the coupling communicator */
    int      _cplCodeRootRankCplComm;   /*!< Root rank associated to the coupled code into the coupling communicator */
    bool      _isCplRank;                /*!< Is a current rank coupled */
  };
}


/**
 * \endcond
 */

#endif //__COMMUNICATION_H__
