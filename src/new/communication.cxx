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

#include "communication.hxx"
#include "coupling.hxx"
#include "pdm_printf.h"


/**
 * \cond
 */

using namespace std;

namespace cwipi {

  /**
   *
   * \brief Constructor.
   *
   */

  Communication::Communication()
    :_localCodeProperties(NULL),
     _cplCodeProperties(NULL),
     _tag(-1),
     _unionGroup(MPI_GROUP_NULL),
     _unionComm(MPI_COMM_NULL),
//     _fvmComm(MPI_COMM_NULL),
     _cplComm(MPI_COMM_NULL),
     _locCodeRootRankUnionComm(-1),
     _cplCodeRootRankUnionComm(-1),
     _isCplRank(false),
     _cplCommCplRanks(NULL),
     _cplCommLocRanks(NULL)
  {
  }

  /**
   *
   * \brief Destructor.
   *
   */

  Communication::~Communication()
  {

    if (_unionComm != MPI_COMM_NULL)
      MPI_Comm_free(&_unionComm);

    if (_cplComm != MPI_COMM_NULL)
      MPI_Comm_free(&_cplComm);

//    if (_fvmComm != MPI_COMM_NULL)
//      MPI_Comm_free(&_fvmComm);
  }

  /**
   *
   * \brief Initialize coupling communicators.
   *
   * \param [in]  localCodeProperties   Local code properties
   * \param [in]  cplCodeProperties     Coupled code properties
   * \param [in]  cplId                 Coupling identifier
   *
   */

  void
  Communication::init
  (
   const CodeProperties &localCodeProperties,
   const CodeProperties &cplCodeProperties,
   const string         &cplId,
   CouplingDB           &cplDB
  )
  {

    _isCplRank = localCodeProperties.isActiveRank();
    if (_cplComm == MPI_COMM_NULL) {
      _cplCodeProperties = &cplCodeProperties;
      _localCodeProperties = &localCodeProperties;

      const int localRootRank = localCodeProperties.rootRankGet();
      const int cplRootRank = cplCodeProperties.rootRankGet();

      const MPI_Comm& globalComm = localCodeProperties.globalCommGet();

      int globalRank;
      MPI_Comm_rank(globalComm, &globalRank);

      if (!localCodeProperties.isActiveRank()) {
        PDM_printf(
           "Warning CWP_Cpl_create : Call CWP_Cpl_create function"
           " on an uncoupled rank (%d) of the '%s' code\n",
            globalRank,
            localCodeProperties.nameGet().c_str());
      }
      else {

        //Build a specific tag through the cplId
        _tag = 0;
        for (size_t i = 0; i < cplId.size(); i++) {
          _tag += cplId[i];
        }

        if (MPI_TAG_UB > 0) {
          _tag = _tag % MPI_TAG_UB;
        }

        // Build the union communicator between the two coupled codes

        if (localCodeProperties.idGet() < cplCodeProperties.idGet()) {
          MPI_Group_union (localCodeProperties.connectableGroupGet(),
                           cplCodeProperties.connectableGroupGet(),
                           &_unionGroup);
        }
        else {
          MPI_Group_union (cplCodeProperties.connectableGroupGet(),
                           localCodeProperties.connectableGroupGet(),
                           &_unionGroup);
        }

        MPI_Comm_create_group(globalComm, _unionGroup, _tag, &_unionComm);

        int mergeInterCommSize;

        MPI_Comm_size(_unionComm, &mergeInterCommSize);

        CWP_Comm_t commType = commTypeGet();

        CWP_Comm_t cplCommType;

        if (globalRank == localRootRank) {
          if (cplCodeProperties.localCodeIs()) {
            if (localRootRank == cplRootRank) {
              Coupling &distCpl = cplDB.couplingGet(cplCodeProperties, cplId);
              cplCommType = distCpl.commTypeGet() ;
            }
            else {
              MPI_Sendrecv (&commType, 1, MPI_INT,
                          cplRootRank, _tag,
                          &cplCommType, 1, MPI_INT,
                           cplRootRank, _tag,
                           globalComm, MPI_STATUS_IGNORE);
            }
          }

          else {
            MPI_Sendrecv (&commType, 1, MPI_INT,
                          cplRootRank, _tag,
                          &cplCommType, 1, MPI_INT,
                           cplRootRank, _tag,
                           globalComm, MPI_STATUS_IGNORE);
          }
        }

        MPI_Request request1;
        MPI_Request request2;

        MPI_Ibcast(&cplCommType, 1, MPI_INT, 0,
                  localCodeProperties.connectableCommGet(), &request1);
        MPI_Wait(&request1, MPI_STATUS_IGNORE);

        if (cplCodeProperties.localCodeIs()) {
          CWP_Comm_t &cplcommType2 = commType;
          MPI_Ibcast(&cplcommType2, 1, MPI_INT, 0,
                    cplCodeProperties.connectableCommGet(), & request2);
          MPI_Wait(&request2, MPI_STATUS_IGNORE);
        }

        //
        // Build the coupling communicator

        _cplCommCreate(cplCommType);
        MPI_Group globalGroup;
        MPI_Comm_group(globalComm, &globalGroup);

        if (_isCplRank) {
          printf("I am coupled\n");
          MPI_Group_translate_ranks (globalGroup, 1, &localRootRank,
                                     _unionGroup, &_locCodeRootRankUnionComm);

          MPI_Group_translate_ranks (globalGroup, 1, &cplRootRank,
                                     _unionGroup, &_cplCodeRootRankUnionComm);
        }

      }
    }
  }


  /**
   *
   * \brief Initialize coupling communicators.
   *
   * \param [in]  cplCodeComm           Coupled code communication
   *
   */

  void
  Communication::init
  (
   Communication &cplCodeComm
  )
  {
    _localCodeProperties    = cplCodeComm._cplCodeProperties;
    _cplCodeProperties      = cplCodeComm._localCodeProperties;

    _tag                    = cplCodeComm._tag;
    _unionGroup             = cplCodeComm._unionGroup;
    _unionComm              = cplCodeComm._unionComm;
    _cplGroup               = cplCodeComm._cplGroup;
    _cplComm                = cplCodeComm._cplComm;
    _locCodeRootRankUnionComm = cplCodeComm._cplCodeRootRankUnionComm;
    _cplCodeRootRankUnionComm = cplCodeComm._locCodeRootRankUnionComm;
    _unionCommCplRanks = new std::vector<int>(*(cplCodeComm._unionCommLocRanks));
    _unionCommLocRanks = new std::vector<int>(*(cplCodeComm._unionCommCplRanks));
    _cplCommCplRanks = new std::vector<int>(*(cplCodeComm._cplCommLocRanks));
    _cplCommLocRanks = new std::vector<int>(*(cplCodeComm._cplCommCplRanks));


    const int localRootRank = _localCodeProperties->rootRankGet();

    const MPI_Comm& globalComm = _localCodeProperties->globalCommGet();

    int currentRank;
    MPI_Comm_rank(globalComm, &currentRank);

    if (this->commTypeGet() != CWP_COMM_PAR_WITH_PART) {
      _isCplRank = localRootRank == currentRank;
    }
    else {
      _isCplRank = true;
    }

  }

  MPI_Comm
  Communication::unionCommGet
  (
  )
  {
    return _unionComm;
  }

  MPI_Comm
  Communication::cplCommGet
  (
  )
  {
    return _cplComm;
  }


  std::vector<int>*
  Communication::unionCommCplRanksGet
  (
  )
  {
    return _unionCommCplRanks;
  }

  std::vector<int>*
  Communication::cplCommCplRanksGet
  (
  )
  {
    return _cplCommCplRanks;
  }

  std::vector<int>*
  Communication::unionCommLocRanksGet
  (
  )
  {
    return _unionCommLocRanks;
  }

  std::vector<int>*
  Communication::cplCommLocRanksGet
  (
  )
  {
    return _cplCommLocRanks;
  }

  int
  Communication::unionCommCplCodeRootRanksGet
  (
  )
  {
     return  _cplCodeRootRankUnionComm;
  }

  int
  Communication::cplCommCplCodeRootRanksGet
  (
  )
  {
     return  _cplCodeRootRankCplComm;
  }

  int
  Communication::unionCommLocCodeRootRanksGet
  (
  )
  {
    return _locCodeRootRankUnionComm;
  }

  int
  Communication::cplCommLocCodeRootRanksGet
  (
  )
  {
    return _locCodeRootRankCplComm;
  }


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
  Communication::iexchGlobalDataBetweenCodesThroughUnionCom
  (
   size_t       s_data,
   int          n_send_data,
   void        *send_data,
   int          n_recv_data,
   void        *recv_data,
   MPI_Request *request
  )
  {


    int unionCommRank;
    MPI_Status status;
    MPI_Comm_rank (_unionComm, &unionCommRank);

    if (unionCommRank ==  _locCodeRootRankUnionComm) {
      if (_cplCodeProperties->localCodeIs()) {
        if (_locCodeRootRankUnionComm == _cplCodeRootRankUnionComm) {
          assert (n_send_data == n_recv_data);
          memcpy (recv_data, send_data, s_data * n_send_data);
        }
        else {
          MPI_Sendrecv (send_data,
                        (int) s_data * n_send_data,
                        MPI_UNSIGNED_CHAR,
                        _locCodeRootRankUnionComm,
                        _tag,
                        recv_data,
                        (int) s_data * n_recv_data,
                        MPI_UNSIGNED_CHAR,
                        _cplCodeRootRankUnionComm,
                        _tag,
                        _unionComm,
                        &status);
        }
      }
      else {
        MPI_Sendrecv (send_data,
                      (int) s_data * n_send_data,
                      MPI_UNSIGNED_CHAR,
                      _locCodeRootRankUnionComm,
                      _tag,
                      recv_data,
                      (int) s_data * n_recv_data,
                      MPI_UNSIGNED_CHAR,
                      _cplCodeRootRankUnionComm,
                      _tag,
                      _unionComm,
                      &status);
      }
    }
   
    // BCast in the intraComm (Use Ibcast to be sure)

    MPI_Ibcast (recv_data,
               (int) s_data * n_recv_data,
               MPI_UNSIGNED_CHAR,
               0,
               _localCodeProperties->connectableCommGet(),
               request);

  }

}

/**
 * \endcond
 */
