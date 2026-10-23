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

#include "communication.hxx"

#include "cwp.h"

#include "coupling.hxx"
#include "coupling_i.hxx"

#include "pdm_printf.h"

#include <mpi.h>

#include <cassert>
#include <cstring>
#include <string>

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  Communication::Communication()
  : _local_code_properties(nullptr),
    _cpl_code_properties(nullptr),
    _tag(-1),
    _union_comm(MPI_COMM_NULL),
    _cpl_comm(MPI_COMM_NULL),
    _loc_code_root_rank_union_comm(-1),
    _cpl_code_root_rank_union_comm(-1)
  {
  }


  Communication::~Communication()
  {
    if (_cpl_comm != _union_comm) {
      if (!_local_code_properties->isActiveRank() ||
          _local_code_properties->idGet() < _cpl_code_properties->idGet()) {
        if (_union_comm != MPI_COMM_NULL) {
          MPI_Comm_free(&_union_comm);
          _union_comm = MPI_COMM_NULL;
        }
      }
    }

    if (!_local_code_properties->isActiveRank() ||
        _local_code_properties->idGet() < _cpl_code_properties->idGet()) {
      if (_cpl_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&_cpl_comm);
        _cpl_comm = MPI_COMM_NULL;
      }
    }
  }


  void
  Communication::init
  (
    const CodeProperties &localCodeProperties,
    const CodeProperties &cplCodeProperties,
    const string         &cplId,
          CouplingDB     &cplDB
  )
  {
    if (_cpl_comm != MPI_COMM_NULL) {
      return;
    }

    _cpl_code_properties   = &cplCodeProperties;
    _local_code_properties = &localCodeProperties;

    const MPI_Comm &globalComm = localCodeProperties.globalCommGet();

    int globalRank;
    MPI_Comm_rank(globalComm, &globalRank);

    if (!localCodeProperties.isActiveRank()) {
      PDM_printf("Warning CWP_Cpl_create: call CWP_Cpl_create function"
                 " on an uncoupled rank (%d) of the '%s' code\n",
                 globalRank,
                 localCodeProperties.nameGet().c_str());
      return;
    }

    const int localRootRank = localCodeProperties.rootRankGet();
    const int cplRootRank   = cplCodeProperties.rootRankGet();

    // Build a specific tag from cplId
    _tag = tagGet(cplId.c_str(), 0);

    // Build union communicator between coupled codes
    MPI_Group unionGroup;
    if (localCodeProperties.idGet() < cplCodeProperties.idGet()) {
      MPI_Group_union(localCodeProperties.connectableGroupGet(),
                      cplCodeProperties.connectableGroupGet(),
                      &unionGroup);
    }
    else {
      MPI_Group_union(cplCodeProperties.connectableGroupGet(),
                      localCodeProperties.connectableGroupGet(),
                      &unionGroup);
    }
    MPI_Comm_create_group(globalComm, unionGroup, _tag, &_union_comm);

    // Get local and coupled codes communication type
    CWP_Comm_t commType = commTypeGet();
    CWP_Comm_t cplCommType;

    if (globalRank == localRootRank) {
      if (cplCodeProperties.localCodeIs()) {
        if (localRootRank == cplRootRank) {
          // Get directly in memory if both codes on current rank and
          // current rank == root rank for both
          Coupling &distCpl = cplDB.couplingGet(cplCodeProperties, cplId);
          cplCommType = distCpl.commTypeGet();
        }
        else {
          // Exchange to coupled code root rank from local code root rank
          // if local code is on local code root rank and coupled code is on
          // current rank != coupled code root rank
          MPI_Sendrecv(&commType, 1, MPI_INT,
                       cplRootRank, _tag,
                       &cplCommType, 1, MPI_INT,
                       cplRootRank, _tag,
                       globalComm, MPI_STATUS_IGNORE);
        }
      }
      else {
        // Exchange to coupled code root rank from local code root rank
        // if local code is on local code root rank and coupled code is not
        // on current rank
        MPI_Sendrecv(&commType, 1, MPI_INT,
                     cplRootRank, _tag,
                     &cplCommType, 1, MPI_INT,
                     cplRootRank, _tag,
                     globalComm, MPI_STATUS_IGNORE);
      }
    }

    // Broadcast cplCommType from local code root rank to all local code connectable ranks
    MPI_Request request;
    MPI_Ibcast(&cplCommType, 1, MPI_INT, 0,
               localCodeProperties.connectableCommGet(), &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);

    // Broadcast commType from coupled code root rank to all coupled code connectable ranks
    // to account for communication structure copy when both codes are on current rank in
    // coupling.cxx
    if (cplCodeProperties.localCodeIs()) {
      MPI_Ibcast(&commType, 1, MPI_INT, 0,
                 cplCodeProperties.connectableCommGet(), &request);
      MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    // Build coupling communicator
    cplCommCreate(cplCommType);

    // Get root ranks in union communicator
    MPI_Group globalGroup;
    MPI_Comm_group(globalComm, &globalGroup);

    MPI_Group_translate_ranks(globalGroup, 1, &localRootRank,
                              unionGroup, &_loc_code_root_rank_union_comm);

    MPI_Group_translate_ranks(globalGroup, 1, &cplRootRank,
                              unionGroup, &_cpl_code_root_rank_union_comm);
    MPI_Group_free(&unionGroup);
  }


  void
  Communication::init
  (
    const Communication &cplCodeComm
  )
  {
    _local_code_properties         = cplCodeComm._cpl_code_properties;
    _cpl_code_properties           = cplCodeComm._local_code_properties;
    _tag                           = cplCodeComm._tag;
    _union_comm                    = cplCodeComm._union_comm;
    _cpl_comm                      = cplCodeComm._cpl_comm;
    _loc_code_root_rank_union_comm = cplCodeComm._cpl_code_root_rank_union_comm;
    _cpl_code_root_rank_union_comm = cplCodeComm._loc_code_root_rank_union_comm;
  }


  uint32_t
  Communication::tagGet
  (
    const string name,
    const int    offset
  )
  {
    MPI_Aint *maxTagTmp;
    int       flag;

    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &maxTagTmp, &flag);
    int maxTag = static_cast<int>(*maxTagTmp);

    uint32_t mpi_tag = 1;
    if (maxTag > 1) {
      mpi_tag = ((_adler32Compute(name.c_str(), name.size()) + offset) % (maxTag - 1)) + 1;
    }

    return mpi_tag;
  }


  void
  Communication::globalDataBetweenCodesThroughUnionCommExch
  (
    const size_t  s_data,
    const int     n_send_data,
    const void   *send_data,
    const int     n_send_data_cpl,
    const void   *send_data_cpl,
    const int     n_recv_data,
          void   *recv_data,
    const int     n_recv_data_cpl,
          void   *recv_data_cpl
  )
  {
    if (_cpl_code_properties->localCodeIs()) {
      if (n_send_data_cpl > 0) {
        assert(send_data_cpl != nullptr);
      }
      if (n_recv_data_cpl > 0) {
        assert(recv_data_cpl != nullptr);
      }
    }
    else {
      assert(send_data_cpl == nullptr);
      assert(recv_data_cpl == nullptr);
    }

    int unionCommRank;
    MPI_Comm_rank(_union_comm, &unionCommRank);

    int codeID    = _local_code_properties->idGet();
    int cplCodeID = _cpl_code_properties->idGet();

    if (unionCommRank == _loc_code_root_rank_union_comm) {
      if (_cpl_code_properties->localCodeIs()) {
        if (_loc_code_root_rank_union_comm == _cpl_code_root_rank_union_comm) {
          if (codeID < cplCodeID) {
            // Copy data directly if both codes are on current rank and
            // current rank == root rank for both and local code has control by ID
            assert(n_send_data == n_recv_data_cpl);
            memcpy(recv_data_cpl, send_data, s_data * n_send_data);

            assert(n_send_data_cpl == n_recv_data);
            memcpy(recv_data, send_data_cpl, s_data * n_send_data_cpl);
          }
        }
        else {
          // Exchange to coupled code root rank from local code root rank
          // if local code is on local code root rank and coupled code is on
          // current rank != coupled code root rank
          MPI_Sendrecv(send_data,
                       static_cast<int>(s_data) * n_send_data,
                       MPI_UNSIGNED_CHAR,
                       _cpl_code_root_rank_union_comm,
                       _tag,
                       recv_data,
                       static_cast<int>(s_data) * n_recv_data,
                       MPI_UNSIGNED_CHAR,
                       _cpl_code_root_rank_union_comm,
                       _tag,
                       _union_comm,
                       MPI_STATUS_IGNORE);
        }
      }
      else {
        // Exchange to coupled code root rank from local code root rank
        // if local code is on local code root rank and coupled code is not
        // on current rank
        MPI_Sendrecv(send_data,
                     static_cast<int>(s_data) * n_send_data,
                     MPI_UNSIGNED_CHAR,
                     _cpl_code_root_rank_union_comm,
                     _tag,
                     recv_data,
                     static_cast<int>(s_data) * n_recv_data,
                     MPI_UNSIGNED_CHAR,
                     _cpl_code_root_rank_union_comm,
                     _tag,
                     _union_comm,
                     MPI_STATUS_IGNORE);
      }
    }
    else {
      if (_cpl_code_properties->localCodeIs()) {
        if (unionCommRank == _cpl_code_root_rank_union_comm) {
          if (codeID < cplCodeID) {
            // Exchange to local code root rank from coupled code root rank
            // if local code is on current rank and current rank != local code root rank
            // and coupled is on current rank and current rank == coupled code root rank
            MPI_Sendrecv(send_data_cpl,
                         static_cast<int>(s_data) * n_send_data_cpl,
                         MPI_UNSIGNED_CHAR,
                         _loc_code_root_rank_union_comm,
                         _tag,
                         recv_data_cpl,
                         static_cast<int>(s_data) * n_recv_data_cpl,
                         MPI_UNSIGNED_CHAR,
                         _loc_code_root_rank_union_comm,
                         _tag,
                         _union_comm,
                         MPI_STATUS_IGNORE);
          }
        }
      }
    }

    if (!_cpl_code_properties->localCodeIs()) {
      // Broadcast from local code root rank to all local code connectable ranks
      // if coupled code not on current rank
      MPI_Bcast(recv_data,
                static_cast<int>(s_data) * n_recv_data,
                MPI_UNSIGNED_CHAR,
                0,
                _local_code_properties->connectableCommGet());
    }
    else {
      if (codeID < cplCodeID) {
        // Broadcast from local code root rank to all local code connectable ranks
        // if coupled code on current rank and local code has control by ID
        MPI_Bcast(recv_data,
                  static_cast<int>(s_data) * n_recv_data,
                  MPI_UNSIGNED_CHAR,
                  0,
                  _local_code_properties->connectableCommGet());

        // Broadcast from coupled code root rank to all coupled code connectable ranks
        // if coupled code on current rank and local code has control by ID
        MPI_Bcast(recv_data_cpl,
                  static_cast<int>(s_data) * n_recv_data_cpl,
                  MPI_UNSIGNED_CHAR,
                  0,
                  _cpl_code_properties->connectableCommGet());
      }
      else {
        // Broadcast from coupled code root rank to all coupled code connectable ranks
        // if coupled code on current rank and coupled code has control by ID
        MPI_Bcast(recv_data_cpl,
                  static_cast<int>(s_data) * n_recv_data_cpl,
                  MPI_UNSIGNED_CHAR,
                  0,
                  _cpl_code_properties->connectableCommGet());
      }
    }
  }


  void
  Communication::globalDataBetweenCodesThroughUnionCommIssend
  (
    const string       global_data_id,
          MPI_Request *data_send_request,
    const size_t       s_entity,
    const int          stride,
    const int          n_entity,
    const void        *send_data
  )
  {
    assert(send_data != nullptr);

    int unionCommRank;
    MPI_Comm_rank(_union_comm, &unionCommRank);

    if (unionCommRank == _loc_code_root_rank_union_comm) {
      if (_cpl_code_properties->localCodeIs()) {
        if (_loc_code_root_rank_union_comm != _cpl_code_root_rank_union_comm) {
          // Send to coupled code root rank from local code root rank
          // if local code is on local code root rank and coupled code is on
          // current rank != coupled code root rank
          uint32_t data_tag = tagGet(global_data_id, 0);
          MPI_Issend(send_data,
                     static_cast<int>(s_entity) * stride * n_entity,
                     MPI_UNSIGNED_CHAR,
                     _cpl_code_root_rank_union_comm,
                     static_cast<int>(data_tag),
                     _union_comm,
                     data_send_request);
        }
      }
      else {
        // Send to coupled code root rank from local code root rank
        // if local code is on local code root rank and coupled code is not
        // on current rank
        uint32_t data_tag = tagGet(global_data_id, 0);
        MPI_Issend(send_data,
                   static_cast<int>(s_entity) * stride * n_entity,
                   MPI_UNSIGNED_CHAR,
                   _cpl_code_root_rank_union_comm,
                   static_cast<int>(data_tag),
                   _union_comm,
                   data_send_request);
      }
    }
  }


  void
  Communication::globalDataBetweenCodesThroughUnionCommIrecv
  (
    const string       global_data_id,
          MPI_Request *data_recv_request,
    const size_t       s_entity,
    const int          stride,
    const int          n_entity,
          void        *recv_data
  )
  {
    assert(recv_data != nullptr);

    int unionCommRank;
    MPI_Comm_rank(_union_comm, &unionCommRank);

    if (unionCommRank == _loc_code_root_rank_union_comm) {
      if (_cpl_code_properties->localCodeIs()) {
        if (_loc_code_root_rank_union_comm != _cpl_code_root_rank_union_comm) {
          // Recv to coupled code root rank from local code root rank
          // if local code is on local code root rank and coupled code is on
          // current rank != coupled code root rank
          uint32_t data_tag = tagGet(global_data_id, 0);
          MPI_Irecv(recv_data,
                    static_cast<int>(s_entity) * stride * n_entity,
                    MPI_UNSIGNED_CHAR,
                    _cpl_code_root_rank_union_comm,
                    static_cast<int>(data_tag),
                    _union_comm,
                    data_recv_request);
        }
      }
      else {
        // Recv to coupled code root rank from local code root rank
        // if local code is on local code root rank and coupled code is not
        // on current rank
        uint32_t data_tag = tagGet(global_data_id, 0);
        MPI_Irecv(recv_data,
                  static_cast<int>(s_entity) * stride * n_entity,
                  MPI_UNSIGNED_CHAR,
                  _cpl_code_root_rank_union_comm,
                  static_cast<int>(data_tag),
                  _union_comm,
                  data_recv_request);
      }
    }
  }


  void
  Communication::globalDataBetweenCodesThroughUnionCommWaitIssend
  (
          MPI_Request  data_send_request,
    const size_t       s_entity,
    const int          stride,
    const int          n_entity,
    const void        *send_data,
          void        *recv_data
  )
  {
    int unionCommRank;
    MPI_Comm_rank(_union_comm, &unionCommRank);

    int codeID    = _local_code_properties->idGet();
    int cplCodeID = _cpl_code_properties->idGet();

    if (unionCommRank == _loc_code_root_rank_union_comm) {
      if (_cpl_code_properties->localCodeIs()) {
        if (_loc_code_root_rank_union_comm == _cpl_code_root_rank_union_comm) {
          if (codeID < cplCodeID) {
            // Copy data directly if both codes are on current rank and
            // current rank == root rank for both and local code has control by ID
            assert(send_data != nullptr);
            memcpy(recv_data, send_data, s_entity * stride * n_entity);
          }
        }
        else {
          // Send to coupled code root rank from local code root rank
          // if local code is on local code root rank and coupled code is on
          // current rank != coupled code root rank
          MPI_Wait(&data_send_request, MPI_STATUS_IGNORE);
        }
      }
      else {
        // Send to coupled code root rank from local code root rank
        // if local code is on local code root rank and coupled code is not
        // on current rank
        MPI_Wait(&data_send_request, MPI_STATUS_IGNORE);
      }
    }
  }


  void
  Communication::globalDataBetweenCodesThroughUnionCommWaitIrecv
  (
          MPI_Request  data_recv_request,
    const size_t       s_entity,
    const int          stride,
    const int          n_entity,
    const void        *send_data,
          void        *recv_data
  )
  {
    int unionCommRank;
    MPI_Comm_rank(_union_comm, &unionCommRank);

    int codeID    = _local_code_properties->idGet();
    int cplCodeID = _cpl_code_properties->idGet();

    if (unionCommRank == _loc_code_root_rank_union_comm) {
      if (_cpl_code_properties->localCodeIs()) {
        if (_loc_code_root_rank_union_comm == _cpl_code_root_rank_union_comm) {
          if (codeID < cplCodeID) {
            // Copy data directly if both codes are on current rank and
            // current rank == root rank for both and local code has control by ID
            assert(send_data != nullptr);
            memcpy(recv_data, send_data, s_entity * stride * n_entity);
          }
        }
        else {
          // Recv to coupled code root rank from local code root rank
          // if local code is on local code root rank and coupled code is on
          // current rank != coupled code root rank
          MPI_Wait(&data_recv_request, MPI_STATUS_IGNORE);
        }
      }
      else {
        // Recv to coupled code root rank from local code root rank
        // if local code is on local code root rank and coupled code is not
        // on current rank
        MPI_Wait(&data_recv_request, MPI_STATUS_IGNORE);
      }
    }

    // Broadcast from local code root rank to all local code connectable ranks
    MPI_Bcast(recv_data,
              static_cast<int>(s_entity) * stride * n_entity,
              MPI_UNSIGNED_CHAR,
              0,
              _local_code_properties->connectableCommGet());
  }


  uint32_t
  Communication::_adler32Compute
  (
    const void   *buf,
    const size_t  buflength
  )
  {
    const uint8_t *buffer = static_cast<const uint8_t *>(buf);

    uint32_t s1 = 1;
    uint32_t s2 = 0;

    for (size_t n = 0; n < buflength; n++) {
      s1 = (s1 + buffer[n]) % 65521;
      s2 = (s2 + s1) % 65521;
    }

    return (s2 << 16) | s1;
  }

}

/**
 * \endcond
 */