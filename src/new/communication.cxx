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
     _isCplRank(false)
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
    if (_cplComm == MPI_COMM_NULL) {
      _cplCodeProperties = &cplCodeProperties;
      _localCodeProperties = &localCodeProperties;

      const int localRootRank = localCodeProperties.rootRankGet();
      const int cplRootRank = cplCodeProperties.rootRankGet();

      const MPI_Comm& globalComm = localCodeProperties.globalCommGet();
      
      int globalRank;
      MPI_Comm_rank(globalComm, &globalRank);
      
      if (!localCodeProperties.isCoupledRank()) {
        PDM_printf(
           "Warning CWP_Cpl_create : Call CWP_Cpl_create function"
           " on an uncoupled rank (%d) of the '%s' code\n",
            globalRank,
            localCodeProperties.nameGet().c_str());
      }

      else {
      
        //Build a specific tag through the cplId
        int tag = 0;
        for (size_t i = 0; i < cplId.size(); i++) {
          tag += cplId[i];
        }

        if (MPI_TAG_UB > 0) {
          tag = tag % MPI_TAG_UB;
        }
    
        //
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
        
        MPI_Comm_create_group(globalComm, _unionGroup, tag, &_unionComm);

        int mergeInterCommSize;

        MPI_Comm_size(_unionComm, &mergeInterCommSize);

        CWP_Comm_t commType = commTypeGet();
        
        CWP_Comm_t cplCommType = CWP_COMM_PAR_WITH_PART;

        if (globalRank == localRootRank) {
          if (cplCodeProperties.localCodeIs()) {
            if (localRootRank == cplRootRank) {
              Coupling &distCpl = cplDB.couplingGet(cplCodeProperties, cplId);
              cplCommType = distCpl.commTypeGet() ;
            }
            else {
              MPI_Sendrecv (&commType, 1, MPI_INT,
                          cplRootRank, tag,
                          &cplCommType, 1, MPI_INT,
                           cplRootRank, tag,
                           globalComm, MPI_STATUS_IGNORE);            
            }
          }
          
          else {
            MPI_Sendrecv (&commType, 1, MPI_INT,
                          cplRootRank, tag,
                          &cplCommType, 1, MPI_INT,
                           cplRootRank, tag,
                           globalComm, MPI_STATUS_IGNORE);
          }
        }
        
        MPI_Bcast(&cplCommType, 1, MPI_INT, 0, 
                  localCodeProperties.connectableCommGet());
        
        if (cplCodeProperties.localCodeIs()) {
          CWP_Comm_t &cplcommType2 = commType;
          MPI_Bcast(&cplcommType2, 1, MPI_INT, 0, 
                    cplCodeProperties.connectableCommGet());          
        }
        
        //
        // Build the coupling communicator

        _cplCommCreate(cplCommType);

        MPI_Group globalGroup;
        MPI_Comm_group(globalComm, &globalGroup);

        if (_isCplRank) {
        
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

  std::vector<int>*
  Communication::unionCommCplRanksGet
  (
  )
  {
    return _unionCommCplRanks;
  }

  std::vector<int>*
  Communication::unionCommLocRanksGet
  (
  )
  {
    return _unionCommLocRanks;
  }

  int
  Communication::unionCommCplCodeRootRanksGet
  (
  )
  {
     return  _cplCodeRootRankUnionComm;
  }

  int
  Communication::unionCommLocCodeRootRanksGet
  (
  )
  {
    return _locCodeRootRankUnionComm;
  }
  

}

