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

#include "commWithoutPart.hxx"
#include "pdm_logging.h"

using namespace std;

namespace cwipi {

  /**
   *
   * \brief Constructor.
   *
   */
  
  CommWithoutPart::CommWithoutPart()
    : Communication::Communication(), 
      _commType(CWP_COMM_PAR_WITHOUT_PART)
  {
  }

  /**
   * \brief Destructor.
   *
   */

  CommWithoutPart::~CommWithoutPart()
  {
  }

  /**
   *
   * \brief Building coupling communicator.
   *
   */

  void 
  CommWithoutPart::_cplCommCreate
  (
   CWP_Comm_t cplCodeCommType
  )
  {
    const int localRootRank = _localCodeProperties->rootRankGet();
    const int cplRootRank   = _cplCodeProperties->rootRankGet();
    
    const MPI_Comm& globalComm = _localCodeProperties->globalCommGet();

    int currentRank;
    MPI_Comm_rank(globalComm, &currentRank);

    //_isCplRank = localRootRank == currentRank;

    vector <int> cplRanks = *(_cplCodeProperties->connectableRanksGet());
    vector <int> locRanks = *(_localCodeProperties->connectableRanksGet());

    _unionCommCplRanks = new std::vector<int>(cplRanks.size());
    _unionCommLocRanks = new std::vector<int>(locRanks.size());


    //Traduction de ces rangs dans _unionCommCplRanks - Contient les rangs (unionComm) correspondants

    MPI_Group globalGroup;
    MPI_Comm_group(_localCodeProperties->globalCommGet(), &globalGroup);
      
    MPI_Group unionGroup;
    MPI_Comm_group(_unionComm, &unionGroup);      

    MPI_Group_translate_ranks(globalGroup, cplRanks.size(), &(cplRanks[0]),
                              unionGroup, &((*_unionCommCplRanks)[0]));

    MPI_Group_translate_ranks(globalGroup, locRanks.size(), &(locRanks[0]),
                              unionGroup , &((*_unionCommLocRanks)[0]));

    _cplCommCplRanks = new std::vector<int>(*_unionCommCplRanks);
    _cplCommLocRanks = new std::vector<int>(*_unionCommLocRanks);
      
    if (cplCodeCommType != CWP_COMM_PAR_WITH_PART) {
      log_trace("cplCodeCommType != CWP_COMM_PAR_WITH_PART\n");

      int cplRanks2[2];
      int gap1 = 0;
      int gap2 = 1;
      
      if (_localCodeProperties->idGet() < _cplCodeProperties->idGet()) {
        gap1 = 1;
        gap2 = 0;
      }
      
      MPI_Group_translate_ranks (globalGroup, 1, &localRootRank, 
                                 unionGroup, cplRanks2 + gap1);

      MPI_Group_translate_ranks (globalGroup, 1, &cplRootRank, 
                                 unionGroup, cplRanks2 + gap2);
      
      MPI_Group_incl(unionGroup, 2, cplRanks2, &_cplGroup);
      
      MPI_Comm_create (_unionComm, _cplGroup, &_cplComm);
      
    }
  
    else {
      log_trace("cplCodeCommType == CWP_COMM_PAR_WITH_PART\n");
      log_trace("cplRanks.size() = %d\n", cplRanks.size());
      log_trace("_localCodeProperties->rootRankGet() = %d\n", _localCodeProperties->rootRankGet());
          
      vector <int> exRanks(cplRanks.size());//-1);
      log_trace("exRanks : %p, size = %d\n", exRanks, exRanks.size());

      int j = 0;
      for (size_t i = 0; i < cplRanks.size(); i++) {
        log_trace("i = %d, cplRanks[i] = %d\n", i, cplRanks[i]);
        if (cplRanks[i] != _localCodeProperties->rootRankGet()) {
          log_trace("  j = %d\n", j);
          exRanks[j++] = cplRanks[i]; 
        }      
      }
      
      vector <int> tExRanks(exRanks.size());
          
      MPI_Group_translate_ranks(globalGroup, exRanks.size(), &(exRanks[0]),
                                unionGroup, &(tExRanks[0]));
      MPI_Group_excl(unionGroup, exRanks.size(), &(tExRanks[0]), &_cplGroup);
      
      MPI_Comm_create(_unionComm, _cplGroup, &_cplComm);
      
    }
  }

  /**
   *
   * \brief Synchronise
   *
   */
/*
  void
  CommWithoutPart::sync
  (
   void *tab, 
   MPI_Datatype mpiType, 
   int tabSize
  )
  {
    MPI_Bcast(tab,
              tabSize,
              mpiType,
              0,
              _localCodeProperties->intraCommGet());

  }*/

}
