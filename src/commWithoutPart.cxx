/*
  This file is part of the CWIPI library. 

  Copyright (C) 2013  ONERA

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

using namespace std;

namespace cwipi {

  /**
   *
   * \brief Constructor.
   *
   */
  
  CommWithoutPart::CommWithoutPart()
    : Communication::Communication(), 
      _commType(CWP_COMM_PAR_WITH_PART)
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

    _isCplRank = localRootRank == currentRank;
    
    if (cplCodeCommType != CWP_COMM_PAR_WITH_PART) {

      MPI_Group globalGroup;
      MPI_Comm_group(globalComm, &globalGroup);

      MPI_Group unionGroup;
      MPI_Comm_group(_unionComm, &unionGroup);      

      int cplRanks[2];
      int gap1 = 0;
      int gap2 = 1;
      
      if (_localCodeProperties->idGet() < _cplCodeProperties->idGet()) {
        gap1 = 1;
        gap2 = 0;
      }
      
      MPI_Group_translate_ranks (globalGroup, 1, &localRootRank, 
                                 unionGroup, cplRanks + gap1);

      MPI_Group_translate_ranks (globalGroup, 1, &cplRootRank, 
                                 unionGroup, cplRanks + gap2);
      
      MPI_Group_incl(unionGroup, 2, cplRanks, &_cplGroup);
      
      MPI_Comm_create (_unionComm, _cplGroup, &_cplComm);
      
    }
    
    else {
      
      const vector <int> &cplRanks = *(_localCodeProperties->connectableRanksGet());
      
      vector <int> exRanks(cplRanks.size()-1);
      
      int j = 0;
      for (int i = 0; i < cplRanks.size(); i++) {
        if (cplRanks[i] != _localCodeProperties->rootRankGet()) {
          exRanks[j++] = cplRanks[i]; 
        }      
      }
      
      vector <int> tExRanks(exRanks.size());
      
      MPI_Group globalGroup;
      MPI_Comm_group(_localCodeProperties->globalCommGet(), &globalGroup);
      
      MPI_Group unionGroup;
      MPI_Comm_group(_unionComm, &unionGroup);      
      
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

  }

}
