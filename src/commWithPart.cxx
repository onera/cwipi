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

#include "commWithPart.hxx"

using namespace std;

namespace cwipi {

  /**
   *
   * \brief Constructor.
   *
   */
  
  CommWithPart::CommWithPart()
    : Communication::Communication(), 
      _commType(CWP_COMM_PAR_WITH_PART)
  {
  }

  /**
   * \brief Destructor.
   *
   */

  CommWithPart::~CommWithPart()
  {
  }

  /**
   *
   * \brief Building coupling communicator.
   *
   */

  void 
  CommWithPart::_cplCommCreate
  (
   CWP_Comm_t cplCodeCommType
  )
  {
    
    _isCplRank = true;
    
    if (cplCodeCommType != CWP_COMM_PAR_WITH_PART) {
      
      const vector <int> &cplRanks = *(_cplCodeProperties->connectableRanksGet());
      
      vector <int> exRanks(cplRanks.size()-1);
      
      int j = 0;
      for (int i = 0; i < cplRanks.size(); i++) {
        if (cplRanks[i] != _cplCodeProperties->rootRankGet()) {
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

    else {

      _cplComm = _unionComm;
      _cplGroup = MPI_Comm_group(_cplComm, &_cplGroup);
      
    }
  }

  /**
   *
   * \brief Synchronise
   *
   */

#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable:869)
#elif defined(__clang__)	
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value" 	
#endif
  void
  CommWithPart::sync
  (
   void *tab, 
   MPI_Datatype mpiType, 
   int tabSize
  )
  {
  }
#if defined(__INTEL_COMPILER)
#pragma warning(pop)
#elif defined(__clang__)	
#pragma clang diagnostic pop
#endif
}
