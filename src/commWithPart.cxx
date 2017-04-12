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

    int mergeInterCommSize;
    MPI_Comm_size(_mergeInterComm, &mergeInterCommSize);
    
    const int localFirstRank = _localCodeProperties->firstRankGet();
    const int localLastRank  = _localCodeProperties->lastRankGet();

    const int nLocalRank = localLastRank
                         - localFirstRank
                         + 1;

    const int cplCodeFirstRank = _cplCodeProperties->firstRankGet();
    const int cplCodeLastRank  = _cplCodeProperties->lastRankGet();
   
    const int nCplCodeRank = cplCodeLastRank
                           - cplCodeFirstRank
                           + 1;

    if (cplCodeCommType != CWP_COMM_PAR_WITH_PART) {

      _cplCodeNRankCplComm = 1;

      //
      // Store coupling active ranks to create the coupling communicator

      int  nRankList = 1 + nLocalRank;
      int *rankList  = new int[nRankList]; 

      if (localFirstRank < cplCodeFirstRank) {
        _cplCodeFirstRankCplComm = nLocalRank;
        for (int i = 0; i < nLocalRank + 1; i++)
          rankList[i] = i;
      }
      else {
        rankList[0] = 0;
        _cplCodeFirstRankCplComm = 0;
        for (int i = 0; i < nLocalRank; i++)
          rankList[i+1] = nCplCodeRank + i;
      }

      //
      // Create the coupling communicator from a new group

      MPI_Group mergeGroup = MPI_GROUP_NULL;
      MPI_Group cplGroup   = MPI_GROUP_NULL;
        
      MPI_Comm_group(_mergeInterComm, &mergeGroup);
        
      MPI_Group_incl(mergeGroup, nRankList, rankList, &cplGroup);
        
      MPI_Comm_create(_mergeInterComm, cplGroup, &_cplComm);
        
      MPI_Group_free(&cplGroup);
      MPI_Group_free(&mergeGroup);
        
      delete [] rankList;
    }

    else {

      _cplComm = _mergeInterComm;
      _cplCodeNRankCplComm = cplCodeLastRank 
                           - cplCodeFirstRank
                           + 1;

      if (localFirstRank < cplCodeFirstRank) {
        _cplCodeFirstRankCplComm = nLocalRank;
      }
      else {
        _cplCodeFirstRankCplComm = 0;
      }

    }

    //
    // Build fvm communicator

    MPI_Comm_dup(_localCodeProperties->intraCommGet(), &_fvmComm);
  }

  /**
   *
   * \brief Synchronise
   *
   */

  void
  CommWithPart::sync
  (
   void *tab, 
   MPI_Datatype mpiType, 
   int tabSize
  )
  {
  }
}
