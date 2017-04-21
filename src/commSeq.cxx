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

#include "commSeq.hxx"

using namespace std;

namespace cwipi {

  /**
   *
   * \brief Constructor.
   *
   */
  
  CommSeq::CommSeq()
    : CommWithoutPart::CommWithoutPart(), 
      _commType(CWP_COMM_SEQ)
  {
  }

  /**
   * \brief Destructor.
   *
   */

  CommSeq::~CommSeq()
  {
  }

  /**
   *
   * \brief Building coupling communicator.
   *
   */

  void 
  CommSeq::_cplCommCreate
  (
   CWP_Comm_t cplCodeCommType
  )
  {
      
//FIXME: Check that CommSeq::_cplCommCreate is same that CommWithoutPart::_cplCommCreate
      
//    int mergeInterCommSize;
//    MPI_Comm_size(_mergeInterComm, &mergeInterCommSize);
//
//    const int localFirstRank = _localCodeProperties->firstRankGet();
//    const int localLastRank  = _localCodeProperties->lastRankGet();
//    
//    const MPI_Comm& globalComm = _localCodeProperties->globalCommGet();
//
//    int currentRank;
//    MPI_Comm_rank(globalComm, &currentRank);
//
//    _isCplRank = localFirstRank == currentRank;
// 
//    const int nLocalRank = localLastRank
//                         - localFirstRank
//                         + 1;
//
//    const int cplCodeFirstRank = _cplCodeProperties->firstRankGet();
//    const int cplCodeLastRank  = _cplCodeProperties->lastRankGet();
//   
//    const int nCplCodeRank = cplCodeLastRank
//                           - cplCodeFirstRank
//                           + 1;
//
//    int  nRankList = 0;
//    int *rankList = NULL;
//    
//    //
//    // Store coupling active ranks to create the coupling communicator
//
//    if (cplCodeCommType == CWP_COMM_PAR_WITH_PART) {
//
//      _cplCodeNRankCplComm = 1;
//
//      nRankList = 1 + nCplCodeRank;
//      rankList  = new int[nRankList]; 
//
//      if (localFirstRank < cplCodeFirstRank) {
//        rankList[0] = 0;
//        _cplCodeFirstRankCplComm = 1;
//        for (int i = 0; i < nCplCodeRank; i++)
//          rankList[i+1] = nLocalRank + i;
//      }
//      else {
//        _cplCodeFirstRankCplComm = 0;
//        for (int i = 0; i < nCplCodeRank; i++)
//          rankList[i] = i;
//        rankList[nCplCodeRank] = nCplCodeRank;
//      }
//
//    }
//
//    else {
//
//      nRankList = 2;
//      rankList[0] = 0;
//          
//      _cplCodeNRankCplComm = 1;
//      if (localFirstRank < cplCodeFirstRank) {
//        rankList[1] = nLocalRank;
//        _cplCodeFirstRankCplComm = 1;
//      }
//      else {
//        rankList[1] = nCplCodeRank;
//        _cplCodeFirstRankCplComm = 0;
//      }
//
//    }
//
//    //
//    // Create the coupling communicator from rankList
//
//    MPI_Group mergeGroup = MPI_GROUP_NULL;
//    MPI_Group cplGroup   = MPI_GROUP_NULL;
//        
//    MPI_Comm_group(_mergeInterComm, &mergeGroup);
//        
//    MPI_Group_incl(mergeGroup, nRankList, rankList, &cplGroup);
//        
//    MPI_Comm_create(_mergeInterComm, cplGroup, &_cplComm);
//        
//    MPI_Group_free(&cplGroup);
//    MPI_Group_free(&mergeGroup);
//        
//    delete [] rankList;
//
//    //
//    // Build the fvm communicator from the local communicator 
//    // (Only the master rank is into fvm communicator)
//
//    int list1 = 0;
//    MPI_Group localGroup = MPI_GROUP_NULL;
//    MPI_Group fvmGroup   = MPI_GROUP_NULL;
//    MPI_Comm localComm   = _localCodeProperties->intraCommGet();
//    
//    MPI_Comm_group(localComm, &localGroup);
//    MPI_Group_incl(localGroup, 1, &list1, &fvmGroup);
//    MPI_Comm_create(localComm,
//                    fvmGroup,
//                    &_fvmComm);
//    MPI_Group_free(&localGroup);
//    MPI_Group_free(&fvmGroup);

  }


#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable:869)
#elif defined(__clang__)	
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value" 	
#endif

  /**
   *
   * \brief Synchronise.
   *
   */

  void
  CommSeq::sync
  (
   void *tab, 
   MPI_Datatype mpiType, 
   int tabSize
  )
  {
    //
    // Nothing to do
    //
  }
#if defined(__INTEL_COMPILER)
#pragma warning(pop)
#elif defined(__clang__)	
#pragma clang diagnostic pop
#endif
}
