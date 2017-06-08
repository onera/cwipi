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

#include "communication.hxx"

using namespace std;

namespace cwipi {

  /**
   *
   * \brief Constructor.
   *
   */

  Communication::Communication()
    :_isCplRank(false), 
     _mergeInterComm(MPI_COMM_NULL),
     _fvmComm(MPI_COMM_NULL),
     _cplComm(MPI_COMM_NULL),
     _cplCodeNRankCplComm(0),
     _cplCodeFirstRankCplComm(0),
     _localCodeProperties(NULL),
     _cplCodeProperties(NULL)
  {
  }

  /**
   *
   * \brief Destructor.
   *
   */

  Communication::~Communication()
  {
    if (_mergeInterComm != MPI_COMM_NULL)
      MPI_Comm_free(&_mergeInterComm);

    if (_cplComm != MPI_COMM_NULL)
      MPI_Comm_free(&_cplComm);

    if (_fvmComm != MPI_COMM_NULL)
      MPI_Comm_free(&_fvmComm);
  }

  /**
   *
   * \brief Initialize coupling communicators.
   *
   * \param [in]  localCodeProperties   Local code properties
   * \param [in]  cplCodeProperties     Coupled code properties
   *
   */

  void 
  Communication::init
  (
   CodeProperties &localCodeProperties, 
   CodeProperties &cplCodeProperties 
  )
  {

    _cplCodeProperties = &cplCodeProperties;
    _localCodeProperties = &localCodeProperties;

    const MPI_Comm& localComm = localCodeProperties.intraCommGet();

    if (_cplComm == MPI_COMM_NULL) {

      const int localRootRank = localCodeProperties.rootRankGet();
      
      const MPI_Comm& globalComm = localCodeProperties.globalCommGet();
      
      int currentRank;
      MPI_Comm_rank(globalComm, &currentRank);

      const int cplCodeRootRank = cplCodeProperties.rootRankGet();
      
      //
      // Define coupling communicator
      

      const int tag = 'm'+'e'+'r'+'g'+'e'+'I'+'n'+'t'+'e'+'r'+'C'+'o'+'m'+'m';

      //
      // Build the inter communicator between the two coupled codes

      MPI_Comm tmpInterComm;
//
//      MPI_Intercomm_create(localComm, 0, globalComm, cplCodeFirstRank, tag, &tmpInterComm);
//
//      //
//      // Merge the inter communicator to obtain an intra communicator
//
//      int sorting;
//      if (localFirstRank < cplCodeFirstRank) {
//        sorting = 0;
//        _cplCodeFirstRankCplComm = nLocalRank ;
//      }
//      else {
//        sorting = 1;
//        _cplCodeFirstRankCplComm = 0 ;
//      }
//
//      MPI_Intercomm_merge(tmpInterComm, sorting, &_mergeInterComm);
//
//      if (tmpInterComm != MPI_COMM_NULL)
//        MPI_Comm_free(&tmpInterComm);
//
//      int mergeInterCommSize;
//
//      MPI_Comm_size(_mergeInterComm, &mergeInterCommSize);
//      
//      CWP_Comm_t* commTypes =
//        new CWP_Comm_t[mergeInterCommSize];
//
//      CWP_Comm_t commType = commTypeGet();
//
//      MPI_Allgather((void*)& commType,
//                    1,
//                    MPI_INT,
//                    commTypes,
//                    1,
//                    MPI_INT,
//                    _mergeInterComm);
//
//      //
//      // Check coupling type
//      
//      int firstRank = 0;
//      int lastRank = 0;
//      
//      if (localFirstRank < cplCodeFirstRank) {
//        firstRank = 0;
//        lastRank = nLocalRank;
//      }
//      else {
//        firstRank = nCplCodeRank;
//        lastRank = nLocalRank + nCplCodeRank;
//      }
//      
//      for (int i = firstRank; i < lastRank; i++)
//        if (commTypes[i] != commType)
//          bftc_error(__FILE__, __LINE__, 0,
//                     "Two different communication types for the '%s' application\n",
//                     localCodeProperties.nameGet().c_str());
//
//      CWP_Comm_t cplCodeCommType;
//      
//      if (localFirstRank < cplCodeFirstRank)
//        cplCodeCommType = commTypes[nLocalRank];
//      else
//        cplCodeCommType= commTypes[0];
//      
//      delete [] commTypes;
//
//      _cplCommCreate(cplCodeCommType);
//
//    }
  }

  //     _isCplRank = localCodeProperties.getFirstRank() == currentRank ||
  //       _couplingType == CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING;
      
  //     //
  //     // Build coupling communicator
      
  //     if (_couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING ||
  //         cplCodeCouplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {
        
  //       int *rankList = new int[coupledCodeCommSize];
  //       int nRankList;
        
  //       if (_couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING &&
  //           cplCodeCouplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {
          
  //         nRankList = 2;
  //         rankList[0] = 0;
          
  //         _cplCodeNRankCplComm = 1;
  //         if (localFirstRank < cplCodeFirstRank) {
  //           rankList[1] = nLocalRank;
  //           _cplCodeFirstRankCplComm = 1;
  //         }
  //         else {
  //           rankList[1] = nCplCodeRank;
  //           _cplCodeFirstRankCplComm = 0;
  //         }
  //       }
        
  //       else if (cplCodeCouplingType == CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {
  //         nRankList = 1 + nCplCodeRank;
  //         _cplCodeNRankCplComm = nCplCodeRank;
  //         if (localFirstRank < cplCodeFirstRank) {
  //           rankList[0] = 0;
  //           _cplCodeFirstRankCplComm = 1;
  //           for (int i = 0; i < nCplCodeRank; i++)
  //             rankList[i+1] = nLocalRank + i;
  //         }
  //         else {
  //           _cplCodeFirstRankCplComm = 0;
  //           for (int i = 0; i < nCplCodeRank; i++)
  //             rankList[i] = i;
  //           rankList[nCplCodeRank] = nCplCodeRank;
  //         }
  //       }
        
  //       else if (_couplingType == CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {
  //         nRankList = 1 + nLocalRank;
  //         _cplCodeNRankCplComm = 1;
  //         if (localFirstRank < cplCodeFirstRank) {
  //           _cplCodeFirstRankCplComm = nLocalRank;
  //           for (int i = 0; i < nLocalRank; i++)
  //             rankList[i] = i;
  //         rankList[nLocalRank] = nLocalRank;
  //         }
          
  //         else {
  //           rankList[0] = 0;
  //           _cplCodeFirstRankCplComm = 0;
  //           for (int i = 0; i < nLocalRank; i++)
  //             rankList[i+1] = nCplCodeRank + i;
  //         }
  //       }

  //       else {
          
  //         bftc_error(__FILE__, __LINE__, 0,
  //                    "Error in 'build coupling communicator'\n");
  //       }
        
  //       MPI_Group mergeGroup = MPI_GROUP_NULL;
  //       MPI_Group couplingGroup = MPI_GROUP_NULL;
        
  //       MPI_Comm_group(_mergeInterComm, &mergeGroup);
        
  //       MPI_Group_incl(mergeGroup, nRankList, rankList, &couplingGroup);
        
  //       MPI_Comm_create(_mergeInterComm, couplingGroup, &_cplComm);
        
  //       MPI_Group_free(&couplingGroup);
  //       MPI_Group_free(&mergeGroup);
        
  //       delete [] rankList;
        
  //     }
  //     else
  //       MPI_Comm_dup(_mergeInterComm, &_cplComm);
  //   }


  //   if (_fvmComm == MPI_COMM_NULL) {
      
  //     if (_couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {
        
  //       int list1 = 0;
  //       MPI_Group localGroup = MPI_GROUP_NULL;
  //       MPI_Group fvmGroup = MPI_GROUP_NULL;
        
  //       MPI_Comm dupLocalComm = MPI_COMM_NULL;
  //       MPI_Comm_dup(localComm, &dupLocalComm);
        
  //       MPI_Comm_group(dupLocalComm, &localGroup);
  //       MPI_Group_incl(localGroup, 1, &list1, &fvmGroup);
  //       MPI_Comm_create(localComm,
  //                       fvmGroup,
  //                       &_fvmComm);
  //       MPI_Group_free(&localGroup);
  //       MPI_Group_free(&fvmGroup);
  //       if (dupLocalComm != MPI_COMM_NULL)
  //         MPI_Comm_free(&dupLocalComm);
  //     }

  //     else
        
  //       MPI_Comm_dup(localComm, &_fvmComm);
      
  //   }
  // }
  }
}

