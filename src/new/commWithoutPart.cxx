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

#include "commWithoutPart.hxx"

#include "cwp.h"

#include <mpi.h>

#include <cstddef>
#include <cstdlib>
#include <vector>

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  CommWithoutPart::CommWithoutPart()
  : Communication::Communication(),
    _commType(CWP_COMM_PAR_WITHOUT_PART)
  {
  }


  CommWithoutPart::~CommWithoutPart()
  {
  }


  void
  CommWithoutPart::cplCommCreate
  (
    CWP_Comm_t cplCodeCommType
  )
  {
    MPI_Group globalGroup;
    MPI_Comm_group(_local_code_properties->globalCommGet(), &globalGroup);

    MPI_Group unionGroup;
    MPI_Comm_group(_union_comm, &unionGroup);

    if (cplCodeCommType != CWP_COMM_PAR_WITH_PART) {

      const int localRootRank = _local_code_properties->rootRankGet();
      const int cplRootRank   = _cpl_code_properties->rootRankGet();

      int cplRanks[2];
      int gap1 = 0;
      int gap2 = 1;
      if (_local_code_properties->idGet() < _cpl_code_properties->idGet()) {
        gap1 = 1;
        gap2 = 0;
      }

      // Translate localRootRank from globalGroup into cplRanks[gap1] in unionGroup
      MPI_Group_translate_ranks(globalGroup, 1, &localRootRank,
                                unionGroup, cplRanks + gap1);

      // Translate cplRootRank from globalGroup into cplRanks[gap2] in unionGroup
      MPI_Group_translate_ranks(globalGroup, 1, &cplRootRank,
                                unionGroup, cplRanks + gap2);

      // Include root ranks of the local and coupled code to obtain the cplComm communicator
      MPI_Group cplGroup;
      if (cplRanks[0] == cplRanks[1]) {
        MPI_Group_incl(unionGroup, 1, &cplRanks[0], &cplGroup);
      }
      else {
        MPI_Group_incl(unionGroup, 2, cplRanks, &cplGroup);
      }
      MPI_Comm_create(_union_comm, cplGroup, &_cpl_comm);
      MPI_Group_free(&cplGroup);

    }
    else {

      vector<int> locRanks = *(_local_code_properties->connectableRanksGet());

      // Translate locRanks from globalGroup into unionCommLocRanks in unionGroup
      vector<int> unionCommLocRanks(locRanks.size());
      MPI_Group_translate_ranks(globalGroup, locRanks.size(), &(locRanks[0]),
                                unionGroup, &(unionCommLocRanks[0]));

      // Exclude the connectable ranks (unionComm) of the local code to obtain the cplComm communicator
      // containing only the root rank of the local code
      const int locRootRankInGlobalComm = _local_code_properties->rootRankGet();
      int locRootRankInUnionComm;
      MPI_Group_translate_ranks(globalGroup, 1, &locRootRankInGlobalComm,
                                unionGroup, &locRootRankInUnionComm);
      int *excludeRanks = static_cast<int *>(malloc(sizeof(int) * (locRanks.size() - 1)));
      int n_exclude = 0;
      for (size_t i = 0; i < locRanks.size(); i++) {
        if (unionCommLocRanks[i] != locRootRankInUnionComm) {
          excludeRanks[n_exclude++] = unionCommLocRanks[i];
        }
      }

      MPI_Group cplGroup;
      MPI_Group_excl(unionGroup, n_exclude, excludeRanks, &cplGroup);
      MPI_Comm_create(_union_comm, cplGroup, &_cpl_comm);
      MPI_Group_free(&cplGroup);
      free(excludeRanks);

    }

    MPI_Group_free(&unionGroup);
    MPI_Group_free(&globalGroup);

  }

}

/**
 * \endcond
 */