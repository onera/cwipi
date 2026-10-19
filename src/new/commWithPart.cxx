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

#include "commWithPart.hxx"

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

  CommWithPart::CommWithPart()
  : Communication::Communication(),
    _commType(CWP_COMM_PAR_WITH_PART)
  {
  }


  CommWithPart::~CommWithPart()
  {
  }


  void
  CommWithPart::cplCommCreate
  (
    const CWP_Comm_t cplCodeCommType
  )
  {
    MPI_Group globalGroup;
    MPI_Comm_group(_local_code_properties->globalCommGet(), &globalGroup);

    if (cplCodeCommType != CWP_COMM_PAR_WITH_PART) {

      MPI_Group unionGroup;
      MPI_Comm_group(_union_comm, &unionGroup);

      vector<int> cplRanks = *(_cpl_code_properties->connectableRanksGet());

      // Translate cplRanks from globalGroup into unionCommCplRanks in unionGroup
      vector<int> unionCommCplRanks(cplRanks.size());
      MPI_Group_translate_ranks(globalGroup, cplRanks.size(), &(cplRanks[0]),
                                unionGroup, &(unionCommCplRanks[0]));

      // Exclude the connectable ranks (unionComm) of the coupled code to obtain the cplComm communicator
      // containing only the root rank of the coupled code
      const int cplRootRankInGlobalComm = _cpl_code_properties->rootRankGet();
      int cplRootRankInUnionComm;
      MPI_Group_translate_ranks(globalGroup, 1, &cplRootRankInGlobalComm,
                                unionGroup, &cplRootRankInUnionComm);
      int *excludeRanks = static_cast<int *>(malloc(sizeof(int) * (cplRanks.size() - 1)));
      int n_exclude = 0;
      for (size_t i = 0; i < cplRanks.size(); i++) {
        if (unionCommCplRanks[i] != cplRootRankInUnionComm) {
          excludeRanks[n_exclude++] = unionCommCplRanks[i];
        }
      }

      MPI_Group cplGroup;
      MPI_Group_excl(unionGroup, n_exclude, excludeRanks, &cplGroup);
      MPI_Comm_create(_union_comm, cplGroup, &_cpl_comm);
      MPI_Group_free(&cplGroup);
      MPI_Group_free(&unionGroup);
      free(excludeRanks);

    }
    else {
      _cpl_comm = _union_comm;
    }

    MPI_Group_free(&globalGroup);

  }

}

/**
 * \endcond
 */