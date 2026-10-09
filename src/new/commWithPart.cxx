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
  CommWithPart::_cplCommCreate
  (
    CWP_Comm_t cplCodeCommType
  )
  {
    MPI_Group globalGroup;
    MPI_Comm_group(_localCodeProperties->globalCommGet(), &globalGroup);

    MPI_Group unionGroup;
    MPI_Comm_group(_unionComm, &unionGroup);

    vector<int> cplRanks = *(_cplCodeProperties->connectableRanksGet());
    vector<int> locRanks = *(_localCodeProperties->connectableRanksGet());

    _unionCommCplRanks = new vector<int>(cplRanks.size());
    _unionCommLocRanks = new vector<int>(locRanks.size());

    // Translate cplRanks from globalGroup into _unionCommCplRanks in unionGroup
    MPI_Group_translate_ranks(globalGroup, cplRanks.size(), &(cplRanks[0]),
                              unionGroup, &((*_unionCommCplRanks)[0]));

    // Translate locRanks from globalGroup into _unionCommLocRanks in unionGroup
    MPI_Group_translate_ranks(globalGroup, locRanks.size(), &(locRanks[0]),
                              unionGroup, &((*_unionCommLocRanks)[0]));

    _cplCommCplRanks = new vector<int>(*_unionCommCplRanks);
    _cplCommLocRanks = new vector<int>(*_unionCommLocRanks);

    if (cplCodeCommType != CWP_COMM_PAR_WITH_PART) {

      // Exclude the connectable ranks (unionComm) of the coupled code to obtain the cplComm communicator
      // containing only the root rank of the coupled code

      const int cplRootRankInGlobalComm = _cplCodeProperties->rootRankGet();
      int cplRootRankInUnionComm;
      MPI_Group_translate_ranks(globalGroup, 1, &cplRootRankInGlobalComm,
                                unionGroup, &cplRootRankInUnionComm);
      int *excludeRanks = static_cast<int *>(malloc(sizeof(int) * (cplRanks.size() - 1)));
      int n_exclude = 0;
      for (size_t i = 0; i < cplRanks.size(); i++) {
        if ((*_unionCommCplRanks)[i] != cplRootRankInUnionComm) {
          excludeRanks[n_exclude++] = (*_unionCommCplRanks)[i];
        }
      }

      MPI_Group cplGroup;
      MPI_Group_excl(unionGroup, n_exclude, excludeRanks, &cplGroup);
      MPI_Comm_create(_unionComm, cplGroup, &_cplComm);
      MPI_Group_free(&cplGroup);
      free(excludeRanks);

    }
    else {
      _cplComm = _unionComm;
    }

    MPI_Group_free(&unionGroup);
    MPI_Group_free(&globalGroup);

  }

}

/**
 * \endcond
 */