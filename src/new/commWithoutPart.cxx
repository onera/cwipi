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
  CommWithoutPart::_cplCommCreate
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

      const int localRootRank = _localCodeProperties->rootRankGet();
      const int cplRootRank   = _cplCodeProperties->rootRankGet();

      int cplRanks2[2];
      int gap1 = 0;
      int gap2 = 1;
      if (_localCodeProperties->idGet() < _cplCodeProperties->idGet()) {
        gap1 = 1;
        gap2 = 0;
      }

      // Translate localRootRank from globalGroup into cplRanks2[gap1] in unionGroup
      MPI_Group_translate_ranks(globalGroup, 1, &localRootRank,
                                unionGroup, cplRanks2 + gap1);

      // Translate cplRootRank from globalGroup into cplRanks2[gap2] in unionGroup
      MPI_Group_translate_ranks(globalGroup, 1, &cplRootRank,
                                unionGroup, cplRanks2 + gap2);

      // Include root ranks of the local and coupled code to obtain the cplComm communicator
      MPI_Group cplGroup;
      if (cplRanks2[0] == cplRanks2[1]) {
        MPI_Group_incl(unionGroup, 1, &cplRanks2[0], &cplGroup);
      }
      else {
        MPI_Group_incl(unionGroup, 2, cplRanks2, &cplGroup);
      }
      MPI_Comm_create(_unionComm, cplGroup, &_cplComm);
      MPI_Group_free(&cplGroup);

    }
    else {

      // Exclude the connectable ranks (unionComm) of the local code to obtain the cplComm communicator
      // containing only the root rank of the local code

      const int locRootRankInGlobalComm = _localCodeProperties->rootRankGet();
      int locRootRankInUnionComm;
      MPI_Group_translate_ranks(globalGroup, 1, &locRootRankInGlobalComm,
                                unionGroup, &locRootRankInUnionComm);
      int *excludeRanks = static_cast<int *>(malloc(sizeof(int) * (locRanks.size() - 1)));
      int n_exclude = 0;
      for (size_t i = 0; i < locRanks.size(); i++) {
        if ((*_unionCommLocRanks)[i] != locRootRankInUnionComm) {
          excludeRanks[n_exclude++] = (*_unionCommLocRanks)[i];
        }
      }

      MPI_Group cplGroup;
      MPI_Group_excl(unionGroup, n_exclude, excludeRanks, &cplGroup);
      MPI_Comm_create(_unionComm, cplGroup, &_cplComm);
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