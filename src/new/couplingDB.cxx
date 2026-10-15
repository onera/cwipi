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

#include "couplingDB.hxx"

#include "codeProperties.hxx"
#include "coupling.hxx"
#include "couplingDB_i.hxx"

#include "pdm_error.h"

#include <map>
#include <string>

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  CouplingDB::CouplingDB()
  : _couplingDB(*(new map<const CodeProperties *, map<string, Coupling *>>()))
  {
  }


  CouplingDB::~CouplingDB()
  {
    typedef map<const CodeProperties *, map<string, Coupling *>>::iterator Iterator;
    typedef map<string, Coupling *>::iterator Iterator2;

    for (Iterator p1 = _couplingDB.begin();
         p1 != _couplingDB.end(); p1++) {
      for (Iterator2 p = p1->second.begin();
           p != p1->second.end(); p++) {
        if (p->second != nullptr)
          delete p->second;
      }
      p1->second.clear();
    }

    _couplingDB.clear();
    delete &_couplingDB;
  }


  void
  CouplingDB::couplingCreate
  (
          CodeProperties       &localCodeProperties,
    const string               &cplId,
          CodeProperties       &coupledCodeProperties,
    const CWP_Interface_t       entities_dim,
    const CWP_Comm_t            commType,
    const CWP_Spatial_interp_t  spatialInterpAlgo,
    const int                   nPart,
    const CWP_Dynamic_mesh_t    movingStatus,
    const CWP_Time_exch_t       recvFreqType
  )
  {
    if (couplingIs(localCodeProperties, cplId)) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' existing coupling\n", cplId.c_str());
    }

    Coupling *newCoupling = new Coupling(cplId,
                                         commType,
                                         localCodeProperties,
                                         coupledCodeProperties,
                                         entities_dim,
                                         spatialInterpAlgo,
                                         nPart,
                                         movingStatus,
                                         recvFreqType,
                                         *this);

    map<string, Coupling *>& codeMap = _couplingDB[&localCodeProperties];
    pair<string, Coupling *> newPair(string(cplId), newCoupling);
    codeMap.insert(newPair);
  }


  void
  CouplingDB::couplingDel
  (
    const CodeProperties &localCodeProperties,
    const string         &cplId
  )
  {
    typedef const map<const CodeProperties *, map<string, Coupling *>>::iterator Iterator;
    Iterator p = _couplingDB.find(&localCodeProperties);

    if (p == _couplingDB.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' coupling not found for '%s' code\n", cplId.c_str(),
                localCodeProperties.nameGet().c_str());
    }
    else {
      typedef map<string, Coupling *>::iterator Iterator2;
      Iterator2 p1 = p->second.find(cplId);

      if (p1 == p->second.end()) {
        PDM_error(__FILE__, __LINE__, 0,
                  "'%s' coupling not found '%s' code\n", cplId.c_str(),
                  localCodeProperties.nameGet().c_str());
      }

      if (p1->second != nullptr) {
        delete p1->second;
      }
      p->second.erase(p1);
    }
  }


  void
  CouplingDB::timeStepBegin
  (
    const CodeProperties &localCodeProperties,
          double          current_time
  )
  {
    typedef const map<const CodeProperties *, map<string, Coupling *>>::iterator Iterator;
    Iterator p = _couplingDB.find(&localCodeProperties);

    if (p == _couplingDB.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' code not found in couplingDB\n",
                localCodeProperties.nameGet().c_str());
    }

    map<string, Coupling *>::iterator itc = p->second.begin();
    while (itc != p->second.end()) {
      itc->second->time_step_beg(current_time);
      itc++;
    }
  }


  void
  CouplingDB::timeStepEnd
  (
    const CodeProperties &localCodeProperties
  )
  {
    typedef const map<const CodeProperties *, map<string, Coupling *>>::iterator Iterator;
    Iterator p = _couplingDB.find(&localCodeProperties);

    if (p == _couplingDB.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' code not found in couplingDB\n",
                localCodeProperties.nameGet().c_str());
    }

    map<string, Coupling *>::iterator itc = p->second.begin();
    while (itc != p->second.end()) {
      itc->second->time_step_end();
      itc++;
    }
  }

}

/**
 * \endcond
 */