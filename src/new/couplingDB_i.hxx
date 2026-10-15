#ifndef __COUPLING_DB_I_H__
#define __COUPLING_DB_I_H__
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

#include "pdm_error.h"

#include <cassert>

/**
 * \cond
 */

namespace cwipi {

  inline Coupling&
  CouplingDB::couplingGet
  (
    const CodeProperties &localCodeProperties,
    const std::string    &cplId
  )
  {
    typedef const std::map<const CodeProperties *, std::map<std::string, Coupling *>>::iterator Iterator;
    Iterator p = _couplingDB.find(&localCodeProperties);

    if (p == _couplingDB.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' code not found in couplingDB\n", localCodeProperties.nameGet().c_str());
    }

    typedef std::map<std::string, Coupling *>::iterator Iterator2;
    Iterator2 p1 = p->second.find(cplId);

    if (p1 == p->second.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' coupling not found '%s' code\n", cplId.c_str(),
                localCodeProperties.nameGet().c_str());
    }

    assert(p1->second != nullptr);
    return *p1->second;
  }


  inline std::map<std::string, Coupling *>&
  CouplingDB::couplingListGet
  (
    const CodeProperties &localCodeProperties
  )
  {
    typedef const std::map<const CodeProperties *, std::map<std::string, Coupling *>>::iterator Iterator;
    Iterator p = _couplingDB.find(&localCodeProperties);

    if (p == _couplingDB.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' code not found in couplingDB\n", localCodeProperties.nameGet().c_str());
    }

    return p->second;
  }


  inline bool
  CouplingDB::couplingIs
  (
    const CodeProperties &localCodeProperties,
    const std::string    &cplId
  )
  {
    bool status = true;
    typedef const std::map<const CodeProperties *, std::map<std::string, Coupling *>>::iterator Iterator;
    Iterator p = _couplingDB.find(&localCodeProperties);

    if (p == _couplingDB.end()) {
      status = false;
    }
    else {
      typedef std::map<std::string, Coupling *>::iterator Iterator2;
      Iterator2 p1 = p->second.find(cplId);

      if (p1 == p->second.end()) {
        status = false;
      }
    }

    return status;
  }

}

/**
 * \endcond
 */

#endif // __COUPLING_DB_I_H__