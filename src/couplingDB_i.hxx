#ifndef __COUPLING_DB_I_H__
#define __COUPLING_DB_I_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

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

#include <cassert>

#include <bftc_error.h>

using namespace std;

namespace cwipi {

  /**
   * \brief Return a coupling object from it identifier
   *
   * \param [in]  cplId              Coupling identifier
   *
   */

  Coupling& 
  CouplingDB::couplingGet
  (
   const string &cplId
  ) 
  {
    const map <string, Coupling * >::iterator p = _couplingDB.find(cplId);
    if (p == _couplingDB.end())
      bftc_error(__FILE__, __LINE__, 0, 
                "'%s' coupling not found \n", cplId.c_str());
    assert( p->second != NULL);
    return *p->second;
  }

}

#endif
