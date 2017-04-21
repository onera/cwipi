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

#include <bftc_printf.h>

#include "singleton.hpp"
#include "couplingDB.hxx"
#include "couplingDB_i.hxx"
#include "codeProperties.hxx"
#include "coupling.hxx"

using namespace std;

namespace cwipi {

  /**
   * \brief Constructor
   *
   */

  CouplingDB::CouplingDB()
    : _couplingDB(*new map <string, Coupling * > ())
  {
  }

  /**
   * \brief Destructor
   *
   */

  CouplingDB::~CouplingDB()
  {
    typedef map <string, Coupling * >::iterator Iterator;
    for (Iterator p = _couplingDB.begin();
         p != _couplingDB.end(); p++) {
      if (p->second != NULL)
        delete p->second;
    }
    _couplingDB.clear();

    delete &_couplingDB;
  }

  /**
   * \brief Building and storage a coupling object 
   *
   * This function creates a coupling object and defines its properties.
   *
   * \param [in]  cplId              Coupling identifier
   * \param [in]  commType           Communication type
   * \param [in]  cplCodeProperties  Coupled code properties
   * \param [in]  geomAlgo           Geometric algorithm
   * \param [in]  supportType        Support type
   * \param [in]  nPart              Number of interface partition 
   * \param [in]  movingStatus       Support moving status
   * \param [in]  recvFreqType       Type of receiving frequency
   *
   */

  void  
  CouplingDB::couplingCreate
  (
   const string                &cplId,
   const CWP_Comm_t            commType,
   const CodeProperties        &localCodeProperties,
   const CodeProperties        &coupledCodeProperties,
   const CWP_Geom_t           geomAlgo,
   const CWP_Support_t   supportType,
   const int                    nPart,
   const CWP_Displacement_t  movingStatus,
   const CWP_Freq_t           recvFreqType
  )
  {

    //
    // Create the new coupling

    Coupling *newCoupling = NULL;
    //TODO: Call new Coupling when the constructor will be available
//    Coupling *newCoupling = new Coupling(cplId,
//                                         commType,
//                                         localCodeProperties,
//                                         coupledCodeProperties,
//                                         geomAlgo,
//                                         supportType,
//                                         nPart,
//                                         movingStatus,
//                                         recvFreqType);

    pair<string, Coupling* >
      newPair(string(cplId), newCoupling);

    pair<map<string, Coupling* >::iterator, bool>
      p = _couplingDB.insert(newPair);

    if (!p.second)
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' existing coupling\n", cplId.c_str());

  }

  /**
   * \brief Deletion a coupling object int the database.
   *
   * \param [in]  cplId              Coupling identifier
   *
   */

  void  
  CouplingDB::couplingDel
  (
   const string &cplId
  )
  {
    const map <string, Coupling * >::iterator p = _couplingDB.find(cplId);
    if (p == _couplingDB.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' coupling not found \n", cplId.c_str());

    if (p->second != NULL)
      delete p->second;

    _couplingDB.erase(p);
  }
}
