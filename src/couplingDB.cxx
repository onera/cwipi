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
    : _couplingDB(*new map < const CodeProperties *, map <string, Coupling * > > ())
  {
  }

  /**
   * \brief Destructor
   *
   */

  CouplingDB::~CouplingDB()
  {
    typedef map < const CodeProperties *, map <string, Coupling * > >::iterator Iterator;
    typedef map < string, Coupling * > ::iterator Iterator2;

    for (Iterator p1 = _couplingDB.begin();
           p1 != _couplingDB.end(); p1++) {
      for (Iterator2 p = p1->second.begin();
           p != p1->second.end(); p++) {
        if (p->second != NULL)
          delete p->second;
      }
      p1->second.clear();
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
   * \param [in]  cplCodeProperties  Coupled code properties
   * \param [in]  commType           Communication type
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
   const CodeProperties        &localCodeProperties,
   const string                &cplId,
   const CodeProperties        &coupledCodeProperties,
   const CWP_Comm_t            commType,
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

//    const map <string, CodeProperties * >::iterator p = 
//      _couplingDB.find(cplId);
    
//    if (p == _couplingDB.end()) {
//      pair<string, Coupling* >
//        newPair(string(cplId), newCoupling);
//
//      _couplingDB.insert(newPair);
//    }
//    
//    else {
//      
//      if (p->second->)
//      
//      bftc_error(__FILE__, __LINE__, 0,
//                "'%s' existing coupling\n", cplId.c_str());
//    }

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
   const CodeProperties &localCodeProperties,
   const string &cplId
  )
  {
//    const map <string, Coupling * >::iterator p = _couplingDB.find(cplId);
//    if (p == _couplingDB.end())
//      bftc_error(__FILE__, __LINE__, 0,
//                "'%s' coupling not found \n", cplId.c_str());
//
//    if (p->second != NULL)
//      delete p->second;
//
//    _couplingDB.erase(p);
  }
}
