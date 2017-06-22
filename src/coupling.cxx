/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011-2017  ONERA

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

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstring>
#include <sstream>

#include <bftc_error.h>
#include <bftc_file.h>

#include <fvmc_parall.h>

#include "coupling.hxx"
#include "coupling_i.hxx"

#include "oldMesh.hxx"
#include "codeProperties.hxx"

#include "solve_ax_b_4.h"
#include "quickSort.h"
#include "cwp.h"

#include "factory.hpp"
#include "geometry.hxx"
#include "communication.hxx"


/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux) &&  !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

using namespace std;

namespace cwipi {

  typedef Factory<Communication, CWP_Comm_t> FC;
  typedef Factory<Geometry, CWP_Geom_t> FG;
  typedef Factory<Support, CWP_Support_t> FS;

  
  Coupling::Coupling
  (
   const string                &cplId,
   const CWP_Comm_t           cplType,
   const CodeProperties        &localCodeProperties,
   const CodeProperties        &coupledCodeProperties,
   const CWP_Geom_t           geomAlgo,
   const CWP_Support_t        supportType,
   const int                    nPart,
   const CWP_Displacement_t  movingStatus,
   const CWP_Freq_t           recvFreqType,
   CouplingDB                 &cplDB          
   )
  :_cplId(cplId),
   _commType(cplType),
   _communication(*(FC::getInstance().CreateObject(cplType))),
   _localCodeProperties(localCodeProperties),
   _coupledCodeProperties(coupledCodeProperties),
   _recvFreqType (recvFreqType),
   _cplDB(cplDB)   
//   _geometry(*(FG::getInstance().CreateObject(geomAlgo))),
{
  
  if (coupledCodeProperties.localCodeIs()) {
    if (cplDB.couplingIs(coupledCodeProperties, cplId)) {
      _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB);
      Coupling &distCpl = cplDB.couplingGet(coupledCodeProperties, cplId);
      distCpl._communication.init(_communication);
    }
  }
  
  else {
    _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB);
  }
    
}


Coupling::~Coupling()
{
#if defined(DEBUG) && 0
  cout << "destroying '" << _name << "' coupling : TODO" << endl;
#endif

}


} // namespace cwipi

