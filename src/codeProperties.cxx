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
#include <mpi.h>

#include "codeProperties.hxx"
#include "bftc_printf.h"

namespace cwipi
{

  /**
   * \brief Constructor.
   *
   * \param [in]  name         Current code name
   * \param [in]  id           Identifier
   * \param [in]  rootRank     Root rank in global communicator
   * \param [in]  isLocal      Is a local code
   * \param [in]  globalComm   MPI communicator containing all processes 
   *                           of all codes
   *
   */

  CodeProperties::CodeProperties
  (
   string &name,
   int    id,   
   int    rootRank,   
   bool   isLocal,   
   const MPI_Comm globalComm
  ): _name(name), _id(id), _isLocal(isLocal),
     _rootRankInGlobalComm(rootRank),  
     _globalComm(globalComm),
     _isCoupledRank(false), 
     _intCtrlParam(*(new map <string, int>())),
     _dblCtrlParam(*(new map <string, double>())),
     _strCtrlParam(*(new map <string, string>()))
  {
    _intraComm = MPI_COMM_NULL;
    _intraCoupledGroup = MPI_GROUP_NULL;
    _intraGroup        = MPI_GROUP_NULL;
    
    int globalCommSize;
    MPI_Comm_size(globalComm, &globalCommSize);
    
  }

  /**
   * \brief Copy constructor.
   *
   * \param [in]  other        other code properties
   *
   */

  CodeProperties::CodeProperties
  (
   const CodeProperties& other
  ): _name(other._name), _id(other._id), _isLocal(other._isLocal),
     _rootRankInGlobalComm(other._rootRankInGlobalComm),  
     _globalComm(other._globalComm), 
     _intraComm(other._intraComm),
     _isCoupledRank(other._isCoupledRank), 
     _intraCoupledGroup(other._intraCoupledGroup),
     _intraGroup(other._intraGroup),
     _intCtrlParam(other._intCtrlParam),
     _dblCtrlParam(other._dblCtrlParam),
     _strCtrlParam(other._strCtrlParam)
  {
  }

  /**
   * \brief Dump properties
   *
   */

  void 
  CodeProperties::dump()
  {
    bftc_printf("'%s' properties\n",_name.c_str());
//    bftc_printf("  - Ranks in global MPI_comm : %i <= ranks <= %i \n",
//               _firstRank,
//               _lastRank);
    bftc_printf("  - Identifier %d :\n", _id);
    bftc_printf("  - Int Control Parameter :\n");

    typedef map <string, int>::iterator Iterator1;
    for (Iterator1 p = _intCtrlParam.begin(); 
                   p != _intCtrlParam.end(); p++)
      bftc_printf("   * '%s' : %i\n", p->first.c_str(), p->second);
    bftc_printf("  - Double Control Parameter :\n");

    typedef map <string, double>::iterator Iterator2;
    for (Iterator2 p = _dblCtrlParam.begin(); 
                   p != _dblCtrlParam.end(); p++)
      bftc_printf("   * '%s' : %12.5e\n", p->first.c_str(), p->second);

    bftc_printf("  - String Control Parameter :\n");

    typedef map <string, string>::iterator Iterator3;
    for (Iterator3 p = _strCtrlParam.begin(); 
                   p != _strCtrlParam.end(); p++)
      bftc_printf("   * '%s' : '%s'\n", p->first.c_str(), p->second.c_str());

    bftc_printf_flush();
  }

  /**
   * \brief Destructor
   *
   */

  CodeProperties::~CodeProperties()
  {
    printf("~CodeProperties\n");
    if (!_intCtrlParam.empty())
      _intCtrlParam.clear();
    delete &_intCtrlParam;
    if (!_dblCtrlParam.empty())
      _dblCtrlParam.clear();
    delete &_dblCtrlParam;
    if (!_strCtrlParam.empty())
      _strCtrlParam.clear();
    delete &_strCtrlParam;
    if (_intraComm != MPI_COMM_NULL)
      MPI_Comm_free(&_intraComm);
    if (_intraGroup != MPI_GROUP_NULL)
      MPI_Group_free(&_intraGroup);
    if (_intraCoupledGroup != MPI_GROUP_NULL)
      MPI_Group_free(&_intraCoupledGroup);
  }
}
