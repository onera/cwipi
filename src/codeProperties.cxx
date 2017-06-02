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
#include <sstream>
#include <iostream>
#include <algorithm>

#include "codeProperties.hxx"
#include "bftc_printf.h"

using namespace std;

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
     _winIntParamIdxName(MPI_WIN_NULL),
     _winIntParamName(MPI_WIN_NULL),
     _winIntParamValue(MPI_WIN_NULL),
     _winIntParamIdxNameData(NULL),
     _winIntParamNameData(NULL),
     _winIntParamValueData(NULL),
     _winDoubleParamIdxName(MPI_WIN_NULL),
     _winDoubleParamName(MPI_WIN_NULL),
     _winDoubleParamValue(MPI_WIN_NULL),
     _winDoubleParamIdxNameData(NULL),
     _winDoubleParamNameData(NULL),
     _winDoubleParamValueData(NULL),
     _winStrParamIdxName(MPI_WIN_NULL),
     _winStrParamName(MPI_WIN_NULL),
     _winStrParamIdxValue(MPI_WIN_NULL),
     _winStrParamValue(MPI_WIN_NULL),
     _winStrParamIdxNameData(NULL),
     _winStrParamNameData(NULL),
     _winStrParamIdxValueData(NULL),
     _winStrParamValueData(NULL)
  {
    _intraComm = MPI_COMM_NULL;
    _intraCoupledGroup = MPI_GROUP_NULL;
    _intraGroup        = MPI_GROUP_NULL;
    
    _winGlobData[0] = 0; // Unlock parameters access
    _winGlobData[1] = 0; // 0 int param
    _winGlobData[2] = 0; // 0 doube param
    _winGlobData[3] = 0; // 0 str param
    
    _winIntParamIdxNameData = (int *) malloc(sizeof(int));
    _winDoubleParamIdxNameData = (int *) malloc(sizeof(int));
    _winStrParamIdxNameData = (int *) malloc(sizeof(int));
    _winStrParamIdxValueData = (int *) malloc(sizeof(int));

    _winIntParamIdxNameData[0] = 0;
    _winDoubleParamIdxNameData[0] = 0;
    _winStrParamIdxNameData[0] = 0;
    _winStrParamIdxValueData[0] = 0;
    
    cout << "toto : " << name << " " << isLocal << endl;
    
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
     _winIntParamIdxName(other._winIntParamIdxName),
     _winIntParamName(other._winIntParamName),
     _winIntParamValue(other._winIntParamValue),
     _winIntParamIdxNameData(other._winIntParamIdxNameData),
     _winIntParamNameData(other._winIntParamNameData),
     _winIntParamValueData(other._winIntParamValueData),
     _winDoubleParamIdxName(other._winDoubleParamIdxName),
     _winDoubleParamName(other._winDoubleParamName),
     _winDoubleParamValue(other._winDoubleParamValue),
     _winDoubleParamIdxNameData(other._winDoubleParamIdxNameData),
     _winDoubleParamNameData(other._winDoubleParamNameData),
     _winDoubleParamValueData(other._winDoubleParamValueData),
     _winStrParamIdxName(other._winStrParamIdxName),
     _winStrParamName(other._winStrParamName),
     _winStrParamIdxValue(other._winStrParamIdxValue),
     _winStrParamValue(other._winStrParamValue),
     _winStrParamIdxNameData(other._winStrParamIdxNameData),
     _winStrParamNameData(other._winStrParamNameData),
     _winStrParamIdxValueData(other._winStrParamIdxValueData),
     _winStrParamValueData(other._winStrParamValueData)

  {
    
    memcpy(_winGlobData, other._winGlobData, 4 *sizeof(int));
    
  }

  /**
   * \brief Dump properties
   *
   */

  void 
  CodeProperties::dump()
  {
    bftc_printf ("'%s' properties\n",_name.c_str());
    bftc_printf ("  - Identifier : %d\n", _id);
    bftc_printf ("  - Is it a local code : %d\n", _isLocal);
    bftc_printf ("  - Root rank in global communicator : %d\n", _rootRankInGlobalComm);
    bftc_printf ("  - Root rank in global communicator : %d\n", _rootRankInGlobalComm);
    bftc_printf ("  - Is local code : %d\n", _isLocal);
    bftc_printf ("  - Is it a coupled rank : %d\n", _isCoupledRank);

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    _updateIntValues ();   
    //_updateDoubleValues ();   
    //_updateStrValues ();   

    bftc_printf ("  - Integer control parameters \n");
    bftc_printf ("     * Number : %d\n", _winGlobData[1]);

    int sParamMax = -1;
    for (int i = 0; i < _winGlobData[1]; i++) {
      int sParam = _winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i];
      sParamMax = max(sParam, sParamMax);
    }

    if (sParamMax > 99) sParamMax == 99;

    char fmtName[24];
    sprintf(fmtName, "     * '%%%d.%ds' : %%d\n",sParamMax, sParamMax);
    char *tmpName = (char *) malloc (sizeof(char) * (sParamMax + 1));
    for (int i = 0; i < _winGlobData[1]; i++) {
      int sParam = _winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i];
      strncpy (tmpName, 
               _winIntParamNameData + _winIntParamIdxNameData[i], 
               min (sParam, sParamMax));
      tmpName[sParam] = '\0';
      
      bftc_printf (fmtName, tmpName, _winIntParamValueData[i]);
      
    }
    free (tmpName);





    MPI_Win_unlock ( _rootRankInGlobalComm, _winGlob);

    
    


    bftc_printf_flush();
  }

  /**
   * \brief Destructor
   *
   */

  CodeProperties::~CodeProperties()
  {
    printf("~CodeProperties %s\n", _name.c_str());
    fflush(stdout);

    if (_winIntParamIdxName != MPI_WIN_NULL) {
      MPI_Win_detach(_winIntParamIdxName, _winIntParamIdxNameData);
      MPI_Win_free(&_winIntParamIdxName);
    }
    if (_winIntParamName != MPI_WIN_NULL) {
      MPI_Win_detach(_winIntParamName, _winIntParamNameData);
      MPI_Win_free(&_winIntParamName);
    }
    if (_winIntParamValue != MPI_WIN_NULL) {
      MPI_Win_detach(_winIntParamValue, _winIntParamValueData);      
      MPI_Win_free(&_winIntParamValue);
    }
    if (_winDoubleParamIdxName != MPI_WIN_NULL) {
      MPI_Win_detach(_winDoubleParamIdxName, _winDoubleParamIdxNameData);
      MPI_Win_free(&_winDoubleParamIdxName);
    }
    if (_winDoubleParamName != MPI_WIN_NULL) {
      MPI_Win_detach(_winDoubleParamName, _winDoubleParamNameData);
      MPI_Win_free(&_winDoubleParamName);
    }
    if (_winDoubleParamValue != MPI_WIN_NULL) {
      MPI_Win_detach(_winDoubleParamValue, _winDoubleParamValueData);
      MPI_Win_free(&_winDoubleParamValue);
    }
    if (_winStrParamIdxName != MPI_WIN_NULL) {
      MPI_Win_detach(_winStrParamIdxName, _winStrParamIdxNameData);
      MPI_Win_free(&_winStrParamIdxName);
    }
    if (_winStrParamName != MPI_WIN_NULL) {
      MPI_Win_detach(_winStrParamName, _winStrParamNameData);
      MPI_Win_free(&_winStrParamName);
    }
    if (_winStrParamIdxValue != MPI_WIN_NULL) {
      MPI_Win_detach(_winStrParamIdxValue, _winStrParamIdxValueData);
      MPI_Win_free(&_winStrParamIdxValue);
    }
    if (_winStrParamValue != MPI_WIN_NULL) {
      MPI_Win_detach(_winStrParamValue, _winStrParamValueData);
      MPI_Win_free(&_winStrParamValue);
    }

    if (_winIntParamIdxNameData != NULL) {
      free (_winIntParamIdxNameData);
      _winIntParamIdxNameData = NULL;
    }

    if (_winIntParamNameData != NULL) {
      free (_winIntParamNameData);
      _winIntParamNameData = NULL;
    }

    if (_winIntParamValueData != NULL) {
      free (_winIntParamValueData);
      _winIntParamValueData = NULL;
    }

    if (_winDoubleParamIdxNameData != NULL) {
      free (_winDoubleParamIdxNameData);
      _winDoubleParamIdxNameData = NULL;
    }

    if (_winDoubleParamNameData != NULL) {
      free (_winDoubleParamNameData);
      _winDoubleParamNameData = NULL;
    }

    if (_winDoubleParamValueData != NULL) {
      free (_winDoubleParamValueData);
      _winDoubleParamValueData = NULL;
    }

    if (_winStrParamIdxNameData != NULL) {
      free (_winStrParamIdxNameData);
      _winStrParamIdxNameData = NULL;
    }

    if (_winStrParamNameData != NULL) {
      free (_winStrParamNameData);
      _winStrParamNameData = NULL;
    }

    if (_winStrParamIdxValueData != NULL) {
      free (_winStrParamIdxValueData);
      _winStrParamIdxValueData = NULL;
    }

    if (_winStrParamValueData != NULL) {
      free (_winStrParamValueData);
      _winStrParamValueData = NULL;
    }
    
    if (_intraComm != MPI_COMM_NULL) {
      MPI_Comm_free(&_intraComm);
    }
    if (_intraGroup != MPI_GROUP_NULL) {
      MPI_Group_free(&_intraGroup);
    }
    if (_intraCoupledGroup != MPI_GROUP_NULL) {
      MPI_Group_free(&_intraCoupledGroup);
    }
  }
}
