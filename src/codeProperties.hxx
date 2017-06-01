#ifndef __CODE_PROPERTIES_H__
#define __CODE_PROPERTIES_H__
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
#include <cassert>
#include <cstring>

#include <map>
#include <vector>
#include <string>
#include <typeinfo>

#include <bftc_error.h>

#include "cwp.h"

using namespace std;

namespace cwipi {

  /** 
   * \class codeProperties
   *        codeProperties.hxx 
   *        "codeProperties.hxx"
   *
   * \brief Code properties management.
   *
   *  This class manages a code properties :
   *  - Local control parameters,
   *  - Distant control parameters,
   *  - MPI communicators
   *  - . 
   * 
   */

  class CodeProperties {

    friend class CodePropertiesDB;

  public:

    /**
     * \brief Return if the current rank is a coupled rank
     *
     * \return  isCoupledRank
     *
     */

    inline bool
    isCoupledRank() const;


    /**
     * \brief Set if the current rank is coupled
     *
     * \param[in] status   Coupled rank status 
     *
     */

    inline void
    isCoupledRankset 
    (
    bool status
    );

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

    CodeProperties
    (
     string         &name,
     int            id,
     int            rootRank,
     bool            isLocal,
     const MPI_Comm  globalComm
    );

    /**
     * \brief Copy constructor.
     *
     * \param [in]  other        other code properties
     *
     */

    CodeProperties
    (
     const CodeProperties& other
    );

    /**
     * \brief Destructor
     *
     */

    virtual ~CodeProperties();

    /**
     * \brief Get code name
     *
     * \return    Name
     *
     */

    inline const string &
    nameGet() const;

    /**
     * \brief Get root rank in global communicator
     *
     * \return    Root rank
     *
     */

    inline int
    rootRankGet() const;

    /**
     * \brief Get identifier
     *
     * \return    Identifier
     *
     */

    inline int
    idGet() const;

    /**
     * \brief Get global MPI communicator
     *
     * \return    Global MPI communicator
     *
     */

    inline const MPI_Comm &
    globalCommGet() const;

    /**
     * \brief Get MPI intra-communicator
     *
     * \return    MPI intra-communicator
     *
     */

    inline const MPI_Comm &
    intraCommGet() const;

    /**
     * \brief Set MPI intra-communicator
     *
     * \param[in] localComm    MPI intra-communicator
     *
     */

    inline void 
    intraCommSet
    (
     MPI_Comm localComm
    );
    
    /**
     * \brief Get MPI Group in global communicator
     *
     * \return    MPI group
     *
     */

    inline const MPI_Group &
    groupGet() const;


    /**
     * \brief Set MPI intra-communicator
     *
     * \param[in] localComm    MPI intra-communicator
     *
     */

    inline void 
    groupSet
    (
     MPI_Group group
    );

    /**
     * \brief Get a integer control parameter value
     *
     * \param[in] name   Control parameter name  
     *
     */

    inline void
    ctrlParamGet
    (
     const string &name,
     int          *value
    );

    /**
     * \brief Get a double control parameter value
     *
     * \param[in] name   Control parameter name  
     *
     */

    inline void
    ctrlParamGet
    (
     const string &name,
     double       *value
    );

    /**
     * \brief Get a string control parameter value
     *
     * \param[in] name   Control parameter name  
     *
     */

    inline void
    ctrlParamGet
    (
     const string &name,
     string       *value
    );

    /**
     * \brief Set an integer control parameter value
     *
     * \param[in] name    Control parameter name  
     * \param[in] value   Value  
     *
     */

    inline void 
    ctrlParamSet
    (
     const string &name, 
     const int     value
    );

    /**
     * \brief Set a double control parameter value
     *
     * \param[in] name    Control parameter name  
     * \param[in] value   Value  
     *
     */

    inline void 
    ctrlParamSet
    (
     const string &name, 
     const double  value
    );


    /**
     * \brief Set a string control parameter value
     *
     * \param[in] name    Control parameter name  
     * \param[in] value   Value  
     *
     */

    inline void 
    ctrlParamSet
    (
     const string &name, 
     const string  value
    );

    /**
     * \brief Add integer control parameter value
     *
     * \param[in] name    Control parameter name  
     * \param[in] value   Value  
     *
     */

    inline void 
    ctrlParamAdd
    (
     const string &name, 
     const int value
    );

    /**
     * \brief Add a double control parameter value
     *
     * \param[in] name    Control parameter name  
     * \param[in] value   Value  
     *
     */

    inline void 
    ctrlParamAdd
    (
     const string &name, 
     const double  value
    );

    /**
     * \brief Add a string control parameter value
     *
     * \param[in] name    Control parameter name  
     * \param[in] value   Value  
     *
     */

    inline void 
    ctrlParamAdd
    (
     const string &name, 
     const string  value
    );

    /**
     * \brief Cancel a control parameter value
     *
     * \param[in] name    Control parameter name  
     * \param[in] value   Value  
     *
     */

    template<typename T>
    void 
    ctrlParamCancel
    (
     const string &name
    );


    /**
     * \brief Return number of parameters
     *
     * \return Number of parameters
     *
     */

    template<typename T>
    int 
    ctrlParamNGet
    (
    );


    /**
     * \brief Return list of parameters
     *
     * \return Number of parameters
     *
     */

    template<typename T>
    void 
    ctrlParamListGet
    (
    int  *nParam, 
    char ***names
    );


    /**
     * \brief  Is a parameter ?
     *
     * \param[in] name
     *
     * \return  1 : true / 0 : false
     *
     */

    template<typename T>
    int
    ctrlParamIs
    (
     const string &name
    );


    /**
     * \brief  Is a local code ?
     *
     * \return  1 : true / 0 : false
     *
     */

    inline bool
    localCodeIs
    (
    ) const;

    
    /**
     * \brief Dump properties
     *
     */

    void 
    dump();

    
    /**
     * \brief Lock access to the control parameters
     *
     */

    inline void
    paramLock();

    
    /**
     * \brief Unlock access to the control parameters
     *
     */

    inline void
    paramUnLock();

    /**
     * \brief set isLocal code
     *
     */

    inline void
    isLocalSet
    (
    bool status
    );


  private:

    /**
     * \brief Default assignment (unavailable)
     *
     */

    CodeProperties &
    operator=
    (
     const CodeProperties &other
     );

  private:
    string    _name;          /*!< Name */
    int       _id;            /*!< Identifier */ 
    bool      _isLocal;       /*!< Is a local code */
    int       _rootRankInGlobalComm; /*!< Root rank 
                                         *   in MPI global communicator */ 
    MPI_Comm  _globalComm;    /*!< MPI global communicator */
    MPI_Comm  _intraComm;     /*!< MPI intra communicator */
    bool      _isCoupledRank;  /*!< Is a coupled rank */
    MPI_Group _intraGroup;     /*!< MPI group in 
                                                    the global communicator */
    MPI_Group _intraCoupledGroup; /*!< coupled MPI group in 
                                                    the global communicator */

    MPI_Win   _winGlob;        /*!< MPI window to store general parameters informations */
    int       _winGlobData[4]; /*!< \ref _winGlob data (defined only on \ref _rootRankInGlobalComm :
                                   *      - Lock Param Status
                                   *      - Number of int parameters
                                   *      - Number of double parameters
                                   *      - Number of string parameters */
    
    MPI_Win   _winIntParamIdxName; /*!< Window to store indexes of int param names 
                                                 * size = Number of int parameters + 1 */
    MPI_Win   _winIntParamName; /*!< Window to store param names 
                                                 * size = \ref _winIntParamIdxName[Number of int parameter] */
    MPI_Win   _winIntParamValue; /*!< Window to store int param values 
                                                 * size = \re Number of int parameters */
    int      *_winIntParamIdxNameData; /* Data of \ref _winIntParamIdxName window */
    char     *_winIntParamNameData; /* Data of \ref _winIntParamName window */
    int      *_winIntParamValueData; /* Data of \ref _winIntParamValue window */

    MPI_Win   _winDoubleParamIdxName;/*!< Window to store indexes of double param names 
                                                   * size = Number of double parameters + 1 */
    MPI_Win   _winDoubleParamName;/*!< Window to store param names 
                                                * size = \ref _winDoubleParamIdxName[Number of double parameter] */
    MPI_Win   _winDoubleParamValue; /*!< Window to store double param values 
                                                  * size = \ref Number of int parameters */
    int      *_winDoubleParamIdxNameData; /* Data of \ref _winDoubleParamIdxName window */
    char     *_winDoubleParamNameData; /* Data of \ref _winDoubleParamName window */
    double   *_winDoubleParamValueData; /* Data of \ref _winDoubleParamValue window */
                                             
    MPI_Win   _winStrParamIdxName; /*!< Window to store indexes of string param names 
                                                   * size = Number of string parameters + 1 */
    MPI_Win   _winStrParamName; /*!< Window to store param names 
                                                * size = \ref _winDoubleParamIdxName[Number of double parameter] */
    MPI_Win   _winStrParamIdxValue; /*!< Window to store indexes of string param values 
                                                   * size = Number of string parameters + 1 */
    MPI_Win   _winStrParamValue; /*!< Window to store string param values 
                                                  * size = \ref Number of string parameters */
    int      *_winStrParamIdxNameData;  /* Data of \ref _winStrParamIdxName window */
    char     *_winStrParamNameData; /* Data of \ref _winStrParamName window */
    int      *_winStrParamIdxValueData; /* Data of \ref _winStrParamIdxValue window */
    char     *_winStrParamValueData; /* Data of \ref _winStrParamValue window */

  };

  
  /**
   * \brief set isLocal code
   *
   */

  void
  CodeProperties::isLocalSet
  (
  bool status
  )
  {
    _isLocal = status;
  }
  
  /**
   * \brief Lock access to the control parameters
   *
   */

  void
  CodeProperties::paramLock()
  {
    _winGlobData[0] = 1;
  }

  
  /**
   * \brief Unlock access to the control parameters
   *
   */

  void
  CodeProperties::paramUnLock()
  {
    _winGlobData[0] = 0;
  }

  
  /**
   * \brief Return if the current rank is a coupled rank
   *
   * \return  isCoupledRank
   *
   */

  bool
  CodeProperties::isCoupledRank() const
  {
    return _isCoupledRank;
  }


  /**
   * \brief Set if the current rank is coupled
   *
   * \param[in] status   Coupled rank status 
   *
   */

  void
  CodeProperties::isCoupledRankset
  (
  bool status
  )
  {
    _isCoupledRank = status;
  }

  
  /**
   * \brief Get code name
   *
   * \return    Name
   *
   */

  const string &
  CodeProperties::nameGet() const
  {
    return _name;
  }

  /**
   * \brief Get root rank in global communicator
   *
   * \return    Root rank
   *
   */

  int
  CodeProperties::rootRankGet() const
  {
    return _rootRankInGlobalComm;
  }

  /**
   * \brief Get identifier
   *
   * \return    Identifier
   *
   */

  int
  CodeProperties::idGet() const
  {
    return _id;
  }

  /**
   * \brief Get global MPI communicator
   *
   * \return    Global MPI communicator
   *
   */

  const MPI_Comm &
  CodeProperties::globalCommGet() const
  {
    return _globalComm;
  }

  /**
   * \brief Get MPI intra-communicator
   *
   * \return    MPI intra-communicator
   *
   */

  const MPI_Comm &
  CodeProperties::intraCommGet() const
  {
    return _intraComm;
  }

  /**
   * \brief Set MPI intra-communicator
   *
   * \param[in] localComm    MPI intra-communicator
   *
   */

  void 
  CodeProperties::intraCommSet
  (
   MPI_Comm intraComm
  )
  {
    _intraComm = intraComm;
  }

  /**
   * \brief Get MPI Group in global communicator
   *
   * \return    MPI group
   *
   */

  const MPI_Group &
  CodeProperties::groupGet() const
  {
    return _intraCoupledGroup;
  }

  /**
   * \brief Set MPI intra-communicator
   *
   * \param[in] localComm    MPI intra-communicator
   *
   */

  void 
  CodeProperties::groupSet
  (
   MPI_Group group
  )
  {
    _intraCoupledGroup = group;
  }

  /**
   * \brief Get a integer control parameter value
   *
   * \param[in] name   Control parameter name  
   *
   */

  void
  CodeProperties::ctrlParamGet
  (
   const string &name,
   int          *value
  )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    int oldLockStatus   = _winGlobData[0];
    int oldNIntParam    = _winGlobData[1];
    int lockStatus      = oldLockStatus;
    int nIntParam       = oldNIntParam;

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    if (rank != _rootRankInGlobalComm) {

      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm, 0, 4, 
                  MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus); 

      nIntParam    = _winGlobData[1];
      
      if (nIntParam > oldNIntParam) {
        _winIntParamIdxNameData = 
          (int *) realloc(_winIntParamIdxNameData, sizeof(int) * (nIntParam + 1));        
        _winIntParamValueData = 
          (int *) realloc(_winIntParamValueData, sizeof(int) * nIntParam);        
      }
      
      int oldSParamNameData = 0;
      if (oldNIntParam > 0) {
        oldSParamNameData = _winIntParamIdxNameData[oldNIntParam];
      }
      
      MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winIntParamIdxName);
      MPI_Get (_winIntParamIdxNameData, nIntParam + 1, 
               MPI_INT, _rootRankInGlobalComm, 0, nIntParam + 1,
               MPI_INT, _winIntParamIdxName);      
      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamIdxName);
   
      if (_winIntParamIdxNameData[nIntParam] > oldSParamNameData) {
        _winIntParamNameData = (char *) realloc(_winIntParamNameData, 
                               sizeof(char) * _winIntParamIdxNameData[nIntParam]);          
      }
      
      MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winIntParamName);
      MPI_Get (_winIntParamNameData, _winIntParamIdxNameData[nIntParam], 
               MPI_INT, _rootRankInGlobalComm, 0, _winIntParamIdxNameData[nIntParam], 
               MPI_INT, _winIntParamName);      
      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamName);

      MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winIntParamValue);
      MPI_Get (_winIntParamValueData, nIntParam, 
               MPI_INT, _rootRankInGlobalComm, 0, nIntParam,
               MPI_INT, _winIntParamValue);      
      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamValue);
      
    }
    
    else {
      if (lockStatus) {
        bftc_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");      
      }
      
    }

    int sName = name.size();
    int found = 0;
    int i;
    for (i = 0; i < nIntParam; i++) {
      int sParam = _winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i];
      if (sName == sParam) {
       found = !strncmp(name.c_str(), 
               _winIntParamNameData + _winIntParamIdxNameData[i], 
               sName);
      }
      if (found) break;
    }

    if (!found) {
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' Unknown parameter on '%s' code\n", name.c_str(), 
                                                         _name.c_str());      
    }
    
    *value = _winIntParamValueData[i];

    MPI_Win_unlock ( _rootRankInGlobalComm, _winGlob);

  }


  /**
   * \brief Get a double control parameter value
   *
   * \param[in] name   Control parameter name  
   *
   */

  void
  CodeProperties::ctrlParamGet
  (
   const string &name,
   double       *value
  )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    int oldLockStatus   = _winGlobData[0];
    int oldNDoubleParam = _winGlobData[2];
    int lockStatus      = oldLockStatus;
    int nDoubleParam    = oldNDoubleParam;

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    if (rank != _rootRankInGlobalComm) {
      
      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm, 0,
                  4, MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus); 

      nDoubleParam = _winGlobData[2];
      
      if (nDoubleParam > oldNDoubleParam) {
        _winDoubleParamIdxNameData = (int *) realloc(_winDoubleParamIdxNameData,
                                              sizeof(int) * (nDoubleParam + 1));        
        _winDoubleParamValueData = (double *) realloc(_winDoubleParamValueData, 
                                                 sizeof(double) * nDoubleParam);        
      }
      
      int oldSParamNameData = 0;
      if (oldNDoubleParam > 0) {
        oldSParamNameData = _winDoubleParamIdxNameData[oldNDoubleParam];
      }
      
      MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winDoubleParamIdxName);
      MPI_Get (_winDoubleParamIdxNameData, nDoubleParam + 1, 
               MPI_INT, _rootRankInGlobalComm, 0, nDoubleParam + 1, MPI_INT, _winDoubleParamIdxName);      
      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamIdxName);
   
      if (_winDoubleParamIdxNameData[nDoubleParam] > oldSParamNameData) {
        _winDoubleParamNameData = (char *) realloc(_winDoubleParamNameData, 
                                       sizeof(char) * _winDoubleParamIdxNameData[nDoubleParam]);          
      }
      
      MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winDoubleParamName);
      MPI_Get (_winDoubleParamNameData, _winDoubleParamIdxNameData[nDoubleParam], 
               MPI_INT, _rootRankInGlobalComm, 0, _winDoubleParamIdxNameData[nDoubleParam], MPI_INT, _winDoubleParamName);      
      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamName);

      MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winDoubleParamValue);
      MPI_Get (_winDoubleParamValueData, nDoubleParam, 
               MPI_INT, _rootRankInGlobalComm, 0, nDoubleParam, MPI_DOUBLE, _winDoubleParamValue);      
      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamValue);
      
    }
    
    else {
      if (lockStatus) {
        bftc_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");      
      }
      
    }

    int sName = name.size();
    int found = 0;
    int i;
    for (i = 0; i < nDoubleParam; i++) {
      int sParam = _winDoubleParamIdxNameData[i+1] - _winDoubleParamIdxNameData[i];
      if (sName == sParam) {
       found = !strncmp(name.c_str(), 
               _winDoubleParamNameData + _winDoubleParamIdxNameData[i], 
               sName);
      }
      if (found) break;
    }

    if (!found) {
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' Unknown parameter on '%s' code\n", name.c_str(), 
                                                         _name.c_str());      
    }
    
    *value = _winDoubleParamValueData[i];

    MPI_Win_unlock ( _rootRankInGlobalComm, _winGlob);
  }

  /**
   * \brief Get a string control parameter value
   *
   * \param[in] name   Control parameter name  
   *
   */

  void
  CodeProperties::ctrlParamGet
  (
   const string &name,
   string       *value
  )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    int oldLockStatus   = _winGlobData[0];
    int oldNStrParam    = _winGlobData[3];
    int lockStatus      = oldLockStatus;
    int nStrParam       = oldNStrParam;

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    if (rank != _rootRankInGlobalComm) {

      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm,
                  0, 4, MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus); 

      nStrParam    = _winGlobData[3];
      
      if (nStrParam > oldNStrParam) {
        _winStrParamIdxNameData = (int *) realloc(_winStrParamIdxNameData, 
                                          sizeof(int) * (nStrParam + 1));        
        _winStrParamIdxValueData = (int *) realloc(_winStrParamValueData, 
                                        sizeof(int) * (nStrParam + 1));        
      }
      
      int oldSParamNameData = 0;
      int oldSParamValueData = 0;
      if (oldNStrParam > 0) {
        oldSParamNameData = _winStrParamIdxNameData[oldNStrParam];
        oldSParamValueData = _winStrParamIdxValueData[oldNStrParam];
      }
      
      MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winStrParamIdxName);
      MPI_Get (_winStrParamIdxNameData, nStrParam + 1, 
               MPI_INT, _rootRankInGlobalComm, 0, 
               nStrParam + 1, MPI_INT, _winStrParamIdxName);      
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamIdxName);
   
      if (_winStrParamIdxNameData[nStrParam] > oldSParamNameData) {
        _winStrParamNameData = (char *) realloc(_winStrParamNameData, 
                                       sizeof(char) * _winStrParamIdxNameData[nStrParam]);          
      }

      MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winStrParamIdxValue);
      MPI_Get (_winStrParamIdxValueData, nStrParam + 1, 
               MPI_INT, _rootRankInGlobalComm, 0, 
               nStrParam + 1, MPI_INT, _winStrParamIdxValue);      
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamIdxValue);
   
      if (_winStrParamIdxValueData[nStrParam] > oldSParamValueData) {
        _winStrParamValueData = (char *) realloc(_winStrParamValueData, 
                                       sizeof(char) * _winStrParamIdxValueData[nStrParam]);          
      }
      
      MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winStrParamName);
      MPI_Get (_winStrParamNameData, _winStrParamIdxNameData[nStrParam], 
               MPI_INT, _rootRankInGlobalComm, 0, 
               _winStrParamIdxNameData[nStrParam], MPI_INT, _winStrParamName);      
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamName);

      MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winStrParamValue);
      MPI_Get (_winStrParamValueData, _winStrParamIdxValueData[nStrParam], 
               MPI_INT, _rootRankInGlobalComm, 0, 
               _winStrParamIdxValueData[nStrParam], MPI_INT, _winStrParamValue);      
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamValue);
      
    }
    
    else {
      if (lockStatus) {
        bftc_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");      
      }
      
    }

    int sName = name.size();
    int found = 0;
    int i;
    for (i = 0; i < nStrParam; i++) {
      int sParam = _winStrParamIdxNameData[i+1] - _winStrParamIdxNameData[i];
      if (sName == sParam) {
       found = !strncmp(name.c_str(), 
               _winStrParamNameData + _winStrParamIdxNameData[i], 
               sName);
      }
      if (found) break;
    }

    if (!found) {
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' Unknown parameter on '%s' code\n", name.c_str(), 
                                                         _name.c_str());      
    }

    int sValue = _winStrParamIdxValueData[i+1] - _winStrParamIdxValueData[i];
    value->resize(sValue);
    strncpy ((char *) value->c_str(),
             _winStrParamValueData + _winStrParamIdxValueData[i], 
             sValue);
    ((char *)value->c_str())[sValue] = '\0';
    
    MPI_Win_unlock ( _rootRankInGlobalComm, _winGlob);
  }

  /**
   * \brief Set an integer control parameter value
   *
   * \param[in] name    Control parameter name  
   * \param[in] value   Value  
   *
   */

  void 
  CodeProperties::ctrlParamSet
  (
   const string &name, 
   const int     value
  )
  {
    if (!_isLocal) {
      bftc_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Set a distant code parameter is not allowed\n", 
                                                   _name.c_str());      
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winIntParamValue);
    
      int sName = name.size();
      int found = 0;
      int i;
      
      int nIntParam    = _winGlobData[1];
      
      for (i = 0; i < nIntParam; i++) {
        int sParam = _winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i];
        if (sName = sParam) {
         found = !strncmp(name.c_str(), 
                 _winIntParamNameData + _winIntParamIdxNameData[i], 
                 sName);
        }
        if (found) break;
      }

      if (!found) {
        bftc_error(__FILE__, __LINE__, 0,
                   "'%s' Unknown parameter on '%s' code\n", name.c_str(), 
                                                           _name.c_str());      
      }

      _winIntParamValueData[i] = value;
    
      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }
  }

  /**
   * \brief Set a double control parameter value
   *
   * \param[in] name    Control parameter name  
   * \param[in] value   Value  
   *
   */

  void 
  CodeProperties::ctrlParamSet
  (
   const string &name, 
   const double  value
  )
  {
    if (!_isLocal) {
      bftc_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Set a distant code parameter is not allowed\n", 
                                                   _name.c_str());      
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winDoubleParamValue);
    
      int sName = name.size();
      int found = 0;
      int i;
      
      int nDoubleParam = _winGlobData[2];
      
      for (i = 0; i < nDoubleParam; i++) {
        int sParam = _winDoubleParamIdxNameData[i+1] - _winDoubleParamIdxNameData[i];
        if (sName = sParam) {
         found = !strncmp(name.c_str(), 
                 _winDoubleParamNameData + _winDoubleParamIdxNameData[i], 
                 sName);
        }
        if (found) break;
      }

      if (!found) {
        bftc_error(__FILE__, __LINE__, 0,
                   "'%s' Unknown parameter on '%s' code\n", name.c_str(), 
                                                           _name.c_str());      
      }

      _winDoubleParamValueData[i] = value;
    
      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }
  }

  /**
   * \brief Set a string control parameter value
   *
   * \param[in] name    Control parameter name  
   * \param[in] value   Value  
   *
   */

  void 
  CodeProperties::ctrlParamSet
  (
   const string &name, 
   const string  value
  )
  {
    if (!_isLocal) {
      bftc_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Set a distant code parameter is not allowed\n", 
                                                   _name.c_str());      
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamValue);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamIdxValue);
    
      int sName = name.size();
      int found = 0;
      int i;
      
      int nStrParam    = _winGlobData[3];
      
      for (i = 0; i < nStrParam; i++) {
        int sParam = _winStrParamIdxNameData[i+1] - _winStrParamIdxNameData[i];
        if (sName = sParam) {
         found = !strncmp(name.c_str(), 
                 _winStrParamNameData + _winStrParamIdxNameData[i], 
                 sName);
        }
        if (found) break;
      }

      if (!found) {
        bftc_error(__FILE__, __LINE__, 0,
                   "'%s' Unknown parameter on '%s' code\n", name.c_str(), 
                                                           _name.c_str());      
      }

      int sValue = _winStrParamIdxValueData[i+1] - _winStrParamIdxValueData[i];
      int gap = value.size() - sValue;
      if (gap > 0) {
        MPI_Win_detach(_winStrParamValue, _winStrParamValueData);
        _winStrParamValueData = (char *) realloc (_winStrParamValueData, 
   (_winStrParamIdxValueData[nStrParam] + value.size() - sValue) * sizeof(char));
        MPI_Aint swin = _winStrParamIdxValueData[nStrParam] + value.size() - sValue;
        MPI_Win_attach (_winStrParamValue, _winStrParamValueData, swin);
      }
      
      if (gap != 0) {
        if (gap > 0) {
          for (int i1 = _winStrParamIdxValueData[nStrParam] - 1; i1 >= _winStrParamIdxValueData[i+1]; i1--) {
            _winStrParamValueData[i1+gap] = _winStrParamValueData[i1]; 
          }
        }
        else {
          for (int i1 = _winStrParamIdxValueData[i+1]; i1 < _winStrParamIdxValueData[nStrParam]; i1++) {
            _winStrParamValueData[i1+gap] = _winStrParamValueData[i1]; 
          }          
        }
        for (int i1 = i+1; i1 < nStrParam; i1++) {
          _winStrParamIdxValueData[i1] += gap; 
        }                
      }

      if (gap < 0) {
       MPI_Win_detach(_winStrParamValue, _winStrParamValueData);
       _winStrParamValueData = (char *) realloc (_winStrParamValueData, 
                                _winStrParamIdxValueData[nStrParam] * sizeof(char));
       MPI_Win_attach (_winStrParamValue, _winStrParamValueData, _winStrParamIdxValueData[nStrParam]);
      }

      strncpy(_winStrParamValueData + _winStrParamIdxValueData[i], 
              value.c_str(), value.size()); 
      
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamIdxValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }
  }

  /**
   * \brief Add an integer control parameter value
   *
   * \param[in] name    Control parameter name  
   * \param[in] value   Value  
   *
   */

  void 
  CodeProperties::ctrlParamAdd
  (
   const string &name, 
   const int     value
  )
  {
    if (!_isLocal) {
      bftc_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Add a distant code parameter is not allowed\n", 
                                                   _name.c_str());      
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winIntParamValue);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winIntParamName);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winIntParamIdxName);
    
      int sName = name.size();
      int found = 0;
      int i;
      
      int nIntParam    = _winGlobData[1];
      
      for (i = 0; i < nIntParam; i++) {
        int sParam = _winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i];
        if (sName = sParam) {
         found = !strncmp(name.c_str(), 
                 _winIntParamNameData + _winIntParamIdxNameData[i], 
                 sName);
        }
        if (found) break;
      }

      if (found) {
        bftc_error(__FILE__, __LINE__, 0,
                   "'%s' is already a parameter off '%s' code\n", name.c_str(), 
                                                                 _name.c_str());      
      }

      if (nIntParam > 0) {
        MPI_Win_detach(_winIntParamValue, _winIntParamValueData);
        MPI_Win_detach(_winIntParamIdxName, _winIntParamIdxNameData);
        MPI_Win_detach(_winIntParamName, _winIntParamNameData);
      }
      
      nIntParam += 1;

      int sWin1 = (nIntParam + 1) * sizeof(int);
      int sWinName1 = (_winIntParamIdxNameData[nIntParam-1] + name.size()) * sizeof(char);
      _winIntParamValueData = (int *) realloc (_winIntParamValueData, sWin1);
      _winIntParamIdxNameData = (int *) realloc (_winIntParamIdxNameData, sWin1);
      _winIntParamNameData = (char *) realloc (_winIntParamNameData, sWinName1);

      MPI_Aint sWin = (MPI_Aint) sWin1;
      MPI_Aint sWinName = (MPI_Aint) sWinName1;

      MPI_Win_attach (_winIntParamValue, _winIntParamValueData, sWin);
      MPI_Win_attach (_winIntParamIdxName, _winIntParamIdxNameData, sWin);
      MPI_Win_attach (_winIntParamName, _winIntParamNameData, sWinName);
      
      _winIntParamIdxNameData[nIntParam] = _winIntParamIdxNameData[nIntParam-1] + name.size();

      strncpy(_winIntParamNameData + _winIntParamIdxNameData[nIntParam-1], 
              name.c_str(), name.size()); 

      _winIntParamValueData[nIntParam-1] = value;
      
      _winGlobData[1] = nIntParam;

      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamIdxName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);

    }      
  }


  /**
   * \brief Add a double control parameter value
   *
   * \param[in] name    Control parameter name  
   * \param[in] value   Value  
   *
   */

  void 
  CodeProperties::ctrlParamAdd
  (
   const string &name, 
   const double  value
  )
  {
    if (!_isLocal) {
      bftc_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Add a distant code parameter is not allowed\n", 
                                                   _name.c_str());      
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winDoubleParamValue);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winDoubleParamName);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winDoubleParamIdxName);
    
      int sName = name.size();
      int found = 0;
      int i;
      
      int nDoubleParam = _winGlobData[2];
      
      for (i = 0; i < nDoubleParam; i++) {
        int sParam = _winDoubleParamIdxNameData[i+1] - _winDoubleParamIdxNameData[i];
        if (sName == sParam) {
         found = !strncmp(name.c_str(), 
                 _winDoubleParamNameData + _winDoubleParamIdxNameData[i], 
                 sName);
        }
        if (found) break;
      }

      if (found) {
        bftc_error(__FILE__, __LINE__, 0,
                   "'%s' is already a parameter off '%s' code\n", name.c_str(), 
                                                                 _name.c_str());      
      }

      if (nDoubleParam > 0) {
        MPI_Win_detach(_winDoubleParamValue, _winDoubleParamValueData);
        MPI_Win_detach(_winDoubleParamIdxName, _winDoubleParamIdxNameData);
        MPI_Win_detach(_winDoubleParamName, _winDoubleParamNameData);
      }
      
      nDoubleParam += 1;

      int sWin1 = (nDoubleParam + 1) * sizeof(int);
      int sWinName1 = (_winDoubleParamIdxNameData[nDoubleParam-1] + name.size()) * sizeof(char);
      _winDoubleParamValueData = (double *) realloc (_winDoubleParamValueData, sWin1);
      _winDoubleParamIdxNameData = (int *) realloc (_winDoubleParamIdxNameData, sWin1);
      _winDoubleParamNameData = (char *) realloc (_winDoubleParamNameData, sWinName1);

      MPI_Aint sWin = (MPI_Aint) sWin1;
      MPI_Aint sWinName = (MPI_Aint) sWinName1;

      MPI_Win_attach (_winDoubleParamValue, _winDoubleParamValueData, sWin);
      MPI_Win_attach (_winDoubleParamIdxName, _winDoubleParamIdxNameData, sWin);
      MPI_Win_attach (_winDoubleParamName, _winDoubleParamNameData, sWinName);
      
      _winDoubleParamIdxNameData[nDoubleParam] = _winDoubleParamIdxNameData[nDoubleParam-1] + name.size();

      strncpy(_winDoubleParamNameData + _winDoubleParamIdxNameData[nDoubleParam-1], 
              name.c_str(), name.size()); 

      _winDoubleParamValueData[nDoubleParam-1] = value;
      
      _winGlobData[2] = nDoubleParam;

      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamIdxName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }
  }

  /**
   * \brief Add a string control parameter value
   *
   * \param[in] name    Control parameter name  
   * \param[in] value   Value  
   *
   */

  void 
  CodeProperties::ctrlParamAdd
  (
   const string &name, 
   const string  value
  )
  {
    
    if (!_isLocal) {
      bftc_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Add a distant code parameter is not allowed\n", 
                                                   _name.c_str());      
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamValue);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamIdxValue);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamName);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamIdxName);
    
      int sName = name.size();
      int found = 0;
      int i;
      
      int nStrParam    = _winGlobData[3];
      
      for (i = 0; i < nStrParam; i++) {
        int sParam = _winStrParamIdxNameData[i+1] - _winStrParamIdxNameData[i];
        if (sName = sParam) {
         found = !strncmp(name.c_str(), 
                 _winStrParamNameData + _winStrParamIdxNameData[i], 
                 sName);
        }
        if (found) break;
      }

      if (found) {
        bftc_error(__FILE__, __LINE__, 0,
                   "'%s' is already a parameter off '%s' code\n", name.c_str(), 
                                                                 _name.c_str());      
      }
      if (nStrParam > 0) {
        MPI_Win_detach(_winStrParamValue, _winStrParamIdxValueData);
        MPI_Win_detach(_winStrParamValue, _winStrParamValueData);
        MPI_Win_detach(_winStrParamIdxName, _winStrParamIdxNameData);
        MPI_Win_detach(_winStrParamName, _winStrParamNameData);
      }

      nStrParam += 1;

      int sWin1 = (nStrParam + 1) * sizeof(int);
      _winStrParamIdxValueData = (int *) realloc (_winStrParamValueData, sWin1);
      _winStrParamIdxNameData = (int *) realloc (_winStrParamIdxNameData, sWin1);

      int sWinName1 = (_winStrParamIdxNameData[nStrParam-1] + name.size()) * sizeof(char);
      _winStrParamNameData = (char *) realloc (_winStrParamNameData, sWinName1);

      int sWinValue1 = (_winStrParamIdxValueData[nStrParam-1] + value.size()) * sizeof(char);
      _winStrParamValueData = (char *) realloc (_winStrParamValueData, sWinValue1);

      MPI_Aint sWin = (MPI_Aint) sWin1;
      MPI_Win_attach (_winStrParamIdxValue, _winStrParamIdxValueData, sWin);
      _winStrParamIdxValueData[nStrParam] = _winStrParamIdxValueData[nStrParam-1] + value.size();

      MPI_Win_attach (_winStrParamIdxName, _winStrParamIdxNameData, sWin);
      _winStrParamIdxNameData[nStrParam] = _winStrParamIdxNameData[nStrParam-1] + name.size();

      MPI_Aint sWinName = (MPI_Aint) sWinName1;
      MPI_Win_attach (_winStrParamName, _winStrParamNameData, sWinName);
      
      MPI_Aint sWinValue = (MPI_Aint) sWinValue1;
      MPI_Win_attach (_winStrParamValue, _winStrParamValueData, sWinValue);

      strncpy(_winStrParamNameData + _winStrParamIdxNameData[nStrParam-1], 
              name.c_str(), name.size()); 

      strncpy(_winStrParamValueData + _winStrParamIdxValueData[nStrParam-1], 
              value.c_str(), value.size()); 
    
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamIdxName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamIdxValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
      
      _winGlobData[3] = nStrParam;

    }
  }

  /**
   * \brief Cancel a control parameter value
   *
   * \param[in] name    Control parameter name  
   * \param[in] value   Value  
   *
   */

  template<typename T>
  void
  CodeProperties::ctrlParamCancel
  (
   const string &name
  )
  {
    if (!_isLocal) {
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' is a distant code. Delete a distant"
                 " code parameter is not allowed\n", 
                                                   _name.c_str());      
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
    
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);

      int nIntParam    = _winGlobData[1];
      int nDoubleParam = _winGlobData[2];
      int nStrParam    = _winGlobData[3];

      MPI_Win  *winTypeParamIdxValue = NULL;
      MPI_Win  *winTypeParamValue = NULL;
      MPI_Win  *winTypeParamIdxName = NULL;
      MPI_Win  *winTypeParamName = NULL;

      int nTypeParam;
      int  *winTypeParamIdxValueData = NULL;
      T *winTypeParamValueData = NULL;
      int  *winTypeParamIdxNameData = NULL;
      char *winTypeParamNameData = NULL;

      if (typeid(T) == typeid(string)) {
        nTypeParam               = nStrParam;
        winTypeParamIdxValue     = &_winStrParamIdxValue;
        winTypeParamValue        = &_winStrParamValue;
        winTypeParamIdxName      = &_winStrParamIdxName;
        winTypeParamName         = &_winStrParamName;
        winTypeParamIdxValueData = _winStrParamIdxValueData;
        winTypeParamValueData    = (T *) _winStrParamValueData;
        winTypeParamIdxNameData  = _winStrParamIdxNameData; 
        winTypeParamNameData     = _winStrParamNameData; 
      }
      else if (typeid(T) == typeid(int)) {
        nTypeParam              = nIntParam;
        winTypeParamValue       = &_winIntParamValue;
        winTypeParamIdxName     = &_winIntParamIdxName;
        winTypeParamName        = &_winIntParamName;
        winTypeParamValueData   = (T *) _winIntParamValueData;
        winTypeParamIdxNameData = _winIntParamIdxNameData; 
        winTypeParamNameData    = _winIntParamNameData; 
      }
      else if (typeid(T) == typeid(double)) {
        nTypeParam              = nDoubleParam;
        winTypeParamValue       = &_winDoubleParamValue;
        winTypeParamIdxName     = &_winDoubleParamIdxName;
        winTypeParamName        = &_winDoubleParamName;
        winTypeParamValueData   = (T *) _winDoubleParamValueData;
        winTypeParamIdxNameData = _winDoubleParamIdxNameData; 
        winTypeParamNameData    = _winDoubleParamNameData; 
      }
      else {
        bftc_error(__FILE__, __LINE__, 0,
                  "Type not taken into account \n");
      }
      
      int i;
      int found;
      
      for (i = 0; i < nTypeParam; i++) {
        int sParam = winTypeParamIdxNameData[i+1] - winTypeParamIdxNameData[i];
        if (name.size() == sParam) {
         found = !strncmp(name.c_str(), 
                 winTypeParamNameData + winTypeParamIdxNameData[i], 
                 name.size());
        }
        if (found) break;
      }

      if (found) {
        if (winTypeParamIdxValueData != NULL) {
          MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, 
                        *winTypeParamIdxValue);
          MPI_Win_detach(*winTypeParamIdxValue, winTypeParamIdxValueData);
        }
        MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, 
                      *winTypeParamValue);
        MPI_Win_detach(*winTypeParamValue, winTypeParamValueData);
        MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, 
                      *winTypeParamName);
        MPI_Win_detach(*winTypeParamName, winTypeParamNameData);
        MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, 
                      *winTypeParamIdxName);
        MPI_Win_detach(*winTypeParamIdxName, winTypeParamIdxNameData);
        
        if (winTypeParamIdxValueData != NULL) {
          int gap = winTypeParamIdxValueData[i+1] - winTypeParamIdxValueData[i];
          for (int i1 = winTypeParamIdxValueData[i]; i1 < winTypeParamIdxValueData[nTypeParam] - gap; i1++) {
            winTypeParamValueData[i1] = winTypeParamValueData[i1+gap]; 
          }
          for (int i1 = i; i1 < nTypeParam - 1; i1++) {
            winTypeParamIdxValueData[i1] = winTypeParamIdxValueData[i1+1] - gap; 
          }                  
        }
        else {
          for (int i1 = i; i1 < nTypeParam - 1; i1++) {
            winTypeParamValueData[i1] = winTypeParamValueData[i1+1]; 
          }        
        }
        
        int gap2 = winTypeParamIdxNameData[i+1] - winTypeParamIdxNameData[i];
        for (int i1 = winTypeParamIdxNameData[i]; i1 < winTypeParamIdxNameData[nTypeParam] - gap2; i1++) {
          winTypeParamNameData[i1] = winTypeParamNameData[i1+gap2]; 
        }        
        
        for (int i1 = i; i1 < nTypeParam - 1; i1++) {
          winTypeParamIdxNameData[i1] = winTypeParamIdxNameData[i1+1] - gap2; 
        }                  
        
        nTypeParam += -1; 

        int sWin1 = (nTypeParam + 1) * sizeof(int);
        if (winTypeParamIdxValueData != NULL) {
          winTypeParamIdxValueData = (int *) realloc (winTypeParamValueData, sWin1);
        }
        winTypeParamIdxNameData = (int *) realloc (winTypeParamIdxNameData, sWin1);

        int sWinName1 = (winTypeParamIdxNameData[nStrParam] + name.size()) * sizeof(char);
        winTypeParamNameData = (char *) realloc (winTypeParamNameData, sWinName1);

        int sWinValue1;
        if (winTypeParamIdxValueData != NULL) {
          sWinValue1 = winTypeParamIdxValueData[nStrParam] * sizeof(T);
        }
        else {
          sWinValue1 = nTypeParam * sizeof(T);          
        }
        
        winTypeParamValueData = (T *) realloc (_winStrParamValueData, sWinValue1);

        MPI_Aint sWin = (MPI_Aint) sWin1;
        if (winTypeParamIdxValueData != NULL) {
          MPI_Win_attach (*winTypeParamIdxValue, winTypeParamIdxValueData, sWin);
        }
        MPI_Win_attach (*winTypeParamIdxName, winTypeParamIdxNameData, sWin);
        winTypeParamIdxNameData[nStrParam+1] = winTypeParamIdxNameData[nStrParam] + name.size();

        MPI_Aint sWinName = (MPI_Aint) sWinName1;
        MPI_Win_attach (*winTypeParamName, winTypeParamNameData, sWinName);

        MPI_Aint sWinValue = (MPI_Aint) sWinValue1;
        MPI_Win_attach (*winTypeParamValue, winTypeParamValueData, sWinValue);
        
        if (typeid(T) == typeid(string)) {
          nStrParam = nTypeParam;
          _winStrParamIdxValueData = winTypeParamIdxValueData;
          _winStrParamValueData    = (char *) winTypeParamValueData;
          _winStrParamIdxNameData  = winTypeParamIdxNameData; 
          _winStrParamNameData     = winTypeParamNameData; 
        }
        else if (typeid(T) == typeid(int)) {
          nIntParam = nTypeParam;
          _winIntParamValueData    = (int *) winTypeParamValueData;
          _winIntParamIdxNameData  = winTypeParamIdxNameData; 
          _winIntParamNameData     = winTypeParamNameData; 
        }
        else if (typeid(T) == typeid(double)) {
          nDoubleParam = nTypeParam;
          _winDoubleParamValueData    = (double *) winTypeParamValueData;
          _winDoubleParamIdxNameData  = winTypeParamIdxNameData; 
          _winDoubleParamNameData     = winTypeParamNameData; 
        }
        else {
          bftc_error(__FILE__, __LINE__, 0,
                    "Type not taken into account \n");
        }

        if (winTypeParamIdxValueData != NULL) {
          MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamIdxValue);
        }
        MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamValue);
        MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamName);
        MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamIdxName);

      } 
      
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }

  }

  /**
   * \brief Return number of parameters
   *
   * \return Number of parameters
   *
   */

  template<typename T>
  int 
  CodeProperties::ctrlParamNGet
  (
   )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    int oldLockStatus   = _winGlobData[0];
    int oldNIntParam    = _winGlobData[1];
    int oldNDoubleParam = _winGlobData[2];
    int oldNStrParam    = _winGlobData[3];
    int lockStatus      = oldLockStatus;
    int nIntParam       = oldNIntParam;
    int nDoubleParam    = oldNDoubleParam;
    int nStrParam       = oldNStrParam;

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    if (rank != _rootRankInGlobalComm) {

      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm,
                  0, 4, MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus); 

      nIntParam    = _winGlobData[1];
      nDoubleParam = _winGlobData[2];
      nStrParam    = _winGlobData[3];
    }
    
    else {
      if (lockStatus) {
        bftc_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");      
      }
      
    }

    MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);

    int nParam;
    if (typeid(T) == typeid(string)) { 
      nParam = nStrParam;
    }
    else if (typeid(T) == typeid(int)) {
      nParam = nIntParam;
    }
    else if (typeid(T) == typeid(double)) {
      nParam = nDoubleParam;
    }
    else {
      bftc_error(__FILE__, __LINE__, 0,
                "Type not taken into account \n");
    }
    
    return nParam;
  }


  /**
   * \brief Return list of parameters
   *
   * \return Number of parameters
   *
   */
  
  template<typename T>
  void  
  CodeProperties::ctrlParamListGet
  (
   int  *nParam, 
   char ***names
  )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    int oldLockStatus   = _winGlobData[0];
    int oldNIntParam    = _winGlobData[1];
    int oldNDoubleParam = _winGlobData[2];
    int oldNStrParam    = _winGlobData[3];
    int lockStatus      = oldLockStatus;
    int nIntParam       = oldNIntParam;
    int nDoubleParam    = oldNDoubleParam;
    int nStrParam       = oldNStrParam;

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    if (rank != _rootRankInGlobalComm) {

      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm,
                  0, 4, MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus); 

      nIntParam    = _winGlobData[1];
      nDoubleParam = _winGlobData[2];
      nStrParam    = _winGlobData[3];

    }
    
    else {
      if (lockStatus) {
        bftc_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");      
      }
      
    }

    MPI_Win  *winTypeParamIdxName = NULL;
    MPI_Win  *winTypeParamName = NULL;

    int nTypeParam;
    int  *winTypeParamIdxNameData = NULL;
    char *winTypeParamNameData = NULL;

    if (typeid(T) == typeid(string)) {
      nTypeParam               = nStrParam;
      winTypeParamIdxName      = &_winStrParamIdxName;
      winTypeParamName         = &_winStrParamName;
      winTypeParamIdxNameData  = _winStrParamIdxNameData; 
      winTypeParamNameData     = _winStrParamNameData; 
    }
    else if (typeid(T) == typeid(int)) {
      nTypeParam              = nIntParam;
      winTypeParamIdxName     = &_winIntParamIdxName;
      winTypeParamName        = &_winIntParamName;
      winTypeParamIdxNameData = _winIntParamIdxNameData; 
      winTypeParamNameData    = _winIntParamNameData; 
    }
    else if (typeid(T) == typeid(double)) {
      nTypeParam              = nDoubleParam;
      winTypeParamIdxName     = &_winDoubleParamIdxName;
      winTypeParamName        = &_winDoubleParamName;
      winTypeParamIdxNameData = _winDoubleParamIdxNameData; 
      winTypeParamNameData    = _winDoubleParamNameData; 
    }
    else {
      bftc_error(__FILE__, __LINE__, 0,
                "Type not taken into account \n");
    }

    MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, 
                  *winTypeParamName);
    MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, 
                  *winTypeParamIdxName);

    *names = (char **) malloc(sizeof(char *) * nTypeParam);

    for (int i = 0; i < nTypeParam; i++) {
      int sParam = winTypeParamIdxNameData[i+1] - winTypeParamIdxNameData[i];
      (*names)[i] = (char *) malloc(sizeof(char) * (sParam + 1));
      char *curName = (*names)[i];

      strncpy(curName,
              winTypeParamNameData + winTypeParamIdxNameData[i], 
              sParam);

      curName[sParam] = '\0';

    }

    *nParam = nTypeParam;

    MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamName);
    MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamIdxName);
    MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);

  }


  /**
   * \brief  Is a parameter ?
   *
   * \param[in] name
   *
   * \return  1 : true / 0 : false
   *
   */

  template<typename T>
  int
  CodeProperties::ctrlParamIs
  (
   const string &name 
  )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    int oldLockStatus   = _winGlobData[0];
    int oldNIntParam    = _winGlobData[1];
    int oldNDoubleParam = _winGlobData[2];
    int oldNStrParam    = _winGlobData[3];
    int lockStatus      = oldLockStatus;
    int nIntParam       = oldNIntParam;
    int nDoubleParam    = oldNDoubleParam;
    int nStrParam       = oldNStrParam;

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    if (rank != _rootRankInGlobalComm) {

      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm,
                  0, 4, MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus); 

      nIntParam    = _winGlobData[1];
      nDoubleParam = _winGlobData[2];
      nStrParam    = _winGlobData[3];

    }
    
    else {
      if (lockStatus) {
        bftc_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");      
      }
      
    }

    MPI_Win  *winTypeParamIdxName = NULL;
    MPI_Win  *winTypeParamName = NULL;

    int  *winTypeParamIdxNameData = NULL;
    char *winTypeParamNameData = NULL;
    
    int nTypeParam;

    if (typeid(T) == typeid(string)) {
      nTypeParam = nStrParam;
      winTypeParamIdxName      = &_winStrParamIdxName;
      winTypeParamName         = &_winStrParamName;
      winTypeParamIdxNameData  = _winStrParamIdxNameData; 
      winTypeParamNameData     = _winStrParamNameData; 
    }
    else if (typeid(T) == typeid(int)) {
      nTypeParam = nIntParam;
      winTypeParamIdxName     = &_winIntParamIdxName;
      winTypeParamName        = &_winIntParamName;
      winTypeParamIdxNameData = _winIntParamIdxNameData; 
      winTypeParamNameData    = _winIntParamNameData; 
    }
    else if (typeid(T) == typeid(double)) {
      nTypeParam = nDoubleParam;
      winTypeParamIdxName     = &_winDoubleParamIdxName;
      winTypeParamName        = &_winDoubleParamName;
      winTypeParamIdxNameData = _winDoubleParamIdxNameData; 
      winTypeParamNameData    = _winDoubleParamNameData; 
    }
    else {
      bftc_error(__FILE__, __LINE__, 0,
                "Type not taken into account \n");
    }

    MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, 
                  *winTypeParamName);
    MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, 
                  *winTypeParamIdxName);

    int sName = name.size();
    int found = 0;
    for (int i = 0; i < nTypeParam; i++) {
      int sParam = winTypeParamIdxNameData[i+1] - winTypeParamIdxNameData[i];
      if (sName == sParam) {
       found = !strncmp(name.c_str(), 
               winTypeParamNameData + winTypeParamIdxNameData[i], 
               sName);
      }
      if (found) break;
    }

    MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamName);
    MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamIdxName);
    MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);

    return found;
  }

  
  /**
   * \brief  Is a local code ?
   *
   * \return  1 : true / 0 : false
   *
   */

  bool
  CodeProperties::localCodeIs
  (
  ) const
  {
    return _isLocal;    
  }

}


#endif /* __CODE_PROPERTIES_H__ */

