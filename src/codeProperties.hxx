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
     * \brief Constructor.
     *
     * \param [in]  name         Current code name
     * \param [in]  globalComm   MPI communicator containing all processes 
     *                           of all codes
     *
     */

    CodeProperties
    (
     string         &name,
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
     * \brief Get the first code rank into the global communicator
     *
     * \return    First code rank into the global communicator
     *
     */

    inline const int &
    firstRankGet() const;

    /**
     * \brief Get the last code rank into the global communicator
     *
     * \return   First code rank into the global communicator
     *
     */

    inline const int &
    lastRankGet() const;

    /**
     * \brief Set the first code rank into the global communicator
     *
     * \param[in] value    First rank number
     *
     */

    inline void 
    firstRankSet
    (
     const int value
    );

    /**
     * \brief Set the last code rank into the global communicator
     *
     * \param[in] value    Last rank number
     *
     */

    inline void 
    lastRankSet
    (
     const int value
    );

    /**
     * \brief Get int control parameters map
     *
     */

    inline void
    ctrlParamGet
    (
     map <string, int> ** params
    ) const;
    
    /**
     * \brief Get double control parameters map
     *
     */

    inline void
    ctrlParamGet
    (
     map <string, double> ** params
    ) const;

    /**
     * \brief Get string control parameters map
     *
     */

    inline void
    ctrlParamGet
    (
     map <string, string> ** params
    ) const;

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
     int          **value
    ) const;

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
     double       **value
    ) const;

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
     string       **value
    ) const;

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
    ) const;


    /**
     * \brief Return list of parameters
     *
     * \return Number of parameters
     *
     */

    template<typename T>
    char ** 
    ctrlParamListGet
    (
    ) const;


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
    ) const;


    /**
     * \brief Dump properties
     *
     */

    void 
    dump();

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
    string                 _name;          /*!< Name */
    MPI_Comm               _globalComm;    /*!< MPI global communicator */
    MPI_Comm               _intraComm;     /*!< MPI intra communicator */
    int                    _firstRank;     /*!< Current code first rank into
                                                MPI global communicator */
    int                    _lastRank;      /*!< Current code last rank into
                                                MPI global communicator */
    map <string, int>    & _intCtrlParam;  /*!< Integer control parameters */ 
    map <string, double> & _dblCtrlParam;  /*!< Double control parameters */
    map <string, string> & _strCtrlParam;  /*!< String control parameters */
  };

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
   * \brief Get the first code rank into the global communicator
   *
   * \return    First code rank into the global communicator
   *
   */

  const int &
  CodeProperties::firstRankGet() const
  {
    return _firstRank;
  }

  /**
   * \brief Get the last code rank into the global communicator
   *
   * \return   First code rank into the global communicator
   *
   */

  const int &
  CodeProperties::lastRankGet() const
  {
    return _lastRank;
  }

  /**
   * \brief Set the first code rank into the global communicator
   *
   * \param[in] value    First rank number
   *
   */

  void 
  CodeProperties::firstRankSet
  (
   const int value
  )
  {
    _firstRank = value;
  }

  /**
   * \brief Set the last code rank into the global communicator
   *
   * \param[in] value    Last rank number
   *
   */

  void 
  CodeProperties::lastRankSet
  (
   const int value
  )
  {
    _lastRank = value;
  }

  /**
   * \brief Get int control parameters map
   *
   */
  
  void
  CodeProperties::ctrlParamGet
  (
   map <string, int> **params
  ) const
  {
    *params = &_intCtrlParam;
  }
  
  /**
   * \brief Get double control parameters map
   *
   */

  void
  CodeProperties::ctrlParamGet
  (
   map <string, double> ** params
  ) const
  {
    *params = &_dblCtrlParam;
  }

  /**
   * \brief Get string control parameters map
   *
   */

  void
  CodeProperties::ctrlParamGet
  (
   map <string, string> ** params
  ) const
  {
    *params = &_strCtrlParam;
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
   int          **value
  ) const
  {
    map <string, int>::iterator p = _intCtrlParam.find(name);
    if (p == _intCtrlParam.end())
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' int control parameter not found \n", name.c_str());
    *value = &(p->second);
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
   double       **value
  ) const
  {
    map <string, double>::iterator p = _dblCtrlParam.find(name);
    if (p == _dblCtrlParam.end())
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' double control parameter not found \n", name.c_str());
    *value = &(p->second);
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
   string       **value
  ) const
  {
    map <string, string>::iterator p = _strCtrlParam.find(name);
    if (p == _strCtrlParam.end())
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' str control parameter not found \n", name.c_str());
    *value = &(p->second);
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
    map <string, int>::iterator p = _intCtrlParam.find(name);
    if (p == _intCtrlParam.end())
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' int control parameter not found \n", name.c_str());
    p->second = value;
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
    map <string, double>::iterator p = _dblCtrlParam.find(name);
    if (p == _dblCtrlParam.end())
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' dbl control parameter not found \n", name.c_str());
    p->second = value;
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
    map <string, string>::iterator p = _strCtrlParam.find(name);
    if (p == _strCtrlParam.end())
      bftc_error(__FILE__, __LINE__, 0,
                   "'%s' string control parameter not found \n", name.c_str());
    p->second = value;
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
    pair<string, int> parameter(name, value);
    pair<map<string, int>::iterator, bool> p = 
      _intCtrlParam.insert(parameter);
    if (!p.second)
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' is already an integer control parameter\n", name.c_str());
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
    pair<string, double> parameter(name, value);
    pair<map<string, double>::iterator, bool> p = 
      _dblCtrlParam.insert(parameter);
    if (!p.second)
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' is already a double control parameter\n", name.c_str());
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
    pair<string, string> parameter(name, value);
    pair<map<string, string>::iterator, bool> p = 
      _strCtrlParam.insert(parameter);
    if (!p.second)
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' is already a string control parameter\n", name.c_str());
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
    if (typeid(T) == typeid(string)) { 
      map <string, string>::iterator p = _strCtrlParam.find(name);
      if (p != _strCtrlParam.end())
        _strCtrlParam.erase(p);
      else
        bftc_error(__FILE__, __LINE__, 0,
                   "'%s' string control parameter not found \n", name.c_str());
    }
    else if (typeid(T) == typeid(int)) {
      map <string, int>::iterator p = _intCtrlParam.find(name);
      if (p != _intCtrlParam.end())
        _intCtrlParam.erase(p);
      else
        bftc_error(__FILE__, __LINE__, 0,
                   "'%s' int control parameter not found \n", name.c_str());
    }
    else if (typeid(T) == typeid(double)) {
      map <string, double>::iterator p = _dblCtrlParam.find(name);
      if (p != _dblCtrlParam.end())
        _dblCtrlParam.erase(p);
      else
        bftc_error(__FILE__, __LINE__, 0,
                   "'%s' double control parameter not found \n", name.c_str());
    }
    else
      bftc_error(__FILE__, __LINE__, 0,
                "Type not taken into account \n");
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
   )  const
  {
    if (typeid(T) == typeid(string)) { 
      return _strCtrlParam.size();
    }
    else if (typeid(T) == typeid(int)) {
      return _intCtrlParam.size();
    }
    else if (typeid(T) == typeid(double)) {
      return _dblCtrlParam.size();
    }
    else {
      bftc_error(__FILE__, __LINE__, 0,
                "Type not taken into account \n");
    }
    return 0;
  }


  /**
   * \brief Return list of parameters
   *
   * \return Number of parameters
   *
   */
  
  template<typename T>
  char ** 
  CodeProperties::ctrlParamListGet
  (
  )  const
  {
    char **names = new char * [ctrlParamNGet <T> ()];
    if (typeid(T) == typeid(string)) { 
      int i = 0;
      for(typename map < string, string >::iterator p = _strCtrlParam.begin(); 
          p != _strCtrlParam.end(); 
          ++p) {
        names[i] = new char [p->first.size() + 1];
        strcpy(names[i], p->first.c_str());
        i += 1;
      }
    } 
    else if (typeid(T) == typeid(int)) {
      int i = 0;
      for(typename map < string, int >::iterator p = _intCtrlParam.begin(); 
          p != _intCtrlParam.end(); 
          ++p) {
        names[i] = new char [p->first.size() + 1];
        strcpy(names[i], p->first.c_str());
        i += 1;
      }
    }
    else if (typeid(T) == typeid(double)) {
      int i = 0;
      for(typename map < string, double >::iterator p = _dblCtrlParam.begin(); 
          p != _dblCtrlParam.end(); 
          ++p) {
        names[i] = new char [p->first.size() + 1];
        strcpy(names[i], p->first.c_str());
        i += 1;
      }
    }
    else
      bftc_error(__FILE__, __LINE__, 0,
                 "Type not taken into account \n");
    return names;
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
  )  const
  {
    if (typeid(T) == typeid(string)) { 
      typename map < string, string >::iterator  p;
      typename map < string, string >::iterator  p_end;
      p     = _strCtrlParam.find(name);
      p_end = _strCtrlParam.end();
      if (p != p_end)
        return 1;
      else
        return 0;
    }
    else if (typeid(T) == typeid(int)) {
      typename map < string, int >::iterator  p;
      typename map < string, int >::iterator  p_end;
      p     = _intCtrlParam.find(name);
      p_end = _intCtrlParam.end();
      if (p != p_end)
        return 1;
      else
        return 0;
    }
    else if (typeid(T) == typeid(double)) {
      typename map < string, double >::iterator  p;
      typename map < string, double >::iterator  p_end;
      p     = _dblCtrlParam.find(name);
      p_end = _dblCtrlParam.end();
      if (p != p_end)
        return 1;
      else
        return 0;
    }
    else
      bftc_error(__FILE__, __LINE__, 0,
                "Type not taken into account \n");
    return 0;
  }
}


#endif /* __CODE_PROPERTIES_H__ */

