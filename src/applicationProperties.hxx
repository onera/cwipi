#ifndef __APPLICATION_PROPERTIES_H__
#define __APPLICATION_PROPERTIES_H__
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

#define MPICH_IGNORE_CXX_SEEK 1

#include <cassert>

#include <map>
#include <string>

#include <mpi.h>
#include <bft_error.h>
#include "applicationPropertiesDataBase.hxx"

namespace cwipi {

  class ApplicationProperties {

    friend void ApplicationPropertiesDataBase::_mergeIntParameters(const std::string &applicationName);
    friend void ApplicationPropertiesDataBase::_mergeDoubleParameters(const std::string &applicationName);
    friend void ApplicationPropertiesDataBase::_mergeStringParameters(const std::string &applicationName);

  public:
    ApplicationProperties(std::string &name,
                          const MPI_Comm globalComm);

    ApplicationProperties(const ApplicationProperties& name);

    virtual ~ApplicationProperties();

    inline const std::string &getName() const;

    inline const MPI_Comm &getGlobalComm() const;

    inline const MPI_Comm &getLocalComm() const;

    inline void setLocalComm(MPI_Comm localComm);

    inline const int &getBeginningRank() const;

    inline const int &getEndRank() const;

    inline void setBeginningRank(const int value);

    inline void setEndRank(const int value);

    inline const int &getIntControlParameter(const std::string &name) const;

    inline const double &getDoubleControlParameter(const std::string &name) const;

    inline const std::string &getStringControlParameter(const std::string &name) const;

    inline void setIntControlParameter(const std::string &name, const int value);

    inline void setDoubleControlParameter(const std::string &name, const double value);

    inline void setStringControlParameter(const std::string &name, const std::string value);

    inline void addIntControlParameter(const std::string &name, const int initialValue);

    inline void addDoubleControlParameter(const std::string &name, const double initialValue);

    inline void addStringControlParameter(const std::string &name, const std::string initialValue);

    inline void eraseIntControlParameter(const std::string &name);

    inline void eraseDoubleControlParameter(const std::string &name);

    inline void eraseStringControlParameter(const std::string &name);

    void dump();

  private:
    ApplicationProperties();
    ApplicationProperties &operator=(const ApplicationProperties &other);

  private:
    std::string  _name;
    MPI_Comm  _globalComm;
    MPI_Comm  _localComm;
    int _beginningRank;
    int _endRank;
    std::map <std::string, int> & _intControlParameters;
    std::map <std::string, double> & _doubleControlParameters;
    std::map <std::string, std::string> &_stringControlParameters;
  };

  const std::string &ApplicationProperties::getName() const
  {
    return _name;
  }

  const MPI_Comm &ApplicationProperties::getGlobalComm() const
  {
    return _globalComm;
  }

  const MPI_Comm &ApplicationProperties::getLocalComm() const
  {
    return _localComm;
  }

  void ApplicationProperties::setLocalComm(MPI_Comm localComm)
  {
    _localComm = localComm;
  }

  const int &ApplicationProperties::getBeginningRank() const
  {
    return _beginningRank;
  }

  const int &ApplicationProperties::getEndRank() const
  {
    return _endRank;
  }

  void ApplicationProperties::setBeginningRank(const int value)
  {
    _beginningRank = value;
  }

  void ApplicationProperties::setEndRank(const int value)
  {
    _endRank = value;
  }

  const int &ApplicationProperties::getIntControlParameter(const std::string &name) const
  {
    const std::map <std::string, int>::iterator p = _intControlParameters.find(name);
    if (p == _intControlParameters.end())
      bft_error(__FILE__, __LINE__, 0,
                "'%s' int control parameter not found \n", name.c_str());
    return p->second;
  }

  const double &ApplicationProperties::getDoubleControlParameter(const std::string &name) const
  {
    const std::map <std::string, double>::iterator p = _doubleControlParameters.find(name);
    if (p == _doubleControlParameters.end())
      bft_error(__FILE__, __LINE__, 0,
                "'%s' double control parameter not found \n", name.c_str());
    return p->second;
  }

  const std::string &ApplicationProperties::getStringControlParameter(const std::string &name) const
  {
    const std::map <std::string, std::string>::iterator p = _stringControlParameters.find(name);
    if (p == _stringControlParameters.end())
      bft_error(__FILE__, __LINE__, 0,
                "'%s' string control parameter not found \n", name.c_str());
    return p->second;
  }

  void ApplicationProperties::setIntControlParameter(const std::string &name, const int value)
  {
    std::map <std::string, int>::iterator p = _intControlParameters.find(name);
    if (p == _intControlParameters.end())
      bft_error(__FILE__, __LINE__, 0,
                "'%s' int control parameter not found \n", name.c_str());
    p->second = value;
  }

  void ApplicationProperties::setDoubleControlParameter(const std::string &name, const double value)
  {
    std::map <std::string, double>::iterator p = _doubleControlParameters.find(name);
    if (p != _doubleControlParameters.end())
      p->second = value;
    else
      bft_error(__FILE__, __LINE__, 0,
                "'%s' double control parameter not found \n", name.c_str());
  }

  void ApplicationProperties::setStringControlParameter(const std::string &name, const std::string value)
  {
    std::map <std::string, std::string>::iterator p = _stringControlParameters.find(name);
    if (p != _stringControlParameters.end())
      p->second = value;
    else
      bft_error(__FILE__, __LINE__, 0,
                "'%s' string control parameter not found \n", name.c_str());
  }

  void ApplicationProperties::addIntControlParameter(const std::string &name, const int initialValue)
  {
    std::pair<std::string, int> parameter(name, initialValue);
    std::pair<std::map<std::string, int>::iterator, bool> p = _intControlParameters.insert(parameter);
    if (!p.second)
      bft_error(__FILE__, __LINE__, 0,
                "'%s' existing int control parameter \n", name.c_str());
  }

  void ApplicationProperties::addDoubleControlParameter(const std::string &name, const double initialValue)
  {
    std::pair<std::string, double> parameter(name, initialValue);
    std::pair<std::map<std::string, double>::iterator, bool> p = _doubleControlParameters.insert(parameter);
    if (!p.second)
      bft_error(__FILE__, __LINE__, 0,
                "'%s' existing double control parameter \n", name.c_str());
  }

  void ApplicationProperties::addStringControlParameter(const std::string &name, const std::string initialValue)
  {
    std::pair<std::string, std::string> parameter(name, initialValue);
    std::pair<std::map<std::string, std::string>::iterator, bool> p = _stringControlParameters.insert(parameter);
    if (!p.second)
      bft_error(__FILE__, __LINE__, 0,
                "'%s' existing string control parameter \n", name.c_str());
  }
  void ApplicationProperties::eraseIntControlParameter(const std::string &name)
  {
    std::map <std::string, int>::iterator p = _intControlParameters.find(name);
    if (p != _intControlParameters.end())
      _intControlParameters.erase(p);
    else
      bft_error(__FILE__, __LINE__, 0,
                "'%s' int control parameter not found \n", name.c_str());
  }

  void ApplicationProperties::eraseDoubleControlParameter(const std::string &name)
  {
    std::map <std::string, double>::iterator p = _doubleControlParameters.find(name);
    if (p != _doubleControlParameters.end())
      _doubleControlParameters.erase(p);
    else
      bft_error(__FILE__, __LINE__, 0,
                "'%s' double control parameter not found \n", name.c_str());
  }

  void ApplicationProperties::eraseStringControlParameter(const std::string &name)
  {
    std::map <std::string, std::string>::iterator p = _stringControlParameters.find(name);
    if (p != _stringControlParameters.end())
      _stringControlParameters.erase(p);
    else
      bft_error(__FILE__, __LINE__, 0,
                "'%s' string control parameter not found \n", name.c_str());
  }
}

#endif /* __APPLICATION_PROPERTIES_H__ */

