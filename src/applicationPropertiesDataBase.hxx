#ifndef __APPLICATION_PROPERTIES_DATA_BASE_H__
#define __APPLICATION_PROPERTIES_DATA_BASE_H__
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

#define MAX(a,b) ((a) > (b) ? (a) : (b))

//Bug mpich2
#define MPICH_IGNORE_CXX_SEEK 1
#include <mpi.h>
#include <map>

#include <bft_printf.h>

#include "singleton.hpp"

namespace cwipi {

  class ApplicationProperties;

  class ApplicationPropertiesDataBase : public Singleton <ApplicationPropertiesDataBase>
  {
    friend class Singleton <ApplicationPropertiesDataBase>;
    friend class ApplicationProperties;

  private:
    ApplicationPropertiesDataBase();
    ApplicationPropertiesDataBase(const ApplicationPropertiesDataBase &other);
    ApplicationPropertiesDataBase & operator=(const ApplicationPropertiesDataBase &other);
    virtual ~ApplicationPropertiesDataBase();

  public:
    MPI_Comm init(const char* name, const MPI_Comm globalComm);

    inline void setPrintfProxy(bft::bft_printf_proxy_t *const callBackPrintf);

    // Access to local MPI properties (synchronisation with a second application)

    inline const MPI_Comm &getLocalComm() const;

    inline const MPI_Comm &getGlobalComm() const;

    inline const int &getBeginningRank() const;

    inline const int &getEndRank() const;

    // Access to distant MPI properties (synchronisation with a second application)

    //inline const MPI_Comm &getDistantLocalComm(const std::string &applicationName);

    inline const int &getDistantBeginningRank(const std::string &applicationName)const;

    inline const int &getDistantEndRank(const std::string &applicationName) const;

    inline const ApplicationProperties &getDistantApplicationProperties(const std::string &applicationName) const;

    inline const ApplicationProperties &getLocalApplicationProperties() const ;

    // Access to local control parameters

    inline void addLocalIntControlParameter(const std::string &name, const int value);

    inline void addLocalDoubleControlParameter(const std::string &name, const double value);

    inline void addLocalStringControlParameter(const std::string &name, const std::string value);

    inline void setLocalIntControlParameter(const std::string &name, const int value);

    inline void setLocalDoubleControlParameter(const std::string &name, const double value);

    inline void setLocalStringControlParameter(const std::string &name, const std::string value);

    inline const int &getLocalIntControlParameter(const std::string &name);

    inline const double &getLocalDoubleControlParameter(const std::string &name);

    inline const std::string &getLocalStringControlParameter(const std::string &name);

    inline void eraseLocalIntControlParameter(const std::string &name);

    inline void eraseLocalDoubleControlParameter(const std::string &name);

    inline void eraseLocalStringControlParameter(const std::string &name);

    // Merge local parameters with distant parameters (synchronisation with a second application)

    void mergeParameters(const std::string &applicationName);

    // Access to distant control parameters

    inline const int &getDistantIntControlParameter(const std::string &applicationName, const std::string &name) const;

    inline const double &getDistantDoubleControlParameter(const std::string &applicationName, const std::string &name) const;

    inline const std::string &getDistantStringControlParameter(const std::string &applicationName, const std::string &name) const;

    void dump();

  private:
    void _mergeIntParameters(const std::string &applicationName);
    void _mergeDoubleParameters(const std::string &applicationName);
    void _mergeStringParameters(const std::string &applicationName);

  private:
    std::map <std::string, ApplicationProperties * > & _distantApplicationPropertiesDataBase;
    ApplicationProperties * _localApplicationProperties;
  };
}

#endif /* __APPLICATION_PROPERTIES_DATA_BASE_H__ */
