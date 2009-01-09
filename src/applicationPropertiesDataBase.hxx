#ifndef __APPLICATION_PROPERTIES_DATA_BASE_H__
#define __APPLICATION_PROPERTIES_DATA_BASE_H__

#define MAX(a,b) ((a) > (b) ? (a) : (b))

//Bug mpich2
#define MPICH_IGNORE_CXX_SEEK 1
#include <mpi.h>
#include <map>

#include <bft_printf.h>

#include "singleton.hpp"

namespace couplings {
  
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
    void init(const char* name, const MPI_Comm &globalComm, MPI_Comm &localComm);

    inline void setPrintfProxy(bft_printf_proxy_t *const callBackPrintf);

    // Access to local MPI properties (synchronisation with a second application)

    inline MPI_Comm &getLocalComm() const;

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

    inline void setLocalIntControlParameter(const std::string &name, const int value);

    inline void setLocalDoubleControlParameter(const std::string &name, const double value);

    inline const int &getLocalIntControlParameter(const std::string &name);

    inline const double &getLocalDoubleControlParameter(const std::string &name);

    inline void eraseLocalIntControlParameter(const std::string &name);

    inline void eraseLocalDoubleControlParameter(const std::string &name);

    // Merge local parameters with distant parameters (synchronisation with a second application)

    void mergeParameters(const std::string &applicationName);

    // Access to distant control parameters

    inline const int &getDistantIntControlParameter(const std::string &applicationName, const std::string &name) const;

    inline const double &getDistantDoubleControlParameter(const std::string &applicationName, const std::string &name) const;

    void dump();

  private:
    void _mergeIntParameters(const std::string &applicationName); 
    void _mergeDoubleParameters(const std::string &applicationName); 

  private:
    std::map <std::string, ApplicationProperties * > & _distantApplicationPropertiesDataBase;
    ApplicationProperties * _localApplicationProperties;

  };
}

#endif /* __APPLICATION_PROPERTIES_DATA_BASE_H__ */
