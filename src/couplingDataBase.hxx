#ifndef __APPLICATION_COUPLING_DATA_BASE_H__
#define __APPLICATION_COUPLING_DATA_BASE_H__

#include <map>
#include <string>

#include <mpi.h>

#include <bft_printf.h>

#include "singleton.hpp"
#include "cwipi.h"

namespace cwipi {
  class Coupling;
  class ApplicationProperties;
  
  class CouplingDataBase : public Singleton <CouplingDataBase>
  {
    friend class Singleton <CouplingDataBase>;
    friend class Coupling;
    
  public:
    void createCoupling(const std::string &name, 
                        const cwipi_coupling_type_t couplingType,
                        const ApplicationProperties& localApplicationProperties,
                        const ApplicationProperties& coupledApplicationProperties,
                        const int entitiesDim,
                        const double tolerance,
                        const cwipi_solver_type_t solverType,
                        const int    outputFrequency,
                        const char  *outputFormat,
                        const char  *outputFormatOption);
    
    void deleteCoupling(const std::string &name);
    
    inline Coupling& getCoupling(const std::string &name);
    
    
  private:
    CouplingDataBase();
    CouplingDataBase(const CouplingDataBase &other);
    CouplingDataBase & operator=(const CouplingDataBase &other);
    virtual ~CouplingDataBase();
    
  private:
    std::map <std::string, Coupling * > & _couplingDataBase;
    MPI_Comm  _fvmComm;

  };
}

#endif
