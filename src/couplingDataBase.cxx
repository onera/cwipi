
#include "couplingDataBase.hxx"
#include "couplingDataBase_i.hxx"
#include "coupling.hxx"

namespace couplings {
  
  CouplingDataBase::CouplingDataBase()
    :  _couplingDataBase(*new std::map <std::string, Coupling * > ())
  {
  }

  CouplingDataBase::~CouplingDataBase()
  {
    typedef std::map <std::string, Coupling * >::iterator Iterator;
    for (Iterator p = _couplingDataBase.begin(); 
         p != _couplingDataBase.end(); p++) {
      if (p->second != NULL)
        delete p->second;
    }
    _couplingDataBase.clear();
    
    delete &_couplingDataBase;
  }

  void  CouplingDataBase::createCoupling(const std::string &name, 
                                         const ApplicationProperties& localApplicationProperties,
                                         const ApplicationProperties& coupledApplicationProperties,
                                         const int entitiesDim,
                                         const int tolerance,
                                         const couplings_solver_type_t solverType,
                                         const int    outputFrequency,
                                         const char  *outputFormat,
                                         const char  *outputFormatOption)
  {

    Coupling *newCoupling = new Coupling(name, 
                                         localApplicationProperties,
                                         coupledApplicationProperties,
                                         entitiesDim,
                                         tolerance,
                                         solverType,
                                         outputFrequency,
                                         outputFormat,
                                         outputFormatOption);

    std::pair<std::string, Coupling* > 
      newPair(std::string(name), newCoupling);
            
    std::pair<std::map<std::string, Coupling* >::iterator, bool> 
      p = _couplingDataBase.insert(newPair);

    if (!p.second)
      bft_error(__FILE__, __LINE__, 0, 
                "'%s' coupling existing\n", name.c_str());
  }

  void  CouplingDataBase::deleteCoupling(const std::string &name)
  {
    const std::map <std::string, Coupling * >::iterator p = _couplingDataBase.find(name);
    if (p == _couplingDataBase.end())
      bft_error(__FILE__, __LINE__, 0, 
                "'%s' coupling not found \n", name.c_str());

    if (p->second != NULL)
      delete p->second;

    _couplingDataBase.erase(p);
  }
}
