#ifndef __APPLICATION_PROPERTIES_DATA_BASE_I_H__
#define __APPLICATION_PROPERTIES_DATA_BASE_I_H__

#include <cassert>

#include "applicationProperties.hxx"


namespace couplings {

  inline void ApplicationPropertiesDataBase::setPrintfProxy(bft_printf_proxy_t *const callBackPrintf)
  {
    if (callBackPrintf != NULL)
      bft_printf_proxy_set(callBackPrintf);
  }
  
  // Access to local MPI properties (synchronisation with a second application)
  
  inline const MPI_Comm &ApplicationPropertiesDataBase::getLocalComm() const
  {
    return _localApplicationProperties->getLocalComm();
  }
  
  inline const MPI_Comm &ApplicationPropertiesDataBase::getGlobalComm() const
  {
    return _localApplicationProperties->getGlobalComm();
  }
  
  inline const int &ApplicationPropertiesDataBase::getBeginningRank() const
  {
    return _localApplicationProperties->getBeginningRank();
  }
 
  inline const int &ApplicationPropertiesDataBase::getEndRank() const
  {
    return _localApplicationProperties->getEndRank();
  }
  
  // Access to distant MPI properties (synchronisation with a second application)

//   inline const MPI_Comm &ApplicationPropertiesDataBase::getDistantLocalComm(const std::string &applicationName)
//   {
//     std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
//     if (p != _distantApplicationPropertiesDataBase.end())
//       return p->second->getLocalComm();
//     else
//       bft_error(__FILE__, __LINE__, 0, 
//                 "'%s' application not found \n", applicationName.c_str());
//   }
  
  inline const int &ApplicationPropertiesDataBase::getDistantBeginningRank(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bft_error(__FILE__, __LINE__, 0, 
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getBeginningRank();
  }
  
  inline const int &ApplicationPropertiesDataBase::getDistantEndRank(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bft_error(__FILE__, __LINE__, 0, 
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getEndRank();
  }
  
  // Access to local control parameters

  inline void ApplicationPropertiesDataBase::addLocalIntControlParameter(const std::string &name, const int initialValue)
  {
    _localApplicationProperties->addIntControlParameter(name, initialValue);
  }
  
  inline void ApplicationPropertiesDataBase::addLocalDoubleControlParameter(const std::string &name, const double initialValue)
  {
    _localApplicationProperties->addDoubleControlParameter(name, initialValue);
  }
  
  inline void ApplicationPropertiesDataBase::setLocalIntControlParameter(const std::string &name, const int value)
  {
    _localApplicationProperties->setIntControlParameter(name, value);
  }
  
  inline void ApplicationPropertiesDataBase::setLocalDoubleControlParameter(const std::string &name, const double value)
  {
    _localApplicationProperties->setDoubleControlParameter(name, value);
  }
  
  inline const int &ApplicationPropertiesDataBase::getLocalIntControlParameter(const std::string &name)
  {
    return _localApplicationProperties->getIntControlParameter(name);
  }
  
  inline const double &ApplicationPropertiesDataBase::getLocalDoubleControlParameter(const std::string &name)
  {
    return _localApplicationProperties->getDoubleControlParameter(name);
  }
  
  inline void ApplicationPropertiesDataBase::eraseLocalIntControlParameter(const std::string &name)
  {
    _localApplicationProperties->eraseIntControlParameter(name);
  }
  
  inline void ApplicationPropertiesDataBase::eraseLocalDoubleControlParameter(const std::string &name)
  {
    _localApplicationProperties->eraseDoubleControlParameter(name);
  }
  
  // Access to distant control parameters
  
  inline const int &ApplicationPropertiesDataBase::getDistantIntControlParameter(const std::string &applicationName, const std::string &name) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bft_error(__FILE__, __LINE__, 0, 
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getIntControlParameter(name);
  }
  
  inline const ApplicationProperties &ApplicationPropertiesDataBase::getDistantApplicationProperties(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bft_error(__FILE__, __LINE__, 0, 
                "'%s' application not found \n", applicationName.c_str());
    assert(p->second != NULL);
    return *(p->second);
  }
  
  inline const ApplicationProperties &ApplicationPropertiesDataBase::getLocalApplicationProperties() const 
  {
    assert(_localApplicationProperties != NULL);
    return *_localApplicationProperties;
  }

  inline const double &ApplicationPropertiesDataBase::getDistantDoubleControlParameter(const std::string &applicationName, const std::string &name) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bft_error(__FILE__, __LINE__, 0, 
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getDoubleControlParameter(name);
  }
  
} // namespace couplings

#endif /* __APPLICATION_PROPERTIES_DATA_BASE_H__ */
