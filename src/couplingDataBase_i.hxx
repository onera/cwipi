#ifndef __APPLICATION_COUPLING_DATA_BASE_I_H__
#define __APPLICATION_COUPLING_DATA_BASE_I_H__

#include <cassert>

#include <bft_error.h>

namespace couplings {

  Coupling& CouplingDataBase::getCoupling(const std::string &name) 
  {
    const std::map <std::string, Coupling * >::iterator p = _couplingDataBase.find(name);
    if (p == _couplingDataBase.end())
      bft_error(__FILE__, __LINE__, 0, 
                "'%s' coupling not found \n", name.c_str());
    assert( p->second != NULL);
    return *p->second;
  }

}

#endif
