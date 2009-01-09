#ifndef __COUPLING_PROPERTIES_I_H__
#define __COUPLING_PROPERTIES_I_H__

namespace couplings {

  void  Coupling::set_interpolation_function(couplings_interpolation_fct_t *fct)
  {
    _interpolationFct = fct;
  }

  const int & Coupling::getNNotlocatedPoint() const
  {
    return _nNotLocatedPoint;
  }

  const int *Coupling::getNotlocatedPoint() const
  {
    if (_nNotLocatedPoint == 0)
      return NULL;
    else
      if (_notLocatedPoint == NULL)
        return fvm_locator_get_interior_list(_fvmLocator);
      else
        return _notLocatedPoint;
  }

} // namespace couplings


#endif //__COUPLING_PROPERTIES_I_H__
