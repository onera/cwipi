#ifndef __COUPLING_PROPERTIES_I_H__
#define __COUPLING_PROPERTIES_I_H__

namespace cwipi {

  void  Coupling::set_interpolation_function(cwipi_interpolation_fct_t *fct)
  {
    _interpolationFct = fct;
  }

  const int * Coupling::getDistantLocation() const
  {
    if (_toLocate)
      bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    return _location;
  }

  const int & Coupling::getNNotlocatedPoint() const
  {
    if (_toLocate)
      bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    return _nNotLocatedPoint;
  }

  const int *Coupling::getNotlocatedPoint() const
  {
    if (_toLocate)
      bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    if (_nNotLocatedPoint == 0)
      return NULL;
    else
      if (_notLocatedPoint == NULL)
        return fvm_locator_get_interior_list(_fvmLocator);
      else
        return _notLocatedPoint;
  }

  int Coupling::getNLocatedPoint() const
  {
    if (_toLocate)
      bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    return _nPointsToLocate - _nNotLocatedPoint;
  }

  const int *Coupling::getLocatedPoint() const
  {
    if (_toLocate)
      bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    return _locatedPoint;
  }

  const int *Coupling::getDistantBarycentricCoordinatesIndex() const
  {
    if (_toLocate)
      bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    return _barycentricCoordinatesIndex;
  }

  const double *Coupling::getDistantBarycentricCoordinates() const
  {
    if (_toLocate)
      bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    return _barycentricCoordinates;
  }

  inline int Coupling::getNDistantPoint() const
  {
    if (_toLocate)
      bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    return _nDistantpoint;
  }


} // namespace cwipi


#endif //__COUPLING_PROPERTIES_I_H__
