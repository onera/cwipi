#ifndef __COUPLING_PROPERTIES_I_H__
#define __COUPLING_PROPERTIES_I_H__

#include "locationToLocalMesh.hxx"
#include "locationToDistantMesh.hxx"

namespace cwipi {

  void  Coupling::set_interpolation_function(cwipi_interpolation_fct_t *fct)
  {
    _interpolationFct = fct;
  }

  const int * Coupling::getDistantLocation() const
  {
    return _locationToLocalMesh->getLocation();
  }

  int Coupling::getNNotlocatedPoint() const
  {
    return _locationToDistantMesh->getNUnlocatedPoint();
  }

  const int *Coupling::getNotlocatedPoint() const
  {
    return _locationToDistantMesh->getUnlocatedPoint();
  }

  int Coupling::getNLocatedPoint() const
  {
    return _locationToDistantMesh->getNLocatedPoint();
  }

  const int *Coupling::getLocatedPoint() const
  {
    return _locationToDistantMesh->getLocatedPoint();
  }

  const int *Coupling::getDistantBarycentricCoordinatesIndex() const
  {
    return _locationToLocalMesh->getBarycentricCoordinatesIndex();
  }

  const double *Coupling::getDistantBarycentricCoordinates() const
  {
    return _locationToLocalMesh->getBarycentricCoordinates();
  }

  int Coupling::getNDistantPoint() const
  {
    return _locationToLocalMesh->getNLocatedDistantPoint();
  }

  ///
  /// \brief Set the type of information to be exchanged at the location
  ///
  ///   @param information
  ///

  void Coupling::setInfo(cwipi_located_point_info_t info)
  {
      return _locationToDistantMesh->setInfo(info);
  }

  ///
  /// \brief Get list of number of vertices of containing element
  ///

  const int *Coupling::getElementContainingNVertex() const
  {
    return _locationToDistantMesh->getElementContainingNVertex();
  }

  ///
  /// \brief Get connectivity of element that contains each located point
  ///

  const int *Coupling::getElementContainingVertex() const
  {
    return _locationToDistantMesh->getElementContainingVertex();
  }

  ///
  /// \brief Get vertices coordinates of the element that contains each located point
  ///

  const double *Coupling::getElementContainingVertexCoords() const
  {
    return _locationToDistantMesh->getElementContainingVertexCoords();
  }

  ///
  /// \brief Get barycentric coordinates for the element that contains each located point
  ///

  const double *Coupling::getElementContainingBarycentricCoordinates() const
  {
    return _locationToDistantMesh->getElementContainingBarycentricCoordinates();
  }

  ///
  /// \brief Get MPI rank that contains each located point
  ///

  const int *Coupling::getElementContainingMPIrank() const
  {
    return _locationToDistantMesh->getElementContainingMPIrank();
  }


  ///
  /// \brief Get distant element that contain located point
  ///

  const int *Coupling::getElementContaining() const
  {
    return _locationToDistantMesh->getElementContaining();
  }

} // name space cwipi


#endif //__COUPLING_PROPERTIES_I_H__
