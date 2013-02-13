#ifndef __COUPLING_PROPERTIES_I_H__
#define __COUPLING_PROPERTIES_I_H__
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

#include "locationToLocalMesh.hxx"
#include "locationToDistantMesh.hxx"

namespace cwipi {


  const float *Coupling::distance() const
  {
    return &(_distance[0]);
  }


  void  Coupling::set_interpolation_function(cwipi_interpolation_fct_t fct)
  {
    _interpolationFct = fct;
  }

  void  Coupling::set_interpolation_function_f(void * fct)
  {
    _interpolationFct_f = fct;
  }

  const int * Coupling::getDistantLocation() const
  {
    return _locationToLocalMesh->getLocation();
  }

  const float * Coupling::getDistantDistance() const
  {
    return _locationToLocalMesh->getDistance();
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
    return &_locationToLocalMesh->getBarycentricCoordinatesIndex()[0];
  }

  const double *Coupling::getDistantBarycentricCoordinates() const
  {
    return &_locationToLocalMesh->getBarycentricCoordinates()[0];
  }

  int Coupling::getNDistantPoint() const
  {
    return _locationToLocalMesh->getNLocatedDistantPoint();
  }

  const double *Coupling::getDistantPointCoordinates() const
  {
    return &_locationToLocalMesh->getPointCoordinates()[0];
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
