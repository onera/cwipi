/*
  This file is part of the CWIPI library.

  Copyright (C) 2012  ONERA

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

#include <spatialInterpIntersection.hxx>

namespace cwipi {
  SpatialInterpIntersection::SpatialInterpIntersection() = default;

  SpatialInterpIntersection::~SpatialInterpIntersection() = default;


  void *SpatialInterpIntersection::interpolate(Field *referenceField) {
      return NULL;
  }

  void SpatialInterpIntersection::weightsCompute() {

  }




  /**
   *
   * \brief Return the number of uncomputed targets
   *
   * \return                Number of uncomputed targets
   *
   */

  int
  SpatialInterpIntersection::nUncomputedTargetsGet(int i_part)  const 
  { 
    return 0;
  }

  /**
   *
   * \brief Return uncomputed targets
   *
   * \return                Uncomputed targets
   *
   */

  const int *
  SpatialInterpIntersection::uncomputedTargetsGet(int i_part)  const 
  {
    return 0;
  }

  /**
   *
   * \brief Return the number of computed targets
   *
   * \return                Number of computed targets
   */

  int
  SpatialInterpIntersection::nComputedTargetsGet(int i_part)  const 
  {
    return 0;
  }

  /**
   *
   * \brief Return computed targets
   *
   *
   * \return                Computed targets
   *
   */

  const int *
  SpatialInterpIntersection::computedTargetsGet(int i_part) const
  {
    return 0;
  }

};