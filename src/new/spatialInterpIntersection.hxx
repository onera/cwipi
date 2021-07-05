#ifndef CWP_SPATIALINTERPINTERSECTION_HXX
#define CWP_SPATIALINTERPINTERSECTION_HXX
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

#include "spatialInterp.hxx"

namespace cwipi {
    class SpatialInterpIntersection : public SpatialInterp {
    public:
        SpatialInterpIntersection();

        ~SpatialInterpIntersection() override;

        void spatialInterpWeightsCompute(CWP_Field_exch_t Texch_t) override;

        void *interpolate(Field *referenceField) override;

        void init(Coupling *coupling, CWP_Dof_location_t pointsCloudLocation, bool slave) override;
    };
}

#endif //CWP_SPATIALINTERPINTERSECTION_HXX
