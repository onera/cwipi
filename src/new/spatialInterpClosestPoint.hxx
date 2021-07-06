#ifndef CWP_SPATIALINTERPCLOSESTPOINT_HXX
#define CWP_SPATIALINTERPCLOSESTPOINT_HXX
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
#include "pdm_closest_points.h"

namespace cwipi {
    class SpatialInterpClosestPoint : public SpatialInterp {
    public:
        SpatialInterpClosestPoint();

        ~SpatialInterpClosestPoint() override;

        void spatialInterpWeightsCompute(CWP_Field_exch_t Texch_t) override;

    private:
        void *interpolate(Field *referenceField) override;

        void init(Coupling *coupling, CWP_Dof_location_t pointsCloudLocation, bool slave) override;

        SpatialInterpClosestPoint *_spatial_interp_cpl{};

        PDM_g_num_t *closest_src_gnum{};
        double *closest_src_dstance{};

    protected:
        PDM_closest_point_t *_id_pdm;
    };
}

#endif //CWP_SPATIALINTERPCLOSESTPOINT_HXX