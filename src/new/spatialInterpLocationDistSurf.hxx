#ifndef CWP_SPATIALINTERPLOCATIONDISTSURF_H
#define CWP_SPATIALINTERPLOCATIONDISTSURF_H
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

#include "spatialInterpLocation.hxx"
#include "pdm_dist_cloud_surf.h"

namespace cwipi {
    class SpatialInterpLocationDistSurf : public SpatialInterpLocation {
    public:
        SpatialInterpLocationDistSurf() = default;

        ~SpatialInterpLocationDistSurf() override = default;

    private:
        void localization_points_cloud_setting() override;

        void localization_surface_setting() override;

        void localization_compute() override;

        void localization_get() override;

        void localization_get_cpl() override;

        void localization_free() override;

    protected:
        PDM_dist_cloud_surf_t *_id_pdm;
    };

/**
 * \endcond
 */

}
#endif //CWP_SPATIALINTERPLOCATIONDISTSURF_H
