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

#include "pdm_closest_points.h"
#include "spatialInterp.hxx"

namespace cwipi {
  class SpatialInterpClosestPoint : public SpatialInterp {
  public:
    SpatialInterpClosestPoint();

    ~SpatialInterpClosestPoint() override;

    void weightsCompute() override;

    private:
        void interpolate(Field *referenceField, double **buffer) override;

        /**
          *
          * \brief Initialization of the SpatialInterp object.
          *
          * \param [in] coupling            Pointer the coupling object.
          * \param [in] pointsCloudLocation Location of the cloud of points.
          * \param [in] coupling            Pointer the coupling object.
          *
          */

        // void init (
        //   Coupling *coupling, 
        //   CWP_Dof_location_t pointsCloudLocation,
        //   CWP_Dof_location_t cplCodeDofLOcation) override;
        void init (
         Coupling                   *coupling,
         CWP_Dof_location_t          localCodeDofLOcation,
         CWP_Dof_location_t          cplCodeDofLOcation,
         SpatialInterpExchDirection  exchDirection) override;

        SpatialInterpClosestPoint *_spatial_interp_cpl;

        PDM_g_num_t *closest_src_gnum;    // do we keep this? (single part?)
        double      *closest_src_dstance; // do we keep this? (single part?)


        int         **_src_to_tgt_idx;
        PDM_g_num_t **_src_to_tgt;
        double      **_src_to_tgt_dist;

    protected:
        PDM_closest_point_t *_id_pdm;
    };
}

#endif //CWP_SPATIALINTERPCLOSESTPOINT_HXX
