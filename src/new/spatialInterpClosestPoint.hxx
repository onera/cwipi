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


    /**
     *
     * \brief Return the number of uncomputed targets
     *
     * \return                Number of uncomputed targets
     *
     */

    int
    nUncomputedTargetsGet(int i_part) const override; 

    /**
     *
     * \brief Return uncomputed targets
     *
     * \return                Uncomputed targets
     *
     */

    const int *
    uncomputedTargetsGet(int i_part) const override;

    /**
     *
     * \brief Return the number of computed targets
     *
     * \return                Number of computed targets
     */

    int
    nComputedTargetsGet(int i_part) const override;

    /**
     *
     * \brief Return computed targets
     *
     *
     * \return                Computed targets
     *
     */

    const int *
    computedTargetsGet(int i_part) const override;




    private:
        void *interpolate(Field *referenceField) override;

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

        SpatialInterpClosestPoint *_spatial_interp_cpl;

        PDM_g_num_t *closest_src_gnum;
        double *closest_src_dstance;

    protected:
        PDM_closest_point_t *_id_pdm;
    };
}

#endif //CWP_SPATIALINTERPCLOSESTPOINT_HXX