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
#include "pdm_mesh_intersection.h"

namespace cwipi {
  class SpatialInterpIntersection : public SpatialInterp {
  public:
    SpatialInterpIntersection();

    ~SpatialInterpIntersection() override;

    void weightsCompute() override;

    void issend     (Field *referenceField) override;
    void waitIssend (Field *referenceField) override;
    void irecv      (Field *referenceField) override;
    void waitIrecv  (Field *referenceField) override;

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

    void init (Coupling                   *coupling,
               CWP_Dof_location_t          localCodeDofLOcation,
               CWP_Dof_location_t          cplCodeDofLOcation,
               SpatialInterpExchDirection  exchDirection) override;

    SpatialInterpIntersection *_spatial_interp_cpl;

    int         **_tgt_to_src_idx;
    PDM_g_num_t **_tgt_to_src_gnum;
    double      **_tgt_to_src_weight;

  protected:
    PDM_mesh_intersection_t *_id_pdm;
  };

}

#endif //CWP_SPATIALINTERPINTERSECTION_HXX
