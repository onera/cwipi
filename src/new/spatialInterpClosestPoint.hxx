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

  const int CWP_CLOSEST_POINTS_N_CLOSEST_PTS = 4;

  class SpatialInterpClosestPoint : public SpatialInterp {
  public:
    SpatialInterpClosestPoint();

    ~SpatialInterpClosestPoint() override;

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

        void init (
         Coupling                   *coupling,
         CWP_Dof_location_t          localCodeDofLOcation,
         CWP_Dof_location_t          cplCodeDofLOcation,
         SpatialInterpExchDirection  exchDirection) override;

        SpatialInterpClosestPoint *_spatial_interp_cpl;

        PDM_g_num_t **_closest_src_gnum;
        double      **_closest_src_dist;

        int         **_tgt_in_src_idx;
        PDM_g_num_t **_tgt_in_src_gnum;
        double      **_tgt_in_src_dist;

        // Exchange of src coordinates for least square interpolation
        // ideally exchange only once (same coords for multiple fields)
        int _coordinates_exchanged;
        const double **_send_coord;         /*!< Coordinates of source points */
        double       **_recv_coord;         /*!< Coordinates of target points */
        int            _send_coord_request; /*!< Send request */
        int            _recv_coord_request; /*!< Recv request */
        uint32_t       _send_coord_adler;   /*!< tag MPI from adler code */
        uint32_t       _recv_coord_adler;   /*!< tag MPI from adler code */
        // std::vector <double **> _send_coord;         /*!< Send buffer  (size = n_field) */
        // std::vector <double **> _recv_coord;         /*!< Recv buffer  (size = n_field) */
        // std::vector <int>       _send_coord_request; /*!< Send request (size = n_field) */
        // std::vector <int>       _recv_coord_request; /*!< Recv request (size = n_field) */
        // std::vector <uint32_t>  _send_coord_adler;   /*!< tag MPI from adler code */
        // std::vector <uint32_t>  _recv_coord_adler;   /*!< tag MPI from adler code */


    protected:
        PDM_closest_point_t *_id_pdm;
    };
}

#endif //CWP_SPATIALINTERPCLOSESTPOINT_HXX
