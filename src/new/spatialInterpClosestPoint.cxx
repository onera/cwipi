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

#include <pdm_closest_points.h>
#include <spatialInterpClosestPoint.hxx>

namespace cwipi {
    SpatialInterpClosestPoint::SpatialInterpClosestPoint() = default;

    SpatialInterpClosestPoint::~SpatialInterpClosestPoint() = default;

    void SpatialInterpClosestPoint::init(Coupling *coupling, CWP_Dof_location_t pointsCloudLocation, bool slave) {
        // SpatialInterp::init(coupling, pointsCloudLocation, slave);

        // interpolation_time = CWP_INTERP_AT_RECV;

        // CouplingDB *cplDB = _cpl->couplingDBGet();
        // string cplId = coupling->IdGet();
        // if (_both_codes_are_local && !_slave) {
        //     Coupling coupling_cpl = cplDB->couplingGet(*_coupledCodeProperties, cplId);
        //     _spatial_interp_cpl = dynamic_cast<SpatialInterpClosestPoint *>( coupling_cpl.spatialInterpGet(_pointsCloudLocation));
        //     _spatial_interp_cpl->_spatial_interp_cpl = this;
        // }

        // int tmp1; // For dummy communications
        // if (!_both_codes_are_local) {
        //     if (_id < _id_cpl) {
        //         MPI_Bcast(&_nb_part_cpl, 1, MPI_INT, _senderRank, _cplComm);
        //         MPI_Bcast(&tmp1, 1, MPI_INT, _senderRank_cpl, _cplComm);
        //     }
        //     else {
        //         MPI_Bcast(&tmp1, 1, MPI_INT, _senderRank_cpl, _cplComm);
        //         MPI_Bcast(&_nb_part_cpl, 1, MPI_INT, _senderRank, _cplComm);
        //     }
        // }
        // else if (!_slave) {
        //     if (_id < _id_cpl) {
        //         MPI_Bcast(&_nb_part_cpl, 1, MPI_INT, _senderRank, _cplComm);
        //         MPI_Bcast(&(_spatial_interp_cpl->_nb_part_cpl), 1, MPI_INT, _senderRank_cpl, _cplComm);
        //     }
        //     else {
        //         MPI_Bcast(&(_spatial_interp_cpl->_nb_part_cpl), 1, MPI_INT, _senderRank_cpl, _cplComm);
        //         MPI_Bcast(&_nb_part_cpl, 1, MPI_INT, _senderRank, _cplComm);
        //     }
        // }
    }

    void SpatialInterpClosestPoint::spatialInterpWeightsCompute(CWP_Field_exch_t Texch_t) {
        // In case of withOutPart the user provided not null data only on the root rank (senderRank).
        if (!_both_codes_are_local) {
            if (_Texch_t == CWP_FIELD_EXCH_RECV && _pointsCloudLocation == CWP_DOF_LOCATION_USER) user_targets_gnum_compute();
        }
        else {
            if (_Texch_t == CWP_FIELD_EXCH_SEND) {
                _spatial_interp_cpl->_Texch_t = CWP_FIELD_EXCH_RECV;
                if (_pointsCloudLocation == CWP_DOF_LOCATION_USER) _spatial_interp_cpl->user_targets_gnum_compute();
            }
        }

        // Get informations about the local and the coupled meshes
        info_mesh();

        _id_pdm = PDM_closest_points_create(_pdm_cplComm, 5, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);
//        cout << cplComm_rank << ": PDM_closest_point created with id " << _id_pdm << endl;

        PDM_closest_points_n_part_cloud_set(_id_pdm, 1, _nb_part);
//        cout << cplComm_rank << ": PDM_closest_point n_part cloud set" << endl;

        for (int i_part = 0 ; i_part < _nb_part ; i_part++)
            PDM_closest_points_src_cloud_set(_id_pdm, i_part, _n_user_targets[i_part], _coords_user_targets[i_part], _gnum_user_targets[i_part]);
//        cout << cplComm_rank << ": PDM_closest_point src cloud set" << endl;

        for (int i_part = 0 ; i_part < _nb_part ; i_part++)
            PDM_closest_points_tgt_cloud_set(_id_pdm, i_part, _n_target[i_part], _coords_target[i_part], _gnum_target[i_part]);

        PDM_closest_points_compute(_id_pdm);
//        cout << cplComm_rank << ": PDM_closest_point computed" << endl;

        PDM_closest_points_dump_times(_id_pdm);
//        cout << cplComm_rank << ": PDM_closest_point times dumped" << _id_pdm << endl;

        for (int i_part = 0 ; i_part < _nb_part ; i_part++) {
            PDM_closest_points_get(_id_pdm, i_part, &closest_src_gnum, &closest_src_dstance);
        }
//        cout << cplComm_rank << ": PDM_closest_point got " << _id_pdm << endl;

        PDM_closest_points_free (_id_pdm);
//        cout << cplComm_rank << ": PDM_closest_point freed " << _id_pdm << endl;
    }

    void *SpatialInterpClosestPoint::interpolate(Field *referenceField) {
        return NULL;
    }
};