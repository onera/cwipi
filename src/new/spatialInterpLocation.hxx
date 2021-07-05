#ifndef SPATIALINTERPLOCATION_H_
#define SPATIALINTERPLOCATION_H_
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

/**
 * \cond
 */

namespace cwipi {
    class SpatialInterpLocation : public SpatialInterp {
    public:

        /**
          *
          * \brief SpatialInterp location constructor.
          *
          */
        SpatialInterpLocation();

        /**
          *
          * \brief SpatialInterp location destructor.
          *
          */
        ~SpatialInterpLocation() override;

        /**
          *
          * \brief Compute of the spatial interpolation weights. Localization and communication
          *        tree building.
          *
          * \param [in] Texch_t    Type of exchange (sending or reception).
          *
          */
        void spatialInterpWeightsCompute(CWP_Field_exch_t Texch_t) override;

        SpatialInterpLocation *_spatial_interp_cpl{};

        // Localization data
        double **_distance{};                 // Distance to the closest element surface by partition
        double **_projected{};                // Projected point coordinates (on the closest element surface)
        CWP_g_num_t **_closest_elt_gnum{};    // Closest element global numbering

    protected:
        /**
          *
          * \brief Interpolation of a point cloud on a reference field.
          *
          * \param [in]   referenceField   Reference field pointer
          *
          */
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
        void init(Coupling *coupling, CWP_Dof_location_t pointsCloudLocation, bool slave) override;

        /***********************************************************
         ***********************************************************
         **                                                       **
         **            Localization object functions              **
         **                                                       **
         ***********************************************************
         ***********************************************************/
        /**
           *
           * \brief Setting of the points cloud for localization.
           *
           * \param [out] id_dist   Localization object identifier.
           *
           */
        virtual void localization_points_cloud_setting();

        /**
          *
          * \brief Setting of the surface mesh and cloud points at
          * null for the localization object in a case of sending
          * code i.e. code which interpolate reference field.
          *
          * \param [out] id_dist   Localization object identifier.
          *
          */
        virtual void localization_null_setting_send();

        /**
          *
          * \brief Setting of the surface mesh and cloud points at
          * null for the localization object in a case of receving
          * code i.e. code which provides cloud points for interpolation.
          *
          * \param [out] id_dist   Localization object identifier.
          *
          */
        virtual void localization_null_setting_recv();

        /**
          *
          * \brief Setting of the surface mesh and cloud points at
          * null for the localization object in a case of receving
          * code i.e. code which provides cloud points for interpolation.
          *
          * \param [out] id_dist   Localization object identifier.
          *
          */
        virtual void localization_surface_setting();

        /**
          *
          * \brief Compute of localization of a points cloud on a surface
          *  mesh through the localization object.
          *
          * \param [int] id_dist   Localization object identifier.
          *
          */
        virtual void localization_compute();

        /**
          *
          * \brief Get localization results from localization object.
          *
          * \param [int] id_dist   Localization object identifier.
          *
          */
        virtual void localization_get();

        /**
          *
          * \brief Get localization results from localization object
          * from the coupled code in the case where the both codes are on
          * the same process.
          *
          * \param [int] id_dist   Localization object identifier.
          *
          */
        virtual void localization_get_cpl();

        virtual void localization_free();

        /***********************************************************
         ***********************************************************
         **                                                       **
         **   Process, partition, num triplet location from       **
         **           global numbering functions                  **
         **                                                       **
         ***********************************************************
         ***********************************************************/
        /**
           *
           * \brief Setting of requested global numbering for process, partition,
           *        num triplet location from global numbering object.
           *
           * \param [in] id_gnum_location  process, partition, num triplet location
           *              from global numbering identifier.
           *
           */
        void triplet_location_request();

        /**
           *
           * \brief Setting of researched global numbering for process, partition,
           *        num triplet location from global numbering object.
           *
           * \param [in] id_gnum_location  rocess, partition, num triplet location
           *              from global numbering identifier.
           *
           */
        void triplet_location_set();

        /**
          *
          * \brief Setting of researched global numbering for process, partition,
          *  num triplet location from global numbering object in a case of sending
          * code i.e. code which interpolates provided cloud points.
          *
          * \param [in] id_gnum_location  rocess, partition, num triplet location
          *              from global numbering identifier.
          *
          */
        void triplet_location_null_send();

        /**
          *
          * \brief Setting of researched global numbering for process, partition,
          *  num triplet location from global numbering object in a case of receving
          * code i.e. code which provides cloud points for interpolation.
          *
          * \param [in] id_gnum_location  rocess, partition, num triplet location
          *              from global numbering identifier.
          *
          */
        void triplet_location_null_recv();

        /**
          *
          * \brief Compute of process, partition, num triplet location
          *        from global numbering object.
          *
          * \param [in] id_gnum_location    Process, partition, num triplet location
          *                                 from global numbering identifier.
          *
          */
        void triplet_location_compute() const;

        /**
          *
          * \brief Get process, partition, num triplet location
          *        the case where the both codes are on
          *        the same process.
          *
          * \param [in] id_gnum_location     Process, partition, num triplet location
          *                                  from global numbering identifier.
          *
          */
        void triplet_location_get();

        /**
          *
          * \brief Get process, partition, num triplet location
          * from the coupled code in the case where the both codes are on
          * the same process.
          *
          * \param [in] id_gnum_location     Process, partition, num triplet location
          *                                  from global numbering identifier.
          */
        void triplet_location_get_cpl();

        /***********************************************************
         ***********************************************************
         **            Communication tree array functions         **
         **                                                       **
         ***********************************************************
         ***********************************************************/
        /**
          *
          * \brief Initialization of the communication tree array
          *        containing localization informations of the coupled
          *        mesh point cloud.
          *
          */
        void initialization_of_receving_communication_tree_array();

        /**
          *
          * \brief Filling of the communication tree array
          *        containing localization informations of the
          *        mesh point cloud.
          *
          */
        void filling_of_sending_communication_tree_array();

        // Triplet global numbering, MPI process, partition results
        int **_target_proc_part_num_idx{};    // Index array of triplet process, partition, numbering for each target
        int **_target_proc_part_num{};        // Array of triplet process, partition, numbering for each target

        // Paradigm structure identifier
        int _id_gnum_location{};              // Identifier for the global numbering to (process,partition,numbering) triplet object of paradigm
    }; //end SpatialInterpLocation

/**
 * \endcond
 */

}
#endif // SPATIALINTERPLOCATION_H_
