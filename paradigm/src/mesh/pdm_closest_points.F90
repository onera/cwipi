!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2020  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

#include "pdm_configf.h"

module pdm_closest_points

  use pdm
  implicit none

  interface

    !>
    !!
    !! \brief Create a structure to look for the closest points of a point cloud
    !! (target cloud) in an other point cloud (source cloud)
    !!
    !! \param [in]   comm           MPI communicator
    !! \param [in]   n_closest      Number of closest source points to find for each
    !!                              target point
    !!
    !! \return     Identifier
    !!
    !!

    function PDM_closest_points_create (fcomm,     &
                                        n_closest, &
                                        owner)     &
                                        result(cls) &
      bind (c, name = 'PDM_closest_points_create_cf')

      use iso_c_binding

      implicit none

      integer(c_int), value :: fComm
      integer(c_int), value :: n_closest
      integer(c_int), value :: owner

      type(c_ptr)        :: cls

    end function PDM_closest_points_create

    !>
    !!
    !! \brief Set the number of partitions of a point cloud
    !!
    !! \param [in]   id                Identifier
    !! \param [in]   n_part_cloud_src  Number of partitions of the source cloud
    !! \param [in]   n_part_cloud_tgt  Number of partitions od the target cloud
    !!
    !!

    subroutine PDM_closest_points_n_part_cloud_set (cls, &
                                                    n_part_cloud_src, &
                                                    n_part_cloud_tgt) &
      bind (c, name = 'PDM_closest_points_n_part_cloud_set')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

      integer(c_int), value :: n_part_cloud_src
      integer(c_int), value :: n_part_cloud_tgt


    end subroutine PDM_closest_points_n_part_cloud_set

    !>
    !!
    !! \brief Set the target point cloud
    !!
    !! \param [in]   id              Identifier
    !! \param [in]   i_part          Index of partition
    !! \param [in]   n_points        Number of points
    !! \param [in]   coords          Point coordinates
    !! \param [in]   gnum            Point global number
    !!
    !!

    subroutine PDM_closest_points_tgt_cloud_set (cls, &
                                                 i_part, &
                                                 n_points, &
                                                 coords, &
                                                 gnum) &
      bind (c, name = 'PDM_closest_points_tgt_cloud_set')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

      integer(c_int), value :: i_part
      integer(c_int), value :: n_points

      type(c_ptr), value    :: coords
      type(c_ptr), value    :: gnum


    end subroutine PDM_closest_points_tgt_cloud_set

    !>
    !!
    !! \brief Set the source point cloud
    !!
    !! \param [in]   id              Identifier
    !! \param [in]   i_part          Index of partition
    !! \param [in]   n_points        Number of points
    !! \param [in]   coords          Point coordinates
    !! \param [in]   gnum            Point global number
    !!
    !!

    subroutine PDM_closest_points_src_cloud_set (cls, &
                                                 i_part, &
                                                 n_points, &
                                                 coords, &
                                                 gnum) &
      bind (c, name = 'PDM_closest_points_src_cloud_set')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

      integer(c_int), value :: i_part
      integer(c_int), value :: n_points

      type(c_ptr), value    :: coords
      type(c_ptr), value    :: gnum


    end subroutine PDM_closest_points_src_cloud_set

    !>
    !!
    !! \brief Look for closest points
    !!
    !! \param [in]   id  Identifier
    !!
    !!

    subroutine PDM_closest_points_compute (cls) &
      bind (c, name = 'PDM_closest_points_compute')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

    end subroutine PDM_closest_points_compute

    !>
    !!
    !! \brief Get mesh distance
    !!
    !! \param [in]   id                    Identifier
    !! \param [in]   i_part_tgt            Index of partition of the cloud
    !! \param [out]  closest_src_g_num     Global number of the closest element (size = n_closest * n_tgt_points)
    !! \param [out]  closest_src_distance  Distance (size = n_closest * n_tgt_points)
    !!
    !!

    subroutine PDM_closest_points_get (cls, &
                                       i_part_tgt, &
                                       closest_src_gnum, &
                                       closest_src_distance) &
      bind (c, name = 'PDM_closest_points_get')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

      integer(c_int), value :: i_part_tgt

      type(c_ptr)      :: closest_src_gnum
      type(c_ptr)      :: closest_src_distance


    end subroutine PDM_closest_points_get

    !>
    !!
    !! \brief Free a distance mesh structure
    !!
    !! \param [in]  id       Identifier
    !! \param [in]  partial  if partial is equal to 0, all data are removed.
    !!                       Otherwise, results are kept.
    !!
    !!

    subroutine PDM_closest_points_free (cls) &
     bind (c, name = 'PDM_closest_points_free')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

    end subroutine PDM_closest_points_free

    !>
    !!
    !! \brief  Dump elapsed an CPU time
    !!
    !! \param [in]  id       Identifier
    !!
    !!

    subroutine PDM_closest_points_dump_times (cls) &
      bind (c, name = 'PDM_closest_points_dump_times')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

    end subroutine PDM_closest_points_dump_times

  end interface

end module pdm_closest_points
