!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2019  ONERA
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

module pdm_overlay

  use pdm
  use iso_c_binding

  implicit none

  !!
  !! \enum PDM_ol_mesh_t
  !! \brief 3 Meshes
  !!
  !!

  integer(c_int), parameter :: PDM_OL_MESH_A=0  !< First mesh to overlay
  integer(c_int), parameter :: PDM_OL_MESH_B=1  !< Second mesh to overlay

  !!
  !! \enum PDM_ol_parameter_t
  !! \brief Parameters for ovelay meshes building
  !!
  !!

  integer(c_int), parameter ::   PDM_OL_CAR_LENGTH_TOL = 0 !< Absolute tolerance for caracteristic length
  integer(c_int), parameter :: PDM_OL_EXTENTS_TOL    = 1   !< Absolute tolerance for extents
  integer(c_int), parameter :: PDM_OL_SAME_PLANE_TOL = 2   !< Absolute tolerance for check if 2 surfaces are the same plane surface

  !!
  !! \enum PDM_ol_mv_t
  !! \brief Type of moving mesh
  !!

  integer(c_int), parameter :: PDM_OL_MV_TRANSFORMATION  = 0 !< Moving with combination of geometric transformations
  integer(c_int), parameter :: PDM_OL_MV_UNKNOWN         = 1 !< Unknown moving type

  !! No iso c binding interface : See PROCF in pdm_overlay.[ch] to find functions to call

  ! interface




  ! end interface

end module pdm_overlay
