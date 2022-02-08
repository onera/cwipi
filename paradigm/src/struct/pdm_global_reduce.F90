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

module pdm_global_reduce

  use pdm

  implicit none

  interface pdm_global_reduce_create ; module procedure &
  pdm_global_reduce_create_
  end interface

  private :: pdm_global_reduce_create_

  interface

  !>
  !!
  !! \brief Create a structure that computes a global reduction
  !!
  !! \param [in]   n_part       Number of local partitions
  !! \param [in]   fcomm        PDM_MPI communicator
  !!
  !! \return     Pointer to \ref PDM_global_reduce object
  !!

  function PDM_global_reduce_create_cf (n_part, &
                                        comm)   &
  result (gre)                                  &
  bind (c, name = 'PDM_global_reduce_create')

  use iso_c_binding
  implicit none

  integer(c_int), value :: n_part
  integer(c_int), value :: comm

  type (c_ptr)          :: gre

end function PDM_global_reduce_create_cf


!>
!!
!! \brief Free a global point reduce structure
!!
!! \param [in]   gre          Pointer to \ref PDM_global_reduce object
!!
!!

subroutine PDM_global_reduce_free (gre) &
  bind (c, name = 'PDM_global_reduce_free')

  use iso_c_binding
  implicit none

  type (c_ptr), value   :: gre

end subroutine pdm_global_reduce_free


!>
!!
!! \brief Set absolute number
!!
!! \param [in]   gre           Pointer to \ref PDM_global_reduce object
!! \param [in]   i_part        Current partition
!! \param [in]   n_pts         Number of points in the partition
!! \param [in]   pts_ln_to_gn  Global ids of points in the partition
!!
!!

subroutine PDM_global_reduce_g_num_set ( &
  gre,                                   &
  i_part,                                &
  n_pts,                                 &
  pts_ln_to_gn)                          &
bind (c, name = 'PDM_global_reduce_g_num_set')

use iso_c_binding
implicit none

type (c_ptr),   value :: gre
integer(c_int), value :: i_part
integer(c_int), value :: n_pts
type (c_ptr),   value :: pts_ln_to_gn

end subroutine PDM_global_reduce_g_num_set


!>
!!
!! \brief Set reduction operation
!!
!! \param [in]   gre          Pointer to \ref PDM_global_reduce object
!! \param [in]   operation    Type of reduction operation
!!

subroutine PDM_global_reduce_operation_set ( &
  gre,                                       &
  operation)                                 &
bind (c, name = 'PDM_global_reduce_operation_set')

use iso_c_binding
implicit none

type (c_ptr),   value :: gre
integer(c_int), value :: operation

end subroutine PDM_global_reduce_operation_set


!>
!!
!! \brief Set local field
!!
!! \param [in]   gre                       Pointer to \ref PDM_global_reduce object
!! \param [in]   i_part                    Current partition
!! \param [in]   stride                    Stride of the field
!! \param [in]   local_field               Local value of field
!!                                         (can be NULL for any other reduction operation)
!! \param [in]   global_reduced_field_ptr  Pointer where global reduced field
!!                                         will be stored after computing
!!

subroutine PDM_global_reduce_field_set ( &
  gre,                                   &
  i_part,                                &
  stride,                                &
  local_field,                           &
  global_reduced_field_ptr)              &
bind (c, name = 'PDM_global_reduce_field_set')

use iso_c_binding
implicit none

type (c_ptr),   value :: gre
integer(c_int), value :: i_part
integer(c_int), value :: stride
type (c_ptr),   value :: local_field
type (c_ptr),   value :: global_reduced_field_ptr

end subroutine PDM_global_reduce_field_set


!>
!!
!! \brief Compute the global reduced field
!!
!! \param [in]   gre     Pointer to \ref PDM_global_reduce object
!!
!!

subroutine PDM_global_reduce_field_compute (gre) &
  bind (c, name = 'PDM_global_reduce_field_compute')

  use iso_c_binding
  implicit none

  type (c_ptr),   value :: gre

end subroutine PDM_global_reduce_field_compute

end interface


contains

  !>
  !!
  !! \brief Create a structure that computes a global reduction
  !!
  !! \param [out]  gre      Pointer to \ref PDM_global_reduce object
  !! \param [in]   n_part   Number of local partitions
  !! \param [in]   f_comm   PDM_MPI communicator
  !!

  subroutine PDM_global_reduce_create_ (gre,    &
                                        n_part, &
                                        f_comm)

  use iso_c_binding
  implicit none

  integer       :: n_part
  integer       :: f_comm
  type (c_ptr)  :: gre

  integer(c_int) :: c_n_part
  integer(c_int) :: c_comm

  c_comm   = PDM_MPI_Comm_f2c(f_comm)
  c_n_part = n_part

  gre = pdm_global_reduce_create_cf(c_n_part, &
                                    c_comm)

end subroutine PDM_global_reduce_create_

end module pdm_global_reduce
