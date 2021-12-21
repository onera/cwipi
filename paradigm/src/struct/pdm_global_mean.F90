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

module pdm_global_mean

  use pdm

  implicit none

  interface

    function pdm_global_mean_create (n_part, fComm) &
                                       result (ptrC) &
         bind (c, name = 'PDM_global_mean_create_cf')

      use iso_c_binding

      implicit none

      integer(c_int), value :: n_part
      integer(c_int), value :: fComm

      type (c_ptr) :: ptrC

    end function pdm_global_mean_create


    subroutine pdm_global_mean_set (ptrC, i_part, n_point, numabs) &
         bind (c, name = 'PDM_global_mean_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: i_part
      integer(c_int), value :: n_point
      type (c_ptr), value :: numabs

      type (c_ptr), value :: ptrC

    end subroutine pdm_global_mean_set

    subroutine pdm_global_mean_field_set (ptrC, i_part, stride, local_field, local_weight, global_mean_field_ptr) &
         bind (c, name = 'PDM_global_mean_field_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: i_part
      integer(c_int), value :: stride
      type (c_ptr), value :: local_field
      type (c_ptr), value :: local_weight
      type (c_ptr), value :: global_mean_field_ptr

      type (c_ptr), value :: ptrC

    end subroutine pdm_global_mean_field_set

    subroutine pdm_global_mean_field_compute (ptrC) &
         bind (c, name = 'PDM_global_mean_field_compute')

      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptrC

    end subroutine pdm_global_mean_field_compute

    subroutine pdm_global_mean_free (ptrC) &
         bind (c, name = 'PDM_global_mean_free')

      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptrC

    end subroutine PDM_global_mean_free

  end interface

end module pdm_global_mean
