!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2022  ONERA
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

module pdm_part_to_part

  use pdm
  use pdm_pointer_array

  implicit none

  integer, parameter :: PDM_PART_TO_PART_DATA_DEF_ORDER_PART1           = 0
  integer, parameter :: PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2  = 1
  integer, parameter :: PDM_PART_TO_PART_DATA_DEF_ORDER_PART2           = 2
  integer, parameter :: PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM = 3


interface

!>
!!
!! \brief Free a part to part structure
!!
!! \param [inout] ptp  Part to part structure
!!
!! \return       NULL
!!

subroutine PDM_part_to_part_free (ptp) &
bind (c, name='PDM_part_to_part_free')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptp

end subroutine PDM_part_to_part_free


!>
!!
!! \brief Get number of partitions
!!
!! \param [in]  ptp       Pointer to \ref PDM_part_to_part_t object
!! \param [out] n_part1   Number of partitions on side 1
!! \param [out] n_part2   Number of partitions on side 2
!!
!!

subroutine PDM_part_to_part_n_part_get (ptp,     &
                                        n_part1, &
                                        n_part2) &
bind (c, name='PDM_part_to_part_n_part_get')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptp
  integer(c_int)     :: n_part1
  integer(c_int)     :: n_part2

end subroutine PDM_part_to_part_n_part_get


!>
!!
!! \brief Wait a partial asynchronus exchange
!!
!! \param [in]  ptp      Part to part structure
!! \param [in]  request  Request
!!
!!

subroutine PDM_part_to_part_iexch_wait (ptp,     &
                                        request) &
bind (c, name='PDM_part_to_part_iexch_wait')
  use iso_c_binding
  implicit none

  type(c_ptr),    value :: ptp
  integer(c_int), value :: request

end subroutine PDM_part_to_part_iexch_wait


!>
!!
!! \brief Wait a partial asynchronus exchange
!!
!! \param [in]  ptp      Part to part structure
!! \param [in]  request  Request
!!
!!

subroutine PDM_part_to_part_reverse_iexch_wait (ptp,     &
                                                request) &
bind (c, name='PDM_part_to_part_reverse_iexch_wait')
  use iso_c_binding
  implicit none

  type(c_ptr),    value :: ptp
  integer(c_int), value :: request

end subroutine PDM_part_to_part_reverse_iexch_wait

end interface



contains


!>
!!
!! \brief Create a partitions to partitions redistribution
!!
!! \param [out]  ptp                Initialized \ref PDM_part_to_part instance
!! \param [in]   gnum_elt1          Element global number (size : \ref n_part1)
!! \param [in]   n_elt1             Local number of elements (size : \ref n_part1)
!! \param [in]   n_part1            Number of partition
!! \param [in]   gnum_elt2          Element global number (size : \ref n_part2)
!! \param [in]   n_elt2             Local number of elements (size : \ref n_part2)
!! \param [in]   n_part2            Number of partition
!! \param [in]   part1_to_part2_idx Index of data to send to gnum2 from gnum1
!!                                  (for each part size : \ref n_elt1+1)
!! \param [in]   part1_to_part2     Data to send to gnum2 from gnum1
!! \param [in]   comm               MPI communicator
!!

subroutine PDM_part_to_part_create (ptp,                &
                                    gnum_elt1,          &
                                    n_elt1,             &
                                    n_part1,            &
                                    gnum_elt2,          &
                                    n_elt2,             &
                                    n_part2,            &
                                    part1_to_part2_idx, &
                                    part1_to_part2,     &
                                    comm)
  use iso_c_binding
  implicit none

  type(c_ptr)                       :: ptp
  type(PDM_pointer_array_t), target :: gnum_elt1
  integer(pdm_l_num_s), pointer     :: n_elt1(:)
  integer, intent(in)               :: n_part1
  type(PDM_pointer_array_t), target :: gnum_elt2
  integer(pdm_l_num_s), pointer     :: n_elt2(:)
  integer, intent(in)               :: n_part2
  type(PDM_pointer_array_t), target :: part1_to_part2_idx
  type(PDM_pointer_array_t), target :: part1_to_part2
  integer, intent(in)               :: comm

  integer(c_int)                    :: c_comm

  interface
    function pdm_part_to_part_create_c (gnum_elt1,          &
                                        n_elt1,             &
                                        n_part1,            &
                                        gnum_elt2,          &
                                        n_elt2,             &
                                        n_part2,            &
                                        part1_to_part2_idx, &
                                        part1_to_part2,     &
                                        comm)               &
    result (ptp)                                            &
    bind(c, name='PDM_part_to_part_create')
      use iso_c_binding
      implicit none

      type(c_ptr)           :: ptp
      type(c_ptr),    value :: gnum_elt1
      type(c_ptr),    value :: n_elt1
      integer(c_int), value :: n_part1
      type(c_ptr),    value :: gnum_elt2
      type(c_ptr),    value :: n_elt2
      integer(c_int), value :: n_part2
      type(c_ptr),    value :: part1_to_part2_idx
      type(c_ptr),    value :: part1_to_part2
      integer(c_int), value :: comm

    end function pdm_part_to_part_create_c
  end interface

  c_comm = PDM_MPI_Comm_f2c(comm)

  ptp = pdm_part_to_part_create_c (c_loc(gnum_elt1%cptr),          &
                                   c_loc(n_elt1),                  &
                                   n_part1,                        &
                                   c_loc(gnum_elt2%cptr),          &
                                   c_loc(n_elt2),                  &
                                   n_part2,                        &
                                   c_loc(part1_to_part2_idx%cptr), &
                                   c_loc(part1_to_part2%cptr),     &
                                   c_comm)

end subroutine PDM_part_to_part_create


!>
!!
!! \brief Initialize a partial asynchronus exchange
!!
!! \param [in]   ptp              Part to part structure
!! \param [in]   k_comm           Kind of MPI communication
!! \param [in]   t_stride         Kind of stride
!! \param [in]   t_part1_data_def Kind of part1 data definition
!! \param [in]   cst_stride       Constant stride
!! \param [in]   part1_stride     Stride of partition 1 data
!! \param [in]   part1_data       Partition 1 data
!! \param [out]  part2_stride     Stride of partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
!! \param [out]  part2_data       Partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
!! \param [out]  request          Request
!!
!!

subroutine PDM_part_to_part_iexch (ptp,              &
                                   k_comm,           &
                                   t_stride,         &
                                   t_part1_data_def, &
                                   cst_stride,       &
                                   part1_stride,     &
                                   part1_data,       &
                                   part2_stride,     &
                                   part2_data,       &
                                   request)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptp
  integer, intent(in)               :: k_comm
  integer, intent(in)               :: t_stride
  integer, intent(in)               :: t_part1_data_def
  integer, intent(in)               :: cst_stride
  type(PDM_pointer_array_t), target :: part1_stride
  type(PDM_pointer_array_t), target :: part1_data
  type(PDM_pointer_array_t), target :: part2_stride
  type(PDM_pointer_array_t), target :: part2_data
  integer, intent(out)              :: request

  integer                           :: n_part1
  integer                           :: n_part2
  type(c_ptr)                       :: c_part2_stride = C_NULL_PTR
  type(c_ptr)                       :: c_part2_data   = C_NULL_PTR
  integer                           :: n_ref
  integer(pdm_l_num_s), pointer     :: ref(:)                 => null()
  integer(pdm_l_num_s), pointer     :: gnum1_come_from_idx(:) => null()
  integer(pdm_g_num_s), pointer     :: gnum1_come_from(:)     => null()
  integer(pdm_l_num_s), pointer     :: stride(:)              => null()
  integer                           :: s_part2_data
  integer                           :: i, j

  integer                           :: s_data

  interface
    subroutine PDM_part_to_part_iexch_c (ptp,              &
                                         k_comm,           &
                                         t_stride,         &
                                         t_part1_data_def, &
                                         cst_stride,       &
                                         s_data,           &
                                         part1_stride,     &
                                         part1_data,       &
                                         part2_stride,     &
                                         part2_data,       &
                                         request)          &
    bind(c, name='PDM_part_to_part_iexch')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: k_comm
      integer(c_int), value :: t_stride
      integer(c_int), value :: t_part1_data_def
      integer(c_int), value :: cst_stride
      integer(c_int), value :: s_data
      type(c_ptr),    value :: part1_stride
      type(c_ptr),    value :: part1_data
      type(c_ptr)           :: part2_stride
      type(c_ptr)           :: part2_data
      integer(c_int)        :: request

    end subroutine PDM_part_to_part_iexch_c
  end interface

  s_data = part1_data%s_data

  call PDM_part_to_part_iexch_c (ptp,                      &
                                 k_comm,                   &
                                 t_stride,                 &
                                 t_part1_data_def,         &
                                 cst_stride,               &
                                 s_data,                   &
                                 c_loc(part1_stride%cptr), &
                                 c_loc(part1_data%cptr),   &
                                 c_part2_stride,           &
                                 c_part2_data,             &
                                 request)

  call PDM_part_to_part_n_part_get (ptp,     &
                                    n_part1, &
                                    n_part2)

  call c_f_pointer(c_part2_stride,    &
                   part2_stride%cptr, &
                   [n_part2])

  call c_f_pointer(c_part2_data,    &
                   part2_data%cptr, &
                   [n_part2])

  allocate(part2_stride%length(n_part2))
  allocate(part2_data%length(n_part2))
  part2_stride%type = PDM_TYPE_INT
  part2_data%type = part1_data%type
  part2_data%s_data = part1_data%s_data

  ! Compute lengths
  do i = 1, n_part2
    call PDM_part_to_part_ref_lnum2_get (ptp,   &
                                         i-1,   &
                                         n_ref, &
                                         ref)

    call PDM_part_to_part_gnum1_come_from_get (ptp,                 &
                                               i-1,                 &
                                               gnum1_come_from_idx, &
                                               gnum1_come_from)

    if (t_stride .eq. PDM_STRIDE_VAR_INTERLACED) then

      part2_stride%length(i) = gnum1_come_from_idx(n_ref+1)

      call PDM_pointer_array_part_get (part2_stride, &
                                       i-1,          &
                                       stride)

      s_part2_data = 0
      do j = 1,gnum1_come_from_idx(n_ref+1)
        s_part2_data = s_part2_data + stride(j)
      end do

    else

      s_part2_data = cst_stride * gnum1_come_from_idx(n_ref+1)

    end if

    part2_data%length(i) = s_part2_data
  end do

end subroutine PDM_part_to_part_iexch


!>
!!
!! \brief Initialize a partial reverse asynchronus exchange
!! \param [in]   ptp              Part to part structure
!! \param [in]   k_comm           Kind of MPI communication
!! \param [in]   t_stride         Kind of stride
!! \param [in]   t_part2_data_def Kind of part2 data definition
!! \param [in]   cst_stride       Constant stride
!! \param [in]   s_data           Data size
!! \param [in]   part2_stride     Stride of partition 1 data (Accordding to t_part2_data_def)
!! \param [in]   part2_data       Partition 1 data (Accordding to t_part2_data_def)
!! \param [out]  part1_stride     Stride of partition 2 data (order given by part1_to_part2)
!! \param [out]  part1_data       Partition 2 data (order given by part1_to_part2)
!! \param [out]  request          Request
!!

subroutine PDM_part_to_part_reverse_iexch (ptp,              &
                                           k_comm,           &
                                           t_stride,         &
                                           t_part2_data_def, &
                                           cst_stride,       &
                                           s_data,           &
                                           part2_stride,     &
                                           part2_data,       &
                                           part1_stride,     &
                                           part1_data,       &
                                           request)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptp
  integer, intent(in)               :: k_comm
  integer, intent(in)               :: t_stride
  integer, intent(in)               :: t_part2_data_def
  integer, intent(in)               :: cst_stride
  integer, intent(in)               :: s_data
  type(PDM_pointer_array_t), target :: part2_stride
  type(PDM_pointer_array_t), target :: part2_data
  type(PDM_pointer_array_t), target :: part1_stride
  type(PDM_pointer_array_t), target :: part1_data
  integer, intent(out)              :: request

  integer                           :: n_part1
  integer                           :: n_part2
  type(c_ptr)                       :: c_part1_stride = C_NULL_PTR
  type(c_ptr)                       :: c_part1_data   = C_NULL_PTR
  integer                           :: n_elt1
  integer(pdm_l_num_s), pointer     :: part1_to_part2_idx(:) => null()
  integer(pdm_g_num_s), pointer     :: part1_to_part2(:)     => null()
  integer(pdm_l_num_s), pointer     :: stride(:)             => null()
  integer                           :: s_part1_data
  integer                           :: i, j

  interface
    subroutine PDM_part_to_part_reverse_iexch_c (ptp,              &
                                                 k_comm,           &
                                                 t_stride,         &
                                                 t_part2_data_def, &
                                                 cst_stride,       &
                                                 s_data,           &
                                                 part2_stride,     &
                                                 part2_data,       &
                                                 part1_stride,     &
                                                 part1_data,       &
                                                 request)          &
    bind (c, name='PDM_part_to_part_reverse_iexch')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: k_comm
      integer(c_int), value :: t_stride
      integer(c_int), value :: t_part2_data_def
      integer(c_int), value :: cst_stride
      integer(c_int), value :: s_data
      type(c_ptr),    value :: part2_stride
      type(c_ptr),    value :: part2_data
      type(c_ptr)           :: part1_stride
      type(c_ptr)           :: part1_data
      integer(c_int)        :: request

    end subroutine PDM_part_to_part_reverse_iexch_c
  end interface

  call PDM_part_to_part_reverse_iexch_c (ptp,                      &
                                         k_comm,                   &
                                         t_stride,                 &
                                         t_part2_data_def,         &
                                         cst_stride,               &
                                         s_data,                   &
                                         c_loc(part2_stride%cptr), &
                                         c_loc(part2_data%cptr),   &
                                         c_part1_stride,           &
                                         c_part1_data,             &
                                         request)

  call PDM_part_to_part_n_part_get (ptp,     &
                                    n_part1, &
                                    n_part2)

  call c_f_pointer(c_part1_stride,    &
                   part1_stride%cptr, &
                   [n_part1])

  call c_f_pointer(c_part1_data,    &
                   part1_data%cptr, &
                   [n_part1])

  allocate(part1_stride%length(n_part1))
  allocate(part1_data%length(n_part1))
  part1_stride%type = PDM_TYPE_INT
  part1_data%type = part2_data%type

  do i = 1, n_part1
    call PDM_part_to_part_part1_to_part2_get (ptp,                &
                                              i-1,                &
                                              n_elt1,             &
                                              part1_to_part2_idx, &
                                              part1_to_part2)

    if (t_stride .eq. PDM_STRIDE_VAR_INTERLACED) then

      part1_stride%length(i) = part1_to_part2_idx(n_elt1+1)

      call PDM_pointer_array_part_get (part1_stride, &
                                       i-1,          &
                                       stride)

      s_part1_data = 0
      do j = 1,part1_to_part2_idx(n_elt1+1)
        s_part1_data = s_part1_data + stride(j)
      end do

    else

      s_part1_data = cst_stride * part1_to_part2_idx(n_elt1+1)

    end if

    part1_data%length(i) = s_part1_data
  end do

end subroutine PDM_part_to_part_reverse_iexch


!>
!!
!! \brief Get referenced gnum2 elements
!!
!! \param [in]   ptp             Pointer to \ref PDM_part_to_part_t object
!! \param [in]   i_part          Id of partition
!! \param [out]  n_ref_lnum2     Number of referenced gnum2
!! \param [out]  ref_lnum2       Referenced gnum2
!!
!!

subroutine PDM_part_to_part_ref_lnum2_get (ptp,         &
                                           i_part,      &
                                           n_ref_lnum2, &
                                           ref_lnum2)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptp
  integer, intent(in)           :: i_part
  integer, intent(out)          :: n_ref_lnum2
  integer(pdm_l_num_s), pointer :: ref_lnum2(:)

  type(c_ptr)                   :: c_ref_lnum2 = C_NULL_PTR

  interface

    subroutine PDM_part_to_part_ref_lnum2_get_c (ptp,         &
                                                 i_part,      &
                                                 n_ref_lnum2, &
                                                 ref_lnum2)   &
    bind (c, name='PDM_part_to_part_ref_lnum2_single_part_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: i_part
      integer(c_int)        :: n_ref_lnum2
      type(c_ptr)           :: ref_lnum2

    end subroutine PDM_part_to_part_ref_lnum2_get_c

  end interface

  call PDM_part_to_part_ref_lnum2_get_c (ptp,         &
                                         i_part,      &
                                         n_ref_lnum2, &
                                         c_ref_lnum2)

  call c_f_pointer(c_ref_lnum2,   &
                   ref_lnum2,     &
                   [n_ref_lnum2])

end subroutine PDM_part_to_part_ref_lnum2_get


!>
!!
!! \brief Get unreferenced gnum2 elements
!!
!! \param [in]   ptp             Pointer to \ref PDM_part_to_part_t object
!! \param [in]   i_part          Id of partition
!! \param [out]  n_unref_lnum2   Number of unreferenced gnum2
!! \param [out]  unref_lnum2     Unreferenced gnum2
!!
!!

subroutine PDM_part_to_part_unref_lnum2_get (ptp,           &
                                             i_part,        &
                                             n_unref_lnum2, &
                                             unref_lnum2)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptp
  integer, intent(in)           :: i_part
  integer, intent(out)          :: n_unref_lnum2
  integer(pdm_l_num_s), pointer :: unref_lnum2(:)

  type(c_ptr)                   :: c_unref_lnum2 = C_NULL_PTR

  interface

    subroutine PDM_part_to_part_unref_lnum2_get_c (ptp,           &
                                                   i_part,        &
                                                   n_unref_lnum2, &
                                                   unref_lnum2)   &
    bind (c, name='PDM_part_to_part_unref_lnum2_single_part_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: i_part
      integer(c_int)        :: n_unref_lnum2
      type(c_ptr)           :: unref_lnum2

    end subroutine PDM_part_to_part_unref_lnum2_get_c

  end interface

  call PDM_part_to_part_unref_lnum2_get_c (ptp,           &
                                           i_part,        &
                                           n_unref_lnum2, &
                                           c_unref_lnum2)

  call c_f_pointer(c_unref_lnum2,   &
                   unref_lnum2,     &
                   [n_unref_lnum2])

end subroutine PDM_part_to_part_unref_lnum2_get


!>
!!
!! \brief Get gnum come from gnum1 for each referenced gnum2
!!
!! \param [in]   ptp                 Pointer to \ref PDM_part_to_part_t object
!! \param [in]   i_part              Id of partition
!! \param [out]  gnum1_come_from_idx Index for gnum1_come_from array
!! \param [out]  gnum1_come_from     Gnum come from gnum1 for each referenced gnum2
!!
!!

subroutine PDM_part_to_part_gnum1_come_from_get (ptp,                 &
                                                 i_part,              &
                                                 gnum1_come_from_idx, &
                                                 gnum1_come_from)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptp
  integer, intent(in)           :: i_part
  integer(pdm_l_num_s), pointer :: gnum1_come_from_idx(:)
  integer(pdm_g_num_s), pointer :: gnum1_come_from(:)

  type(c_ptr)                   :: c_gnum1_come_from_idx = C_NULL_PTR
  type(c_ptr)                   :: c_gnum1_come_from     = C_NULL_PTR
  integer                       :: n_ref
  integer(pdm_l_num_s), pointer :: ref(:) => null()

  interface
    subroutine PDM_part_to_part_gnum1_come_from_get_c (ptp,                 &
                                                       i_part,              &
                                                       gnum1_come_from_idx, &
                                                       gnum1_come_from)     &
    bind (c, name='PDM_part_to_part_gnum1_come_from_single_part_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: i_part
      type(c_ptr)           :: gnum1_come_from_idx
      type(c_ptr)           :: gnum1_come_from

    end subroutine PDM_part_to_part_gnum1_come_from_get_c
  end interface

  call PDM_part_to_part_ref_lnum2_get (ptp,    &
                                       i_part, &
                                       n_ref,  &
                                       ref)

  call PDM_part_to_part_gnum1_come_from_get_c (ptp,                   &
                                               i_part,                &
                                               c_gnum1_come_from_idx, &
                                               c_gnum1_come_from)

  call c_f_pointer(c_gnum1_come_from_idx, &
                   gnum1_come_from_idx,   &
                   [n_ref + 1])

  call c_f_pointer(c_gnum1_come_from, &
                   gnum1_come_from,   &
                   [gnum1_come_from_idx(n_ref+1)])

end subroutine PDM_part_to_part_gnum1_come_from_get


!>
!!
!! \brief
!!
!! \param [in]   ptp                 Pointer to \ref PDM_part_to_part_t object
!! \param [in]   i_part              Id of partition
!! \param [out]  part1_to_part2_idx  Index for part1_to_part2 array
!! \param [out]  part1_to_part2      Gnum2 for each element of part1
!!
!!

subroutine PDM_part_to_part_part1_to_part2_get (ptp,                &
                                                i_part,             &
                                                n_elt1,             &
                                                part1_to_part2_idx, &
                                                part1_to_part2)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptp
  integer, intent(in)           :: i_part
  integer, intent(out)          :: n_elt1
  integer(pdm_l_num_s), pointer :: part1_to_part2_idx(:)
  integer(pdm_g_num_s), pointer :: part1_to_part2(:)

  type(c_ptr)                   :: c_part1_to_part2_idx = C_NULL_PTR
  type(c_ptr)                   :: c_part1_to_part2     = C_NULL_PTR

interface
  subroutine PDM_part_to_part_part1_to_part2_get_c (ptp,                &
                                                    i_part,             &
                                                    n_elt1,             &
                                                    part1_to_part2_idx, &
                                                    part1_to_part2)     &
  bind (c, name='PDM_part_to_part_part1_to_part2_single_part_get')
    use iso_c_binding
    implicit none

    type(c_ptr),    value :: ptp
    integer(c_int), value :: i_part
    integer(c_int)        :: n_elt1
    type(c_ptr)           :: part1_to_part2_idx
    type(c_ptr)           :: part1_to_part2

  end subroutine PDM_part_to_part_part1_to_part2_get_c
end interface

  call PDM_part_to_part_part1_to_part2_get_c (ptp,                  &
                                              i_part,               &
                                              n_elt1,               &
                                              c_part1_to_part2_idx, &
                                              c_part1_to_part2)

  call c_f_pointer(c_part1_to_part2_idx, &
                   part1_to_part2_idx,   &
                   [n_elt1 + 1])

  call c_f_pointer(c_part1_to_part2, &
                   part1_to_part2,   &
                   [part1_to_part2_idx(n_elt1+1)])

end subroutine PDM_part_to_part_part1_to_part2_get


end module pdm_part_to_part
