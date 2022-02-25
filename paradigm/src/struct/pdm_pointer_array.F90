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

module pdm_pointer_array

  use pdm
  use pdm_fortran
  use iso_c_binding

  implicit none

  integer, parameter :: PDM_TYPE_INT    = 0
  integer, parameter :: PDM_TYPE_G_NUM  = 1
  integer, parameter :: PDM_TYPE_DOUBLE = 2




  type PDM_pointer_array_t

    integer                       :: type = -1
    type(c_ptr),          pointer :: cptr(:)   => null()
    integer(pdm_l_num_s), pointer :: length(:) => null()

  end type PDM_pointer_array_t




  interface PDM_pointer_array_part_set
    module procedure PDM_pointer_array_part_set_int
#ifdef PDM_LONG_G_NUM
    module procedure PDM_pointer_array_part_set_g_num
#endif
    module procedure PDM_pointer_array_part_set_double
  end interface

  interface PDM_pointer_array_part_get
    module procedure PDM_pointer_array_part_get_int
#ifdef PDM_LONG_G_NUM
    module procedure PDM_pointer_array_part_get_g_num
#endif
    module procedure PDM_pointer_array_part_get_double
  end interface


  contains


  !>
  !! \brief Initialize a \ref PDM_pointer_array_t object
  !!
  !! \param [out]  pa      \ref PDM_pointer_array_t object
  !! \param [in]   n_part  Number of partitions
  !! \param [in]   type    Data type of pointers
  !!

  subroutine PDM_pointer_array_create (pa,     &
                                       n_part, &
                                       type)
    use iso_c_binding
    implicit none

    type(PDM_pointer_array_t), target  :: pa
    integer, intent(in)                :: n_part
    integer, intent(in)                :: type

    integer                            :: i

    pa%type = type
    allocate(pa%cptr(n_part))
    allocate(pa%length(n_part))

    do i = 1, n_part
      pa%cptr(i)   = C_NULL_PTR
      pa%length(i) = 0
    end do

  end subroutine PDM_pointer_array_create



  !>
  !! \brief Free a \ref PDM_pointer_array_t object
  !!
  !! \param [in, out]  pa      \ref PDM_pointer_array_t object
  !!

  subroutine PDM_pointer_array_free (pa)
    use iso_c_binding
    implicit none

    type(PDM_pointer_array_t)  :: pa

    if (associated(pa%cptr)) then
      deallocate(pa%cptr)
    end if

    if (associated(pa%length)) then
      deallocate(pa%length)
    end if

  end subroutine PDM_pointer_array_free


  !>
  !! \brief Free a \ref PDM_pointer_array_t object
  !!
  !! \param [in, out]  pa      \ref PDM_pointer_array_t object
  !!

  subroutine PDM_pointer_array_free_from_c (pa)
    use iso_c_binding
    implicit none

    type(PDM_pointer_array_t)  :: pa

    integer                    :: i

    if (associated(pa%cptr)) then
      do i = 1, size(pa%cptr)
        call pdm_fortran_free_c(pa%cptr(i))
      end do
      call pdm_fortran_free_c(c_loc(pa%cptr))
    end if

    if (associated(pa%length)) then
      deallocate(pa%length)
    end if

  end subroutine PDM_pointer_array_free_from_c



  !>
  !! \brief Set a partition from a Fortran pointer
  !!
  !! \param [in]  pa         Array of \ref PDM_pointer_array_t
  !! \param [in]  i_part     Id of partition
  !! \param [in]  pointer_f  Pointer to an integer array
  !!

  subroutine PDM_pointer_array_part_set_int (pa,        &
                                             i_part,    &
                                             pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_pointer_array_t), target  :: pa
    integer, intent(in)                :: i_part
    integer(pdm_l_num_s),      pointer :: pointer_f(:)

    if (pa%type .ne. PDM_TYPE_INT) then
      print *, "PDM_pointer_array_part_set_int : wrong type"
      stop
    end if

    if (i_part .ge. size(pa%cptr)) then
      print *, "PDM_pointer_array_part_set_int : wrong i_part"
      stop
    end if

    pa%cptr(i_part+1)   = c_loc(pointer_f)
    pa%length(i_part+1) = size(pointer_f)

  end subroutine PDM_pointer_array_part_set_int


  !>
  !! \brief Set a partition from a Fortran pointer
  !!
  !! \param [in]  pa         Array of \ref PDM_pointer_array_t
  !! \param [in]  i_part     Id of partition
  !! \param [in]  pointer_f  Pointer to a g_num array
  !!

  subroutine PDM_pointer_array_part_set_g_num (pa,        &
                                               i_part,    &
                                               pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_pointer_array_t), target  :: pa
    integer, intent(in)                :: i_part
    integer(pdm_g_num_s),      pointer :: pointer_f(:)

    if (pa%type .ne. PDM_TYPE_G_NUM) then
      print *, "PDM_pointer_array_part_set_g_num : wrong type"
      stop
    end if

    if (i_part .ge. size(pa%cptr)) then
      print *, "PDM_pointer_array_part_set_g_num : wrong i_part"
      stop
    end if

    pa%cptr(i_part+1)   = c_loc(pointer_f)
    pa%length(i_part+1) = size(pointer_f)

  end subroutine PDM_pointer_array_part_set_g_num


  !>
  !! \brief Set a partition from a Fortran pointer
  !!
  !! \param [in]  pa         Array of \ref PDM_pointer_array_t
  !! \param [in]  i_part     Id of partition
  !! \param [in]  pointer_f  Pointer to a double array
  !!

  subroutine PDM_pointer_array_part_set_double (pa,        &
                                                i_part,    &
                                                pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_pointer_array_t), target  :: pa
    integer, intent(in)                :: i_part
    double precision,          pointer :: pointer_f(:)

    if (pa%type .ne. PDM_TYPE_DOUBLE) then
      print *, "PDM_pointer_array_part_set_double : wrong type"
      stop
    end if

    if (i_part .ge. size(pa%cptr)) then
      print *, "PDM_pointer_array_part_set_double : wrong i_part"
      stop
    end if

    pa%cptr(i_part+1)   = c_loc(pointer_f)
    pa%length(i_part+1) = size(pointer_f)

  end subroutine PDM_pointer_array_part_set_double


  !>
  !! \brief Set a partition from a C pointer
  !!
  !! \param [in]  pa         Array of \ref PDM_pointer_array_t
  !! \param [in]  i_part     Id of partition
  !! \param [in]  pointer_c  C pointer
  !!

  subroutine PDM_pointer_array_part_set_from_cptr (pa,        &
                                                   i_part,    &
                                                   pointer_c, &
                                                   length)
    use iso_c_binding
    implicit none

    type(PDM_pointer_array_t), target  :: pa
    integer, intent(in)                :: i_part
    type(c_ptr), value                 :: pointer_c
    integer, intent(in)                :: length

    if (i_part .ge. size(pa%cptr)) then
      print *, "PDM_pointer_array_part_set_int : wrong i_part"
      stop
    end if

    pa%cptr(i_part+1)   = pointer_c
    pa%length(i_part+1) = length

  end subroutine PDM_pointer_array_part_set_from_cptr


  !>
  !! \brief Get a partition
  !!
  !! Maps a Fortran pointer onto a C pointer
  !!
  !! \param [in]       pa         Array of \ref PDM_pointer_array_t
  !! \param [in]       i_part     Id of partition
  !! \param [in, out]  pointer_f  Pointer to an integer array
  !!

  subroutine PDM_pointer_array_part_get_int (pa,        &
                                             i_part,    &
                                             pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_pointer_array_t), target  :: pa
    integer, intent(in)                :: i_part
    integer(pdm_l_num_s),      pointer :: pointer_f(:)

    if (pa%type .ne. PDM_TYPE_INT) then
      print *, "PDM_pointer_array_part_set_int : wrong type"
      stop
    end if

    if (i_part .ge. size(pa%cptr)) then
      print *, "PDM_pointer_array_part_set_int : wrong i_part"
      stop
    end if


    call c_f_pointer(pa%cptr(i_part+1),     &
                     pointer_f,             &
                     [pa%length(i_part+1)])

  end subroutine PDM_pointer_array_part_get_int


  !>
  !! \brief Get a partition
  !!
  !! Maps a Fortran pointer onto a C pointer
  !!
  !! \param [in]       pa         Array of \ref PDM_pointer_array_t
  !! \param [in]       i_part     Id of partition
  !! \param [in, out]  pointer_f  Pointer to a g_num array
  !!

#ifdef PDM_LONG_G_NUM
  subroutine PDM_pointer_array_part_get_g_num (pa,        &
                                               i_part,    &
                                               pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_pointer_array_t), target  :: pa
    integer, intent(in)                :: i_part
    integer(pdm_g_num_s),      pointer :: pointer_f(:)

    if (pa%type .ne. PDM_TYPE_G_NUM) then
      print *, "PDM_pointer_array_part_set_int : wrong type"
      stop
    end if

    if (i_part .ge. size(pa%cptr)) then
      print *, "PDM_pointer_array_part_set_int : wrong i_part"
      stop
    end if


    call c_f_pointer(pa%cptr(i_part+1),     &
                     pointer_f,             &
                     [pa%length(i_part+1)])

  end subroutine PDM_pointer_array_part_get_g_num
#endif


  !>
  !! \brief Get a partition
  !!
  !! Maps a Fortran pointer onto a C pointer
  !!
  !! \param [in]       pa         Array of \ref PDM_pointer_array_t
  !! \param [in]       i_part     Id of partition
  !! \param [in, out]  pointer_f  Pointer to a double array
  !!

  subroutine PDM_pointer_array_part_get_double (pa,        &
                                                i_part,    &
                                                pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_pointer_array_t), target  :: pa
    integer, intent(in)                :: i_part
    double precision,          pointer :: pointer_f(:)

    if (pa%type .ne. PDM_TYPE_double) then
      print *, "PDM_pointer_array_part_set_double : wrong type"
      stop
    end if

    if (i_part .ge. size(pa%cptr)) then
      print *, "PDM_pointer_array_part_set_double : wrong i_part"
      stop
    end if


    call c_f_pointer(pa%cptr(i_part+1),     &
                     pointer_f,             &
                     [pa%length(i_part+1)])

  end subroutine PDM_pointer_array_part_get_double





end module pdm_pointer_array
