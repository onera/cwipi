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

program testf

  use pdm
  use mpi
  use pdm_block_to_part
  use iso_c_binding
  use pdm_fortran

  implicit none

  integer :: code
  integer :: i_rank
  integer :: n_rank

  integer (kind = pdm_g_num_s), pointer :: block_distrib_index(:)
  type(c_ptr)                           :: cptr_block_distrib_index
  integer (kind = pdm_g_num_s), pointer :: gnum_elt(:)
  type(c_ptr), pointer                  :: cptr_gnum_elt(:)
  type(c_ptr)                           :: cptr_cptr_gnum_elt

  integer(c_int), parameter :: fComm = MPI_COMM_WORLD

  type(c_ptr)              :: btp

  integer(c_int), pointer   :: n_elt(:)
  type(c_ptr)              :: cptr_n_elt

  integer(c_int), parameter :: n_part = 1

  integer(c_int), parameter :: n_elt_case = 5

  integer(c_int), parameter :: t_stride = PDM_STRIDE_VAR

  integer(c_int), pointer :: block_stride(:)
  type(c_ptr)             :: cptr_block_stride !  c_ptr of C array containing block stride
                                            ! (NULL for this case)
  double precision, pointer :: block_data(:)
  type(c_ptr)               :: cptr_block_data  ! c_ptr of C array containing block data

  integer(c_int), pointer :: part_stride(:) ! Fortran array containing data for part 0
                                            ! (Unused for this case with a constant stride)
  type(c_ptr), pointer :: cptr_part_stride(:) ! Fortran array of c_ptr containing c_loc (part_stride)
                                              ! for each partition (for this case n_part = 1)
                                              ! (Unused for this case with a constant stride)
  type(c_ptr)          :: cptr_cptr_part_stride ! c_loc about cptr_part_stride

  integer(c_size_t), parameter :: s_data = 8 ! size of type of the exchanged data
                                             ! (for this case exchange data are 'real*8')

  double precision, pointer :: part_data(:) ! Fortran array containing data for part 0
  type(c_ptr), pointer  :: cptr_part_data(:) ! Fortran array of c_ptr containing c_loc (part_data)
                                            ! for each partition (for this case n_part = 1)
  type(c_ptr)          :: cptr_cptr_part_data ! c_loc about cptr_part_data

  !
  ! Init
  !

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, i_rank, code)
  call mpi_comm_size(mpi_comm_world, n_rank, code)

  if (n_rank .ne. 1) then
    print *,'Error : 2 MPI processes are mandatory'
    call mpi_finalize(code)
    stop
  end if

  allocate(block_distrib_index(2))
  block_distrib_index(1) = 0
  block_distrib_index(2) = 5

  cptr_block_distrib_index = c_loc(block_distrib_index)

  allocate(cptr_gnum_elt(n_part))
  allocate(gnum_elt(n_elt_case))
  gnum_elt(1) = 3
  gnum_elt(2) = 5
  gnum_elt(3) = 2
  gnum_elt(4) = 4
  gnum_elt(5) = 1

  cptr_gnum_elt(1) = c_loc(gnum_elt)
  cptr_cptr_gnum_elt = c_loc(cptr_gnum_elt)

  allocate(n_elt(n_part))
  n_elt(1) = n_elt_case
  cptr_n_elt = c_loc(n_elt)

  btp = PDM_block_to_part_create (cptr_block_distrib_index, &
                                  cptr_cptr_gnum_elt, &
                                  cptr_n_elt, &
                                  n_part, &
                                  fcomm)

  allocate(block_stride(n_elt_case))
  block_stride(1) = 1
  block_stride(2) = 1
  block_stride(3) = 1
  block_stride(4) = 1
  block_stride(5) = 1

  cptr_block_stride = c_loc(block_stride)

  allocate(block_data(n_elt_case))

  block_data(1) = 10.d0
  block_data(2) = 20.d0
  block_data(3) = 30.d0
  block_data(4) = 40.d0
  block_data(5) = 50.d0

  cptr_block_data = c_loc(block_data)

  call PDM_block_to_part_exch2 (btp, &
                                s_data, &
                                t_stride, &
                                cptr_block_stride, &
                                cptr_block_data, &
                                cptr_cptr_part_stride, &
                                cptr_cptr_part_data)

  call c_f_pointer(cptr_cptr_part_stride, cptr_part_stride, [n_part])
  call c_f_pointer(cptr_cptr_part_data, cptr_part_data, [n_part])

  call c_f_pointer(cptr_part_stride(1), part_stride, [n_elt_case])
  call c_f_pointer(cptr_part_data(1), part_data, [n_elt_case])

  btp = PDM_block_to_part_free (btp)

  print*, "part_stride", part_stride(1), part_stride(2), part_stride(3), part_stride(4), part_stride(5)
  print*, "part_data", part_data(1), part_data(2), part_data(3), part_data(4), part_data(5)

  deallocate(n_elt)
  deallocate(gnum_elt)
  deallocate(cptr_gnum_elt)

  deallocate(block_data)
  deallocate(block_stride)
  deallocate(block_distrib_index)

  call pdm_fortran_free_c (cptr_part_data(1))
  call pdm_fortran_free_c (cptr_cptr_part_data)

  call pdm_fortran_free_c (cptr_part_stride(1))
  call pdm_fortran_free_c (cptr_cptr_part_stride)

  call mpi_finalize(code)

end program testf
