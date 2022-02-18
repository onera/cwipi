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

module pdm_block_to_part

  use pdm

  implicit none

  interface


    !>
    !!
    !! \brief Create a block to partitions redistribution
    !!
    !! \param [in]   blockDistribIdx Block distribution (size : \ref size of \ref comm + 1)
    !! \param [in]   gnum_elt        Element global number (size : \ref n_part)
    !! \param [in]   n_elt           Local number of elements (size : \ref n_part)
    !! \param [in]   n_part          Number of partition
    !! \param [in]   comm            MPI communicator
    !!
    !! \return   Initialized \ref PDM_block_to_part instance
    !!
    !!

    function PDM_block_to_part_create (blockDistribIdx, &
                                       gnum_elt, &
                                       n_elt, &
                                       n_part, &
                                       fcomm) &
                                       result (btp) &
      bind (c, name = 'PDM_block_to_part_create_cf')

      use iso_c_binding

      implicit none

      type (c_ptr), value :: blockDistribIdx
      type (c_ptr), value :: gnum_elt
      type (c_ptr), value :: n_elt
      integer (c_int), value :: n_part
      integer (c_int), value :: fcomm

      type (c_ptr) :: btp

      !> A l'appel :
      !>   integer (kind = pdm_g_num_s), pointer  :: n_elt(n_rank + 1) Index de distribution des blocs
      !>   type(c_ptr), pointer     :: gnum_elt(n_part)
      !>   allocate(gnum_elt(n_part) ) Tableau de d'adresses C contenant l'adresse de chaque tableau de numerotation
      !>                               globale de chaque partition
      !>   do i=1,n_part
      !>     gnum_elt(i)=c_loc(gnum_elt_fortran_i_part)  (Attention les tableaux gnum_elt_fortran_i_part contiennent
      !<                                                  des entiers de taille pdm_g_num_s)
      !>   enddo
      !>   integer (c_int), pointer  :: n_elt(n_part)
      !>

    end function PDM_block_to_part_create

    !>
    !!
    !! \brief Exchange data from blocks to partitions (part_stride and part_data are allocated by user)
    !!
    !! Prefer PDM_block_to_part_exch2 which allocates the results (part_stride and part_data)
    !!
    !! \param [in]   btp          Block to part structure
    !! \param [in]   s_data       Data size
    !! \param [in]   t_stride     Stride type
    !! \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
    !!                            Constant stride for \ref PDM_STRIDE_VAR
    !! \param [in]   block_data   Block data
    !! \param [out]  part_stride  Partition stride
    !! \param [out]  part_data    Partition data
    !!
    !!

    subroutine PDM_block_to_part_exch (btp, &
                                       s_data, &
                                       t_stride, &
                                       block_stride, &
                                       block_data, &
                                       part_stride, &
                                       part_data) &
      bind (c, name = 'PDM_block_to_part_exch')

      use iso_c_binding

      implicit none
      type (c_ptr), value       :: btp
      integer (c_size_t), value :: s_data
      integer (c_int), value    :: t_stride !PDM_STRIDE_CST or PDM_STRIDE_VAR
      type (c_ptr), value       :: block_stride ! Array of stride if variable otherwize
      type (c_ptr), value       :: block_data
      type (c_ptr), value       :: part_stride
      type (c_ptr), value       :: part_data

      !> A l'appel : part_stride n'est pas utilise si PDM_STRIDE_CST. A allouer comme gnum_elt dans PDM_block_to_part_create
      !>  (tableau d'adresses C)
      !> part_data a allouer comme gnum_elt dans PDM_block_to_part_create
      !>  (tableau d'adresses C)

    end subroutine PDM_block_to_part_exch


    !>
    !!
    !! \brief Initialize an exchange
    !! (part_stride and part_data are allocated in function)
    !!
    !! \param [in]   btp          Block to part structure
    !! \param [in]   s_data       Data size
    !! \param [in]   t_stride     Stride type
    !! \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
    !!                            Constant stride for \ref PDM_STRIDE_VAR
    !! \param [in]   block_data   Block data
    !! \param [out]  part_stride  Partition stride
    !! \param [out]  part_data    Partition data
    !!
    !!

    subroutine PDM_block_to_part_exch2 (btp, &
                                        s_data, &
                                        t_stride, &
                                        block_stride, &
                                        block_data, &
                                        part_stride, &
                                        part_data) &
      bind (c, name = 'PDM_block_to_part_exch2')

      use iso_c_binding

      implicit none
      type (c_ptr), value       :: btp
      integer (c_size_t), value :: s_data
      integer (c_int), value    :: t_stride !PDM_STRIDE_CST or PDM_STRIDE_VAR
      type (c_ptr), value       :: block_stride ! Array of stride if variable otherwize
      type (c_ptr), value       :: block_data
      type (c_ptr)              :: part_stride
      type (c_ptr)              :: part_data

    end subroutine PDM_block_to_part_exch2


    !>
    !!
    !! \brief Free a block to part structure
    !!
    !! \param [inout] btp  Block to part structure
    !!
    !! \return       NULL
    !!

    function PDM_block_to_part_free (btp)&
                                    result(ptr_null)&
      bind (c, name = 'PDM_block_to_part_free')
      use iso_c_binding

      implicit none

      type (c_ptr), value :: btp
      type (c_ptr) :: ptr_null

    end function PDM_block_to_part_free

    !>
    !!
    !! \brief Return index in the block for a gnum
    !!
    !! \param [in] ptb         Part to block structure
    !! \param [in] gNum        Global number
    !!
    !! \return  Index
    !!

    function PDM_block_to_part_gnum_idx_get (btp, &
                                             gnum)&
                                             result(index)&
      bind (c, name = 'PDM_block_to_part_gnum_idx_get')
      use iso_c_binding
      use pdm

      implicit none

      type (c_ptr), value :: btp
#ifdef PDM_LONG_G_NUM
      integer (c_long),  value :: gnum
#else
      integer (c_int),  value :: gnum
#endif
      integer (c_int)     :: index

    end function PDM_block_to_part_gnum_idx_get

  end interface

end module pdm_block_to_part
