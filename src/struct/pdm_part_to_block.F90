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

module pdm_part_to_block

  use pdm

  implicit none

  integer, parameter :: PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC          = 0  ! Distribute block on all processes
  integer, parameter :: PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE = 1  ! Distribute block on one processe pere node
  integer, parameter :: PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE      = 2  ! Distribute block on part of nodes

  integer, parameter :: PDM_PART_TO_BLOCK_POST_NOTHING              = 0    ! No post processing
  integer, parameter :: PDM_PART_TO_BLOCK_POST_CLEANUP              = 1  ! Cleanup multi-elements
  integer, parameter :: PDM_PART_TO_BLOCK_POST_MERGE                = 2  ! Merge multi-elements


  interface

    !> \brief Create a partitioning to block redistribution
    !!
    !! \param [in]   t_distrib       Distribution type
    !! \param [in]   t_post          Post processing type
    !! \param [in]   partActiveNode  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
    !! \param [in]   gnum_elt        Element global number
    !! \param [in]   weight          Weight of elements (or NULL)
    !! \param [in]   n_elt           Local number of elements
    !! \param [in]   n_part          Number of partition
    !! \param [in]   comm            MPI communicator
    !!
    !! \return   Initialized cs_part_to_block

    function PDM_part_to_block_create (t_distrib, &
                                       t_post, &
                                       part_active_node,&
                                       gnum_elt,&
                                       weight,&
                                       n_elt,&
                                       n_part,&
                                       fcomm) &
                                       result(ptrC) &
      bind (c, name = 'PDM_part_to_block_create_cf')

      use iso_c_binding

      implicit none

      integer(c_int), value :: t_distrib ! PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                         ! PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE,
                                         ! PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE
      integer(c_int), value :: t_post ! PDM_PART_TO_BLOCK_POST_NOTHING
                                      ! PDM_PART_TO_BLOCK_POST_CLEANUP
                                      ! PDM_PART_TO_BLOCK_POST_MERGE

      real (c_double), value :: part_active_node

      type (c_ptr), value :: gnum_elt
      type (c_ptr), value :: weight
      type (c_ptr), value :: n_elt
      integer(c_int), value :: n_part
      integer(c_int), value :: fComm
      type (c_ptr) :: ptrC

      !> La routine qui appelle cette interface doit disposer des modules
      !> use iso_c_binding
      !> use pdm_part_to_block

      !> La routine qui appelle cette interface doit de la variable
      !> type(c_ptr)              :: PDM_part_to_block_t

      !> Les petits trucs:
      !> type (c_ptr), value :: gnum_elt
      !>   passer tableau local2Global(:)-> gnum_elt
      !>   integer                  :: iCell
      !>   integer    , pointer     :: local2Global(:)
      !>   type(c_ptr), pointer     :: gnum_elt(:,:)=>null()
      !>
      !>   allocate( gnum_elt(1:size(local2Global),1:1) )
      !>   do iCell=1,size(local2Global)
      !>     gnum_elt(iCell,1)=c_loc( local2Global(iCell))
      !>   enddo
      !>
      !> type (c_ptr), value :: n_elt
      !>   passer nCellDomain -> n_elt
      !>   type (c_ptr), value :: n_elt
      !>   integer                 :: nCellDomain
      !>   integer(c_int), pointer :: n_elt(:)
      !>
      !>   allocate(n_elt(1:1))
      !>   n_elt(1)=int(nCellDomain,kind=c_int)
      !>
      !> Appel fonction PDM_part_to_block_create
      !>   PDM_part_to_block_t = PDM_part_to_block_create( &
      !>   &    t_distrib=int(0,kind=c_int)               ,& !> 0
      !>   &    t_post=PDM_PART_TO_BLOCK_POST_CLEANUP     ,& !> 0
      !>   &    part_active_node=real(0,kind=c_double)    ,& !> 0
      !>   &    gnum_elt=c_loc(gnum_elt)                  ,& !> localToGlobal
      !>   &    weight=C_NULL_PTR                         ,& !> null
      !>   &    n_elt=c_loc(n_elt)                        ,& !> nombre element de ma partition
      !>   &    n_part=int(1,kind=c_int)                  ,& !> 1
      !>   &    fcomm=commSpace                            ) !> communicateur fortran


    end function PDM_part_to_block_create

    function PDM_part_to_block_create2 (t_distrib, &
                                        t_post, &
                                        part_active_node,&
                                        gnum_elt,&
                                        data_distrib_index,&
                                        weight,&
                                        n_elt,&
                                        n_part,&
                                        fcomm) &
                                        result(ptrC) &
      bind (c, name = 'PDM_part_to_block_create2_cf')

      use iso_c_binding

      implicit none

      integer(c_int), value :: t_distrib ! PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                         ! PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE,
                                         ! PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE
      integer(c_int), value :: t_post ! PDM_PART_TO_BLOCK_POST_NOTHING
                                      ! PDM_PART_TO_BLOCK_POST_CLEANUP
                                      ! PDM_PART_TO_BLOCK_POST_MERGE

      real (c_double), value :: part_active_node

      type (c_ptr), value :: gnum_elt
      type (c_ptr), value :: data_distrib_index
      type (c_ptr), value :: weight
      type (c_ptr), value :: n_elt
      integer(c_int), value :: n_part
      integer(c_int), value :: fComm
      type (c_ptr) :: ptrC

    end function PDM_part_to_block_create2

    !> \brief Return number of active ranks
    !!
    !! \param [in]   ptb          Part to block structure
    !!
    !!\return Number of active ranks

    function PDM_part_to_block_n_active_ranks_get (ptb)&
         result(n_active_ranks)&
         bind (c, name = 'PDM_part_to_block_n_active_ranks_get')
      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptb
      integer(c_int) :: n_active_ranks

    end function PDM_part_to_block_n_active_ranks_get

    !> \brief Return if current rank is active
    !!
    !! \param [in]   ptb          Part to block structure
    !!
    !! \return  if current rank is active

    function PDM_part_to_block_is_active_rank (ptb)&
         result(is_active_rank)&
         bind (c, name = 'PDM_part_to_block_is_active_rank')
      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptb
      integer(c_int) :: is_active_rank

    end function PDM_part_to_block_is_active_rank

    !> \brief Return active ranks
    !!
    !! \param [in]   ptb          Part to block structure
    !!
    !! \return  active ranks

    function PDM_part_to_block_active_ranks_get (ptb)&
         result(ptr_active_ranks)&
         bind (c, name = 'PDM_part_to_block_active_ranks_get')
      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptb
      type (c_ptr) :: ptr_active_ranks

    end function PDM_part_to_block_active_ranks_get

    !> \brief Return number of element in the current process
    !!
    !! \param [in]   ptb          Part to block structure
    !!
    !! \return Number of element in the current process

    function PDM_part_to_block_n_elt_block_get (ptb)&
         result(n_elt_block)&
         bind (c, name = 'PDM_part_to_block_n_elt_block_get')
      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptb
      integer(c_int) :: n_elt_block

    end function PDM_part_to_block_n_elt_block_get

    !> \brief Return global numbers of element in the current process
    !!
    !! \param [in]   ptb          Part to block structure
    !!
    !! \return  Global numbers

    function PDM_part_to_block_block_gnum_get (ptb)&
         result(ptr_block_gnum)&
         bind (c, name = 'PDM_part_to_block_block_gnum_get')
      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptb
      type (c_ptr) :: ptr_block_gnum

    end function PDM_part_to_block_block_gnum_get

    !> \brief Initialize a data exchange
    !!
    !! \param [in]   ptb          Part to block structure
    !! \param [in]   s_data       Data size
    !! \param [in]   t_stride     Stride type
    !! \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
    !! \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
    !!
    !! \return       Size of highest block

    function PDM_part_to_block_exch (ptb, s_data, t_stride, cst_stride, part_stride, part_data, &
                                     block_stride, block_data)&
         result(size_highest_block)&
         bind (c, name = 'PDM_part_to_block_exch')
      use iso_c_binding

      implicit none

      type (c_ptr), value       :: ptb
      integer (c_size_t), value :: s_data
      integer (c_int), value    :: t_stride !PDM_STRIDE_CST or PDM_STRIDE_VAR
      integer (c_int), value    :: cst_stride ! stride value if PDM_STRIDE_CST is selected
      type (c_ptr), value       :: part_stride
      type (c_ptr), value       :: part_data
      type (c_ptr)              :: block_stride
      type (c_ptr)              :: block_data

      integer(c_int)            :: size_highest_block

    end function PDM_part_to_block_exch

    !> \brief Free a part to block structure
    !!
    !! \param [inout] ptb         Part to block structure
    !!
    !! \return       NULL

    function PDM_part_to_block_free (ptb)&
         result(ptr_null)&
         bind (c, name = 'PDM_part_to_block_free')
      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptb
      type (c_ptr) :: ptr_null

    end function PDM_part_to_block_free

    !> \brief Return block distribution index
    !!
    !! \param [in] ptb         Part to block structure
    !!
    !! \return  Distribution (size = communicator size + 1)

    function PDM_part_to_block_distrib_index_get (ptb)&
         result(ptr_distrib_index)&
         bind (c, name = 'PDM_part_to_block_distrib_index_get')
      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptb
      type (c_ptr) :: ptr_distrib_index

    end function PDM_part_to_block_distrib_index_get

    !> \brief Return processus destination
    !!
    !! \param [in] ptb         Part to block structure
    !!
    !! \return  Destination (size = sum of partition elements)

    function PDM_part_to_block_destination_get(ptb)&
         result(ptr_block_dest)&
         bind (c, name = 'PDM_part_to_block_destination_get')
      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptb
      type (c_ptr) :: ptr_block_dest

    end function PDM_part_to_block_destination_get

  end interface

end module pdm_part_to_block
