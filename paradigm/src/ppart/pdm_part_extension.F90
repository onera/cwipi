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

module pdm_part_extension

  use pdm

  implicit none

  integer, parameter :: PDM_EXTEND_FROM_FACE = 0
  integer, parameter :: PDM_EXTEND_FROM_EDGE = 1
  integer, parameter :: PDM_EXTEND_FROM_VTX  = 2

interface

!>
!!
!! \brief Compute a part extension structure
!!
!! \param [in]   part_ext          PDM_part_extension_t
!!
!!

subroutine PDM_part_extension_compute (part_ext) &
bind (c, name='PDM_part_extension_compute')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: part_ext

end subroutine PDM_part_extension_compute



!>
!!
!! \brief Free a part extension structure
!!
!! \param [in]   part_ext          PDM_part_extension_t
!!
!!

subroutine PDM_part_extension_free (part_ext) &
bind (c, name='PDM_part_extension_free')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: part_ext

end subroutine PDM_part_extension_free

end interface


contains


!>
!!
!! \brief Initialize a part extension structure
!!
!! \param [out]   part_ext          Pointer to a new \ref PDM_part_extension_t object
!! \param [in]    n_part            Number of partitions
!! \param [in]    extend_type       Type of extension
!! \param [in]    depth             Depth of extension
!! \param [in]    comm              MPI communicator
!! \param [in]    owner             Ownership
!!
!!

subroutine PDM_part_extension_create (part_ext,    &
                                      n_domain,    &
                                      n_part,      &
                                      extend_type, &
                                      depth,       &
                                      comm,        &
                                      owner)
  use iso_c_binding
  implicit none

  type(c_ptr)                   :: part_ext
  integer, intent(in)           :: n_domain
  integer(pdm_l_num_s), pointer :: n_part(:)
  integer, intent(in)           :: extend_type
  integer, intent(in)           :: depth
  integer, intent(in)           :: comm
  integer, intent(in)           :: owner

  integer(c_int)                :: c_comm

  interface
    function PDM_part_extension_create_c (n_domain,    &
                                          n_part,      &
                                          extend_type, &
                                          depth,       &
                                          comm,        &
                                          owner)       &
    result (part_ext)                                  &
    bind (c, name='PDM_part_extension_create')
      use iso_c_binding
      implicit none

      integer(c_int), value :: n_domain
      type(c_ptr),    value :: n_part
      integer(c_int), value :: extend_type
      integer(c_int), value :: depth
      integer(c_int), value :: comm
      integer(c_int), value :: owner
      type(c_ptr)           :: part_ext

    end function PDM_part_extension_create_c
  end interface

  c_comm = PDM_MPI_Comm_f2c(comm)

  part_ext = PDM_part_extension_create_c (n_domain,      &
                                          c_loc(n_part), &
                                          extend_type,   &
                                          depth,         &
                                          c_comm,        &
                                          owner)

end subroutine PDM_part_extension_create


!>
!! \brief Set
!!
!! \param [in]  part_ext                  Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_domain                  Id of the domain
!! \param [in]  i_part                    Id of the partition
!! \param [in]  n_cell                    Number of cells
!! \param [in]  n_face                    Number of faces
!! \param [in]  n_face_part_bound         Number of partition boundary faces
!! \param [in]  n_face_group              Number of face groups
!! \param [in]  n_edge                    Number of edges
!! \param [in]  n_vtx                     Number of vertices
!! \param [in]  cell_face_idx             Cell-face connectivity index (size = \ref n_cell + 1)
!! \param [in]  cell_face                 Cell-face connectivity (size = \ref cell_face_idx(\ref n_cell + 1))
!! \param [in]  face_cell                 Face-cell connectivity (size = 2 * \ref n_face)
!! \param [in]  face_edge_idx             Face-edge connectivity index (size = \ref n_face + 1)
!! \param [in]  face_edge                 Face-edge connectivity (size = \ref face_edge_idx(\ref n_face + 1))
!! \param [in]  face_vtx_idx              Face-vertex connectivity index (size = \ref n_face + 1)
!! \param [in]  face_vtx                  Face-vertex connectivity (size = \ref face_vtx_idx(\ref n_face + 1))
!! \param [in]  edge_vtx                  Edge-vertex connectivity (size = 2 * \ref n_edge)
!! \param [in]  face_bound_idx
!! \param [in]  face_bound
!! \param [in]  face_join_idx
!! \param [in]  face_join
!! \param [in]  face_part_bound_proc_idx
!! \param [in]  face_part_bound_part_idx
!! \param [in]  face_part_bound
!! \param [in]  vtx_part_bound_proc_idx
!! \param [in]  vtx_part_bound_part_idx
!! \param [in]  vtx_part_bound
!! \param [in]  cell_ln_to_gn             Cell global ids (size = \ref n_cell)
!! \param [in]  face_ln_to_gn             Face global ids (size = \ref n_face)
!! \param [in]  edge_ln_to_gn             Edge global ids (size = \ref n_edge)
!! \param [in]  vtx_ln_to_gn              Vertex global ids (size = \ref n_vtx)
!! \param [in]  face_group_ln_to_gn
!! \param [in]  vtx_coord                 Vertex coordinates (size = 3 * \ref n_vtx)
!!

subroutine PDM_part_extension_set_part (part_ext,                 &
                                        i_domain,                 &
                                        i_part,                   &
                                        n_cell,                   &
                                        n_face,                   &
                                        n_face_part_bound,        &
                                        n_face_group,             &
                                        n_edge,                   &
                                        n_vtx,                    &
                                        cell_face_idx,            &
                                        cell_face,                &
                                        face_cell,                &
                                        face_edge_idx,            &
                                        face_edge,                &
                                        face_vtx_idx,             &
                                        face_vtx,                 &
                                        edge_vtx,                 &
                                        face_bound_idx,           &
                                        face_bound,               &
                                        face_join_idx,            &
                                        face_join,                &
                                        face_part_bound_proc_idx, &
                                        face_part_bound_part_idx, &
                                        face_part_bound,          &
                                        vtx_part_bound_proc_idx,  &
                                        vtx_part_bound_part_idx,  &
                                        vtx_part_bound,           &
                                        cell_ln_to_gn,            &
                                        face_ln_to_gn,            &
                                        edge_ln_to_gn,            &
                                        vtx_ln_to_gn,             &
                                        face_group_ln_to_gn,      &
                                        vtx_coord)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: part_ext
  integer, intent(in)           :: i_domain
  integer, intent(in)           :: i_part
  integer, intent(in)           :: n_cell
  integer, intent(in)           :: n_face
  integer, intent(in)           :: n_face_part_bound
  integer, intent(in)           :: n_face_group
  integer, intent(in)           :: n_edge
  integer, intent(in)           :: n_vtx
  integer(pdm_l_num_s), pointer :: cell_face_idx(:)
  integer(pdm_l_num_s), pointer :: cell_face(:)
  integer(pdm_l_num_s), pointer :: face_cell(:)
  integer(pdm_l_num_s), pointer :: face_edge_idx(:)
  integer(pdm_l_num_s), pointer :: face_edge(:)
  integer(pdm_l_num_s), pointer :: face_vtx_idx(:)
  integer(pdm_l_num_s), pointer :: face_vtx(:)
  integer(pdm_l_num_s), pointer :: edge_vtx(:)
  integer(pdm_l_num_s), pointer :: face_bound_idx(:)
  integer(pdm_l_num_s), pointer :: face_bound(:)
  integer(pdm_l_num_s), pointer :: face_join_idx(:)
  integer(pdm_l_num_s), pointer :: face_join(:)
  integer(pdm_l_num_s), pointer :: face_part_bound_proc_idx(:)
  integer(pdm_l_num_s), pointer :: face_part_bound_part_idx(:)
  integer(pdm_l_num_s), pointer :: face_part_bound(:)
  integer(pdm_l_num_s), pointer :: vtx_part_bound_proc_idx(:)
  integer(pdm_l_num_s), pointer :: vtx_part_bound_part_idx(:)
  integer(pdm_l_num_s), pointer :: vtx_part_bound(:)
  integer(pdm_g_num_s), pointer :: cell_ln_to_gn(:)
  integer(pdm_g_num_s), pointer :: face_ln_to_gn(:)
  integer(pdm_g_num_s), pointer :: edge_ln_to_gn(:)
  integer(pdm_g_num_s), pointer :: vtx_ln_to_gn(:)
  integer(pdm_g_num_s), pointer :: face_group_ln_to_gn(:)
  double precision,     pointer :: vtx_coord(:,:)

  interface
    subroutine PDM_part_extension_set_part_c (part_ext,                 &
                                              i_domain,                 &
                                              i_part,                   &
                                              n_cell,                   &
                                              n_face,                   &
                                              n_face_part_bound,        &
                                              n_face_group,             &
                                              n_edge,                   &
                                              n_vtx,                    &
                                              cell_face_idx,            &
                                              cell_face,                &
                                              face_cell,                &
                                              face_edge_idx,            &
                                              face_edge,                &
                                              face_vtx_idx,             &
                                              face_vtx,                 &
                                              edge_vtx,                 &
                                              face_bound_idx,           &
                                              face_bound,               &
                                              face_join_idx,            &
                                              face_join,                &
                                              face_part_bound_proc_idx, &
                                              face_part_bound_part_idx, &
                                              face_part_bound,          &
                                              vtx_part_bound_proc_idx,  &
                                              vtx_part_bound_part_idx,  &
                                              vtx_part_bound,           &
                                              cell_ln_to_gn,            &
                                              face_ln_to_gn,            &
                                              edge_ln_to_gn,            &
                                              vtx_ln_to_gn,             &
                                              face_group_ln_to_gn,      &
                                              vtx_coord)                &
    bind (c, name='PDM_part_extension_set_part')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: part_ext
      integer(c_int), value :: i_domain
      integer(c_int), value :: i_part
      integer(c_int), value :: n_cell
      integer(c_int), value :: n_face
      integer(c_int), value :: n_face_part_bound
      integer(c_int), value :: n_face_group
      integer(c_int), value :: n_edge
      integer(c_int), value :: n_vtx
      type(c_ptr),    value :: cell_face_idx
      type(c_ptr),    value :: cell_face
      type(c_ptr),    value :: face_cell
      type(c_ptr),    value :: face_edge_idx
      type(c_ptr),    value :: face_edge
      type(c_ptr),    value :: face_vtx_idx
      type(c_ptr),    value :: face_vtx
      type(c_ptr),    value :: edge_vtx
      type(c_ptr),    value :: face_bound_idx
      type(c_ptr),    value :: face_bound
      type(c_ptr),    value :: face_join_idx
      type(c_ptr),    value :: face_join
      type(c_ptr),    value :: face_part_bound_proc_idx
      type(c_ptr),    value :: face_part_bound_part_idx
      type(c_ptr),    value :: face_part_bound
      type(c_ptr),    value :: vtx_part_bound_proc_idx
      type(c_ptr),    value :: vtx_part_bound_part_idx
      type(c_ptr),    value :: vtx_part_bound
      type(c_ptr),    value :: cell_ln_to_gn
      type(c_ptr),    value :: face_ln_to_gn
      type(c_ptr),    value :: edge_ln_to_gn
      type(c_ptr),    value :: vtx_ln_to_gn
      type(c_ptr),    value :: face_group_ln_to_gn
      type(c_ptr),    value :: vtx_coord

    end subroutine PDM_part_extension_set_part_c
  end interface

  call PDM_part_extension_set_part_c (part_ext,                        &
                                      i_domain,                        &
                                      i_part,                          &
                                      n_cell,                          &
                                      n_face,                          &
                                      n_face_part_bound,               &
                                      n_face_group,                    &
                                      n_edge,                          &
                                      n_vtx,                           &
                                      c_loc(cell_face_idx),            &
                                      c_loc(cell_face),                &
                                      c_loc(face_cell),                &
                                      c_loc(face_edge_idx),            &
                                      c_loc(face_edge),                &
                                      c_loc(face_vtx_idx),             &
                                      c_loc(face_vtx),                 &
                                      c_loc(edge_vtx),                 &
                                      c_loc(face_bound_idx),           &
                                      c_loc(face_bound),               &
                                      c_loc(face_join_idx),            &
                                      c_loc(face_join),                &
                                      c_loc(face_part_bound_proc_idx), &
                                      c_loc(face_part_bound_part_idx), &
                                      c_loc(face_part_bound),          &
                                      c_loc(vtx_part_bound_proc_idx),  &
                                      c_loc(vtx_part_bound_part_idx),  &
                                      c_loc(vtx_part_bound),           &
                                      c_loc(cell_ln_to_gn),            &
                                      c_loc(face_ln_to_gn),            &
                                      c_loc(edge_ln_to_gn),            &
                                      c_loc(vtx_ln_to_gn),             &
                                      c_loc(face_group_ln_to_gn),      &
                                      c_loc(vtx_coord))

end subroutine PDM_part_extension_set_part


!>
!!
!! \brief Get connectivity
!!
!! \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_domain     Id of current domain
!! \param [in]  i_part       Id of current partition
!! \param [in]  mesh_entity  Type of mesh entity
!! \param [out] n_elt         Number of elements
!! \param [out] connect      Entity->group graph (size = \ref connect_idx[\ref n_elt])
!! \param [out] connect_idx  Index for entity->group graph (size = \ref n_elt + 1)
!!
!!

subroutine PDM_part_extension_connectivity_get (part_ext,          &
                                                i_domain,          &
                                                i_part,            &
                                                connectivity_type, &
                                                n_elt,             &
                                                connect,           &
                                                connect_idx)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: part_ext
  integer, intent(in)           :: i_domain
  integer, intent(in)           :: i_part
  integer, intent(in)           :: connectivity_type
  integer, intent(out)          :: n_elt
  integer(pdm_l_num_s), pointer :: connect(:)
  integer(pdm_l_num_s), pointer :: connect_idx(:)

  type(c_ptr)                   :: c_connect     = C_NULL_PTR
  type(c_ptr)                   :: c_connect_idx = C_NULL_PTR

  interface
    function PDM_part_extension_connectivity_get_c (part_ext,          &
                                                    i_domain,          &
                                                    i_part,            &
                                                    connectivity_type, &
                                                    connect,           &
                                                    connect_idx)       &
    result (n_elt)                                                     &
    bind (c, name='PDM_part_extension_connectivity_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: part_ext
      integer(c_int), value :: i_domain
      integer(c_int), value :: i_part
      integer(c_int), value :: connectivity_type
      type(c_ptr)           :: connect
      type(c_ptr)           :: connect_idx
      integer(c_int)        :: n_elt

    end function PDM_part_extension_connectivity_get_c
  end interface

  n_elt =  PDM_part_extension_connectivity_get_c (part_ext,          &
                                                  i_domain,          &
                                                  i_part,            &
                                                  connectivity_type, &
                                                  c_connect,         &
                                                  c_connect_idx)

  call c_f_pointer(c_connect_idx, &
                   connect_idx,   &
                   [n_elt+1])

  call c_f_pointer(c_connect, &
                   connect,   &
                   [connect_idx(n_elt + 1)])

end subroutine PDM_part_extension_connectivity_get


!>
!!
!! \brief Get global ids
!!
!! \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_domain     Id of current domain
!! \param [in]  i_part       Id of current partition
!! \param [in]  mesh_entity  Type of mesh entity
!! \param [out] n_elt        Number of elements
!! \param [out] ln_to_gn     Global ids (size = \ref n_elt)
!!

subroutine PDM_part_extension_ln_to_gn_get (part_ext,    &
                                            i_domain,    &
                                            i_part,      &
                                            mesh_entity, &
                                            n_elt,       &
                                            ln_to_gn)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: part_ext
  integer, intent(in)           :: i_domain
  integer, intent(in)           :: i_part
  integer, intent(in)           :: mesh_entity
  integer, intent(out)          :: n_elt
  integer(pdm_g_num_s), pointer :: ln_to_gn(:)

  type(c_ptr)                   :: c_ln_to_gn = C_NULL_PTR

  interface
    function PDM_part_extension_ln_to_gn_get_c (part_ext,    &
                                                i_domain,    &
                                                i_part,      &
                                                mesh_entity, &
                                                ln_to_gn)    &
    result (n_elt)                                            &
    bind (c, name='PDM_part_extension_ln_to_gn_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: part_ext
      integer(c_int), value :: i_domain
      integer(c_int), value :: i_part
      integer(c_int), value :: mesh_entity
      type(c_ptr)           :: ln_to_gn
      integer(c_int)        :: n_elt

    end function PDM_part_extension_ln_to_gn_get_c
  end interface

  n_elt =  PDM_part_extension_ln_to_gn_get_c (part_ext,    &
                                              i_domain,    &
                                              i_part,      &
                                              mesh_entity, &
                                              c_ln_to_gn)

  call c_f_pointer(c_ln_to_gn, &
                   ln_to_gn,   &
                   [n_elt])

end subroutine PDM_part_extension_ln_to_gn_get


!>
!!
!! \brief Get groups
!!
!! \param [in]  part_ext       Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_domain       Id of current domain
!! \param [in]  i_part         Id of current partition
!! \param [in]  mesh_entity    Type of mesh entity
!! \param [out] n_elt          Number of elements
!! \param [out] elt_group      Entity->group graph (size = \ref elt_group_idx[\ref n_elt])
!! \param [out] elt_group_idx  Index for entity->group graph (size = \ref n_elt + 1)
!! \param [out] ln_to_gn       Global ids (size = \ref elt_group_idx[\ref n_elt])
!!

subroutine PDM_part_extension_group_get (part_ext,      &
                                         i_domain,      &
                                         i_part,        &
                                         n_elt,         &
                                         elt_group,     &
                                         elt_group_idx, &
                                         ln_to_gn)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: part_ext
  integer, intent(in)           :: i_domain
  integer, intent(in)           :: i_part
  integer, intent(out)          :: n_elt
  integer(pdm_l_num_s), pointer :: elt_group(:)
  integer(pdm_l_num_s), pointer :: elt_group_idx(:)
  integer(pdm_g_num_s), pointer :: ln_to_gn(:)

  type(c_ptr)                   :: c_elt_group     = C_NULL_PTR
  type(c_ptr)                   :: c_elt_group_idx = C_NULL_PTR
  type(c_ptr)                   :: c_ln_to_gn      = C_NULL_PTR

  interface
    function PDM_part_extension_group_get_c (part_ext,      &
                                             i_domain,      &
                                             i_part,        &
                                             elt_group,     &
                                             elt_group_idx, &
                                             ln_to_gn)      &
    result (n_elt)                                          &
    bind (c, name='PDM_part_extension_group_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: part_ext
      integer(c_int), value :: i_domain
      integer(c_int), value :: i_part
      type(c_ptr)           :: elt_group
      type(c_ptr)           :: elt_group_idx
      type(c_ptr)           :: ln_to_gn
      integer(c_int)        :: n_elt

    end function PDM_part_extension_group_get_c
  end interface

  n_elt =  PDM_part_extension_group_get_c (part_ext,        &
                                           i_domain,        &
                                           i_part,          &
                                           c_elt_group,     &
                                           c_elt_group_idx, &
                                           c_ln_to_gn)

  call c_f_pointer(c_elt_group_idx, &
                   elt_group_idx,   &
                   [n_elt+1])

  call c_f_pointer(c_elt_group, &
                   elt_group,   &
                   [elt_group_idx(n_elt + 1)])

  call c_f_pointer(c_ln_to_gn, &
                   ln_to_gn,   &
                   [elt_group_idx(n_elt + 1)])

end subroutine PDM_part_extension_group_get



!>
!!
!! \brief Get vertex coordinates
!!
!! \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_domain     Id of current domain
!! \param [in]  i_part       Id of current partition
!! \param [out] n_vtx        Number of vertices
!! \param [out] vtx_coord    Vertex coordinates (size = \ref n_vtx * 3)
!!

subroutine PDM_part_extension_coord_get (part_ext,    &
                                         i_domain,    &
                                         i_part,      &
                                         n_vtx,       &
                                         vtx_coord)
  use iso_c_binding
  implicit none

  type(c_ptr), value        :: part_ext
  integer, intent(in)       :: i_domain
  integer, intent(in)       :: i_part
  integer, intent(out)      :: n_vtx
  double precision, pointer :: vtx_coord(:,:)

  type(c_ptr)               :: c_vtx_coord = C_NULL_PTR

  interface
    function PDM_part_extension_coord_get_c (part_ext,    &
                                             i_domain,    &
                                             i_part,      &
                                             vtx_coord)   &
    result (n_vtx)                                        &
    bind (c, name='PDM_part_extension_coord_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: part_ext
      integer(c_int), value :: i_domain
      integer(c_int), value :: i_part
      type(c_ptr)           :: vtx_coord
      integer(c_int)        :: n_vtx

    end function PDM_part_extension_coord_get_c
  end interface

  n_vtx =  PDM_part_extension_coord_get_c (part_ext,    &
                                           i_domain,    &
                                           i_part,      &
                                           c_vtx_coord)

  call c_f_pointer(c_vtx_coord, &
                   vtx_coord,   &
                   [3,n_vtx])

end subroutine PDM_part_extension_coord_get


!>
!! \brief Set a domain interface
!!
!! \param [in]  part_ext        Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_interface     Id of the interface
!! \param [in]  interface_kind  Kind of the interface
!! \param [in]  n_pair          Number of pairs
!! \param [in]  interface_ids   Global ids (size = 2 * \ref n_pair)
!! \param [in]  interface_dom   Ids of the domains (size = 2 * \ref n_pair or 2)
!!

! subroutine PDM_part_extension_domain_interface_set (part_ext,       &
!                                                     i_interface,    &
!                                                     interface_kind, &
!                                                     n_pair,         &
!                                                     interface_ids,  &
!                                                     interface_dom)
!   use iso_c_binding
!   implicit none

!   type(c_ptr), value            :: part_ext
!   integer                       :: i_interface
!   integer                       :: interface_kind
!   integer                       :: n_pair
!   integer(pdm_g_num_s), pointer :: interface_ids(:)
!   integer(pdm_l_num_s), pointer :: interface_dom(:)

!   type(c_ptr)                   :: c_interface_ids = C_NULL_PTR
!   type(c_ptr)                   :: c_interface_dom = C_NULL_PTR

!   interface
!     subroutine PDM_part_extension_domain_interface_set_c (part_ext,       &
!                                                           i_interface,    &
!                                                           interface_kind, &
!                                                           n_pair,         &
!                                                           interface_ids,  &
!                                                           interface_dom)  &
!     bind (c, name='PDM_part_extension_domain_interface_set')
!       use iso_c_binding
!       implicit none

!       type(c_ptr),    value :: part_ext
!       integer(c_int), value :: i_interface
!       integer(c_int), value :: interface_kind
!       integer(c_int), value :: n_pair
!       type(c_ptr),    value :: interface_ids
!       type(c_ptr),    value :: interface_dom

!     end subroutine PDM_part_extension_domain_interface_set_c
!   end interface

!   c_interface_ids = c_loc(interface_ids)
!   c_interface_dom = c_loc(interface_dom)

!   call PDM_part_extension_domain_interface_set_c (part_ext,         &
!                                                   i_interface,      &
!                                                   interface_kind,   &
!                                                   n_pair,           &
!                                                   c_interface_ids,  &
!                                                   c_interface_dom)

! end subroutine PDM_part_extension_domain_interface_set

!>
!! \brief Set the translation vector for a periodic interface
!!
!! \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_interface  Id of the interface
!! \param [in]  vect         Translation vector (size = 3)
!!

! subroutine PDM_part_extension_domain_interface_translation_set (part_ext,    &
!                                                                 i_interface, &
!                                                                 vect)
!   use iso_c_binding
!   implicit none

!   type(c_ptr), value        :: part_ext
!   integer                   :: i_interface
!   double precision, pointer :: vect(:)

!   type(c_ptr)               :: c_vect = C_NULL_PTR

!   interface
!     subroutine PDM_part_extension_domain_interface_translation_set_c (part_ext,    &
!                                                                       i_interface, &
!                                                                       vect)        &
!     bind (c, name='PDM_part_extension_domain_interface_translation_set')
!       use iso_c_binding
!       implicit none

!       type(c_ptr),    value :: part_ext
!       integer(c_int), value :: i_interface
!       type(c_ptr),    value :: vect

!     subroutine PDM_part_extension_domain_interface_translation_set_c
!   end interface

!   c_vect = c_loc(vect)

!   call PDM_part_extension_domain_interface_translation_set_c (part_ext,    &
!                                                               i_interface, &
!                                                               c_vect)

! end subroutine PDM_part_extension_domain_interface_translation_set


!>
!! \brief Set the rotation for a periodic interface
!!
!! \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_interface  Id of the interface
!! \param [in]  direction    Rotation axis (size = 3)
!! \param [in]  center       Rotation center (size = 3)
!! \param [in]  angle        Rotation angle
!!

! subroutine PDM_part_extension_domain_interface_rotation_set (part_ext,    &
!                                                              i_interface, &
!                                                              direction,   &
!                                                              center,      &
!                                                              angle)
!   use iso_c_binding
!   implicit none

!   type(c_ptr), value        :: part_ext
!   integer                   :: i_interface
!   double precision, pointer :: direction(:)
!   double precision, pointer :: center(:)
!   double precision          :: angle

!   type(c_ptr)               :: c_direction = C_NULL_PTR
!   type(c_ptr)               :: c_center    = C_NULL_PTR

!   interface
!     subroutine PDM_part_extension_domain_interface_rotation_set_c (part_ext,    &
!                                                                    i_interface, &
!                                                                    direction,   &
!                                                                    center,      &
!                                                                    angle)       &
!     bind (c, name='PDM_part_extension_domain_interface_rotation_set')
!       use iso_c_binding
!       implicit none

!       type(c_ptr),    value :: part_ext
!       integer(c_int), value :: i_interface
!       type(c_ptr),    value :: direction
!       type(c_ptr),    value :: center
!       real(c_double), value :: vect

!     subroutine PDM_part_extension_domain_interface_rotation_set_c
!   end interface

!   c_direction = c_loc(direction)
!   c_center    = c_loc(center)

!   call PDM_part_extension_domain_interface_rotation_set_c (part_ext,    &
!                                                            i_interface, &
!                                                            c_direction, &
!                                                            c_center,    &
!                                                            angle)

! end subroutine PDM_part_extension_domain_interface_rotation_set



end module pdm_part_extension
