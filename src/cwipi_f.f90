!-----------------------------------------------------------------------------
! This file is part of the CWIPI library. 
!
! Copyright (C) 2011  ONERA
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

module cwipi
  implicit none

  !
  ! Parameters
  ! ----------
  !
  ! cwipi_nature_t
  ! TODO:  cwipi_nature a supprimer ?
  integer, parameter :: cwipi_nature_element_center = 0
  integer, parameter :: cwipi_nature_node = 1
  !
  ! cwipi_coupling_type_t
  integer, parameter :: cwipi_cpl_parallel_with_part = 0
  integer, parameter :: cwipi_cpl_parallel_without_part = 1
  integer, parameter :: cwipi_cpl_sequential = 2
  !
  ! cwipi_type_t
  integer, parameter :: cwipi_type_float = 0
  integer, parameter :: cwipi_type_double = 1
  !
  ! cwipi_interpolation_t
  integer, parameter :: cwipi_interpolation_standard = 0
  integer, parameter :: cwipi_interpolation_user = 1
  !
  ! mesh type
  integer, parameter :: cwipi_static_mesh = 0
  integer, parameter :: cwipi_mobile_mesh = 1
  !
  ! solver type
  integer, parameter :: cwipi_solver_cell_center = 0
  integer, parameter :: cwipi_solver_cell_vertex = 1
  !
  ! exchange status
  integer, parameter :: cwipi_exchange_ok = 0
  integer, parameter :: cwipi_exchange_bad_receiving = 1

  !
  ! Logical unit for listing
  integer, save :: ifile

  !
  ! Public interfaces
  interface cwipi_exchange_f ; module procedure &
    cwipi_exch_without_user_itp_f_, &
    cwipi_exch_with_user_itp_f_
  end interface

  interface cwipi_send_f     ; module procedure  &
    cwipi_send_without_user_itp_f_, &
    cwipi_send_with_user_itp_f_
  end interface

  interface cwipi_issend_f     ; module procedure  &
    cwipi_issend_without_user_itp_f_, &
    cwipi_issend_with_user_itp_f_
  end interface

  interface cwipi_wait_issend_f     ; module procedure  &
    cwipi_wait_issend_f_
  end interface

  interface cwipi_wait_irecv_f     ; module procedure  &
    cwipi_wait_irecv_f_
  end interface

  interface cwipi_init_f                  ; module procedure &
    cwipi_init_f_
  end interface

  interface cwipi_set_output_listing_f    ; module procedure &
    cwipi_set_output_listing_f_
  end interface

  interface cwipi_add_loc_int_ctrl_param_f ; module procedure &
    cwipi_add_loc_int_ctrl_param_f_
  end interface

  interface cwipi_add_loc_dbl_ctrl_param_f ; module procedure &
    cwipi_add_loc_dbl_ctrl_param_f_
  end interface

  interface cwipi_add_loc_str_ctrl_param_f ; module procedure &
    cwipi_add_loc_str_ctrl_param_f_
  end interface

  interface cwipi_set_loc_int_ctrl_param_f ; module procedure &
    cwipi_set_loc_int_ctrl_param_f_
  end interface

  interface cwipi_set_loc_dbl_ctrl_param_f ; module procedure &
    cwipi_set_loc_dbl_ctrl_param_f_
  end interface

  interface cwipi_set_loc_str_ctrl_param_f ; module procedure &
    cwipi_set_loc_str_ctrl_param_f_
  end interface

  interface cwipi_get_loc_int_ctrl_param_f ; module procedure &
    cwipi_get_loc_int_ctrl_param_f_
  end interface

  interface cwipi_get_loc_dbl_ctrl_param_f ; module procedure &
    cwipi_get_loc_dbl_ctrl_param_f_
  end interface

  interface cwipi_get_loc_str_ctrl_param_f ; module procedure &
    cwipi_get_loc_str_ctrl_param_f_
  end interface

  interface cwipi_del_loc_int_ctrl_param_f ; module procedure &
    cwipi_del_loc_int_ctrl_param_f_
  end interface

  interface cwipi_del_loc_dbl_ctrl_param_f ; module procedure &
    cwipi_del_loc_dbl_ctrl_param_f_
  end interface

  interface cwipi_del_loc_str_ctrl_param_f ; module procedure &
    cwipi_del_loc_str_ctrl_param_f_
  end interface

  interface cwipi_get_dis_int_ctrl_param_f ; module procedure &
    cwipi_get_dis_int_ctrl_param_f_
  end interface

  interface cwipi_get_dis_dbl_ctrl_param_f ; module procedure &
    cwipi_get_dis_dbl_ctrl_param_f_
  end interface

  interface cwipi_get_dis_str_ctrl_param_f ; module procedure &
    cwipi_get_dis_str_ctrl_param_f_
  end interface

  interface cwipi_get_n_located_dist_pts_f ; module procedure &
    cwipi_get_n_located_dist_pts_f_
  end interface

  interface cwipi_synch_ctrl_param_f       ; module procedure &
    cwipi_synch_ctrl_param_f_
  end interface

  interface cwipi_create_coupling_f        ; module procedure &
   cwipi_create_coupling_f_
  end interface

  interface cwipi_set_points_to_locate_f   ; module procedure &
   cwipi_set_points_to_locate_f_
  end interface

  interface cwipi_define_mesh_f            ; module procedure &
   cwipi_define_mesh_f_
  end interface

  interface cwipi_add_polyhedra_f          ; module procedure &
    cwipi_add_polyhedra_f_
  end interface

  interface cwipi_locate_f                 ; module procedure &
    cwipi_locate_f_
  end interface

  interface cwipi_update_location_f                 ; module procedure &
    cwipi_update_location_f_
  end interface

  interface cwipi_get_bary_coord_f         ; module procedure &
    cwipi_get_bary_coord_f_
  end interface

  interface cwipi_get_bary_coord_idx_f     ; module procedure &
    cwipi_get_bary_coord_idx_f_
  end interface

  interface cwipi_get_location_f           ; module procedure &
    cwipi_get_location_f_
  end interface

  interface cwipi_get_distance_f           ; module procedure &
    cwipi_get_distance_f_
  end interface

  interface cwipi_get_coord_f           ; module procedure &
    cwipi_get_coord_f_
  end interface

  interface cwipi_receive_f                ; module procedure &
    cwipi_receive_f_
  end interface

  interface cwipi_ireceive_f                ; module procedure &
    cwipi_ireceive_f_
  end interface

  interface cwipi_delete_coupling_f        ; module procedure &
    cwipi_delete_coupling_f_
  end interface

  interface cwipi_get_not_located_pts_f ; module procedure &
    cwipi_get_not_located_pts_f_
  end interface

  interface cwipi_get_n_not_located_pts_f ; module procedure &
    cwipi_get_n_not_located_pts_f_
  end interface

  interface cwipi_get_n_located_pts_f ; module procedure &
    cwipi_get_n_located_pts_f_
  end interface

  interface cwipi_get_located_pts_f ; module procedure &
    cwipi_get_located_pts_f_
  end interface

  interface  cwipi_dump_appli_properties_f ; module procedure &
      cwipi_dump_appli_properties_f_
  end interface

  interface cwipi_get_elt_cont_f ; module procedure &
      cwipi_get_elt_cont_f_
  end interface

  interface  cwipi_get_elt_cont_n_vtx_f ; module procedure &
    cwipi_get_elt_cont_n_vtx_f_
  end interface

  interface  cwipi_get_elt_cont_vtx_f ; module procedure &
    cwipi_get_elt_cont_vtx_f_
  end interface

  interface  cwipi_get_elt_cont_vtx_coo_f ; module procedure &
    cwipi_get_elt_cont_vtx_coo_f_
  end interface

  interface  cwipi_get_elt_cont_bar_coo_f ; module procedure &
    cwipi_get_elt_cont_bar_coo_f_
  end interface

  interface  cwipi_get_elt_cont_MPI_rank_f ; module procedure &
    cwipi_get_elt_cont_MPI_rank_f_
  end interface

  interface  cwipi_exch_cellvtxfd_eltcont_f ; module procedure &
    cwipi_exch_cellvtxfd_eltcont_f_
  end interface

  interface  cwipi_send_cellvtxfd_eltcont_f ; module procedure &
    cwipi_send_cellvtxfd_eltcont_f_
  end interface

  interface  cwipi_recv_cellvtxfd_eltcont_f ; module procedure &
    cwipi_recv_cellvtxfd_eltcont_f_
  end interface

  interface  cwipi_exch_cellcenfd_eltcont_f ; module procedure &
    cwipi_exch_cellcenfd_eltcont_f_
  end interface

  interface  cwipi_send_cellcenfd_eltcont_f ; module procedure &
    cwipi_send_cellcenfd_eltcont_f_
  end interface

  interface  cwipi_recv_cellcenfd_eltcont_f ; module procedure &
    cwipi_recv_cellcenfd_eltcont_f_
  end interface

  interface  cwipi_finalize_f; module procedure &
    cwipi_finalize_f_
  end interface

  interface  cwipi_set_info_f; module procedure &
    cwipi_set_info_f_
  end interface

  interface cwipi_dist_located_pts_get_f; module procedure &
    cwipi_dist_located_pts_get_f_  
  end interface cwipi_dist_located_pts_get_f

  !
  ! Private

  private :: cwipi_init_f_,                   &
             cwipi_set_output_listing_f_,     &
             cwipi_exch_without_user_itp_f_,  &
             cwipi_exch_with_user_itp_f_,     &
             cwipi_send_without_user_itp_f_,  &
             cwipi_send_with_user_itp_f_,     &
             cwipi_add_loc_int_ctrl_param_f_, &
             cwipi_add_loc_dbl_ctrl_param_f_, &
             cwipi_add_loc_str_ctrl_param_f_, &
             cwipi_set_loc_int_ctrl_param_f_, &
             cwipi_set_loc_dbl_ctrl_param_f_, &
             cwipi_set_loc_str_ctrl_param_f_, &
             cwipi_get_loc_int_ctrl_param_f_, &
             cwipi_get_loc_dbl_ctrl_param_f_, &
             cwipi_get_loc_str_ctrl_param_f_, &
             cwipi_del_loc_int_ctrl_param_f_, &
             cwipi_del_loc_dbl_ctrl_param_f_, &
             cwipi_del_loc_str_ctrl_param_f_, &
             cwipi_get_dis_int_ctrl_param_f_, &
             cwipi_get_dis_dbl_ctrl_param_f_, &
             cwipi_get_dis_str_ctrl_param_f_, &
             cwipi_synch_ctrl_param_f_,       &
             cwipi_create_coupling_f_,        &
             cwipi_set_points_to_locate_f_,   &
             cwipi_define_mesh_f_,            &
             cwipi_add_polyhedra_f_,          &
             cwipi_locate_f_,                 &
             cwipi_update_location_f_,        &
             cwipi_receive_f_,                &
             cwipi_delete_coupling_f_,        &
             cwipi_dump_appli_properties_f_,  &
             cwipi_finalize_f_,               &
             cwipi_get_location_f_,           &
             cwipi_get_distance_f_,           &
             cwipi_get_bary_coord_f_,         &
             cwipi_get_bary_coord_idx_f_,     &
             cwipi_get_coord_f_,              &
             cwipi_get_not_located_pts_f_,    &
             cwipi_get_elt_cont_f_,           &
             cwipi_get_elt_cont_n_vtx_f_,     &
             cwipi_get_elt_cont_vtx_f_,       &
             cwipi_get_elt_cont_vtx_coo_f_,   &
             cwipi_get_elt_cont_bar_coo_f_,   &
             cwipi_get_elt_cont_MPI_rank_f_,  &
             cwipi_exch_cellvtxfd_eltcont_f_, &
             cwipi_exch_cellcenfd_eltcont_f_, &
             cwipi_dist_located_pts_get_f_,   &
             cwipi_set_info_f_

contains

 subroutine cwipi_finalize_f_

    implicit none

    call cwipi_finalize_cf

  end subroutine cwipi_finalize_f_


  subroutine cwipi_dump_appli_properties_f_

    implicit none

    call cwipi_dump_appli_properties_cf

  end subroutine cwipi_dump_appli_properties_f_


!
!*******************************************************************************
!
! cwipi_init_f_
!
!  Initialize the cwipi library.
!  Redirect outputs in a file (Standard output with output_listing = NULL or
!  output_logical_unit = -1)
!  Create the current communicator application from 'common_comm'.
!
!  parameters:
!    globalComm    <-- Common MPI communicator
!    appliName     <-- Current application name
!    appliComm     --> Internal MPI communicator for the current
!                      application
!
!  It is a synchronization point between all applications
!
!*******************************************************************************
!

  subroutine cwipi_init_f_ (globalComm, appliName, appliComm)

    implicit none

    integer :: globalcomm, applicomm
    character (len = *) :: appliname

    integer :: l1

    l1 = len(appliname)

    call cwipi_init_cf (globalcomm, appliname, l1, applicomm)

  end subroutine cwipi_init_f_

!
!*******************************************************************************
!
!  Set up the file used for the output listing
!
!  parameters:
!    output_listing      <-- Output listing file
!
!*******************************************************************************
!

  subroutine cwipi_set_output_listing_f_ (outputUnit)

  implicit none

  integer :: outputUnit

  ifile =  outputUnit

  call cwipi_set_output_listing_cf

  end subroutine cwipi_set_output_listing_f_

!
!********************************************************************************
! cwipi_add_loc_int_ctrl_param_f
!
! Add a integer control parameter
!
! parameters
!    name           <-- parameter name
!    initial_value  <-- initial value
!
!********************************************************************************
!

  subroutine cwipi_add_loc_int_ctrl_param_f_ (name, initialvalue)

    implicit none

    character (len = *) :: name
    integer :: initialvalue

    integer :: l

    l = len(name)

    call cwipi_add_loc_int_ctrl_param_cf (name, l, initialvalue)

  end subroutine cwipi_add_loc_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_add_loc_dbl_ctrl_param_f
!
! Add a double control parameter
!
! parameters
!    name           <-- parameter name
!    initial_value  <-- initial value
!********************************************************************************
!

  subroutine cwipi_add_loc_dbl_ctrl_param_f_ (name, initialvalue)

    implicit none

    character (len = *) ::name
    double precision :: initialvalue

    integer :: l

    l = len(name)

    call cwipi_add_loc_dbl_ctrl_param_cf (name, l, initialvalue)

  end subroutine cwipi_add_loc_dbl_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_add_loc_str_ctrl_param_f
!
! Add a double control parameter
!
! parameters
!    name           <-- parameter name
!    initial_value  <-- initial value
!********************************************************************************
!

  subroutine cwipi_add_loc_str_ctrl_param_f_(name, initialvalue)

    implicit none

    character (len = *) :: name
    character (len = *) :: initialvalue

    integer :: l1, l2

    l1 = len(name)
    l2 = len(initialvalue)

    call cwipi_add_loc_str_ctrl_param_cf(name, l1, initialvalue, l2)

  end subroutine cwipi_add_loc_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_set_loc_int_ctrl_param_f
!
!  Set a integer control parameter
!
!  parameters
!     name           <-- parameter name
!     value          <-- value
!
!********************************************************************************
!

  subroutine cwipi_set_loc_int_ctrl_param_f_(name, initialvalue)

    implicit none

    character (len = *) ::name
    integer :: initialvalue

    integer :: l

    l = len(name)

    call cwipi_set_loc_int_ctrl_param_cf (name, l, initialvalue)

  end subroutine cwipi_set_loc_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_set_loc_dbl_ctrl_param_f
!
! Set a double control parameter
!
! parameters
!    name           <-- parameter name
!    value          <-- value
!
!
!********************************************************************************
!

  subroutine cwipi_set_loc_dbl_ctrl_param_f_ (name, initialvalue)

    implicit none

    character (len = *) :: name
    double precision :: initialvalue

    integer :: l

    l = len(name)

    call cwipi_set_loc_dbl_ctrl_param_cf (name, l, initialvalue)

  end subroutine cwipi_set_loc_dbl_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_set_loc_str_ctrl_param_f
!
! Set a double control parameter
!
! parameters
!    name           <-- parameter name
!    value          <-- value
!
!
!********************************************************************************
!

  subroutine cwipi_set_loc_str_ctrl_param_f_ (name, initialvalue)

    implicit none

    character (len = *) :: name
    character (len = *) :: initialvalue

    integer :: l1, l2

    l1 = len(name)
    l2 = len(initialvalue)

    call cwipi_set_loc_str_ctrl_param_cf (name, l1, initialvalue, l2)

  end subroutine cwipi_set_loc_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_loc_int_ctrl_param_f
!
! Get a integer control parameter of the current application
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_loc_int_ctrl_param_f_ (name, value)

    implicit none

    character (len = *) :: name
    integer ::value

    integer :: l

    l = len(name)

    call cwipi_get_loc_int_ctrl_param_cf (name, l, value)

  end subroutine cwipi_get_loc_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_loc_dbl_ctrl_param_f
!
! Get a double control parameter of the current application
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_loc_dbl_ctrl_param_f_ (name, value)

    implicit none

    character (len = *) :: name
    double precision :: value

    integer :: l

    l = len(name)

    call cwipi_get_loc_dbl_ctrl_param_cf (name, l, value)

  end subroutine cwipi_get_loc_dbl_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_loc_str_ctrl_param_f
!
! Get a double control parameter of the current application
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_loc_str_ctrl_param_f_(name, value_str_f)

    implicit none

    character (len = *) :: name
    character (len = *) :: value_str_f

    integer :: l1, l2

    l1 = len(name)
    l2 = len(value_str_f)

    call cwipi_get_loc_str_ctrl_param_cf(name, l1, value_str_f, l2)

  end subroutine cwipi_get_loc_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_del_loc_int_ctrl_param_f
!
! Delete a current application int parameter
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_del_loc_int_ctrl_param_f_(name)

    implicit none

    character (len = *) :: name
    integer l

    l = len(name)

    call cwipi_del_loc_int_ctrl_param_cf(name, l)

  end subroutine cwipi_del_loc_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_del_loc_dbl_ctrl_param_f
!
! Delete a current application double parameter
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_del_loc_dbl_ctrl_param_f_ (name)

    implicit none

    character (len = *) :: name
    integer :: l

    l = len(name)

    call cwipi_del_loc_dbl_ctrl_param_cf (name, l)

  end subroutine cwipi_del_loc_dbl_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_del_loc_str_ctrl_param_f
!
! Delete a current application double parameter
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_del_loc_str_ctrl_param_f_ (name)

    implicit none

    character (len = *) :: name
    integer :: l

    l = len(name)

    call cwipi_del_loc_str_ctrl_param_cf (name, l)

  end subroutine cwipi_del_loc_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_dis_int_ctrl_param_f
!
! Get a integer control parameter of a other application
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_dis_int_ctrl_param_f_ (appliName, &
                                                             paramName, &
                                                             value)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    integer :: value

    integer :: l1, l2

    l1 = len(appliName)
    l2 = len(paramName)

    call cwipi_get_dis_int_ctrl_param_cf (appliName, &
                                                         l1, &
                                                         paramName, &
                                                         l2, &
                                                         value)

  end subroutine cwipi_get_dis_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_dis_dbl_ctrl_param_f
!
! Get a double control parameter of a other application
!
! parameters
!    application_name    <-- application name
!    name                <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_dis_dbl_ctrl_param_f_  (appliName, &
                                                               paramName, &
                                                               value)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    double precision :: value

    integer l1, l2

    l1 = len(appliName)
    l2 = len(paramName)

    call cwipi_get_dis_dbl_ctrl_param_cf (appliName, &
                                                            l1, &
                                                            paramName, &
                                                            l2, &
                                                            value)

  end subroutine cwipi_get_dis_dbl_ctrl_param_f_


!
!********************************************************************************
!
! cwipi_get_dis_dbl_ctrl_param_f
!
! Get a double control parameter of a other application
!
! parameters
!    application_name    <-- application name
!    name                <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_dis_str_ctrl_param_f_(appliName, &
                                             paramName, &
                                             value_str_f)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    character (len = *) :: value_str_f

    integer :: l1, l2, l3

    l1 = len(appliName)
    l2 = len(paramName)
    l3 = len(value_str_f)

    call cwipi_get_dis_str_ctrl_param_cf (appliName, &
                                          l1, &
                                          paramName, &
                                          l2, &
                                          value_str_f, &
                                          l3)

  end subroutine cwipi_get_dis_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_synchronize_control_parameter_f
!
! Synchronize local control parameters with an other application.
!  It is a synchronization point with this second application
!
! parameters
!    appliName           <-- application name
!
!********************************************************************************
!

  subroutine cwipi_synch_ctrl_param_f_ (appliName)

    implicit none

    character (len = *) :: appliName

    integer l

    l = len(appliName)

    call cwipi_synch_ctrl_param_cf (appliName, l)

  end subroutine cwipi_synch_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_dump_application_properties_f (define into cwipi_cf.hxx)
!
! Dump application properties
!
!********************************************************************************
!

!
!********************************************************************************
!
! cwipi_create_coupling_f
!
! Create a coupling object
!
! parameters:
!   couplingName            <-- Coupling identifier
!   couplingType            <-- Coupling type
!   cplAppli                <-- Coupled application name
!   entitiesDim             <-- Mesh entities dimension (1, 2 or 3)
!   tolerance               <-- Geometric tolerance to locate
!   meshT                   <-- CWIPI_STATIC_MESH
!                               CWIPI_MOBILE_MESH (not implemented yet)
!   solverT                 <-- CWIPI_SOLVER_CELL_CENTER
!                               CWIPI_SOLVER_CELL_VERTEX
!   outputFreq              <-- Output frequency
!   outputFmt               <-- Output format to visualize exchanged fields
!                               on the coupled mesh. Choice between :
!                                 - "EnSight Gold"
!                                 - "MED_ficher"
!                                 - "CGNS"
!   outputFmtOpt            <-- Output options
!                             text                output text files
!                             binary              output binary files (default)
!                             big_endian          force binary files
!                                                 to big-endian
!                             discard_polygons    do not output polygons
!                                                 or related values
!                             discard_polyhedra   do not output polyhedra
!                                                 or related values
!                             divide_polygons     tesselate polygons
!                                                 with triangles
!                             divide_polyhedra    tesselate polyhedra
!                                                 with tetrahedra and pyramids
!                                                 (adding a vertex near
!                                                 each polyhedron's center)
!
!********************************************************************************
!

  subroutine cwipi_create_coupling_f_ (couplingName, &
                                           couplingType, &
                                           cplAppli, &
                                           entitiesDim, &
                                           tolerance, &
                                           meshT, &
                                           solvert, &
                                           outputfreq, &
                                           outputfmt, &
                                           outputfmtopt)

    implicit none

    character (len = *) :: couplingName, cplAppli
    integer :: entitiesDim, meshT, solverT, couplingType
    double precision :: tolerance
    integer :: outputFreq
    character (len = *) :: outputFmt, outputFmtOpt

    integer :: lCouplingName, lCplAppli
    integer :: lOutputFmt, lOutputFmtOpt

    lCouplingName = len(couplingName)
    lCplAppli     = len(cplAppli)
    lOutputFmt    = len(outputFmt)
    lOutputFmtOpt = len(outputFmtOpt)

    call cwipi_create_coupling_cf(couplingName, &
                                  lCouplingName, &
                                  couplingType, &
                                  cplAppli, &
                                  lCplAppli, &
                                  entitiesDim, &
                                  tolerance, &
                                  meshT, &
                                  solverT, &
                                  outputFreq, &
                                  outputFmt, &
                                  lOutputFmt, &
                                  outputFmtOpt, &
                                  lOutputFmtOpt)
    
  end subroutine cwipi_create_coupling_f_

!
!********************************************************************************
!
! cwipi_set_points_to_locate_f
!
! Set points to locate. This function must be called if the points to locate
! do not correspond to :
!        - vertices for CELL_VERTEX nature
!        - cell center for CELL_CENTER nature
!
! parameters:
!   couplingName       <-- coupling identifier
!   nPts               <-- number of points to locate
!   coords             <-- coordinates of points to locate (enterlaced)
!
!********************************************************************************
!

  subroutine cwipi_set_points_to_locate_f_ (couplingName, &
                                                nPts, &
                                                coords)

    implicit none

    character (len = *) :: couplingName

    integer :: nPts
    double precision, dimension(3 * npts) :: coords

    integer :: lCouplingName

    lCouplingName  = len(couplingName)

    call cwipi_set_points_to_locate_cf(couplingName, &
                                           lCouplingName, &
                                           nPts, &
                                           coords)

  end subroutine cwipi_set_points_to_locate_f_

!
!********************************************************************************
!
! cwipi_define_mesh_f
!
!
! Define the support mesh for a coupling. The connectivity is sorted if
! necessary.
!
! Order definition :
!    1D : edges
!    2D : triangles, quadrangles, polygons
!    3D : tetrahedra, pyramids, prism, hexaedra
!
! Local connectivity for the following element type :
!
!  - edge :
!
!   1 x-------x 2
!
!  - triangle :
!
!   1 x-------x 3
!      \     /
!       \   /
!        \ /
!         x 2
!
!  - quadrangle :
!
!      4 x-------x 3
!       /       /
!      /       /
!   1 x-------x2
!
!   - tetrahedra :
!
!         x 4
!        /|\
!       / | \
!      /  |  \
!   1 x- -|- -x 3
!      \  |  /
!       \ | /
!        \|/
!         x 2
!
!   - pyramid :
!
!          5 x
!           /|\
!          //| \
!         // |  \
!      4 x/--|---x 3
!       //   |  /
!      //    | /
!   1 x-------x 2
!
!  - prism :
!
!   4 x-------x 6
!     |\     /|
!     | \   / |
!   1 x- \-/ -x 3
!      \ 5x  /
!       \ | /
!        \|/
!         x 2
!
!  -  hexahedra :
!
!      8 x-------x 7
!       /|      /|
!      / |     / |
!   5 x-------x6 |
!     | 4x----|--x 3
!     | /     | /
!     |/      |/
!   1 x-------x 2
!
!
!
! parameters:
!   couplingName       <-- coupling name
!   nVertex            <-- number of vertices
!   nElts              <-- number of elements
!   coords             <-- vertex interlaced coordinates
!   connecIndex        <-- element -> vertices index (O to n-1)
!                          size: n_elements + 1
!   connec             <-- element -> vertex connectivity (1 to n)
!                          size: connectivity_index[n_elements]
!
!********************************************************************************
!

  subroutine cwipi_define_mesh_f_ (couplingName, &
                                       nVertex, &
                                       nElts, &
                                       coords, &
                                       connecIndex, &
                                       connec)

    implicit none

    character (len = *) :: couplingName
    integer :: lCouplingName

    integer :: nElts, nVertex
    integer, dimension (nelts+1) :: connecIndex(nelts+1)

    integer, dimension (*) ::  connec
    double precision, dimension(3 * nVertex) :: coords

    lCouplingName    = len(couplingName)

    call cwipi_define_mesh_cf(couplingName, &
                                  lCouplingName, &
                                  nVertex, &
                                  nElts, &
                                  coords, &
                                  connecindex, &
                                  connec)

  end subroutine cwipi_define_mesh_f_

!
!********************************************************************************
!
! cwipi_add_polyhedra_f
!
! parameters:
!   couplingName       <-- Coupling name
!   nElts              <-- Number of elements
!   cellToFaceIdx      <-- Cell -> faces connectivity index (0 to n-1)
!                          size: nElts + 1
!   cellToFace         <-- Cell -> faces connectivity (1 to n)                         
!                          size: cellToFaceIdx(nElts)
!   nFaces             <-- Number of faces
!   faceConnecIdx      <-- Face -> vertex connectivity index (O to n-1)
!                          size: nFaces + 1
!   faceConnec         <-- Face -> vertex connectivity (1 to n)
!                          size: faceConnecIdx[nFaces]
!*******************************************************************************
!

  subroutine cwipi_add_polyhedra_f_ (couplingName, &
                                     nElts, &
                                     cellToFaceIdx, &
                                     cellToFace, &
                                     nFaces, &
                                     faceConnecIdx, &
                                     faceConnec)

    implicit none

    character (len = *) :: couplingName
    integer :: lCouplingname, nElts, nFaces
    integer, dimension(nelts) :: cellToFaceIdx
    integer, dimension(*) :: cellToFace, faceConnecIdx, faceConnec

    lCouplingName = len(couplingName)

    call cwipi_add_polyhedra_cf (couplingName, &
                                 lCouplingName, &
                                 nElts, &
                                 cellToFaceIdx, &
                                 cellToFace, &
                                 nFaces, &
                                 faceConnecIdx, &
                                 faceConnec)

  end subroutine cwipi_add_polyhedra_f_

!
!********************************************************************************
!
! cwipi_locate_f
!
! Location completion.
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine cwipi_locate_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname

    lCouplingName = len(couplingName)

    call cwipi_locate_cf(couplingName, lCouplingName)
  end subroutine cwipi_locate_f_

!
!********************************************************************************
!
! cwipi_dist_located_pts_get_f_
!
! Return located points to interface
!
! parameters
!   couplingName         <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine cwipi_dist_located_pts_get_f_(couplingName, distance)
    implicit none

    character (len = *) :: couplingName
    real(kind=4) :: distance(*)
    integer :: lcouplingname

    lCouplingName = len(couplingName)

    call cwipi_dist_located_pts_get_cf(couplingName, lCouplingName, distance)
  end subroutine cwipi_dist_located_pts_get_f_
!
!********************************************************************************
!
! cwipi_locate_f
!
! Location completion.
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine cwipi_update_location_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname

    lCouplingName = len(couplingName)

    call cwipi_update_location_cf(couplingName, lCouplingName)
  end subroutine cwipi_update_location_f_


!
!********************************************************************************
!
! cwipi_get_location_f
!
! Get located points location
!
! parameters
!   couplingName         <-- Coupling identifier
!   location             <-- Get located points location
!
!*******************************************************************************
!

  subroutine cwipi_get_location_f_(couplingName, location)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname
    integer, dimension(*) :: location

    lCouplingName = len(couplingName)

    call cwipi_get_distant_location_cf (couplingName, lCouplingName, location)
  end subroutine cwipi_get_location_f_

!
!********************************************************************************
!
! cwipi_get_location_f
!
! Get located points location
!
! parameters
!   couplingName         <-- Coupling identifier
!   location             <-- Get located points location
!
!*******************************************************************************
!

  subroutine cwipi_get_distance_f_(couplingName, distance)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname
    real*4, dimension(*) :: distance

    lCouplingName = len(couplingName)

    call cwipi_get_distant_distance_cf (couplingName, lCouplingName, distance)
  end subroutine cwipi_get_distance_f_

!
!********************************************************************************
!
! cwipi_get_coord
!
! Get distant points coordinates
!
! parameters
!   couplingName      <-- Coupling identifier
!   coordinates       <-- Get distant point coordinates
! 
!*******************************************************************************
!

  subroutine cwipi_get_coord_f_(couplingName, &
                                coordinates)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname
    double precision, dimension(*) :: coordinates

    lCouplingName = len(couplingName)

    call cwipi_get_dis_coord_cf(couplingName, &
                                lCouplingName, &
                                coordinates)
  end subroutine cwipi_get_coord_f_

!
!********************************************************************************
!
! cwipi_get_bary_coord_idx_f
!
! Get located points barycentric coordinates index
!
! parameters
!   couplingName                 <-- Coupling identifier
!   barycentricCoordinatesIndex  <-- Get located points barycentric coordinates
!                                    index
!*******************************************************************************
!

  subroutine cwipi_get_bary_coord_idx_f_(couplingName, &
                                         barycentricCoordinatesIndex)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname
    integer, dimension(*) :: barycentricCoordinatesIndex

    lCouplingName = len(couplingName)

    call cwipi_get_dis_bary_coord_idx_cf (couplingName, &
                                          lCouplingName, &
                                          barycentricCoordinatesIndex)
  end subroutine cwipi_get_bary_coord_idx_f_

!
!********************************************************************************
!
! cwipi_get_bary_coord_f
!
! Get located points barycentric coordinates
!
! parameters
!   couplingName              <-- Coupling identifier
!   barycentricCoordinates   <-- Get located points barycentric coordinates
!
!*******************************************************************************
!

  subroutine cwipi_get_bary_coord_f_(couplingName, &
                                     barycentricCoordinates)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname
    double precision, dimension(*) :: barycentricCoordinates

    lCouplingName = len(couplingName)

    call cwipi_get_dis_bary_coord_cf (couplingName, &
                                      lCouplingName, &
                                      barycentricCoordinates)
  end subroutine cwipi_get_bary_coord_f_

!
!********************************************************************************
!
! cwipi_exch_without_user_itp_f_
!
! Exchange data with the coupled application.
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field (NULL -> no sending)
!   receivingFieldName   <-- Receiving field name
!   receivingField       --> Receiving field
!   nNotLocatedPoints    --> Number of not located points
!   status               --> Cwipi exchange status
!
!********************************************************************************
!

  subroutine cwipi_exch_without_user_itp_f_ (couplingName, &
                                             exchangeName, &
                                             stride, &
                                             nStep, &
                                             timeValue, &
                                             sendingFieldName, &
                                             sendingField, &
                                             receivingFieldName, &
                                             receivingField, &
                                             nNotLocatedPoints, &
                                             status)

    implicit none

    character (len = *) :: couplingName, exchangeName, sendingFieldName
    character (len = *) :: receivingFieldName
    integer :: stride, nStep, status
    integer :: nNotLocatedPoints
    double precision :: timeValue
    double precision, dimension(*) :: sendingField, receivingField

    integer :: lCouplingName, lExchangeName, lSendingFieldName
    integer :: lReceivingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)
    lReceivingFieldName = len(receivingFieldName)

    call cwipi_exchange_cf(couplingName, &
                           lCouplingName, &
                           exchangeName, &
                           lExchangeName, &
                           stride, &
                           nStep, &
                           timeValue, &
                           sendingFieldName, &
                           lSendingFieldName, &
                           sendingField, &
                           receivingFieldName, &
                           lReceivingFieldName, &
                           receivingField, &
                           nNotLocatedPoints, &
                           status)

  end subroutine cwipi_exch_without_user_itp_f_

!
!********************************************************************************
!
! cwipi_issend_without_user_itp_f
!
! Send data to the coupled application (only send)
! Non-blocking communication 
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   tag                  <-- Exchange tag
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field
!   request              -->Sending request
!
!********************************************************************************
!

  subroutine cwipi_issend_without_user_itp_f_(couplingName, &
                                              exchangeName, &
                                              tag, &
                                              stride, &
                                              nStep, &
                                              timeValue, &
                                              sendingFieldName, &
                                              sendingField, &
                                              request)

    implicit none

    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer :: stride, nStep, status, tag, request
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call cwipi_issend_cf(couplingName, &
                         lCouplingName, &
                         exchangeName, &
                         lExchangeName, &
                         tag, &
                         stride, &
                         nStep, &
                         timeValue, &
                         sendingFieldName, &
                         lSendingFieldName, &
                         sendingField, &
                         request)

  end subroutine cwipi_issend_without_user_itp_f_

 
!********************************************************************************
!
! cwipi_issend_with_user_itp_f
!
! Exchange data with the coupled application (only send)
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field (NULL -> no sending)
!   ptInterpolationFct   <-- Callback for interpolation
!   request              --> Exchange request
!
!********************************************************************************
!
  subroutine cwipi_issend_with_user_itp_f_ (couplingName, &
                                            exchangeName, &
                                            tag, &
                                            stride, &
                                            nStep, &
                                            timeValue, &
                                            sendingFieldName, &
                                            sendingField, &
                                            ptInterpolationFct, &
                                            request)

    implicit none

    interface
       subroutine  ptInterpolationFct(entitiesDim, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer :: entitiesDim
         integer :: nLocalVertex
         integer :: nLocalElement
         integer :: nLocalPolyhedra
         integer :: nDistantPoint
         double precision, dimension(*) :: localCoordinates
         integer, dimension(*) :: localConnectivityIndex
         integer, dimension(*) :: localConnectivity
         integer, dimension(*) :: localPolyFaceIndex
         integer, dimension(*) :: localPolyCellToFaceConnec
         integer, dimension(*) :: localPolyFaceConnecIdx
         integer, dimension(*) :: localPolyFaceConnec
         double precision, dimension(*) :: disPtsCoordinates
         integer, dimension(*) :: disPtsLocation
         real*4, dimension(*) :: disPtsDistance
         integer, dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         integer :: stride
         integer :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine ptInterpolationFct
    end interface
    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer :: stride, nStep, tag, request
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call cwipi_issend_with_user_itp_cf(couplingName, &
                                       lCouplingName, &
                                       exchangeName, &
                                       lExchangeName, &
                                       tag, &
                                       stride, &
                                       nStep, &
                                       timeValue, &
                                       sendingFieldName, &
                                       lSendingFieldName, &
                                       sendingField, &
                                       ptInterpolationFct, &
                                       request)

  end subroutine cwipi_issend_with_user_itp_f_


!
!********************************************************************************
!
! cwipi_send_without_user_itp_f
!
! Exchange data with the coupled application (only send)
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field (NULL -> no sending)
!   status               --> Cwipi exchange status
!
!********************************************************************************
!

  subroutine cwipi_send_without_user_itp_f_ (couplingName, &
                                             exchangeName, &
                                             stride, &
                                             nStep, &
                                             timeValue, &
                                             sendingFieldName, &
                                             sendingField, &
                                             status)

    implicit none

    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer :: stride, nStep, status
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call cwipi_send_cf(couplingName, &
                       lCouplingName, &
                       exchangeName, &
                       lExchangeName, &
                       stride, &
                       nStep, &
                       timeValue, &
                       sendingFieldName, &
                       lSendingFieldName, &
                       sendingField, &
                       status)

  end subroutine cwipi_send_without_user_itp_f_

 
!********************************************************************************
!
! cwipi_send_with_user_itp_f
!
! Exchange data with the coupled application (only send)
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field (NULL -> no sending)
!   ptInterpolationFct   <-- Callback for interpolation
!   status               --> Cwipi exchange status
!
!********************************************************************************
!
  subroutine cwipi_send_with_user_itp_f_ (couplingName, &
                                          exchangeName, &
                                          stride, &
                                          nStep, &
                                          timeValue, &
                                          sendingFieldName, &
                                          sendingField, &
                                          ptInterpolationFct, &
                                          status)

    implicit none

    interface
       subroutine  ptInterpolationFct(entitiesDim, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer :: entitiesDim
         integer :: nLocalVertex
         integer :: nLocalElement
         integer :: nLocalPolyhedra
         integer :: nDistantPoint
         double precision, dimension(*) :: localCoordinates
         integer, dimension(*) :: localConnectivityIndex
         integer, dimension(*) :: localConnectivity
         integer, dimension(*) :: localPolyFaceIndex
         integer, dimension(*) :: localPolyCellToFaceConnec
         integer, dimension(*) :: localPolyFaceConnecIdx
         integer, dimension(*) :: localPolyFaceConnec
         double precision, dimension(*) :: disPtsCoordinates
         integer, dimension(*) :: disPtsLocation
         real*4, dimension(*) :: disPtsDistance
         integer, dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         integer :: stride
         integer :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine ptInterpolationFct
    end interface
    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer :: stride, nStep, status
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call cwipi_send_with_user_itp_cf(couplingName, &
                                        lCouplingName, &
                                        exchangeName, &
                                        lExchangeName, &
                                        stride, &
                                        nStep, &
                                        timeValue, &
                                        sendingFieldName, &
                                        lSendingFieldName, &
                                        sendingField, &
                                        ptInterpolationFct, &
                                        status)

  end subroutine cwipi_send_with_user_itp_f_

!
!********************************************************************************
!
! cwipi_receive_f
!
! Exchange data with the coupled application. (only receive)
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   receivingFieldName   <-- Receiving field name
!   receivingField       --> Receiving field
!   nNotLocatedPoints    --> Number of not located points
!   status               --> Cwipi exchange status
!
!********************************************************************************
!

  subroutine cwipi_receive_f_ (couplingName, &
                                   exchangeName, &
                                   stride, &
                                   nStep, &
                                   timeValue, &
                                   receivingFieldName, &
                                   receivingField, &
                                   nNotLocatedPoints, &
                                   status)

    implicit none

    character (len = *) :: couplingName, exchangeName, receivingFieldName
    integer :: stride, nStep, status
    integer :: nNotlocatedPoints
    double precision :: timeValue
    double precision, dimension(*) :: receivingField

    integer :: lCouplingName, lExchangeName, lReceivingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lReceivingFieldName = len(receivingFieldName)

    call cwipi_receive_cf(couplingName, &
                              lCouplingName, &
                              exchangeName, &
                              lExchangeName, &
                              stride, &
                              nStep, &
                              timeValue, &
                              receivingFieldName, &
                              lReceivingFieldName, &
                              receivingField, &
                              nNotLocatedPoints, &
                              status)

  end subroutine cwipi_receive_f_


!
!********************************************************************************
!
! cwipi_ireceive_f
!
! Receive data from the coupled application.
! Non-blocking communication. receivingField is fuly updated 
! after cwipi_wait_irecv calling
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   tag                  <-- Exchange tag
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   receivingFieldName   <-- Receiving field name
!   receivingField       <-- Receiving field
!   request              --> Exchange request
!
!********************************************************************************
!

  subroutine cwipi_ireceive_f_ (couplingName, &
                                exchangeName, &
                                tag, &
                                stride, &
                                nStep, &
                                timeValue, &
                                receivingFieldName, &
                                receivingField, &
                                request)

    implicit none

    character (len = *) :: couplingName, exchangeName, receivingFieldName
    integer :: stride, nStep, status, tag, request
    double precision :: timeValue
    double precision, dimension(*) :: receivingField

    integer :: lCouplingName, lExchangeName, lReceivingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lReceivingFieldName = len(receivingFieldName)

    call cwipi_ireceive_cf(couplingName, &
                           lCouplingName, &
                           exchangeName, &
                           lExchangeName, &
                           tag, &
                           stride, &
                           nStep, &
                           timeValue, &
                           receivingFieldName, &
                           lReceivingFieldName, &
                           receivingField, &
                           request)

  end subroutine cwipi_ireceive_f_


!
!********************************************************************************
!
! cwipi_wait_issend
!
! Wait for cwipi_issend
!
! parameters
!   couplingName         <-- Coupling identifier
!   request              <-- Exchange request
!
!********************************************************************************
!

  subroutine cwipi_wait_issend_f_(couplingName, &
                                  request)

    implicit none

    character (len = *) :: couplingName
    integer :: request
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_wait_issend_cf(couplingName, lCouplingName, request)

  end subroutine cwipi_wait_issend_f_


!
!********************************************************************************
!
! cwipi_wait_irecv
!
! Wait for cwipi_irecv
!
! parameters
!   couplingName         <-- Coupling identifier
!   request              <-- Exchange request
!
!********************************************************************************
!

  subroutine cwipi_wait_irecv_f_(couplingName, &
                                 request)

    implicit none

    character (len = *) :: couplingName
    integer :: request
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_wait_irecv_cf(couplingName, lCouplingName, request)

  end subroutine cwipi_wait_irecv_f_


!
!********************************************************************************
!
! cwipi_exchange_f_with_user_itp_f
!
! Exchange data with the coupled application.
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field (NULL -> no sending)
!   receivingFieldName   <-- Receiving field name
!   receivingField       --> Receiving field
!   ptInterpolationFct   <-- Callback for interpolation
!   nNotLocatedPoints    --> Number of not located points
!   status               --> Cwipi exchange status
!
!********************************************************************************
!

  subroutine cwipi_exch_with_user_itp_f_ (couplingName, &
                                          exchangeName, &
                                          exchangeDim, &
                                          nStep, &
                                          timeValue, &
                                          sendingFieldName, &
                                          sendingField, &
                                          receivingFieldName, &
                                          receivingField, &
                                          ptInterpolationFct, &
                                          nNotLocatedPoints, &
                                          status)

    implicit none

    interface
       subroutine  ptInterpolationFct(entitiesDim, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer :: entitiesDim
         integer :: nLocalVertex
         integer :: nLocalElement
         integer :: nLocalPolyhedra
         integer :: nDistantPoint
         double precision, dimension(*) :: localCoordinates
         integer, dimension(*) :: localConnectivityIndex
         integer, dimension(*) :: localConnectivity
         integer, dimension(*) :: localPolyFaceIndex
         integer, dimension(*) :: localPolyCellToFaceConnec
         integer, dimension(*) :: localPolyFaceConnecIdx
         integer, dimension(*) :: localPolyFaceConnec
         double precision, dimension(*) :: disPtsCoordinates
         integer, dimension(*) :: disPtsLocation
         real*4, dimension(*) :: disPtsDistance
         integer, dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         integer :: stride
         integer :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine ptInterpolationFct
    end interface

    character (len = *) :: couplingName, exchangeName, sendingFieldName
    character (len = *) :: receivingFieldName
    integer :: exchangeDim, nStep, status
    integer :: nnotlocatedpoints
    double precision :: timeValue
    double precision, dimension(*) ::  sendingField, receivingField

    integer :: lCouplingName, lExchangeName, lSendingFieldName
    integer :: lReceivingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)
    lReceivingFieldName = len(receivingFieldName)

    call cwipi_exch_with_user_itp_cf(couplingName, &
                                     lCouplingName, &
                                     exchangeName, &
                                     lExchangeName, &
                                     exchangeDim, &
                                     nStep, &
                                     timeValue, &
                                     sendingFieldName, &
                                     lSendingFieldName, &
                                     sendingField, &
                                     receivingFieldName, &
                                     lReceivingFieldName, &
                                     receivingField, &
                                     ptInterpolationFct, &
                                     nnotlocatedpoints, &
                                     status)

  end subroutine cwipi_exch_with_user_itp_f_

!
!********************************************************************************
!
! cwipi_delete_coupling_f
!
! Delete a coupling
!
! parameters
!   couplingName         <-- Coupling identifier
!
!
!********************************************************************************
!

  subroutine cwipi_delete_coupling_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_delete_coupling_cf(couplingName, lCouplingName)

  end subroutine cwipi_delete_coupling_f_

!
!********************************************************************************
!
! cwipi_finalize()  (define into cwipi_cf.hxx)
!
! Finalize cwipi. This is a synchronization point between all applications
!
!********************************************************************************
!

!
!********************************************************************************
!
! cwipi_get_not_located_pts_f
!
! Get located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   notLocatedPoints     --> Not located points (1 to n)
!
!********************************************************************************
!
  subroutine cwipi_get_not_located_pts_f_(couplingName, notLocatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer, dimension(*) :: notLocatedPoints
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_not_located_pts_cf(couplingName, lCouplingName, notLocatedPoints)
  end subroutine cwipi_get_not_located_pts_f_

!
!********************************************************************************
!
! cwipi_get_n_located_dist_pts_f
!
! Get located points
!
! parameters
!   couplingName             <-- Coupling identifier
!   nLocatedDistantPoints     --> Number oflocated distan points
!
!********************************************************************************
!
  subroutine cwipi_get_n_located_dist_pts_f_(couplingName, nLocatedDistantPoints)

    implicit none

    character (len = *) :: couplingName
    integer :: nLocatedDistantPoints
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_n_located_dist_pts_cf(couplingName, &
                                                   lCouplingName, &
                                                   nLocatedDistantPoints)
  end subroutine cwipi_get_n_located_dist_pts_f_



!
!********************************************************************************
!
! Get located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   locatedPoints        --> Located points (1 to n)
!
!********************************************************************************
!

  subroutine cwipi_get_located_pts_f_(couplingName, locatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer, dimension(*) :: locatedPoints
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_located_pts_cf(couplingName, lCouplingName, locatedPoints)
  end subroutine cwipi_get_located_pts_f_

!
!********************************************************************************
!
! Get number of located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   locatedPoints        --> Located points (1 to n)
!
!********************************************************************************
!

  subroutine cwipi_get_n_located_pts_f_(couplingName, nLocatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer :: nLocatedPoints
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_n_located_pts_cf(couplingName, &
                                         lCouplingName, &
                                         nLocatedPoints)

  end subroutine cwipi_get_n_located_pts_f_

!
!********************************************************************************
!
! Get number of not located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   locatedPoints        --> Located points (1 to n)
!
!********************************************************************************
!

  subroutine cwipi_get_n_not_located_pts_f_(couplingName, nNotLocatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer :: nNotLocatedPoints
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_n_not_located_pts_cf(couplingName, &
                                        lCouplingName, &
                                        nNotLocatedPoints)

  end subroutine cwipi_get_n_not_located_pts_f_

!
!********************************************************************************
!
! Get distant elements that contain located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   elements             --> Element that contain located points
!
!********************************************************************************
!

  subroutine cwipi_get_elt_cont_f_ (couplingName, elements)
    implicit none

    character (len = *) :: couplingName
    integer, dimension(*) :: elements
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_elt_cont_cf (couplingName, lCouplingName, elements)

  end subroutine cwipi_get_elt_cont_f_

!
!********************************************************************************
!
! Get number of vertices of distant elements that contain located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   n_vertices           --> Number of vertices of element that contain
!                            located pointlocatedPoints        --> Located points (1 to n)
!
!********************************************************************************
!

  subroutine cwipi_get_elt_cont_n_vtx_f_ (couplingName, n_vertices)
    implicit none

    character (len = *) :: couplingName
    integer, dimension(*) :: n_vertices
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_elt_cont_n_vtx_cf(couplingname, lCouplingName, n_vertices)

  end subroutine cwipi_get_elt_cont_n_vtx_f_

!
!********************************************************************************
!
! Get vertices id of distant elements that contain located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   vertices             --> Vertices id
!
!********************************************************************************
!

  subroutine cwipi_get_elt_cont_vtx_f_ (couplingName, vertices)
    implicit none

    character (len = *) :: couplingName
    integer, dimension(*) :: vertices
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_elt_cont_n_vtx_cf(couplingName, lCouplingName, vertices)

  end subroutine cwipi_get_elt_cont_vtx_f_

!
!********************************************************************************
!
! Get vertices coordinates of distant elements that contain located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   coordinates          --> Vertices coordinates
!
!********************************************************************************
!

  subroutine cwipi_get_elt_cont_vtx_coo_f_ (couplingName, coordinates)
    implicit none

    character (len = *) :: couplingName
    double precision, dimension(*) :: coordinates
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_elt_cont_vtx_coo_cf(couplingName, lCouplingName, coordinates)

  end subroutine cwipi_get_elt_cont_vtx_coo_f_

!
!********************************************************************************
!
! Get barycentric coords in distant elements for located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   coordinates          --> Barycentric coordinates
!
!********************************************************************************
!

  subroutine cwipi_get_elt_cont_bar_coo_f_(couplingName, coordinates)
    implicit none

    character (len = *) :: couplingName
    double precision, dimension(*) :: coordinates
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_elt_cont_bar_coo_cf(couplingName, lCouplingName, coordinates)

  end subroutine cwipi_get_elt_cont_bar_coo_f_

!
!********************************************************************************
!
! For each located point get the MPI rank of distant element
!
! parameters
!   couplingName         <-- Coupling identifier
!   MPIranks             --> MPI ranks that contains located point
!
!********************************************************************************
!

  subroutine cwipi_get_elt_cont_MPI_rank_f_(couplingName, MPIrank)
    implicit none

    character (len = *) :: couplingName
    integer, dimension(*) :: MPIrank
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_elt_cont_MPI_rank_cf(couplingName, lCouplingName, MPIrank)

  end subroutine cwipi_get_elt_cont_MPI_rank_f_

!
!********************************************************************************
!
! Exchange Fields on vertices of element containing each located point
!
! parameters
!   couplingName         <-- Coupling identifier
!   sendingField         <-- Field defined on local mesh vertices
!   receivingField       --> Field defined on vertices of distant
!                             elements that contain each located point
!   stride               <-- Number of field component
!
!********************************************************************************
!

  subroutine cwipi_exch_cellvtxfd_eltcont_f_(couplingName, &
                                              sendingField, &
                                              receivingField, &
                                              stride)
    implicit none

    character (len = *) :: couplingName
    double precision, dimension(*) :: sendingField
    double precision, dimension(*) :: receivingField
    integer :: stride
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_exch_cellvtxfd_eltcont_cf(couplingName, &
                                         lCouplingName, &
                                         sendingField, &
                                         receivingField, &
                                         stride)
  end subroutine cwipi_exch_cellvtxfd_eltcont_f_


  subroutine cwipi_send_cellvtxfd_eltcont_f_(couplingName, &
                                              sendingField, &
                                              stride)
    implicit none

    character (len = *) :: couplingName
    double precision, dimension(*) :: sendingField
    integer :: stride
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_send_cellvtxfd_eltcont_cf(couplingName, &
                                         lCouplingName, &
                                         sendingField, &
                                         stride)
  end subroutine cwipi_send_cellvtxfd_eltcont_f_

  subroutine cwipi_recv_cellvtxfd_eltcont_f_(couplingName, &
                                              receivingField, &
                                              stride)
    implicit none

    character (len = *) :: couplingName
    double precision, dimension(*) :: receivingField
    integer :: stride
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_recv_cellvtxfd_eltcont_cf(couplingName, &
                                         lCouplingName, &
                                         receivingField, &
                                         stride)
  end subroutine cwipi_recv_cellvtxfd_eltcont_f_
!
!********************************************************************************
!
! Get number of not located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   sendingField         <-- Field defined on local mesh vertices
!   receivingField       --> Field defined on vertices of distant
!                            elements that contain each located point
!   stride               <-- Number of field component
!
!********************************************************************************
!

  subroutine cwipi_exch_cellcenfd_eltcont_f_(couplingName, &
                                              sendingField, &
                                              receivingField, &
                                              stride)

    implicit none

    character (len = *) :: couplingName
    double precision, dimension(*) :: sendingField
    double precision, dimension(*) :: receivingField
    integer :: stride
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_exch_cellcenfd_eltcont_cf(couplingName, &
                                         lCouplingName, &
                                         sendingField, &
                                         receivingField, &
                                         stride)
  end subroutine cwipi_exch_cellcenfd_eltcont_f_

  subroutine cwipi_send_cellcenfd_eltcont_f_(couplingName, &
                                              sendingField, &
                                              stride)
    implicit none

    character (len = *) :: couplingName
    double precision, dimension(*) :: sendingField
    integer :: stride
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_send_cellcenfd_eltcont_cf(couplingName, &
                                         lCouplingName, &
                                         sendingField, &
                                         stride)
  end subroutine cwipi_send_cellcenfd_eltcont_f_

  subroutine cwipi_recv_cellcenfd_eltcont_f_(couplingName, &
                                              receivingField, &
                                              stride)
    implicit none

    character (len = *) :: couplingName
    double precision, dimension(*) :: receivingField
    integer :: stride
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_recv_cellcenfd_eltcont_cf(couplingName, &
                                         lCouplingName, &
                                         receivingField, &
                                         stride)
  end subroutine cwipi_recv_cellcenfd_eltcont_f_
!
!********************************************************************************
!
! Set coupling info
!
! parameters
!   couplingName         <-- Coupling identifier
!   sendingField         <-- Field defined on local mesh vertices
!   info                 <-- Coupling info
!
!********************************************************************************
!

 subroutine cwipi_set_info_f_(couplingName, info)

    implicit none

    character (len = *) :: couplingName
    integer :: info
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_set_info_cf(couplingName, &
                           lCouplingName, &
                           info)

  end subroutine cwipi_set_info_f_

end module cwipi
