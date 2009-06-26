module couplings
  implicit none

  !
  ! Parameters
  ! ----------
  !
  ! couplings_nature_t
  integer, parameter :: couplings_nature_element_center = 0
  integer, parameter :: couplings_nature_node = 1
  !
  ! couplings_type_t
  integer, parameter :: couplings_type_float = 0
  integer, parameter :: couplings_type_double = 1
  !
  ! couplings_interpolation_t
  integer, parameter :: couplings_interpolation_standard = 0
  integer, parameter :: couplings_interpolation_user = 1
  !
  ! couplings_mpi_ranks_for_coupling
  integer, parameter :: couplings_mpi_ranks_all_ranks = 0
  integer, parameter :: couplings_mpi_ranks_only_master = 1
  !
  ! mesh type
  integer, parameter :: couplings_static_mesh = 0
  integer, parameter :: couplings_mobile_mesh = 1
  !
  ! solver type
  integer, parameter :: couplings_solver_cell_center = 0
  integer, parameter :: couplings_solver_cell_vertex = 1
  !
  ! exchange status
  integer, parameter :: couplings_exchange_ok = 0
  integer, parameter :: couplings_exchange_bad_receiving = 1

  !
  ! Logical unit for listing
  integer, save :: ifile

  !
  ! Public interfaces
  interface couplings_exchange_f ; module procedure &
    couplings_exchange_without_user_interpolation_f_, &
    couplings_exchange_with_user_interpolation_f_
  end interface
  interface couplings_send_f     ; module procedure  &
    couplings_send_without_user_interpolation_f_, &
    couplings_send_with_user_interpolation_f_
  end interface

  interface couplings_init_f ; module procedure couplings_init_f_ ; end interface

  interface couplings_set_output_listing_f ; module procedure &
    couplings_set_output_listing_f_
  end interface

  interface couplings_add_local_int_control_parameter_f ; module procedure &
    couplings_add_local_int_control_parameter_f_
  end interface

  interface couplings_add_local_double_control_parameter_f ; module procedure &
    couplings_add_local_double_control_parameter_f_
  end interface

  interface couplings_set_local_int_control_parameter_f ; module procedure &
    couplings_set_local_int_control_parameter_f_
  end interface

  interface couplings_set_local_double_control_parameter_f ; module procedure &
    couplings_set_local_double_control_parameter_f_
  end interface

  interface couplings_get_local_int_control_parameter_f ; module procedure &
    couplings_get_local_int_control_parameter_f_
  end interface

  interface couplings_get_local_double_control_parameter_f ; module procedure &
    couplings_get_local_double_control_parameter_f_
  end interface

  interface couplings_delete_local_int_control_parameter_f ; module procedure &
    couplings_delete_local_int_control_parameter_f_
  end interface

  interface couplings_delete_local_double_control_parameter_f ; module procedure &
    couplings_delete_local_double_control_parameter_f_
  end interface

  interface couplings_get_distant_int_control_parameter_f ; module procedure &
    couplings_get_distant_int_control_parameter_f_
  end interface

  interface couplings_get_distant_double_control_parameter_f ; module procedure &
    couplings_get_distant_double_control_parameter_f_
  end interface

  interface couplings_get_n_located_distant_points_f ; module procedure &
    couplings_get_n_located_distant_points_f_
  end interface
  interface couplings_synchronize_control_parameter_f ; module procedure couplings_synchronize_control_parameter_f_ ; end interface
  interface couplings_create_coupling_f               ; module procedure couplings_create_coupling_f_               ; end interface
  interface couplings_set_points_to_locate_f          ; module procedure couplings_set_points_to_locate_f_          ; end interface
  interface couplings_define_mesh_f                   ; module procedure couplings_define_mesh_f_                   ; end interface
  interface couplings_add_polyhedra_f                 ; module procedure couplings_add_polyhedra_f_                 ; end interface
  interface couplings_locate_f                        ; module procedure couplings_locate_f_                        ; end interface
  interface couplings_get_barycentric_coordinates_f   ; module procedure couplings_get_barycentric_coordinates_f_   ; end interface
  interface couplings_get_barycentric_coordinates_index_f; module procedure &
    couplings_get_barycentric_coordinates_index_f_
  end interface
  interface couplings_get_location_f                  ; module procedure couplings_get_location_f_                  ; end interface
  interface couplings_receive_f                       ; module procedure couplings_receive_f_                       ; end interface
  interface couplings_delete_coupling_f               ; module procedure couplings_delete_coupling_f_               ; end interface
  interface couplings_get_not_located_points_f        ; module procedure couplings_get_not_located_points_f_        ; end interface
  interface couplings_get_n_not_located_points_f ; module procedure &
    couplings_get_n_not_located_points_f_
  end interface
  interface couplings_get_n_located_points_f ; module procedure &
    couplings_get_n_located_points_f_
  end interface
  !
  ! Private

  private :: couplings_init_f_,                                 &
             couplings_set_output_listing_f_,                   &
             couplings_exchange_without_user_interpolation_f_,  &
             couplings_exchange_with_user_interpolation_f_,     &
             couplings_send_without_user_interpolation_f_,      &
             couplings_send_with_user_interpolation_f_,         &
             couplings_add_local_int_control_parameter_f_,      &
             couplings_add_local_double_control_parameter_f_,   &
             couplings_set_local_int_control_parameter_f_,      &
             couplings_set_local_double_control_parameter_f_,   &
             couplings_get_local_int_control_parameter_f_,      &
             couplings_get_local_double_control_parameter_f_,   &
             couplings_delete_local_int_control_parameter_f_,   &
             couplings_delete_local_double_control_parameter_f_,&
             couplings_get_distant_int_control_parameter_f_,    &
             couplings_get_distant_double_control_parameter_f_, &
             couplings_synchronize_control_parameter_f_,        &
             couplings_create_coupling_f_,                      &
             couplings_set_points_to_locate_f_,                 &
             couplings_define_mesh_f_,                          &
             couplings_add_polyhedra_f_,                        &
             couplings_locate_f_,                               &
             couplings_receive_f_,                              &
             couplings_delete_coupling_f_,                      &
             couplings_get_not_located_points_f_

contains

!
!*******************************************************************************
!
! couplings_init_f_
!
!  Initialize the couplings library.
!  Redirect outputs in a file (Standard output with output_listing = NULL or
!  output_logical_unit = -1)
!  Create the current communicator application from 'common_comm'.
!
!  parameters:
!    globalComm    <-- Common MPI communicator
!    outputunit    <-- Output listing logical unit
!    appliName     <-- Current application name
!    appliComm     --> Internal MPI communicator for the current
!                      application
!
!  It is a synchronization point between all applications
!
!*******************************************************************************
!

  subroutine couplings_init_f_ (globalComm, outputUnit, appliName, mpiRanks, appliComm)

    implicit none

    integer :: globalcomm, outputUnit, mpiRanks, applicomm
    character (len = *) :: appliname

    integer :: l1
    integer :: i

    l1 = len(appliname)
    ifile =  outputUnit

    call couplings_init_cf (globalcomm, outputunit, appliname, l1, mpiRanks,  applicomm)

  end subroutine couplings_init_f_

!
!*******************************************************************************
!
!  Set up the file used for the output listing
!
!  parameters:
!    output_listing      <-- Output listing file (C function)
!
!*******************************************************************************
!

  subroutine couplings_set_output_listing_f_ (outputUnit)

  implicit none

  integer :: outputUnit

  ifile =  outputUnit

  call couplings_set_output_listing_cf (outputUnit)

  end subroutine couplings_set_output_listing_f_

!
!********************************************************************************
! couplings_add_local_int_control_parameter_f
!
! Add a integer control parameter
!
! parameters
!    name           <-- parameter name
!    initial_value  <-- initial value
!
!********************************************************************************
!

  subroutine couplings_add_local_int_control_parameter_f_ (name, initialvalue)

    implicit none

    character (len = *) :: name
    integer :: initialvalue

    integer :: l

    l = len(name)

    call couplings_add_local_int_control_parameter_cf (name, l, initialvalue)

  end subroutine couplings_add_local_int_control_parameter_f_

!
!********************************************************************************
!
! couplings_add_local_double_control_parameter_f
!
! Add a double control parameter
!
! parameters
!    name           <-- parameter name
!    initial_value  <-- initial value
!********************************************************************************
!

  subroutine couplings_add_local_double_control_parameter_f_ (name, initialvalue)

    implicit none

    character (len = *) ::name
    double precision :: initialvalue

    integer :: l

    l = len(name)

    call couplings_add_local_double_control_parameter_cf (name, l, initialvalue)

  end subroutine couplings_add_local_double_control_parameter_f_

!
!********************************************************************************
!
! couplings_set_local_int_control_parameter_f
!
!  Set a integer control parameter
!
!  parameters
!     name           <-- parameter name
!     value          <-- value
!
!********************************************************************************
!

  subroutine couplings_set_local_int_control_parameter_f_(name, initialvalue)

    implicit none

    character (len = *) ::name
    integer :: initialvalue

    integer :: l

    l = len(name)

    call couplings_set_local_int_control_parameter_cf (name, l, initialvalue)

  end subroutine couplings_set_local_int_control_parameter_f_

!
!********************************************************************************
!
! couplings_set_local_double_control_parameter_f
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

  subroutine couplings_set_local_double_control_parameter_f_ (name, initialvalue)

    implicit none

    character (len = *) :: name
    double precision :: initialvalue

    integer :: l

    l = len(name)

    call couplings_set_local_double_control_parameter_cf (name, l, initialvalue)

  end subroutine couplings_set_local_double_control_parameter_f_

!
!********************************************************************************
!
! couplings_get_local_int_control_parameter_f
!
! Get a integer control parameter of the current application
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine couplings_get_local_int_control_parameter_f_ (name, value)

    implicit none

    character (len = *) :: name
    integer ::value

    integer :: l

    l = len(name)

    call couplings_get_local_int_control_parameter_cf (name, l, value)

  end subroutine couplings_get_local_int_control_parameter_f_

!
!********************************************************************************
!
! couplings_get_local_double_control_parameter_f
!
! Get a double control parameter of the current application
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine couplings_get_local_double_control_parameter_f_ (name, value)

    implicit none

    character (len = *) :: name
    double precision :: value

    integer :: l

    l = len(name)

    call couplings_get_local_double_control_parameter_cf (name, l, value)

  end subroutine couplings_get_local_double_control_parameter_f_

!
!********************************************************************************
!
! couplings_delete_local_int_control_parameter_f
!
! Delete a current application int parameter
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine couplings_delete_local_int_control_parameter_f_ (name)

    implicit none

    character (len = *) :: name
    integer l

    l = len(name)

    call couplings_delete_local_int_control_parameter_cf (name, l)

  end subroutine couplings_delete_local_int_control_parameter_f_

!
!********************************************************************************
!
! couplings_delete_local_double_control_parameter_f
!
! Delete a current application double parameter
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine couplings_delete_local_double_control_parameter_f_ (name)

    implicit none

    character (len = *) :: name
    integer :: l

    l = len(name)

    call couplings_delete_local_double_control_parameter_cf (name, l)

  end subroutine couplings_delete_local_double_control_parameter_f_

!
!********************************************************************************
!
! couplings_get_distant_int_control_parameter_f
!
! Get a integer control parameter of a other application
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine couplings_get_distant_int_control_parameter_f_ (appliName, &
                                                             paramName, &
                                                             value)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    integer :: value

    integer :: l1, l2

    l1 = len(appliName)
    l2 = len(paramName)

    call couplings_get_distant_int_control_parameter_cf (appliName, &
                                                         l1, &
                                                         paramName, &
                                                         l2, &
                                                         value)

  end subroutine couplings_get_distant_int_control_parameter_f_

!
!********************************************************************************
!
! couplings_get_distant_double_control_parameter_f
!
! Get a double control parameter of a other application
!
! parameters
!    application_name    <-- application name
!    name                <-- parameter name
!
!********************************************************************************
!

  subroutine couplings_get_distant_double_control_parameter_f_  (appliName, &
                                                               paramName, &
                                                               value)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    double precision :: value

    integer l1, l2

    l1 = len(appliName)
    l2 = len(paramName)

    call couplings_get_distant_double_control_parameter_cf (appliName, &
                                                            l1, &
                                                            paramName, &
                                                            l2, &
                                                            value)

  end subroutine couplings_get_distant_double_control_parameter_f_

!
!********************************************************************************
!
! couplings_synchronize_control_parameter_f
!
! Synchronize local control parameters with an other application.
!  It is a synchronization point with this second application
!
! parameters
!    appliName           <-- application name
!
!********************************************************************************
!

  subroutine couplings_synchronize_control_parameter_f_ (appliName)

    implicit none

    character (len = *) :: appliName

    integer l

    l = len(appliName)

    call couplings_synchronize_control_parameter_cf (appliName, l)

  end subroutine couplings_synchronize_control_parameter_f_

!
!********************************************************************************
!
! couplings_dump_application_properties_f (define into couplings_cf.hxx)
!
! Dump application properties
!
!********************************************************************************
!

!
!********************************************************************************
!
! couplings_create_coupling_f
!
! Create a coupling object
!
! parameters:
!   couplingName            <-- Coupling identifier
!   cplAppli                <-- Coupled application name
!   entitiesDim             <-- Mesh entities dimension (1, 2 or 3)
!   tolerance               <-- Geometric tolerance to locate
!   meshT                   <-- COUPLINGS_STATIC_MESH
!                               COUPLINGS_MOBILE_MESH (not implemented yet)
!   solverT                 <-- COUPLINGS_SOLVER_CELL_CENTER
!                               COUPLINGS_SOLVER_CELL_VERTEX
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

  subroutine couplings_create_coupling_f_ (couplingName, &
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
    integer :: entitiesDim, meshT, solverT
    double precision :: tolerance
    integer :: outputFreq
    character (len = *) :: outputFmt, outputFmtOpt

    integer :: lCouplingName, lCplAppli
    integer :: lOutputFmt, lOutputFmtOpt

    lCouplingName = len(couplingName)
    lCplAppli     = len(cplAppli)
    lOutputFmt    = len(outputFmt)
    lOutputFmtOpt = len(outputFmtOpt)

    call couplings_create_coupling_cf(couplingName, &
                                      lCouplingName, &
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

  end subroutine couplings_create_coupling_f_

!
!********************************************************************************
!
! couplings_set_points_to_locate_f
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

  subroutine couplings_set_points_to_locate_f_ (couplingName, &
                                                nPts, &
                                                coords)

    implicit none

    character (len = *) :: couplingName

    integer :: nPts
    double precision, dimension(3 * npts) :: coords

    integer :: lCouplingName

    lCouplingName  = len(couplingName)

    call couplings_set_points_to_locate_cf(couplingName, &
                                           lCouplingName, &
                                           nPts, &
                                           coords)

  end subroutine couplings_set_points_to_locate_f_

!
!********************************************************************************
!
! couplings_define_mesh_f
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
!  -  hexaedra :
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
!   connecIndex        <-> element -> vertices index (O to n-1)
!                          size: n_elements + 1
!                          (out : stored connectivity_index)
!   connec             <-> element -> vertex connectivity (1 to n)
!                          size: connectivity_index[n_elements]
!                          (out : stored connectivity)
!
!********************************************************************************
!

  subroutine couplings_define_mesh_f_ (couplingName, &
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

    call couplings_define_mesh_cf(couplingName, &
                                  lCouplingName, &
                                  nVertex, &
                                  nElts, &
                                  coords, &
                                  connecindex, &
                                  connec)

  end subroutine couplings_define_mesh_f_

!
!********************************************************************************
!
! couplings_add_polyhedra_f
!
!*******************************************************************************
!

  subroutine couplings_add_polyhedra_f_ (couplingName, &
                                         nElts, &
                                         faceIdx, &
                                         cellToFace, &
                                         faceConnecIdx, &
                                         faceConnec)

    implicit none

    character (len = *) :: couplingName
    integer :: lCouplingname, nElts
    integer, dimension(nelts) :: faceIdx
    integer, dimension(*) :: cellToFace, faceConnecIdx, faceConnec

    lCouplingName = len(couplingName)

    call couplings_add_polyhedra_cf (couplingName, &
                                     lCouplingName, &
                                     nElts, &
                                     faceIdx, &
                                     cellToFace, &
                                     faceConnecIdx, &
                                     faceConnec)

  end subroutine couplings_add_polyhedra_f_

!
!********************************************************************************
!
! couplings_locate_f
!
! Location completion.
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine couplings_locate_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname

    lCouplingName = len(couplingName)

    call couplings_locate_cf(couplingName, lCouplingName)
  end subroutine couplings_locate_f_

!
!********************************************************************************
!
! couplings_get_location_f
!
! Get located points location
!
! parameters
!   couplingName         <-- Coupling identifier
!   location             <-- Get located points location
!
!*******************************************************************************
!

  subroutine couplings_get_location_f_(couplingName, location)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname
    integer, dimension(*) :: location

    lCouplingName = len(couplingName)

    call couplings_get_distant_location_cf (couplingName, lCouplingName, location)
  end subroutine couplings_get_location_f_

!
!********************************************************************************
!
! couplings_get_barycentric_coordinates_index_f
!
! Get located points barycentric coordinates index
!
! parameters
!   couplingName                 <-- Coupling identifier
!   barycentricCoordinatesIndex  <-- Get located points barycentric coordinates
!                                    index
!*******************************************************************************
!

  subroutine couplings_get_barycentric_coordinates_index_f_(couplingName, &
                                                            barycentricCoordinatesIndex)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname
    integer, dimension(*) :: barycentricCoordinatesIndex

    lCouplingName = len(couplingName)

    call couplings_get_distant_barycentric_coordinates_index_cf (couplingName, &
                                                         lCouplingName, &
                                                         barycentricCoordinatesIndex)
  end subroutine couplings_get_barycentric_coordinates_index_f_

!
!********************************************************************************
!
! couplings_get_barycentric_coordinates_f
!
! Get located points barycentric coordinates
!
! parameters
!   couplingName              <-- Coupling identifier
!   barycentricCoordinates   <-- Get located points barycentric coordinates
!
!*******************************************************************************
!

  subroutine couplings_get_barycentric_coordinates_f_(couplingName, &
                                                      barycentricCoordinates)

    implicit none

    character (len = *) :: couplingName
    integer :: lcouplingname
    double precision, dimension(*) :: barycentricCoordinates

    lCouplingName = len(couplingName)

    call couplings_get_distant_barycentric_coordinates_cf (couplingName, &
                                                   lCouplingName, &
                                                   barycentricCoordinates)
  end subroutine couplings_get_barycentric_coordinates_f_

!
!********************************************************************************
!
! couplings_exchange_without_user_interpolation_f_
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
!   status               --> Couplings exchange status
!
!********************************************************************************
!

  subroutine couplings_exchange_without_user_interpolation_f_ (couplingName, &
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

    call couplings_exchange_cf(couplingName, &
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

  end subroutine couplings_exchange_without_user_interpolation_f_

!
!********************************************************************************
!
! couplings_send_without_user_interpolation_f
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
!   status               --> Couplings exchange status
!
!********************************************************************************
!

  subroutine couplings_send_without_user_interpolation_f_ (couplingName, &
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

    call couplings_send_cf(couplingName, &
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

  end subroutine couplings_send_without_user_interpolation_f_

!
!********************************************************************************
!
! couplings_send_with_user_interpolation_f
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
!   status               --> Couplings exchange status
!
!********************************************************************************
!

  subroutine couplings_send_with_user_interpolation_f_ (couplingName, &
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
                                      stored, &
                                      localParentEltNum, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyhedraFaceIndex, &
                                      localPolyhedraCellToFaceConnectivity, &
                                      localPolyhedraFaceConnectivityIndex, &
                                      localPolyhedraFaceConnectivity, &
                                      distantPointsCoordinates, &
                                      distantPointsLocation, &
                                      distantPointsBarycentricCoordinatesIndex, &
                                      distantPointsBarycentricCoordinates, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer :: entitiesDim
         integer :: nLocalVertex
         integer :: nLocalElement
         integer :: nLocalPolyhedra
         integer :: nDistantPoint
         integer :: stored
         double precision, dimension(*) :: localCoordinates
         integer, dimension(*) :: localParentEltNum
         integer, dimension(*) :: localConnectivityIndex
         integer, dimension(*) :: localConnectivity
         integer, dimension(*) :: localPolyhedraFaceIndex
         integer, dimension(*) :: localPolyhedraCellToFaceConnectivity
         integer, dimension(*) :: localPolyhedraFaceConnectivityIndex
         integer, dimension(*) :: localPolyhedraFaceConnectivity
         double precision, dimension(*) :: distantPointsCoordinates
         integer, dimension(*) :: distantPointsLocation
         integer, dimension(*) :: distantPointsBarycentricCoordinatesIndex
         double precision, dimension(*) :: distantPointsBarycentricCoordinates
         integer :: stride
         integer :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine ptInterpolationFct
    end interface
    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer :: stride, nStep, status
    integer :: nNotLocatedPoints
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call couplings_send_with_user_interpolation_cf(couplingName, &
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

  end subroutine couplings_send_with_user_interpolation_f_

!
!********************************************************************************
!
! couplings_receive_f
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
!   status               --> Couplings exchange status
!
!********************************************************************************
!

  subroutine couplings_receive_f_ (couplingName, &
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

    call couplings_receive_cf(couplingName, &
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

  end subroutine couplings_receive_f_

!
!********************************************************************************
!
! couplings_exchange_f_with_user_interpolation_f
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
!   status               --> Couplings exchange status
!
!********************************************************************************
!

  subroutine couplings_exchange_with_user_interpolation_f_ (couplingName, &
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
                                      stored, &
                                      localParentEltNum, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyhedraFaceIndex, &
                                      localPolyhedraCellToFaceConnectivity, &
                                      localPolyhedraFaceConnectivityIndex, &
                                      localPolyhedraFaceConnectivity, &
                                      distantPointsCoordinates, &
                                      distantPointsLocation, &
                                      distantPointsBarycentricCoordinatesIndex, &
                                      distantPointsBarycentricCoordinates, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)

         integer :: entitiesDim
         integer :: nLocalVertex
         integer :: nLocalElement
         integer :: nLocalPolyhedra
         integer :: nDistantPoint
         integer :: stored
         double precision, dimension(*) :: localCoordinates
         integer, dimension(*) :: localParentEltNum
         integer, dimension(*) :: localConnectivityIndex
         integer, dimension(*) :: localConnectivity
         integer, dimension(*) :: localPolyhedraFaceIndex
         integer, dimension(*) :: localPolyhedraCellToFaceConnectivity
         integer, dimension(*) :: localPolyhedraFaceConnectivityIndex
         integer, dimension(*) :: localPolyhedraFaceConnectivity
         double precision, dimension(*) :: distantPointsCoordinates
         integer, dimension(*) :: distantPointsLocation
         integer, dimension(*) :: distantPointsBarycentricCoordinatesIndex
         double precision, dimension(*) :: distantPointsBarycentricCoordinates
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

    call couplings_exchange_with_user_interpolation_cf(couplingName, &
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

  end subroutine couplings_exchange_with_user_interpolation_f_

!
!********************************************************************************
!
! couplings_delete_coupling_f
!
! Delete a coupling
!
! parameters
!   couplingName         <-- Coupling identifier
!
!
!********************************************************************************
!

  subroutine couplings_delete_coupling_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call couplings_delete_coupling_cf(couplingName, lCouplingName)

  end subroutine couplings_delete_coupling_f_

!
!********************************************************************************
!
! couplings_finalize()  (define into couplings_cf.hxx)
!
! Finalize couplings. This is a synchronization point between all applications
!
!********************************************************************************
!

!
!********************************************************************************
!
! couplings_get_not_located_points_f
!
! Get located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   notLocatedPoints     --> Not located points (1 to n)
!
!********************************************************************************
!
  subroutine couplings_get_not_located_points_f_(couplingName, notLocatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer, dimension(*) :: notLocatedPoints
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call couplings_get_not_located_points_cf(couplingName, lCouplingName, notLocatedPoints)
  end subroutine couplings_get_not_located_points_f_

!
!********************************************************************************
!
! couplings_get_n_located_distant_points_f
!
! Get located points
!
! parameters
!   couplingName             <-- Coupling identifier
!   nLocatedDistantPoints     --> Number oflocated distan points
!
!********************************************************************************
!
  subroutine couplings_get_n_located_distant_points_f_(couplingName, nLocatedDistantPoints)

    implicit none

    character (len = *) :: couplingName
    integer :: nLocatedDistantPoints
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call couplings_get_n_located_distant_points_cf(couplingName, &
                                                   lCouplingName, &
                                                   nLocatedDistantPoints)
  end subroutine couplings_get_n_located_distant_points_f_



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

  subroutine couplings_get_located_points_f_(couplingName, locatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer, dimension(*) :: locatedPoints
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call couplings_get_located_points_cf(couplingName, lCouplingName, locatedPoints)
  end subroutine couplings_get_located_points_f_

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

  subroutine couplings_get_n_located_points_f_(couplingName, nLocatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer :: nLocatedPoints
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call couplings_get_n_located_points_cf(couplingName, &
                                         lCouplingName, &
                                         nLocatedPoints)

  end subroutine couplings_get_n_located_points_f_

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

  subroutine couplings_get_n_not_located_points_f_(couplingName, nNotLocatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer :: nNotLocatedPoints
    integer :: lCouplingName

    lCouplingName       = len(couplingName)

    call couplings_get_n_not_located_points_cf(couplingName, &
                                               lCouplingName, &
                                               nNotLocatedPoints)

  end subroutine couplings_get_n_not_located_points_f_

end module couplings
