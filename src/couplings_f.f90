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
  ! couplings_dimension_t
  integer, parameter :: couplings_dimension_scalar = 0
  integer, parameter :: couplings_dimension_interlaced_vector = 1
  !
  ! couplings_type_t
  integer, parameter :: couplings_type_float = 0
  integer, parameter :: couplings_type_double = 1
  !
  ! couplings_interpolation_t
  integer, parameter :: couplings_interpolation_standard = 0
  integer, parameter :: couplings_interpolation_user = 1
  !
  ! couplings_not_located_point_treatment_t
  integer, parameter :: couplings_not_located_point_treatment_standard = 0
  integer, parameter :: couplings_not_located_point_treatnent_user = 1
  !
  ! mesh type   
  integer, parameter :: couplings_static_mesh = 0
  integer, parameter :: couplings_mobile_mesh = 1 
  !
  ! solver type
  integer, parameter :: couplings_solver_cell_center = 0
  integer, parameter :: couplings_solver_cell_vertex = 1 

  !
  ! Commons
  ! -------
  !
  ! Logical unit for listing
  integer, save :: ifile

  interface couplings_exchange_f
     module procedure couplings_exchange_without_user_interpolation_f, &
                      couplings_exchange_with_user_interpolation_f
  end interface

  interface couplings_send_f
     module procedure couplings_send_without_user_interpolation_f, &
                      couplings_send_with_user_interpolation_f
  end interface

contains

!
!*******************************************************************************
!
! couplings_init_f 
!
!*******************************************************************************
!
  
  subroutine couplings_init_f (globalcomm, outputunit, appliname, applicomm)
    
    implicit none
    
    integer :: globalcomm, outputUnit, applicomm
    character (len = *) :: appliname

    integer :: l1
    integer :: i
    
    l1 = len(appliname)
    ifile =  outputUnit
    
    call couplings_init_cf (globalcomm, outputunit, appliname, l1,  applicomm)
    
  end subroutine couplings_init_f

!
!********************************************************************************
!
! couplings_add_local_int_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_add_local_int_control_parameter_f (name, initialvalue)

    implicit none

    character (len = *) :: name
    integer :: initialvalue

    integer :: l

    l = len(name)

    call couplings_add_local_int_control_parameter_cf (name, l, initialvalue) 

  end subroutine couplings_add_local_int_control_parameter_f

!
!********************************************************************************
!
! couplings_add_local_double_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_add_local_double_control_parameter_f (name, initialvalue)

    implicit none

    character (len = *) ::name
    double precision :: initialvalue

    integer :: l

    l = len(name)

    call couplings_add_local_double_control_parameter_cf (name, l, initialvalue) 

  end subroutine couplings_add_local_double_control_parameter_f

!
!********************************************************************************
!
! couplings_set_local_int_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_set_local_int_control_parameter_f(name, initialvalue)
        
    implicit none

    character (len = *) ::name
    integer :: initialvalue

    integer :: l

    l = len(name)

    call couplings_set_local_int_control_parameter_cf (name, l, initialvalue) 

  end subroutine couplings_set_local_int_control_parameter_f

!
!********************************************************************************
!
! couplings_set_local_double_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_set_local_double_control_parameter_f (name, initialvalue)

    implicit none

    character (len = *) :: name
    double precision :: initialvalue

    integer :: l

    l = len(name)
    
    call couplings_set_local_double_control_parameter_cf (name, l, initialvalue) 

  end subroutine couplings_set_local_double_control_parameter_f

!
!********************************************************************************
!
! couplings_get_local_int_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_get_local_int_control_parameter_f (name, value)

    implicit none
    
    character (len = *) :: name
    integer ::value

    integer :: l
    
    l = len(name)
    
    call couplings_get_local_int_control_parameter_cf (name, l, value) 
    
  end subroutine couplings_get_local_int_control_parameter_f

!
!********************************************************************************
!
! couplings_get_local_double_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_get_local_double_control_parameter_f (name, value)

    implicit none
    
    character (len = *) :: name
    double precision :: value
    
    integer :: l

    l = len(name)

    call couplings_get_local_double_control_parameter_cf (name, l, value) 

  end subroutine couplings_get_local_double_control_parameter_f

!
!********************************************************************************
!
! couplings_delete_local_int_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_delete_local_int_control_parameter_f (name)

    implicit none

    character (len = *) :: name
    integer l

    l = len(name)

    call couplings_delete_local_int_control_parameter_cf (name, l) 

  end subroutine couplings_delete_local_int_control_parameter_f

!
!********************************************************************************
!
! couplings_delete_local_double_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_delete_local_double_control_parameter_f (name)

    implicit none

    character (len = *) :: name
    integer :: l

    l = len(name)

    call couplings_delete_local_double_control_parameter_cf (name, l) 

  end subroutine couplings_delete_local_double_control_parameter_f

!
!********************************************************************************
!
! couplings_get_distant_int_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_get_distant_int_control_parameter_f (appliName, &
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
    
  end subroutine couplings_get_distant_int_control_parameter_f

!
!********************************************************************************
!
! couplings_get_distant_double_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_get_distant_double_control_parameter_f (appliName, &
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

  end subroutine couplings_get_distant_double_control_parameter_f

!
!********************************************************************************
!
! couplings_synchronise_control_parameter_f
!
!********************************************************************************
!

  subroutine couplings_synchronise_control_parameter_f (appliName)

    implicit none

    character (len = *) :: appliname
      
    integer l
    
    l = len(appliname)

    call couplings_synchronise_control_parameter_cf (appliName, l)

  end subroutine couplings_synchronise_control_parameter_f
 
!
!********************************************************************************
!
! couplings_create_coupling_f
!
!********************************************************************************
!

  subroutine couplings_create_coupling_f (couplingname, &
                                          cplappli, &
                                          dim, &
                                          tolerance, &
                                          mesht, &
                                          solvert, &
                                          outputfreq, &
                                          outputfmt, &
                                          outputfmtopt)

    implicit none

    character (len = *) :: couplingname, cplappli
    integer :: dim, mesht, solvert
    double precision :: tolerance
    integer :: outputfreq
    character (len = *) :: outputfmt, outputfmtopt

    integer :: lcouplingname, lcplappli
    integer :: loutputfmt, loutputfmtopt

    lcouplingname = len(couplingname)
    lcplappli     = len(cplappli)
    loutputfmt    = len(outputfmt)
    loutputfmtopt = len(outputfmtopt)

    call couplings_create_coupling_cf(couplingname, &
                                      lcouplingname, &
                                      cplappli, &
                                      lcplappli, &
                                      dim, &
                                      tolerance, &
                                      mesht, &
                                      solvert, &
                                      outputfreq, &
                                      outputfmt, &
                                      loutputfmt, &
                                      outputfmtopt, &
                                      loutputfmtopt)

  end subroutine couplings_create_coupling_f

!
!********************************************************************************
!
! couplings_set_points_to_locate_f
!
!********************************************************************************
!

  subroutine couplings_set_points_to_locate_f (couplingname, &
                                               npts, &
                                               coords)

    implicit none
    
    character (len = *) :: couplingname
      
    integer :: npts
    double precision, dimension(3 * npts) :: coords
    
    integer :: lcouplingname

    lcouplingname  = len(couplingname)

    call couplings_set_points_to_locate_cf(couplingname, & 
                                           lcouplingname, & 
                                           npts, &
                                           coords)
    
  end subroutine couplings_set_points_to_locate_f

!
!********************************************************************************
!
! couplings_define_mesh_f
!
!********************************************************************************
!

  subroutine couplings_define_mesh_f (couplingName, &
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

  end subroutine couplings_define_mesh_f

!
!********************************************************************************
!
! couplings_add_polyhedra_f
!
!*******************************************************************************
!

  subroutine couplings_add_polyhedra_f (couplingname, &
                                        nelts, &
                                        faceidx, &
                                        celltoface, &
                                        faceconnecidx, &
                                        faceconnec)

    implicit none

    character (len = *) :: couplingname
    integer :: lcouplingname, nelts
    integer, dimension(nelts) :: faceidx
    integer, dimension(*) :: celltoface, faceconnecidx, faceconnec

    lcouplingname = len(couplingname)

    call couplings_add_polyhedra_cf (couplingname, &
                                     lcouplingname, &
                                     nelts, &
                                     faceidx, &
                                     celltoface, &
                                     faceconnecidx, &
                                     faceconnec)

  end subroutine couplings_add_polyhedra_f

!
!********************************************************************************
!
! couplings_exchange_f
!
!********************************************************************************
!

  subroutine couplings_exchange_without_user_interpolation_f (couplingname, &
                                   exchangename, &
                                   exchangedim, &
                                   nstep, &
                                   timevalue, &
                                   sendingfieldname, &
                                   sendingfield, &
                                   receivingfieldname, &
                                   receivingfield, &
                                   status)
    
    implicit none

    character (len = *) :: couplingname, exchangename, sendingfieldname
    character (len = *) :: receivingfieldname
    integer :: exchangedim, nstep, status
    double precision :: timevalue
    double precision, dimension(*) :: sendingfield, receivingfield
    
    integer :: lcouplingname, lexchangename, lsendingfieldname
    integer :: lreceivingfieldname

    lcouplingname       = len(couplingname)
    lexchangename       = len(exchangename)
    lsendingfieldname   = len(sendingfieldname)
    lreceivingfieldname = len(receivingfieldname)
    
    call couplings_exchange_cf(couplingname, &
                               lcouplingname, &
                               exchangename, &
                               lexchangename, &
                               exchangedim, &
                               nstep, &
                               timevalue, &
                               sendingfieldname, &
                               lsendingfieldname, &
                               sendingfield, &
                               receivingfieldname, &
                               lreceivingfieldname, &
                               receivingfield, &
                               status)
    
  end subroutine couplings_exchange_without_user_interpolation_f

!
!********************************************************************************
!
! couplings_send_f
!
!********************************************************************************
!

  subroutine couplings_send_without_user_interpolation_f (couplingname, &
                               exchangename, &
                               exchangedim, &
                               nstep, &
                               timevalue, &
                               sendingfieldname, &
                               sendingfield, &
                               status)
    
    implicit none

    character (len = *) :: couplingname, exchangename, sendingfieldname
    integer :: exchangedim, nstep, status
    double precision :: timevalue
    double precision, dimension(*) :: sendingfield 
    
    integer :: lcouplingname, lexchangename, lsendingfieldname
 
    lcouplingname       = len(couplingname)
    lexchangename       = len(exchangename)
    lsendingfieldname   = len(sendingfieldname)
     
    call couplings_send_cf(couplingname, &
                           lcouplingname, &
                           exchangename, &
                           lexchangename, &
                           exchangedim, &
                           nstep, &
                           timevalue, &
                           sendingfieldname, &
                           lsendingfieldname, &
                           sendingfield, &
                           status)
    
  end subroutine couplings_send_without_user_interpolation_f


!
!********************************************************************************
!
! couplings_send_with_user_interpolation_f
!
!********************************************************************************
!

  subroutine couplings_send_with_user_interpolation_f (couplingname, &
                               exchangename, &
                               exchangedim, &
                               nstep, &
                               timevalue, &
                               sendingfieldname, &
                               sendingfield, &
                               ptInterpolationFct, &
                               status)
    
    implicit none

    interface
       subroutine  ptInterpolationFct(entities_dim, &
                                      n_local_vertex, &
                                      n_local_element, &
                                      n_local_polhyedra, &
                                      n_distant_point, &
                                      local_coordinates, &
                                      stored, &
                                      local_parent_vertex_num, &
                                      local_connectivity_index, &
                                      local_connectivity, &
                                      local_polyhedra_face_index, &
                                      local_polyhedra_cell_to_face_connectivity, &
                                      local_polyhedra_face_connectivity_index, &
                                      local_polyhedra_face_connectivity, &
                                      distant_points_coordinates, &
                                      distant_points_location, &
                                      distant_points_barycentric_coordinates_index, &
                                      distant_points_barycentric_coordinates, &
                                      data_dimension, &
                                      solver_type, &
                                      local_field, &
                                      distant_field)
         integer :: entities_dim
         integer :: n_local_vertex
         integer :: n_local_element
         integer :: n_local_polhyedra
         integer :: n_distant_point
         integer :: stored
         double precision, dimension(*) :: local_coordinates
         integer, dimension(*) :: local_parent_vertex_num
         integer, dimension(*) :: local_connectivity_index
         integer, dimension(*) :: local_connectivity
         integer, dimension(*) :: local_polyhedra_face_index
         integer, dimension(*) :: local_polyhedra_cell_to_face_connectivity
         integer, dimension(*) :: local_polyhedra_face_connectivity_index
         integer, dimension(*) :: local_polyhedra_face_connectivity
         double precision, dimension(*) :: distant_points_coordinates
         integer, dimension(*) :: distant_points_location
         integer, dimension(*) :: distant_points_barycentric_coordinates_index
         double precision, dimension(*) :: distant_points_barycentric_coordinates
         integer :: data_dimension
         integer :: solver_type
         double precision, dimension(*) :: local_field
         double precision, dimension(*) :: distant_field
       end subroutine ptInterpolationFct
    end interface
    character (len = *) :: couplingname, exchangename, sendingfieldname
    integer :: exchangedim, nstep, status
    double precision :: timevalue
    double precision, dimension(*) :: sendingfield
    
    integer :: lcouplingname, lexchangename, lsendingfieldname
 
    lcouplingname       = len(couplingname)
    lexchangename       = len(exchangename)
    lsendingfieldname   = len(sendingfieldname)
     
    call couplings_send_with_user_interpolation_cf(couplingname, &
                                                   lcouplingname, &
                                                   exchangename, &
                                                   lexchangename, &
                                                   exchangedim, &
                                                   nstep, &
                                                   timevalue, &
                                                   sendingfieldname, &
                                                   lsendingfieldname, &
                                                   sendingfield, &
                                                   ptInterpolationFct, &
                                                   status)
    
  end subroutine couplings_send_with_user_interpolation_f

!
!********************************************************************************
!
! couplings_receive_f
!
!********************************************************************************
!

  subroutine couplings_receive_f (couplingname, &
                                  exchangename, &
                                  exchangedim, &
                                  nstep, &
                                  timevalue, &
                                  receivingfieldname, &
                                  receivingfield, &
                                  status)
    
    implicit none

    character (len = *) :: couplingname, exchangename, receivingfieldname
    integer :: exchangedim, nstep, status
    double precision :: timevalue
    double precision, dimension(*) :: receivingfield
    
    integer :: lcouplingname, lexchangename, lreceivingfieldname
 
    lcouplingname       = len(couplingname)
    lexchangename       = len(exchangename)
    lreceivingfieldname   = len(receivingfieldname)
     
    call couplings_receive_cf(couplingname, &
                              lcouplingname, &
                              exchangename, &
                              lexchangename, &
                              exchangedim, &
                              nstep, &
                              timevalue, &
                              receivingfieldname, &
                              lreceivingfieldname, &
                              receivingfield, &
                              status)
    
  end subroutine couplings_receive_f

!
!********************************************************************************
!
! couplings_exchange_f_with_user_interpolation_f
!
!********************************************************************************
!

  subroutine couplings_exchange_with_user_interpolation_f (couplingName, &
                                                           exchangeName, &
                                                           exchangeDim, &
                                                           nStep, &
                                                           timeValue, &
                                                           sendingFieldName, &
                                                           sendingField, &
                                                           receivingFieldName, &
                                                           receivingField, &
                                                           ptInterpolationFct, &
                                                           status)

    implicit none
    
    interface
       subroutine  ptInterpolationFct(entities_dim, &
                                      n_local_vertex, &
                                      n_local_element, &
                                      n_local_polhyedra, &
                                      n_distant_point, &
                                      local_coordinates, &
                                      stored, &
                                      local_parent_vertex_num, &
                                      local_connectivity_index, &
                                      local_connectivity, &
                                      local_polyhedra_face_index, &
                                      local_polyhedra_cell_to_face_connectivity, &
                                      local_polyhedra_face_connectivity_index, &
                                      local_polyhedra_face_connectivity, &
                                      distant_points_coordinates, &
                                      distant_points_location, &
                                      distant_points_barycentric_coordinates_index, &
                                      distant_points_barycentric_coordinates, &
                                      data_dimension, &
                                      solver_type, &
                                      local_field, &
                                      distant_field)
         integer :: entities_dim
         integer :: n_local_vertex
         integer :: n_local_element
         integer :: n_local_polhyedra
         integer :: n_distant_point
         integer :: stored
         double precision, dimension(*) :: local_coordinates
         integer, dimension(*) :: local_parent_vertex_num
         integer, dimension(*) :: local_connectivity_index
         integer, dimension(*) :: local_connectivity
         integer, dimension(*) :: local_polyhedra_face_index
         integer, dimension(*) :: local_polyhedra_cell_to_face_connectivity
         integer, dimension(*) :: local_polyhedra_face_connectivity_index
         integer, dimension(*) :: local_polyhedra_face_connectivity
         double precision, dimension(*) :: distant_points_coordinates
         integer, dimension(*) :: distant_points_location
         integer, dimension(*) :: distant_points_barycentric_coordinates_index
         double precision, dimension(*) :: distant_points_barycentric_coordinates
         integer :: data_dimension
         integer :: solver_type
         double precision, dimension(*) :: local_field
         double precision, dimension(*) :: distant_field
       end subroutine ptInterpolationFct
    end interface

    character (len = *) :: couplingName, exchangeName, sendingFieldName
    character (len = *) :: receivingFieldName
    integer :: exchangeDim, nStep, status
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
                                                       status)

  end subroutine couplings_exchange_with_user_interpolation_f

!
!********************************************************************************
!
! couplings_delete_coupling_f
!
!********************************************************************************
!

  subroutine couplings_delete_coupling_f(couplingName)

    implicit none

    character (len = *) :: couplingName 
    integer :: lCouplingName

    lCouplingName       = len(couplingName)
    
    call couplings_delete_coupling_cf(couplingName, lCouplingName)
    
  end subroutine couplings_delete_coupling_f

end module couplings
