
subroutine printStatus(iiunit, status)
  use couplings
  implicit none
  integer :: status
  integer :: iiunit

  select case (status)
     case(couplings_exchange_ok)
        write(iiunit,*) "Exchange ok"
     case(couplings_exchange_bad_receiving)
        write(iiunit,*) "no or bad receiving"
     case default
        write(iiunit,*) "Unknown receiving status"
        stop
  end select

end subroutine printStatus

!
! ------------------------------------------------------------------------------
! Fonction d'interpolation bidon, juste pour voir si c'est bien pris en compte
! ------------------------------------------------------------------------------
!

subroutine  interpolationbidon_f(entities_dim, &
                                 n_local_vertex, &
                                 n_local_element, &
                                 n_local_polhyedra, &
                                 n_distant_point, &
                                 local_coordinates, &
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
  double precision, dimension(*) :: local_coordinates
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

  integer :: i

  do i = 1, n_distant_point
     distant_field(i) = i
  enddo

end subroutine interpolationbidon_f

!
! -----------------------------------------------------------------------------
! Programme de tests
! -----------------------------------------------------------------------------
!

program testf
 
  use mpi
  use couplings
 
  implicit none

  interface
     subroutine interpolationbidon_f(entities_dim, &
                                     n_local_vertex, &
                                     n_local_element, &
                                     n_local_polhyedra, &
                                     n_distant_point, &
                                     local_coordinates, &
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
       double precision, dimension(*) :: local_coordinates
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
     end subroutine interpolationbidon_f

  end interface

  integer, allocatable, dimension(:) :: location
  integer, allocatable, dimension(:) :: baryCooIdx
  double precision, allocatable, dimension(:) :: baryCoo
  integer :: nLocatedPoints
  integer :: nNotLocatedPoints
  integer :: nDistantPoints

  integer :: localcom, localGroup, p1Group, p1Comm
  integer :: irank, currentRank
  character (len = 4) :: proc
  integer :: code
  integer :: iiunit
  integer :: ivalue
  double precision :: dvalue
      
  double precision, parameter :: xmin = -100.d0
  double precision, parameter :: xmax =  100.d0
  double precision, parameter :: ymin = -100.d0
  double precision, parameter :: ymax =  100.d0
  integer, parameter  :: nx   = 24    
  integer, parameter  :: ny   = 24
  integer, parameter  :: initrandom = 2

  integer nvertex, nelts, lconnecindex

  integer, parameter :: nvertexm = 4000 
  integer, parameter :: neltsm = 4000
  integer, parameter :: lconnecindexm = 12000
  integer, parameter :: nptstolocate = 21
  double precision, dimension(3*nptstolocate) :: coordstolocate
    
  double precision, dimension(3*nvertexm) :: coords
  integer, dimension(neltsm+1) :: connecindex
  integer, dimension(lconnecindexm) :: connec
  
  double precision, dimension(3*nvertexm) :: values, localvalues

  integer status
  
  integer i, order
  
  integer :: vpar = 10
  character (len = 6) :: cpar = "niterf"

  integer :: stride = 1
  integer, dimension(1) :: rl

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, irank, code)


!
! -----------------------------------------
! Initialisation de l'interface de couplage
! -----------------------------------------
!
  
  call couplings_init_f (mpi_comm_world, & 
                         "CodeFortran", & 
                         localcom)


!
! ------------------------------------------------
! Creation du fichier de sortie listing
! (en parallele, un fichier listing par processus)
! ------------------------------------------------ 
!

  call mpi_comm_rank(localcom, currentRank, code)


  write(proc,'(i4.4)') currentRank
  iiunit = 9
  open(unit=iiunit, file='listing_test2D_1_c2_'//proc, &
       form='formatted', status='unknown')

  call couplings_set_output_listing_f(iiunit)

  write(iiunit,*)
  write(iiunit,*) "dump apres initialisation"
  write(iiunit,*) "-------------------------"
  write(iiunit,*)

  call couplings_dump_application_properties_f

!
! -------------------------------
! Test des parametres de controle
! -------------------------------
!

!
! Ajout de parametres de controle

!  call couplings_add_local_int_control_parameter_f("niterf", 10)
  call couplings_add_local_int_control_parameter_f(cpar, vpar)

  call couplings_add_local_double_control_parameter_f("physicaltimef", 1.123d0)

  write(iiunit,*)
  write(iiunit,*) "dump apres ajout de parametres"
  write(iiunit,*) "------------------------------"
  write(iiunit,*)
  call couplings_dump_application_properties_f
!
! Modification des parametres de controle
  call couplings_get_local_int_control_parameter_f("niterf", ivalue)

  ivalue = ivalue + 1

  call couplings_set_local_int_control_parameter_f("niterf", ivalue)

  call couplings_get_local_double_control_parameter_f("physicaltimef", dvalue)

  dvalue = dvalue + 0.1d0

  call couplings_set_local_double_control_parameter_f("physicaltimef", dvalue)

  write(iiunit,*)
  write(iiunit,*) "dump apres modification de parametres"
  write(iiunit,*) "-------------------------------------"
  write(iiunit,*)
  call couplings_dump_application_properties_f
!
! Echange des parametres de controle

  call couplings_synchronize_control_parameter_f("CodeC")

  write(iiunit,*)
  write(iiunit,*) "dump apres synchronisation"
  write(iiunit,*) "--------------------------"
  write(iiunit,*)
  call couplings_dump_application_properties_f
!
! Suppression des parametres de controle
  call couplings_delete_local_int_control_parameter_f("niterf")
  call couplings_delete_local_double_control_parameter_f("physicaltimef")

  write(iiunit,*)
  write(iiunit,*) "dump apres suppression des parametres"
  write(iiunit,*) "-------------------------------------"
  write(iiunit,*)
  call couplings_dump_application_properties_f
! 
! ------------------------
! Test couplage p1 <-> p1
! ------------------------
! 
! Construction de "l'objet" couplage

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 1"
  write(iiunit, *)

  call couplings_create_coupling_f("test2D_1", & 
                                   couplings_cpl_parallel_with_part,&
                                   "CodeC", &          
                                   2,     & ! Dimension des entites geometriques                 
                                   0.1d0, & ! Tolerance geometrique                        
                                   couplings_static_mesh, &       
                                   couplings_solver_cell_vertex, &
                                   1, &                          
                                   "Ensight Gold",&              
                                   "text")                     

! 
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm  
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call couplings_define_mesh_f("test2D_1", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

! 
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo
  
  stride = 1
  call couplings_exchange_f ("test2D_1", &
                             "echange1", &
                             stride, & 
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call couplings_delete_coupling_f("test2D_1")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

! 
! -------------------------------------
! Test couplage P1 -> P0 puis P0 -> P1 
! -------------------------------------
!
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 2"
  write(iiunit, *)

  call couplings_create_coupling_f("test2D_2", & 
                                   couplings_cpl_parallel_with_part,&
                                   "CodeC", &          
                                   2,     & ! Dimension des entites geometriques                 
                                   0.1d0, & ! Tolerance geometrique                        
                                   couplings_static_mesh, &       
                                   couplings_solver_cell_vertex, &
                                   1, &                          
                                   "Ensight Gold",&              
                                   "text")                     

! 
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm  
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call couplings_define_mesh_f("test2D_2", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

! 
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo
  
  stride = 1
  call couplings_send_f ("test2D_2", &
                         "echange1", &
                         stride, & 
                         1, &
                         0.1d0, &
                         "cooY", &
                         values, &
                         status)
  call printStatus(iiunit, status)

  call couplings_receive_f ("test2D_2", &
                            "echange2", &
                            stride, & 
                            1, &
                            0.1d0, &
                            "cooYY", &
                            values, &
                            nNotLocatedPoints, &
                            status)
  call printStatus(iiunit, status)


!
! Suppression de l'objet couplage "couplingcellvertex"

  call couplings_delete_coupling_f("test2D_2");
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

! 
! -------------------------------------
! Test de definition des points 
! a interpoler
! -------------------------------------
!
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 3"
  write(iiunit, *)

  call couplings_create_coupling_f("test2D_3", & 
                                   couplings_cpl_parallel_with_part,&
                                   "CodeC", &          
                                   2,     & ! Dimension des entites geometriques                 
                                   0.1d0, & ! Tolerance geometrique                        
                                   couplings_static_mesh, &       
                                   couplings_solver_cell_vertex, &
                                   1, &                          
                                   "Ensight Gold",&              
                                   "text")                     

! 
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm  
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call couplings_define_mesh_f("test2D_3", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

! 
! Definition des points a localiser

  coordstolocate(1) = -75.d0
  coordstolocate(2) = -75.d0
  coordstolocate(3) = 0.d0

  coordstolocate(4) = -75.d0
  coordstolocate(5) = -50.d0
  coordstolocate(6) = 0.d0

  coordstolocate(7) = -75.d0
  coordstolocate(8) = -25.d0
  coordstolocate(9) = 0.d0

  coordstolocate(10) = -75.d0
  coordstolocate(11) = 0.d0
  coordstolocate(12) = 0.d0

  coordstolocate(13) = -75.d0
  coordstolocate(14) = 25.d0
  coordstolocate(15) = 0.d0

  coordstolocate(16) = -75.d0
  coordstolocate(17) = 50.d0
  coordstolocate(18) = 0.d0

  coordstolocate(19) = -75.d0
  coordstolocate(20) = 75.d0
  coordstolocate(21) = 0.d0

  coordstolocate(22) = -25.d0
  coordstolocate(23) = -75.d0
  coordstolocate(24) = 0.d0

  coordstolocate(25) = -25.d0
  coordstolocate(26) = -50.d0
  coordstolocate(27) = 0.d0

  coordstolocate(28) = -25.d0
  coordstolocate(29) = -25.d0
  coordstolocate(30) = 0.d0

  coordstolocate(31) = -25.d0
  coordstolocate(32) = 0.d0
  coordstolocate(33) = 0.d0

  coordstolocate(34) = -25.d0
  coordstolocate(35) = 25.d0
  coordstolocate(36) = 0.d0

  coordstolocate(37) = -25.d0
  coordstolocate(38) = 50.d0
  coordstolocate(39) = 0.d0

  coordstolocate(40) = -25.d0
  coordstolocate(41) = 75.d0
  coordstolocate(42) = 0.d0

  coordstolocate(43) = 25.d0
  coordstolocate(44) = -75.d0
  coordstolocate(45) = 0.d0

  coordstolocate(46) = 25.d0
  coordstolocate(47) = -50.d0
  coordstolocate(48) = 0.d0

  coordstolocate(49) = 25.d0
  coordstolocate(50) = -25.d0
  coordstolocate(51) = 0.d0

  coordstolocate(52) = 25.d0
  coordstolocate(53) = 0.d0
  coordstolocate(54) = 0.d0

  coordstolocate(55) = 25.d0
  coordstolocate(56) = 25.d0
  coordstolocate(57) = 0.d0

  coordstolocate(58) = 25.d0
  coordstolocate(59) = 50.d0
  coordstolocate(60) = 0.d0

  coordstolocate(61) = 25.d0
  coordstolocate(62) = 75.d0
  coordstolocate(63) = 0.d0

  call couplings_set_points_to_locate_f ("test2D_3", &
                                         nptstolocate, &
                                         coordstolocate)

! 
! Envoi de la coordonnee Y a codeC
! Reception de la coordonnee Y provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo
  
  stride = 1
  call couplings_exchange_f ("test2D_3", &
                             "echange1", &
                             stride, & 
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)
  write(iiunit,*) "valeurs recues test2D_3"
  write(iiunit,*) (localvalues(i),i=1,nptstolocate)

!
! Suppression de l'objet couplage 

  call couplings_delete_coupling_f("test2D_3");
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

! 
! -------------------------------------
! Test de definition d'une fonction 
! d'interpolation en fortran (callback)
! -------------------------------------
!
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 4"
  write(iiunit, *)

  call couplings_create_coupling_f("test2D_4", & 
                                   couplings_cpl_parallel_with_part,&
                                   "CodeC", &          
                                   2,     & ! Dimension des entites geometriques                 
                                   0.1d0, & ! Tolerance geometrique                        
                                   couplings_static_mesh, &       
                                   couplings_solver_cell_vertex, &
                                   1, &                          
                                   "Ensight Gold",&              
                                   "text")                     

  


! 
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm  
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call couplings_define_mesh_f("test2D_4", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

! 
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo
  
  stride = 1
  call couplings_exchange_f ("test2D_4", &
                             "echange1", &
                             stride, & 
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             interpolationbidon_f, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call couplings_delete_coupling_f("test2D_4")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

! 
! -------------------------------------
! test de la transmission d'un vecteur
! -------------------------------------
!
! 
! Construction de "l'objet" couplage
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 5"
  write(iiunit, *)

  call couplings_create_coupling_f("test2D_5", & 
                                   couplings_cpl_parallel_with_part,&
                                   "CodeC", &          
                                   2,     & ! Dimension des entites geometriques                 
                                   0.1d0, & ! Tolerance geometrique                        
                                   couplings_static_mesh, &       
                                   couplings_solver_cell_vertex, &
                                   1, &                          
                                   "Ensight Gold",&              
                                   "text")                     

! 
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm  
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call couplings_define_mesh_f("test2D_5", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

! 
! Envoi de la coory a codec
! Reception de la coorx provenant de codec
  
  stride = 3
  call couplings_exchange_f ("test2D_5", &
                             "echange1", &
                             stride, & 
                             1, &
                             0.1d0, &
                             "cooy", &
                             coords, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call couplings_delete_coupling_f("test2D_5")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

!
! -------------------------------------
! test des sorties d'erreur
! ------------------------------------- 
!

! 
! Construction de "l'objet" couplage
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 6"
  write(iiunit, *)

  call couplings_create_coupling_f("test2D_6", & 
                                   couplings_cpl_parallel_with_part,&
                                   "CodeC", &          
                                   2,     & ! Dimension des entites geometriques                 
                                   0.1d0, & ! Tolerance geometrique                        
                                   couplings_static_mesh, &       
                                   couplings_solver_cell_vertex, &
                                   1, &                          
                                   "Ensight Gold",&              
                                   "text")                     

! 
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm  
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call couplings_define_mesh_f("test2D_6", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

! 
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo
  
  stride = 1
  call couplings_exchange_f ("test2D_6", &
                             "echange1", &
                             stride, & 
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)

  !
  ! Reception mais pas d'envoi alors que code_C attend quelque chose
  !

  stride = 1
  call couplings_receive_f ("test2D_6", &
                            "echange2", &
                            stride, & 
                            1, &
                            0.1d0, &
                            "cooYY", &
                            values, &
                            nNotLocatedPoints, &
                            status)
  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call couplings_delete_coupling_f("test2D_6")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

! 
! -------------------------------------
! Test simple localisation
! -------------------------------------
!

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 7"
  write(iiunit, *)

  call couplings_create_coupling_f("test2D_7", & 
                                   couplings_cpl_parallel_with_part,&
                                   "CodeC", &          
                                   2,     & ! Dimension des entites geometriques                 
                                   0.1d0, & ! Tolerance geometrique                        
                                   couplings_static_mesh, &       
                                   couplings_solver_cell_vertex, &
                                   1, &                          
                                   "Ensight Gold",&              
                                   "text")                     

! 
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm  
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call couplings_define_mesh_f("test2D_7", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

  call couplings_locate_f("test2D_7")

  call couplings_get_n_located_distant_points_f("test2D_7", nDistantPoints)

  allocate(location(nDistantPoints))
  allocate(baryCooIdx(nDistantPoints+1))

  call couplings_get_location_f("test2D_7", location)

  write(iiunit,*) "location",(location(i),i=1,nDistantPoints)

  call couplings_get_barycentric_coordinates_index_f("test2D_7", baryCooIdx)

  allocate(baryCoo(baryCooIdx(nDistantPoints+1)))

  call couplings_get_barycentric_coordinates_f("test2D_7", baryCoo)

  deallocate(location)
  deallocate(baryCooIdx)
  deallocate(baryCoo)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call couplings_delete_coupling_f("test2D_7")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

! 
! ------------------------
! Test couplage p1 <-> p1
! ------------------------
! 
! Construction de "l'objet" couplage

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 8"
  write(iiunit, *)

  call couplings_create_coupling_f("test2D_8", & 
                                   couplings_cpl_parallel_with_part,&
                                   "CodeC", &          
                                   2,     & ! Dimension des entites geometriques                 
                                   0.1d0, & ! Tolerance geometrique                        
                                   couplings_static_mesh, &       
                                   couplings_solver_cell_vertex, &
                                   1, &                          
                                   "Ensight Gold",&              
                                   "text")                     


! 
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm  
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call couplings_define_mesh_f("test2D_8", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)



! 
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo
  
  stride = 1
  call couplings_exchange_f ("test2D_8", &
                             "echange1", &
                             stride, & 
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)


!!$  call couplings_create_coupling_f("test2D_9", & 
!!$                                   couplings_cpl_parallel_without_part,&
!!$                                   "CodeC", &          
!!$                                   2,     & ! Dimension des entites geometriques                 
!!$                                   0.1d0, & ! Tolerance geometrique                        
!!$                                   couplings_static_mesh, &       
!!$                                   couplings_solver_cell_vertex, &
!!$                                   1, &                          
!!$                                   "Ensight Gold",&              
!!$                                   "text")                     

  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

 ! call couplings_delete_coupling_f("test2D_8")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

! 
! ------------------------
! Test couplage p1 <-> p1
! ------------------------
! 
! Construction de "l'objet" couplage

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 9"
  write(iiunit, *)

  call couplings_create_coupling_f("test2D_9", & 
                                   couplings_cpl_parallel_without_part,&
                                   "CodeC", &          
                                   2,     & ! Dimension des entites geometriques                 
                                   0.1d0, & ! Tolerance geometrique                        
                                   couplings_static_mesh, &       
                                   couplings_solver_cell_vertex, &
                                   1, &                          
                                   "Ensight Gold",&              
                                   "text")                     

! 
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm  
  lconnecindex = lconnecindexm
  order = 1

  call mpi_comm_group(localcom, localGroup, code)
    
  rl(1) = 0
  call mpi_group_incl(localGroup, 1, rl, p1Group, code)
  call mpi_comm_create(localcom, p1Group, p1Comm, code)

  if (currentRank == 0) then

     call creemaillagepolygone2d_f (order, &
          p1Comm, &
          xmin, &
          xmax, &
          ymin, &
          ymax, &
          initrandom, &
          nx, &
          ny, &
          nvertex, &
          coords, &
          nelts, &
          lconnecindex, &
          connecindex, &
          connec)

     call couplings_define_mesh_f("test2D_9", &
                                  nvertex, &
                                  nelts, &
                                  coords, &
                                  connecindex, &
                                  connec)
  endif

! 
! Envoi de la coory a codec
! Reception de la coorx provenant de codec
  stride = 1

  if (currentRank == 0) then
     do i = 1, nvertex
        values(i) = coords(3*(i-1) + 2)
     enddo
  
     call couplings_exchange_f ("test2D_9", &
                                "echange1", &
                                stride, & 
                                1, &
                                0.1d0, &
                                "cooy", &
                                values, &
                                "coox", &
                                localvalues, &
                                nNotLocatedPoints, &
                                status)
  else

     call couplings_receive_f ("test2D_9", &
                               "echange1", &
                               stride, & 
                               1, &
                               0.1d0, &
                               "coox", &
                               values, &
                               nNotLocatedPoints, &
                               status)

  endif

  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call couplings_delete_coupling_f("test2D_9")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

! 
! ------------------------
! Test couplage p1 <-> p1
! ------------------------
! 
! Construction de "l'objet" couplage

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 10"
  write(iiunit, *)

  call couplings_create_coupling_f("test2D_10", & 
                                   couplings_cpl_parallel_without_part,&
                                   "CodeC", &          
                                   2,     & ! Dimension des entites geometriques                 
                                   0.1d0, & ! Tolerance geometrique                        
                                   couplings_static_mesh, &       
                                   couplings_solver_cell_vertex, &
                                   1, &                          
                                   "Ensight Gold",&              
                                   "text")                     

! 
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm  
  lconnecindex = lconnecindexm
  order = 1

  call mpi_comm_group(localcom, localGroup, code)
    
  rl(1) = 0
  call mpi_group_incl(localGroup, 1, rl, p1Group, code)
  call mpi_comm_create(localcom, p1Group, p1Comm, code)

  if (currentRank == 0) then

     call creemaillagepolygone2d_f (order, &
          p1Comm, &
          xmin, &
          xmax, &
          ymin, &
          ymax, &
          initrandom, &
          nx, &
          ny, &
          nvertex, &
          coords, &
          nelts, &
          lconnecindex, &
          connecindex, &
          connec)

     call couplings_define_mesh_f("test2D_10", &
                                  nvertex, &
                                  nelts, &
                                  coords, &
                                  connecindex, &
                                  connec)
  endif

! 
! Envoi de la coory a codec
! Reception de la coorx provenant de codec
  stride = 1

  if (currentRank == 0) then
     do i = 1, nvertex
        values(i) = coords(3*(i-1) + 2)
     enddo
  
     call couplings_exchange_f ("test2D_10", &
                                "echange1", &
                                stride, & 
                                1, &
                                0.1d0, &
                                "cooy", &
                                values, &
                                "coox", &
                                localvalues, &
                                nNotLocatedPoints, &
                                status)
  else

     call couplings_receive_f ("test2D_10", &
                               "echange1", &
                               stride, & 
                               1, &
                               0.1d0, &
                               "coox", &
                               values, &
                               nNotLocatedPoints, &
                               status)

  endif

  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call couplings_delete_coupling_f("test2D_10")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

  call couplings_finalize_f()
      
end program testf
