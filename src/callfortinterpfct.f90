subroutine callfortinterpfct(entities_dim, &
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
                             stride, &
                             solver_type, &
                             local_field, &
                             distant_field, &
                             ptInterpolationFct )

  implicit none

  interface
     subroutine  ptInterpolationFct(entities_dim, &
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
                                    stride, &
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
       integer :: stride
       integer :: solver_type
       double precision, dimension(*) :: local_field
       double precision, dimension(*) :: distant_field
     end subroutine ptInterpolationFct
  end interface

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
  integer :: stride
  integer :: solver_type
  double precision, dimension(*) :: local_field
  double precision, dimension(*) :: distant_field


  call ptInterpolationFct (entities_dim, &
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
                           stride, &
                           solver_type, &
                           local_field, &
                           distant_field)

end subroutine callfortinterpfct

