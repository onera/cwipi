subroutine callfortinterpfct(entities_dim, &
                             n_local_vertex, &
                             n_local_element, &
                             n_local_polhyedra, &
                             n_distant_point, &
                             local_coordinates, &
                             local_connectivity_index, &
                             local_connectivity, &
                             local_poly_face_index, &
                             local_poly_cell2face_connec, &
                             local_poly_face_connec_idx, &
                             local_poly_face_connec, &
                             dist_pts_coord, &
                             dist_pts_location, &
                             dist_pts_barycentric_coord_idx, &
                             dist_pts_barycentric_coord, &
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
                                    local_poly_face_index, &
                                    local_poly_cell2face_connec, &
                                    local_poly_face_connec_idx, &
                                    local_poly_face_connec, &
                                    dist_pts_coord, &
                                    dist_pts_location, &
                                    dist_pts_barycentric_coord_idx, &
                                    dist_pts_barycentric_coord, &
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
       integer, dimension(*) :: local_poly_face_index
       integer, dimension(*) :: local_poly_cell2face_connec
       integer, dimension(*) :: local_poly_face_connec_idx
       integer, dimension(*) :: local_poly_face_connec
       double precision, dimension(*) :: dist_pts_coord
       integer, dimension(*) :: dist_pts_location
       integer, dimension(*) :: dist_pts_barycentric_coord_idx
       double precision, dimension(*) :: dist_pts_barycentric_coord
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
  integer, dimension(*) :: local_poly_face_index
  integer, dimension(*) :: local_poly_cell2face_connec
  integer, dimension(*) :: local_poly_face_connec_idx
  integer, dimension(*) :: local_poly_face_connec
  double precision, dimension(*) :: dist_pts_coord
  integer, dimension(*) :: dist_pts_location
  integer, dimension(*) :: dist_pts_barycentric_coord_idx
  double precision, dimension(*) :: dist_pts_barycentric_coord
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
                           local_poly_face_index, &
                           local_poly_cell2face_connec, &
                           local_poly_face_connec_idx, &
                           local_poly_face_connec, &
                           dist_pts_coord, &
                           dist_pts_location, &
                           dist_pts_barycentric_coord_idx, &
                           dist_pts_barycentric_coord, &
                           stride, &
                           solver_type, &
                           local_field, &
                           distant_field)

end subroutine callfortinterpfct

