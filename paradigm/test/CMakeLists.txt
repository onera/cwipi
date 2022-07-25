# --------------------------
# Construction des Tests
# --------------------------

if (MPI_C_COMPILER)
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

if (MPI_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

if (MPI_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

unset(LINK_LIBRARIES)

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)
  if (PDM_ENABLE_STATIC)
    if (PDM_ENABLE_EXTENSION_PDMA)
      list(APPEND LINK_LIBRARIES
         pdmaf_static)
    endif()
    list(APPEND LINK_LIBRARIES pdmf_static)
    if (PDM_ENABLE_EXTENSION_PDMA)
      list(APPEND LINK_LIBRARIES
         pdmaf_static)
    endif()
  else ()
    if (PDM_ENABLE_EXTENSION_PDMA)
      list(APPEND LINK_LIBRARIES
         pdmaf_shared)
    endif()
    list(APPEND LINK_LIBRARIES pdmf_shared)
    if (PDM_ENABLE_EXTENSION_PDMA)
      list(APPEND LINK_LIBRARIES
         pdmaf_shared)
    endif()
  endif()
endif()

if (PDM_ENABLE_STATIC)
 if (PDM_ENABLE_EXTENSION_PDMA)
  list(APPEND LINK_LIBRARIES
       pdma_static)
 endif()
 list(APPEND LINK_LIBRARIES
      pdm_static
      pdm_mpi_static
      ${COMMON_LIBRARIES})
 if (PDM_ENABLE_EXTENSION_PDMA)
  list(APPEND LINK_LIBRARIES
       pdma_static)
 endif()
else ()
 list(APPEND LINK_LIBRARIES
      pdm_shared
      pdm_mpi_shared ${COMMON_LIBRARIES})
 if (PDM_ENABLE_EXTENSION_PDMA)
  list(APPEND LINK_LIBRARIES
       pdma_shared)
 endif()
endif()

if (OPENMP_FOUND)
  list(APPEND LINK_LIBRARIES ${OpenMP_C_FLAGS})
endif()

if (PDM_ENABLE_CUDA)
  list(APPEND LINK_LIBRARIES pdm_gpu)
endif()

get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)

file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/meshes DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
if (PDM_ENABLE_EXTENSION_PDMA)
  file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/meshes DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../extensions/paradigma/test/")
endif()

# Test pdm_t_part_cube
# --------------

include(TestCreate)

if(PDM_ENABLE_UNIT_TEST)
  file(GLOB_RECURSE __test_files
       CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../src/*.test.cpp")
  test_cpp_unit_create(pdm_t_unit_test ${MPIEXEC_MAX_NUMPROCS} SOURCES ${__test_files} )
endif()

if (PDM_HAVE_PARMETIS OR PDM_HAVE_PTSCOTCH)
  test_c_create (pdm_t_part_dcube ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_part_dcube_perfo_btp ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_part_dcube_nodal ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_dmesh_nodal_to_dmesh 2)
  test_c_create (pdm_t_part_extension 2)
  test_c_create (pdm_t_partitioning_dcube ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_partitioning_dcube_nodal 2)
  test_c_create (pdm_t_multipart ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_global_mean ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_global_reduce ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_writer_ensight_3d ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_writer_ensight_moving_mesh_3d ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_plane_agglo 2)
  test_c_create (pdm_t_cellface_orient ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_gnum ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_gnum_location ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_gnum_with_merge ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_dist ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_dist_strip 2)
  test_c_create (pdm_t_dist_rocket 2)
  test_c_create (pdm_t_dist_axitransbump 2)
  test_c_create (pdm_t_knn_cube ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_create_edge_gnum 1)
  test_c_create (pdm_t_mesh_location_dcube ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_part1_to_selected_part2_from_mesh_location ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_mesh_location_disjoint_cubes ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_mesh_location_poly_surf ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_plane_meshes_intersect 2)
  test_c_create (pdm_t_dsurf_mesh_gen 2)
  test_c_create (pdm_t_polysurf_extrusion 2)
  test_c_create (pdm_t_dcube_nodal_gen 2)
  test_c_create (pdm_t_extract_part_from_part 1)
  test_c_create (pdm_t_compute_part_edges 1)
  test_c_create (pdm_t_extract_part 1)
  test_c_create (pdm_t_extract_part_nodal 1)
  test_c_create (pdm_t_extract_part_from_target 1)
  test_c_create (pdm_t_extract_surf_part 1)
  test_c_create (pdm_t_edge_upwind_and_downwind 1)
  test_c_create (pdm_t_overlay_from_vol_mesh 1)
endif()
test_c_create (pdm_t_cellface_orient_poly 1)
test_c_create (pdm_t_mesh_adapt 1)
test_c_create (pdm_t_ptb ${MPIEXEC_MAX_NUMPROCS})
test_c_create (pdm_t_par_octree ${MPIEXEC_MAX_NUMPROCS})
test_c_create (pdm_t_closest_points_small 1)
test_c_create (pdm_t_closest_points ${MPIEXEC_MAX_NUMPROCS})
test_c_create (pdm_t_distant_neighbor 1)
test_c_create (pdm_t_distant_neighbor_nomatch 1)
test_c_create (pdm_t_order 1)
test_c_create (pdm_t_gnum_from_hash_values 1)
test_c_create (pdm_t_plane_intersect_unit 2)
test_c_create (pdm_t_plane_meshes_intersect_unit 1)
test_c_create (pdm_t_intersect_line_box 1)
test_c_create (pdm_t_intersect_volume_box 1)
test_c_create (pdm_t_linear_programming_box_box 1)
test_c_create (pdm_t_box_tree_intersect_volume_box 1)
test_c_create (pdm_t_dbbtree_intersect_volume_box 1)
test_c_create (pdm_t_volume_intersect_back_mesh 1)
test_c_create (pdm_t_mesh_extract 1)
test_c_create (pdm_t_connect_zone 1)
test_c_create (pdm_t_part_to_block_geom 1)
if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)
  test_fortran_create (pdm_t_closest_points_f 2)
  test_fortran_create (pdm_t_mesh_location_f 1)
  test_fortran_create (pdm_t_part_to_block_f 1)
  test_fortran_create (pdm_t_block_to_part_f 1)
  test_fortran_create (pdm_t_dcube_gen_f 1)
  test_fortran_create (pdm_t_plane_intersect_unit_f 1)
  test_fortran_create (pdm_t_global_reduce_f 2)
  test_fortran_create (pdm_t_mpi_comm_f2c 2)
  test_fortran_create (pdm_t_dist_f 2)
  test_fortran_create (pdm_t_gnum_f 2)
  test_fortran_create (pdm_t_part_f 2)
  test_fortran_create (pdm_t_writer_f 2)
  test_fortran_create (pdm_t_mesh_check_f 1)
  test_fortran_create (pdm_t_part_to_part_f 2)
  test_fortran_create (pdm_t_gnum_location_f 2)
  test_fortran_create (pdm_t_part_extension_f 2)
  test_fortran_create (pdm_t_io  2)
endif()

test_c_create (pdm_t_closest_points_1_point 1)
test_c_create (pdm_t_dist_ridge 1)
# test_c_create (pdm_t_multi_block_merge 1)

test_c_create (pdm_t_ho_ordering 1)

#TODO: Not fully implemented 
#test_c_create (pdm_t_win_shared 1)

#test_c_create (test_mpi_morton 2)

test_c_create (pdm_t_point_in_polygon 1)
test_c_create (pdm_t_line_intersection 1)
test_c_create (pdm_t_mesh_check 1)
test_c_create (pdm_t_join_domains 2)
test_c_create (pdm_t_part_to_part 2)
test_c_create (pdm_t_ho_location 1)
test_c_create (pdm_t_subdivision_bezier_triangle 1)
test_c_create (pdm_t_projection_back_mesh 1)

test_c_create (pdm_t_ho_bezier_intersect 1)

test_c_create (pdm_t_sphere_vol_gen 1)
test_c_create (pdm_t_separate_groups 1)

add_subdirectory(tutorial)
