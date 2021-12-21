#include "doctest/extensions/doctest_mpi.h"
#include <array>
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_dcube_nodal_gen.h"

MPI_TEST_CASE("[1p] dcube_nodal_gen HEXA ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected_vol = {
    {1, 2, 11, 10, 4, 5, 14, 13, 2, 3, 12, 11, 5, 6, 15, 14, 4, 5, 14, 13, 7, 8, 17, 16, 5, 6, 15, 14, 8, 9, 18, 17, 10, 11, 20, 19, 13, 14, 23, 22, 11, 12, 21, 20, 14, 15, 24, 23, 13, 14, 23, 22, 16, 17, 26, 25, 14, 15, 24, 23, 17, 18, 27, 26} // HEXA
  };

  std::vector<std::vector<PDM_g_num_t>> connec_expected_surf = {
    {1, 10, 13, 4, 4, 13, 16, 7, 10, 19, 22, 13, 13, 22, 25, 16},     // QUAD IMIN
    {3, 12, 15, 6, 6, 15, 18, 9, 12, 21, 24, 15, 15, 24, 27, 18},     // QUAD IMAX
    {1, 2, 11, 10, 2, 3, 12, 11, 10, 11, 20, 19, 11, 12, 21, 20},     // QUAD JMIN
    {7, 8, 17, 16, 8, 9, 18, 17, 16, 17, 26, 25, 17, 18, 27, 26},     // QUAD JMAX
    {1, 4, 5, 2, 2, 5, 6, 3, 4, 7, 8, 5, 5, 8, 9, 6},                 // QUAD KMIN
    {19, 22, 23, 20, 20, 23, 24, 21, 22, 25, 26, 23, 23, 26, 27, 24}  // QUAD KMAX
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_vol = {
    {0, 8} // HEXA
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_surf = {
    {0, 4}, // QUAD IMIN
    {0, 4}, // QUAD IMAX
    {0, 4}, // QUAD JMIN
    {0, 4}, // QUAD JMAX
    {0, 4}, // QUAD KMIN
    {0, 4}  // QUAD KMAX
  };

  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_vol = {PDM_MESH_NODAL_HEXA8};

  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_surf = {PDM_MESH_NODAL_QUAD4, PDM_MESH_NODAL_QUAD4,
                                                                   PDM_MESH_NODAL_QUAD4, PDM_MESH_NODAL_QUAD4,
                                                                   PDM_MESH_NODAL_QUAD4, PDM_MESH_NODAL_QUAD4};

  int n_vtx_seg = 3;

  /*PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_HEXA8,
                                                      PDM_OWNERSHIP_KEEP);*/
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_HEXA8,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_ordering_set (dcube, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs ==  8);
  CHECK( n_face_abs == -1);
  CHECK( n_edge_abs == -1);
  CHECK( n_vtx_abs  == 27);

  int n_section_vol = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  CHECK( n_section_vol == 1); // HEXA + 6 * QUAD

  int* sections_id_vol = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  for(int i_section = 0; i_section < n_section_vol; ++i_section) {

    int id_section = sections_id_vol[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_vol[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_vol[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_vol [i_section], n_vtx_per_elmt * dn_elmt);
  }

  int n_section_surf = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  CHECK( n_section_surf == 6); // HEXA + 6 * QUAD

  int* sections_id_surf = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  for(int i_section = 0; i_section < n_section_surf; ++i_section) {

    int id_section = sections_id_surf[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_surf[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_surf[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_surf [i_section], n_vtx_per_elmt * dn_elmt);
  }

  PDM_dcube_nodal_gen_free(dcube);
}

MPI_TEST_CASE("[1p] dcube_nodal_gen PRISM ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected_vol = {
    {1, 5, 4, 10, 14, 13, 1, 2, 5, 10, 11, 14, 2, 6, 5, 11, 15, 14, 2, 3, 6, 11, 12, 15, 4, 8, 7, 13, 17, 16, 4, 5, 8, 13, 14, 17, 5, 9, 8, 14, 18, 17, 5, 6, 9, 14, 15, 18, 10, 14, 13, 19, 23, 22, 10, 11, 14, 19, 20, 23, 11, 15, 14, 20, 24, 23, 11, 12, 15, 20, 21, 24, 13, 17, 16, 22, 26, 25, 13, 14, 17, 22, 23, 26, 14, 18, 17, 23, 27, 26, 14, 15, 18, 23, 24, 27} // HEXA
  };

  std::vector<std::vector<PDM_g_num_t>> connec_expected_surf = {
    {1, 10, 13, 4, 4, 13, 16, 7, 10, 19, 22, 13, 13, 22, 25, 16}, // QUAD IMIN
    {3, 12, 15, 6, 6, 15, 18, 9, 12, 21, 24, 15, 15, 24, 27, 18}, // QUAD IMAX
    {1, 2, 11, 10, 2, 3, 12, 11, 10, 11, 20, 19, 11, 12, 21, 20}, // QUAD JMIN
    {7, 8, 17, 16, 8, 9, 18, 17, 16, 17, 26, 25, 17, 18, 27, 26}, // QUAD JMAX
    {1, 4, 5, 1, 5, 2, 2, 5, 6, 2, 6, 3, 4, 7, 8, 4, 8, 5, 5, 8, 9, 5, 9, 6}, // TRI KMIN
    {19, 22, 23, 19, 23, 20, 20, 23, 24, 20, 24, 21, 22, 25, 26, 22, 26, 23, 23, 26, 27, 23, 27, 24}  // TRI KMAX
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_vol = {
    {0, 16}, // PRISM
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_surf = {
    {0,  4}, // QUAD IMIN
    {0,  4}, // QUAD IMAX
    {0,  4}, // QUAD JMIN
    {0,  4}, // QUAD JMAX
    {0,  8}, // TRI KMIN
    {0,  8}  // TRI KMAX
  };
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_vol = {PDM_MESH_NODAL_PRISM6};
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_surf = {PDM_MESH_NODAL_QUAD4, PDM_MESH_NODAL_QUAD4,
                                                                   PDM_MESH_NODAL_QUAD4, PDM_MESH_NODAL_QUAD4,
                                                                   PDM_MESH_NODAL_TRIA3, PDM_MESH_NODAL_TRIA3};

  int n_vtx_seg = 3;

  /*PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_PRISM6,
                                                      PDM_OWNERSHIP_KEEP);*/
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_PRISM6,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_ordering_set (dcube, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs == 16);
  CHECK( n_face_abs == -1);
  CHECK( n_edge_abs == -1);
  CHECK( n_vtx_abs  == 27);
  int n_section_vol = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  CHECK( n_section_vol == 1); // HEXA + 6 * QUAD

  int* sections_id_vol = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  for(int i_section = 0; i_section < n_section_vol; ++i_section) {

    int id_section = sections_id_vol[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_vol[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_vol[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_vol [i_section], n_vtx_per_elmt * dn_elmt);
  }

  int n_section_surf = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  CHECK( n_section_surf == 6); // HEXA + 6 * QUAD

  int* sections_id_surf = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  for(int i_section = 0; i_section < n_section_surf; ++i_section) {

    int id_section = sections_id_surf[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_surf[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_surf[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_surf [i_section], n_vtx_per_elmt * dn_elmt);
  }

  PDM_dcube_nodal_gen_free(dcube);
}



MPI_TEST_CASE("[1p] dcube_nodal_gen TETRA ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected_vol = {
    {1, 2, 4, 10, 2, 5, 4, 14, 4, 10, 14, 13, 2, 10, 11, 14, 2, 4, 10, 14, 2, 6, 5, 14, 2, 3, 6, 12, 6, 12, 15, 14, 12, 11, 14, 2, 12, 6, 2, 14, 4, 8, 7, 16, 4, 5, 8, 14, 8, 14, 17, 16, 14, 13, 16, 4, 14, 8, 4, 16, 5, 6, 8, 14, 6, 9, 8, 18, 8, 14, 18, 17, 6, 14, 15, 18, 6, 8, 14, 18, 10, 14, 13, 22, 10, 11, 14, 20, 14, 20, 23, 22, 20, 19, 22, 10, 20, 14, 10, 22, 11, 12, 14, 20, 12, 15, 14, 24, 14, 20, 24, 23, 12, 20, 21, 24, 12, 14, 20, 24, 13, 14, 16, 22, 14, 17, 16, 26, 16, 22, 26, 25, 14, 22, 23, 26, 14, 16, 22, 26, 14, 18, 17, 26, 14, 15, 18, 24, 18, 24, 27, 26, 24, 23, 26, 14, 24, 18, 14, 26}
  };

  std::vector<std::vector<PDM_g_num_t>> connec_expected_surf = {
    {1, 4, 10, 4, 13, 10, 4, 7, 16, 4, 16, 13, 10, 13, 22, 10, 22, 19, 13, 16, 22, 16, 25, 22}, // TRI IMIN
    {3, 6, 12, 6, 15, 12, 6, 9, 18, 6, 18, 15, 12, 15, 24, 12, 24, 21, 15, 18, 24, 18, 27, 24}, // TRI IMAX
    {1, 2, 10, 2, 11, 10, 2, 3, 12, 2, 12, 11, 10, 11, 20, 10, 20, 19, 11, 12, 20, 12, 21, 20}, // TRI JMIN
    {7, 8, 16, 8, 17, 16, 8, 9, 18, 8, 18, 17, 16, 17, 26, 16, 26, 25, 17, 18, 26, 18, 27, 26}, // TRI JMAX
    {1, 4, 2, 4, 5, 2, 2, 5, 6, 2, 6, 3, 4, 7, 8, 4, 8, 5, 5, 8, 6, 8, 9, 6}, // TRI KMIN
    {19, 22, 20, 22, 23, 20, 20, 23, 24, 20, 24, 21, 22, 25, 26, 22, 26, 23, 23, 26, 24, 26, 27, 24}  // TRI KMAX
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_vol = {
    {0, 40}, // PRISM
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_surf = {
    {0,  8}, // QUAD IMIN
    {0,  8}, // QUAD IMAX
    {0,  8}, // QUAD JMIN
    {0,  8}, // QUAD JMAX
    {0,  8}, // TRI KMIN
    {0,  8}  // TRI KMAX
  };
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_vol = {PDM_MESH_NODAL_TETRA4};

  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_surf = {PDM_MESH_NODAL_TRIA3, PDM_MESH_NODAL_TRIA3,
                                                                   PDM_MESH_NODAL_TRIA3, PDM_MESH_NODAL_TRIA3,
                                                                   PDM_MESH_NODAL_TRIA3, PDM_MESH_NODAL_TRIA3};


  int n_vtx_seg = 3;

  /*PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_TETRA4,
                                                      PDM_OWNERSHIP_KEEP);*/
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_TETRA4,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_ordering_set (dcube, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs == 40);
  CHECK( n_face_abs == -1);
  CHECK( n_edge_abs == -1);
  CHECK( n_vtx_abs  == 27);

  int n_section_vol = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  CHECK( n_section_vol == 1); // HEXA + 6 * QUAD

  int* sections_id_vol = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  for(int i_section = 0; i_section < n_section_vol; ++i_section) {

    int id_section = sections_id_vol[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_vol[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_vol[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_vol [i_section], n_vtx_per_elmt * dn_elmt);
  }

  int n_section_surf = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  CHECK( n_section_surf == 6); // HEXA + 6 * QUAD

  int* sections_id_surf = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  for(int i_section = 0; i_section < n_section_surf; ++i_section) {

    int id_section = sections_id_surf[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_surf[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_surf[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_surf [i_section], n_vtx_per_elmt * dn_elmt);
  }


  PDM_dcube_nodal_gen_free(dcube);
}


MPI_TEST_CASE("[1p] dcube_nodal_gen QUAD ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected_surf = {
    {1,2,5,4,2,3,6,5,4,5,8,7,5,6,9,8}
  };

  std::vector<std::vector<PDM_g_num_t>> connec_expected_ridge = {
    {4,1,7,4}, // TRI IMIN
    {6,3,9,6}, // TRI IMAX
    {2,1,3,2}, // TRI JMIN
    {8,7,9,8}, // TRI JMAX
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_surf = {
    {0, 4} // QUAD
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_ridge = {
    {0, 2}, // BAR IMIN
    {0, 2}, // BAR IMAX
    {0, 2}, // BAR JMIN
    {0, 2}, // BAR JMAX
  };

  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_surf  = {PDM_MESH_NODAL_QUAD4};
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_ridge = {PDM_MESH_NODAL_BAR2, PDM_MESH_NODAL_BAR2,
                                                                    PDM_MESH_NODAL_BAR2, PDM_MESH_NODAL_BAR2};


  int n_vtx_seg = 3;

  /*PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_QUAD4,
                                                      PDM_OWNERSHIP_KEEP);*/
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_QUAD4,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_ordering_set (dcube, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs ==  4);
  CHECK( n_face_abs == -1);
  CHECK( n_edge_abs == -1);
  CHECK( n_vtx_abs  ==  9);

  int n_section_surf = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  CHECK( n_section_surf == 1); // HEXA + 6 * QUAD

  int* sections_id_surf = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  for(int i_section = 0; i_section < n_section_surf; ++i_section) {

    int id_section = sections_id_surf[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_surf[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_surf[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_surf [i_section], n_vtx_per_elmt * dn_elmt);
  }

  int n_section_ridge = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE);

  CHECK( n_section_ridge == 4); // HEXA + 6 * QUAD

  int* sections_id_ridge = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE);

  for(int i_section = 0; i_section < n_section_ridge; ++i_section) {

    int id_section = sections_id_ridge[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_ridge[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_ridge[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_ridge [i_section], n_vtx_per_elmt * dn_elmt);
  }

  PDM_dcube_nodal_gen_free(dcube);
}


MPI_TEST_CASE("[1p] dcube_nodal_gen TRI ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected_surf = {
    {1, 2, 4, 2, 5, 4, 2, 3, 5, 3, 6, 5, 4, 5, 7, 5, 8, 7, 5, 6, 8, 6, 9, 8 }
  };

  std::vector<std::vector<PDM_g_num_t>> connec_expected_ridge = {
    {4,1,7,4}, // TRI IMIN
    {6,3,9,6}, // TRI IMAX
    {2,1,3,2}, // TRI JMIN
    {8,7,9,8}, // TRI JMAX
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_surf = {
    {0, 8}, // TRI
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_ridge = {
    {0, 2}, // BAR IMIN
    {0, 2}, // BAR IMAX
    {0, 2}, // BAR JMIN
    {0, 2}, // BAR JMAX
  };
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_surf = {PDM_MESH_NODAL_TRIA3};
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_ridge = {PDM_MESH_NODAL_BAR2, PDM_MESH_NODAL_BAR2,
                                                                    PDM_MESH_NODAL_BAR2, PDM_MESH_NODAL_BAR2};


  int n_vtx_seg = 3;

  /*PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_TRIA3,
                                                      PDM_OWNERSHIP_KEEP);*/
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_TRIA3,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_ordering_set (dcube, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs ==  8);
  CHECK( n_face_abs == -1);
  CHECK( n_edge_abs == -1);
  CHECK( n_vtx_abs  ==  9);

  int n_section_surf = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  CHECK( n_section_surf == 1); // HEXA + 6 * QUAD

  int* sections_id_surf = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  for(int i_section = 0; i_section < n_section_surf; ++i_section) {

    int id_section = sections_id_surf[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_surf[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_surf[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_surf [i_section], n_vtx_per_elmt * dn_elmt);
  }

  int n_section_ridge = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE);

  CHECK( n_section_ridge == 4); // HEXA + 6 * QUAD

  int* sections_id_ridge = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE);

  for(int i_section = 0; i_section < n_section_ridge; ++i_section) {

    int id_section = sections_id_ridge[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_ridge[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_ridge[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_ridge [i_section], n_vtx_per_elmt * dn_elmt);
  }

  PDM_dcube_nodal_gen_free(dcube);
}
