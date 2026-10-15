#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "cwp.h"

#include <string.h>
#include <vector>

static double _max_error = 1e-6;

#define CHECK_EPS(array, array_expected, size) \
  for(int i = 0; i < size; ++i){ CHECK(array[i] == doctest::Approx(array_expected[i]).epsilon(_max_error)); }

static void
_gen_mesh
(
  const MPI_Comm             comm,           // MPI communicator
  const int                  n_part,         // Number of partitions
  const CWPT_split_dual_t    part_method,    // Partitioning method
  const CWP_g_num_t          n_vtx_seg,      // Number of vertices on square side
        int                **pn_face,        // Number of faces
        int                **pn_vtx,         // Number of vertices
        int               ***pface_vtx_idx,  // Face vertex connectivity index
        int               ***pface_vtx,      // Face vertex connectivity
        double            ***pvtx_coord,     // Vertex coordinates
        CWP_g_num_t       ***pface_ln_to_gn, // Face global numbering
        CWP_g_num_t       ***pvtx_ln_to_gn   // Vertex global numbering
)
{
  int          *pn_edge        = nullptr;
  int         **pface_edge     = nullptr;
  int         **pedge_vtx      = nullptr;
  CWP_g_num_t **pedge_ln_to_gn = nullptr;
  CWPT_generate_mesh_rectangle_ngon(comm,
                                    static_cast<CWPT_Mesh_nodal_elt_t>(CWP_BLOCK_FACE_QUAD4),
                                    0.,
                                    0.,
                                    0.,
                                    1.,
                                    1.,
                                    n_vtx_seg,
                                    n_vtx_seg,
                                    n_part,
                                    part_method,
                                    0,
                                    pn_vtx,
                                    &pn_edge,
                                    pn_face,
                                    pvtx_coord,
                                    &pedge_vtx,
                                    pface_vtx_idx,
                                    &pface_edge,
                                    pface_vtx,
                                    pvtx_ln_to_gn,
                                    &pedge_ln_to_gn,
                                    pface_ln_to_gn);

  for (int i = 0; i < n_part; i++) {
    free(pface_edge    [i]);
    free(pedge_vtx     [i]);
    free(pedge_ln_to_gn[i]);
  }
  free(pn_edge);
  free(pface_edge    );
  free(pedge_vtx     );
  free(pedge_ln_to_gn);
}

static void
_get_values
(
  MPI_Comm   comm,       // MPI communicator
  CWP_Comm_t part1,      // Type of partitioning support for code 1
  CWP_Comm_t part2,      // Type of partitioning support for code 2
  bool       share_ranks // Share ranks for codes 1 and 2
)
{
  int i_rank;
  MPI_Comm_rank(comm, &i_rank);
  int n_rank;
  MPI_Comm_size(comm, &n_rank);
  assert(n_rank <= 2);

  std::vector<const char *> code_names;
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  // 1 rank = 2 codes on current rank
  if (n_rank == 1) {
    code_names = {"code1", "code2"};
  }
  // 2 ranks = dispatch codes to share ranks if share_ranks == true
  else {
    if (share_ranks) {
      code_names = {"code1", "code2"};
    }
    else {
      if (i_rank == 1) {
        code_names.push_back("code1");
      }
      else { // i_rank == 0
        code_names.push_back("code2");
      }
    }
  }
  int n_code = static_cast<int>(code_names.size());

  MPI_Comm *intra_comm = (MPI_Comm *) malloc(sizeof(MPI_Comm) * n_code);
  CWP_Init(comm,
           n_code,
           code_names.data(),
           is_active_rank,
           intra_comm);

  // Create coupling
  int        *n_part    = (int *) malloc(sizeof(int) * n_code);
  int         n_vtx_seg = 9;
  const char *cpl_name  = "partDataTest";

  for (int i_code = 0; i_code < n_code; i_code++) {

    const char *coupled_code_name;
    CWP_Comm_t  part;

    if (strcmp(code_names[i_code], "code1") == 0) {
      coupled_code_name = "code2";
      part           = part1;
      n_part[i_code] = 1;
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      coupled_code_name = "code1";
      part           = part2;
      n_part[i_code] = 2;
    }

    CWP_Cpl_create(code_names[i_code],
                   cpl_name,
                   coupled_code_name,
                   CWP_INTERFACE_SURFACE,
                   part,
                   CWP_SPATIAL_INTERP_FROM_IDENTITY,
                   n_part    [i_code],
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

  }

  // Define mesh
  int          **pn_face        = (int          **) malloc(sizeof(int          *) * n_code);
  int          **pn_vtx         = (int          **) malloc(sizeof(int          *) * n_code);
  int         ***pface_vtx_idx  = (int         ***) malloc(sizeof(int         **) * n_code);
  int         ***pface_vtx      = (int         ***) malloc(sizeof(int         **) * n_code);
  double      ***pvtx_coord     = (double      ***) malloc(sizeof(double      **) * n_code);
  CWP_g_num_t ***pface_ln_to_gn = (CWP_g_num_t ***) malloc(sizeof(CWP_g_num_t **) * n_code);
  CWP_g_num_t ***pvtx_ln_to_gn  = (CWP_g_num_t ***) malloc(sizeof(CWP_g_num_t **) * n_code);

  for (int i_code = 0; i_code < n_code; i_code++) {

    // Generate mesh on intra_comm communicator or local rank communicator
    // if no partitioning is required
    MPI_Comm mesh_comm;
    if ((part1 == CWP_COMM_PAR_WITHOUT_PART && strcmp(code_names[i_code], "code1") == 0) ||
        (part2 == CWP_COMM_PAR_WITHOUT_PART && strcmp(code_names[i_code], "code2") == 0)) {
      int color;
      MPI_Comm_rank(intra_comm[i_code], &color);
      MPI_Comm_split(intra_comm[i_code], color, color, &mesh_comm);
    }
    else {
      MPI_Comm_dup(intra_comm[i_code], &mesh_comm);
    }
    _gen_mesh(mesh_comm,
              n_part         [i_code],
              CWPT_SPLIT_DUAL_WITH_HILBERT,
              n_vtx_seg,
              &pn_face       [i_code],
              &pn_vtx        [i_code],
              &pface_vtx_idx [i_code],
              &pface_vtx     [i_code],
              &pvtx_coord    [i_code],
              &pface_ln_to_gn[i_code],
              &pvtx_ln_to_gn [i_code]);
    MPI_Comm_free(&mesh_comm);

  }

  // Allocate exchange buffers and fill send buffers
  double      ***send_val1 = (double      ***) malloc(sizeof(double      **) * n_code);
  double      ***recv_val1 = (double      ***) malloc(sizeof(double      **) * n_code);
  double      ***send_val2 = (double      ***) malloc(sizeof(double      **) * n_code);
  double      ***recv_val2 = (double      ***) malloc(sizeof(double      **) * n_code);
  CWP_g_num_t ***send_val3 = (CWP_g_num_t ***) malloc(sizeof(CWP_g_num_t **) * n_code);
  CWP_g_num_t ***recv_val3 = (CWP_g_num_t ***) malloc(sizeof(CWP_g_num_t **) * n_code);
  CWP_g_num_t ***send_val4 = (CWP_g_num_t ***) malloc(sizeof(CWP_g_num_t **) * n_code);
  CWP_g_num_t ***recv_val4 = (CWP_g_num_t ***) malloc(sizeof(CWP_g_num_t **) * n_code);

  for (int i_code = 0; i_code < n_code; i_code++) {

    send_val1[i_code] = (double      **) malloc(sizeof(double      *) * n_part[i_code]);
    recv_val1[i_code] = (double      **) malloc(sizeof(double      *) * n_part[i_code]);
    send_val2[i_code] = (double      **) malloc(sizeof(double      *) * n_part[i_code]);
    recv_val2[i_code] = (double      **) malloc(sizeof(double      *) * n_part[i_code]);
    send_val3[i_code] = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * n_part[i_code]);
    recv_val3[i_code] = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * n_part[i_code]);
    send_val4[i_code] = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * n_part[i_code]);
    recv_val4[i_code] = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * n_part[i_code]);

    for (int i_part = 0; i_part < n_part[i_code]; i_part++) {

      send_val1[i_code][i_part] = pvtx_coord    [i_code][i_part];
      send_val2[i_code][i_part] = (double      *) malloc(sizeof(double     ) *     pn_vtx [i_code][i_part]);
      send_val3[i_code][i_part] = pface_ln_to_gn[i_code][i_part];
      send_val4[i_code][i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * 2 * pn_face[i_code][i_part]);
      recv_val1[i_code][i_part] = (double      *) malloc(sizeof(double     ) * 3 * pn_vtx [i_code][i_part]);
      recv_val2[i_code][i_part] = (double      *) malloc(sizeof(double     ) *     pn_vtx [i_code][i_part]);
      recv_val3[i_code][i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) *     pn_face[i_code][i_part]);
      recv_val4[i_code][i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * 2 * pn_face[i_code][i_part]);

      for (int i_vtx = 0; i_vtx < pn_vtx[i_code][i_part]; i_vtx++) {
        send_val2[i_code][i_part][i_vtx] = static_cast<double>(pvtx_ln_to_gn[i_code][i_part][i_vtx]);
      }
      for (int i_face = 0; i_face < pn_face[i_code][i_part]; i_face++) {
        send_val4[i_code][i_part][2 * i_face    ] = 2 * pface_ln_to_gn[i_code][i_part][i_face];
        send_val4[i_code][i_part][2 * i_face + 1] = 3 * pface_ln_to_gn[i_code][i_part][i_face];
      }

    }

  }

  // Create partData objects
  const char *part_data_name1 = "part_data_vtx";
  const char *part_data_name2 = "part_data_face";

  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      CWP_Part_data_create(code_names   [i_code],
                           cpl_name,
                           part_data_name1,
                           CWP_PARTDATA_SEND,
                           pvtx_ln_to_gn[i_code],
                           pn_vtx       [i_code],
                           n_part       [i_code]);

      CWP_Part_data_create(code_names    [i_code],
                           cpl_name,
                           part_data_name2,
                           CWP_PARTDATA_RECV,
                           pface_ln_to_gn[i_code],
                           pn_face       [i_code],
                           n_part        [i_code]);
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      CWP_Part_data_create(code_names   [i_code],
                           cpl_name,
                           part_data_name1,
                           CWP_PARTDATA_RECV,
                           pvtx_ln_to_gn[i_code],
                           pn_vtx       [i_code],
                           n_part       [i_code]);

      CWP_Part_data_create(code_names    [i_code],
                           cpl_name,
                           part_data_name2,
                           CWP_PARTDATA_SEND,
                           pface_ln_to_gn[i_code],
                           pn_face       [i_code],
                           n_part        [i_code]);
    }
  }

  // Send/recv data with matching exch_id (used to compute communication MPI tag)
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      CWP_Part_data_issend(code_names[i_code],
                           cpl_name,
                           part_data_name1,
                           0,
                           sizeof(double),
                           3,
                 (void **) send_val1 [i_code]);

      CWP_Part_data_issend(code_names[i_code],
                           cpl_name,
                           part_data_name1,
                           1,
                           sizeof(double),
                           1,
                 (void **) send_val2 [i_code]);

      CWP_Part_data_irecv (code_names[i_code],
                           cpl_name,
                           part_data_name2,
                           2,
                           sizeof(CWP_g_num_t),
                           1,
                 (void **) recv_val3 [i_code]);

      CWP_Part_data_irecv (code_names[i_code],
                           cpl_name,
                           part_data_name2,
                           3,
                           sizeof(CWP_g_num_t),
                           2,
                 (void **) recv_val4 [i_code]);
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      CWP_Part_data_irecv (code_names[i_code],
                           cpl_name,
                           part_data_name1,
                           0,
                           sizeof(double),
                           3,
                 (void **) recv_val1 [i_code]);

      CWP_Part_data_irecv (code_names[i_code],
                           cpl_name,
                           part_data_name1,
                           1,
                           sizeof(double),
                           1,
                 (void **) recv_val2 [i_code]);

      CWP_Part_data_issend(code_names[i_code],
                           cpl_name,
                           part_data_name2,
                           2,
                           sizeof(CWP_g_num_t),
                           1,
                 (void **) send_val3 [i_code]);

      CWP_Part_data_issend(code_names[i_code],
                           cpl_name,
                           part_data_name2,
                           3,
                           sizeof(CWP_g_num_t),
                           2,
                 (void **) send_val4 [i_code]);
    }
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      CWP_Part_data_wait_issend(code_names[i_code],
                                cpl_name,
                                part_data_name1,
                                0);

      CWP_Part_data_wait_issend(code_names[i_code],
                                cpl_name,
                                part_data_name1,
                                1);

      CWP_Part_data_wait_irecv (code_names[i_code],
                                cpl_name,
                                part_data_name2,
                                2);

      CWP_Part_data_wait_irecv (code_names[i_code],
                                cpl_name,
                                part_data_name2,
                                3);
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      CWP_Part_data_wait_irecv (code_names[i_code],
                                cpl_name,
                                part_data_name1,
                                0);

      CWP_Part_data_wait_irecv (code_names[i_code],
                                cpl_name,
                                part_data_name1,
                                1);

      CWP_Part_data_wait_issend(code_names[i_code],
                                cpl_name,
                                part_data_name2,
                                2);

      CWP_Part_data_wait_issend(code_names[i_code],
                                cpl_name,
                                part_data_name2,
                                3);
    }
  }

  // Check consistency
  for (int i_code = 0; i_code < n_code; i_code++) {
    for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
      if (strcmp(code_names[i_code], "code1") == 0) {
        CHECK_EPS(recv_val3[i_code][i_part], send_val3[i_code][i_part],     pn_face[i_code][i_part]);
        CHECK_EPS(recv_val4[i_code][i_part], send_val4[i_code][i_part], 2 * pn_face[i_code][i_part]);
      }

      if (strcmp(code_names[i_code], "code2") == 0) {
        CHECK_EPS(recv_val1[i_code][i_part], send_val1[i_code][i_part], 3 * pn_vtx[i_code][i_part]);
        CHECK_EPS(recv_val2[i_code][i_part], send_val2[i_code][i_part],     pn_vtx[i_code][i_part]);
      }
    }
  }

  // Free all
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      CWP_Part_data_del(code_names[i_code],
                        cpl_name,
                        part_data_name1);

      CWP_Part_data_del(code_names[i_code],
                        cpl_name,
                        part_data_name2);
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      CWP_Part_data_del(code_names[i_code],
                        cpl_name,
                        part_data_name1);

      CWP_Part_data_del(code_names[i_code],
                        cpl_name,
                        part_data_name2);
    }

    CWP_Cpl_del(code_names[i_code], cpl_name);
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
      free(pface_vtx_idx [i_code][i_part]);
      free(pface_vtx     [i_code][i_part]);
      free(pvtx_coord    [i_code][i_part]);
      free(pface_ln_to_gn[i_code][i_part]);
      free(pvtx_ln_to_gn [i_code][i_part]);
      free(recv_val1     [i_code][i_part]);
      free(send_val2     [i_code][i_part]);
      free(recv_val2     [i_code][i_part]);
      free(recv_val3     [i_code][i_part]);
      free(send_val4     [i_code][i_part]);
      free(recv_val4     [i_code][i_part]);
    }
    free(pn_face       [i_code]);
    free(pn_vtx        [i_code]);
    free(pface_vtx_idx [i_code]);
    free(pface_vtx     [i_code]);
    free(pvtx_coord    [i_code]);
    free(pface_ln_to_gn[i_code]);
    free(pvtx_ln_to_gn [i_code]);
    free(send_val1     [i_code]);
    free(recv_val1     [i_code]);
    free(send_val2     [i_code]);
    free(recv_val2     [i_code]);
    free(send_val3     [i_code]);
    free(recv_val3     [i_code]);
    free(send_val4     [i_code]);
    free(recv_val4     [i_code]);
  }
  free(pn_face);
  free(pn_vtx);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pvtx_coord);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(n_part);
  free(send_val1);
  free(recv_val1);
  free(send_val2);
  free(recv_val2);
  free(send_val3);
  free(recv_val3);
  free(send_val4);
  free(recv_val4);

  free(intra_comm);

  CWP_Finalize();
}

MPI_TEST_CASE("[partData] - 1p", 1) {

  SUBCASE("WP - WP - no-shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITH_PART,
                CWP_COMM_PAR_WITH_PART,
                false);
  }

  SUBCASE("WoP - WoP - no-shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITHOUT_PART,
                CWP_COMM_PAR_WITHOUT_PART,
                false);
  }

  SUBCASE("WP - WoP - no-shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITH_PART,
                CWP_COMM_PAR_WITHOUT_PART,
                false);
  }

  SUBCASE("WoP - WP - no-shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITHOUT_PART,
                CWP_COMM_PAR_WITH_PART,
                false);
  }

}

MPI_TEST_CASE("[partData] - 2p", 2) {

  SUBCASE("WP - WP - no-shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITH_PART,
                CWP_COMM_PAR_WITH_PART,
                false);
  }

  SUBCASE("WP - WoP - no-shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITH_PART,
                CWP_COMM_PAR_WITHOUT_PART,
                false);
  }

  SUBCASE("WoP - WoP - no-shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITHOUT_PART,
                CWP_COMM_PAR_WITHOUT_PART,
                false);
  }

  SUBCASE("WoP - WP - no-shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITHOUT_PART,
                CWP_COMM_PAR_WITH_PART,
                false);
  }

  SUBCASE("WP - WP - shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITH_PART,
                CWP_COMM_PAR_WITH_PART,
                true);
  }

  SUBCASE("WP - WoP - shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITH_PART,
                CWP_COMM_PAR_WITHOUT_PART,
                true);
  }

  SUBCASE("WoP - WP - shared") {
    _get_values(test_comm,
                CWP_COMM_PAR_WITHOUT_PART,
                CWP_COMM_PAR_WITH_PART,
                true);
  }

}