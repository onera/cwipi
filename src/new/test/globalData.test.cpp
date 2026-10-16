#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "cwp.h"

#include <string.h>
#include <vector>

static double _max_error = 1e-6;

#define CHECK_EPS(array, array_expected, size) \
  for(int i = 0; i < size; ++i){ CHECK(array[i] == doctest::Approx(array_expected[i]).epsilon(_max_error)); }

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
  int        *n_part   = (int *) malloc(sizeof(int) * n_code);
  const char *cpl_name = "globalDataTest";

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

  // Allocate buffers
  const char *global_data_name1 = "lapin";
  size_t  s_entity1  = sizeof(double);
  int     stride1    = 2;
  int     n_entity1  = 1;
  double *send_data1 = (double *) malloc(s_entity1 * stride1 * n_entity1);
  double *recv_data1 = (double *) malloc(s_entity1 * stride1 * n_entity1);

  const char *global_data_name2 = "capybara";
  size_t  s_entity2  = sizeof(double);
  int     stride2    = 2;
  int     n_entity2  = 3;
  double *send_data2 = (double *) malloc(s_entity2 * stride2 * n_entity2);
  double *recv_data2 = (double *) malloc(s_entity2 * stride2 * n_entity2);

  // Fill buffers
  send_data1[0] = 1991.07;
  send_data1[1] = 2003.06;

  send_data2[0] = 1988.05;
  send_data2[1] = 1986.03;
  send_data2[2] = 1954.05;
  send_data2[3] = 1961.03;
  send_data2[4] = 1981.12;
  send_data2[5] = 1974.02;

  // Send/recv data
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      CWP_Global_data_issend(code_names[i_code],
                             cpl_name,
                             global_data_name1,
                             s_entity1,
                             stride1,
                             n_entity1,
                    (void *) send_data1);

      CWP_Global_data_irecv(code_names[i_code],
                            cpl_name,
                            global_data_name2,
                            s_entity2,
                            stride2,
                            n_entity2,
                   (void *) recv_data2);
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      CWP_Global_data_issend(code_names[i_code],
                             cpl_name,
                             global_data_name2,
                             s_entity2,
                             stride2,
                             n_entity2,
                    (void *) send_data2);

      CWP_Global_data_irecv(code_names[i_code],
                            cpl_name,
                            global_data_name1,
                            s_entity1,
                            stride1,
                            n_entity1,
                   (void *) recv_data1);
    }
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      CWP_Global_data_wait_issend(code_names[i_code],
                                  cpl_name,
                                  global_data_name1);

      CWP_Global_data_wait_irecv(code_names[i_code],
                                 cpl_name,
                                 global_data_name2);
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      CWP_Global_data_wait_issend(code_names[i_code],
                                  cpl_name,
                                  global_data_name2);

      CWP_Global_data_wait_irecv(code_names[i_code],
                                 cpl_name,
                                 global_data_name1);
    }
  }

  // Check consistency
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      CHECK_EPS(recv_data2, send_data2, stride2 * n_entity1);
    }
    if (strcmp(code_names[i_code], "code2") == 0) {
      CHECK_EPS(recv_data1, send_data1, stride1 * n_entity1);
    }
  }

  // Free all
  for (int i_code = 0; i_code < n_code; i_code++) {
    CWP_Cpl_del(code_names[i_code], cpl_name);
  }

  free(send_data1);
  free(recv_data1);
  free(send_data2);
  free(recv_data2);
  free(n_part);

  free(intra_comm);

  CWP_Finalize();
}

MPI_TEST_CASE("[globalData] - 1p", 1) {

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

MPI_TEST_CASE("[globalData] - 2p", 2) {

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