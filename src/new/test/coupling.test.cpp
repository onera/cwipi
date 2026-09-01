#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "cwp.h"

MPI_TEST_CASE("[coupling] - 2p - simple ", 2) {

  int i_rank;
  MPI_Comm_rank(test_comm, &i_rank);

  int n_code = 0;
  const char **code_names = NULL;
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  if (i_rank == 0) {
    n_code = 1;
    code_names = (const char **) malloc(sizeof(char *) * n_code);
    code_names[0] = "code1";
  }
  else if (i_rank == 1) {
    n_code = 1;
    code_names = (const char **) malloc(sizeof(char *) * n_code);
    code_names[0] = "code2";
  }

  MPI_Comm *localComm = (MPI_Comm *) malloc(sizeof(MPI_Comm) * n_code);

  CWP_Init(test_comm,
           n_code,
           (const char **) code_names,
           is_active_rank,
           localComm);

  free(localComm);
  free(code_names);

  CWP_Finalize();

}
