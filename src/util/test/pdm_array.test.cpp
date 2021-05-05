#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"

#include "pdm_array.h"

MPI_TEST_CASE("[1p] _PDM_array_zeros", 1) {
  int *array = PDM_array_zeros_int(5);
  for (int i = 0; i < 5; i++)
    CHECK(array[i] == 0);
  free(array);
}

MPI_TEST_CASE("[1p] _PDM_array_const", 1) {
  int *array = PDM_array_const_int(7, 42);
  for (int i = 0; i < 7; i++)
    CHECK(array[i] == 42);
  free(array);

  PDM_g_num_t *array_gnum = PDM_array_const_gnum(10, 33550336);
  for (int i = 0; i < 10; i++)
    CHECK(array_gnum[i] == 33550336);
  free(array_gnum);
}

MPI_TEST_CASE("[1p] _PDM_array_reset", 1) {
  int *array = (int *) malloc(5*sizeof(int));
  for (int i = 0; i < 5; i++)
    array[i] == i;

  PDM_array_reset_int(array, 5, 0);
  for (int i = 0; i < 5; i++)
    CHECK(array[i] == 0);
  PDM_array_reset_int(array, 5, -999);
  for (int i = 0; i < 5; i++)
    CHECK(array[i] == -999);

  free(array);

  PDM_g_num_t *array_gnum = (PDM_g_num_t *) malloc(5*sizeof(PDM_g_num_t));
  for (int i = 0; i < 5; i++)
    array_gnum[i] == i;

  PDM_array_reset_gnum(array_gnum, 5, 0);
  for (int i = 0; i < 5; i++)
    CHECK(array_gnum[i] == 0);
  free(array_gnum);
}

MPI_TEST_CASE("[1p] _PDM_array_new_idx_from_sizes", 1) {
  int size_array[] = {5, 5, 2, 5};
  int *idx_array = PDM_array_new_idx_from_sizes_int(size_array, 4);
  int expected_idx_array[] = {0, 5, 10, 12, 17};
  CHECK_EQ_C_ARRAY(idx_array, expected_idx_array, 4+1);
  free(idx_array);

  PDM_g_num_t *idx_array_gnum = PDM_array_new_idx_from_sizes_gnum(size_array, 4);
  PDM_g_num_t expected_idx_array_gnum[] = {0, 5, 10, 12, 17};
  CHECK_EQ_C_ARRAY(idx_array_gnum, expected_idx_array_gnum, 4+1);
  free(idx_array_gnum);

  int empty_array[] = {};
  idx_array = PDM_array_new_idx_from_sizes_int(empty_array, 0);
  CHECK(idx_array[0] == 0);
  free(idx_array);
}

MPI_TEST_CASE("[1p] _PDM_array_idx_from_sizes", 1) {
  int size_array[] = {5, 5, 2, 5};
  int idx_array[] = {-1,-1,-1,-1,-1};
  PDM_array_idx_from_sizes_int(size_array, 4, idx_array);
  int expected_idx_array[] = {0, 5, 10, 12, 17};
  CHECK_EQ_C_ARRAY(idx_array, expected_idx_array, 4+1);

  PDM_g_num_t idx_array_gnum[] = {-1,-1,-1,-1,-1};
  PDM_array_idx_from_sizes_gnum(size_array, 4, idx_array_gnum);
  int expected_idx_array_gnum[] = {0, 5, 10, 12, 17};
  CHECK_EQ_C_ARRAY(idx_array_gnum, expected_idx_array_gnum, 4+1);
}