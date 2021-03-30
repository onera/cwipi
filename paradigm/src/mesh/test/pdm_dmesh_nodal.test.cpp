#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"

/*
 *  Use case
 *
 *       3           4           2
 *       +-----------+-----------+
 *       |           |           |
 *       |           |           |
 *       |           |           |
 *       +-----------+-----------+
 *       5           1           6
 *
 *  A l'issu de l'algorithme on doit identifer 7 edges -->
 */

// PDM_g_num_t quad_connec[8] = {1, 4, 3, 5,  // First QUAD
//                               2, 4, 1, 6}; // Second QUAD

// For PDM_parent_elmt_find
// PDM_g_num_t line_connect[12] = {5, 1, 1, 6, 6, 2, 2, 4, 4, 3, 3, 5};


// MPI_TEST_CASE("[1p] dmesh_nodal 2D ",1) {
// }

// MPI_TEST_CASE("[2p] dmesh_nodal 2D ",2) {
// }
