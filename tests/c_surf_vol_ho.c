#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cwipi.h"


static void
_gen_mesh_2d
(
  const int         order,
  int              *out_n_vtx,
  int              *out_n_elt,
  double          **out_vtx_coord,
  int             **out_connec_idx,
  int             **out_connec,
  int             **out_ijk_grid,
  cwipi_element_t  *out_elt_type,
  int              *out_n_vtx_per_elt
)
{
  // IJK grid
  int n_vtx_per_elt = (order+1)*(order+1);

  int *ijk_grid = (int *) malloc(sizeof(int) * n_vtx_per_elt * 2);
  int idx = 0;
  for (int j = 0; j <= order; j++) {
    for (int i = 0; i <= order; i++) {
      ijk_grid[idx++] = i;
      ijk_grid[idx++] = j;
    }
  }

  int n_elt_x = 5;
  int n_elt_y = 5;

  // Vertices
  int n_vtx = (n_elt_x*order + 1) * (n_elt_y*order + 1);
  double *vtx_coord = (double *) malloc(sizeof(double) * n_vtx * 3);

  double step_x = 1./(n_elt_x*order);
  double step_y = 1./(n_elt_y*order);

  idx = 0;
  for (int j = 0; j <= n_elt_y*order; j++) {
    for (int i = 0; i <= n_elt_x*order; i++) {
      vtx_coord[idx++] = i*step_x - 0.5;
      vtx_coord[idx++] = j*step_y - 0.5;
      vtx_coord[idx++] = 0;
    }
  }

  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    vtx_coord[3*i_vtx+2] = 0.5*(vtx_coord[3*i_vtx]*vtx_coord[3*i_vtx] + vtx_coord[3*i_vtx+1]*vtx_coord[3*i_vtx+1]) - 0.25;
  }

  // Elements
  int n_elt = n_elt_x*n_elt_y;
  int *connec_idx = (int *) malloc(sizeof(int) * (n_elt+1));
  int *connec     = (int *) malloc(sizeof(int) * n_elt * n_vtx_per_elt);

  connec_idx[0] = 0;
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    connec_idx[i_elt+1] = connec_idx[i_elt] + n_vtx_per_elt;
  }

  idx = 0;
  for (int j = 0; j < n_elt_y; j++) {
    for (int i = 0; i < n_elt_x; i++) {
      for (int jj = 0; jj <= order; jj++) {
        for (int ii = 0; ii <= order; ii++) {
          connec[idx++] = 1 + (j*order + jj)*(n_elt_x*order + 1) + i*order + ii;
        }
      }
    }
  }

  // Outputs
  *out_n_vtx         = n_vtx;
  *out_n_elt         = n_elt;
  *out_vtx_coord     = vtx_coord;
  *out_connec_idx    = connec_idx;
  *out_connec        = connec;
  *out_ijk_grid      = ijk_grid;
  *out_elt_type      = CWIPI_FACE_QUADHO;
  *out_n_vtx_per_elt = n_vtx_per_elt;
}



static void
_gen_mesh_3d
(
  const int         order,
  int              *out_n_vtx,
  int              *out_n_elt,
  double          **out_vtx_coord,
  int             **out_connec_idx,
  int             **out_connec,
  int             **out_ijk_grid,
  cwipi_element_t  *out_elt_type,
  int              *out_n_vtx_per_elt
)
{
  // IJK grid
  int n_vtx_per_elt = (order+1)*(order+1)*(order+1);

  int *ijk_grid = (int *) malloc(sizeof(int) * n_vtx_per_elt * 3);
  int idx = 0;
  for (int k = 0; k <= order; k++) {
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order; i++) {
        ijk_grid[idx++] = i;
        ijk_grid[idx++] = j;
        ijk_grid[idx++] = k;
      }
    }
  }

  int n_elt_x = 3;
  int n_elt_y = 3;
  int n_elt_z = 3;

  // Vertices
  int n_vtx = (n_elt_x*order + 1) * (n_elt_y*order + 1) * (n_elt_z*order + 1);
  double *vtx_coord = (double *) malloc(sizeof(double) * n_vtx * 3);

  double step_x = 1./(n_elt_x*order);
  double step_y = 1./(n_elt_y*order);
  double step_z = 1./(n_elt_z*order);

  idx = 0;
  for (int k = 0; k <= n_elt_z*order; k++) {
    for (int j = 0; j <= n_elt_y*order; j++) {
      for (int i = 0; i <= n_elt_x*order; i++) {
        vtx_coord[idx++] = i*step_x - 0.5;
        vtx_coord[idx++] = j*step_y - 0.5;
        vtx_coord[idx++] = k*step_z - 0.5;
      }
    }
  }

  // Elements
  int n_elt = n_elt_x*n_elt_y*n_elt_z;
  int *connec_idx = (int *) malloc(sizeof(int) * (n_elt+1));
  int *connec     = (int *) malloc(sizeof(int) * n_elt * n_vtx_per_elt);

  connec_idx[0] = 0;
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    connec_idx[i_elt+1] = connec_idx[i_elt] + n_vtx_per_elt;
  }

  idx = 0;
  for (int k = 0; k < n_elt_z; k++) {
    for (int j = 0; j < n_elt_y; j++) {
      for (int i = 0; i < n_elt_x; i++) {
        for (int kk = 0; kk <= order; kk++) {
          for (int jj = 0; jj <= order; jj++) {
            for (int ii = 0; ii <= order; ii++) {
              connec[idx++] = 1 +
                (k*order + kk)*(n_elt_x*order + 1)*(n_elt_y*order+1) +
                (j*order + jj)*(n_elt_x*order + 1)                   +
                 i*order + ii;
            }
          }
        }
      }
    }
  }

  // Outputs
  *out_n_vtx         = n_vtx;
  *out_n_elt         = n_elt;
  *out_vtx_coord     = vtx_coord;
  *out_connec_idx    = connec_idx;
  *out_connec        = connec;
  *out_ijk_grid      = ijk_grid;
  *out_elt_type      = CWIPI_CELL_HEXAHO;
  *out_n_vtx_per_elt = n_vtx_per_elt;
}


/**
 * One-way coupling with two high-order meshes:
 * code1 (surface, order3)  <--  code2 (volume, order 2)
 */
int
main
(
  int   argc,
  char *argv[]
)
{
  // Initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  if (n_rank != 2) {
    if (i_rank == 0) {
      printf("Error: This test can only be run with 2 MPI tasks\n");
    }
    exit(EXIT_FAILURE);
  }


  const char *code_name;
  const char *cpl_code_name;

  if (i_rank < n_rank/2) {
    code_name     = "code1";
    cpl_code_name = "code2";
  }
  else {
    code_name     = "code2";
    cpl_code_name = "code1";
  }

  // Initialize CWIPI
  MPI_Comm intra_comm;
  cwipi_init(MPI_COMM_WORLD,
             code_name,
             &intra_comm);

  // Generate mesh
  int              dim;
  int              n_vtx;
  int              n_elt;
  int              order;
  double          *vtx_coord;
  int             *connec_idx;
  int             *connec;
  int             *ijk_grid;
  cwipi_element_t  elt_type;
  int              n_vtx_per_elt;
  if (i_rank < n_rank/2) {
    dim   = 2;
    order = 3;
    _gen_mesh_2d(order,
                 &n_vtx,
                 &n_elt,
                 &vtx_coord,
                 &connec_idx,
                 &connec,
                 &ijk_grid,
                 &elt_type,
                 &n_vtx_per_elt);
  }
  else {
    dim   = 3;
    order = 2;
    _gen_mesh_3d(order,
                 &n_vtx,
                 &n_elt,
                 &vtx_coord,
                 &connec_idx,
                 &connec,
                 &ijk_grid,
                 &elt_type,
                 &n_vtx_per_elt);
  }

  // Create coupling environment
  cwipi_create_coupling("surf_vol_ho",                             // Coupling id
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                        cpl_code_name,                             // Coupled application id
                        dim,                                       // Geometric entities dimension
                        0.1,                                       // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        CWIPI_SOLVER_CELL_VERTEX,                  // Solver type
                        0,                                         // Postprocessing frequency
                        "EnSight Gold",                            // Postprocessing format
                        "text");                                   // Postprocessing option

  // Set interface mesh
  cwipi_ho_define_mesh("surf_vol_ho",
                       n_vtx,
                       n_elt,
                       order,
                       vtx_coord,
                       connec_idx,
                       connec);

  cwipi_ho_ordering_from_IJK_set("surf_vol_ho",
                                 elt_type,
                                 n_vtx_per_elt,
                                 ijk_grid);


  // Exchange field
  double *send_coord_x = NULL;
  double *recv_coord_x = NULL;

  if (i_rank < n_rank/2) {
    recv_coord_x = (double *) malloc(sizeof(double) * n_vtx);
  }
  else {
    send_coord_x = (double *) malloc(sizeof(double) * n_vtx);
    for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
      send_coord_x[i_vtx] = vtx_coord[3*i_vtx];
    }
  }

  int n_unlocated = 0;
  cwipi_exchange("surf_vol_ho",
                 "exchange",
                 1,
                 1,
                 0.,
                 "coord_x",
                 send_coord_x,
                 "coord_x",
                 recv_coord_x,
                 &n_unlocated);

  int error = 0;
  if (i_rank < n_rank/2) {
    if (n_unlocated > 0) {
      printf("All points should have been located (n_unlocated = %d / %d)\n", n_unlocated, n_vtx);
      exit(EXIT_FAILURE);
    }
    double max_err = 0;
    for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
      double err = fabs(recv_coord_x[i_vtx] - vtx_coord[3*i_vtx]);
      if (err > max_err) {
        max_err = err;
      }
    }
    printf("max_err = %e\n", max_err);
    error = (max_err > 1e-14);
  }


  // Finalize CWIPI
  cwipi_finalize();

  // Free memory
  free(vtx_coord);
  free(connec_idx);
  free(connec);
  free(ijk_grid);
  free(send_coord_x);
  free(recv_coord_x);

  // Finalize MPI
  MPI_Finalize();

  return error;
}