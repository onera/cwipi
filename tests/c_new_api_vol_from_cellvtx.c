/*
  This file is part of the CWIPI library.

  Copyright (C) 2024  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cwipi.h"
#include "cwp.h"
#include "cwp_priv.h"


/*----------------------------------------------------------------------
 *
 * Main : volume coupling test : c_new_api_vol_from_cellvtx
 *
 *---------------------------------------------------------------------*/
int main
(
  int   argc,
  char *argv[]
)
{
  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  assert(n_rank == 2);

  /* Initialize CWIPI */
  int n_code = 1;
  const char *code_name;
  const char *coupled_code_name;

  int code_id = (i_rank%2);

  if (code_id == 0) {
    code_name = "code1";
    coupled_code_name = "code2";
  }
  else {
    code_name = "code2";
    coupled_code_name = "code1";
  }

  CWP_Status_t is_active_rank = CWP_STATUS_ON;
  MPI_Comm     intra_comm;

  CWP_Init(comm,
           n_code,
           &code_name,
           is_active_rank,
           &intra_comm);


  /* Create Coupling */
  int n_part = 1;
  const char *cpl_name = "c_new_api_vol_from_cellvtx";
  CWP_Cpl_create(code_name,
                 cpl_name,
                 coupled_code_name,
                 CWP_INTERFACE_VOLUME,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  CWP_Visu_set(code_name,
               cpl_name,
               1,
               CWP_VISU_FORMAT_ENSIGHT,
               "text");

  /* Define interface mesh */
  const int n_vtx = 27;
  double vtx_coord[27*3];
  int idx = 0;
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        vtx_coord[idx++] = i - 1;
        vtx_coord[idx++] = j - 1;
        vtx_coord[idx++] = k - 1;
      }
    }
  }

  if (code_id == 0) {
    // Rotate 90° around z-axis
    for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
      double x = vtx_coord[3*i_vtx];
      vtx_coord[3*i_vtx  ] = -vtx_coord[3*i_vtx+1];
      vtx_coord[3*i_vtx+1] = x;
    }
  }

  int n_cell = 22;
  int cell_vtx_idx[23] = {0, 8, 16, 22, 28, 34, 40, 45, 50, 55, 60, 65, 70, 74, 78, 82, 86, 90, 94, 98, 102, 106, 110};
  int cell_vtx[110] = {
    1, 2, 5, 4, 10, 11, 14, 13,
    4, 5, 8, 7, 13, 14, 17, 16,
    2, 3, 5, 11, 12, 14,
    3, 6, 5, 12, 15, 14,
    5, 6, 9, 14, 15, 18,
    5, 9, 8, 14, 18, 17,
    10, 11, 14, 13, 20,
    13, 14, 23, 22, 20,
    10, 13, 22, 19, 20,
    13, 14, 17, 16, 26,
    13, 22, 23, 14, 26,
    13, 16, 25, 22, 26,
    11, 12, 14, 20,
    14, 23, 20, 24,
    14, 15, 24, 12,
    14, 12, 24, 20,
    12, 24, 20, 21,
    14, 15, 18, 24,
    14, 23, 24, 26,
    14, 18, 17, 26,
    14, 18, 26, 24,
    18, 26, 24, 27
  };

  CWP_Mesh_interf_vtx_set(code_name,
                          cpl_name,
                          0,
                          n_vtx,
                          vtx_coord,
                          NULL);

  CWP_Mesh_interf_from_cellvtx_set(code_name,
                                   cpl_name,
                                   0,
                                   n_cell,
                                   cell_vtx_idx,
                                   cell_vtx,
                                   NULL);

  CWP_Mesh_interf_finalize(code_name, cpl_name);


  /* Create and set fields */
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name = "field";

  double *send_val = malloc(sizeof(double) * n_vtx);
  double *recv_val = malloc(sizeof(double) * n_vtx);

  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    send_val[i_vtx] = vtx_coord[3*i_vtx] + vtx_coord[3*i_vtx+1] + vtx_coord[3*i_vtx+2];
  }

  CWP_Field_create(code_name,
                   cpl_name,
                   field_name,
                   CWP_DOUBLE,
                   CWP_FIELD_STORAGE_INTERLACED,
                   1,
                   CWP_DOF_LOCATION_NODE,
                   CWP_FIELD_EXCH_SENDRECV,
                   visu_status);

  CWP_Time_step_beg(code_name,
                    0.0);

  CWP_Field_data_set(code_name,
                     cpl_name,
                     field_name,
                     0,
                     CWP_FIELD_MAP_SOURCE,
                     send_val);

  CWP_Field_data_set(code_name,
                     cpl_name,
                     field_name,
                     0,
                     CWP_FIELD_MAP_TARGET,
                     recv_val);

  /* Compute spatial interpolation weights */
  CWP_Spatial_interp_weights_compute(code_name, cpl_name);

  /* Exchange interpolated field */
  CWP_Field_issend(code_name, cpl_name, field_name);
  CWP_Field_irecv (code_name, cpl_name, field_name);

  CWP_Field_wait_issend(code_name, cpl_name, field_name);
  CWP_Field_wait_irecv (code_name, cpl_name, field_name);

  CWP_Time_step_end(code_name);

  /* Check interpolated field */
  double max_err = 0;
  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    double err = fabs(send_val[i_vtx] - recv_val[i_vtx]);
    if (err > max_err) {
      max_err = err;
    }
  }

  /* Free memory */
  free(send_val);
  free(recv_val);

  /* Finalize CWIPI */
  CWP_Finalize();


  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("End\n");
    fflush(stdout);
  }

  /* Finalize MPI */
  MPI_Finalize();

  if (max_err > 1e-14) {
    printf("Interpolation error is too high : %e\n", max_err);
    fflush(stdout);
    return EXIT_FAILURE;
   }
   else {
    return EXIT_SUCCESS;
   }
}