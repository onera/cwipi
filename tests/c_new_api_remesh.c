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
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#include "cwp.h"
#include "cwp_priv.h"

#include "cwp_priv.h"

static
void
_eval_field1
(
  const int     n,
  const double *coord,
  const double  time,
        double *field_value
)
{
  for (int i = 0; i < n; i++) {
    field_value[i] = cos( sqrt(coord[3*i]*coord[3*i] + coord[3*i+1]*coord[3*i+1]) - 0.1*time );
  }
}


static
void
_eval_field2
(
  const int     n_elt,
  const int    *elt_vtx_idx,
  const int    *elt_vtx,
  const double *coord,
  const double  time,
        double *field_value
)
{
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    field_value[i_elt] = 0.;

    for (int i = elt_vtx_idx[i_elt]; i < elt_vtx_idx[i_elt+1]; i++) {
      double x = coord[3*(elt_vtx[i] - 1)];
      field_value[i_elt] += cos(x + 0.1*time);
    }

    field_value[i_elt] /= elt_vtx_idx[i_elt+1] - elt_vtx_idx[i_elt];
  }
}



int
main
(
 int   argc,
 char *argv[]
)
{
  // Initialize MPI
  int i_rank;
  int n_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);


  int i_code = i_rank%2;

  // Initialize CWIPI
  const char   *code_name         = NULL;
  const char   *coupled_code_name = NULL;
  CWP_Status_t  is_active_rank    = CWP_STATUS_ON;
  int           n_part            = 1;

  if (i_code == 0) {
    code_name         = "code1";
    coupled_code_name = "code2";
  }
  else {
    code_name         = "code2";
    coupled_code_name = "code1";
  }


  MPI_Comm intra_comm;

  CWP_Init(comm,
           1,
           &code_name,
           is_active_rank,
           &intra_comm);

  int i_rank_intra;
  MPI_Comm_rank(intra_comm, &i_rank_intra);


  // Create coupling
  const char *cpl_name = "c_new_api_remesh";
  CWP_Cpl_create(code_name,
                 cpl_name,
                 coupled_code_name,
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 n_part,
                 CWP_DYNAMIC_MESH_VARIABLE,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  CWP_Visu_set(code_name, cpl_name, 1, CWP_VISU_FORMAT_ENSIGHT, "text");


  // Initial mesh
  int          n_vtx       = 0;
  int          n_elt       = 0;
  double      *coords      = NULL;
  int         *elt_vtx_idx = NULL;
  int         *elt_vtx     = NULL;
  CWP_g_num_t  n_vtx_seg   = 10 - 3*i_code;
  CWPT_generate_mesh_rectangle_simplified(intra_comm,
                                          n_vtx_seg,
                                          &n_vtx,
                                          &n_elt,
                                          &coords,
                                          &elt_vtx_idx,
                                          &elt_vtx);

  CWP_Mesh_interf_vtx_set(code_name,
                          cpl_name,
                          0,
                          n_vtx,
                          coords,
                          NULL);

  int block_id = CWP_Mesh_interf_block_add(code_name,
                                           cpl_name,
                                           CWP_BLOCK_FACE_TRIA3);

  CWP_Mesh_interf_block_std_set(code_name,
                                cpl_name,
                                0,
                                block_id,
                                n_elt,
                                elt_vtx,
                                NULL);

  CWP_Mesh_interf_finalize(code_name, cpl_name);



  // Create fields
  const char *field_name1 = "field1";
  const char *field_name2 = "field2";

  double *send_data1 = NULL;
  double *send_data2 = NULL;
  double *recv_data1 = NULL;
  double *recv_data2 = NULL;

  if (i_code == 0) {
    CWP_Field_create(code_name,
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SEND,
                     CWP_STATUS_ON);

    CWP_Field_create(code_name,
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     CWP_STATUS_ON);

    send_data1 = malloc(sizeof(double) * n_vtx);
    recv_data2 = malloc(sizeof(double) * n_vtx);

  }
  else {
    CWP_Field_create(code_name,
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_RECV,
                     CWP_STATUS_ON);

    CWP_Field_create(code_name,
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     CWP_STATUS_ON);

    recv_data1 = malloc(sizeof(double) * n_elt);
    send_data2 = malloc(sizeof(double) * n_elt);
  }



  // Time loop
  int n_time_steps = 30;
  double time = 0.0;


  int freq_remesh = 5 + 3*i_code;

  for (int i_step = 0; i_step < n_time_steps; i_step++) {

    if (i_rank == 0) {
      printf("Step %d\n", i_step);
      fflush(stdout);
    }

    // Begind time step
    CWP_Time_step_beg(code_name,
                      time);
    time += 1.;

    // Remesh if required
    if ((i_step+1)%freq_remesh == 0) {

      if (i_rank_intra == 0) {
        printf("  Remesh %s\n", code_name);
        fflush(stdout);
      }

      CWP_Mesh_interf_del(code_name, cpl_name);

      free(coords     );
      free(elt_vtx_idx);
      free(elt_vtx    );

      n_vtx_seg += 2 + i_code;
      CWPT_generate_mesh_rectangle_simplified(intra_comm,
                                              n_vtx_seg,
                                              &n_vtx,
                                              &n_elt,
                                              &coords,
                                              &elt_vtx_idx,
                                              &elt_vtx);

      CWP_Mesh_interf_vtx_set(code_name,
                              cpl_name,
                              0,
                              n_vtx,
                              coords,
                              NULL);

      block_id = CWP_Mesh_interf_block_add(code_name,
                                           cpl_name,
                                           CWP_BLOCK_FACE_TRIA3);

      CWP_Mesh_interf_block_std_set(code_name,
                                    cpl_name,
                                    0,
                                    block_id,
                                    n_elt,
                                    elt_vtx,
                                    NULL);

      CWP_Mesh_interf_finalize(code_name, cpl_name);

      if (i_code == 0) {
        send_data1 = realloc(send_data1, sizeof(double) * n_vtx);
        recv_data2 = realloc(recv_data2, sizeof(double) * n_vtx);
      }
      else {
        recv_data1 = realloc(recv_data1, sizeof(double) * n_elt);
        send_data2 = realloc(send_data2, sizeof(double) * n_elt);
      }
    }

    // Set field value pointers
    if (i_code == 0) {
      _eval_field1(n_vtx, coords, time, send_data1);

      CWP_Field_data_set(code_name,
                         cpl_name,
                         field_name1,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_data1);

      CWP_Field_data_set(code_name,
                         cpl_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_data2);
    }
    else {
      _eval_field2(n_elt, elt_vtx_idx, elt_vtx, coords, time, send_data2);

      CWP_Field_data_set(code_name,
                         cpl_name,
                         field_name1,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_data1);

      CWP_Field_data_set(code_name,
                         cpl_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_data2);
    }


    // Compute spatial interpolation weights
    CWP_Spatial_interp_weights_compute(code_name, cpl_name);


    // Exchange fields
    if (i_code == 0) {
      CWP_Field_issend(code_name, cpl_name, field_name1);
      CWP_Field_irecv (code_name, cpl_name, field_name2);
    }
    else {
      CWP_Field_irecv (code_name, cpl_name, field_name1);
      CWP_Field_issend(code_name, cpl_name, field_name2);
    }

    // Overlap exchanges with computations if possible...

    if (i_code == 0) {
      CWP_Field_wait_issend(code_name, cpl_name, field_name1);
      CWP_Field_wait_irecv (code_name, cpl_name, field_name2);
    }
    else {
      CWP_Field_wait_irecv (code_name, cpl_name, field_name1);
      CWP_Field_wait_issend(code_name, cpl_name, field_name2);
    }

    // End time step
    CWP_Time_step_end(code_name);
  }


  // Delete mesh
  CWP_Mesh_interf_del(code_name, cpl_name);

  // Delete coupling
  CWP_Cpl_del(code_name, cpl_name);

  // Finalize
  CWP_Finalize();

  // Free memory
  free(coords     );
  free(elt_vtx_idx);
  free(elt_vtx    );
  if (i_code == 0) {
    free(send_data1);
    free(recv_data2);
  }
  else {
    free(recv_data1);
    free(send_data2);
  }

  if (i_rank == 0) {
    printf("The End :)\n");
    fflush(stdout);
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}
