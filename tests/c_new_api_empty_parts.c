/*
  This file is part of the CWIPI library.

  Copyright (C) 2024 ONERA

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

#include "cwp.h"
#include "cwp_priv.h"


static void
_generate_mesh
(
 MPI_Comm      intra_comm,
 int           n_vtx_seg,
 int          *n_vtx,
 double      **vtx_coord,
 CWP_g_num_t **vtx_ln_to_gn,
 int          *n_elt,
 int         **connec,
 CWP_g_num_t **elt_ln_to_gn
 )
{
  int i_rank;
  MPI_Comm_rank(intra_comm, &i_rank);

  if (i_rank == 0) {
    *n_vtx = n_vtx_seg*n_vtx_seg;
    *n_elt = (n_vtx_seg-1)*(n_vtx_seg-1);

    *vtx_coord    = (double      *) malloc(sizeof(double     ) * (*n_vtx) * 3);
    *vtx_ln_to_gn = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * (*n_vtx));
    *connec       = (int         *) malloc(sizeof(int        ) * (*n_elt) * 4);
    *elt_ln_to_gn = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * (*n_elt));

    double step = 1./(n_vtx_seg - 1);

    int k = 0;
    for (int j = 0; j < n_vtx_seg; j++) {
      for (int i = 0; i < n_vtx_seg; i++) {
        (*vtx_coord)[3*k  ] = i*step;
        (*vtx_coord)[3*k+1] = j*step;
        (*vtx_coord)[3*k+2] = 0;
        (*vtx_ln_to_gn)[k] = k + 1;
        k++;
      }
    }

    k = 0;
    for (int j = 0; j < n_vtx_seg-1; j++) {
      for (int i = 0; i < n_vtx_seg-1; i++) {
        (*connec)[4*k  ] = n_vtx_seg*j     + i     + 1;
        (*connec)[4*k+1] = n_vtx_seg*j     + (i+1) + 1;
        (*connec)[4*k+2] = n_vtx_seg*(j+1) + (i+1) + 1;
        (*connec)[4*k+3] = n_vtx_seg*(j+1) + i     + 1;
        (*elt_ln_to_gn)[k] = k + 1;
        k++;
      }
    }
  }
  else {
    *n_vtx = 0;
    *n_elt = 0;

    *vtx_coord    = NULL;
    *vtx_ln_to_gn = NULL;
    *connec       = NULL;
    *elt_ln_to_gn = NULL;
  }
}


/*----------------------------------------------------------------------
 *
 * Main
 *
 *---------------------------------------------------------------------*/

int
main
(
 int argc, 
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

  int is_code1 = (i_rank%2 == 0);


  const char *code_name;
  const char *coupled_code_name;
  if (is_code1) {
    code_name         = "code1";
    coupled_code_name = "code2";
  }
  else {
    code_name         = "code2";
    coupled_code_name = "code1";
  }


  // Initialize CWIPI
  int n_code         = 1;
  int is_active_rank = 1;
  MPI_Comm intra_comm;
  CWP_Init(comm,
           n_code,
           (const char **) &code_name,
           is_active_rank,
           &intra_comm);


  // Create coupling environment
  int n_part = 1;
  const char *coupling_name = "c_new_api_empty_parts";
  CWP_Cpl_create(code_name,
                 coupling_name,
                 coupled_code_name,
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  CWP_Visu_set(code_name,
               coupling_name,
               1,
               CWP_VISU_FORMAT_ENSIGHT,
               "text");

  // Define interface mesh
  int          n_vtx        = 0;
  double      *vtx_coord    = NULL;
  CWP_g_num_t *vtx_ln_to_gn = NULL;
  int          n_elt        = 0;
  int         *connec       = NULL;
  CWP_g_num_t *elt_ln_to_gn = NULL;

  _generate_mesh(intra_comm,
                 3 + is_code1,
                 &n_vtx,
                 &vtx_coord,
                 &vtx_ln_to_gn,
                 &n_elt,
                 &connec,
                 &elt_ln_to_gn);

  printf("rank %d has %d vtx and %d elt\n", i_rank, n_vtx, n_elt);
  fflush(stdout);

  CWP_Mesh_interf_vtx_set(code_name,
                          coupling_name,
                          0,
                          n_vtx,
                          vtx_coord,
                          vtx_ln_to_gn);


  int block_id = CWP_Mesh_interf_block_add(code_name,
                                           coupling_name,
                                           CWP_BLOCK_FACE_QUAD4);

  CWP_Mesh_interf_block_std_set(code_name,
                                coupling_name,
                                0,
                                block_id,
                                n_elt,
                                connec,
                                elt_ln_to_gn);

  CWP_Mesh_interf_finalize(code_name,
                           coupling_name);

  // Define fields
  const char *field_name = "field";

  double *send_data = malloc(sizeof(double) * n_elt);
  double *recv_data = malloc(sizeof(double) * n_vtx);

  if (is_code1) {
    for (int i = 0; i < n_elt; i++) {
      send_data[i] = elt_ln_to_gn[i];
    }

    CWP_Field_create(code_name,
                     coupling_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     CWP_STATUS_ON);

    CWP_Field_data_set(code_name,
                       coupling_name,
                       field_name,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_data);
  }
  else {
    CWP_Field_create(code_name,
                     coupling_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     CWP_STATUS_ON);

    CWP_Field_data_set(code_name,
                       coupling_name,
                       field_name,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_data);
  }


  CWP_Time_step_beg(code_name,
                    0.0);

  // Compute spatial interpolation weights
  CWP_Spatial_interp_weights_compute(code_name,
                                     coupling_name);


  // Exchange interpolated fields
  if (is_code1) {
    CWP_Field_issend(code_name,
                     coupling_name,
                     field_name);
  }
  else {
    CWP_Field_irecv(code_name,
                    coupling_name,
                    field_name);
  }

  if (is_code1) {
    CWP_Field_wait_issend(code_name,
                          coupling_name,
                          field_name);
  }
  else {
    CWP_Field_wait_irecv(code_name,
                         coupling_name,
                         field_name);
  }

  CWP_Time_step_end(code_name);

  if (i_rank == 0) {
    printf("End :)\n");
    fflush(stdout);
  }

  // Free memory
  if (vtx_coord != NULL) {
    free(vtx_coord);
  }
  if (vtx_ln_to_gn != NULL) {
    free(vtx_ln_to_gn);
  }
  if (connec != NULL) {
    free(connec);
  }
  if (elt_ln_to_gn != NULL) {
    free(elt_ln_to_gn);
  }
  free(send_data);
  free(recv_data);

  // Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();


  return EXIT_SUCCESS;
}
