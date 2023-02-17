/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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

#include "cwp.h"
#include "cwp_priv.h"

static void
_userInterpolation
(
 const char           *local_code_name,
 const char           *cpl_id,
 const char           *field_id,
 CWP_Spatial_interp_t  spartial_interp_algorithm,
 double               *buffer_in,
 double               *buffer_out
)
{
  int **elt_pts_inside_idx = CWP_Interp_location_elt_pts_inside_idx_get(local_code_name,
                                                                        cpl_id,
                                                                        field_id,
                                                                        spartial_interp_algorithm);

  CWP_g_num_t **points_gnum = CWP_Interp_location_points_gnum_get(local_code_name,
                                                                  cpl_id,
                                                                  field_id,
                                                                  spartial_interp_algorithm);

  double **points_coords = CWP_Interp_location_points_coords_get(local_code_name,
                                                                  cpl_id,
                                                                  field_id,
                                                                  spartial_interp_algorithm);

  double **points_uvw = CWP_Interp_location_points_uvw_get(local_code_name,
                                                           cpl_id,
                                                           field_id,
                                                           spartial_interp_algorithm);

  double **points_dist2 = CWP_Interp_location_points_dist2_get(local_code_name,
                                                               cpl_id,
                                                               field_id,
                                                               spartial_interp_algorithm);

  double **points_projected_coords = CWP_Interp_location_points_projected_coords_get(local_code_name,
                                                                                     cpl_id,
                                                                                     field_id,
                                                                                     spartial_interp_algorithm);

  int **cell_vtx_idx = CWP_Interp_location_cell_vtx_idx_get(local_code_name,
                                                            cpl_id,
                                                            field_id,
                                                            spartial_interp_algorithm);

  int **cell_vtx = CWP_Interp_location_cell_vtx_get(local_code_name,
                                                    cpl_id,
                                                    field_id,
                                                    spartial_interp_algorithm);

  buffer_out[0] = buffer_in[0];
  buffer_out[1] = buffer_in[1];
}

/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // init cwipi
  int n_part = 1;
  int n_code = 1;
  const char  **code_name         = malloc(sizeof(char *) * n_code);
  const char  **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t *is_active_rank    = malloc(sizeof(CWP_Status_t) * n_code);
  double       *time_init         = malloc(sizeof(double) * n_code);

  int code_id;
  if (i_rank % 2 == 0) {
    code_id = 1;
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
  }
  else {
    code_id = 2;
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
  }

  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           time_init,
           intra_comm);

  // create coupling
  const char *coupling_name = "coupling";
  CWP_Spatial_interp_t loc_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 loc_method,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  // set mesh
  // int block_id = CWP_Mesh_interf_block_add(code_name[0],
  //                                          coupling_name,
  //                                          CWP_BLOCK_FACE_TRIA3);
  int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                           coupling_name,
                                           CWP_BLOCK_FACE_POLY);

  int n_vtx = 4;
  double *vtx_coord = malloc(sizeof(double) * 3 * n_vtx);
  vtx_coord[0]  = 0.0;
  vtx_coord[1]  = 0.0;
  vtx_coord[2]  = 0.0;
  vtx_coord[3]  = 1.0;
  vtx_coord[4]  = 0.0;
  vtx_coord[5]  = 0.0;
  vtx_coord[6]  = 0.0;
  vtx_coord[7]  = 1.0;
  vtx_coord[8]  = 0.0;
  vtx_coord[9]  = 1.0;
  vtx_coord[10] = 1.0;
  vtx_coord[11] = 0.0;
  CWP_Mesh_interf_vtx_set(code_name[0],
                          coupling_name,
                          0,
                          n_vtx,
                          vtx_coord,
                          NULL);

  int n_face = 2;
  int *face_vtx_idx = malloc(sizeof(int) * (n_face + 1));
  face_vtx_idx[0] = 0;
  face_vtx_idx[1] = 3;
  face_vtx_idx[2] = 6;
  int *face_vtx     = malloc(sizeof(int) * 3 * n_face);
  face_vtx[0] = 1;
  face_vtx[1] = 2;
  face_vtx[2] = 3;
  face_vtx[3] = 2;
  face_vtx[4] = 4;
  face_vtx[5] = 3;
  CWP_g_num_t *gnum_elt = malloc(sizeof(CWP_g_num_t) * n_face);
  gnum_elt[0] = 1;
  gnum_elt[1] = 2;
  // CWP_Mesh_interf_block_std_set(code_name[0],
  //                                  coupling_name,
  //                                  0,
  //                                  block_id,
  //                                  n_face,
  //                                  face_vtx,
  //                                  gnum_elt);
  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   coupling_name,
                                   0,
                                   block_id,
                                   n_face,
                                   face_vtx_idx,
                                   face_vtx,
                                   gnum_elt);

  CWP_Mesh_interf_finalize(code_name[0], coupling_name);

  MPI_Barrier(MPI_COMM_WORLD);

  // create field
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name = "Maine Coon";

  if (code_id == 1) {
    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);
  } else {
    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);
  }

  double *send_buff = malloc(sizeof(double) * n_face);
  double *recv_buff = malloc(sizeof(double) * n_face);
  if (code_id == 1) {
    send_buff[0] = gnum_elt[0];
    send_buff[1] = gnum_elt[1];
    CWP_Field_data_set(code_name[0],
                       coupling_name,
                       field_name,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_buff);
  } else {
    CWP_Field_data_set(code_name[0],
                       coupling_name,
                       field_name,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_buff);
  }

  CWP_Interp_function_set(code_name[0], coupling_name, field_name, _userInterpolation);

  MPI_Barrier(MPI_COMM_WORLD);

  CWP_Spatial_interp_weights_compute(code_name[0], coupling_name);

  MPI_Barrier(MPI_COMM_WORLD);

  // exchange field
  if (code_id == 1) {
    CWP_Field_issend(code_name[0], coupling_name, field_name);
  } else {
    CWP_Field_irecv(code_name[0], coupling_name, field_name);
  }

  if (code_id == 1) {
    CWP_Field_wait_issend(code_name[0], coupling_name, field_name);
    printf("send_buff[0] : %f / send_buff[1] : %f\n", send_buff[0], send_buff[1]);
    fflush(stdout);
  } else {
    CWP_Field_wait_irecv(code_name[0], coupling_name, field_name);
    printf("recv_buff[0] : %f / recv_buff[1] : %f\n", recv_buff[0], recv_buff[1]);
    fflush(stdout);
  }

  // del
  CWP_Mesh_interf_del(code_name[0], coupling_name);
  CWP_Cpl_del(code_name[0], coupling_name);

  // finalize cwipi
  CWP_Finalize();

  // free
  free(code_name);
  free(coupled_code_name);
  free(is_active_rank);
  free(time_init);
  free(vtx_coord);
  free(face_vtx_idx);
  free(face_vtx);
  free(send_buff);
  free(recv_buff);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
