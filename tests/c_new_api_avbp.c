/*
  This file is part of the CWIPI library.

  Copyright (C) 2017  ONERA

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

#include "cwp.h"

int main(int argc, char *argv[]) {
  int mode;
  if (argc == 1) {
    mode = 0;
  }
  else if (argc == 2) {
    if (strcmp(argv[1], "-0") == 0) mode = 0;
    else if (strcmp(argv[1], "-1") == 0) mode = 1;
    else {
      printf("Unknown argument value\n");
      exit(0);
    }
  }
  else {
    mode = -1;
    printf("Wrong number of arguments\n");
    exit(0);
  }

  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm commcwipi = MPI_COMM_WORLD;
  MPI_Comm_rank(commcwipi, &rank);
  MPI_Comm_size(commcwipi, &comm_world_size);

  const char *coupling_name = "cpl12";

  const char *code_name, *code_coupled_name;
  if (rank == 0) {
    code_name = "avbp1";
    code_coupled_name = "avbp2";
  }
  else if (rank == 1) {
    code_name = "avbp2";
    code_coupled_name = "avbp1";
  }
  else {
    printf("Invalid rank\n");
    fflush(stdout);
    exit(0);
  }

  int nb_codes = 1;
  const char **code_names = (const char **) malloc(sizeof(char *) * nb_codes);
  double *time_init = (double *) malloc(sizeof(double) * nb_codes);
  CWP_Status_t *statuses = (CWP_Status_t *) malloc(sizeof(CWP_Status_t) * nb_codes);
  MPI_Comm *cwp_comms = (MPI_Comm *) malloc(sizeof(MPI_Comm) * nb_codes);

  code_names[0] = code_name;
  time_init[0] = 0.;
  statuses[0] = CWP_STATUS_ON;

  CWP_Init(commcwipi, nb_codes, (const char **) code_names, statuses, time_init, cwp_comms);
  printf("%d --- CWIPI initialised\n", rank);
  fflush(stdout);

  CWP_Cpl_create(code_name, coupling_name, code_coupled_name,
                 CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITH_PART, CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
  printf("%d --- Coupling created\n", rank);
  fflush(stdout);

  CWP_Visu_set(code_name, coupling_name, 1, CWP_VISU_FORMAT_ENSIGHT, "text");

  int nnode_cpl_s, ncell_cpl_s, nnode_cpl_d;
  int n_vtx_per_cell = 8;

  char fn[100];
  sprintf(fn, "avbp/periodic/avbp_dimensions_%d.out", rank);
  FILE *fil = fopen(fn, "r");
  fscanf(fil, "%d", &nnode_cpl_s);
  fscanf(fil, "%d", &ncell_cpl_s);
  fscanf(fil, "%d", &nnode_cpl_d);
  fclose(fil);

  // Limit the mesh size
//  ncell_cpl_s = 100;

  int n_var, stride;
  if (mode == 0) {
    n_var = 6;
    stride = 4;
  }
  else {
    n_var = 1;
    stride = 1;
  }

  double *xyz_s = (double *) malloc(sizeof(double) * (3 * nnode_cpl_s));
  int *connind_s = (int *) malloc(sizeof(int) * (ncell_cpl_s + 1));
  int *elems_s = (int *) malloc(sizeof(int) * (n_vtx_per_cell * ncell_cpl_s));
  double *xyz_dest = (double *) malloc(sizeof(double) * (3 * nnode_cpl_d));

  double *sfields = (double *) malloc(sizeof(double) * (n_var * nnode_cpl_s));
  double *rfields = (double *) malloc(sizeof(double) * (n_var * nnode_cpl_d));

  // xyz_s
  sprintf(fn, "avbp/periodic/avbp_xyz_s_%d.out", rank);
  fil = fopen(fn, "r");
  for (int i = 0 ; i < 3 * nnode_cpl_s ; i++) {
    fscanf(fil, "%lf", &xyz_s[i]);
  }
  fclose(fil);
  // connind_s
  sprintf(fn, "avbp/periodic/avbp_connind_s_%d.out", rank);
  fil = fopen(fn, "r");
  for (int i = 0 ; i < 3 * ncell_cpl_s + 1 ; i++) {
    fscanf(fil, "%d", &connind_s[i]);
  }
  fclose(fil);
  // elems_s
  sprintf(fn, "avbp/periodic/avbp_elems_s_%d.out", rank);
  fil = fopen(fn, "r");
  for (int i = 0 ; i < 3 * n_vtx_per_cell * ncell_cpl_s ; i++) {
    fscanf(fil, "%d", &elems_s[i]);
  }
  fclose(fil);
  // xyz_dest
  sprintf(fn, "avbp/periodic/avbp_xyz_dest_%d.out", rank);
  fil = fopen(fn, "r");
  for (int i = 0 ; i < 3 * nnode_cpl_d ; i++) {
    fscanf(fil, "%lf", &xyz_dest[i]);
  }
  fclose(fil);

  // Fields
  for (int i = 0 ; i < n_var * nnode_cpl_d ; i++) {
    rfields[i] = 0.;
  }
  if (mode == 0) {
    // avbp_sfields
    sprintf(fn, "avbp/periodic/avbp_sfields_%d.out", rank);
    fil = fopen(fn, "r");
    for (int i = 0 ; i < n_var * nnode_cpl_s ; i++) {
      fscanf(fil, "%lf", &sfields[i]);
    }
    fclose(fil);
  }
  else {
    // Test valeurs non uniformes entrelacees
    for (int i = 0 ; i < n_var * nnode_cpl_s ; i++) {
      sfields[i] = 2. * (i + 1);
    }
  }
  printf("%d --- Data read\n", rank);
  fflush(stdout);

  CWP_Mesh_interf_vtx_set(code_name, coupling_name, 0, nnode_cpl_s, xyz_s, NULL);

//  int block_id = CWP_Mesh_interf_block_add(code_name, coupling_name, CWP_BLOCK_FACE_QUAD4);
//  CWP_Mesh_interf_block_std_set(code_name, coupling_name, 0, block_id, ncell_cpl_s, elems_s, NULL);

  int n_tetra = 0, n_prism = 0;
  int n_pts_per_elt;
  for (int i = 0 ; i < ncell_cpl_s ; ++i) {
    n_pts_per_elt = connind_s[i + 1] - connind_s[i];
    switch (n_pts_per_elt) {
      case 4:
        n_tetra++;
        break;
      case 6:
        n_prism++;
        break;
      default:
        printf("Invalid n_pts_per_elt (count) = %d\n", n_pts_per_elt);
        fflush(stdout);
    }
  }

  printf("\tn_tetra: %d\n", n_tetra);
  printf("\tn_prism: %d\n", n_prism);
  printf("\tn_elts: %d\n", n_tetra + n_prism);
  fflush(stdout);

  int *elems_tetra = (int *) malloc(4 * n_tetra * sizeof(int));
  int *elems_prism = (int *) malloc(6 * n_prism * sizeof(int));

  int tetra_order_avbp[4] = {0, 1, 2, 3};
  int tetra_order_cwp[4] = {0, 1, 3, 2};
  int tetra_order_grid[4] = {0, 2, 1, 3};

  int prism_order_avbp[6] = {0, 1, 2, 3, 4, 5};
  int prism_order_cwp[6] = {0, 3, 4, 1, 5, 2};
  int prism_order_grid[6] = {5, 0, 3, 4, 1, 2};

  int *tetra_order1 = tetra_order_avbp;
  int *prism_order1 = prism_order_avbp;

  int *tetra_order2 = tetra_order_avbp;
  int *prism_order2 = prism_order_avbp;

  int tetra_counter = 0, prism_counter = 0, elt_counter = 0;
  for (int i = 0 ; i < ncell_cpl_s ; ++i) {
    n_pts_per_elt = connind_s[i + 1] - connind_s[i];
    switch (n_pts_per_elt) {
      case 4:
        for (int j = 0 ; j < n_pts_per_elt ; ++j) {
          elems_tetra[tetra_counter + tetra_order1[j]] = elems_s[elt_counter + tetra_order2[j]];
        }
        tetra_counter += n_pts_per_elt;
        break;
      case 6:
        for (int j = 0 ; j < n_pts_per_elt ; ++j) {
          elems_prism[prism_counter + prism_order1[j]] = elems_s[elt_counter + prism_order2[j]];
        }
        prism_counter += n_pts_per_elt;
        break;
      default:
        printf("Invalid n_pts_per_elt (fill) = %d\n", n_pts_per_elt);
        fflush(stdout);
    }
    elt_counter += n_pts_per_elt;
  }

  int block_tetra = CWP_Mesh_interf_block_add(code_name, coupling_name, CWP_BLOCK_CELL_TETRA4);
  CWP_Mesh_interf_block_std_set(code_name, coupling_name, 0, block_tetra, n_tetra, elems_tetra, NULL);
//  int block_prism = CWP_Mesh_interf_block_add(code_name, coupling_name, CWP_BLOCK_CELL_PRISM6);
//  CWP_Mesh_interf_block_std_set(code_name, coupling_name, 0, block_prism, n_prism, elems_prism, NULL);

  CWP_Mesh_interf_finalize(code_name, coupling_name);
  printf("%d --- Geometry set\n", rank);

  const char *send_field_name, *recv_field_name;
  if (rank == 0) {
    send_field_name = "code1_to_code2";
    recv_field_name = "code2_to_code1";
  }
  else {
    send_field_name = "code2_to_code1";
    recv_field_name = "code1_to_code2";
  }

  CWP_Field_create(code_name, coupling_name, send_field_name, CWP_DOUBLE, CWP_FIELD_STORAGE_BLOCK, stride,
                   CWP_DOF_LOCATION_NODE, CWP_FIELD_EXCH_SEND, CWP_STATUS_OFF);
  CWP_Field_data_set(code_name, coupling_name, send_field_name, 0, CWP_FIELD_MAP_SOURCE, sfields);

  CWP_Field_create(code_name, coupling_name, recv_field_name, CWP_DOUBLE, CWP_FIELD_STORAGE_BLOCK, stride,
                   CWP_DOF_LOCATION_USER, CWP_FIELD_EXCH_RECV, CWP_STATUS_OFF);
  CWP_Field_data_set(code_name, coupling_name, recv_field_name, 0, CWP_FIELD_MAP_TARGET, rfields);

  printf("%d --- Fields created\n", rank);

  CWP_User_tgt_pts_set(code_name, coupling_name, 0, nnode_cpl_d, xyz_dest, NULL);
  printf("%d --- Points to locate set\n", rank);

  CWP_Spatial_interp_weights_compute(code_name, coupling_name);
  printf("%d --- Localisation done\n", rank);

  int n_not_located_points;
  int n_located_points;
  int n_involved_srcs;
//  const int *located_points;
  const int *not_located_points;
//  const int *involved_srcs;

  n_not_located_points = CWP_N_uncomputed_tgts_get(code_name, coupling_name, recv_field_name, 0);
  n_located_points = CWP_N_computed_tgts_get(code_name, coupling_name, recv_field_name, 0);
  n_involved_srcs = CWP_N_involved_srcs_get(code_name, coupling_name, send_field_name, 0);

//  located_points = CWP_Computed_tgts_get(code_name, coupling_name, recv_field_name, 0);
  not_located_points = CWP_Uncomputed_tgts_get(code_name, coupling_name, recv_field_name, 0);
//  involved_srcs = CWP_Involved_srcs_get(code_name, coupling_name, send_field_name, 0);

  printf("%d --- n_not_located_points %d\n", rank, n_not_located_points);
  printf("%d --- n_located_points %d\n", rank, n_located_points);
  printf("%d --- n_distant_located_points %d\n", rank, n_involved_srcs);
  printf("not located points :");
  for (int i = 0 ; i < n_not_located_points ; ++i) printf("%d ", not_located_points[i]);
//  printf("\nlocated points: ");
//  for (int i = 0 ; i < n_located_points ; ++i) printf("%d ", located_points[i]);
//  printf("\ninvolved source points");
//  for (int i = 0 ; i < n_involved_srcs ; ++i) printf("%d ", involved_srcs[i]);
  printf("\n");

  // Coordinates of not located points
  if (n_not_located_points != 0) {
    printf("%d --- not located points coordinates:\n", rank);
    for (int i = 0 ; i < n_not_located_points ; ++i) {
      printf("\t%f %f %f\n", xyz_dest[3 * not_located_points[i]], xyz_dest[3 * not_located_points[i] + 1], xyz_dest[3 * not_located_points[i] + 2]);
    }
    printf("\n");
  }

  CWP_Field_irecv(code_name, coupling_name, recv_field_name);
  CWP_Field_issend(code_name, coupling_name, send_field_name);
  CWP_Field_wait_irecv(code_name, coupling_name, recv_field_name);
  CWP_Field_wait_issend(code_name, coupling_name, send_field_name);

  printf("%d --- Exchanges completed\n", rank);

  CWP_Mesh_interf_del(code_name, coupling_name);
  printf("%d --- Interface deleted\n", rank);

  CWP_Cpl_del(code_name, coupling_name);
  printf("%d --- Coupling deleted\n", rank);

  CWP_Finalize();
  printf("%d --- CWIPI finalized\n", rank);

  free(code_names);
  free(time_init);
  free(statuses);
  free(cwp_comms);
  free(xyz_s);
  free(connind_s);
  free(elems_s);
  free(xyz_dest);
  free(sfields);
  free(rfields);
}
