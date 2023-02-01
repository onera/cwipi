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
#include <unistd.h>

#include "cwp.h"
#include "cwp_priv.h"
#include "grid_mesh.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_io.h"
#include "pdm_distrib.h"
#include "pdm_printf.h"
#include "pdm_partitioning_algorithm.h"
#include "client_server/client.h"

#include "cwp_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CWP_HEADER_SIZE    32

/*=============================================================================
 * Util functions
 *============================================================================*/

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -c     Filename of the server configuration file.\n\n"
     "  -h     This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
 int            argc,
 char         **argv,
 char         **config  // filename for server ip adresses + ports
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *config = argv[i];
      }
    }

    else
      _usage(EXIT_FAILURE);
    i++;

  }

}

/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main
(
 int argc,
 char *argv[]
)
{
  // default
  char *config     = NULL;

  _read_args(argc,
             argv,
             &config);

  if (config == NULL) {
    config = (char *) "../bin/cwp_config_srv.txt";
  }

  // mpi
  int rank;
  int commWorldSize;

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &commWorldSize);

  // launch server
  char launch_server[99];
  sprintf(launch_server, "mpirun -n %d ../bin/server_main &", commWorldSize);
  system(launch_server);

  while (access(config, R_OK) != 0) {
    printf("HERE\n");
    // wait
  }
  sleep(5);

  // Read args from command line
  int nVertexSeg = 10;
  double randLevel = 0.4;

  // Initialization
  int n_code = 1;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t *is_active_rank = malloc(sizeof(CWP_Status_t) * n_code);
  double *time_init = malloc(sizeof(double) * n_code);

  is_active_rank[0] = CWP_STATUS_ON;
  time_init[0] = 0.;
  if (rank % 2 == 0) {
    printf("%d - Working for code1\n", rank);
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
  }

  else {
    printf("%d - Working for code2\n", rank);
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
  }

  CWP_client_Init(comm,
                  config,
                  n_code,
                  (const char **) code_name,
                  is_active_rank,
                  time_init);

  printf("%d - Create coupling\n", rank);
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1";
  int nb_part = 1;
  CWP_client_Cpl_create(code_name[0],                                          // Code name
                 cpl_name,                                              // Coupling id
                 coupled_code_name[0],                                  // Coupled application id
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,                                // Coupling type
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                 nb_part,                                               // Number of partitions
                 CWP_DYNAMIC_MESH_STATIC,                               // Mesh type
                 CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency

  printf("%d - Set visu\n", rank);
  CWP_client_Visu_set(code_name[0],            // Code name
               cpl_name,                // Coupling id
               1,                       // Postprocessing frequency
               CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
               "text");                 // Postprocessing option

  // Mesh definition
  printf("%d - Create mesh\n", rank);

  int nVertex;
  double **coords = NULL;
  int nElts;
  int **eltsConnecPointer = NULL;
  int **eltsConnec = NULL;

  // Domain bounds
  const double xmin = -10;
  const double xmax = 10;
  const double ymin = -10;
  const double ymax = 10;

  nVertex = nVertexSeg * nVertexSeg;
  nElts = (nVertexSeg - 1) * (nVertexSeg - 1);

  coords            = (double **) malloc(sizeof(double *) * n_code);
  eltsConnecPointer = (int    **) malloc(sizeof(int    *) * n_code);
  eltsConnec        = (int    **) malloc(sizeof(int    *) * n_code);

  coords[0]            = (double *) malloc(sizeof(double) * 3 * nVertex);
  eltsConnecPointer[0] = (int *) malloc(sizeof(int) * (nElts + 1));
  eltsConnec[0]        = (int *) malloc(sizeof(int) * 4 * nElts);
  randLevel = 0.4;

  // Create local comm to code using comm_split

  PDM_MPI_Comm LocalComm;
  int color;
  if (rank % 2 == 0) {
    color = 1;
  } else {
    color = 0;
  }
  PDM_MPI_Comm_split(PDM_MPI_mpi_2_pdm_mpi_comm(&comm), color, 0, &LocalComm);
  int size = -1;
  PDM_MPI_Comm_size(LocalComm, &size);

  srand(time(NULL));

  if (strcmp(code_name[0], "code1") == 0) {
    grid_mesh(xmin,
              xmax,
              ymin,
              ymax,
              randLevel,
              nVertexSeg,
              (int) sqrt(size),
              coords[0],
              eltsConnecPointer[0],
              eltsConnec[0],
              * (MPI_Comm *) PDM_MPI_2_mpi_comm(LocalComm));
  }

  randLevel = 0.2;
  if (strcmp(code_name[0], "code2") == 0) {
    grid_mesh(xmin,
              xmax,
              ymin,
              ymax,
              randLevel,
              nVertexSeg,
              (int) sqrt(size),
              coords[0],
              eltsConnecPointer[0],
              eltsConnec[0],
              * (MPI_Comm *) PDM_MPI_2_mpi_comm(LocalComm));
  }

  printf("%d - Number of vertex   : %d\n", rank, nVertex);
  printf("%d - Number of elements : %d\n", rank, nElts);

  CWP_g_num_t *global_num_vtx = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * nVertex);
  for (int i = 0; i < nVertex; i++) {
    global_num_vtx[i] = i + 1;
  }

  CWP_client_Mesh_interf_vtx_set(code_name[0], cpl_name, 0, nVertex, coords[0], global_num_vtx);

  int block_id = CWP_client_Mesh_interf_block_add(code_name[0], cpl_name, CWP_BLOCK_FACE_POLY);

  CWP_g_num_t *global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * nElts);
  for (int i = 0; i < nElts; i++) {
    global_num[i] = i + 1;
  }

  int GETnElts = -1;
  int *GETeltsConnecPointer  = NULL;
  int *GETeltsConnec         = NULL;
  CWP_g_num_t *GETglobal_num = NULL;

  CWP_client_Mesh_interf_f_poly_block_set(code_name[0],
                                          cpl_name,
                                          0,
                                          block_id,
                                          nElts,
                                          eltsConnecPointer[0],
                                          eltsConnec[0],
                                          global_num);

  CWP_client_Mesh_interf_f_poly_block_get(code_name[0],
                                          cpl_name,
                                          0,
                                          block_id,
                                          &GETnElts,
                                          &GETeltsConnecPointer,
                                          &GETeltsConnec,
                                          &GETglobal_num);

  printf("GETnElts == nElts: %d\n", GETnElts == nElts);
  int equal = -1;
  for (int i = 0; i < (nElts + 1); i++) {
    equal = (GETeltsConnecPointer[i] == eltsConnecPointer[0][i]);
    if (equal == 0) {
      break;
    }
  }
  printf("eltsConnecPointer equal: %d\n", equal);
  equal = -1;
  for (int i = 0; i < 4 * nElts; i++) {
    equal = (GETeltsConnec[i] == eltsConnec[0][i]);
    if (equal == 0) {
      break;
    }
  }
  printf("eltsConnec equal: %d\n", equal);
  equal = -1;
  for (int i = 0; i < nElts; i++) {
    equal = (global_num[i] == GETglobal_num[i]);
    if (equal == 0) {
      break;
    }
  }
  printf("global_num equal: %d\n", equal);

  CWP_client_Mesh_interf_finalize(code_name[0], cpl_name);

  printf("%d - Exchange Code1 <-> Code2\n", rank);

  double **sendValues = (double **) malloc(sizeof(double *) * n_code);
  double **recvValues = (double **) malloc(sizeof(double *) * n_code);

  if (strcmp(code_name[0], "code1") == 0) {
    sendValues[0] = (double *) malloc(sizeof(double) * nVertex);
    recvValues[0] = (double *) malloc(sizeof(double) * nElts);
    for (int i = 0 ; i < nVertex ; i++) {
      sendValues[0][i] = coords[0][3 * i];
    }
  }
  else {
    sendValues[0] = (double *) malloc(sizeof(double) * nElts);
    recvValues[0] = (double *) malloc(sizeof(double) * nVertex);
    for (int i = 0 ; i < nElts ; i++) {
      sendValues[0][i] = rank;
    }
  }

  // Exchange
  // field1: code1 -> code2
  // field2: code2 -> code1
  const char *field_name1 = "cooX";
  const char *field_name2 = "rank";

  CWP_Status_t visu_status = CWP_STATUS_ON;
  printf("%d - Defining fields\n", rank);
  if (strcmp(code_name[0], "code1") == 0) {
    CWP_client_Field_create(code_name[0],
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLEAVED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);
    CWP_client_Field_create(code_name[0],
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLEAVED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);

    CWP_client_Field_data_set(code_name[0], cpl_name, field_name1, 0, CWP_FIELD_MAP_SOURCE, nVertex, sendValues[0]);
    CWP_client_Field_data_set(code_name[0], cpl_name, field_name2, 0, CWP_FIELD_MAP_TARGET, nElts, recvValues[0]);
  }
  
  else if (strcmp(code_name[0], "code2") == 0) {
    CWP_client_Field_create(code_name[0],
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLEAVED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);
    CWP_client_Field_create(code_name[0],
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLEAVED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_client_Field_data_set(code_name[0], cpl_name, field_name2, 0, CWP_FIELD_MAP_SOURCE, nElts, sendValues[0]);
    CWP_client_Field_data_set(code_name[0], cpl_name, field_name1, 0, CWP_FIELD_MAP_TARGET, nVertex, recvValues[0]);
  }

  printf("%d - Before compute\n", rank);
  // CWP_next_recv_time_set(code_name[0],
  //                        cpl_name,
  //                        0.);
  CWP_client_Spatial_interp_weights_compute(code_name[0], cpl_name);

  int n_uncomputed = 0;

  if (strcmp(code_name[0], "code1") == 0) {
    n_uncomputed = CWP_client_N_uncomputed_tgts_get (code_name[0], cpl_name, field_name2, 0);
  }

  else if (strcmp(code_name[0], "code2") == 0) {
    n_uncomputed = CWP_client_N_uncomputed_tgts_get (code_name[0], cpl_name, field_name1, 0);
  }

  printf("%d - After compute %d\n", rank, n_uncomputed);

  if (strcmp(code_name[0], "code1") == 0) {
    CWP_client_Field_issend(code_name[0], cpl_name, field_name1);
    CWP_client_Field_irecv(code_name[0], cpl_name, field_name2);
  }
  else if (strcmp(code_name[0], "code2") == 0) {
    CWP_client_Field_irecv(code_name[0], cpl_name, field_name1);
    CWP_client_Field_issend(code_name[0], cpl_name, field_name2);
  }

  free(recvValues[0]);
  recvValues[0] = NULL;

  if (strcmp(code_name[0], "code1") == 0) {
    CWP_client_Field_wait_issend(code_name[0], cpl_name, field_name1);
    CWP_client_Field_wait_irecv(code_name[0], cpl_name, field_name2, &recvValues[0]);
  }
  else if (strcmp(code_name[0], "code2") == 0) {
    CWP_client_Field_wait_irecv(code_name[0], cpl_name, field_name1, &recvValues[0]);
    CWP_client_Field_wait_issend(code_name[0], cpl_name, field_name2);
  }

  // Mesh_interf_from_faceedge_set
  int          **part_pface_edge_idx = NULL;
  int          **part_pface_edge     = NULL;
  int           *part_pn_edge        = NULL;
  int          **part_pedge_vtx      = NULL;
  PDM_g_num_t  **part_pedge_ln_to_gn = NULL;

  PDM_g_num_t   *pface_ln_to_gn = malloc(sizeof(PDM_g_num_t) * nElts);
  PDM_g_num_t   *pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t) * nVertex);

  PDM_g_num_t *distrib_elt = PDM_compute_entity_distribution(PDM_MPI_mpi_2_pdm_mpi_comm(&comm), nElts);
  PDM_g_num_t *distrib_vtx = PDM_compute_entity_distribution(PDM_MPI_mpi_2_pdm_mpi_comm(&comm), nVertex);

  for (int i = 0; i < nElts; i++) {
    pface_ln_to_gn[i] = distrib_elt[rank] + i + 1;
  }

  for (int i = 0; i < nVertex; i++) {
    pvtx_ln_to_gn[i] = distrib_vtx[rank] + i + 1;
  }

  PDM_compute_face_edge_from_face_vtx(PDM_MPI_mpi_2_pdm_mpi_comm(&comm),
                                      1,
                                      &nElts,
                                      &nVertex,
                                      &eltsConnecPointer[0],
                                      &eltsConnec[0],
                                      &pface_ln_to_gn,
                                      &pvtx_ln_to_gn,
                                      &part_pface_edge_idx,
                                      &part_pface_edge,
                                      &part_pn_edge,
                                      &part_pedge_vtx,
                                      &part_pedge_ln_to_gn);

  int *edge_vtx_idx = malloc(sizeof(int) * (part_pn_edge[0] + 1));
  edge_vtx_idx[0] = 0;
  for (int i = 1; i < part_pn_edge[0] + 1; i++) {
    edge_vtx_idx[i] = edge_vtx_idx[i-1] + 2;
  }

  CWP_client_Mesh_interf_from_faceedge_set(code_name[0],
                                           cpl_name,
                                           0,
                                           nElts,
                                           part_pface_edge_idx[0],
                                           part_pface_edge[0],
                                           part_pn_edge[0],
                                           edge_vtx_idx,
                                           part_pedge_vtx[0],
                                           global_num);

  printf("%d - Delete mesh\n", rank);
  CWP_client_Mesh_interf_del(code_name[0], cpl_name);

  printf("%d - Delete coupling\n", rank);
  CWP_client_Cpl_del(code_name[0], cpl_name);

  // Freeing memory
  free(coords[0]           );
  free(eltsConnecPointer[0]);
  free(eltsConnec[0]       );

  free(coords);
  free(eltsConnecPointer);
  free(eltsConnec       );

  free(sendValues[0]);
  free(sendValues);

  free(code_name);
  free(coupled_code_name);
  free(is_active_rank);
  free(time_init);

  free(global_num);
  free(global_num_vtx);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn );
  free(distrib_elt);
  free(distrib_vtx);
  free(edge_vtx_idx );

  free(part_pface_edge_idx);
  free(part_pface_edge    );
  free(part_pn_edge       );
  free(part_pedge_vtx     );
  free(part_pedge_ln_to_gn);


  // Finalize
  CWP_client_Finalize();

  MPI_Finalize();
  return EXIT_SUCCESS;
}
