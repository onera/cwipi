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
#include "grid_mesh.h"


/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_display_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -n     <level>  Number of vertices in band width.\n\n"
         "  -rand  <level>  Random level ( > 0 and < 0.4) \n\n"
         "  -h             this message.\n\n");
  exit(exit_code);
}


/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nVertex             <-- Number of vertices in bandwidth
 *   randLevel           <-- Random level
 *---------------------------------------------------------------------*/

static void
_read_args(int argc, char **argv, int *nVertex, double *randLevel) {
  int i = 1;

  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {_display_usage(EXIT_SUCCESS);}
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {_display_usage(EXIT_FAILURE);}
      else {*nVertex = atoi(argv[i]);}
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      i++;
      if (i >= argc) {_display_usage(EXIT_FAILURE);}
      else {*randLevel = atof(argv[i]);}
    }
    i++;
  }
}


/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank;
  int commWorldSize;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

  // Read args from command line
  int nVertexSeg = 10;
  double randLevel = 0.4;

  _read_args(argc, argv, &nVertexSeg, &randLevel);

  // Initialization
  int n_code = 1;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t *is_active_rank = malloc(sizeof(CWP_Status_t) * n_code);
  double *time_init = malloc(sizeof(double) * n_code);
  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
  MPI_Comm *connectableLocalComm = malloc(sizeof(MPI_Comm) * n_code);
  int *connectableLocalCommSize = malloc(sizeof(int) * n_code);

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

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           time_init,
           intra_comm);

  printf("%d - Create coupling\n", rank);
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1";
  int nb_part = 1;
  CWP_Cpl_create(code_name[0],                                          // Code name
                 cpl_name,                                              // Coupling id
                 coupled_code_name[0],                                  // Coupled application id
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,                                // Coupling type
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                 nb_part,                                               // Number of partitions
                 CWP_DYNAMIC_MESH_STATIC,                               // Mesh type
                 CWP_TIME_EXCH_CPL_TIME_STEP);                          // Postprocessing frequency

  printf("%d - Set visu\n", rank);
  CWP_Visu_set(code_name[0],            // Code name
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

  coords = (double **) malloc(sizeof(double *) * n_code);
  eltsConnecPointer = (int **) malloc(sizeof(int *) * n_code);
  eltsConnec = (int **) malloc(sizeof(int *) * n_code);

  coords[0] = (double *) malloc(sizeof(double) * 3 * nVertex);
  eltsConnecPointer[0] = (int *) malloc(sizeof(int) * (nElts + 1));
  eltsConnec[0] = (int *) malloc(sizeof(int) * 4 * nElts);
  randLevel = 0.4;

  connectableLocalComm[0] = CWP_Connectable_comm_get((char *) code_name[0]);
  MPI_Comm_size(connectableLocalComm[0], &connectableLocalCommSize[0]);

  srand(time(NULL));

  if (strcmp(code_name[0], "code1") == 0) {
    grid_mesh(xmin,
              xmax,
              ymin,
              ymax,
              randLevel,
              nVertexSeg,
              (int) sqrt(connectableLocalCommSize[0]),
              coords[0],
              eltsConnecPointer[0],
              eltsConnec[0],
              connectableLocalComm[0]);
  }

  randLevel = 0.2;
  if (strcmp(code_name[0], "code2") == 0) {
    grid_mesh(xmin,
              xmax,
              ymin,
              ymax,
              randLevel,
              nVertexSeg,
              (int) sqrt(connectableLocalCommSize[0]),
              coords[0],
              eltsConnecPointer[0],
              eltsConnec[0],
              connectableLocalComm[0]);
  }

  printf("%d - Number of vertex   : %d\n", rank, nVertex);
  printf("%d - Number of elements : %d\n", rank, nElts);

  CWP_Mesh_interf_vtx_set(code_name[0], cpl_name, 0, nVertex, coords[0], NULL);

  int block_id = CWP_Mesh_interf_block_add(code_name[0], cpl_name, CWP_BLOCK_FACE_POLY);

  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   cpl_name,
                                   0,
                                   block_id,
                                   nElts,
                                   eltsConnecPointer[0],
                                   eltsConnec[0],
                                   NULL);

  CWP_Mesh_interf_finalize(code_name[0], cpl_name);

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
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_BLOCK,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_BLOCK,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);

    CWP_Field_data_set(code_name[0], cpl_name, field_name1, 0, CWP_FIELD_MAP_SOURCE, sendValues[0]);
    CWP_Field_data_set(code_name[0], cpl_name, field_name2, 0, CWP_FIELD_MAP_TARGET, recvValues[0]);
  }
  
  else if (strcmp(code_name[0], "code2") == 0) {
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_BLOCK,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_BLOCK,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_Field_data_set(code_name[0], cpl_name, field_name2, 0, CWP_FIELD_MAP_SOURCE, sendValues[0]);
    CWP_Field_data_set(code_name[0], cpl_name, field_name1, 0, CWP_FIELD_MAP_TARGET, recvValues[0]);
  }

  printf("%d - Before compute\n", rank);
  CWP_Spatial_interp_weights_compute(code_name[0], cpl_name);

  int n_uncomputed = 0;

  if (strcmp(code_name[0], "code1") == 0) {
    n_uncomputed = CWP_N_uncomputed_tgts_get (code_name[0], cpl_name, field_name2, 0);
  }

  else if (strcmp(code_name[0], "code2") == 0) {
    n_uncomputed = CWP_N_uncomputed_tgts_get (code_name[0], cpl_name, field_name1, 0);
  }

  printf("%d - After compute %d\n", rank, n_uncomputed);

  if (strcmp(code_name[0], "code1") == 0) {
    CWP_Field_issend(code_name[0], cpl_name, field_name1);
    CWP_Field_irecv(code_name[0], cpl_name, field_name2);
  }
  else if (strcmp(code_name[0], "code2") == 0) {
    CWP_Field_irecv(code_name[0], cpl_name, field_name1);
    CWP_Field_issend(code_name[0], cpl_name, field_name2);
  }

  if (strcmp(code_name[0], "code1") == 0) {
    CWP_Field_wait_issend(code_name[0], cpl_name, field_name1);
    CWP_Field_wait_irecv(code_name[0], cpl_name, field_name2);
  }
  else if (strcmp(code_name[0], "code2") == 0) {
    CWP_Field_wait_irecv(code_name[0], cpl_name, field_name1);
    CWP_Field_wait_issend(code_name[0], cpl_name, field_name2);
  }

  printf("%d - Delete mesh\n", rank);
  CWP_Mesh_interf_del(code_name[0], cpl_name);

  printf("%d - Delete coupling\n", rank);
  CWP_Cpl_del(code_name[0], cpl_name);

  // Freeing memory
  free(coords);
  free(sendValues);
  free(recvValues);

  // Finalize
  CWP_Finalize();
  MPI_Finalize();
  return EXIT_SUCCESS;
}
