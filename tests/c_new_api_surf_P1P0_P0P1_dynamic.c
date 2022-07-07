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
    if (strcmp(argv[i], "-h") == 0) {
      _display_usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *nVertex = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      i++;
      if (i >= argc) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *randLevel = atof(argv[i]);
      }
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

  srand(rank + time(0));

  int n_partition = 0;
  while (1.5 * (double) pow(n_partition, 2) < commWorldSize) {
    n_partition++;
  }
  int n2 = (int) (1.5 * pow(n_partition, 2));

  if (n2 != commWorldSize) {
    if (rank == 0) {
      printf("Not executed : only available if the number of processus in the form of '1.5 * n^2' \n");
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  // Read args from command line
  int nVertexSeg = 30;
  double randLevel = 0.4;
  _read_args(argc, argv, &nVertexSeg, &randLevel);

  // Initialization
  int n_code = 0;
  const char **code_name = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;
  const char **coupled_code_name = NULL;

  int rankCode1[commWorldSize];
  int rankCode2[commWorldSize];

  // Random distribution of codes on processes.
  if (rank == 0) {
    time_t t;
    // Intializes random number generator
    srand((unsigned) time(&t));
    int size_code1_domain = (int) (commWorldSize * 2.0 / 3.0);
    int size_code2_domain = (int) (commWorldSize * 2.0 / 3.0);

    if ((int) sqrt(size_code1_domain) != (int) sqrt(size_code1_domain)) {
      size_code1_domain = (int) pow(sqrt(size_code1_domain), 2.0);
    }

    if ((int) sqrt(size_code2_domain) != (int) sqrt(size_code2_domain)) {
      size_code2_domain = (int) pow((int) sqrt(size_code2_domain), 2.0);
    }

    for (int i = 0 ; i < commWorldSize ; i++) {
      rankCode1[i] = 0;
      rankCode2[i] = 0;
    }
    for (int i = 0 ; i < size_code1_domain ; i++) {
      int tmp = i; //rand() % commWorldSize;
      while (rankCode1[tmp] == 1) {
        tmp = rand() % commWorldSize;
      }
      rankCode1[tmp] = 1;
    }

    for (int i = 0 ; i < size_code2_domain ; i++) {
      int tmp = size_code1_domain / 2 + i; //rand() % commWorldSize;
      while (rankCode2[tmp] == 1) {
        tmp = rand() % commWorldSize;
      }
      rankCode2[tmp] = 1;
    }

    int ind = 0;
    for (int i = 0 ; i < commWorldSize ; i++) {
      if (rankCode2[i] || rankCode1[i]) {
//        nranks[ind] = i;
        ind++;
      }
    }
  }

  int rankCode1Rcv = -1;
  int rankCode2Rcv = -1;

  MPI_Scatter(&rankCode1, 1, MPI_INT, &rankCode1Rcv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(&rankCode2, 1, MPI_INT, &rankCode2Rcv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(&rankCode2, 1, MPI_INT, &rankCode2Rcv, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rankCode1Rcv && rankCode2Rcv) {
    printf("%d - Working for code1 and code2\n", rank);
  }
  else {
    if (rankCode1Rcv) {
      printf("%d - Working for code1\n", rank);
    }
    if (rankCode2Rcv) {
      printf("%d - Working for code2\n", rank);
    }
  }

  if (rankCode1Rcv) {
    n_code = 1;
    code_name = malloc(sizeof(char *) * n_code);
    coupled_code_name = malloc(sizeof(char *) * n_code);
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }

  if (rankCode2Rcv) {
    n_code = 1;
    code_name = malloc(sizeof(char *) * n_code);
    coupled_code_name = malloc(sizeof(char *) * n_code);
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }

  if (rankCode1Rcv && rankCode2Rcv) {
    n_code = 2;
    code_name = malloc(sizeof(char *) * n_code);
    coupled_code_name = malloc(sizeof(char *) * n_code);
    code_name[0] = "code1";
    code_name[1] = "code2";
    coupled_code_name[0] = "code2";
    coupled_code_name[1] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
  }

  times_init = malloc(sizeof(double) * n_code);
  for (int i = 0 ; i < n_code ; i++) {
    times_init[i] = 0;
  }

  MPI_Comm *localComm = malloc(sizeof(MPI_Comm) * n_code);
  CWP_Init(MPI_COMM_WORLD,
           n_code,
           code_name,
           is_coupled_rank,
           times_init,
           localComm);

  // Output redirection
  int currentRank[n_code];
  int localCommSize[n_code];

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (localComm[i_code] != MPI_COMM_NULL) {
      MPI_Comm_rank(localComm[i_code], &currentRank[i_code]);
    }
    if (localComm[i_code] != MPI_COMM_NULL) {
      MPI_Comm_size(localComm[i_code], &localCommSize[i_code]);
    }
  }

  // Coupling creation
  printf("%d - Create coupling\n", rank);
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1_dynamic";
  int nb_part = 1;
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Cpl_create(code_name[i_code],                                     // Code name
                   cpl_name,                                              // Coupling id
                   coupled_code_name[i_code],                             // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,                                // Coupling type
                   CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                   nb_part,                                               // Partition number
                   CWP_DYNAMIC_MESH_DEFORMABLE,                           // Mesh displacement type
                   CWP_TIME_EXCH_CPL_TIME_STEP);                          // Postprocessing frequency
  }

  printf("%d - Set visu\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Visu_set(code_name[i_code],       // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

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

  srand(time(NULL));

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    coords[i_code] = (double *) malloc(sizeof(double) * 3 * nVertex);
    eltsConnecPointer[i_code] = (int *) malloc(sizeof(int) * (nElts + 1));
    eltsConnec[i_code] = (int *) malloc(sizeof(int) * 4 * nElts);
    grid_mesh(xmin,
              xmax,
              ymin,
              ymax,
              randLevel,
              nVertexSeg,
              (int) sqrt(localCommSize[i_code]),
              coords[i_code],
              eltsConnecPointer[i_code],
              eltsConnec[i_code],
              localComm[i_code]);

    carre2rond(xmin, xmax, ymin, ymax, coords[i_code], nVertex);
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Mesh_interf_vtx_set(code_name[i_code], cpl_name, 0, nVertex, coords[i_code], NULL);

    int block_id = CWP_Mesh_interf_block_add(code_name[i_code], cpl_name, CWP_BLOCK_FACE_POLY);

    CWP_Mesh_interf_f_poly_block_set(code_name[i_code],
                                     cpl_name,
                                     0,
                                     block_id,
                                     nElts,
                                     eltsConnecPointer[i_code],
                                     eltsConnec[i_code],
                                     NULL);
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Mesh_interf_finalize(code_name[i_code], cpl_name);
  }

  printf("%d - Exchange Code1 <-> Code2\n", rank);

  double **sendValues = (double **) malloc(sizeof(double *) * n_code);
  double **recvValues = (double **) malloc(sizeof(double *) * n_code);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (strcmp(code_name[i_code], "code1") == 0) {
      sendValues[i_code] = (double *) malloc(sizeof(double) * 3 * nVertex);
      recvValues[i_code] = (double *) malloc(sizeof(double) * nElts);
      for (int i = 0 ; i < nVertex ; i++) {
        sendValues[i_code][3 * i] = coords[i_code][3 * i];
        sendValues[i_code][3 * i + 1] = coords[i_code][3 * i + 1];
        sendValues[i_code][3 * i + 2] = coords[i_code][3 * i + 2];
      }
    }
    else {
      sendValues[i_code] = (double *) malloc(sizeof(double) * nElts);
      recvValues[i_code] = (double *) malloc(sizeof(double) * 3 * nVertex);
      for (int i = 0 ; i < nElts ; i++) {
        sendValues[i_code][i] = rank;
      }
    }
  }

  // Exchange
  // field1: code1 -> code2
  // field2: code2 -> code1
  const char *fieldName1 = "cooX_t0";;
  const char *fieldName2 = "code2_elt_rank";;

  CWP_Status_t visu_status = CWP_STATUS_ON;

  printf("%d - Defining fields\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (strcmp(code_name[i_code], "code1") == 0) {
      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       3,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       fieldName2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         fieldName1,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         sendValues[i_code]);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         fieldName2,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recvValues[i_code]);
    }
    else {
      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       3,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       fieldName2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         fieldName2,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         sendValues[i_code]);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         fieldName1,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         recvValues[i_code]);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double recv_time = 0.;
  for (int il_nb_ite = 0 ; il_nb_ite < 1 ; il_nb_ite++) {
    recv_time += 1.0;

    // Mesh rotation and new localisation
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      if (strcmp(code_name[i_code], "code2") == 0) {
        mesh_rotate(coords[i_code], nVertex, 3 * recv_time);
      }
      if (strcmp(code_name[i_code], "code1") == 0) {
        mesh_rotate(coords[i_code], nVertex, recv_time);
      }

      CWP_Spatial_interp_weights_compute(code_name[i_code], cpl_name);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      if (strcmp(code_name[i_code], "code1") == 0) {
        CWP_Field_issend(code_name[i_code], cpl_name, fieldName1);
        CWP_Field_irecv(code_name[i_code], cpl_name, fieldName2);
      }
      else {
        CWP_Field_irecv(code_name[i_code], cpl_name, fieldName1);
        CWP_Field_issend(code_name[i_code], cpl_name, fieldName2);
      }
    }

    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      if (strcmp(code_name[i_code], "code1") == 0) {
        CWP_Field_wait_issend(code_name[i_code], cpl_name, fieldName1);
        CWP_Field_wait_irecv(code_name[i_code], cpl_name, fieldName2);
      }
      else {
        CWP_Field_wait_irecv(code_name[i_code], cpl_name, fieldName1);
        CWP_Field_wait_issend(code_name[i_code], cpl_name, fieldName2);
      }
    }
  }

  printf("%d - Delete mesh\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Mesh_interf_del(code_name[i_code], cpl_name);
  }

  printf("%d - Delete coupling\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Cpl_del(code_name[i_code], cpl_name);
  }

  // Freeing memory
  free(coords);
  free(sendValues);
  free(recvValues);

  // Finalize
  CWP_Finalize();
  MPI_Finalize();
  return EXIT_SUCCESS;
}
