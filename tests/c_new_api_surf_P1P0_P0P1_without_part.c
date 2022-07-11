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
#include <time.h>
#include <math.h>

#include "cwp.h"
#include "grid_mesh.h"
#include "cwp_priv.h"


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

  // Parse and check command line
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


static void
_userInterpolation(const int interface_type, const int n_src_vtcs, const int n_src_std_elts,
                   const int n_src_poly, const int n_tgt_pts, const double src_vtcs_coords[],
                   const CWP_g_num_t src_global_elts_num[], const CWP_g_num_t src_global_vtcs_num[],
                   const int src_connec_idx[], const int src_connec[],
                   const int src_poly_cell_face_idx[], const int src_poly_cell_face_connec[],
                   const int src_poly_face_vtx_idx[], const int src_poly_face_vtx_connec[],
                   const double tgt_pts_coords[], const int tgt_pts_target_location[],
                   const double tgt_pts_dist[], const int tgt_pts_bary_coords_idx[],
                   const double tgt_pts_bary_coords[], const int stride,
                   const CWP_Dof_location_t src_field_location, const void *src_field,
                   const CWP_Dof_location_t tgt_field_location, void *tgt_field) {
  CWP_UNUSED (interface_type);
  CWP_UNUSED (n_src_vtcs);
  CWP_UNUSED (n_src_std_elts);
  CWP_UNUSED (n_src_poly);
  CWP_UNUSED (src_vtcs_coords);
  CWP_UNUSED (src_global_elts_num);
  CWP_UNUSED (src_global_vtcs_num);
  CWP_UNUSED (src_connec_idx);
  CWP_UNUSED (src_connec);
  CWP_UNUSED(src_poly_cell_face_idx);
  CWP_UNUSED(src_poly_cell_face_connec);
  CWP_UNUSED(src_poly_face_vtx_idx);
  CWP_UNUSED(src_poly_face_vtx_connec);
  CWP_UNUSED (tgt_pts_coords);
  CWP_UNUSED (tgt_pts_target_location);
  CWP_UNUSED (tgt_pts_dist);
  CWP_UNUSED (tgt_pts_bary_coords_idx);
  CWP_UNUSED (tgt_pts_bary_coords);
  CWP_UNUSED (stride);
  CWP_UNUSED (tgt_field_location);

  if (src_field_location == CWP_DOF_LOCATION_NODE) {
    for (int i = 0 ; i < n_tgt_pts ; i++) {

//      int ielt = tgt_pts_target_location[i];

      //      int ivertex[4];
      //      ivertex[0] = src_connec[src_connec_idx[ielt]];
      //      ivertex[1] = src_connec[src_connec_idx[ielt] + 1];
      //      ivertex[2] = src_connec[src_connec_idx[ielt] + 2];
      //      ivertex[3] = src_connec[src_connec_idx[ielt] + 3];

//      double *shapef = (double *) malloc(4 * n_src_std_elts * sizeof(double));
//      double *shapef_elt = shapef + 4 * ielt;

      //  printf("n_tgt_pts USER_INTERP %i\n",n_tgt_pts);
      ((double *) tgt_field)[i] = 0.0;
      for (int k = 0 ; k < 4 ; k++) {
        ((double *) tgt_field)[i] = ((double *) src_field)[0];//ivertex[k]];
      }

//      free(shapef);
    }
  }
  else {
    printf("Error in userInterpolation : bad solver_type\n");
    exit(EXIT_FAILURE);
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

  if (pow(n_partition, 2) > commWorldSize) {
    n_partition--;
  }

  int size_code1 = (int) (pow(n_partition, 2));
  int size_code2 = (int) (pow(n_partition, 2));

  // Read args from command line
  int nVertexSeg = 100;
  double randLevel = 0.4;

  _read_args(argc, argv, &nVertexSeg, &randLevel);

  // Initialization
  int n_code = 0;
  const char **codeName = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;

  const char **codeCoupledName = NULL;

  int rankCode1[commWorldSize];
  int rankCode2[commWorldSize];

  // Random distribution of codes on processes.
  if (rank == 0) {
    // Intializes random number generator
    srand((unsigned) time(NULL));

    int size_code1_domain = size_code1;
    int size_code2_domain = size_code2;
    printf("size_code1 %i\n", size_code1);
    for (int i = 0 ; i < commWorldSize ; i++) {
      rankCode1[i] = 0;
      rankCode2[i] = 0;
    }
    for (int i = 0 ; i < size_code1_domain ; i++) {
      int tmp = rand() % commWorldSize;
      while (rankCode1[tmp] == 1) {
        tmp = rand() % commWorldSize;
      }
      rankCode1[tmp] = 1;
    }

    for (int i = 0 ; i < size_code2_domain ; i++) {
      int tmp = rand() % commWorldSize;
      while (rankCode2[tmp] == 1) {
        tmp = rand() % commWorldSize;
      }
      rankCode2[tmp] = 1;
    }

    for (int i = 0 ; i < commWorldSize ; i++) {
      if (rankCode1[i] == 0 && rankCode2[i] == 0) {
        int t = rand() % commWorldSize;
        while (!(rankCode1[t] == 1 && rankCode2[t] == 1)) {
          t = rand() % commWorldSize;
        }
        int h12 = rand() % 2;
        if (h12 == 0) {
          rankCode1[t] = 0;
          rankCode1[i] = 1;
        }
        else {
          rankCode2[t] = 0;
          rankCode2[i] = 1;
        }
      }
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
    printf("I am process %i and I own the code1 and the code2.\n", rank);
  }
  else {
    if (rankCode1Rcv) {
      printf("maxI am process %i and I own the code1.\n", rank);
    }
    if (rankCode2Rcv) {
      printf("maxI am process %i and I own the code2.\n", rank);
    }
  }

  int *nbPartSeg = NULL;
  int *nbPart = NULL;

  if (rankCode1Rcv) {
    n_code = 1;
    codeName = malloc(sizeof(char *) * n_code);
    codeCoupledName = malloc(sizeof(char *) * n_code);
    codeName[0] = "code1";
    codeCoupledName[0] = "code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
    nbPartSeg = (int *) malloc(sizeof(int) * n_code);
    nbPartSeg[0] = 1;
    nbPart = (int *) malloc(sizeof(int) * n_code);
    nbPart[0] = nbPartSeg[0] * nbPartSeg[0];
  }

  if (rankCode2Rcv) {
    n_code = 1;
    codeName = malloc(sizeof(char *) * n_code);
    codeCoupledName = malloc(sizeof(char *) * n_code);
    codeName[0] = "code2";
    codeCoupledName[0] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
    nbPartSeg = (int *) malloc(sizeof(int) * n_code);
    nbPartSeg[0] = 1;
    nbPart = (int *) malloc(sizeof(int) * n_code);
    nbPart[0] = nbPartSeg[0] * nbPartSeg[0];
  }

  if (rankCode1Rcv && rankCode2Rcv) {
    n_code = 2;
    codeName = malloc(sizeof(char *) * n_code);
    codeCoupledName = malloc(sizeof(char *) * n_code);
    codeName[0] = "code1";
    codeName[1] = "code2";
    codeCoupledName[0] = "code2";
    codeCoupledName[1] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
    nbPartSeg = (int *) malloc(sizeof(int) * n_code);
    nbPartSeg[0] = 1;
    nbPartSeg[1] = 1;
    nbPart = (int *) malloc(sizeof(int) * n_code);
    nbPart[0] = nbPartSeg[0] * nbPartSeg[0];
    nbPart[1] = nbPartSeg[1] * nbPartSeg[1];
  }

  times_init = malloc(sizeof(double) * n_code);

  for (int i = 0 ; i < n_code ; i++) {
    times_init[i] = 0;
  }

  MPI_Comm *localComm = malloc(sizeof(MPI_Comm) * n_code);
  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) codeName,
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

  const char *cpl_name = "c_new_api_surf_P1P0_P0P1_without_part";
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    printf("        Create coupling %i code %i\n", rank, i_code);
    CWP_Cpl_create(codeName[i_code],                                      // Code name
                   cpl_name,                                              // Coupling id
                   codeCoupledName[i_code],                               // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITHOUT_PART,                             // Coupling type
                   CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                   nbPart[i_code],                                        // Partition number
                   CWP_DYNAMIC_MESH_STATIC,                               // Mesh type
                   CWP_TIME_EXCH_USER_CONTROLLED);                           // frequency
    printf("   After     Create coupling %i code %i\n", rank, i_code);
  }

  printf("        Set visu\n");
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Visu_set(codeName[i_code],        // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

  // Mesh definition
  printf("        Create mesh %i\n", rank);

  int **nVertex = NULL;            // Number of vertex
  double ***coords = NULL;         // Vertex coordinates
  int **nElts = NULL;              // Number of elements
  int ***eltsConnecPointer = NULL; // Connectivity index
  int ***eltsConnec = NULL;        // Connectivity

  // Domain bounds
  const double xmin = -10;
  const double xmax = 10;
  const double ymin = -10;
  const double ymax = 10;

  coords = (double ***) malloc(sizeof(double **) * n_code);
  eltsConnecPointer = (int ***) malloc(sizeof(int **) * n_code);
  eltsConnec = (int ***) malloc(sizeof(int **) * n_code);
  nVertex = (int **) malloc(sizeof(int *) * n_code);
  nElts = (int **) malloc(sizeof(int *) * n_code);

  srand(time(NULL));

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    coords[i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    eltsConnecPointer[i_code] = (int **) malloc(sizeof(int *) * nbPart[i_code]);
    eltsConnec[i_code] = (int **) malloc(sizeof(int *) * nbPart[i_code]);
    nVertex[i_code] = (int *) malloc(sizeof(int) * nbPart[i_code]);
    nElts[i_code] = (int *) malloc(sizeof(int) * nbPart[i_code]);
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    randLevel = 0.4;
    double xSegPart = (xmax - xmin) / (double) nbPartSeg[i_code];
    double ySegPart = (ymax - ymin) / (double) nbPartSeg[i_code];
    for (int u = 0 ; u < nbPartSeg[i_code] ; u++) {
      for (int v = 0 ; v < nbPartSeg[i_code] ; v++) {
        int i_part = nbPartSeg[i_code] * u + v;
        int nVertexSegPart = nVertexSeg / nbPartSeg[i_code];
        nVertex[i_code][i_part] = nVertexSegPart * nVertexSegPart;
        nElts[i_code][i_part] = (nVertexSegPart - 1) * (nVertexSegPart - 1);

        coords[i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nVertex[i_code][i_part]);
        eltsConnecPointer[i_code][i_part] = (int *) malloc(sizeof(int) * (nElts[i_code][i_part] + 1));
        eltsConnec[i_code][i_part] = (int *) malloc(sizeof(int) * 4 * nElts[i_code][i_part]);
        if (strcmp(codeName[i_code], "code1") == 0) {
          double xminPart = xmin + xSegPart * v;
          double yminPart = ymin + ySegPart * u;
          double xmaxPart = xminPart + xSegPart;
          double ymaxPart = yminPart + ySegPart;
          printf("INFO rank %i i_part %i x %f %f y %f %f nbPart %i test %i nVertex %i \n",
                 rank,
                 i_part,
                 xminPart,
                 yminPart,
                 xmaxPart,
                 ymaxPart,
                 nbPart[i_code],
                 nVertexSeg / nbPartSeg[i_code],
                 nVertex[i_code][i_part]);

          grid_mesh(xminPart,
                    xmaxPart,
                    yminPart,
                    ymaxPart,
                    randLevel,
                    nVertexSegPart,
                    (int) sqrt(localCommSize[i_code]),
                    coords[i_code][i_part],
                    eltsConnecPointer[i_code][i_part],
                    eltsConnec[i_code][i_part],
                    localComm[i_code]);
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    randLevel = 0.2;
    double xSegPart = (xmax - xmin) / (double) nbPartSeg[i_code];
    double ySegPart = (ymax - ymin) / (double) nbPartSeg[i_code];
    for (int u = 0 ; u < nbPartSeg[i_code] ; u++) {
      for (int v = 0 ; v < nbPartSeg[i_code] ; v++) {
        int i_part = nbPartSeg[i_code] * u + v;
        int nVertexSegPart = nVertexSeg / nbPartSeg[i_code];
        if (strcmp(codeName[i_code], "code2") == 0) {
          double xminPart = xmin + xSegPart * v;
          double yminPart = ymin + ySegPart * u;
          double xmaxPart = xminPart + xSegPart;
          double ymaxPart = yminPart + ySegPart;
          printf("INFO rank %i i_part %i xmin %f ymin %f xmax %f ymax %f nbPart %i test %i nVertex %i \n",
                 rank,
                 i_part,
                 xminPart,
                 yminPart,
                 xmaxPart,
                 ymaxPart,
                 nbPart[i_code],
                 nVertexSeg / nbPartSeg[i_code],
                 nVertex[i_code][i_part]);

          grid_mesh(xminPart,
                    xmaxPart,
                    yminPart,
                    ymaxPart,
                    randLevel,
                    nVertexSegPart,
                    (int) sqrt(localCommSize[i_code]),
                    coords[i_code][i_part],
                    eltsConnecPointer[i_code][i_part],
                    eltsConnec[i_code][i_part],
                    localComm[i_code]);
        }
      }
    }
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    for (int i_part = 0 ; i_part < nbPart[i_code] ; i_part++) {
      CWP_Mesh_interf_vtx_set(codeName[i_code],
                              cpl_name,
                              i_part,
                              nVertex[i_code][i_part],
                              coords[i_code][i_part],
                              NULL);
    }
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    int block_id = CWP_Mesh_interf_block_add(codeName[i_code], cpl_name, CWP_BLOCK_FACE_POLY);

    for (int i_part = 0 ; i_part < nbPart[i_code] ; i_part++) {
      CWP_Mesh_interf_f_poly_block_set(codeName[i_code],
                                       cpl_name,
                                       i_part,
                                       block_id,
                                       nElts[i_code][i_part],
                                       eltsConnecPointer[i_code][i_part],
                                       eltsConnec[i_code][i_part],
                                       NULL);
    }
  }

  printf("       Finalize %i\n", rank);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Mesh_interf_finalize(codeName[i_code], cpl_name);
  }

  // Fields exchange
  printf("        Exchange Code1 <-> Code2 %i\n", rank);

  double ***sendValues = (double ***) malloc(sizeof(double *) * n_code);
  double ***recvValues = (double ***) malloc(sizeof(double *) * n_code);
  double ***Values2Vertex = (double ***) malloc(sizeof(double *) * n_code);
  double ***recvValuesUser = (double ***) malloc(sizeof(double *) * n_code);

  int nbPoints = 3;
  int **nbPointsUser = (int **) malloc(sizeof(int *) * n_code);

  double ***coordsPointsUser = (double ***) malloc(sizeof(double **) * n_code);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    sendValues[i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    nbPointsUser[i_code] = (int *) malloc(sizeof(int) * nbPart[i_code]);
    recvValues[i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    Values2Vertex[i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    recvValuesUser[i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    coordsPointsUser[i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    for (int i_part = 0 ; i_part < nbPart[i_code] ; i_part++) {

      if (strcmp(codeName[i_code], "code1") == 0) {
        sendValues[i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);
        recvValues[i_code][i_part] = (double *) malloc(sizeof(double) * nElts[i_code][i_part]);
        Values2Vertex[i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);
        nbPointsUser[i_code][i_part] = 0;
        recvValuesUser[i_code][i_part] = (double *) malloc(sizeof(double) * nbPointsUser[i_code][i_part]);
        coordsPointsUser[i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nbPointsUser[i_code][i_part]);
        for (int i = 0 ; i < nVertex[i_code][i_part] ; i++) {
          sendValues[i_code][i_part][i] = rank;
        }
      }
      else {
        sendValues[i_code][i_part] = (double *) malloc(sizeof(double) * nElts[i_code][i_part]);
        recvValues[i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nVertex[i_code][i_part]);
        Values2Vertex[i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);
        nbPointsUser[i_code][i_part] = nbPoints;
        recvValuesUser[i_code][i_part] = (double *) malloc(sizeof(double) * nbPointsUser[i_code][i_part]);
        coordsPointsUser[i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nbPointsUser[i_code][i_part]);
        for (int i = 0 ; i < nElts[i_code][i_part] ; i++) {
          sendValues[i_code][i_part][i] = rank;
        }
        for (int i = 0 ; i < nVertex[i_code][i_part] ; i++) {
          Values2Vertex[i_code][i_part][i] = i;
        }

        for (int i = 0 ; i < nbPointsUser[i_code][i_part] ; i++) {
          coordsPointsUser[i_code][i_part][3 * i] = 0.0 + 0.01 * i_part + 0.001 * i;
          coordsPointsUser[i_code][i_part][3 * i + 1] = rank * 0.1 + 0.01 * i_part + 0.001 * i;
          coordsPointsUser[i_code][i_part][3 * i + 2] = 0.0;
        }
      }
    }
  }

  // Define fields
  int code1I = 1;
  int code2I = 1;
  int code3I = 1;
  int code4I = 1;
  int code5I = 1;

  const char *fieldName1;
  const char *fieldName2;
  const char *fieldName3;
  fieldName1 = "cooX";
  fieldName2 = "rank";
  fieldName3 = "userField";
  const char *fieldName4 = "userInterpolation";
  const char *fieldName5 = "sendRecvField";

  printf("        Field %i\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (strcmp(codeName[i_code], "code1") == 0) {
      if (code1I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName1,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_NODE,
                         CWP_FIELD_EXCH_SEND,
                         CWP_STATUS_ON);
      }
      if (code2I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName2,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_CELL_CENTER,
                         CWP_FIELD_EXCH_RECV,
                         CWP_STATUS_ON);
      }
      if (code3I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName3,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_NODE,
                         CWP_FIELD_EXCH_SEND,
                         CWP_STATUS_OFF);
      }
      if (code4I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName4,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_NODE,
                         CWP_FIELD_EXCH_SEND,
                         CWP_STATUS_ON);
      }
      if (code5I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName5,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_NODE,
                         CWP_FIELD_EXCH_RECV,
                         CWP_STATUS_OFF);
      }

      for (int i_part = 0 ; i_part < nbPart[i_code] ; i_part++) {
        if (code1I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName1,
                             i_part,
                             CWP_FIELD_MAP_SOURCE,
                             sendValues[i_code][i_part]);
        }
        if (code2I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName2,
                             i_part,
                             CWP_FIELD_MAP_TARGET,
                             recvValues[i_code][i_part]);
        }
        if (code3I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName3,
                             i_part,
                             CWP_FIELD_MAP_SOURCE,
                             sendValues[i_code][i_part]);
        }
        if (code4I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName4,
                             i_part,
                             CWP_FIELD_MAP_SOURCE,
                             sendValues[i_code][i_part]);
        }
        if (code5I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName5,
                             i_part,
                             CWP_FIELD_MAP_TARGET,
                             Values2Vertex[i_code][i_part]);
        }
      }

      if (code4I == 1) {
        CWP_Interp_from_location_set(codeName[i_code], cpl_name, fieldName4, _userInterpolation);
      }
    }
    else {
      if (code1I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName1,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_NODE,
                         CWP_FIELD_EXCH_RECV,
                         CWP_STATUS_ON);
      }
      if (code2I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName2,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_CELL_CENTER,
                         CWP_FIELD_EXCH_SEND,
                         CWP_STATUS_ON);
      }
      if (code3I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName3,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_USER,
                         CWP_FIELD_EXCH_RECV,
                         CWP_STATUS_OFF);
      }
      if (code4I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName4,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_NODE,
                         CWP_FIELD_EXCH_RECV,
                         CWP_STATUS_ON);
      }
      if (code5I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName5,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_NODE,
                         CWP_FIELD_EXCH_SEND,
                         CWP_STATUS_ON);
      }

      for (int i_part = 0 ; i_part < nbPart[i_code] ; i_part++) {
        if (code1I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName1,
                             i_part,
                             CWP_FIELD_MAP_TARGET,
                             Values2Vertex[i_code][i_part]);
        }
        if (code2I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName2,
                             i_part,
                             CWP_FIELD_MAP_SOURCE,
                             sendValues[i_code][i_part]);
        }
        if (code3I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName3,
                             i_part,
                             CWP_FIELD_MAP_TARGET,
                             recvValuesUser[i_code][i_part]);
        }
        if (code4I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName4,
                             i_part,
                             CWP_FIELD_MAP_TARGET,
                             recvValues[i_code][i_part]);
        }
        if (code5I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName5,
                             i_part,
                             CWP_FIELD_MAP_SOURCE,
                             Values2Vertex[i_code][i_part]);
        }
        if (code3I == 1) {
          CWP_User_tgt_pts_set(codeName[i_code],
                               cpl_name,
                               i_part,
                               nbPointsUser[i_code][i_part],
                               coordsPointsUser[i_code][i_part],
                               NULL);
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  printf("Before compute\n");
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Spatial_interp_weights_compute(codeName[i_code], cpl_name);
  }
  printf("After compute\n");

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (strcmp(codeName[i_code], "code1") == 0) {
      if (code1I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName1);
      }
      if (code1I == 1) {
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName1);
      }
    }
    else {
      if (code1I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName1);
      }
      if (code1I == 1) {
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName1);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (strcmp(codeName[i_code], "code1") == 0) {
      if (code2I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName2);
      }
      if (code2I == 1) {
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName2);
      }
    }
    else {
      if (code2I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName2);
      }
      if (code2I == 1) {
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName2);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (strcmp(codeName[i_code], "code1") == 0) {
      if (code3I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName3);
      }
      if (code3I == 1) {
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName3);
      }
    }
    else {
      if (code3I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName3);
      }
      if (code3I == 1) {
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName3);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (strcmp(codeName[i_code], "code1") == 0) {
      if (code4I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName4);
      }
      if (code4I == 1) {
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName4);
      }
    }
    else {
      if (code4I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName4);
      }
      if (code4I == 1) {
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName4);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (strcmp(codeName[i_code], "code1") == 0) {
      if (code5I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName5);
      }
      if (code5I == 1) {
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName5);
      }
    }
    else {
      if (code5I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName5);
      }
      if (code5I == 1) {
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName5);
      }
    }
  }

  // Coupling deletion
  printf("%d - Delete mesh\n", rank);
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Mesh_interf_del(codeName[i_code], cpl_name);
  }

  if (rank == 0) {
    printf("%d - Delete coupling\n", rank);
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Cpl_del(codeName[i_code], cpl_name);
  }

  // Freeing memory
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    for (int i_part = 0 ; i_part < nbPart[i_code] ; i_part++) {
      free(coords[i_code][i_part]);
      free(sendValues[i_code][i_part]);
      free(recvValues[i_code][i_part]);
      free(Values2Vertex[i_code][i_part]);
    }
    free(coords[i_code]);
    free(sendValues[i_code]);
    free(recvValues[i_code]);
    free(Values2Vertex[i_code]);
  }

  free(coords);
  free(sendValues);
  free(Values2Vertex);
  free(recvValues);

  // Finalize
  CWP_Finalize();
  MPI_Finalize();

  return EXIT_SUCCESS;
}
