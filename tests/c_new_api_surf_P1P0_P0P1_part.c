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
_userInterpolation
(
  const int                  interface_type,
  const char                *code_name,
  const int                  src_n_block,
  const CWP_Block_t          src_blocks_type[],
  const int                  src_i_part,
  const int                  src_n_vtx,
  const double               src_vtx_coords[],
  const CWP_g_num_t          src_vtx_global_num[],
  const int                  src_n_elts,
  const int                  src_i_block[],
  const int                  src_elt_in_block[],
  const int                  src_elt_vtx_idx[],
  const int                  src_elt_vtx[],
  const CWP_g_num_t          src_elts_global_num[],
  const int                  tgt_n_pts,
  const int                  tgt_pts_elt_idx[],
  const double               tgt_pts_coords[],
  const double               tgt_pts_dist[],
  const double               tgt_pts_uvw[],
  const int                  tgt_pts_weights_idx[],
  const double               tgt_pts_weights[],
  const int                  stride,
  const CWP_Dof_location_t   src_field_dof_location,
  const void                *src_field,
  void                      *tgt_field
)
{
  CWP_UNUSED(interface_type);
  CWP_UNUSED(code_name);
  CWP_UNUSED(src_n_block);
  CWP_UNUSED(src_blocks_type);
  CWP_UNUSED(src_i_part);
  CWP_UNUSED(src_n_vtx);
  CWP_UNUSED(src_vtx_coords);
  CWP_UNUSED(src_vtx_global_num);
  CWP_UNUSED(src_n_elts);
  CWP_UNUSED(src_i_block);
  CWP_UNUSED(src_elt_in_block);
  CWP_UNUSED(src_elt_vtx_idx);
  CWP_UNUSED(src_elt_vtx);
  CWP_UNUSED(src_elts_global_num);
  CWP_UNUSED(tgt_n_pts);
  CWP_UNUSED(tgt_pts_elt_idx);
  CWP_UNUSED(tgt_pts_coords);
  CWP_UNUSED(tgt_pts_dist);
  CWP_UNUSED(tgt_pts_uvw);
  CWP_UNUSED(tgt_pts_weights_idx);
  CWP_UNUSED(tgt_pts_weights);
  CWP_UNUSED(stride);
  CWP_UNUSED(src_field_dof_location);
  CWP_UNUSED(src_field);
  CWP_UNUSED(tgt_field);

  // Compute shapef
  // if (src_field_location == CWP_DOF_LOCATION_NODE) {
  //   for (int i = 0 ; i < n_tgt_pts ; i++) {
  //     int ielt = tgt_pts_target_location[i];
  //     double *shapef = (double *) malloc(4 * n_src_std_elts * sizeof(double));
  //     double *shapef_elt = shapef + 4 * ielt;
  //     CWP_UNUSED (shapef_elt);

  //     // Interpolation
  //     ((double *) tgt_field)[i] = 0.0;
  //     for (int k = 0 ; k < 4 ; k++) {
  //       ((double *) tgt_field)[i] = ((double *) src_field)[0];//ivertex[k]];
  //     }
  //     free(shapef);
  //   }
  // }
  // else {
  //   printf("Error in userInterpolation : bad solver_type\n");
  //   exit(EXIT_FAILURE);
  // }
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

  srand(rank);// + time(0));

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

  if (rank == 0) {
    time_t t;
    // Intializes random number generator
    srand((unsigned) time(&t));
    int size_code1_domain = size_code1;
    int size_code2_domain = size_code2;
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
        int t1 = rand() % commWorldSize;
        while (!(rankCode1[t1] == 1 && rankCode2[t1] == 1)) {
          t1 = rand() % commWorldSize;
        }
        int h12 = rand() % 2;
        if (h12 == 0) {
          rankCode1[t1] = 0;
          rankCode1[i] = 1;
        }
        else {
          rankCode2[t1] = 0;
          rankCode2[i] = 1;
        }
      }
    }

    int ind = 0;
    for (int i = 0 ; i < commWorldSize ; i++) {
      if (rankCode2[i] || rankCode1[i]) {
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
    nbPartSeg[0] = 2;
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
    nbPartSeg[1] = 2;
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

  // Coupling creation
  printf("%d - Create coupling\n", rank);
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1_part";
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Cpl_create(codeName[i_code],                                      // Code name
                   cpl_name,                                              // Coupling id
                   codeCoupledName[i_code],                               // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,                                // Coupling type
                   CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                   nbPart[i_code],                                        // Partition number
                   CWP_DYNAMIC_MESH_STATIC,                               // Mesh type
                   CWP_TIME_EXCH_USER_CONTROLLED);                           // Postprocessing frequency
  }

  printf("%d - Set visu\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Visu_set(codeName[i_code],        // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

  // Mesh definition
  printf("%d - Create mesh\n", rank);

  int     **nVertex           = NULL;
  double ***coords            = NULL;
  int     **nElts             = NULL;
  int    ***eltsConnecPointer = NULL;
  int    ***eltsConnec        = NULL;

  // Domain bounds
  const double xmin = -10;
  const double xmax =  10;
  const double ymin = -10;
  const double ymax =  10;

  coords            = (double ***) malloc(sizeof(double **) * n_code);
  eltsConnecPointer = (int    ***) malloc(sizeof(int    **) * n_code);
  eltsConnec        = (int    ***) malloc(sizeof(int    **) * n_code);
  nVertex           = (int     **) malloc(sizeof(int    * ) * n_code);
  nElts             = (int     **) malloc(sizeof(int    * ) * n_code);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    coords           [i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    eltsConnecPointer[i_code] = (int    **) malloc(sizeof(int    *) * nbPart[i_code]);
    eltsConnec       [i_code] = (int    **) malloc(sizeof(int    *) * nbPart[i_code]);
    nVertex          [i_code] = (int     *) malloc(sizeof(int     ) * nbPart[i_code]);
    nElts            [i_code] = (int     *) malloc(sizeof(int     ) * nbPart[i_code]);
  }

  srand(time(NULL));

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    randLevel = 0.4;
    double xSegPart = (xmax - xmin) / (double) nbPartSeg[i_code];
    double ySegPart = (ymax - ymin) / (double) nbPartSeg[i_code];
    for (int u = 0 ; u < nbPartSeg[i_code] ; u++) {
      for (int v = 0 ; v < nbPartSeg[i_code] ; v++) {
        int i_part = nbPartSeg[i_code] * u + v;
        int nVertexSegPart = nVertexSeg / nbPartSeg[i_code];
        nVertex[i_code][i_part] =  nVertexSegPart      *  nVertexSegPart;
        nElts  [i_code][i_part] = (nVertexSegPart - 1) * (nVertexSegPart - 1);

        coords           [i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nVertex[i_code][i_part]   );
        eltsConnecPointer[i_code][i_part] = (int    *) malloc(sizeof(int   ) *    (nElts[i_code][i_part] + 1));
        eltsConnec       [i_code][i_part] = (int    *) malloc(sizeof(int   ) * 4 * nElts[i_code][i_part]     );
        if (!strcmp(codeName[i_code], "code1")) {
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
        if (!strcmp(codeName[i_code], "code2")) {
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

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Mesh_interf_finalize(codeName[i_code], cpl_name);
  }

  printf("%d - Exchange Code1 <-> Code2\n", rank);

  double ***sendValues     = (double ***) malloc(sizeof(double **) * n_code);
  double ***sendValues2    = (double ***) malloc(sizeof(double **) * n_code);
  double ***sendValues3    = (double ***) malloc(sizeof(double **) * n_code);
  double ***recvValues     = (double ***) malloc(sizeof(double **) * n_code);
  double ***recvValues2    = (double ***) malloc(sizeof(double **) * n_code);
  double ***recvValues3    = (double ***) malloc(sizeof(double **) * n_code);
  double ***Values2Vertex  = (double ***) malloc(sizeof(double **) * n_code);
  double ***recvValuesUser = (double ***) malloc(sizeof(double **) * n_code);

  int nbPoints = 3;
  int **nbPointsUser = (int **) malloc(sizeof(int *) * n_code);

  double ***coordsPointsUser = (double ***) malloc(sizeof(double **) * n_code);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    sendValues [i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    sendValues2[i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    sendValues3[i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);

    nbPointsUser    [i_code] = (int     *) malloc(sizeof(int     ) * nbPart[i_code]);
    recvValues      [i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    recvValues2     [i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    recvValues3     [i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    Values2Vertex   [i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    recvValuesUser  [i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    coordsPointsUser[i_code] = (double **) malloc(sizeof(double *) * nbPart[i_code]);
    for (int i_part = 0 ; i_part < nbPart[i_code] ; i_part++) {

      if (!strcmp(codeName[i_code], "code1")) {
        sendValues [i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);
        sendValues2[i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);
        sendValues3[i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);

        recvValues      [i_code][i_part] = (double *) malloc(sizeof(double) * nElts[i_code][i_part]);
        Values2Vertex   [i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);
        nbPointsUser    [i_code][i_part] = 0;
        recvValuesUser  [i_code][i_part] = (double *) malloc(sizeof(double) * nbPointsUser[i_code][i_part]);
        coordsPointsUser[i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nbPointsUser[i_code][i_part]);
        for (int i = 0 ; i < nVertex[i_code][i_part] ; i++) {
          sendValues[i_code][i_part][i] = coords[i_code][i_part][3 * i];
        }
      }
      else {
        sendValues [i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nElts[i_code][i_part]);
        sendValues2[i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nElts[i_code][i_part]);
        sendValues3[i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nElts[i_code][i_part]);

        recvValues      [i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nVertex[i_code][i_part]);
        recvValues2     [i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nVertex[i_code][i_part]);
        recvValues3     [i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nVertex[i_code][i_part]);
        Values2Vertex   [i_code][i_part] = (double *) malloc(sizeof(double) *     nVertex[i_code][i_part]);
        nbPointsUser    [i_code][i_part] = nbPoints;
        recvValuesUser  [i_code][i_part] = (double *) malloc(sizeof(double) *     nbPointsUser[i_code][i_part]);
        coordsPointsUser[i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nbPointsUser[i_code][i_part]);
        for (int i = 0 ; i < nElts[i_code][i_part] ; i++) {
          sendValues [i_code][i_part][3 * i] = rank;
          sendValues2[i_code][i_part][3 * i] = i_part;
          sendValues3[i_code][i_part][3 * i] = i;
        }
        for (int i = 0 ; i < nVertex[i_code][i_part] ; i++) {
          Values2Vertex[i_code][i_part][i] = i;//coords[3 * i];
        }

        for (int i = 0 ; i < nbPointsUser[i_code][i_part] ; i++) {
          coordsPointsUser[i_code][i_part][3 * i    ] = 0.0 + 0.01 * i_part + 0.001 * i;
          coordsPointsUser[i_code][i_part][3 * i + 1] = rank * 0.1 + 0.01 * i_part + 0.001 * i;
          coordsPointsUser[i_code][i_part][3 * i + 2] = 0.0;
        }
      }
    }
  }

  // Exchange
  int code1I = 1;
  int code2I = 1;
  int code3I = 1;
  int code4I = 1;
  int code5I = 1;
  int code6I = 1;
  int code7I = 1;

  const char *fieldName1 = "cooX";;
  const char *fieldName2 = "rank";;
  const char *fieldName3 = "userField";;
  const char *fieldName4 = "userInterpolation";
  const char *fieldName5 = "sendRecvField";
  const char *fieldName6 = "part";
  const char *fieldName7 = "num";

  printf("%d - Defining fields\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (!strcmp(codeName[i_code], "code1")) {
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
      if (code6I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName6,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_CELL_CENTER,
                         CWP_FIELD_EXCH_SEND,
                         CWP_STATUS_ON);
      }
      if (code7I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName7,
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
        if (code6I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName6,
                             i_part,
                             CWP_FIELD_MAP_SOURCE,
                             sendValues2[i_code][i_part]);
        }
        if (code7I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName7,
                             i_part,
                             CWP_FIELD_MAP_SOURCE,
                             sendValues3[i_code][i_part]);
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
        CWP_Interp_function_set(codeName[i_code], cpl_name, fieldName4, _userInterpolation);
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

      if (code6I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName6,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         1,
                         CWP_DOF_LOCATION_NODE,
                         CWP_FIELD_EXCH_RECV,
                         CWP_STATUS_ON);
      }

      if (code7I == 1) {
        CWP_Field_create(codeName[i_code],
                         cpl_name,
                         fieldName7,
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
                             recvValues[i_code][i_part]);
        }
        if (code6I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName6,
                             i_part,
                             CWP_FIELD_MAP_TARGET,
                             recvValues2[i_code][i_part]);
        }
        if (code7I == 1) {
          CWP_Field_data_set(codeName[i_code],
                             cpl_name,
                             fieldName7,
                             i_part,
                             CWP_FIELD_MAP_TARGET,
                             recvValues3[i_code][i_part]);
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
          CWP_User_tgt_pts_set(codeName[i_code],
                               cpl_name,
                               i_part,
                               nbPointsUser[i_code][i_part],
                               coordsPointsUser[i_code][i_part],
                               NULL);
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
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  printf("%d - Before compute\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Spatial_interp_weights_compute(codeName[i_code], cpl_name);
  }
  printf("%d - After compute\n", rank);

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (!strcmp(codeName[i_code], "code1")) {
      if (code1I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName1);
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName1);
      }

      if (code6I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName6);
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName6);
      }

      if (code7I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName7);
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName7);
      }

    }
    else {
      if (code1I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName1);
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName1);
      }

      if (code6I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName6);
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName6);
      }

      if (code7I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName7);
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName7);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (!strcmp(codeName[i_code], "code1")) {
      if (code2I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName2);
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName2);
      }
    }
    else {
      if (code2I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName2);
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName2);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (!strcmp(codeName[i_code], "code1")) {
      if (code3I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName3);
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName3);
      }
    }
    else {
      if (code3I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName3);
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName3);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (!strcmp(codeName[i_code], "code1")) {
      if (code4I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName4);
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName4);
      }
    }
    else {
      if (code4I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName4);
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName4);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (!strcmp(codeName[i_code], "code1")) {
      if (code5I == 1) {
        CWP_Field_irecv(codeName[i_code], cpl_name, fieldName5);
        CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName5);
      }
    }
    else {
      if (code5I == 1) {
        CWP_Field_issend(codeName[i_code], cpl_name, fieldName5);
        CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName5);
      }
    }
  }

  printf("%d - Delete mesh\n", rank);
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Mesh_interf_del(codeName[i_code], cpl_name);
  }

  printf("%d - Delete coupling\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Cpl_del(codeName[i_code], cpl_name);
  }

  // Freeing memory
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    for (int i_part = 0 ; i_part < nbPart[i_code] ; i_part++) {
      free(sendValues   [i_code][i_part]);
      free(recvValues   [i_code][i_part]);
      free(sendValues2  [i_code][i_part]);
      if (strcmp(codeName[i_code], "code1")) {
        free(recvValues2  [i_code][i_part]);
        free(recvValues3  [i_code][i_part]);
      }
      free(sendValues3  [i_code][i_part]);
      free(Values2Vertex[i_code][i_part]);
      free(recvValuesUser[i_code][i_part]);

      free(coords           [i_code][i_part]);
      free(eltsConnecPointer[i_code][i_part]);
      free(eltsConnec       [i_code][i_part]);

      free(coordsPointsUser[i_code][i_part]);
    }
    free(sendValues   [i_code]);
    free(recvValues   [i_code]);
    free(sendValues2  [i_code]);
    free(recvValues2  [i_code]);
    free(sendValues3  [i_code]);
    free(recvValues3  [i_code]);
    free(Values2Vertex[i_code]);
    free(recvValuesUser[i_code]);

    free(coords           [i_code]);
    free(eltsConnecPointer[i_code]);
    free(eltsConnec       [i_code]);
    free(nVertex          [i_code]);
    free(nElts            [i_code]);

    free(nbPointsUser    [i_code]);
    free(coordsPointsUser[i_code]);
  }

  free(sendValues);
  free(recvValues);
  free(sendValues2);
  free(recvValues2);
  free(sendValues3);
  free(recvValues3);
  free(Values2Vertex);
  free(recvValuesUser);

  free(coords);
  free(eltsConnecPointer);
  free(eltsConnec);
  free(nVertex);
  free(nElts);

  free(times_init);
  free(codeName);
  free(codeCoupledName);
  free(is_coupled_rank);
  free(nbPartSeg);
  free(nbPart);
  free(localComm);

  free(nbPointsUser);
  free(coordsPointsUser);

  // Finalize
  CWP_Finalize();
  MPI_Finalize();
  return EXIT_SUCCESS;
}
