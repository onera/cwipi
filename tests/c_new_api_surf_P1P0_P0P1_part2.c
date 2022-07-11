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
#include <assert.h>
#include <math.h>
#include <time.h>

#include "cwp.h"
#include "cwp_priv.h"
#include "pdm_timer.h"


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
         "  -nx     <val>    Global number of vertices in the side of domaine.\n\n"
         "  -part   <val>    Part of active ranks in the coupling.\n\n"
         "  -s      <val>    Size of domain.\n\n"
         "  -h               this message.\n\n");
  exit(exit_code);
}


/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nx             --> Global number of vertices in the side of domaine
 *   part           --> Part of active ranks in the coupling
 *   s              --> Size of domain
 *
 *---------------------------------------------------------------------*/

static void
_read_args(int argc, char **argv, CWP_g_num_t *nx, double *part, double *s, int *n_compute) {
  int i = 1;

  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {_display_usage(EXIT_SUCCESS);}
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc + 1) {_display_usage(EXIT_FAILURE);}
      else {*nx = atoi(argv[i]);}
    }
    else if (strcmp(argv[i], "-part") == 0) {
      i++;
      if (i >= argc) {_display_usage(EXIT_FAILURE);}
      else {*part = atof(argv[i]);}
    }
    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc) {_display_usage(EXIT_FAILURE);}
      else {*s = atof(argv[i]);}
    }
    else if (strcmp(argv[i], "-nc") == 0) {
      i++;
      if (i >= argc) {_display_usage(EXIT_FAILURE);}
      else {*n_compute = atoi(argv[i]);}
    }
    i++;
  }
}


static int
_random_global_int(void) {
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int resultInt = -1;
  int *randInt = NULL;
  if (rank == 0) {
    srand(time(NULL));
    randInt = (int *) malloc(sizeof(int) * size);
    for (int i = 0 ; i < size ; i++) {
      int test;
      int ok;
      do {
        test = rand() % size;
        ok = 1;
        for (int j = 0 ; j < i ; j++) {
          if (test == randInt[j]) {
            ok = 0;
            break;
          }
        }
      }
      while (!ok);
      randInt[i] = test;
    }
  }
  MPI_Scatter(randInt, 1, MPI_INT, &resultInt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return resultInt;
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

  // Read args from command line
  CWP_g_num_t nx = 10;
  double part = 1.;
  double s = 10.;
  int n_compute = 1;
  const double dev_limit = 0.05;

  _read_args(argc, argv, &nx, &part, &s, &n_compute);

  // Coupled MPI process fraction for each code
  double *cpl_frac = (double *) malloc(sizeof(double) * 2);
  cpl_frac[0] = part;
  cpl_frac[1] = part;

  // Non null interface mesh fraction for each code
  double *non_null_mesh_frac = (double *) malloc(sizeof(double) * 2);
  non_null_mesh_frac[0] = 1.0;
  non_null_mesh_frac[1] = 1.0;

  // Init + create coupling
  assert (commWorldSize >= 2);

  int n_code;
  double prop1 = 1. / 2.;
  double prop2 = 1. / 2.;
  double prop12 = 0. / 3.;

  int partial_covering_mpi_domains = 1;
  if (partial_covering_mpi_domains == 1) {
    prop1 = 1. / 3.;
    prop2 = 1. / 3.;
    prop12 = 1. / 3.;
  }

  int randomGlobalInt = _random_global_int();

  if (randomGlobalInt < (double) commWorldSize * prop1) {n_code = 1;}
  else if (randomGlobalInt < (double) commWorldSize * (prop1 + prop12)) {n_code = 2;}
  else {n_code = 1;}

  CWP_g_num_t *nxCode = (CWP_g_num_t *) malloc(n_code * sizeof(CWP_g_num_t));

  const char **codeName = (const char **) malloc(sizeof(const char *) * n_code);
  const char **codeCoupledName = (const char **) malloc(sizeof(const char *) * n_code);
  int *codeId = (int *) malloc(sizeof(int) * n_code);
  MPI_Comm *localComm = malloc(sizeof(MPI_Comm) * n_code);
  int *localRank = (int *) malloc(sizeof(int) * n_code);
  CWP_Status_t *is_coupled_rank = (CWP_Status_t *) malloc(sizeof(CWP_Status_t) * n_code);
  double *time_init = (double *) malloc(sizeof(double) * n_code);
  int *nb_part = (int *) malloc(sizeof(int) * n_code);

  if (randomGlobalInt < (double) commWorldSize * prop1) {
    codeName[0] = "code1";
    codeId[0] = 1;
    codeCoupledName[0] = "code2";

    if (randomGlobalInt < (double) commWorldSize * prop1 * cpl_frac[0]) {is_coupled_rank[0] = CWP_STATUS_ON;}
    else {is_coupled_rank[0] = CWP_STATUS_OFF;}

    time_init[0] = 0.0;
    nxCode[0] = nx;
    nb_part[0] = 2;
  }
  else if (randomGlobalInt < (double) commWorldSize * (prop1 + prop12)) {
    codeName[0] = "code1";
    codeId[0] = 1;
    codeCoupledName[0] = "code2";

    if (randomGlobalInt < (double) commWorldSize * (prop1 + cpl_frac[0] * prop12)) {is_coupled_rank[0] = CWP_STATUS_ON;}
    else {is_coupled_rank[0] = CWP_STATUS_OFF;}

    time_init[0] = 0.0;
    nxCode[0] = nx;
    nb_part[0] = 2;

    codeName[1] = "code2";
    codeId[1] = 2;
    codeCoupledName[1] = "code1";
    if (randomGlobalInt < (double) commWorldSize * (prop1 + cpl_frac[1] * prop12)) {is_coupled_rank[1] = CWP_STATUS_ON;}
    else {is_coupled_rank[1] = CWP_STATUS_OFF;}

    time_init[1] = 0.0;
    nxCode[1] = (CWP_g_num_t) (2.0 * (double) nx);
    nb_part[1] = 4;
  }
  else if (randomGlobalInt < (double) commWorldSize * (prop1 + prop12 + prop2)) {
    codeName[0] = "code2";
    codeId[0] = 2;
    codeCoupledName[0] = "code1";

    if (randomGlobalInt < (double) commWorldSize * (prop1 + prop12 + cpl_frac[1] * prop2)) {is_coupled_rank[0] = CWP_STATUS_ON;}
    else {is_coupled_rank[0] = CWP_STATUS_OFF;}

    time_init[0] = 0.0;
    nxCode[0] = (CWP_g_num_t) (2.0 * (double) nx);
    nb_part[0] = 4;
  }

  CWP_Init(MPI_COMM_WORLD, n_code, codeName, is_coupled_rank, time_init, localComm);

  printf("%d - Create coupling\n", rank);
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1_part2";
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (is_coupled_rank[i_code] == CWP_STATUS_ON) {
      MPI_Comm_rank(localComm[i_code], &(localRank[i_code]));
      assert(localComm[i_code] != MPI_COMM_NULL);

      CWP_Cpl_create(codeName[i_code],                                      // Code name
                     cpl_name,                                              // Coupling id
                     codeCoupledName[i_code],                               // Coupled application id
                     CWP_INTERFACE_SURFACE,
                     CWP_COMM_PAR_WITH_PART,                                // Coupling type
                     CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                     nb_part[i_code],                                       // Partition number
                     CWP_DYNAMIC_MESH_STATIC,                               // Mesh type
                     CWP_TIME_EXCH_USER_CONTROLLED);                          // Postprocessing frequency
    }
  }

  PDM_timer_t *timer = PDM_timer_create();
  PDM_timer_t *timer2 = PDM_timer_create();
  PDM_timer_init(timer);

  int n_int = 1;
  double compute_time[n_compute];
  double compute_exch_time[n_int];

  double ***sendValues = (double ***) malloc(sizeof(double **) * n_code);
  double ***recvValues = (double ***) malloc(sizeof(double **) * n_code);

  MPI_Barrier(MPI_COMM_WORLD);

  printf("%d - Set visu\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (is_coupled_rank[i_code] == CWP_STATUS_ON) {
      CWP_Visu_set(codeName[i_code],        // Code name
                   cpl_name,                // Coupling id
                   1,                       // Postprocessing frequency
                   CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                   "text");                 // Postprocessing option
      MPI_Comm connectableComm = CWP_Connectable_comm_get((char*) codeName[i_code]);

      CWP_surf_gen_init((char*) codeName[i_code],
                        (int) nxCode[i_code],
                        (int) nxCode[i_code],
                        nb_part[i_code],
                        &connectableComm,
                        non_null_mesh_frac[i_code],
                        s,
                        (double) codeId[i_code]);
    }
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (is_coupled_rank[i_code] == CWP_STATUS_ON) {
      CWP_surf_gen_compute((char*) codeName[i_code]);
    }
  }

  int **eltsConnecPolyIndex = NULL;
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (is_coupled_rank[i_code] == CWP_STATUS_ON) {
      sendValues[i_code] = (double **) malloc(sizeof(double *) * nb_part[i_code]);
      recvValues[i_code] = (double **) malloc(sizeof(double *) * nb_part[i_code]);

      double **coords = (double **) malloc(sizeof(double *) * nb_part[i_code]);

      CWP_g_num_t **vtxGnum = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * nb_part[i_code]);
      int *nVtx = (int *) malloc(sizeof(int) * nb_part[i_code]);

      int *n_tri = (int *) malloc(sizeof(int) * nb_part[i_code]);
      int **eltsConnecTri = (int **) malloc(sizeof(int *) * nb_part[i_code]);
      CWP_g_num_t **eltsGnumTri = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * nb_part[i_code]);

      int *n_quad = (int *) malloc(sizeof(int) * nb_part[i_code]);
      int **eltsConnecQuad = (int **) malloc(sizeof(int *) * nb_part[i_code]);
      CWP_g_num_t **eltsGnumQuad = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * nb_part[i_code]);

      int *n_poly2d = (int *) malloc(sizeof(int) * nb_part[i_code]);
      eltsConnecPolyIndex = (int **) malloc(sizeof(int *) * nb_part[i_code]);
      int **eltsConnecPoly = (int **) malloc(sizeof(int *) * nb_part[i_code]);

      int *n_faces = (int *) malloc(sizeof(int) * nb_part[i_code]);
      int **faceEdgeIdx = (int **) malloc(sizeof(int *) * nb_part[i_code]);
      int **faceEdge = (int **) malloc(sizeof(int *) * nb_part[i_code]);

      int *n_edges = (int *) malloc(sizeof(int) * nb_part[i_code]);
      int **edgeVtxIdx = (int **) malloc(sizeof(int *) * nb_part[i_code]);
      int **edgeVtx = (int **) malloc(sizeof(int *) * nb_part[i_code]);

      CWP_g_num_t **faceLNToGN = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * nb_part[i_code]);
      CWP_g_num_t **eltsGnumPoly = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * nb_part[i_code]);

      int gnum_compute = 0;
      int face_edge = 1;

      for (int i_part = 0 ; i_part < nb_part[i_code] ; i_part++) {
        nVtx[i_part] = 0;
        coords[i_part] = NULL;
        vtxGnum[i_part] = NULL;

        int nElts = 0;

        n_tri[i_part] = 0;
        eltsConnecTri[i_part] = NULL;
        eltsGnumTri[i_part] = NULL;

        n_quad[i_part] = 0;
        eltsConnecQuad[i_part] = NULL;
        eltsGnumQuad[i_part] = NULL;

        n_poly2d[i_part] = 0;
        eltsConnecPolyIndex[i_part] = NULL;
        eltsConnecPoly[i_part] = NULL;
        eltsGnumPoly[i_part] = NULL;

        CWP_surf_gen_by_block_get((char*) codeName[i_code],
                                  i_part,
                                  &nVtx[i_part],
                                  &coords[i_part],
                                  &vtxGnum[i_part],
                                  &nElts,
                                  &n_tri[i_part],
                                  &eltsConnecTri[i_part],
                                  &eltsGnumTri[i_part],
                                  &n_quad[i_part],
                                  &eltsConnecQuad[i_part],
                                  &eltsGnumQuad[i_part],
                                  &n_poly2d[i_part],
                                  &eltsConnecPolyIndex[i_part],
                                  &eltsConnecPoly[i_part],
                                  &eltsGnumPoly[i_part]);

        if (face_edge == 1) {
          CWP_surf_face_edge_get((char*) codeName[i_code],
                                 i_part,
                                 &nVtx[i_part],
                                 &coords[i_part],
                                 &vtxGnum[i_part],
                                 &n_faces[i_part],
                                 &faceEdgeIdx[i_part],
                                 &faceEdge[i_part],
                                 &n_edges[i_part],
                                 &edgeVtxIdx[i_part],
                                 &edgeVtx[i_part],
                                 &faceLNToGN[i_part]);
        }

        if (gnum_compute == 0) {
          vtxGnum[i_part] = NULL;
          eltsGnumTri[i_part] = NULL;
          eltsGnumQuad[i_part] = NULL;
          eltsGnumPoly[i_part] = NULL;
        }

        CWP_Mesh_interf_vtx_set(codeName[i_code],
                                cpl_name,
                                i_part,
                                nVtx[i_part],
                                coords[i_part],
                                vtxGnum[i_part]);

        if (face_edge == 1) {
          CWP_Mesh_interf_from_faceedge_set(codeName[i_code],
                                            cpl_name,
                                            i_part,
                                            n_faces[i_part],
                                            faceEdgeIdx[i_part],
                                            faceEdge[i_part],
                                            n_edges[i_part],
                                            edgeVtxIdx[i_part],
                                            edgeVtx[i_part],
                                            faceLNToGN[i_part]);
        }
        else {
          int block_id_tri = CWP_Mesh_interf_block_add(codeName[i_code],
                                                       cpl_name,
                                                       CWP_BLOCK_FACE_TRIA3);

          CWP_Mesh_interf_block_std_set(codeName[i_code],
                                        cpl_name,
                                        i_part,
                                        block_id_tri,
                                        n_tri[i_part],
                                        eltsConnecTri[i_part],
                                        eltsGnumTri[i_part]);

          int block_id_quad = CWP_Mesh_interf_block_add(codeName[i_code],
                                                        cpl_name,
                                                        CWP_BLOCK_FACE_QUAD4);

          CWP_Mesh_interf_block_std_set(codeName[i_code],
                                        cpl_name,
                                        i_part,
                                        block_id_quad,
                                        n_quad[i_part],
                                        eltsConnecQuad[i_part],
                                        eltsGnumQuad[i_part]);

          int block_id_poly = CWP_Mesh_interf_block_add(codeName[i_code],
                                                        cpl_name,
                                                        CWP_BLOCK_FACE_POLY);

          CWP_Mesh_interf_f_poly_block_set(codeName[i_code],
                                           cpl_name,
                                           i_part,
                                           block_id_poly,
                                           n_poly2d[i_part],
                                           eltsConnecPolyIndex[i_part],
                                           eltsConnecPoly[i_part],
                                           eltsGnumPoly[i_part]);
        }

        if (codeId[i_code] == 1) {
          sendValues[i_code][i_part] = (double *) malloc(sizeof(double) * nVtx[i_part]);
          recvValues[i_code][i_part] = (double *) malloc(sizeof(double) * nElts);
          for (int i = 0 ; i < nVtx[i_part] ; i++) {
            sendValues[i_code][i_part][i] = coords[i_part][3 * i];
          }
        }
        else {
          double *field_tri;
          double *field_quad;
          double *field_poly2d;
          CWP_surf_gen_tri_field_get((char*) codeName[i_code], i_part, &field_tri);
          CWP_surf_gen_quad_field_get((char*) codeName[i_code], i_part, &field_quad);
          CWP_surf_gen_poly_field_get((char*) codeName[i_code], i_part, &field_poly2d);
          sendValues[i_code][i_part] = (double *) malloc(sizeof(double) * nElts);
          recvValues[i_code][i_part] = (double *) malloc(sizeof(double) * nVtx[i_part]);

          int ind = 0;
          for (int i = 0 ; i < n_tri[i_part] ; i++) {
            sendValues[i_code][i_part][ind] = field_tri[i];
            ind++;
          }
          for (int i = 0 ; i < n_quad[i_part] ; i++) {
            sendValues[i_code][i_part][ind] = field_quad[i];
            ind++;
          }
          for (int i = 0 ; i < n_poly2d[i_part] ; i++) {
            sendValues[i_code][i_part][ind] = field_poly2d[i];
            ind++;
          }
        }
      }
    }
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (is_coupled_rank[i_code] == CWP_STATUS_ON) {
      CWP_Mesh_interf_finalize(codeName[i_code], cpl_name);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  printf("%d - Exchange Code1 <-> Code2\n", rank);

  // Exchange
  // field1: code1 -> code2
  // field2: code2 -> code1
  const char *fieldName1 = "cooX";
  const char *fieldName2 = "rank";

  CWP_Status_t visu_status = CWP_STATUS_ON;
  printf("%d - Defining fields\n", rank);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (strcmp(codeName[i_code], "code1") == 0 && is_coupled_rank[i_code] == CWP_STATUS_ON) {
      CWP_Field_create(codeName[i_code],
                       cpl_name,
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_create(codeName[i_code],
                       cpl_name,
                       fieldName2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      for (int i_part = 0 ; i_part < nb_part[i_code] ; i_part++) {
        CWP_Field_data_set(codeName[i_code],
                           cpl_name,
                           fieldName1,
                           i_part,
                           CWP_FIELD_MAP_SOURCE,
                           sendValues[i_code][i_part]);
        CWP_Field_data_set(codeName[i_code],
                           cpl_name,
                           fieldName2,
                           i_part,
                           CWP_FIELD_MAP_TARGET,
                           recvValues[i_code][i_part]);
      }
    }
    else if (strcmp(codeName[i_code], "code2") == 0 && is_coupled_rank[i_code] == CWP_STATUS_ON) {
      CWP_Field_create(codeName[i_code],
                       cpl_name,
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_create(codeName[i_code],
                       cpl_name,
                       fieldName2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      for (int i_part = 0 ; i_part < nb_part[i_code] ; i_part++) {
        CWP_Field_data_set(codeName[i_code],
                           cpl_name,
                           fieldName1,
                           i_part,
                           CWP_FIELD_MAP_TARGET,
                           recvValues[i_code][i_part]);
        CWP_Field_data_set(codeName[i_code],
                           cpl_name,
                           fieldName2,
                           i_part,
                           CWP_FIELD_MAP_SOURCE,
                           sendValues[i_code][i_part]);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  PDM_timer_init(timer);
  PDM_timer_resume(timer);
  PDM_timer_init(timer);
  PDM_timer_init(timer2);

  PDM_timer_resume(timer);
  double mean = 0.0;
  double mean2 = 0.0;
  double std_dev = 0.0;
  int n_it = 0;

  for (int i = 0 ; i < n_compute ; i++) {
    PDM_timer_init(timer2);
    PDM_timer_resume(timer2);
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      if (is_coupled_rank[i_code] == CWP_STATUS_ON) {
        CWP_Spatial_interp_weights_compute(codeName[i_code], cpl_name);
      }
    }

    PDM_timer_hang_on(timer2);
    compute_time[i] = PDM_timer_elapsed(timer2);

    mean += compute_time[i];
    MPI_Allreduce(&mean, &mean2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    mean2 = mean2 / ((double) (i + 1) * (double) commWorldSize);

    std_dev = 0.0;
    for (int h = 0 ; h <= i ; h++) {
      std_dev += pow((compute_time[h] - mean2) / mean2, 2);
    }
    std_dev = sqrt(std_dev) / (double) (i + 1);

    double std2;
    MPI_Allreduce(&std_dev, &std2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    std_dev = std2 / (double) commWorldSize;

    n_it = i;
    if (i > 3 && std_dev < dev_limit) {
      i = n_compute + 1;
    }

    if (localRank[0] == 0) {printf("Survey exchange %i %5.4e\n", i, std_dev);}
  }

  PDM_timer_hang_on(timer);

  MPI_Barrier(MPI_COMM_WORLD);

  if (localRank[0] == 0) {
    printf("New localization time %5.4e codeName[0] %s deviation %5.4e nb_it %i \n",
           mean2,
           codeName[0],
           std_dev,
           n_it);
  }

  PDM_timer_init(timer);
  PDM_timer_resume(timer);

  mean = 0.0;
  std_dev = 0.0;

  MPI_Barrier(MPI_COMM_WORLD);


  for (int i = 0 ; i < n_int ; i++) {
    PDM_timer_init(timer2);
    PDM_timer_resume(timer2);

    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      if (is_coupled_rank[i_code] == CWP_STATUS_ON) {
        if (strcmp(codeName[i_code], "code1") == 0) {
          CWP_Field_issend(codeName[i_code], cpl_name, fieldName1);
          CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName1);
          CWP_Field_irecv(codeName[i_code], cpl_name, fieldName2);
          CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName2);
        }
        else if (strcmp(codeName[i_code], "code2") == 0) {
          CWP_Field_irecv(codeName[i_code], cpl_name, fieldName1);
          CWP_Field_wait_irecv(codeName[i_code], cpl_name, fieldName1);
          CWP_Field_issend(codeName[i_code], cpl_name, fieldName2);
          CWP_Field_wait_issend(codeName[i_code], cpl_name, fieldName2);
        }
      }
    }

    PDM_timer_hang_on(timer2);
    compute_exch_time[i] = PDM_timer_elapsed(timer2);

    mean += compute_exch_time[i];
    MPI_Allreduce(&mean, &mean2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    mean2 = mean2 / ((double) (i + 1) * (double) commWorldSize);

    std_dev = 0.0;
    for (int h = 0 ; h <= i ; h++) {
      std_dev += pow((compute_exch_time[h] - mean2) / mean2, 2);
    }
    std_dev = sqrt(std_dev) / (double) (i + 1);

    double std2;
    MPI_Allreduce(&std_dev, &std2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    std_dev = std2 / (double) commWorldSize;

    if (i > 3 && std_dev < dev_limit) {
      i = n_int + 1;
    }
    if (localRank[0] == 0) {printf("Survey exchange %i %5.4e\n", i, std_dev);}
    MPI_Barrier(MPI_COMM_WORLD);
  }

  PDM_timer_hang_on(timer);

  if (localRank[0] == 0) {
    printf("New exchange time for %i iterations %5.4e s codeName[0] %s deviation %5.4e\n",
           n_int,
           mean2,
           codeName[0],
           std_dev);
  }


  MPI_Barrier(MPI_COMM_WORLD);

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (is_coupled_rank[i_code] == CWP_STATUS_ON) {
      CWP_Mesh_interf_del(codeName[i_code], cpl_name);
      CWP_Cpl_del(codeName[i_code], cpl_name);
      free(sendValues[i_code]);
      free(recvValues[i_code]);
    }
  }

  free(sendValues);
  free(recvValues);

  // Finalize
  CWP_Finalize();
  MPI_Finalize();
  return EXIT_SUCCESS;
}
