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
#include <assert.h>

#include "cwp.h"

/*----------------------------------------------------------------------
 *
 * Main : linear coupling test
 *
 *---------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
  FILE *outputFile;

  MPI_Init(&argc, &argv);

  int rank, comm_world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  char *srcName = (char *) malloc(sizeof(char) * (strlen(__FILE__) + 1));
  strcpy(srcName, __FILE__);
  char *srcBaseName = NULL;
  srcBaseName = strrchr(srcName, '.');
  if (srcBaseName != NULL) *srcBaseName = '\0';
  srcBaseName = NULL;
  srcBaseName = strrchr(srcName, '/');
  if (srcBaseName != NULL) srcBaseName += 1;
  else srcBaseName = srcName;

  if (rank == 0) printf("\nSTART: %s\n", srcBaseName);

  // Initialization
  int n_code_name = 0;
  const char **codeNames = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;

  if (rank == 0) {
    n_code_name = 2;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code1";
    codeNames[1] = "code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
  }
  else if (rank == 1) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }
  else if (rank == 2) {
    n_code_name = 4;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code1";
    codeNames[1] = "code2";
    codeNames[2] = "code3";
    codeNames[3] = "code4";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
    is_coupled_rank[2] = CWP_STATUS_ON;
    is_coupled_rank[3] = CWP_STATUS_ON;
  }
  else if (rank == 3) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code3";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }
  else if (rank == 4) {
    n_code_name = 2;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code3";
    codeNames[1] = "code4";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
  }
  else if (rank == 5) {
    n_code_name = 2;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code1";
    codeNames[1] = "code3";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
  }
  else if (rank == 6) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }
  else if (rank == 7) {
    n_code_name = 3;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code1";
    codeNames[1] = "code2";
    codeNames[2] = "code3";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
    is_coupled_rank[2] = CWP_STATUS_ON;
  }
  else if (rank == 8) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code4";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }
  else if (rank == 9) {
    n_code_name = 2;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code2";
    codeNames[1] = "code3";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
  }

  const char *fileName = NULL;
  if (rank == 0) fileName = "c_new_api_0000.txt";
  else if (rank == 1) fileName = "c_new_api_0001.txt";
  else if (rank == 2) fileName = "c_new_api_0002.txt";
  else if (rank == 3) fileName = "c_new_api_0003.txt";
  else if (rank == 4) fileName = "c_new_api_0004.txt";
  else if (rank == 5) fileName = "c_new_api_0005.txt";
  else if (rank == 6) fileName = "c_new_api_0006.txt";
  else if (rank == 7) fileName = "c_new_api_0007.txt";
  else if (rank == 8) fileName = "c_new_api_0008.txt";
  else if (rank == 9) fileName = "c_new_api_0009.txt";
  outputFile = fopen(fileName, "w");

  times_init = malloc(sizeof(double) * n_code_name);

  for (int i = 0 ; i < n_code_name ; i++) times_init[i] = 0;

  MPI_Comm *localComm = malloc(sizeof(MPI_Comm) * n_code_name);
  CWP_Init(MPI_COMM_WORLD, n_code_name, (const char **) codeNames, is_coupled_rank, times_init, localComm);

  // Output redirection
  int currentRank;
  int localCommSize;

  for (int i = 0 ; i < n_code_name ; i++) {
    MPI_Comm_rank(localComm[i], &currentRank);
    MPI_Comm_size(localComm[i], &localCommSize);
  }

  // Finalize
  if (rank == 0 || rank == 1 || rank == 2 || rank == 5 || rank == 7) {
    int toto = 111;
    CWP_Param_lock("code1");
    //  MPI_Barrier (MPI_COMM_WORLD);
    CWP_Param_add("code1", "toto", CWP_INT, &toto);
    const char *A = "Bonjour !";
    CWP_Param_add("code1", "toto2", CWP_CHAR, &A);
    CWP_Param_unlock("code1");
  }
  MPI_Barrier(MPI_COMM_WORLD);

  int titi;
  CWP_Param_get("code1", "toto", CWP_INT, &titi);

  char *titi2;
  CWP_Param_get("code1", "toto2", CWP_CHAR, &titi2);

  free(titi2);
  assert(titi == 111);

  char cpl_id1[] = "cpl1_code1_code2";
  char cpl_id2[] = "cpl2_code1_code3";
  char cpl_id3[] = "cpl3_code2_code3";
  char cpl_id4[] = "cpl4_code4_code3";
  char cpl_id5[] = "cpl5_code1_code4";
  char cpl_id6[] = "cpl6_code2_code4";

  CWP_Spatial_interp_t interp_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;

  // cpl1
  if (rank == 0 || rank == 1 || rank == 2 || rank == 5 || rank == 7) {
    int v = -1;
    if (rank == 0) {
      v = 11;
    }
    MPI_Bcast(&v, 1, MPI_INT, 0, localComm[0]);
    printf("code 1 v : %d\n", v);
    fflush(stdout);
    CWP_Cpl_create("code1", cpl_id1, "code2", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
  }
  if (rank == 0 || rank == 2 || rank == 6 || rank == 7 || rank == 9) {
    CWP_Cpl_create("code2", cpl_id1, "code1", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
    int v = -2;
    if (rank == 0) {
      v = 21;
    }
    MPI_Comm intraComm;
    if (rank == 0 || rank == 2 || rank == 2 || rank == 7) {
      intraComm = localComm[1];
    }
    else if (rank == 6 || rank == 9) {
      intraComm = localComm[0];
    }
    MPI_Bcast(&v, 1, MPI_INT, 0, intraComm);
    printf("code 2 v : %d\n", v);
    fflush(stdout);
  }
  printf("Fin premier couplage\n");
  MPI_Barrier(MPI_COMM_WORLD);

  // cpl2
  if (rank == 0 || rank == 1 || rank == 2 || rank == 5 || rank == 7) {
    CWP_Cpl_create("code1", cpl_id2, "code3", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
  }
  if (rank == 2 || rank == 3 || rank == 4 || rank == 5 || rank == 7 || rank == 9) {
    CWP_Cpl_create("code3", cpl_id2, "code1", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
    int v = -2;
    if (rank == 2) {
      v = 31;
    }
    MPI_Comm intraComm;
    if (rank == 2 || rank == 7) {
      intraComm = localComm[2];
    }
    else if (rank == 3 || rank == 4) {
      intraComm = localComm[0];
    }
    else if (rank == 5 || rank == 9) {
      intraComm = localComm[1];
    }
    MPI_Bcast(&v, 1, MPI_INT, 0, intraComm);
    printf("code 3 v : %d\n", v);
    fflush(stdout);
  }

  // cpl3
  if (rank == 0 || rank == 2 || rank == 6 || rank == 7 || rank == 9) {
    CWP_Cpl_create("code2", cpl_id3, "code3", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
  }
  if (rank == 2 || rank == 3 || rank == 4 || rank == 5 || rank == 7 || rank == 9) {
    CWP_Cpl_create("code3", cpl_id3, "code2", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
  }

  // cpl4
  if (rank == 2 || rank == 4 || rank == 8) {
    CWP_Cpl_create("code4", cpl_id4, "code3", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
    int v = -2;
    if (rank == 2) {
      v = 41;
    }
    MPI_Comm intraComm;
    if (rank == 2) {
      intraComm = localComm[3];
    }
    else if (rank == 4) {
      intraComm = localComm[1];
    }
    else if (rank == 8) {
      intraComm = localComm[0];
    }
    MPI_Bcast(&v, 1, MPI_INT, 0, intraComm);
    printf("code 4 v : %d\n", v);
    fflush(stdout);
  }

  if (rank == 2 || rank == 3 || rank == 4 || rank == 5 || rank == 7 || rank == 9) {
    CWP_Cpl_create("code3", cpl_id4, "code4", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
  }

  // cpl5
  if (rank == 0 || rank == 1 || rank == 2 || rank == 5 || rank == 7) {
    CWP_Cpl_create("code1", cpl_id5, "code4", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
  }
  if (rank == 2 || rank == 4 || rank == 8) {
    CWP_Cpl_create("code4", cpl_id5, "code1", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
  }

  // cpl6
  if (rank == 0 || rank == 2 || rank == 6 || rank == 7 || rank == 9) {
    CWP_Cpl_create("code2", cpl_id6, "code4", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
  }
  if (rank == 2 || rank == 4 || rank == 8) {
    CWP_Cpl_create("code4", cpl_id6, "code2", CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
  }

  printf("All done for rank %d\n", rank);

  CWP_Finalize();
  MPI_Finalize();

  free(srcName);
  free(localComm);
  free(codeNames);
  free(is_coupled_rank);
  free(times_init);
  fclose(outputFile);

  return 0;
}
