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
#include <stddef.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#include "cwp.h"

/*----------------------------------------------------------------------
 *                                                                     
 * Main : linear coupling test                                         
 *
 *---------------------------------------------------------------------*/
 
int main
(
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
)
{

  FILE *outputFile;

  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  char *srcName = (char *) malloc (sizeof(char) * (strlen(__FILE__) + 1));
  strcpy(srcName, __FILE__);
  char *srcBaseName = NULL;
  srcBaseName = strrchr(srcName, '.');
  if (srcBaseName != NULL)
    *srcBaseName = '\0';
  srcBaseName = NULL;
  srcBaseName = strrchr(srcName, '/');
  if (srcBaseName != NULL)
    srcBaseName += 1;
  else
    srcBaseName = srcName;

  if (rank == 0) 
    printf("\nSTART: %s\n", srcBaseName);


  /* Initialization
   * -------------- */

  int n_code_name = 0;
  char **codeNames = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;

  if (rank == 0) {
    n_code_name = 2;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] ="code1";
    codeNames[1] ="code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
    is_coupled_rank[2] = CWP_STATUS_ON;
    is_coupled_rank[3] = CWP_STATUS_ON;
  }
  else if (rank == 1) {
    n_code_name = 2;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] ="code1";
    codeNames[1] ="code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
    is_coupled_rank[2] = CWP_STATUS_ON;
    is_coupled_rank[3] = CWP_STATUS_ON;
  }
  else if (rank == 2) {
    n_code_name = 2;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] ="code1";
    codeNames[1] ="code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
    is_coupled_rank[2] = CWP_STATUS_ON;
    is_coupled_rank[3] = CWP_STATUS_ON;
  }
  else if (rank == 3) {
    n_code_name = 2;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] ="code1";
    codeNames[1] ="code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
    is_coupled_rank[2] = CWP_STATUS_ON;
    is_coupled_rank[3] = CWP_STATUS_ON;
  }

  

   
  char* fileName = NULL;
  if (rank == 0) 
    fileName="c_new_api_0000.txt";
  else if (rank == 1)
    fileName="c_new_api_0001.txt";
  else if (rank == 2)
    fileName="c_new_api_0002.txt";
  else if (rank == 3)
    fileName="c_new_api_0003.txt";

  outputFile = fopen(fileName,"w");

  times_init = malloc(sizeof(double) * n_code_name);

  //CWP_Output_file_set (outputFile);

  for (int i = 0; i < n_code_name; i++) {
    times_init[i] = 0; 
  }
  
  MPI_Comm *localComm = malloc(sizeof(MPI_Comm)*n_code_name);
  CWP_Init(MPI_COMM_WORLD,
           n_code_name,
           (const char **) codeNames,
           is_coupled_rank,
           times_init,
           localComm);


  /* Output redirection
   * ------------------ */

  int currentRank;
  int localCommSize;

  for (int i = 0; i < n_code_name; i++ ) {
    MPI_Comm_rank(localComm[i], &currentRank);
    MPI_Comm_size(localComm[i], &localCommSize);
  }

  /* Finalize
   * -------- *//*
  
  if (rank == 0 || rank == 1 || rank == 2 || rank == 5 || rank == 7) {
    int toto = 111;
    CWP_Param_lock ("code1");    
    //  MPI_Barrier (MPI_COMM_WORLD);
    CWP_Param_add ("code1", "toto", CWP_INT, &toto);
    char *A = "Bonjour !";
    CWP_Param_add ("code1", "toto2", CWP_CHAR, &A);
    CWP_Param_unlock ("code1");    
  }
  //  else {
 
  MPI_Barrier (MPI_COMM_WORLD);
    
  // }

  int titi;
  CWP_Param_get ("code1", "toto", CWP_INT, &titi);
  
  char *titi2;
  CWP_Param_get ("code1", "toto2", CWP_CHAR, &titi2);
  
  free (titi2);
  assert(titi == 111);

  CWP_Properties_dump ();
*/   
  char cpl_id1[] = "cpl_code1_code2";

  // cpl1
  
 if (rank == 0 || rank == 1 || rank == 2 || rank == 3 ) {
    CWP_Cpl_create ("code1", cpl_id1, "code2", CWP_COMM_PAR_WITH_PART,
                    CWP_GEOM_LOCATION, 1,
                    CWP_DISPLACEMENT_STATIC, CWP_FREQ_CPL_TIME_STEP);


  }
   

  printf("pass4\n");
  fflush(stdout);
  
  CWP_Status_t visu_status;
  
  if (rank == 0 || rank == 1 || rank == 2 || rank == 3  ) {
  CWP_Field_create (cpl_id1,"code1","unchamp",CWP_DOUBLE,CWP_FIELD_STORAGE_BLOCK,3,
                     CWP_FIELD_VALUE_CELL_POINT,
                     CWP_FIELD_EXCH_SENDRECV,
                     visu_status);
  }

  int n_component;
  n_component=CWP_Field_n_component_get("code1",cpl_id1,"unchamp");
  printf("n_component: %i\n",n_component);
  
  
  CWP_Field_value_t field_location = CWP_Field_location_get("code1",cpl_id1,"unchamp");
  printf("field_location: %d\n",field_location);
  
  
  CWP_Type_t field_Type = CWP_Field_type_get("code1",cpl_id1,"unchamp");
  printf("field_Type: %d\n",field_Type);
  
  
  CWP_Field_storage_t storage_type = CWP_Field_storage_get("code1",cpl_id1,"unchamp");
  printf("storage_type: %d\n",storage_type);


  CWP_Field_del("code1",cpl_id1,"unchamp");


  CWP_Finalize();

  MPI_Finalize();

  free (srcName);
  free (localComm);
  free (codeNames);
  free (is_coupled_rank);
  free (times_init);
  fclose (outputFile);

  return 0;
}
