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
#include <math.h>
#include "grid_mesh.h"
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

  
  int n_partition = 0;
  const int two = 2;
  while(two * pow(n_partition, two) < comm_world_size) n_partition++;

  int n2 = (int) (two * pow(n_partition, two));

  if (n2 != comm_world_size) {
    if (rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  
  
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

  if (rank == 0 ) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] ="code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
  }
   
  if (rank == 1 ) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] ="code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
  }
  
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
    printf("Size of localComm[%i]=%i et rang du proc=%i.\n",i,localCommSize,currentRank );
  }

  /* Finalize
   * -------- */
 
  char cpl_id1[] = "cpl_code1_code2";

  printf("Coupling creation %s\n",codeNames[0]);
  if(rank==0) {

    CWP_Cpl_create ("code1", cpl_id1, "code2", CWP_COMM_PAR_WITHOUT_PART,
                    CWP_GEOM_LOCATION, 1,
                    CWP_DISPLACEMENT_STATIC, CWP_FREQ_CPL_TIME_STEP);
                    
     CWP_Visu_set("code1", cpl_id1,1.25,Ensight,"binary");
  
     printf("Visu Set\n");
                      
    }

   if(rank==1) {
    CWP_Cpl_create ("code2", cpl_id1, "code1", CWP_COMM_PAR_WITHOUT_PART,
                    CWP_GEOM_LOCATION, 1,
                    CWP_DISPLACEMENT_STATIC, CWP_FREQ_CPL_TIME_STEP);
                    
    CWP_Visu_set("code2", cpl_id1,1.25,Ensight,"binary");
                    
    }   
       

  printf("Coupling created %s\n",codeNames[0]);
  
  
  


  
  
  
 
  int nVertex=3;                

  //  coords = (double *) malloc(sizeof(double) * 3 * nVertex );

  int nVertexSeg    = 10 ;
  double randLevel  = 0.2;
  int postFreq      = -1 ;
  int t_com         = 0  ;
  int tostdout      = 0  ;
  int randFromClock = 0  ;
  int *eltsConnecPointer = NULL; // Connectivity index
  int *eltsConnec = NULL;        // Connectivity
  double *coords  = NULL; 
  int     nElts;

  const double xmin = -10;
  const double xmax =  10;
  const double ymin = -10;
  const double ymax =  10;
  
  nVertex = nVertexSeg * nVertexSeg;
  nElts = (nVertexSeg - 1) * (nVertexSeg - 1);

  coords = (double *) malloc(sizeof(double) * 3 * nVertex );
  eltsConnecPointer = (int *) malloc(sizeof(int) * (nElts + 1));
  eltsConnec = (int *) malloc(sizeof(int) * 4 * nElts);
  
  printf("Mesh grid creation\n");
  grid_mesh(xmin, 
            xmax, 
            ymin, 
            ymax, 
            randLevel,
            nVertexSeg,
            n_partition, 
            coords, 
            eltsConnecPointer,
            eltsConnec,
            localComm[0]); 

  char* code_name = (char*)malloc(sizeof(char)*5);
  if(rank==0) {
    code_name = "code1";
  }
  if(rank==1) {
    code_name = "code2";
  }
  
  printf("vtx_set\n");
  CWP_Mesh_interf_vtx_set(code_name, cpl_id1,0,nVertex,coords,NULL);
 
  printf("Standard Block Add\n");
    
  int block_id = CWP_Mesh_interf_block_add(code_name, cpl_id1,CWP_BLOCK_FACE_QUAD4);
    
  printf("Standard block set\n");
  CWP_Mesh_interf_block_std_set(code_name, cpl_id1,0,block_id,nElts,eltsConnec,NULL);

  printf("Create a Field\n");
  CWP_Status_t visu_status = CWP_STATUS_OFF;
  
  
  if(rank==0){

    CWP_Field_create (code_name,cpl_id1,"unchamp",CWP_DOUBLE,CWP_FIELD_STORAGE_BLOCK,1,
                     CWP_FIELD_VALUE_CELL_POINT,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    double data[nElts];
    for(int i=0;i<nElts;i++)
       {data[i]=(double)i;
        printf("data[i] %f\n",data[i]);
       }
    CWP_Field_data_set(code_name,cpl_id1,"unchamp",0,data);
  
  }
  
  if(rank==1){

    CWP_Field_create (code_name,cpl_id1,"unchamp",CWP_DOUBLE,CWP_FIELD_STORAGE_BLOCK,1,
                     CWP_FIELD_VALUE_CELL_POINT,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);       
  }

  int* n_uncomputed_tgt;
  CWP_Geom_compute(code_name,cpl_id1,CWP_FIELD_VALUE_CELL_POINT, n_uncomputed_tgt);
  
  if(rank==0){
    CWP_Issend (code_name,cpl_id1,"unchamp");
    CWP_Wait_issend (code_name,cpl_id1,"unchamp");
  }

  if(rank==1){
    CWP_Irecv (code_name,cpl_id1,"unchamp");   
    CWP_Wait_irecv (code_name,cpl_id1,"unchamp");
  }

  
  
  printf("Interface Mesh deletion\n");
  CWP_Mesh_interf_del(code_name, cpl_id1);
 
 

  fflush(stdout);

  CWP_Finalize();

  MPI_Finalize();

  free (srcName);
  free (localComm);
  free (codeNames);
  free (is_coupled_rank);
  free (times_init);

  return 0;
}
