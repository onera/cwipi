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
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cwp.h"
#include "grid_mesh.h"


/*----------------------------------------------------------------------
 *                                                                     
 * Dump status exchange                                                
 *                                                                     
 * parameters:
 *   status              <-- Exchange status           
 *---------------------------------------------------------------------*/

static void _dumpStatus(FILE* outputFile, CWP_Status_t status)
{
  switch(status) {
  case CWP_STATUS_ON :
    fprintf(outputFile, "Exchange Ok\n");
    break;
  case CWP_STATUS_OFF :
    fprintf(outputFile, "Bad receiving\n");
    break;
  default :
    printf("Error : bad exchange status\n");
    exit(1);
  }
}

/*----------------------------------------------------------------------
 *                                                                     
 * Dump not located points                                             
 *                                                                     
 * parameters:
 *   coupling_id         <-- Coupling id               
 *   nNotLocatedPoints   <-- Number of not located points
 *---------------------------------------------------------------------*/

static void _dumpNotLocatedPoints(FILE* outputFile,
                                  const char *coupling_id,
                                  const int nNotLocatedPoints)
{
  if ( nNotLocatedPoints > 0) {
    fprintf(outputFile, "Not located points :\n");
    const int* notLocatedPoints = cwipi_get_not_located_points(coupling_id);
    for(int i = 0; i < nNotLocatedPoints; i++)
     fprintf(outputFile, "%i ", notLocatedPoints[i]);
    fprintf(outputFile, "\n");
  }
}

/*----------------------------------------------------------------------
 *                                                                     
 * Display usage                                             
 *                                                                     
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_usage(int exit_code)
{
  printf
    ("\n"
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
_read_args(int            argc,
           char         **argv,
           int          *nVertex,
           double       *randLevel)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *nVertex = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *randLevel = atof(argv[i]);
    }
    i++;
  }
}

/*----------------------------------------------------------------------
 *                                                                     
 * Main : surface coupling test : P1P0_P0P1 
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
  int commWorldSize;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

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

  srand(rank+time(0));

  int n_partition = 0;
  while(1.5 * (double)pow(n_partition, 2) < commWorldSize) n_partition++;

  int n2 = (int)(1.5 * pow(n_partition, 2));

  if (n2 != commWorldSize) {
    if (rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '1.5 * n^2' \n");
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  /* Read args from command line
   * --------------------------- */

  int nVertexSeg = 30;
  double randLevel = 0.4;

  _read_args(argc, argv, &nVertexSeg, &randLevel);

  /* Initialization
   * -------------- */

  int nb_part = 1; 

  int n_code_name = 0;
  char **codeName = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;

  char **codeCoupledName = NULL;
  int codeId;

  int rankCode1[ commWorldSize];
  int rankCode2[ commWorldSize];


  //Random distribution of codes on processes.
  int nsize = commWorldSize;   
  int nranks[nsize];
  
 
  
  if(rank==0){
    time_t t;
    /* Intializes random number generator */
    srand((unsigned) time(&t));
    int size_code1_domain = (double)commWorldSize * (2.0/3.0);
    int size_code2_domain = (double)commWorldSize * (2.0/3.0);
   
    if(sqrt(size_code1_domain)!=(int)sqrt(size_code1_domain)){
      size_code1_domain = pow((int)sqrt(size_code1_domain),2.0);
    }

    if(sqrt(size_code2_domain)!=(int)sqrt(size_code2_domain)){
      size_code2_domain = pow((int)sqrt(size_code2_domain),2.0);
    }
   

    for(int i=0; i < commWorldSize; i++) {
      rankCode1[i] = 0;
      rankCode2[i] = 0;
    }
    for(int i=0; i < size_code1_domain; i++){
      int tmp = i;//rand() % commWorldSize;   
      while(rankCode1[tmp] == 1){
        tmp = rand() % commWorldSize;
      }
      rankCode1[tmp]=1; 
    }

    for(int i=0; i < size_code2_domain; i++){
      int tmp = size_code1_domain/2 + i;//rand() % commWorldSize;   
      while(rankCode2[tmp] == 1){
        tmp = rand() % commWorldSize;
      }
      rankCode2[tmp]=1; 
    }
 
 /*
    rankCode1[0]=1; 
    rankCode1[1]=1; 
    rankCode1[2]=1;
    rankCode1[3]=1;  
    rankCode1[4]=0; 
    rankCode1[5]=0; 
  
    rankCode2[0]=0; 
    rankCode2[1]=0; 
    rankCode2[2]=1;
    rankCode2[3]=1;  
    rankCode2[4]=1; 
    rankCode2[5]=1; 
  
    */

    int ind=0;
    for(int i=0; i < commWorldSize; i++) {
      if(rankCode2[i] || rankCode1[i]){
       nranks[ind]=i;
       ind++;
      }
    }    
  } // if rank==0


  int nsizercv=nsize;   
  
  int rankCode1Rcv = -1;
  int rankCode2Rcv = -1;

  MPI_Scatter(&rankCode1    ,1 ,MPI_INT,
              &rankCode1Rcv,1 ,MPI_INT,
              0,
              MPI_COMM_WORLD
             );
    
  MPI_Scatter(&rankCode2    ,1 ,MPI_INT,
              &rankCode2Rcv,1 ,MPI_INT,
              0,
              MPI_COMM_WORLD
             );

  MPI_Scatter(&rankCode2    ,1 ,MPI_INT,
              &rankCode2Rcv,1 ,MPI_INT,
              0,
              MPI_COMM_WORLD
             );
   


  if(rankCode1Rcv && rankCode2Rcv) {
    printf("ZZI am process %i and I own the code1 and the code2.\n",rank);
  }
  else {
    if(rankCode1Rcv){
      printf("ZZI am process %i and I own the code1.\n",rank);
    }
    if(rankCode2Rcv){
      printf("ZZI am process %i and I own the code2.\n",rank);
    }    
  } 
  
  
  if (rankCode1Rcv) {
    n_code_name = 1;
    codeName = malloc(sizeof(char*)*n_code_name);
    codeCoupledName = malloc(sizeof(char*)*n_code_name);
    codeName[0] = "code1";
    codeCoupledName[0] = "code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;    
  }
  
  if (rankCode2Rcv) {
    n_code_name = 1;
    codeName = malloc(sizeof(char*)*n_code_name);
    codeCoupledName = malloc(sizeof(char*)*n_code_name);
    codeName[0] = "code2";
    codeCoupledName[0] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;    
  }

  if (rankCode1Rcv && rankCode2Rcv) {
    n_code_name = 2;
    codeName        = malloc(sizeof(char*)*n_code_name);
    codeCoupledName = malloc(sizeof(char*)*n_code_name);
    codeName[0] = "code1";
    codeName[1] = "code2";
    codeCoupledName[0] = "code2";
    codeCoupledName[1] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;    
    is_coupled_rank[1] = CWP_STATUS_ON;  
  }  
 
  char* fileName = (char *) malloc(sizeof(char) * 35);
  sprintf(fileName,"c_surf_coupling_P1P0_P0P1_%4.4d.txt",rank);
  
  outputFile = fopen(fileName,"w");

  free(fileName);
  CWP_Output_file_set (outputFile);


  times_init = malloc(sizeof(double) * n_code_name);

  for (int i = 0; i < n_code_name; i++) {
    times_init[i] = 0; 
  }
  
  
  MPI_Comm *localComm = malloc(sizeof(MPI_Comm)*n_code_name);
  CWP_Init(MPI_COMM_WORLD,
           n_code_name,
           (const char **) codeName,
           is_coupled_rank,
           times_init,
           localComm);
  
  /* Output redirection
   * ------------------ */

  int currentRank[n_code_name];
  int localCommSize[n_code_name];

  for(int i_code = 0; i_code < n_code_name; i_code++) { 
    if(localComm[i_code]!=MPI_COMM_NULL) MPI_Comm_rank(localComm[i_code], &currentRank[i_code]);
    if(localComm[i_code]!=MPI_COMM_NULL) MPI_Comm_size(localComm[i_code], &localCommSize[i_code]);
  }
  
  fprintf(outputFile, "  Surface coupling test : P1P0_P0P1 with polygon\n");
  fprintf(outputFile, "\n");

  fprintf(outputFile, "\nDump after initialization\n");
  fprintf(outputFile, "-------------------------\n");


  /* Coupling creation
   * ----------------- */
  if (rank == 0)
    printf("        Create coupling \n");


  for(int i_code = 0; i_code < n_code_name; i_code++) { 
    CWP_Cpl_create (codeName[i_code],             // Code name
                    "c_new_api_surf_cpl_P1P0_P0P1_dynamic",  // Coupling id
                    codeCoupledName[i_code],         // Coupled application id
                    CWP_COMM_PAR_WITH_PART,  // Coupling type
                    CWP_SPATIAL_INTERP_FROM_LOCATION,       // Solver type
                    nb_part,                 // Partition number
                    CWP_DISPLACEMENT_DEFORMABLE, // Mesh displacement type
                    CWP_FREQ_CPL_TIME_STEP); // Postprocessing frequency
  }
  
  if (rank == 0)
    printf("        Set visu\n");


  for(int i_code = 0; i_code < n_code_name; i_code++) {                     
    CWP_Visu_set(codeName[i_code], // Code name
                 "c_new_api_surf_cpl_P1P0_P0P1_dynamic",     // Coupling id
                 1,           // Postprocessing frequency
                 Ensight,     // Postprocessing format
                 "text");     // Postprocessing option
  }

  /* Mesh definition
   * --------------- */

  if (rank == 0)
    printf("        Create mesh\n");
  
  int nVertex = 0;               // Number of vertex
  double **coords = NULL;         // Vertex coordinates
  int nElts = 0;                 // Number of elements
  int **eltsConnecPointer = NULL; // Connectivity index
  int **eltsConnec = NULL;        // Connectivity
  
  /* Domain bounds */

  const double xmin = -10;
  const double xmax =  10;
  const double ymin = -10;
  const double ymax =  10;

  nVertex = nVertexSeg * nVertexSeg;
  nElts = (nVertexSeg - 1) * (nVertexSeg - 1);
  
  coords = (double **) malloc(sizeof(double*) * n_code_name );
  eltsConnecPointer = (int **) malloc(sizeof(int*) * n_code_name);
  eltsConnec = (int **) malloc(sizeof(int*) * n_code_name);

  srand (time(NULL));

  for(int i_code = 0; i_code < n_code_name; i_code++) { 
    coords           [i_code] = (double *) malloc(sizeof(double) * 3 * nVertex );
    eltsConnecPointer[i_code] = (int *)    malloc(sizeof(int) * (nElts + 1));
    eltsConnec       [i_code] = (int *)    malloc(sizeof(int) * 4 * nElts);
    grid_mesh(xmin, 
            xmax, 
            ymin, 
            ymax, 
            randLevel,
            nVertexSeg,
            sqrt(localCommSize[i_code]), 
            coords           [i_code], 
            eltsConnecPointer[i_code],
            eltsConnec       [i_code],
            localComm        [i_code]);
    
    carre2rond(xmin, 
	       xmax, 
	       ymin, 
	       ymax,
	       coords[i_code], 
	       nVertex);

  }
 
  //MPI_Barrier(MPI_COMM_WORLD); 


  fprintf(outputFile, "   Number of vertex   : %i\n", nVertex);
  fprintf(outputFile, "   Number of elements : %i\n", nElts);

  for(int i_code = 0; i_code < n_code_name; i_code++) {        
    CWP_Mesh_interf_vtx_set(codeName[i_code],             //Code name
                            "c_new_api_surf_cpl_P1P0_P0P1_dynamic",  // Coupling id
                            0,
                            nVertex,
                            coords[i_code],
                            NULL);  

    int block_id = CWP_Mesh_interf_block_add(codeName[i_code],             // Code name
                                           "c_new_api_surf_cpl_P1P0_P0P1_dynamic",  // Coupling id
                                           CWP_BLOCK_FACE_POLY);
                       
    CWP_Mesh_interf_f_poly_block_set(codeName[i_code],             //Code name
                                     "c_new_api_surf_cpl_P1P0_P0P1_dynamic",  // Coupling id
                                     0,
                                     block_id,
                                     nElts,
                                     eltsConnecPointer[i_code],
                                     eltsConnec       [i_code],
                                     NULL);

  }

  
  for(int i_code = 0; i_code < n_code_name; i_code++) {                                      
    CWP_Mesh_interf_finalize (codeName[i_code],
                              "c_new_api_surf_cpl_P1P0_P0P1_dynamic"  // Coupling id
                             );
  }
  
  /* Fields exchange
   *     - Proc 0 : Send X coordinates
   *                Recv Y coordinates
   *     - Proc 1 : Send Y coordinates
   *                Recv X coordinates
   * --------------------------------- */

 // if (rank == 0)
  printf("        Exchange Code1 <-> Code2 %i\n",rank);

  double **sendValues = (double **) malloc(sizeof(double*) * n_code_name);
  double **recvValues = (double **) malloc(sizeof(double*) * n_code_name);
  for(int i_code = 0; i_code < n_code_name; i_code++) {   
    if (codeName[i_code] == "code1") {
      sendValues[i_code] = (double *) malloc(sizeof(double) * 3* nVertex);
      recvValues[i_code] = (double *) malloc(sizeof(double) * nElts); 
      for (int i = 0; i <nVertex; i++) {
        sendValues[i_code][3*i] = coords[i_code][3 * i];
        sendValues[i_code][3*i+1] = coords[i_code][3 * i+1];
        sendValues[i_code][3*i+2] = coords[i_code][3 * i+2];
      }
    }
    else {
      sendValues[i_code] = (double *) malloc(sizeof(double) * nElts);
      recvValues[i_code] = (double *) malloc(sizeof(double) * 3*nVertex); 
      for (int i = 0; i <nElts; i++) {
        sendValues[i_code][i] = rank;//i;//coords[3 * i];
      }  
    }
  
  }
  /* Define fields to send (X coordinate or Y coordinate) */

  /* Exchange */

  int nNotLocatedPoints = 0;
  char *fieldName1;
  char *fieldName2;
  
  fieldName1 = "cooX_t0";
  fieldName2 = "code2_elt_rank";

  CWP_Status_t visu_status = CWP_STATUS_OFF; 
  
  printf("        Field %i\n",rank);
  for(int i_code = 0; i_code < n_code_name; i_code++) {   
    if (codeName[i_code] == "code1") {

     CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_dynamic",
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       3,
                       CWP_FIELD_VALUE_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

     CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_dynamic",
                       fieldName2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_FIELD_VALUE_CELL_POINT,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName1,0, sendValues[i_code]); 
      CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName2,0, recvValues[i_code]);
//    
    
 //   _dumpStatus(outputFile, status);
//    _dumpNotLocatedPoints(outputFile, "c_new_api_surf_cpl_P1P0_P0P1_dynamic", nNotLocatedPoints);

    }
    else {
  
     CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_dynamic",
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       3,
                       CWP_FIELD_VALUE_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
                     
      CWP_Field_create (codeName[i_code],
                     "c_new_api_surf_cpl_P1P0_P0P1_dynamic",
                     fieldName2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_BLOCK,
                     1,
                     CWP_FIELD_VALUE_CELL_POINT,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName2,0, sendValues[i_code]); 
    CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName1,0, recvValues[i_code]);

  //  _dumpStatus(outputFile, status);
  //  _dumpNotLocatedPoints(outputFile, "c_new_api_surf_cpl_P1P0_P0P1_dynamic", nNotLocatedPoints);

    }

  }

   
  MPI_Barrier(MPI_COMM_WORLD);
  double recv_time = 0.;
  for (int il_nb_ite = 0; il_nb_ite < 10; il_nb_ite++) {
    recv_time += 1.0;
 
    /* Mesh rotation and new localisation - Mapping Compute */


    for(int i_code = 0; i_code < n_code_name; i_code++) {  
      if(codeName[i_code]=="code2") mesh_rotate(coords[i_code], nVertex, 3*recv_time); 
      if(codeName[i_code]=="code1") mesh_rotate(coords[i_code], nVertex, recv_time); 
      CWP_next_recv_time_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",recv_time);
    }
    
    CWP_Spatial_interp_weights_compute("c_new_api_surf_cpl_P1P0_P0P1_dynamic");
          
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i_code = 0; i_code < n_code_name; i_code++) {    
      
      if (codeName[i_code] == "code1") {
	CWP_Issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName1);   
	CWP_Irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName2);  
      }
      else {
	CWP_Irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName1);   
	CWP_Issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName2);    
      }
    }

    for(int i_code = 0; i_code < n_code_name; i_code++) {       
      if (codeName[i_code] == "code1") {
	CWP_Wait_issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName1);  
	CWP_Wait_irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName2);    
      }
      else {
	CWP_Wait_irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName1);       
	CWP_Wait_issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic",fieldName2);   
      }    
    }
  }
  /* Coupling deletion
   * ----------------- */

  //if (rank == 0)
    printf("        Delete mesh\n",rank);
 
  for(int i_code = 0; i_code < n_code_name; i_code++) 
     CWP_Mesh_interf_del(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic");

  if (rank == 0)
    printf("        Delete coupling\n");

  for(int i_code = 0; i_code < n_code_name; i_code++)    
    CWP_Cpl_del(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_dynamic");

  /* Freeing memory
   * -------------- */

  free(coords);

  free(sendValues);   
  free(recvValues);
  free(srcName);

  /* Finalize
   * -------- */

  CWP_Finalize();

  MPI_Finalize();

  fclose(outputFile);

  return EXIT_SUCCESS;
}



