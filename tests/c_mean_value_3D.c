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
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <mpi.h>

#include "cwipi.h"
#include "grid_mesh.h"


/*----------------------------------------------------------------------
 *                                                                     
 * Dump status exchange                                                
 *                                                                     
 * parameters:
 *   status              <-- Exchange status           
 *---------------------------------------------------------------------*/

static void _dumpStatus(FILE* outputFile, cwipi_exchange_status_t status)
{
  switch(status) {
  case CWIPI_EXCHANGE_OK :
    fprintf(outputFile, "Exchange Ok\n");
    break;
  case CWIPI_EXCHANGE_BAD_RECEIVING :
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
 * Main : surface coupling test : P1P1 
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
  char *srcBaseName = strrchr(srcName, '.');
  *srcBaseName = '\0';
  srcBaseName = strrchr(srcName, '/') + 1;

  if (rank == 0)
    printf("\nSTART: %s\n", srcBaseName);

  srand(rank+time(0));

  int n_partition = 0;
  while(2 * pow(n_partition, 2) < commWorldSize) n_partition++;

  int n2 = 2 * pow(n_partition, 2);

  if (n2 != commWorldSize) {
    if (rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  /* Read args from command line
   * --------------------------- */

  int nVertexSeg = 10;
  double randLevel = 0.1;

  _read_args(argc, argv, &nVertexSeg, &randLevel);

  /* Initialization
   * -------------- */

  char *codeName;
  int codeId;
  char *codeCoupledName;

  int third_size; 

  if (rank < commWorldSize / 2) {
    codeName = "code1";
    codeId = 1;
    codeCoupledName = "code2";
  }
  else {
    codeName = "code2";
    codeId = 2;
    codeCoupledName = "code1";
  }

  char* fileName = (char *) malloc(sizeof(char) * 25);
  sprintf(fileName,"c_mean_value_3D_%4.4d.txt",rank);

  outputFile = fopen(fileName,"w");

  free(fileName);

  cwipi_set_output_listing(outputFile);

  MPI_Comm localComm;
  cwipi_init(MPI_COMM_WORLD,
             codeName ,
             &localComm);

  /* Output redirection
   * ------------------ */

  int currentRank;
  int localCommSize;

  MPI_Comm_rank(localComm, &currentRank);
  MPI_Comm_size(localComm, &localCommSize);

  fprintf(outputFile, "  Mean value test\n");
  fprintf(outputFile, "\n");

  fprintf(outputFile, "\nDump after initialization\n");
  fprintf(outputFile, "-------------------------\n");
  cwipi_dump_application_properties();

  if (rank == 0)
    printf("        Create coupling\n");
  
  cwipi_solver_type_t solver_type;
  
  solver_type = CWIPI_SOLVER_CELL_VERTEX;
  
  /* Coupling creation
   * ----------------- */

  cwipi_create_coupling("c_mean_value_3D",                                // Coupling id
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                        codeCoupledName,                           // Coupled application id
                        3,                                         // Geometric entities dimension
                        0.1,                                       // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        solver_type,                               // Solver type
                        1,                                         // Postprocessing frequency
                        "EnSight Gold",                            // Postprocessing format
                        "text");                                   // Postprocessing option
  
  /* Mesh definition
   * --------------- */

  if (rank == 0)
    printf("        Create mesh\n");
  /* tetraedre */
  int nVertex = 7;               // Number of vertex
  double *coords = (double *) malloc(nVertex * 3 * sizeof(double));
  int nElts = 1;                 // Number of elements
  int *eltsConnecPointer = (int *) malloc(2 * sizeof(int));
  int *eltsConnec = (int *) malloc(4 * sizeof(int));

  coords[0] = 0.; 
  coords[1] = 0.; 
  coords[2] = 0.; 
  
  coords[3] = 0.; 
  coords[4] = 1.; 
  coords[5] = 0.; 
     
  coords[6] = 1./2; 
  coords[7] = 1./2; 
  coords[8] = 1.;

  coords[9] = 0.; 
  coords[10] = 0.; 
  coords[11] = 1.;
  
  coords[12] = 1; 
  coords[13] = 0; 
  coords[14] = 0; 
    
  coords[15] =  1; 
  coords[16] = 1; 
  coords[17] = 0;
  
  
  coords[18] = 1./2; 
  coords[19] = 1./2; 
  coords[20] = -1;

  /*coords[0] = 0; 
  coords[1] = 0; 
  coords[2] = 0; 
  
  coords[3] = 1; 
  coords[4] = 0; 
  coords[5] = 0; 
  
  coords[6] = 0; 
  coords[7] = 1; 
  coords[8] = 0; 
  
  coords[9] =  0; 
  coords[10] = 0; 
  coords[11] = 1;*/
  
  eltsConnecPointer[0] = 0; 
  eltsConnecPointer[1] = 4; 
  
  eltsConnec[0] = 1; 
  eltsConnec[1] = 4; 
  eltsConnec[2] = 3; 
  eltsConnec[3] = 2; 


  cwipi_define_mesh("c_mean_value_3D",
                    nVertex,
                    1,
                    coords,
                    eltsConnecPointer,
                    eltsConnec);




  nVertex = 6;               // Number of vertex
  /*octaedre*/

  int *face_index = (int*) malloc(nElts * sizeof(int));
  int *cell_to_face_connectivity = (int*) malloc((8+1) * sizeof(int));
  int *face_connectivity = (int*) malloc(3 * 8 * sizeof(int));
  int *face_connectivity_index = (int*) malloc((8+1) * sizeof(int));
 
 

  face_index[0] = 0;
  face_index[1] = 8;

  for(int i = 0; i < 9 ; i++)
    cell_to_face_connectivity[i] = i+1;

  for (int i = 0; i < 9 ; i++)
    face_connectivity_index[i] = 3*i;
     
  face_connectivity[0]  = 1;
  face_connectivity[1]  = 3;
  face_connectivity[2]  = 2;

  face_connectivity[3]  = 2;
  face_connectivity[4]  = 3;
  face_connectivity[5]  = 6;

  face_connectivity[6]  = 6;
  face_connectivity[7]  = 3;
  face_connectivity[8]  = 5;

  face_connectivity[9]  = 5;
  face_connectivity[10] = 3;
  face_connectivity[11] = 1;

  face_connectivity[12] = 1;
  face_connectivity[13] = 2;
  face_connectivity[14] = 7;

  face_connectivity[15] = 2;
  face_connectivity[16] = 6;
  face_connectivity[17] = 7;

  face_connectivity[18] = 6;
  face_connectivity[19] = 5;
  face_connectivity[20] = 7;

  face_connectivity[21] = 5;
  face_connectivity[22] = 1;
  face_connectivity[23] = 7;

  
  cwipi_add_polyhedra("c_mean_value_3D",
                     nElts,
                     face_index,
                     cell_to_face_connectivity,
                     face_connectivity_index,
                     face_connectivity);
  


  //ajouter des polyedres par cwipi_add_polyedra ...


  /* Points to locate
   * --------------- */

  int nPts = 7;               // Number of vertex
  double *coordsPts = (double *) malloc(nPts * 3 * sizeof(double));

  coordsPts[0] = 0.;
  coordsPts[1] = 0.;
  coordsPts[2] = 0.;

  coordsPts[3] = 1.;
  coordsPts[4] = 0.;
  coordsPts[5] = 0.;

  coordsPts[6] = 0.;
  coordsPts[7] = 1./2;
  coordsPts[8] = 1./2;
  
  coordsPts[9] = 0.;
  coordsPts[10] = 1./4;
  coordsPts[11] = 1./4;

  coordsPts[12] = 1./3;
  coordsPts[13] = 1./2;
  coordsPts[14] = -1./10;

  coordsPts[15] = 0.;
  coordsPts[16] = 0.;
  coordsPts[17] = 1.;

  coordsPts[18] = 1./8;
  coordsPts[19] = 1./3;
  coordsPts[20] = 1./2;

  //cwipi_set_points_to_locate ("c_mean_value_3D", nPts, coordsPts);

  double* value = (double *) malloc (nVertex * sizeof(double));
  double* lvalue = (double *) malloc (nVertex * sizeof(double));

  int nNotLocatedPoints;

  for (int i = 0; i < nVertex; i++) 
    value[i] = 1.2;

  cwipi_exchange_status_t status = cwipi_exchange("c_mean_value_3D",
                                                  "ech",
                                                  1,
                                                  1,     // n_step
                                                  0.1,   // physical_time
                                                  "val1",
                                                  value,
                                                  "val2",
                                                  lvalue,
                                                  &nNotLocatedPoints); 

  

  //cwipi_locate("c_mean_value_3D");
  

  cwipi_delete_coupling("c_mean_value_3D");
  
  /* Finalize
   * -------- */

  free(coords);
  free(coordsPts);
  free(eltsConnec);
  free(eltsConnecPointer);

  cwipi_finalize();

  MPI_Finalize();

  fclose(outputFile);

  return EXIT_SUCCESS;
}
