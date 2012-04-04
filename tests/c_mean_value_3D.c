/*
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
                        1e-6,                                       // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        solver_type,                               // Solver type
                        1,                                         // Postprocessing frequency
                        "EnSight Gold",                            // Postprocessing format
                        "text");                                   // Postprocessing option
  
  /* Mesh definition
   * --------------- */

  if (rank == 0)
    printf("        Create mesh\n");


  /**** Coord tetraedre *****/


  /*int nVertex = 4;               // Number of vertex
  double *coords = (double *) malloc(nVertex * 3 * sizeof(double));

  coords[0] = 0.; 
  coords[1] = 0.; 
  coords[2] = 0.; 
  
  coords[3] = 1.; 
  coords[4] = 0.; 
  coords[5] = 0.; 
     
  coords[6] = 0.; 
  coords[7] = 1.; 
  coords[8] = 0.;

  coords[9] = 0.; 
  coords[10] = 0.; 
  coords[11] = 1.;

  /*  coord tetraedre + octaedre */

  /* int nVertex = 7;               // Number of vertex
  double *coords = (double *) malloc(nVertex * 3 * sizeof(double));

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
  
  /*  coord tetraedre + dodecaedre */

  /* double X = 1.;
  double Z = (1./2)*(1 + sqrt(5));

  int nVertex = 13;               // Number of vertex
  double *coords = (double *) malloc(nVertex * 3 * sizeof(double));

  coords[0] = 0;
  coords[1] =  Z/3 + 1*X ;
  coords[2] =  (X + 2*Z)/3 + 1.5*X;

  coords[3] = -X;
  coords[4] =  0;
  coords[5] =  Z;

  coords[6] =  X;
  coords[7] =  0;
  coords[8] =  Z;

  coords[9] = -X;
  coords[10] =  0;
  coords[11] = -Z;

  coords[12] =  X;
  coords[13] =  0;
  coords[14] = -Z;

  coords[15] =  0;
  coords[16] =  Z;
  coords[17] =  X;

  coords[18] =  0;
  coords[19] =  Z;
  coords[20] = -X;

  coords[21] =  0;
  coords[22] = -Z;
  coords[23] =  X;

  coords[24] =  0;
  coords[25] = -Z;
  coords[26] = -X;

  coords[27] =  Z;
  coords[28] =  X;
  coords[29] =  0;

  coords[30] = -Z;
  coords[31] =  X;
  coords[32] =  0;

  coords[33] =  Z;
  coords[34] = -X;
  coords[35] =  0;

  coords[36] = -Z;
  coords[37] = -X;
  coords[38] =  0;

  /*********  coord tetraedre + octaedre non convexe ********/

  /* int nVertex = 8;               // Number of vertex
  double *coords = (double *) malloc(nVertex * 3 * sizeof(double));

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
  coords[20] = 1./2;


  /*********  coord tetraedre + 2 octaedres ********/

  /*int nVertex = 9;               // Number of vertex
  double *coords = (double *) malloc(nVertex * 3 * sizeof(double));

  coords[0] = 0.; 
  coords[1] = 0.; 
  coords[2] = 0.; 
  
  coords[3] = 0.; 
  coords[4] = 1.; 
  coords[5] = 0.; 
     
  coords[6] = 1.; 
  coords[7] = 1.; 
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
  
  
  coords[18] = 1.; 
  coords[19] = 1.; 
  coords[20] = -1.;

  coords[21] = 2.; 
  coords[22] = 0.; 
  coords[23] = 0.;

  coords[24] = 2.; 
  coords[25] = 1.; 
  coords[26] = 0.;

  /******* tetra + pyramide *******/

  int nVertex = 14;               // Number of vertex
  double *coords = (double *) malloc(nVertex * 3 * sizeof(double));

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
  
  coords[12] = 1.; 
  coords[13] = 0.; 
  coords[14] = 0.; 
    
  coords[15] = 1.; 
  coords[16] = 1.; 
  coords[17] = 0.;

  coords[18] = 0.; 
  coords[19] = 0.; 
  coords[20] = -1.;

  coords[21] = 0.; 
  coords[22] = 1.; 
  coords[23] = -1.;

  coords[24] = 1.5; 
  coords[25] = 0.; 
  coords[26] = -0.5;
  
  coords[27] = 0.5; 
  coords[28] = 0.; 
  coords[29] = -1.5;

  coords[30] = 0.4; 
  coords[31] = 0.8; 
  coords[32] = -1.2;

  coords[33] = 1.8; 
  coords[34] = 1.3; 
  coords[35] = 0.2;

  coords[36] = 1.5; 
  coords[37] = 1.; 
  coords[38] = 1.7;

  coords[39] = 1.1; 
  coords[40] = 0.3; 
  coords[41] = 0.1;


  /******tetraedre + pyramide*******/

  int nElts = 4;                 // Number of elements
  int *eltsConnecPointer = (int *) malloc((nElts+1) * sizeof(int));
  int *eltsConnec = (int *) malloc( 23* sizeof(int));


  eltsConnecPointer[0] = 0; 
  eltsConnecPointer[1] = 4; 
  eltsConnecPointer[2] = 9;
  eltsConnecPointer[3] = 15;
  eltsConnecPointer[4] = 23;


  eltsConnec[0] = 1; 
  eltsConnec[1] = 3; 
  eltsConnec[2] = 2; 
  eltsConnec[3] = 4; 

  eltsConnec[4] = 1; 
  eltsConnec[5] = 5; 
  eltsConnec[6] = 6; 
  eltsConnec[7] = 2; 
  eltsConnec[8] = 3; 

  eltsConnec[9] =  1; 
  eltsConnec[10] = 5; 
  eltsConnec[11] = 7; 
  eltsConnec[12] = 2; 
  eltsConnec[13] = 6; 
  eltsConnec[14] = 8;

  eltsConnec[15] = 5; 
  eltsConnec[16] = 7; 
  eltsConnec[17] = 8; 
  eltsConnec[18] = 6; 
  eltsConnec[19] = 9; 
  eltsConnec[20] = 10;
  eltsConnec[21] = 11;
  eltsConnec[22] = 12;


  cwipi_define_mesh("c_mean_value_3D",
                    nVertex,
                    nElts,
                    coords,
                    eltsConnecPointer,
                    eltsConnec);


  /************ forme quelconque ( au moins 1 face avec > 4 sommets) **************/
  
 /*  int nbr_face = 4; */
  
 /*  int *face_index = (int*) malloc(2 * sizeof(int)); */
 /*  int *cell_to_face_connectivity = (int*) malloc(nbr_face* sizeof(int)); */
 /*  int *face_connectivity = (int*) malloc(16 * sizeof(int)); */
 /*  int *face_connectivity_index = (int*) malloc((nbr_face+1) * sizeof(int)); */

 /*  face_index[0] = 0; */
 /*  face_index[1] = 4; */

 /*  for(int i = 0; i < 4 ; i++) */
 /*    cell_to_face_connectivity[i] = i+1; */

 /*  face_connectivity_index[0] = 0; */
 /*  face_connectivity_index[1] = 5;   */
 /*  face_connectivity_index[2] = 9; */
 /*  face_connectivity_index[3] = 13; */
 /*  face_connectivity_index[4] = 16; */
    
 /*  face_connectivity[0]  = 5; */
 /*  face_connectivity[1]  = 3; */
 /*  face_connectivity[2]  = 13; */
 /*  face_connectivity[3]  = 12; */
 /*  face_connectivity[4]  = 9; */

 /*  face_connectivity[5]  = 9; */
 /*  face_connectivity[6]  = 12; */
 /*  face_connectivity[7]  = 6; */
 /*  face_connectivity[8]  = 5; */

 /*  face_connectivity[9]  = 6; */
 /*  face_connectivity[10]  = 12; */
 /*  face_connectivity[11]  = 13; */
 /*  face_connectivity[12]  = 3; */

 /*  face_connectivity[13]  = 5; */
 /*  face_connectivity[14]  = 6; */
 /*  face_connectivity[15]  = 3; */

 /* cwipi_add_polyhedra("c_mean_value_3D", */
 /*                     1, */
 /*                     face_index, */
 /*                     cell_to_face_connectivity, */
 /*                     4, */
 /*                     face_connectivity_index, */
 /*                     face_connectivity); */

  /************octaedre**************/

  
  /*int nbr_face = 8;
  
  int *face_index = (int*) malloc(2 * sizeof(int));
  int *cell_to_face_connectivity = (int*) malloc(nbr_face* sizeof(int));
  int *face_connectivity = (int*) malloc(3 * nbr_face * sizeof(int));
  int *face_connectivity_index = (int*) malloc((nbr_face+1) * sizeof(int));
 
 

  face_index[0] = 0;
  face_index[1] = 8;

  for(int i = 0; i < 8 ; i++)
    cell_to_face_connectivity[i] = i+1;

  for (int i = 0; i < 9 ; i++)
    face_connectivity_index[i] = 3*i;
     
  face_connectivity[0]  = 1;
  face_connectivity[1]  = 2;
  face_connectivity[2]  = 3;

  face_connectivity[3]  = 2;
  face_connectivity[4]  = 6;
  face_connectivity[5]  = 3;

  face_connectivity[6]  = 6;
  face_connectivity[7]  = 5;
  face_connectivity[8]  = 3;

  face_connectivity[9]  = 5;
  face_connectivity[10] = 1;
  face_connectivity[11] = 3;

  face_connectivity[12] = 1;
  face_connectivity[13] = 7;
  face_connectivity[14] = 2;

  face_connectivity[15] = 2;
  face_connectivity[16] = 7;
  face_connectivity[17] = 6;

  face_connectivity[18] = 6;
  face_connectivity[19] = 7;
  face_connectivity[20] = 5;

  face_connectivity[21] = 5;
  face_connectivity[22] = 7;
  face_connectivity[23] = 1;

  
  cwipi_add_polyhedra("c_mean_value_3D",
                     nElts,
                     face_index,
                     cell_to_face_connectivity,
                     face_connectivity_index,
                     face_connectivity);

  /************deux octaedres**************/


  /*int nbr_face = 14;
  
  int *face_index = (int*) malloc(3 * sizeof(int));
  int *cell_to_face_connectivity = (int*) malloc(16* sizeof(int));
  int *face_connectivity = (int*) malloc(3 * nbr_face * sizeof(int));
  int *face_connectivity_index = (int*) malloc((nbr_face+1) * sizeof(int));
 
 

  face_index[0] = 0;
  face_index[1] = 8;
  face_index[2] = 16;


  for(int i = 0; i < 8 ; i++)
    cell_to_face_connectivity[i] = i+1;

  cell_to_face_connectivity[8] = -3;
  cell_to_face_connectivity[9] = 9;
  cell_to_face_connectivity[10] = 10;
  cell_to_face_connectivity[11] = 11;
  cell_to_face_connectivity[12] = -7;
  cell_to_face_connectivity[13] = 12;
  cell_to_face_connectivity[14] = 13;
  cell_to_face_connectivity[15] = 14;



  for (int i = 0; i < nbr_face + 1 ; i++)
    face_connectivity_index[i] = 3*i;
     
  face_connectivity[0]  = 1;
  face_connectivity[1]  = 2;
  face_connectivity[2]  = 3;

  face_connectivity[3]  = 2;
  face_connectivity[4]  = 6;
  face_connectivity[5]  = 3;

  face_connectivity[6]  = 6;
  face_connectivity[7]  = 5;
  face_connectivity[8]  = 3;

  face_connectivity[9]  = 5;
  face_connectivity[10] = 1;
  face_connectivity[11] = 3;

  face_connectivity[12] = 1;
  face_connectivity[13] = 7;
  face_connectivity[14] = 2;

  face_connectivity[15] = 2;
  face_connectivity[16] = 7;
  face_connectivity[17] = 6;

  face_connectivity[18] = 6;
  face_connectivity[19] = 7;
  face_connectivity[20] = 5;

  face_connectivity[21] = 5;
  face_connectivity[22] = 7;
  face_connectivity[23] = 1;

  face_connectivity[24] = 6;
  face_connectivity[25] = 9;
  face_connectivity[26] = 3;

  face_connectivity[27] = 9;
  face_connectivity[28] = 8;
  face_connectivity[29] = 3;

  face_connectivity[30] = 8;
  face_connectivity[31] = 5;
  face_connectivity[32] = 3;

  face_connectivity[33] = 6;
  face_connectivity[34] = 7;
  face_connectivity[35] = 9;

  face_connectivity[36] = 9;
  face_connectivity[37] = 7;
  face_connectivity[38] = 8;

  face_connectivity[39] = 8;
  face_connectivity[40] = 7;
  face_connectivity[41] = 5;

  
  cwipi_add_polyhedra("c_mean_value_3D",
                     2,
                     face_index,
                     cell_to_face_connectivity,
                     face_connectivity_index,
                     face_connectivity);

  /**********Dodecaedre**********/

  /* int nbr_face = 20;
  
  int *face_index = (int*) malloc(2 * sizeof(int));
  int *cell_to_face_connectivity = (int*) malloc(nbr_face*sizeof(int));
  int *face_connectivity = (int*) malloc(3 * nbr_face * sizeof(int));
  int *face_connectivity_index = (int*) malloc((nbr_face+1) * sizeof(int));
  
  
  
  face_index[0] = 0;
  face_index[1] = nbr_face;

  for(int i = 0; i < nbr_face ; i++)
    cell_to_face_connectivity[i] = i+1;

  for (int i = 0; i < nbr_face+1 ; i++)
    face_connectivity_index[i] = 3*i;
     
  face_connectivity[0]  = 2;
  face_connectivity[1]  = 6;
  face_connectivity[2]  = 3;

  face_connectivity[3]  = 2;
  face_connectivity[4]  = 11;
  face_connectivity[5]  = 6;

  face_connectivity[6]  = 11;
  face_connectivity[7]  = 7;
  face_connectivity[8]  = 6;

  face_connectivity[9]  = 6;
  face_connectivity[10]  = 7;
  face_connectivity[11]  = 10;

  face_connectivity[12]  = 6;
  face_connectivity[13]  = 10;
  face_connectivity[14]  = 3;

  face_connectivity[15]  = 10;
  face_connectivity[16]  = 12;
  face_connectivity[17]  = 3;

  face_connectivity[18]  = 10;
  face_connectivity[19]  = 5;
  face_connectivity[20]  = 12;

  face_connectivity[21]  = 7;
  face_connectivity[22]  = 5;
  face_connectivity[23]  = 10;

  face_connectivity[24]  = 7;
  face_connectivity[25]  = 4;
  face_connectivity[26]  = 5;

  face_connectivity[27]  = 4;
  face_connectivity[28]  = 9;
  face_connectivity[29]  = 5;

  face_connectivity[30]  = 9;
  face_connectivity[31]  = 12;
  face_connectivity[32]  = 5;

  face_connectivity[33]  = 9;
  face_connectivity[34]  = 8;
  face_connectivity[35]  = 12;

  face_connectivity[36]  = 9;
  face_connectivity[37]  = 13;
  face_connectivity[38]  = 8;

  face_connectivity[39]  = 13;
  face_connectivity[40]  = 2;
  face_connectivity[41]  = 8;

  face_connectivity[42]  = 2;
  face_connectivity[43]  = 3;
  face_connectivity[44]  = 8;

  face_connectivity[45]  = 8;
  face_connectivity[46]  = 3;
  face_connectivity[47]  = 12;

  face_connectivity[48]  = 11;
  face_connectivity[49]  = 2;
  face_connectivity[50]  = 13;

  face_connectivity[51]  = 11;
  face_connectivity[52]  = 13;
  face_connectivity[53]  = 4;
                                                 
  face_connectivity[54]  = 11;
  face_connectivity[55]  = 4;
  face_connectivity[56]  = 7;

  face_connectivity[57]  = 9;
  face_connectivity[58]  = 4;
  face_connectivity[59]  = 13;

  cwipi_add_polyhedra("c_mean_value_3D",
                     nElts,
                     face_index,
                     cell_to_face_connectivity,
                     face_connectivity_index,
                     face_connectivity);


  /*********** Pyramide **********/

  /*int nbr_face = 5;
  
  int *face_index = (int*) malloc(2 * sizeof(int));
  int *cell_to_face_connectivity = (int*) malloc(nbr_face* sizeof(int));
  int *face_connectivity = (int*) malloc(16 * sizeof(int));
  int *face_connectivity_index = (int*) malloc((nbr_face+1) * sizeof(int));
 
 

  face_index[0] = 0;
  face_index[1] = 5;

  for(int i = 0; i < nbr_face ; i++)
    cell_to_face_connectivity[i] = i+1;

  face_connectivity_index[0] = 0;
  face_connectivity_index[1] = 4;

  for (int i = 2; i < nbr_face+1 ; i++)
    face_connectivity_index[i] = 3*i+1;
     
  face_connectivity[0]  = 1;
  face_connectivity[1]  = 5;
  face_connectivity[2]  = 6;
  face_connectivity[3]  = 2;

  face_connectivity[4]  = 5;
  face_connectivity[5]  = 1;
  face_connectivity[6]  = 3;

  face_connectivity[7]  = 6;
  face_connectivity[8]  = 5;
  face_connectivity[9]  = 3;

  face_connectivity[10] = 2;
  face_connectivity[11] = 6;
  face_connectivity[12] = 3;

  face_connectivity[13] = 1;
  face_connectivity[14] = 2;
  face_connectivity[15] = 3;


   cwipi_add_polyhedra("c_mean_value_3D",
                     1,
                     face_index,
                     cell_to_face_connectivity,
                     face_connectivity_index,
                     face_connectivity);*/


  //ajouter des polyedres par cwipi_add_polyedra ...


  /* Points to locate
   * --------------- */

  int nPts = 19;               // Number of vertex
  double *coordsPts = (double *) malloc(nPts * 3 * sizeof(double));

  coordsPts[0] = 0.;
  coordsPts[1] = 0.;
  coordsPts[2] = 0.;
  
  /**** exemple tetraedre ****/

  /*coordsPts[3] = 1.;
  coordsPts[4] = 0.;
  coordsPts[5] = 0.;

  coordsPts[6] = 0.;
  coordsPts[7] = 1.;
  coordsPts[8] = 0.;
  
  coordsPts[9] =  0.;
  coordsPts[10] = 0.;
  coordsPts[11] = 1.;

  coordsPts[12] = 1./7;
  coordsPts[13] = 1./8;
  coordsPts[14] = 1./4;

  coordsPts[15] = 1./4;
  coordsPts[16] = 1./4;
  coordsPts[17] = 1./4;

  /***** exemple dodecaedre ****/

  /*coordsPts[3] = X;
  coordsPts[4] = 0.;
  coordsPts[5] = Z;

  coordsPts[6] = 0.;
  coordsPts[7] = -Z;
  coordsPts[8] = X   ;
  
  coordsPts[9] =  0.3;
  coordsPts[10] = 0.9;
  coordsPts[11] = 0.7;

  coordsPts[12] = 0.6;
  coordsPts[13] = -0.2;
  coordsPts[14] = -0.1;

  coordsPts[15] = 0.;
  coordsPts[16] = 0.;
  coordsPts[17] = 1.;

  coordsPts[18] = -1./8;
  coordsPts[19] = -1./3;
  coordsPts[20] = -1./2;

  coordsPts[21]= 0;
  coordsPts[22] =  Z/3 + X ;
  coordsPts[23] =  (X + 2*Z)/3 + 1.5*X;

  coordsPts[24]= 0;
  coordsPts[25] =  4*Z/9 + X/3 ;
  coordsPts[26] =  2./3 * Z + X /3 + coords[2]/3;*/

  /***** exemple avec deux octaedre ****/

  /* coordsPts[3] = 3./2;
  coordsPts[4] = 1./2;
  coordsPts[5] = 0.;

  coordsPts[6] = 3./2;
  coordsPts[7] = 1./2;
  coordsPts[8] = -1./4;
  
  coordsPts[9] =  1.2;
  coordsPts[10] = 0.9;
  coordsPts[11] = 3./4;

  coordsPts[12] = 1.;
  coordsPts[13] = 1.;
  coordsPts[14] = 1.;

  coordsPts[15] = 1.;
  coordsPts[16] = 1.;
  coordsPts[17] = -1.;

  coordsPts[18] = 2.;
  coordsPts[19] = 0.;
  coordsPts[20] = 0.;

  coordsPts[21] = 2.;
  coordsPts[22] = 1.;
  coordsPts[23] = 0.;

  coordsPts[24] = 2.;
  coordsPts[25] = 1.;
  coordsPts[26] = 1.;

  coordsPts[27] = 3./2-1./4;
  coordsPts[28] = 1./2+1./5;
  coordsPts[29] = -1./3;

  /**** octaedre non convexe ***/

  /* coordsPts[3] =  1./2;
  coordsPts[4] = 1./2;
  coordsPts[5] = 1./2;

  coordsPts[6] = 1./2;
  coordsPts[7] = 1./2;
  coordsPts[8] = 3./4;
  
  coordsPts[9] =  1./2;
  coordsPts[10] = 1./2;
  coordsPts[11] = 1.;

  coordsPts[12] = 1./2;
  coordsPts[13] = 1./2;
  coordsPts[14] = 0.;

  coordsPts[15] = 1./2;
  coordsPts[16] = 1./2;
  coordsPts[17] = 1./4;

  /****** Pyramide *****/

  coordsPts[3] =  0.;
  coordsPts[4] = 1.;
  coordsPts[5] = 0.;

  coordsPts[6] = 1.;
  coordsPts[7] = 0.;
  coordsPts[8] = 0.;
  
  coordsPts[9] =  1.;
  coordsPts[10] = 0.;
  coordsPts[11] = -2.;

  coordsPts[12] = 1./4;
  coordsPts[13] = 1./4;
  coordsPts[14] = 1./10;

  coordsPts[15] = 0.;
  coordsPts[16] = 1.;
  coordsPts[17] = -2.;

  coordsPts[18] = 1./2;
  coordsPts[19] = 1./2;
  coordsPts[20] = 1./2;

  coordsPts[21] = 1./2;
  coordsPts[22] = 1./3;
  coordsPts[23] = 1./4;

  coordsPts[24] = 1./2;
  coordsPts[25] = 1./3;
  coordsPts[26] = -1./4;

  coordsPts[27] = 1./2;
  coordsPts[28] = 1./4;
  coordsPts[29] = -3./2;

  coordsPts[30] = 1.8;
  coordsPts[31] = 1.3;
  coordsPts[32] = 0.2;

  coordsPts[33] = 1.5;
  coordsPts[34] = 0.;
  coordsPts[35] = -0.5;

  coordsPts[36] = 0.4;
  coordsPts[37] = 0.8;
  coordsPts[38] = -1.2;

  coordsPts[39] = 0.5;
  coordsPts[40] = 0.;
  coordsPts[41] = -1.5;

  coordsPts[42] = 0.9;
  coordsPts[43] = 0.6;
  coordsPts[44] = -0.7;

  coordsPts[45] = 1.8; 
  coordsPts[46] = 1.3; 
  coordsPts[47] = 0.2;

  coordsPts[48] = 1.5; 
  coordsPts[49] = 1. ; 
  coordsPts[50] = 1.7;
  
  coordsPts[51] = 1.65; 
  coordsPts[52] = 0.65 ;
  coordsPts[53] = -0.15;

  coordsPts[54] = 1.5; 
  coordsPts[55] = 1. ;
  coordsPts[56] = 1.5;

   cwipi_set_points_to_locate ("c_mean_value_3D", nPts, coordsPts);


  /*double* value = (double *) malloc (nVertex * sizeof(double));
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
                                                  &nNotLocatedPoints); */
  

   cwipi_locate("c_mean_value_3D");
       
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
