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
#include <assert.h>

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

/* static void _dumpStatus(FILE* outputFile, cwipi_exchange_status_t status) */
/* { */
/*   switch(status) { */
/*   case CWIPI_EXCHANGE_OK : */
/*     fprintf(outputFile, "Exchange Ok\n"); */
/*     break; */
/*   case CWIPI_EXCHANGE_BAD_RECEIVING : */
/*     fprintf(outputFile, "Bad receiving\n"); */
/*     break; */
/*   default : */
/*     printf("Error : bad exchange status\n"); */
/*     exit(1); */
/*   } */
/* } */

/*----------------------------------------------------------------------
 *
 *
 *
 * parameters:
 *   status              <-- Exchange status
 *---------------------------------------------------------------------*/

static double _f(double x, double y, double z)
{
  return 2*x*x + z*z - 3*x*z + z - x + 2. + 3*z;
}

static double _y(double x)
{
  return x*x + 2*x -1;
}

static double _z(double x)
{
  return x*x + 2;
}


/*----------------------------------------------------------------------
 *
 * Read mesh
 *
 * parameters:
 *   f                   <-- Mesh file
 *   dimension           --> Dimension
 *   nvertex             <-- number of vertices
 *   nElements           <-- number of elements
 *   nConnecVertex       <-- size of connectivity
 *   coords              --> vertices coordinates
 *   connecPointer       --> connectivity index
 *   connec              --> connectivity
 *---------------------------------------------------------------------*/

static int _read_mesh(FILE *f,
                      int *format,
                      int *dimension,
                      int *nVertex,
                      int *nElt,
                      double *coords,
                      int *eltsConnecPointer,
                      int *eltsConnec)
{


  int r;
  int nConnecVertex;
  int _format;
  int _dimension;
  int _nVertex;
  int _nElt;
  int *un, loop = 0;
  char key[256];



  while (loop == 0){
    r = fscanf(f, "%s",key);
    printf("key = %s\n", key);
  switch (key[0]) {
    case 'M':
      r = fscanf(f, "%d",format);
      _format = *format;
      printf("format = %d, r = %i\n", _format, r);
      break;

    case 'D':
      r = fscanf(f, "%d",dimension);
      _dimension = *dimension;
      printf("dimension = %d, r = %i\n", _dimension, r);
      break;

    case 'V':
      r = fscanf(f, "%d",nVertex);
      _nVertex = *nVertex;
      printf("nVertex = %d, r = %i\n", _nVertex, r);
      coords = (double *) malloc(sizeof(double) * 3 * _nVertex );
      for (int i = 0; i < _nVertex; i++) {

        for (int j = 0; j < 3; j++) {
        r = fscanf(f, "%lf",coords + i * 3 + j);
        }
      }
      break;

      case 'E':
        r = fscanf(f, "%d",nElt);
        _nElt = *nElt;
        printf("nElt = %d, r = %i\n", _nElt, r);
        nConnecVertex = _nElt * 3;
        eltsConnec = (int *) malloc(sizeof(int) * nConnecVertex);
        for (int i = 0; i < nConnecVertex; i++) {

          for (int j = 0; j < 3; j++) {
          r = fscanf(f, "%d",eltsConnec + i * 3 + j);
          }
          r = fscanf(f, "%d",un);
        }
        break;

      case 'F':
        loop = 1;
        break;
  };
}

  eltsConnecPointer = (int *) malloc(sizeof(int) * (_nElt + 1));

  for (int i = 0; i < _nElt; i++) {
    eltsConnecPointer[i] = 3*i;
  }
  eltsConnecPointer[_nElt] = nConnecVertex;

  return 1;
}
/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

/* static void */
/* _usage(int exit_code) */
/* { */
/*   printf */
/*     ("\n" */
/*      "  Usage: \n\n" */
/*      "  -n     <level>  Number of vertices in band width.\n\n" */
/*      "  -rand  <level>  Random level ( > 0 and < 0.4) \n\n" */
/*      "  -visu           Ensight outputs \n\n" */
/*      "  -a              Unlocking communication \n\n" */
/*      "  -stdout         Standard output \n\n" */
/*      "  -h             this message.\n\n"); */

/*   exit(exit_code); */
/* } */

/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nVertex             <-- Number of vertices in bandwidth
 *   randLevel           <-- Random level
 *---------------------------------------------------------------------*/

/* static void */
/* _read_args(int            argc, */
/*            char         **argv, */
/*            int          *nVertex, */
/*            double       *randLevel, */
/* 	   int          *randFromClock, */
/* 	   int          *postFreq, */
/* 	   int          *t_com, */
/* 	   int          *tostdout) */

/* { */
/*   int i = 1; */

/*   /\* Parse and check command line *\/ */

/*   while (i < argc) { */

/*     if (strcmp(argv[i], "-h") == 0) */
/*       _usage(EXIT_SUCCESS); */

/*     else if (strcmp(argv[i], "-n") == 0) { */
/*       i++; */
/*       if (i >= argc) */
/*         _usage(EXIT_FAILURE); */
/*       else */
/*         *nVertex = atoi(argv[i]); */
/*     } */
/*     else if (strcmp(argv[i], "-rand") == 0) { */
/*       i++; */
/*       if (i >= argc) */
/*         _usage(EXIT_FAILURE); */
/*       else */
/*         *randLevel = atof(argv[i]); */
/*     } */
/*     else if (strcmp(argv[i], "-randFromClock") == 0) { */
/*       *randFromClock = 1; */
/*     } */
/*     else if (strcmp(argv[i], "-a") == 0) { */
/*       *t_com = 1; */
/*     } */
/*     else if (strcmp(argv[i], "-visu") == 0) { */
/*       *postFreq = 1; */
/*     } */
/*     else if (strcmp(argv[i], "-stdout") == 0) { */
/*       *tostdout = 1; */
/*     } */
/*     else */
/*       _usage(EXIT_FAILURE); */
/*     i++; */
/*   } */
/* } */


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
  FILE* meshFile;

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

  srand(rank + time(0));

  int n_partition = 0;
  const int two = 2;
  while(two * pow(n_partition, two) < commWorldSize) n_partition++;

  if (two != commWorldSize) {
    if (rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  /* Initialization
   * -------------- */

  char *codeName;
  int codeId;
  char *codeCoupledName;

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

  char* fileName = (char *) malloc(sizeof(char) * 37);
  sprintf(fileName,"c_linear_location_hexaP2_%4.4d.txt",rank);

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

  fprintf(outputFile, "  Volume coupling test : location in hexahedron P2\n");
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

  const int postFreq = -1;

  cwipi_create_coupling("c_volumic_cpl_location_hexaP2",            // Coupling id
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                        codeCoupledName,                           // Coupled application id
                        3,                                         // Geometric entities dimension
                        0.1,                                       // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        solver_type,                               // Solver type
                        postFreq,                                  // Postprocessing frequency
                        "EnSight Gold",                            // Postprocessing format
                        "text");                                   // Postprocessing option

  /* Mesh definition
   * --------------- */

  if (rank == 0)
    printf("        Create mesh\n");

  int format = 0;
  int dimension = 0;
  int nVertex = 0;               // Number of vertex
  double *coords = NULL;         // Vertex coordinates
  int nElts = 0;                 // Number of elements
  int *eltsConnecPointer = NULL; // Connectivity index
  int *eltsConnec = NULL;        // Connectivity

  /* Domain bounds */

  const double xmin =  0.0;
  const double xmax =  1.0;
  const double ymin =  0.0;
  const double ymax =  1.0;
  const double zmin =  0.0;
  const double zmax =  1.0;

  nVertex = 27;
  nElts = 1;

/*
 meshFile = fopen("meshes/hexahedronp2.mesh", "r");

  assert (meshFile != NULL);
//  _read_mesh(meshFile, &format, &dimension, &nVertex, &nElts, coords, eltsConnecPointer, eltsConnec);
int r;
int nConnecVertex;
int _format;
int _dimension;
int _nVertex;
int _nElts;
int *un, loop = 0;
char key[40];



while (loop == 0){
  r = fscanf(meshFile, "%s",key);
  printf("key = %s\n", key);
switch (key[0]) {
  case 'M':
    r = fscanf(meshFile, "%d",&format);
    _format = format;
    printf("format = %d, r = %i\n", _format, r);
    break;

  case 'D':
    r = fscanf(meshFile, "%d",&dimension);
    _dimension = dimension;
    printf("dimension = %d, r = %i\n", _dimension, r);
    break;

  case 'V':
    r = fscanf(meshFile, "%d",&nVertex);
    _nVertex = nVertex;
    printf("nVertex = %d, r = %i\n", _nVertex, r);
    coords = (double *) malloc(sizeof(double) * 3 * _nVertex );
    for (int i = 0; i < _nVertex; i++) {

      r = fscanf(meshFile, "%lf, %lf, %lf, %d",coords + i * 3, coords + i * 3 + 1, coords + i * 3 + 2, un);
    }
    break;

    case 'T':
      //switch (key[1]) {
        //case 'd':
          r = fscanf(meshFile, "%d",&nElts);
          _nElts = nElts;
          printf("nElts = %d, r = %i\n", _nElts, r);
          nConnecVertex = _nElts * 3;
          eltsConnec = (int *) malloc(sizeof(int) * 11);
          for (int i = 0; i < 11; i++) {

              r = fscanf(meshFile, "%d",eltsConnec + i );
          }
          r = fscanf(meshFile, "%d",un);
          break;
        //case 'n':
    case 'E':
          loop = 1;
          break;

    }
}

eltsConnecPointer = (int *) malloc(sizeof(int) * (_nElts + 1));

for (int i = 0; i < _nElts; i++) {
  eltsConnecPointer[i] = 3*i;
}
eltsConnecPointer[_nElts] = nConnecVertex;
*/


  coords = (double *) malloc(sizeof(double) * 3 * nVertex );
  eltsConnecPointer = (int *) malloc(sizeof(int) * (nElts + 1));
  eltsConnec = (int *) malloc(sizeof(int) * nVertex);

  eltsConnecPointer[0] = 0;
  eltsConnecPointer[1] = 27;

  eltsConnec[0]  = 1;
  eltsConnec[1]  = 2;
  eltsConnec[2]  = 3;
  eltsConnec[3]  = 4;
  eltsConnec[4]  = 5;
  eltsConnec[5]  = 6;
  eltsConnec[6]  = 7;
  eltsConnec[7]  = 8;
  eltsConnec[8]  = 9;
  eltsConnec[9]  = 10;
  eltsConnec[10] = 11;
  eltsConnec[11] = 12;
  eltsConnec[12] = 13;
  eltsConnec[13] = 14;
  eltsConnec[14] = 15;
  eltsConnec[15] = 16;
  eltsConnec[16] = 17;
  eltsConnec[17] = 18;
  eltsConnec[18] = 19;
  eltsConnec[19] = 20;
  eltsConnec[20] = 21;
  eltsConnec[21] = 22;
  eltsConnec[22] = 23;
  eltsConnec[23] = 24;
  eltsConnec[24] = 25;
  eltsConnec[25] = 26;
  eltsConnec[26] = 27;

  coords[0] = xmin;
  coords[1] = ymin;
  coords[2] = zmin;

  coords[3] = xmax;
  coords[4] = ymin;
  coords[5] = zmin;

  coords[6] = xmin;
  coords[7] = ymax;
  coords[8] = zmin;

  coords[9]  = xmax;
  coords[10] = ymax;
  coords[11] = zmin;

  coords[12] = xmin;
  coords[13] = ymin;
  coords[14] = zmax;

  coords[15] = xmax;
  coords[16] = ymin;
  coords[17] = zmax;

  coords[18] = xmin;
  coords[19] = ymax;
  coords[20] = zmax;

  coords[21] = xmax;
  coords[22] = ymax;
  coords[23] = zmax;

  coords[24] = (xmin + xmax) / 2;
  coords[25] = ymin + 0.1;
  coords[26] = zmin + 0.1;


  coords[27] = xmin + 0.1;
  coords[28] = (ymin + ymax) / 2;
  coords[29] = zmin + 0.1;

  coords[30] = (xmin + xmax) / 2;
  coords[31] = (ymin + ymax) / 2;
  coords[32] = zmin + 0.2;

  coords[33] = xmax - 0.1;
  coords[34] = (ymin + ymax) / 2;
  coords[35] = zmin + 0.1;

  coords[36] = (xmin + xmax) / 2;
  coords[37] = ymax - 0.1;
  coords[38] = zmin + 0.1;

  coords[39] = xmin + 0.1;
  coords[40] = ymin + 0.1;
  coords[41] = (zmin + zmax) / 2;

  coords[42] = (xmin + xmax) / 2;
  coords[43] = ymin + 0.2;
  coords[44] = (zmin + zmax) / 2;

  coords[45] = xmax - 0.1;
  coords[46] = ymin + 0.1;
  coords[47] = (zmin + zmax) / 2;

  coords[48] = xmin + 0.2;
  coords[49] = (ymin + ymax) / 2;
  coords[50] = (zmin + zmax) / 2;

  coords[51] = (xmin + xmax) / 2;
  coords[52] = (ymin + ymax) / 2;
  coords[53] = (zmin + zmax) / 2;

  coords[54] = xmax - 0.2;
  coords[55] = (ymin + ymax) / 2;
  coords[56] = (zmin + zmax) / 2;

  coords[57] = xmin + 0.1;
  coords[58] = ymax - 0.1;
  coords[59] = (zmin + zmax) / 2;

  coords[60] = (xmin + xmax) / 2;
  coords[61] = ymax - 0.2;
  coords[62] = (zmin + zmax) / 2;

  coords[63] = xmax - 0.1;
  coords[64] = ymax - 0.1;
  coords[65] = (zmin + zmax) / 2;

  coords[66] = (xmin + xmax) / 2;
  coords[67] = ymin + 0.1;
  coords[68] = zmax - 0.1;

  coords[69] = xmin + 0.1;
  coords[70] = (ymin + ymax) / 2;
  coords[71] = zmax - 0.1;

  coords[72] = (xmin + xmax) / 2;
  coords[73] = (ymin + ymax) / 2;
  coords[74] = zmax - 0.2;

  coords[75] = xmax - 0.1;
  coords[76] = (ymin + ymax) / 2;
  coords[77] = zmax - 0.1;

  coords[78] = (xmin + xmax) / 2;
  coords[79] = ymax - 0.1;
  coords[80] = zmax - 0.1;


  fprintf(outputFile, "   Number of vertex   : %i\n", nVertex);
  fprintf(outputFile, "   Number of elements : %i\n", nElts);

  const int order = 2;

  cwipi_ho_define_mesh("c_volumic_cpl_location_hexaP2",
                       nVertex,
                       nElts,
                       order,
                       coords,
                       eltsConnecPointer,
                       eltsConnec);



  const int n_node = 27;

  int *ijk = malloc(sizeof(int)*3*n_node);

  ijk[ 0] = 0;
  ijk[ 1] = 0;
  ijk[ 2] = 0;

  ijk[ 3] = 2;
  ijk[ 4] = 0;
  ijk[ 5] = 0;

  ijk[ 6] = 0;
  ijk[ 7] = 2;
  ijk[ 8] = 0;

  ijk[ 9] = 2;
  ijk[10] = 2;
  ijk[11] = 0;

  ijk[12] = 0;
  ijk[13] = 0;
  ijk[14] = 2;

  ijk[15] = 2;
  ijk[16] = 0;
  ijk[17] = 2;

  ijk[18] = 0;
  ijk[19] = 2;
  ijk[20] = 2;

  ijk[21] = 2;
  ijk[22] = 2;
  ijk[23] = 2;

  ijk[24] = 1;
  ijk[25] = 0;
  ijk[26] = 0;

  ijk[27] = 0;
  ijk[28] = 1;
  ijk[29] = 0;

  ijk[30] = 1;
  ijk[31] = 1;
  ijk[32] = 0;

  ijk[33] = 2;
  ijk[34] = 1;
  ijk[35] = 0;

  ijk[36] = 1;
  ijk[37] = 2;
  ijk[38] = 0;

  ijk[39] = 0;
  ijk[40] = 0;
  ijk[41] = 1;

  ijk[42] = 1;
  ijk[43] = 0;
  ijk[44] = 1;

  ijk[45] = 2;
  ijk[46] = 0;
  ijk[47] = 1;

  ijk[48] = 0;
  ijk[49] = 1;
  ijk[50] = 1;

  ijk[51] = 1;
  ijk[52] = 1;
  ijk[53] = 1;

  ijk[54] = 2;
  ijk[55] = 1;
  ijk[56] = 1;

  ijk[57] = 0;
  ijk[58] = 2;
  ijk[59] = 1;

  ijk[60] = 1;
  ijk[61] = 2;
  ijk[62] = 1;

  ijk[63] = 2;
  ijk[64] = 2;
  ijk[65] = 1;

  ijk[66] = 1;
  ijk[67] = 0;
  ijk[68] = 2;

  ijk[69] = 0;
  ijk[70] = 1;
  ijk[71] = 2;

  ijk[72] = 1;
  ijk[73] = 1;
  ijk[74] = 2;

  ijk[75] = 2;
  ijk[76] = 1;
  ijk[77] = 2;

  ijk[78] = 1;
  ijk[79] = 2;
  ijk[80] = 2;

  cwipi_ho_ordering_from_IJK_set ("c_volumic_cpl_location_hexaP2",
                                  CWIPI_CELL_HEXAHO,
                                  n_node,
                                  ijk);



  int n_pts_to_locate = 1;

  double *pts_to_locate = (double *) malloc(sizeof(double) * 3 * n_pts_to_locate);

  pts_to_locate[0]  = 0.5;
  pts_to_locate[1]  = 0.5;
  pts_to_locate[2]  = 0.1;



  for (int i = 0; i < n_pts_to_locate; i++) {
    printf("%12.5e %12.5e %12.5e\n",  pts_to_locate[3*i], pts_to_locate[3*i+1], pts_to_locate[3*i+2]);
  }

  cwipi_set_points_to_locate ("c_volumic_cpl_location_hexaP2",
                              n_pts_to_locate,
                              pts_to_locate);

  /* Fields exchange
   *     - Proc 0 : Send X coordinates
   *                Recv Y coordinates
   *     - Proc 1 : Send Y coordinates
   *                Recv X coordinates
   * --------------------------------- */

  if (rank == 0)
    printf("        Exchange Code1 <-> Code2\n");

  double *sendValues = NULL;
  double *recvValues = NULL;

  sendValues = (double *) malloc(sizeof(double) * nVertex);
  recvValues = (double *) malloc(sizeof(double) * n_pts_to_locate);

  /* Define fields to send (X coordinate or Y coordinate) */

  for (int i = 0; i < nVertex; i++) {
    sendValues[i] = _f(coords[3 * i], coords[3 * i+1], coords[3 * i+2]);
  }

  /* Exchange */

  int nNotLocatedPoints = 0;
  char *sendValuesName;
  char *recvValuesName;

  sendValuesName = "_fs";
  recvValuesName = "_fr";


    printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");

  cwipi_locate("c_volumic_cpl_location_hexaP2");




  nNotLocatedPoints = cwipi_get_n_not_located_points("c_volumic_cpl_location_hexaP2");
  if (nNotLocatedPoints > 0) {
    printf("--- Error --- : %d not located points found\n", nNotLocatedPoints);
    exit(1);
  }

  int sRequest, rRequest;
  int tag = 1;

  cwipi_irecv("c_volumic_cpl_location_hexaP2",
              "ech",
              tag,
              1,
              1,
              0.1,
              recvValuesName,
              recvValues,
              &rRequest);


  cwipi_issend("c_volumic_cpl_location_hexaP2",
               "ech",
               tag,
               1,
               1,
               0.1,
               sendValuesName,
               sendValues,
               &sRequest);

  cwipi_wait_irecv("c_volumic_cpl_location_hexaP2", rRequest);
  cwipi_wait_issend("c_volumic_cpl_location_hexaP2", sRequest);



  /* Coupling deletion
   * ----------------- */

  if (rank == 0)
    printf("        Delete coupling\n");

  cwipi_delete_coupling("c_volumic_cpl_location_hexaP2");



  /* Free
   * ---- */

  free(coords);
  free(eltsConnecPointer);
  free(eltsConnec);
  free(sendValues);
  free(recvValues);
  free(srcName);

  /* Finalize
   * -------- */

  cwipi_finalize();

  MPI_Finalize();




  fclose(outputFile);
//fclose(meshFile);

  return EXIT_SUCCESS;
}
