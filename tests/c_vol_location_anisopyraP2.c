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
#include "fvmc_ho_location.h"
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
  return 2*x*x + z*z -0.4*y*y + x*y - 3*x*z + 0.1*y*z + z - x + 2. + 3*z;
}

static double _x(double y)
{
  return 0.2*y*(y-1);
}

static double _y(double x)
{
  return 0.2*x*(x-1);
}

static double _z(double x)
{
  return 0.005*x*x + 0.75*x + 1;
}

/*------------------------------------------------------------------------
 *  Random
 *------------------------------------------------------------------------*/

static double frand_a_b(double a, double b){
    return (( rand()/(double)RAND_MAX ) * (b-a) + a);
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
  sprintf(fileName,"c_linear_location_pyraP2_%4.4d.txt",rank);

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

  fprintf(outputFile, "  Volume coupling test : location in pyrahedron P2\n");
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

  cwipi_create_coupling("c_volumic_cpl_location_pyraP2",            // Coupling id
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

  srand(time(NULL));
  double dila = frand_a_b(0.0001, 0.01);
  printf("dila = %12.5e\n", dila);

  const double xmin =  0.0;
  const double xmax =  5.0;
  const double ymin =  0.0;
  const double ymax =  5.0;
  const double zmin =  0.0;
  const double zmax =  1.0 * dila;

  nVertex = 14;
  nElts = 1;


  coords = (double *) malloc(sizeof(double) * 3 * nVertex );
  eltsConnecPointer = (int *) malloc(sizeof(int) * (nElts + 1));
  eltsConnec = (int *) malloc(sizeof(int) * nVertex);

  eltsConnecPointer[0] = 0;
  eltsConnecPointer[1] = 14;

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

  coords[12] = (xmin + xmax) / 2;
  coords[13] = (ymin + ymax) / 2;
  coords[14] = zmax;

  coords[15] = (xmin + xmax) / 2;
  coords[16] = _y((xmin + xmax) / 2) + ymin;
  coords[17] = zmin;

  coords[18] = _x((ymin + ymax) / 2) + xmin;
  coords[19] = (ymin + ymax) / 2;
  coords[20] = zmin;

  coords[21] = (xmin + xmax) / 2;
  coords[22] = (ymin + ymax) / 2;
  coords[23] = zmin;

  coords[24] = _x((ymin + ymax) / 2) + xmax;
  coords[25] = (ymin + ymax) / 2;
  coords[26] = zmin;

  coords[27] = (xmin + xmax) / 2;
  coords[28] = _y((xmin + xmax) / 2) + ymax;
  coords[29] = zmin;

  coords[30] = (coords[0] + coords[12]) / 2;
  coords[31] = (coords[1] + coords[13]) / 2;
  coords[32] = (coords[2] + coords[14]) / 2;

  coords[33] = (coords[3] + coords[12]) / 2;
  coords[34] = (coords[4] + coords[13]) / 2;
  coords[35] = (coords[5] + coords[14]) / 2;

  coords[36] = (coords[6] + coords[12]) / 2;
  coords[37] = (coords[7] + coords[13]) / 2;
  coords[38] = (coords[8] + coords[14]) / 2;

  coords[39] = (coords[9]  + coords[12]) / 2;
  coords[40] = (coords[10] + coords[13]) / 2;
  coords[41] = (coords[11] + coords[14]) / 2;


  if (rank == 0) {
    for (int g = 0; g < nVertex; g++) {
      printf("%12.5e, %12.5e, %12.5e \n", coords[3*g], coords[3*g+1], coords[3*g+2]);
    }
  }


  fprintf(outputFile, "   Number of vertex   : %i\n", nVertex);
  fprintf(outputFile, "   Number of elements : %i\n", nElts);

  const int order = 2;

  cwipi_ho_define_mesh("c_volumic_cpl_location_pyraP2",
                       nVertex,
                       nElts,
                       order,
                       coords,
                       eltsConnecPointer,
                       eltsConnec);

printf("\n" );

  const int n_node = 14;

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

  ijk[15] = 1;
  ijk[16] = 0;
  ijk[17] = 0;

  ijk[18] = 0;
  ijk[19] = 1;
  ijk[20] = 0;

  ijk[21] = 1;
  ijk[22] = 1;
  ijk[23] = 0;

  ijk[24] = 2;
  ijk[25] = 1;
  ijk[26] = 0;

  ijk[27] = 1;
  ijk[28] = 2;
  ijk[29] = 0;

  ijk[30] = 0;
  ijk[31] = 0;
  ijk[32] = 1;

  ijk[33] = 1;
  ijk[34] = 0;
  ijk[35] = 1;

  ijk[36] = 0;
  ijk[37] = 1;
  ijk[38] = 1;

  ijk[39] = 1;
  ijk[40] = 1;
  ijk[41] = 1;


  cwipi_ho_ordering_from_IJK_set ("c_volumic_cpl_location_pyraP2",
                                  CWIPI_CELL_PYRAMHO,
                                  n_node,
                                  ijk);



  int n_pts_to_locate = 14;

  double *pts_to_locate = (double *) malloc(sizeof(double) * 3 * n_pts_to_locate);

  for (int i=0; i<n_pts_to_locate; i++) {
    pts_to_locate[3*i]  = frand_a_b(xmin - (xmax-xmin)*0.1, xmax + (xmax-xmin)*0.1);
    pts_to_locate[3*i+1]  = _y(pts_to_locate[3*i]) + frand_a_b(ymin - (ymax-ymin)*0.1, ymax + (ymax-ymin)*0.1);
    pts_to_locate[3*i+2]  = _z(pts_to_locate[3*i]) + frand_a_b(zmin - (zmax-zmin)*0.1, zmax + (zmax-zmin)*0.1);
  }



  if (rank == 0) {
    for (int i = 0; i < n_pts_to_locate; i++) {
      printf("%12.5e %12.5e %12.5e\n",  pts_to_locate[3*i], pts_to_locate[3*i+1], pts_to_locate[3*i+2]);
    }
  }

  cwipi_set_points_to_locate ("c_volumic_cpl_location_pyraP2",
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

  // Define fields to send (X coordinate or Y coordinate)

  for (int i = 0; i < nVertex; i++) {
    sendValues[i] = _f(coords[3 * i], coords[3 * i+1], coords[3 * i+2]);
  }

  // Exchange

  int nNotLocatedPoints = 0;
  char *sendValuesName;
  char *recvValuesName;

  sendValuesName = "_fs";
  recvValuesName = "_fr";


    printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
  cwipi_locate("c_volumic_cpl_location_pyraP2");




  nNotLocatedPoints = cwipi_get_n_not_located_points("c_volumic_cpl_location_pyraP2");
  if (nNotLocatedPoints > 0) {
    printf("--- Error --- : %d not located points found\n", nNotLocatedPoints);
    exit(1);
  }

  int sRequest, rRequest;
  int tag = 1;

  cwipi_irecv("c_volumic_cpl_location_pyraP2",
              "ech",
              tag,
              1,
              1,
              0.1,
              recvValuesName,
              recvValues,
              &rRequest);


  cwipi_issend("c_volumic_cpl_location_pyraP2",
               "ech",
               tag,
               1,
               1,
               0.1,
               sendValuesName,
               sendValues,
               &sRequest);

  cwipi_wait_irecv("c_volumic_cpl_location_pyraP2", rRequest);
  cwipi_wait_issend("c_volumic_cpl_location_pyraP2", sRequest);



  // Coupling deletion

  if (rank == 0)
    printf("        Delete coupling\n");

  cwipi_delete_coupling("c_volumic_cpl_location_pyraP2");



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
