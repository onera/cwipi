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
  return 2*y*x + z*z - 3*x*z + z - x + 2. + 3*z;
}


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
  sprintf(fileName,"c_linear_location_tetraP2_%4.4d.txt",rank);

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

  fprintf(outputFile, "  Volume coupling test : location in tetrahedron P2\n");
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

  cwipi_create_coupling("c_volumic_cpl_location_tetraP2",            // Coupling id
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                        codeCoupledName,                           // Coupled application id
                        3,                                         // Geometric entities dimension
                        2,                                       // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        solver_type,                               // Solver type
                        postFreq,                                  // Postprocessing frequency
                        "EnSight Gold",                            // Postprocessing format
                        "text");                                   // Postprocessing option


  cwipi_ho_options_set("c_volumic_cpl_location_tetraP2",
                      "opt_bbox_step",
                      "-1");
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

  nVertex = 10;
  nElts = 1;


  coords = (double *) malloc(sizeof(double) * 3 * nVertex );
  eltsConnecPointer = (int *) malloc(sizeof(int) * (nElts + 1));
  eltsConnec = (int *) malloc(sizeof(int) * nVertex);

  eltsConnecPointer[0] = 0;
  eltsConnecPointer[1] = 10;

  eltsConnec[0] = 1;
  eltsConnec[1] = 2;
  eltsConnec[2] = 3;
  eltsConnec[3] = 4;
  eltsConnec[4] = 5;
  eltsConnec[5] = 6;
  eltsConnec[6] = 7;
  eltsConnec[7] = 8;
  eltsConnec[8] = 9;
  eltsConnec[9] = 10;

/*
  coords[0] = -4.903926402015830e-01;//xmin;
  coords[1] =  9.754516100822440e-02;//ymin;
  coords[2] =  0.000000000000000e+00;//zmin;

  coords[3] = -4.789922246525480e-01;//xmax;
  coords[4] =  1.392031233369340e-01;//ymin;
  coords[5] =  3.448099731221840e-02;//zmin;

  coords[6] = -4.615116344320200e-01;//xmin;
  coords[7] =  1.797641367398650e-01;//ymax;
  coords[8] =  6.849720013297110e-02;//zmin;

  coords[9]  = -4.205148508640080e-01;//xmin;
  coords[10] =  6.032984012355020e-02;//ymin;
  coords[11] =  3.039338835464140e-02;//zmax;

  coords[12] = -4.060743479792260e-01;//(xmin + xmax) / 2;
  coords[13] =  1.014393279893710e-01;//ymin;
  coords[14] =  6.464198842112701e-02;//zmin;

  coords[15] = -3.506370615264330e-01;//(xmin + xmax) / 2;
  coords[16] =  2.311451923887600e-02;//(ymin + xmax) / 2;
  coords[17] =  6.078677670928280e-02;//zmin;

  coords[18] = -4.174275540611020e-01;//xmin;
  coords[19] =  1.117825698303370e-01;//(ymin + ymax) / 2;
  coords[20] = -1.540283024855780e-02;//zmin;

  coords[21] = -4.029870511763210e-01;//xmin;
  coords[22] =  1.528920576961570e-01;//ymin;
  coords[23] =  1.884576981792770e-02;//(zmin + zmax) / 2;

  coords[24] = -3.475497647235270e-01;//(xmin + xmax) / 2;
  coords[25] =  7.456724894566230e-02;//ymin;
  coords[26] =  1.499055810608360e-02;//(zmin + zmax) / 2;

  coords[27] = -3.444624679206210e-01;//xmin;
  coords[28] =  1.260199786524490e-01;//(ymin + ymax) / 2;
  coords[29] = -3.080566049711560e-02;//(zmin + zmax) / 2;
*/



coords[0] = 9.030342459149850e-02;
coords[1] = 2.423433238303810e-01;
coords[2] = 4.279193906589080e-01;

coords[3] = 3.513641864412240e-02;
coords[4] = 1.996536901834650e-01;
coords[5] = 3.618983441488240e-01;

coords[6] = -2.003058730325370e-02;
coords[7] = 1.569640565365490e-01;
coords[8] = 2.958772976387400e-01;

coords[9] = 4.549623479940860e-02;
coords[10] = 2.612693950859740e-01;
coords[11] = 4.238730892737660e-01;

coords[12] = -1.001529365162680e-02;
coords[13] = 2.173745865234740e-01;
coords[14] = 3.558060518948060e-01;

coords[15] = 0.000000000000000e+00;
coords[16] = 2.777851165103990e-01;
coords[17] = 4.157348061508730e-01;

 coords[18] = 8.776630908293880e-02;
 coords[19] = 2.488799519208400e-01;
 coords[20] = 3.484222722132320e-01;

 coords[21] = 3.259930313556270e-02;
 coords[22] = 2.061903182739240e-01;
 coords[23] = 2.824012257031480e-01;

 coords[24] = 4.261459678718950e-02;
 coords[25] = 2.666008482608490e-01;
 coords[26] = 3.423299799592140e-01;

 coords[27] = 8.522919357437909e-02;
 coords[28] = 2.554165800112990e-01;
 coords[29] = 2.689251537675560e-01;

  fprintf(outputFile, "   Number of vertex   : %i\n", nVertex);
  fprintf(outputFile, "   Number of elements : %i\n", nElts);

  const int order = 2;

  cwipi_ho_define_mesh("c_volumic_cpl_location_tetraP2",
                       nVertex,
                       nElts,
                       order,
                       coords,
                       eltsConnecPointer,
                       eltsConnec);



  const int n_node = 10;

  int *ijk = malloc(sizeof(int)*3*n_node);

  ijk[ 0] = 0;
  ijk[ 1] = 0;
  ijk[ 2] = 0;

  ijk[ 3] = 1;//2;
  ijk[ 4] = 0;
  ijk[ 5] = 0;

  ijk[ 6] = 2;//0;
  ijk[ 7] = 0;//2;
  ijk[ 8] = 0;

  ijk[ 9] = 0;
  ijk[10] = 1;//0;
  ijk[11] = 0;//2;

  ijk[12] = 1;
  ijk[13] = 1;//0;
  ijk[14] = 0;

  ijk[15] = 0;//1;
  ijk[16] = 2;//1;
  ijk[17] = 0;

  ijk[18] = 0;
  ijk[19] = 0;//1;
  ijk[20] = 1;//0;

  ijk[21] = 1;//0;
  ijk[22] = 0;
  ijk[23] = 1;

  ijk[24] = 0;//1;
  ijk[25] = 1;//0;
  ijk[26] = 1;

  ijk[27] = 0;
  ijk[28] = 0;//1;
  ijk[29] = 2;//1;

  cwipi_ho_ordering_from_IJK_set ("c_volumic_cpl_location_tetraP2",
                                  CWIPI_CELL_TETRAHO,
                                  n_node,
                                  ijk);



  int n_pts_to_locate = 1;

  double *pts_to_locate = (double *) malloc(sizeof(double) * 3 * n_pts_to_locate);


  pts_to_locate[0]  = -8.407698607790286e-04;
  pts_to_locate[1]  =  1.964277933020044e-01;
  pts_to_locate[2]  =  3.223073602567958e-01;

/*  for (int i = 0; i < n_pts_to_locate; i++){
    pts_to_locate[3*i]    =  frand_a_b(xmin,0.3);
    pts_to_locate[3*i+1]  =  frand_a_b(ymin,0.3);
    pts_to_locate[3*i+2]  =  frand_a_b(zmin,0.3);
  }*/


  for (int i = 0; i < n_pts_to_locate; i++) {
    printf("point : %12.15e %12.15e %12.15e\n",  pts_to_locate[3*i], pts_to_locate[3*i+1], pts_to_locate[3*i+2]);
  }

  cwipi_set_points_to_locate ("c_volumic_cpl_location_tetraP2",
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



  cwipi_locate("c_volumic_cpl_location_tetraP2");




  nNotLocatedPoints = cwipi_get_n_not_located_points("c_volumic_cpl_location_tetraP2");
  if (nNotLocatedPoints > 0) {
    printf("--- Error --- : %d not located points found\n", nNotLocatedPoints);
    exit(1);
  }

  int sRequest, rRequest;
  int tag = 1;

  cwipi_irecv("c_volumic_cpl_location_tetraP2",
              "ech",
              tag,
              1,
              1,
              0.1,
              recvValuesName,
              recvValues,
              &rRequest);


  cwipi_issend("c_volumic_cpl_location_tetraP2",
               "ech",
               tag,
               1,
               1,
               0.1,
               sendValuesName,
               sendValues,
               &sRequest);

  cwipi_wait_irecv("c_volumic_cpl_location_tetraP2", rRequest);
  cwipi_wait_issend("c_volumic_cpl_location_tetraP2", sRequest);



  /* Coupling deletion
   * ----------------- */

  if (rank == 0)
    printf("        Delete coupling\n");

  cwipi_delete_coupling("c_volumic_cpl_location_tetraP2");



  /* Check barycentric coordinates */

  if (rank == 0)
    printf("        Check results\n");

  double *res = (double *) malloc(sizeof(double) *  n_pts_to_locate);

  for (int i = 0; i < n_pts_to_locate; i++) {
    res[i] = _f(pts_to_locate[3*i], pts_to_locate[3*i+1], pts_to_locate[3*i+2]);
  }

  double err;

  for (int i = 0; i < n_pts_to_locate; i++) {
    err = fabs(recvValues[i] - res[i]);
    //    if (err > 1e-6) {
    printf ("[%d] err %d : %12.5e %12.5e %12.5e\n", codeId, i, err, recvValues[i], res[i]);
      // }
  }

  double err_max;
  MPI_Allreduce(&err, &err_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if (err_max >= 1e-6) {
    if (rank == 0) {
      printf("        !!! Error = %12.5e\n", err_max);
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }




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
