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
     "  -visu           Ensight outputs \n\n"
     "  -a              Unlocking communication \n\n"
     "  -stdout         Standard output \n\n"
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
           double       *randLevel,
	   int          *randFromClock,
	   int          *postFreq,
	   int          *t_com,
	   int          *tostdout)

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
    else if (strcmp(argv[i], "-randFromClock") == 0) {
      *randFromClock = 1;
    }
    else if (strcmp(argv[i], "-a") == 0) {
      *t_com = 1;
    }
    else if (strcmp(argv[i], "-visu") == 0) {
      *postFreq = 1;
    }
    else if (strcmp(argv[i], "-stdout") == 0) {
      *tostdout = 1;
    }
    else
      _usage(EXIT_FAILURE);
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

  int n2 = (int) (two * pow(n_partition, two));

  if (n2 != commWorldSize) {
    if (rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  /* Read args from command line
   * --------------------------- */

/*   int nVertexSeg = 10; */
/*   double randLevel = 0.4; */

  int nVertexSeg = 10;
  double randLevel = 0.4;
  int postFreq     = -1;
  int t_com        = 0;
  int tostdout     = 0;
  int randFromClock = 0;

  _read_args(argc, argv, &nVertexSeg, &randLevel,
             &randFromClock, &postFreq, &t_com, &tostdout);

  if (randFromClock == 1)
    srand(rank + time(0));
  else
    srand(rank + 1);

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

  char* fileName = (char *) malloc(sizeof(char) * 30);
  sprintf(fileName,"c_surf_coupling_P1P1_%4.4d.txt",rank);

  if (tostdout)
    outputFile = stdout;
  else
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

  fprintf(outputFile, "  Surface coupling test : P1P1 with polygon\n");
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

  cwipi_create_coupling("c_surf_cpl_P1P1",                                // Coupling id
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                        codeCoupledName,                           // Coupled application id
                        2,                                         // Geometric entities dimension
                        1e-3,                                       // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        solver_type,                               // Solver type
                        postFreq,                                         // Postprocessing frequency
                        "EnSight Gold",                            // Postprocessing format
                        "text");                                   // Postprocessing option
  
  /* Mesh definition
   * --------------- */

  if (rank == 0)
    printf("        Create mesh\n");

  int nVertex = 0;               // Number of vertex
  double *coords = NULL;         // Vertex coordinates
  int nElts = 0;                 // Number of elements
  int *eltsConnecPointer = NULL; // Connectivity index
  int *eltsConnec = NULL;        // Connectivity
  
  /* Domain bounds */

  const double xmin = -10;
  const double xmax =  10;
  const double ymin = -10;
  const double ymax =  10;

  nVertex = nVertexSeg * nVertexSeg;
  nElts = (nVertexSeg - 1) * (nVertexSeg - 1);

  coords = (double *) malloc(sizeof(double) * 3 * nVertex );
  eltsConnecPointer = (int *) malloc(sizeof(int) * (nElts + 1));
  eltsConnec = (int *) malloc(sizeof(int) * 4 * nElts);
  
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
            localComm); 


  fprintf(outputFile, "   Number of vertex   : %i\n", nVertex);
  fprintf(outputFile, "   Number of elements : %i\n", nElts);

  cwipi_define_mesh("c_surf_cpl_P1P1",
                    nVertex,
                    nElts,
                    coords,
                    eltsConnecPointer,
                    eltsConnec);

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
  recvValues = (double *) malloc(sizeof(double) * nVertex);

  /* Define fields to send (X coordinate or Y coordinate) */

  for (int i = 0; i < nVertex; i++) {
    if (codeId == 1)
      sendValues[i] = coords[3 * i];
    else
      sendValues[i] = coords[3 * i + 1];
  }

  /* Exchange */

  int nNotLocatedPoints = 0;
  char *sendValuesName;
  char *recvValuesName;
  if (codeId == 1) {
    sendValuesName = "cooX";
    recvValuesName = "cooY";
  }
  else{
    sendValuesName = "cooY";
    recvValuesName = "cooX";
  }

  cwipi_locate("c_surf_cpl_P1P1");

  nNotLocatedPoints = cwipi_get_n_not_located_points("c_surf_cpl_P1P1");
  if (nNotLocatedPoints > 0) {
    printf("--- Error --- : %d not located points found\n", nNotLocatedPoints);
    exit(1);
  }

  if (t_com == 0) {
    cwipi_exchange_status_t status = cwipi_exchange("c_surf_cpl_P1P1",
						    "ech",
						    1,
						    1,     // n_step
						    0.1,   // physical_time
						    sendValuesName,
						    sendValues,
						    recvValuesName,
						    recvValues,
						    &nNotLocatedPoints);
    _dumpStatus(outputFile, status);
  }

  else if (t_com == 1) {

    int sRequest, rRequest;
    int tag = 1;

    cwipi_irecv("c_surf_cpl_P1P1",
		"ech",
		tag,
		1,
		1,
		0.1,
		recvValuesName,
		recvValues,
		&rRequest);
  

    cwipi_issend("c_surf_cpl_P1P1",
		 "ech",
		 tag,
		 1,
		 1,
		 0.1,
		 sendValuesName,
		 sendValues,
		 &sRequest);

    cwipi_wait_irecv("c_surf_cpl_P1P1", rRequest);
    cwipi_wait_issend("c_surf_cpl_P1P1", sRequest);

  }


  /* Coupling deletion
   * ----------------- */

  if (rank == 0)
    printf("        Delete coupling\n");

  cwipi_delete_coupling("c_surf_cpl_P1P1");

  /* Check barycentric coordinates */

  if (rank == 0)
    printf("        Check results\n");    

  double err;
  if (codeId == 1)
    err = fabs(recvValues[0] - coords[3 * 0 + 1]);
  else
    err = fabs(recvValues[0] - coords[3 * 0    ]);
 
  for (int i = 0; i < nVertex; i++) {
    if (codeId == 1) {
      err = ((fabs(recvValues[i] - coords[3 * i + 1])) < (err) ? (err) : 
             (fabs(recvValues[i] - coords[3 * i + 1])));
      if (fabs(recvValues[i] - coords[3 * i + 1]) > 1e-6) {
        printf ("[%d] err %d : %12.5e\n", codeId, i, err);
      }
    }
    else {
      err = ((fabs(recvValues[i] - coords[3 * i    ])) < (err) ? (err) : 
             (fabs(recvValues[i] - coords[3 * i    ])));
      if (fabs(recvValues[i] - coords[3 * i    ]) > 1e-6) {
        printf ("[%d] err %d : %12.5e %12.5e\n", codeId, i, err, fabs(recvValues[i] - coords[3 * i    ]));
      }
    }
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

  if (!tostdout)
    fclose(outputFile);

  return EXIT_SUCCESS;
}