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

#include <mpi.h>

#include "creeMaillagePolygone2D.h"
#include "cwipi.h"

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
 * Main : surface coupling test : P1P1 with polygon                            
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
  char *srcBaseName = strrchr(srcName, '.');
  *srcBaseName = '\0';
  srcBaseName = strrchr(srcName, '/') + 1;

  if (rank == 0)
    printf("\nSTART: %s\n", srcBaseName);

  if (comm_world_size != 2) {
    if (rank == 0)
      printf("      Not executed : only available for 2 processus\n");
    MPI_Finalize();
    return 0;
  }

  /* Initialization
   * -------------- */

  char *codeName;
  char *codeCoupledName;

  if (rank == 0) {
    codeName="code1";
    codeCoupledName="code2";
  }
  else {
    codeName="code2";
    codeCoupledName="code1";
  }

  char* fileName = NULL;
  if (rank == 0) 
    fileName = "c_surf_coupling_P1P1_polygon_0.txt";
  else
    fileName = "c_surf_coupling_P1P1_polygon_1.txt";

  outputFile = fopen(fileName,"w");

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

  cwipi_create_coupling("c_surf_cpl_P1P1_polygon",                 // Coupling id
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                        codeCoupledName,                           // Coupled application id
                        2,                                         // Geometric entities dimension
                        1,                                         // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        solver_type,                               // Solver type
                        1,                                         // Postprocessing frequency
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
  
  const double xmin = -1;
  const double xmax =  1;
  const double ymin = -1;
  const double ymax =  1;
  int    nx;
  int    ny;

  if (rank == 0) {
    nx = 20;
    ny = 16;
  }
  else {
    nx = 56;
    ny = 60;
  } 

  const int   order = 1;

  creeMaillagePolygone2D(order,
                         localComm,
                         xmin,
                         xmax,
                         ymin,
                         ymax,
                         rank+1,
                         nx,
                         ny,
                         &nVertex,
                         &coords,
                         &nElts,
                         &eltsConnecPointer,
                         &eltsConnec);
  
  fprintf(outputFile, "   Number of vertex : %i\n", nVertex);
  fprintf(outputFile, "   Number of elements : %i\n", nElts);

  cwipi_define_mesh("c_surf_cpl_P1P1_polygon",
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
    printf("        Exchange Proc 0 <-> Proc 1\n");
  
  double *sendValues = NULL;
  double *recvValues = NULL;
  
  sendValues = (double *) malloc(sizeof(double) * nVertex);
  recvValues = (double *) malloc(sizeof(double) * nVertex);

  for (int i = 0; i < nVertex; i++) {
    if (rank == 0)
      sendValues[i] = coords[3*i+1];
    else
      sendValues[i] = coords[3*i+1];
  }

  int nNotLocatedPoints = 0;
  cwipi_exchange_status_t status = cwipi_exchange("c_surf_cpl_P1P1_polygon",
                                                  "echange1",
                                                  1,
                                                  1,     // n_step
                                                  0.1,   // physical_time
                                                  "cooY",
                                                  sendValues,
                                                  "cooY",
                                                  recvValues,
                                                  &nNotLocatedPoints);

  _dumpStatus(outputFile, status);
  _dumpNotLocatedPoints(outputFile, "c_surf_cpl_P1P1_polygon", nNotLocatedPoints);

  /* Coupling deletion
   * ----------------- */

  if (rank == 0)
    printf("        Delete coupling\n");

  cwipi_delete_coupling("c_surf_cpl_P1P1_polygon");

  /* Freeing memory
   * -------------- */

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

  return 0;
}
