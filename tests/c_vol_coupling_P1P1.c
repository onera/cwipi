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

/**
 * Simple 3D test
 *
 */

// **************
// ** Includes **
// **************

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "cwipi.h"
#include "cwipi_config.h"

// *************************
// ** Static functions    **
// *************************

/*----------------------------------------------------------------------
 *
 * Dump status exchange
 *
 * parameters:
 *   status              <-- Exchange status
 *---------------------------------------------------------------------*/

static void _dumpStatus(FILE *outputFile, cwipi_exchange_status_t status)
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

static void _dumpNotLocatedPoints(FILE *outputFile,
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
 * Read mesh dimension
 *
 * parameters:
 *   f                   <-- Mesh file
 *   dimension           --> Dimension
 *   nvertex             --> number of vertices
 *   nElements           --> number of elements
 *   nConnecVertex       --> size of connectivity
 *---------------------------------------------------------------------*/

static int _read_mesh_dim(FILE *f,
                          int *dimension,
                          int *nVertex,
                          int *nElements,
                          int *nConnecVertex)
{
  int r;
  r = fscanf( f, "%d %d %d %d", dimension, nVertex, nElements, nConnecVertex );
  if (r == EOF)
    return 0;
  else return 1;
}

/*----------------------------------------------------------------------
 *
 * Read mesh dimension
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
                      int dimension,
                      int nVertex,
                      int nElements,
                      double *coords,
                      int *connecPointer,
                      int *connec)
{
  int i, j, r;

  // Read coordinates
  for (i = 0; i < nVertex; i++) {
    for (j = 0; j < dimension; j++) {
      r = fscanf(f, "%lf", coords + i * dimension + j);
      if (r == EOF)
        return EXIT_FAILURE;
    }
  }

  // Read connectivity index
  for (i = 0; i <= nElements; i++ ) {
    r = fscanf(f, "%d", connecPointer + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read connectivity
  for (i = 0; i < connecPointer[nElements]; i++ ) {
    r = fscanf(f, "%d", &connec[i]);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

// ******************
// ** Main program **
// ******************

int main( int argc, char* argv[] ) {

  // Declarations
  int rank, nb_procs, localRank;
  int comm_world_size;
  char* fileOutput = NULL;
  FILE* outputFile;
  FILE* meshFile;
  MPI_Comm localComm;

  // System initializations
  MPI_Init(&argc, &argv);

  // Initialisation of the variables
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

  if (comm_world_size != 2) {
    if (rank == 0)
      printf("      Not executed : only available for 2 processus\n");
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  if (rank == 0)
    meshFile = fopen(CWP_MESH_DIR"test3D_1_c1.mesh", "r");
  else
    meshFile = fopen(CWP_MESH_DIR"test3D_1_c2.mesh", "r");

  fileOutput = (char *) malloc((strlen("c_vol_cpl_P1P1_c1_") + 4 + 1 + 4) * sizeof(char));
  sprintf(fileOutput, "c_vol_cpl_P1P1_c1_%4.4d.txt", rank);
  outputFile = fopen(fileOutput, "w");
  free(fileOutput);

  cwipi_set_output_listing( outputFile );

  /* Initializations
   * --------------- */

  if (rank == 0)
    cwipi_init(MPI_COMM_WORLD,
               "codeC1",
               &localComm);
  else
    cwipi_init(MPI_COMM_WORLD,
               "codeC2",
               &localComm);

  MPI_Comm_rank(localComm, &localRank);
  MPI_Comm_size(localComm, &nb_procs);

  if (nb_procs > 1) {
    printf("Test3D_1_c1 is not a parallel application\n");
    exit(EXIT_FAILURE);
  }


  fprintf(outputFile, "\nDump after initialization\n");
  fprintf(outputFile, "---------------------------\n");
  cwipi_dump_application_properties();

  /* -----------------------
   * Test coupling P1 <-> P1
   * ----------------------- */

  {

    /* Initialization of the coupling */
    fprintf(outputFile, "Test 1 : Test coupling P1 <-> P1\n");
    fprintf(outputFile, "\n");

    if (rank == 0)
      printf("        Create coupling\n");

    if (rank == 0)
      cwipi_create_coupling("c_vol_cpl_P1P1",                   // Name of the coupling
                            CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                            "codeC2",              // Coupled code
                            3,                            // Dimension of the geometry
                            0.1,                          // Geometrical epsilon
                            CWIPI_STATIC_MESH,        // Static mesh
                            CWIPI_SOLVER_CELL_VERTEX, // Type of the fields
                            1,                            // Post-processing frequency
                            "EnSight Gold",               // Post-processing format
                            "text");                      // Post-processing options
    else
      cwipi_create_coupling("c_vol_cpl_P1P1",                   // Name of the coupling
                            CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                            "codeC1",              // Coupled code
                            3,                            // Dimension of the geometry
                            0.1,                          // Geometrical epsilon
                            CWIPI_STATIC_MESH,        // Static mesh
                            CWIPI_SOLVER_CELL_VERTEX, // Type of the fields
                            1,                            // Post-processing frequency
                            "EnSight Gold",               // Post-processing format
                            "text");                      // Post-processing options

    /* Building of the local mesh */

    int dimension = 0;             // Dimension of the space
    int nVertex = 0;               // Number of points in the mesh
    double* coords = NULL;         // Coordinates of the points
    int nElements = 0;             // Number of cells
    int nConnecVertex = 0;         // Number of cell vertices
    int* eltsConnecPointer = NULL; // Connectivity pointer
    int* eltsConnec = NULL;        // Connectivity of the cells
    double* values = NULL;         // Received field
    double* localValues = NULL;    // Sent field
    int nNotLocatedPoints;         // Number of points out of the mesh

    if  (rank == 0)
      printf("        Read mesh\n");

    assert (meshFile != NULL);
    _read_mesh_dim( meshFile, &dimension, &nVertex, &nElements, &nConnecVertex );

    coords = (double *) malloc(dimension * nVertex * sizeof(double));
    eltsConnecPointer = (int *) malloc((nElements + 1) * sizeof(int));
    eltsConnec = (int *) malloc(nConnecVertex * sizeof(int));

    _read_mesh( meshFile, dimension, nVertex, nElements, coords, eltsConnecPointer, eltsConnec );
    fclose(meshFile);

    if  (rank == 0)
      printf("        Define mesh\n");

    cwipi_define_mesh("c_vol_cpl_P1P1",
                      nVertex,
                      nElements,
                      coords,
                      eltsConnecPointer,
                      eltsConnec);

    /* Sending of the coordinate X
       Receiving of the coordinate Y*/

    values = (double *) malloc(nVertex * sizeof(double));

    for (int i = 0; i < nVertex; i++)
      values[i] = coords[3*i];

    localValues = (double *) malloc(nVertex * sizeof(double));

    cwipi_locate("c_vol_cpl_P1P1");

    if (rank == 0)
      printf("        Exchange\n");


    cwipi_exchange_status_t status = cwipi_exchange("c_vol_cpl_P1P1",
                                                    "echange1",
                                                    1,
                                                    1,     // n_step
                                                    0.1,   // physical_time
                                                    "cooX",
                                                    values,
                                                    "cooY",
                                                    localValues,
                                                    &nNotLocatedPoints);
    _dumpStatus(outputFile, status);
    _dumpNotLocatedPoints(outputFile, "c_vol_cpl_P1P1", nNotLocatedPoints);

    /* Deletion of the coupling object */

    if (rank == 0)
      printf("        Delete coupling\n");

    cwipi_delete_coupling("c_vol_cpl_P1P1");

    /* Check barycentric coordinates */

    if (rank == 0)
      printf("        Check results\n");

    double err = fabs(localValues[0] - values[0]);

    for (int i = 1; i < nVertex; i++) {
      err = ((fabs(localValues[i] - values[i])) < (err) ? (err) :
             (fabs(localValues[i] - values[i])));

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

    free(coords);
    free(eltsConnec);
    free(eltsConnecPointer);
    free(values);
    free(localValues);

    fprintf(outputFile, "--------------------------------------------------------\n");
  }

  /* End of the MPI communications */
  /* ----------------------------- */

  cwipi_finalize();

  fclose(outputFile);

  free (srcName);

  MPI_Finalize();

  return EXIT_SUCCESS;
}