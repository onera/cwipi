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

#include "cwipi.h"

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
                          int *nFace, 
                          int *nElt,
                          int *lFaceConnec,
                          int *lCellConnec)
 
{
  int r;
  r = fscanf(f, "%d %d %d %d %d %d", 
             dimension, 
             nVertex, 
             nFace, 
             nElt,
             lFaceConnec,
             lCellConnec);
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
                      int nFace,
                      int nElt,
                      int lFaceConnec,
                      int lCellConnec,
                      double *coords, 
                      int *faceVertexIdx, 
                      int *faceVertex, 
                      int *cellFaceIdx, 
                      int *cellFace) 
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

  // Read face -> vertex connectivity index
  for (i = 0; i < nFace + 1; i++ ) {
    r = fscanf(f, "%d", faceVertexIdx + i);
    if (r == EOF) 
      return EXIT_FAILURE;
  }

  // Read face -> vertex connectivity
  for (i = 0; i < lFaceConnec; i++ ) {
    r = fscanf(f, "%d", faceVertex + i);
    if (r == EOF) 
      return EXIT_FAILURE;
  }

  // Read cell -> face connectivity index
  for (i = 0; i < nElt + 1; i++ ) {
    r = fscanf(f, "%d", cellFaceIdx + i);
    if (r == EOF) 
      return EXIT_FAILURE;
  }

  // Read cell -> face connectivity
  for (i = 0; i < lCellConnec; i++ ) {
    r = fscanf(f, "%d", cellFace + i);
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
    meshFile = fopen("meshes/mesh_poly_d1", "r");
  else
    meshFile = fopen("meshes/mesh_poly_d2", "r");

  fileOutput = (char *) malloc((strlen("c_vol_poly_cpl_P1P1_") + 4 + 1 + 4) * sizeof(char));
  sprintf(fileOutput, "c_vol_poly_cpl_P1P1_%4.4d.txt", rank);
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
      cwipi_create_coupling("c_vol_cpl_poly_P1P1",                   // Name of the coupling
                            CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                            "codeC2",              // Coupled code
                            3,                            // Dimension of the geometry
                            0.7,                          // Geometrical epsilon
                            CWIPI_STATIC_MESH,        // Static mesh
                            CWIPI_SOLVER_CELL_VERTEX, // Type of the fields
                            1,                            // Post-processing frequency
                            "EnSight Gold",               // Post-processing format
                            "text");                      // Post-processing options
    else
      cwipi_create_coupling("c_vol_cpl_poly_P1P1",                   // Name of the coupling
                            CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                            "codeC1",              // Coupled code
                            3,                            // Dimension of the geometry
                            0.7,                          // Geometrical epsilon
                            CWIPI_STATIC_MESH,        // Static mesh
                            CWIPI_SOLVER_CELL_VERTEX, // Type of the fields
                            1,                            // Post-processing frequency
                            "EnSight Gold",               // Post-processing format
                            "text");                      // Post-processing options
    
    /* Building of the local mesh */
    
    int dimension = 0;             // Dimension of the space
    int nVertex = 0;               // Number of points in the mesh
    int nFace = 0;                 // Number of face
    int nElements = 0;             // Number of cells
    int lFaceConnec = 0;
    int lCellConnec = 0;

    double* coords = NULL;         // Coordinates of the points
    int nConnecVertex = 0;         // Number of cell vertices
    double* values = NULL;         // Received field
    double* localValues = NULL;    // Sent field
    int nNotLocatedPoints;         // Number of points out of the mesh
    int *faceVertexIdx = NULL;
    int *faceVertex    = NULL;
    int *cellFaceIdx   = NULL;
    int *cellFace      = NULL;

    if  (rank == 0)
      printf("        Read mesh\n");

    
    int unique = 0;

    if (!unique) {

      _read_mesh_dim (meshFile, &dimension, &nVertex, &nFace, &nElements, &lFaceConnec, &lCellConnec);
      printf("dims : %i %i %i %i %i %i\n",dimension, nVertex, nFace, nElements, lFaceConnec, lCellConnec);
    
    }

    else {
      nVertex     = 16;
      nFace       = 10;
      nElements   = 1;
      lFaceConnec = 48;
      lCellConnec = nFace;
      dimension   = 3;
    }

    coords        = (double *) malloc(dimension * nVertex * sizeof(double));
    faceVertexIdx = (int *) malloc((nFace + 1) * sizeof(int));
    faceVertex    = (int *) malloc(lFaceConnec * sizeof(int));
    cellFaceIdx   = (int *) malloc((nElements + 1) * sizeof(int));
    cellFace      = (int *) malloc(lCellConnec * sizeof(int));

    if (unique) {

      int imax = 244;
      int *indirec = (int *) malloc((imax  ) * sizeof(int));
      
      for (int k = 0; k < imax ; k++) {
        indirec[k] = 0;
      }
 
      indirec[62] = 7;
      indirec[63] = 3;
      indirec[64] = 4;
      indirec[71] = 5;
      indirec[75] = 8;
      indirec[136] = 2;
      indirec[137] = 9;
      indirec[138] = 10;
      indirec[141] = 16;
      indirec[142] = 14;
      indirec[181] = 6;
      indirec[182] = 15;
      indirec[222] = 1;
      indirec[224] = 13;
      indirec[243] = 11;
      indirec[242] = 12;
    
      int k1 = 0;
      cellFaceIdx[0] = 0;
      cellFaceIdx[1] = nFace;
      
      cellFace[k1++] =  (k1+1);
      cellFace[k1++] =  (k1+1);
      cellFace[k1++] = -(k1+1);
      cellFace[k1++] = -(k1+1);
      cellFace[k1++] = -(k1+1);
      cellFace[k1++] = -(k1+1);
      cellFace[k1++] = -(k1+1);
      cellFace[k1++] = -(k1+1);
      cellFace[k1++] = -(k1+1);
      cellFace[k1++] = -(k1+1);

      printf("toto %d %d\n", cellFace[0], cellFace[8]);
      
      faceVertexIdx[0] = 0;
      
      k1 = 0;
      faceVertex[k1++] = indirec[222];
      faceVertex[k1++] = indirec[137];
      faceVertex[k1++] = indirec[138];
      faceVertex[k1++] = indirec[243];
      faceVertex[k1++] = indirec[242];
      faceVertex[k1++] = indirec[224];
      faceVertexIdx[1] = 6;
      
      faceVertex[k1++] = indirec[224];
      faceVertex[k1++] = indirec[182];
      faceVertex[k1++] = indirec[142];
      faceVertex[k1++] = indirec[141];
      faceVertex[k1++] = indirec[222];
      faceVertexIdx[2] = 11; 
      
      faceVertex[k1++] = indirec[63];
      faceVertex[k1++] = indirec[62];
      faceVertex[k1++] = indirec[75];
      faceVertex[k1++] = indirec[71];
      faceVertex[k1++] = indirec[64];
      faceVertexIdx[3] = 16;
      
      faceVertex[k1++] = indirec[62];
      faceVertex[k1++] = indirec[63];
      faceVertex[k1++] = indirec[136];
      faceVertex[k1++] = indirec[137];
      faceVertex[k1++] = indirec[138];
      faceVertexIdx[4] = 21;
      
      faceVertex[k1++] = indirec[142];
      faceVertex[k1++] = indirec[141];
      faceVertex[k1++] = indirec[136];
      faceVertex[k1++] = indirec[63];
      faceVertex[k1++] = indirec[64];
      faceVertexIdx[5] = 26;
      
      faceVertex[k1++] = indirec[182];
      faceVertex[k1++] = indirec[142];
      faceVertex[k1++] = indirec[64];
      faceVertex[k1++] = indirec[71];
      faceVertex[k1++] = indirec[181];
      faceVertexIdx[6] = 31;
      
      faceVertex[k1++] = indirec[141];
      faceVertex[k1++] = indirec[222];
      faceVertex[k1++] = indirec[137];
      faceVertex[k1++] = indirec[136];
      faceVertexIdx[7] = 35;
      
      faceVertex[k1++] = indirec[181];
      faceVertex[k1++] = indirec[71];
      faceVertex[k1++] = indirec[75];
      faceVertex[k1++] = indirec[243];
      faceVertex[k1++] = indirec[242];
      faceVertexIdx[8] = 40;
      
      faceVertex[k1++] = indirec[138];
      faceVertex[k1++] = indirec[243];
      faceVertex[k1++] = indirec[75];
      faceVertex[k1++] = indirec[62];
      faceVertexIdx[9] = 44;
      
      faceVertex[k1++] = indirec[242];
      faceVertex[k1++] = indirec[224];
      faceVertex[k1++] = indirec[182];
      faceVertex[k1++] = indirec[181];
      faceVertexIdx[10] = 48;
      
      k1 = 0;
      coords[k1++] =1.00000e+00;
      coords[k1++] =    1.00000e+00;
      coords[k1++] =    5.00000e-01;

      coords[k1++] =8.62620e-01;
      coords[k1++] =    8.64922e-01;
      coords[k1++] =    5.00920e-01;

      coords[k1++] =6.80400e-01;
      coords[k1++] =    7.01313e-01;
      coords[k1++] =    2.42025e-01;

      coords[k1++] =5.43250e-01;
      coords[k1++] =    8.41293e-01;
      coords[k1++] =   -3.95872e-02;

      coords[k1++] =6.75290e-01;
      coords[k1++] =    6.95194e-01;
      coords[k1++] =   -3.15702e-01;

      coords[k1++] =8.60060e-01;
      coords[k1++] =    8.62214e-01;
      coords[k1++] =   -5.05754e-01;

      coords[k1++] =7.93540e-01;
      coords[k1++] =    4.32371e-01;
      coords[k1++] =    9.53616e-02;

      coords[k1++] =8.06580e-01;
      coords[k1++] =    4.26460e-01;
      coords[k1++] =   -1.71143e-01;

      coords[k1++] =1.00000e+00;
      coords[k1++] =    7.92893e-01;
      coords[k1++] =    5.00000e-01;

      coords[k1++] =1.00000e+00;
      coords[k1++] =    5.00000e-01;
      coords[k1++] =    2.07107e-01;

      coords[k1++] =1.00000e+00;
      coords[k1++] =    5.00000e-01;
      coords[k1++] =   -2.07107e-01;

      coords[k1++] =1.00000e+00;
      coords[k1++] =    7.92893e-01;
      coords[k1++] =   -5.00000e-01;

      coords[k1++] =1.00000e+00;
      coords[k1++] =    1.00000e+00;
      coords[k1++] =   -5.00000e-01;

      coords[k1++] =7.07110e-01;
      coords[k1++] =    1.00000e+00;
      coords[k1++] =    0.00000e+00;

      coords[k1++] =7.92890e-01;
      coords[k1++] =    1.00000e+00;
      coords[k1++] =   -5.00000e-01;

      coords[k1++] =7.92890e-01;
      coords[k1++] =    1.00000e+00;
      coords[k1++] =    5.00000e-01;


    }

    else {
      _read_mesh (meshFile,
                  dimension,
                  nVertex,
                  nFace,
                  nElements,
                  lFaceConnec,
                  lCellConnec,
                  coords,
                  faceVertexIdx,
                  faceVertex,
                  cellFaceIdx,
                  cellFace);
    
      fclose(meshFile);
    }

    if (!unique) {

      for (int i = 0; i < 3*nVertex; i++) {
        coords[i] = coords[i] - 0.1;
      }
      
      for (int i = 0; i < 3*nVertex; i++) {
        coords[i] = 10 * coords[i];
      }
   
      const double dila = 1.1;
      
      if (rank == 0) {
        for (int i = 0; i < 3*nVertex; i++) {
        coords[i] = dila * coords[i];
        }
      }

    }

    if  (rank == 0)
      printf("        Define mesh\n");

    int t[2];
    t[0] = 0;
    t[1] = 0;
    cwipi_define_mesh("c_vol_cpl_poly_P1P1",
                      nVertex,
                      0,
                      coords,
                      t,
                      t);

    cwipi_add_polyhedra("c_vol_cpl_poly_P1P1",
                        nElements,
                        cellFaceIdx,
                        cellFace,
                        nFace,
                        faceVertexIdx,
                        faceVertex);

    int nPts = 1;
    double *coordsPts        = (double *) malloc(3 * nPts * sizeof(double));
    coordsPts[3*0  ]   = 1.01945e+00;
    coordsPts[3*0+1  ] = 1.10000e+00;
    coordsPts[3*0+2  ] = 2.75000e-01;

    if (unique) {
      cwipi_set_points_to_locate ("c_vol_cpl_poly_P1P1", nPts, coordsPts);
    }
    else {
      nPts = nVertex;
    }

    /* Sending of the coordinate X
       Receiving of the coordinate Y*/

    values = (double *) malloc(nVertex * sizeof(double));
   
    for (int i = 0; i < nVertex; i++) 
      values[i] = coords[3*i];
    
    localValues = (double *) malloc(nPts * sizeof(double));

    if (rank == 0)
      printf("        Exchange\n");

    cwipi_exchange_status_t status;

    status = cwipi_exchange("c_vol_cpl_poly_P1P1",
                            "echange1",
                            1,
                            1,     // n_step
                            0.1,   // physical_time
                            "cooX",
                            values,
                            "cooX",
                            localValues,
                            &nNotLocatedPoints);

    
    for (int i = 0; i < nVertex; i++) 
      localValues[i] = fabs(localValues[i] - values[i]);

    status = cwipi_exchange("c_vol_cpl_poly_P1P1",
                            "echange1",
                            1,
                            1,     // n_step
                            0.1,   // physical_time
                            "cooX2",
                            localValues,
                            "cooX2",
                            values,
                            &nNotLocatedPoints);



    _dumpStatus(outputFile, status);
    _dumpNotLocatedPoints(outputFile, "c_vol_cpl_poly_P1P1", nNotLocatedPoints);

    /* Deletion of the coupling object */

    if (rank == 0)
      printf("        Delete coupling\n");

    cwipi_delete_coupling("c_vol_cpl_poly_P1P1");

    /* Check barycentric coordinates */
    if (rank == 0)
      printf("        Check results\n");    

    if (!unique) {
      double err = fabs(localValues[0] - values[0]);
      
      for (int i = 1; i < nVertex; i++) {
        err = ((fabs(localValues[i] - values[i])) < (err) ? (err) :
               (fabs(localValues[i] - values[i])));
        
      }

      if (err >= 1e-6) {
        if (rank == 0) {
          printf("        !!! Error = %12.5e\n", err);
          return EXIT_FAILURE;
        }
      }
    }

    free(coords);
    free(faceVertexIdx);
    free(faceVertex);
    free(cellFaceIdx);
    free(cellFace);
    free(values);
    free(localValues);

    fprintf(outputFile, "--------------------------------------------------------\n");
  }

  /* End of the MPI communications */
  /* ----------------------------- */

  cwipi_finalize();

  fclose(outputFile);

  MPI_Finalize();

  return EXIT_SUCCESS;
}
