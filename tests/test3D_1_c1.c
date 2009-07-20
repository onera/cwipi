/**
 * Simple 3D test of the good working of CWIPI
 *
 */


// **************
// ** Includes **
// **************

// Standards
#include <stdio.h>
//#include <stddef.h>
#include <string.h>

// MPI
#include <mpi.h>

// BFT for monitoring the memory
#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

// The CWIPI
#include "cwipi.h"


// *************************
// ** Auxiliary functions **
// *************************

static void _dumpStatus(cwipi_exchange_status_t status)
{
  switch(status) {
  case CWIPI_EXCHANGE_OK :
    bft_printf("Exchange Ok\n");
    break;
  case CWIPI_EXCHANGE_BAD_RECEIVING :
    bft_printf("Bad receiving\n");
    break;
  default :
    bft_error(__FILE__, __LINE__, 0,"bad exchange status\n");
  }
}

static void _dumpNotLocatedPoints(const char *coupling_id,
                                  const int nNotLocatedPoints)
{
  if ( nNotLocatedPoints > 0) {
    bft_printf("Points non localises :\n");
    const int* notLocatedPoints = cwipi_get_not_located_points(coupling_id);
    for(int i = 0; i < nNotLocatedPoints; i++)
      bft_printf("%i ", notLocatedPoints[i]);
    bft_printf("\n");
  }
}

int read_mesh_dim( FILE* f, int* dimension, int* nVertex, int* nElements, int* nConnecVertex ) {

  int r;
  r = fscanf( f, "%d %d %d %d", dimension, nVertex, nElements, nConnecVertex );
  if( r == EOF ) return 0;
  else return 1;
}

int read_mesh( FILE* f, int dimension, int nVertex, int nElements, double* coords, int* connecPointer, int* connec ) {

  int i, j, r;

  // Reading of the coordinates
  for( i=0; i < nVertex; i++ ) {
    for( j=0; j<dimension; j++ ) {
      r = fscanf( f, "%lf", coords + i*dimension+j);
      if( r == EOF ) return 0;
    }
  }

  // Reading of the connectivity pointer
  for( i=0; i <= nElements; i++ ) {
    r = fscanf( f, "%d", connecPointer + i );
    if( r == EOF ) return 0;
  }

  // Reading of the connectivity
  for( i=0; i<connecPointer[nElements]; i++ ) {
    r = fscanf( f, "%d", &connec[i] );
    if( r == EOF ) return 0;
  }

  return 1;
}


// ***************
// ** Constants **
// ***************


// ******************
// ** Main program **
// ******************

int main( int argc, char* argv[] ) {

  // Declarations
  int rank, nb_procs, localRank;
  char* fileOutput = NULL;
  FILE* outputFile;
  FILE* meshFile;
  MPI_Comm localComm;

  // System initializations
  bft_mem_init( "logmem.txt" );
  MPI_Init( &argc, &argv );

  // Initialisation of the variables
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  //meshFile = fopen( "meshes/test3D_1_c2.mesh", "r" );

  //TODO: corriger le maillage meshes/test3D_1_c2.mesh
  meshFile = fopen( "meshes/test3D_1_c1.mesh", "r" );


  /* Initializations
   * --------------- */

  cwipi_init( MPI_COMM_WORLD,
                  "codeC1",
                  &localComm    );

  MPI_Comm_rank( localComm, &localRank );
  MPI_Comm_size( localComm, &nb_procs );

  if (nb_procs > 1)
    bft_error(__FILE__, __LINE__, 0,"Test3D_1_c1 is not a parallel application\n");


  BFT_MALLOC(fileOutput, strlen("listing_test3D_1_c1_")+ 4 + 1, char );
  sprintf( fileOutput, "listing_test3D_1_c1_%4.4d", localRank );
  outputFile = fopen( fileOutput, "w" );
  BFT_FREE( fileOutput );

  cwipi_set_output_listing( outputFile );

  bft_printf("\nDump after initialization\n");
  bft_printf("---------------------------\n");
  cwipi_dump_application_properties();

  /* -----------------------
   * Test coupling P1 <-> P1
   * ----------------------- */

  {

    /* Initialization of the coupling */
    bft_printf("Test 1 : Test coupling P1 <-> P1\n");
    bft_printf("\n");

    if  (localRank == 0)
      printf("Test 1 : Test couplage P1 <-> P1\n");

    if  (localRank == 0)
      printf("         Create coupling\n");

    cwipi_create_coupling("test3D_1",                   // Name of the coupling
                              CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                              "codeC2",              // Coupled code
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

    if  (localRank == 0)
      printf("         Read mesh\n");

    read_mesh_dim( meshFile, &dimension, &nVertex, &nElements, &nConnecVertex );
    BFT_MALLOC(coords, dimension*nVertex, double);
    BFT_MALLOC(eltsConnecPointer, nElements+1, int);
    BFT_MALLOC(eltsConnec, nConnecVertex, int);
    read_mesh( meshFile, dimension, nVertex, nElements, coords, eltsConnecPointer, eltsConnec );
    fclose(meshFile);


    if  (localRank == 0)
      printf("         Define mesh\n");

    cwipi_define_mesh("test3D_1",
                          nVertex,
                          nElements,
                          coords,
                          eltsConnecPointer,
                          eltsConnec);

    /* Sending of the coordinate X
       Receiving of the coordinate Y*/

    BFT_MALLOC(values, nVertex, double);
    for( int i = 0; i < nVertex; i++ ) values[i] = coords[3*i]*coords[3*i];
    BFT_MALLOC(localValues, nVertex, double);

    if  (localRank == 0)
      printf("         Exchange\n");

   cwipi_exchange_status_t status = cwipi_exchange("test3D_1",
                                                            "echange1",
                                                            1,
                                                            1,     // n_step
                                                            0.1,   // physical_time
                                                            "cooX",
                                                            values,
                                                            "cooY",
                                                            localValues,
                                                            &nNotLocatedPoints);
    _dumpStatus(status);
    _dumpNotLocatedPoints("test3D_1", nNotLocatedPoints);

    /* Deletion of the coupling object */

    if  (localRank == 0)
      printf("         Delete coupling\n");

    cwipi_delete_coupling("test3D_1");

    /* Rest in peace */

    if (coords != NULL)            
      BFT_FREE(coords);

    if (eltsConnec != NULL)        
      BFT_FREE(eltsConnec);

    if (eltsConnecPointer != NULL)
      BFT_FREE(eltsConnecPointer);

    if (values != NULL)           
      BFT_FREE(values);

    if (localValues != NULL)      
      BFT_FREE(localValues);

    bft_printf("--------------------------------------------------------\n");
  }

  /* End of the MPI communications */
  /* ----------------------------- */

  cwipi_finalize();

  bft_mem_end();
  fclose(outputFile);
  return 0;
}
