#include <stdlib.h>
#include <mpi.h>

#include "grid_mesh.h"

#define ABS(a)     ((a) <  0  ? -(a) : (a))

/*----------------------------------------------------------------------
 *                                                                     
 * Return a random number in [-1, 1]                                    
 *                                                                     
 * parameters:
 *   coupling_id         <-- Coupling id               
 *   nNotLocatedPoints   <-- Number of not located points
 *---------------------------------------------------------------------*/

static double _random01()
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / ABS(rsigna - rsignb);
  double resultat =   sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}


/*----------------------------------------------------------------------
 *                                                                     
 * Dump status exchange                                                
 *                                                                     
 * parameters:
 *   xmin                <-- grid xmin (global) 
 *   xmax                <-- grid xmax (global)
 *   ymin                <-- grid ymin (global)
 *   ymax                <-- grid ymax (global)
 *   randLevel           <-- random level
 *   nVertexSeg          <-- number of vertices in X and Y
 *   nPartitionSeg       <-- number of partitions in X and Y
 *   nVertex             --> number of vertices of the local partition
 *   coords              --> coordinates of vertices of the local partition
 *   nElts               --> number of elements of the local partition
 *   eltsConnecPointer   --> connectivity index of the local partition 
 *   eltsConnec          --> connectivity of the local partition
 *   localComm           <-- MPI comm of the global grid
 *---------------------------------------------------------------------*/

void  grid_mesh(double xmin, 
                double xmax, 
                double ymin, 
                double ymax, 
                double randLevel,
                int nVertexSeg,
                int nPartitionSeg, 
                int *nVertex,
                double **coords, 
                int *nElts,
                int **eltsConnecPointer,
                int **eltsConnec,
                MPI_Comm localComm)
{

  int currentRank;
  int localCommSize;

  MPI_Comm_rank(localComm, &currentRank);
  MPI_Comm_size(localComm, &localCommSize);

  /* Bandwidth */

  double lX = (xmax - xmin) / nPartitionSeg; 
  double lY = (ymax - ymin) / nPartitionSeg;

  /* Compute local partition bounds with random level */

  int nBoundVerticesSeg = nPartitionSeg + 1;
  int nBoundVertices = nBoundVerticesSeg * nBoundVerticesSeg;
  double *boundRanks = (double *) malloc(3 * sizeof(double) * nBoundVertices);

  if (currentRank == 0) {
    
    for (int j = 0; j < nBoundVerticesSeg; j++) {
      for (int i = 0; i < nBoundVerticesSeg; i++) {
        boundRanks[3 * (j * nBoundVerticesSeg + i)] = xmin + i * lX;
        boundRanks[3 * (j * nBoundVerticesSeg + i) + 1] = ymin + j * lY;
        boundRanks[3 * (j * nBoundVerticesSeg + i) + 2] = 0.;
        if (j != 0 && j != (nBoundVerticesSeg - 1))
          boundRanks[3 * (j * nBoundVerticesSeg + i) + 1] += _random01() * randLevel * lY;  
        if (i != 0 && i != (nBoundVerticesSeg - 1))
          boundRanks[3 * (j * nBoundVerticesSeg + i)] += _random01() * randLevel * lX;  
      }
    }
  }

  MPI_Bcast(boundRanks, 3 * nBoundVertices, MPI_DOUBLE, 0, localComm);

  double *boundRank = (double *) malloc(3 * sizeof(double) * 4);

  int p1, p2, p3, p4;

  int ii = currentRank % nPartitionSeg;
  int jj = currentRank / nPartitionSeg;

  p1 = (nBoundVerticesSeg * jj)       + ii;
  p2 = (nBoundVerticesSeg * jj)       + ii + 1;
  p3 = (nBoundVerticesSeg * (jj + 1)) + ii + 1;
  p4 = (nBoundVerticesSeg * (jj + 1)) + ii;

  boundRank[0 * 3 + 0] = boundRanks[3 * p1    ]; 
  boundRank[0 * 3 + 1] = boundRanks[3 * p1 + 1];
  boundRank[0 * 3 + 2] = boundRanks[3 * p1 + 2];

  boundRank[1 * 3 + 0] = boundRanks[3 * p2    ];
  boundRank[1 * 3 + 1] = boundRanks[3 * p2 + 1];
  boundRank[1 * 3 + 2] = boundRanks[3 * p2 + 2];

  boundRank[2 * 3 + 0] = boundRanks[3 * p3    ];
  boundRank[2 * 3 + 1] = boundRanks[3 * p3 + 1];
  boundRank[2 * 3 + 2] = boundRanks[3 * p3 + 2];

  boundRank[3 * 3 + 0] = boundRanks[3 * p4    ];
  boundRank[3 * 3 + 1] = boundRanks[3 * p4 + 1];
  boundRank[3 * 3 + 2] = boundRanks[3 * p4 + 2];

  free(boundRanks);

  /* Number of vertices and elements in the partition */

  *nVertex = nVertexSeg * nVertexSeg;
  *nElts   = (nVertexSeg - 1) * (nVertexSeg - 1);

  /* Define coordinates */

  *coords = (double *) malloc(sizeof(double) * 3 * (*nVertex) );
  double *const _coords = *coords;

  double deltaU = 2.0/(nVertexSeg - 1);
  double deltaV = 2.0/(nVertexSeg - 1);
  double u = -1;
  double v = -1;
  for (int j = 0; j < nVertexSeg; j++) {
    for (int i = 0; i < nVertexSeg; i++) {
      double randU = u;
      double randV = v;
      if ((i != 0) && (j != 0) && (j != nVertexSeg - 1) && (i != nVertexSeg - 1)) {
        randU +=  _random01() * randLevel * deltaU;
        randV +=  _random01() * randLevel * deltaV;
      }
      
      _coords[3 * (j * nVertexSeg + i) + 0] = 
        0.25 * ((1 - randU - randV + randU * randV) * boundRank[0 * 3 + 0] +
                (1 + randU - randV - randU * randV) * boundRank[1 * 3 + 0] +
                (1 + randU + randV + randU * randV) * boundRank[2 * 3 + 0] +
                (1 - randU + randV - randU * randV) * boundRank[3 * 3 + 0] );
      
      _coords[3 * (j * nVertexSeg + i) + 1] = 
        0.25 * ((1 - randU - randV + randU * randV) * boundRank[0 * 3 + 1] +
                (1 + randU - randV - randU * randV) * boundRank[1 * 3 + 1] +
                (1 + randU + randV + randU * randV) * boundRank[2 * 3 + 1] +
                (1 - randU + randV - randU * randV) * boundRank[3 * 3 + 1] );

      _coords[3 * (j * nVertexSeg + i) + 2] = 0.;

      u += deltaU;
    }
    v += deltaV;
    u = -1;
  }

  free(boundRank);

  /* Define connectivity */

  *eltsConnecPointer = (int *) malloc(sizeof(int) * (*nElts + 1));
  int *const _eltsConnecPointer = *eltsConnecPointer;

  *eltsConnec = (int *) malloc(sizeof(int) * 4 * (*nElts));
  int *const _eltsConnec = *eltsConnec;

  _eltsConnecPointer[0] = 0;
  for (int i = 1; i < *nElts + 1; i++) 
    _eltsConnecPointer[i] = _eltsConnecPointer[i-1] + 4;

  int k = 0;
  for (int j = 0; j < (nVertexSeg - 1); j++) {
    for (int i = 0; i < (nVertexSeg - 1); i++) {
      _eltsConnec[4 * k]     =       j * nVertexSeg + i     + 1;
      _eltsConnec[4 * k + 1] =       j * nVertexSeg + i + 1 + 1;
      _eltsConnec[4 * k + 2] = (j + 1) * nVertexSeg + i + 1 + 1;
      _eltsConnec[4 * k + 3] = (j + 1) * nVertexSeg + i     + 1;
      k++;
    }
  }
}
