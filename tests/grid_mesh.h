#ifndef __GRID_MESH_H__
#define __GRID_MESH_H__

#include <mpi.h>

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
                MPI_Comm localComm);
#endif
