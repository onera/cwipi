/*
  This file is part of the CWIPI library. 

  Copyright (C) 2017  ONERA

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
#include <assert.h>
#include <math.h>
#include "grid_mesh.h"
#include <mpi.h>

#include "cwp.h"

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



/*----------------------------------------------------------------------
 *                                                                     
 * Main : linear coupling test                                         
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
  
  FILE* meshFile;
  meshFile = fopen("meshes/mesh_poly_d1", "r");



  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  
  int n_partition = 0;
  const int two = 2;
  while(two * pow(n_partition, two) < comm_world_size) n_partition++;

  int n2 = (int) (two * pow(n_partition, two));

  if (n2 != comm_world_size) {
    if (rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

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


  /* Initialization
   * -------------- */

  int n_code_name = 0;
  char **codeNames = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;

  if (rank == 0 || rank==1 || rank == 2 || rank == 3) {
    n_code_name = 2;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] ="code1";
    codeNames[1] ="code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
  }
   
  char* fileName = NULL;
  if (rank == 0) 
    fileName="c_new_api_0000.txt";
  else if (rank == 1)
    fileName="c_new_api_0001.txt";
  else if (rank == 2)
    fileName="c_new_api_0002.txt";
  else if (rank == 3)
    fileName="c_new_api_0003.txt";

  outputFile = fopen(fileName,"w");

  times_init = malloc(sizeof(double) * n_code_name);

  //CWP_Output_file_set (outputFile);

  for (int i = 0; i < n_code_name; i++) {
    times_init[i] = 0; 
  }
  
  MPI_Comm *localComm = malloc(sizeof(MPI_Comm)*n_code_name);
  CWP_Init(MPI_COMM_WORLD,
           n_code_name,
           (const char **) codeNames,
           is_coupled_rank,
           times_init,
           localComm);


  /* Output redirection
   * ------------------ */

  int currentRank;
  int localCommSize;

  for (int i = 0; i < n_code_name; i++ ) {
    MPI_Comm_rank(localComm[i], &currentRank);
    MPI_Comm_size(localComm[i], &localCommSize);
    printf("Size of localComm[%i]=%i et rang du proc=%i.\n",i,localCommSize,currentRank );
  }

  /* Finalize
   * -------- */
 
  char cpl_id1[] = "cpl_code1_code2";

  printf("Coupling creation\n");
  if (rank == 0 || rank == 1 || rank == 2 || rank == 3) {
    CWP_Cpl_create ("code1", cpl_id1, "code2", CWP_COMM_PAR_WITH_PART,
                    CWP_GEOM_LOCATION, 1,
                    CWP_DISPLACEMENT_STATIC, CWP_FREQ_CPL_TIME_STEP);
  printf("Coupling created\n");
              
    /* Building of the local mesh */
    
    int dimension = 0;             // Dimension of the space
    int nVertex = 0;               // Number of points in the mesh
    int nFace = 0;                 // Number of face
    int nElements = 0;             // Number of cells
    int lFaceConnec = 0;
    int lCellConnec = 0;

    double* coords = NULL;         // Coordinates of the points
    int *faceVertexIdx = NULL;
    int *faceVertex    = NULL;
    int *cellFaceIdx   = NULL;
    int *cellFace      = NULL;

    if  (rank == 0)
      printf("        Read mesh\n");

    
    _read_mesh_dim (meshFile, &dimension, &nVertex, &nFace, &nElements, &lFaceConnec, &lCellConnec);

    coords        = (double *) malloc(dimension * nVertex * sizeof(double));
    faceVertexIdx = (int *) malloc((nFace + 1) * sizeof(int));
    faceVertex    = (int *) malloc(lFaceConnec * sizeof(int));
    cellFaceIdx   = (int *) malloc((nElements + 1) * sizeof(int));
    cellFace      = (int *) malloc(lCellConnec * sizeof(int));

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

    printf("vtx_set\n");
    CWP_Mesh_interf_vtx_set("code1", cpl_id1,0,nVertex,coords,NULL);
 
    printf("3D Cell Polyhedra Block Add\n");
    CWP_Mesh_interf_c_poly_block_add("code1",
                                     cpl_id1,
                                     0,
                                     nElements,
                                     cellFaceIdx,
                                     cellFace,
                                     nFace,
                                     faceVertexIdx,
                                     faceVertex,
                                     NULL);
      

    printf("Finalize Interface Mesh (end_set)\n");
    CWP_Mesh_interf_end_set("code1", cpl_id1);
    
    printf("Interface Mesh deletion\n");
    CWP_Mesh_interf_del("code1", cpl_id1);
        printf("Interface Mesh deleted\n");
  }
   

  fflush(stdout);

  CWP_Finalize();

  MPI_Finalize();

  free (srcName);
  free (localComm);
  free (codeNames);
  free (is_coupled_rank);
  free (times_init);
  fclose (outputFile);
  printf("fin test\n");
  fflush(stdout);
  return 0;
}
