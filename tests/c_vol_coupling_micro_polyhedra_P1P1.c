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
#include <float.h>

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
 *   nb_Vertex             --> number of vertices
 *   nElements           --> number of elements
 *   nConnecVertex       --> size of connectivity
 *---------------------------------------------------------------------*/

static int _read_mesh_dim(FILE *f, 
                          int *dimension, 
                          int *nb_Vertex, 
                          int *nFace, 
                          int *nElt,
                          int *lFaceConnec,
                          int *lCellConnec)
 
{
  int r;
  r = fscanf(f, "%d %d %d %d %d %d", 
             dimension, 
             nb_Vertex, 
             nFace, 
             nElt,
             lFaceConnec,
             lCellConnec);
  if (r == EOF)
    return 0;
  else return 1;
}




void _goto(FILE *f,char* word) 
{ char test[1000];
  int r;
  while(strcmp(test,word)!=0) {
      r = fscanf(f, "%s",test);
      if (r == EOF) 
        return EXIT_FAILURE;
      
      printf("%s\n",&test[0]);
  }
}


typedef struct elType elType;
/* Declare the struct with integer members x, y */
struct elType {
   int    nNodes;
   char*    descri;
};



/*----------------------------------------------------------------------
 *                                                                     
 * Read mesh dimension                                             
 *                                                                     
 * parameters:
 *   f                   <-- Mesh file                 
 *   dimension           --> Dimension                   
 *   nb_Vertex             <-- number of vertices
 *   nElements           <-- number of elements
 *   nConnecVertex       <-- size of connectivity
 *   coords              --> vertices coordinates
 *   connecPointer       --> connectivity index  
 *   connec              --> connectivity
 *---------------------------------------------------------------------*/


static int _read_mesh(FILE *f, 
                      int* nb_Vertex, 
                      int* nb_Elts,
                      int* nBlock,
                      int** nElBlock,
                      int** typeBlock,
                      int*** connec,
                      double **coords) {

  elType* elementType = (elType*)malloc(sizeof(elType*)*100);
  


elementType[1] = (elType){2,"line"};
elementType[2] = (elType){3,"triangle"};
elementType[3] = (elType){4,"quadrangle"};
elementType[4] = (elType){4,"tetrahedron"};
elementType[5] = (elType){8,"hexahedron"};
elementType[6] = (elType){6,"prism"};
elementType[7] = (elType){5,"pyramid"};
elementType[8] = (elType){3,"second order line (2 nodes associated with the vertices and 1 with the edge)"};
elementType[9] = (elType){6,"second order triangle (3 nodes associated with the vertices and 3 with the edges)"};
elementType[10] = (elType){9,"second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)"};
elementType[11] = (elType){10,"second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)"};
elementType[12] = (elType){27,"second order hexahedron (8 nodes associated with the vertices, 12 with the edges 6 with the faces and 1 with the volume)"};
elementType[13] = (elType){18,"second order prism (6 nodes associated with the vertices], 9 with the edges and 3 with the quadrangular faces)"};
elementType[14] = (elType){14,"second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)"};
elementType[15] = (elType){1,"point"};
elementType[16] = (elType){8,"second order quadrangle (4 nodes associated with the vertices and 4 with the edges)"};
elementType[17] = (elType){20,"second order hexahedron (8 nodes associated with the vertices and 12 with the edges)"};
elementType[18] = (elType){15,"second order prism (6 nodes associated with the vertices and 9 with the edges)"};
elementType[19] = (elType){13,"second order pyramid (5 nodes associated with the vertices and 8 with the edges)"};
elementType[20] = (elType){9,"third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)"};
elementType[21] = (elType){10,"third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)"};
elementType[22] = (elType){12,"fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)"};
elementType[23] = (elType){15,"fourth order triangle (3 nodes associated with the vertices, 9 with the edges 3 with the face)"};
elementType[24] = (elType){15,"fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)"};
elementType[25] = (elType){21,"fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges 6 with the face)"};
elementType[26] = (elType){4,"third order edge (2 nodes associated with the vertices 2 internal to the edge)"};
elementType[27] = (elType){5,"fourth order edge (2 nodes associated with the vertices 3 internal to the edge)"};
elementType[28] = (elType){6,"fifth order edge (2 nodes associated with the vertices 4 internal to the edge)"};
elementType[29] = (elType){20,"third order tetrahedron (4 nodes associated with the vertices 12 with the edges 4 with the faces)"};
elementType[30] = (elType){35,"fourth order tetrahedron (4 nodes associated with the vertices 18 with the edges 12 with the faces 1 in the volume)"};
elementType[31] = (elType){56,"fifth order tetrahedron (4 nodes associated with the vertices 24 with the edges 24 with the faces 4 in the volume)"};
elementType[92] = (elType){64,"third order hexahedron (8 nodes associated with the vertices 24 with the edges 24 with the faces 8 in the volume)"};
elementType[93] = (elType){125,"fourth order hexahedron (8 nodes associated with the vertices 36 with the edges 54 with the faces 27 in the volume)"};

  int i, j, r;


  char test[1000];
  
  _goto(f,"$Nodes");
   // _goto(f,"$EndNodes");
  int nv,nEl;
  r = fscanf(f, "%i",nBlock);
  r = fscanf(f, "%i",nb_Vertex);
  printf("%i\n",*nb_Vertex);

  int dimension =3;

  *coords = (double *) malloc(dimension * (*nb_Vertex) * sizeof(double));

  for(int block = 1; block<(*nBlock)+1;block++) {
    r = fscanf(f, "%i",&nv);  
    printf("nv %i\n",nv);
    r = fscanf(f, "%i",&nv);  
    printf("nv %i\n",nv);
    r = fscanf(f, "%i",&nv);
    r = fscanf(f, "%i",&nv);
    printf("nv %i\n",nv);  
 //   r = fscanf(f, "%i",&nv);  

    for (int i = 0; i < nv; i++) {
      r = fscanf(f, "%i",&nEl);  
      nEl=nEl-1;
      r = fscanf(f, "%lf",(*coords)+3*nEl);
      r = fscanf(f, "%lf",(*coords)+3*nEl+1);
      r = fscanf(f, "%lf",(*coords)+3*nEl+2);       
     // printf("nEl %i block %i nv %i x %f y %f z %f\n",nEl,block,nv,(*coords)[3*nEl],(*coords)[3*nEl+1],(*coords)[3*nEl+2]);
      if (r == EOF) 
        return EXIT_FAILURE;
    }
  }
  
  int IelType,nEl2,poub;

  
  _goto(f,"$Elements");
  r = fscanf(f, "%i",&poub);
  r = fscanf(f, "%i",nb_Elts);
  printf("%i\n",*nb_Elts);
  *connec = (int**)malloc((*nBlock)*sizeof(int*));
  *nElBlock = (int*) malloc(sizeof(int)*(*nBlock));
  *typeBlock = (int*) malloc(sizeof(int)*(*nBlock));

  for(int block = 0; block<(*nBlock);block++) {
    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);
    r = fscanf(f, "%i",&IelType);  
    r = fscanf(f, "%i",&nv);
    (*nElBlock)[block] = nv;

    (*typeBlock)[block] = IelType;

    int size_el;

    size_el = elementType[IelType].nNodes;

    (*connec)[block]  = (int *) malloc(nv * size_el * sizeof(int));
    
    for (int i = 0; i < nv; i++) {
      r = fscanf(f, "%i",&nEl2); 
      for(int jv=0;jv<size_el;jv++) { 
        r = fscanf(f, "%i",(*connec)[block]+size_el*i+jv);
       // printf("nEl2 %i block %i nv %i connec %i size_el %i\n",nEl2,block,nv,(*connec)[block][size_el*i+jv],size_el);
      }
      if (r == EOF) 
        return EXIT_FAILURE;
    }
   // free(connec);
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

  meshFile = fopen("carre01.msh", "r");


  printf("        Read mesh\n");


  int nb_Vertex,nb_Elts,nBlock;
  int* nElBlock=NULL;
  int* typeBlock= NULL;
  int** connec = NULL;
  double* coords = NULL;
  
  _read_mesh(meshFile,
             &nb_Vertex, 
             &nb_Elts,
             &nBlock,
             &nElBlock,
             &typeBlock,
             &connec,
             &coords);

    printf("coords %i\n",coords[1]);
    fclose(meshFile);

  
  return EXIT_SUCCESS;
}
