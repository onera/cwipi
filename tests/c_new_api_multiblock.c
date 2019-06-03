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







void _goto(FILE *f,char* word) 
{ char test[1000];
  int r;
  while(strcmp(test,word)!=0) {
      r = fscanf(f, "%s",test);
      if (r == EOF) 
        return EXIT_FAILURE;
      
 //     printf("%s\n",&test[0]);
  }
}

void _generate_gmsh_mesh(char* geofile,int localCommSize,int order){
  char s[1000]; //= (char*) malloc(1000*sizeof(char));
  int len=0;
  len = strlen(geofile);
 // for(len = 0; geofile[len] != '\0'; ++len);
  
  char filename[len-4];//(char*)malloc(sizeof(char)*(len-4));
  for (int i=0;i<len-4;i++) filename[i]=geofile[i];


  printf("UUU %i\n",len);
  if(localCommSize>1) 
    sprintf(s, "gmsh %s -2 -format msh -order %i  -part %i -part_split -o %s.msh",geofile,order,localCommSize,filename);
  else
    sprintf(s, "gmsh %s -2 -format msh -order %i  -o %s.msh",geofile,order,filename);    

  printf("sssss %s\n",s);  
  system(s);
  //free(filename);
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
 *   nVtx             <-- number of vertices
 *   nElements           <-- number of elements
 *   nConnecVertex       <-- size of connectivity
 *   coords              --> vertices coordinates
 *   connecPointer       --> connectivity index  
 *   connec              --> connectivity
 *---------------------------------------------------------------------*/
 static int _read_mesh_dim(FILE *f, 
                      int* nVtx, 
                      int* nb_Elts,
                      int* nBlock,
                      int** nElBlock,
                      int** typeBlock) {

  elType elementType[100];// = (elType*)malloc(sizeof(elType*)*100);
  
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
  double poubd;
  r = fscanf(f, "%i",nBlock);
  r = fscanf(f, "%i",nVtx);
  
  printf("nVtx %i\n",*nVtx);

  int dimension =3;

  for(int block = 1; block<(*nBlock)+1;block++) {
    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);
    r = fscanf(f, "%i",&nv);

    for (int i = 0; i < nv; i++) {
      r = fscanf(f, "%i",&nEl);  
      nEl=nEl-1;
      r = fscanf(f, "%lf",&poubd);
      r = fscanf(f, "%lf",&poubd);
      r = fscanf(f, "%lf",&poubd);       

      if (r == EOF) 
        return EXIT_FAILURE;
    }
  }

  int IelType,nEl2,poub;

  _goto(f,"$Elements");
  r = fscanf(f, "%i",nBlock);
  r = fscanf(f, "%i",nb_Elts);
  printf("nb_Elts %i\n",*nb_Elts);
  printf("nBlock %i\n",*nBlock);

  *nElBlock = (int*)malloc(sizeof(int)*(*nBlock));
  *typeBlock = (int*)malloc(sizeof(int)*(*nBlock));

  int block1 =0;
  for( int block=0;block<(*nBlock);block++) {
    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);
    r = fscanf(f, "%i",&IelType);  
    r = fscanf(f, "%i",&nv);
    //To use with Paraview
    if(IelType == 1 || IelType == 15) {/**nBlock=*nBlock-1;*//*block=block-1;*/*nb_Elts=*nb_Elts-nv;
       for(int s=0;s<nv*(1+elementType[IelType].nNodes);s++) {
         r = fscanf(f, "%i",&poub);
       }
    }
    else {
    (*nElBlock)[block1] = nv;
    (*typeBlock)[block1] = IelType;

    int size_el;
    size_el = elementType[IelType].nNodes;   

    for (int i = 0; i < nv; i++) {
      r = fscanf(f, "%i",&nEl2); 
      for(int jv=0;jv<size_el;jv++) { 
        r = fscanf(f, "%i",&poub);
      }

      if (r == EOF) 
        return EXIT_FAILURE;
    }
    block1++;
    }
   // free(connec);
  }  
  *nBlock = block1;
    printf("nb_Elts %i\n",*nb_Elts);
  return EXIT_SUCCESS;
}


int tabSearch2(int value,int** gnum,int gnum_size) {
 /* Déclarations */

 int POS;   /* position de la valeur */
 int I;     /* indice courant */
 int INF, MIL, SUP; /* limites du champ de recherche */
 
 /* Initialisation des limites du domaine de recherche */
 INF=0;
 SUP=gnum_size-1;
 /* Recherche de la position de la valeur */
 POS=INF;
 while (POS<SUP)
        {
          if(value == (*gnum)[POS]) break;
          POS++;
        }
  if(POS==SUP+1) POS=-1;
  return POS;
}



int tabSearch(int value,int** gnum,int gnum_size) {
 /* Déclarations */

 int POS;   /* position de la valeur */
 int I;     /* indice courant */
 int INF, MIL, SUP; /* limites du champ de recherche */
 
 /* Initialisation des limites du domaine de recherche */
 INF=0;
 SUP=gnum_size-1;
 /* Recherche de la position de la valeur */
 POS=-1;
 while ((INF<=SUP) && (POS==-1))
        {
         MIL=(SUP+INF)/2;
         //printf("value %i size %i MIL %i\n",value,MIL,gnum_size);
         if (value < (*gnum)[MIL])
               SUP=MIL-1;
         else if (value > (*gnum)[MIL])
               INF=MIL+1;
         else
               POS=MIL;
        }
 
  /* Edition du résultat */
 if (POS==-1) {
     printf("La valueeur recherchée ne se trouve pas "
            "dans le tableau. %i %i\n",value,gnum_size);
     return -1;
 }
 else {
     return POS;
 }
}


void tricroissant( int* a, int b) {
    int ind_min=0;
    int i=ind_min;
    int x=ind_min;
    int j=ind_min;
 
    for(i=ind_min;i<b;i++)
    {
        for(j=ind_min+1;j<b;j++)
        {
            if(a[i]<a[j])
            {
                x=a[i];
                a[i]=a[j];
                a[j]=x;
                j--;
                }
 
        }
 
        }
 
    x=a[ind_min];
    for(i=ind_min;i<b;i++)
    a[i]=a[i+1];
    a[b-1]=x;
 
}







/*----------------------------------------------------------------------
 *                                                                     
 * Read mesh dimension                                             
 *                                                                     
 * parameters:
 *   f                   <-- Mesh file                 
 *   dimension           --> Dimension                   
 *   nVtx             <-- number of vertices
 *   nElements           <-- number of elements
 *   nConnecVertex       <-- size of connectivity
 *   coords              --> vertices coordinates
 *   connecPointer       --> connectivity index  
 *   connec              --> connectivity
 *---------------------------------------------------------------------*/
 static int _read_mesh(FILE *f, 
                      int* nVtx, 
                      int* nb_Elts,
                      int* nBlock,
                      int** nElBlock,
                      int** typeBlock,
                      int** connec,
                      double *coords,
                      int *gnum_coords) {

  elType elementType [100];//= (elType*)malloc(sizeof(elType*)*100);
  
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
  r = fscanf(f, "%i",nVtx);
  printf("nVtx1 %i\n",*nVtx);

  int dimension =3;

  int i_el=0;
  for(int block = 1; block<(*nBlock)+1;block++) {
    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);
    r = fscanf(f, "%i",&nv);

    for (int i = 0; i < nv; i++) {
      //if(i_el>=0) {
      r = fscanf(f, "%i",gnum_coords+i_el);  
      r = fscanf(f, "%lf",coords+3*i_el);
      r = fscanf(f, "%lf",coords+3*i_el+1);
      r = fscanf(f, "%lf",coords+3*i_el+2);       
      printf("gnum_coords %i block %i nv %i x %f y %f z %f\n",
      gnum_coords[i_el],block,nv,coords[3*i_el],coords[3*i_el+1],coords[3*i_el+2]);
    //  }
      if (r == EOF) 
        return EXIT_FAILURE;
      i_el++;
    }
  }
//while(1==1){}
  int IelType,nEl2,poub;
  
  _goto(f,"$Elements");
  r = fscanf(f, "%i",nBlock);
  r = fscanf(f, "%i",nb_Elts);
  //printf("nb_Elts %i\n",*nb_Elts);
  //printf("nBlock %i\n",*nBlock);
    
  int** toto = *connec;
  
  int iop = *nb_Elts;
  
  int block1 = 0;

  for( int block=0;block<(*nBlock);block++) {
  

    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);
    r = fscanf(f, "%i",&IelType);  
    r = fscanf(f, "%i",&nv);
    //To use with Paraview
    if(IelType != 2 && IelType != 3) {*nb_Elts=*nb_Elts-nv;
       for(int s=0;s<nv*(1+elementType[IelType].nNodes);s++) {
         r = fscanf(f, "%i",&poub);
       }
    }
    else {
    (*nElBlock)[block1] = nv;

    (*typeBlock)[block1] = IelType;

    int size_el;
    if(IelType!=2 && IelType!=3) printf("IelType %i nv %i\n",IelType,nv);
    size_el = elementType[IelType].nNodes;   
    //  printf("IelType %i nv %i size_el %i block %i block1 %i nBlock %i\n",IelType,nv,size_el,block,block1,*nBlock); 
    for (int i = 0; i < nv; i++) {
      r = fscanf(f, "%i",&nEl2); 
      for(int jv=0;jv<size_el;jv++) { 
        r = fscanf(f, "%i",connec[block1]+size_el*i+jv);
        
        connec[block1][size_el*i+jv] =  1+tabSearch2(connec[block1][size_el*i+jv],&gnum_coords,*nVtx);
        if(connec[block1][size_el*i+jv]<-1 || connec[block1][size_el*i+jv] >100000)
           printf("Alerte wrong connectivity\n");
        //printf("connect %i \n",connec[block1][size_el*i+jv]);
        //printf("%i nEl2 %i block %i nv %i connec %i %i %i size_el %i IelType %i\n",
        //*nb_Elts,nEl2,block,nv,i,jv,connec[block1][size_el*i+jv],size_el,IelType);
      }
      
          
      if (r == EOF) 
        return EXIT_FAILURE;
    }
    block1++;
    }
   // free(connec);
  }  
  *nBlock = block1;
  //printf("nb_Elts %i *nBlock %i\n",*nb_Elts,*nBlock);

  return EXIT_SUCCESS;
}


int sizeForType(int type) {
  switch(type){
    case 1 : return 2;
    case 2 : return 3;
    case 3 : return 4;
    case 15: return 1;
    default: return 0;
  
  }
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
  

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  int n_partition = 0;
  const int two = 2;
  while(two * pow(n_partition, two) < comm_world_size) n_partition++;

  int n2 = (int) (two * pow(n_partition, two));




  /* Initialization
   * -------------- */

  int n_code_name = 0;
  char **codeNames = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;

  if (rank < 3 ) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] ="code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * 1);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }
   
  if( rank > 4) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] ="code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * 1);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }

  if(rank>=3 && rank<=4) {
    n_code_name = 2;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] ="code1";
    codeNames[1] ="code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
  }

  char* code_name = (char*)malloc(sizeof(char)*5);
  
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

  int currentRankA[n_code_name];
  int localCommSize[n_code_name];

  for (int i = 0; i < n_code_name; i++ ) {
    MPI_Comm_rank(localComm[i], &currentRankA[i]);
    MPI_Comm_size(localComm[i], &localCommSize[i]);
  }

 
  char cpl_id1[] = "cpl_code1_code2";

  int nb_part = 10;

    int** nElts;
    int** nBlock;
    int* nBlockOld;
    int*** typeBlock;
    int** nVtx;
    int**** eltsConnec;
    double*** coords;
    int**    gnum_coord;
    int*** nElBlock;

  char*** geofile =(char***)malloc(sizeof(char**)*n_code_name);

/**************************
  Loop on Code Names
*************************/

    nElts = (int**)malloc(sizeof(int*)*n_code_name);
    nVtx  = (int**)malloc(sizeof(int*)*n_code_name);
    eltsConnec = (int****)malloc(sizeof(int***)*n_code_name);
    nBlock = (int**)malloc(sizeof(int*)*n_code_name);
    nBlockOld = (int*)malloc(sizeof(int)*nb_part);
    typeBlock = (int***)malloc(sizeof(int**)*n_code_name);

    
    coords = (double***)malloc(sizeof(double**)*n_code_name);
    gnum_coord = (int**)malloc(sizeof(int*)*nb_part);
    nElBlock = (int***)malloc(sizeof(int**)*n_code_name);


  for (int i = 0; i < n_code_name; i++ ) {

   
   int localComm_size = localCommSize[i];
   code_name = codeNames[i];
   int currentRank = currentRankA[i];

   nElts[i] = (int*)malloc(sizeof(int)*nb_part);
   nVtx [i] = (int*)malloc(sizeof(int)*nb_part);
   eltsConnec[i] = (int***)malloc(sizeof(int**)*nb_part);
   coords[i] = (double**)malloc(sizeof(double*)*nb_part);  
   nElBlock[i] = (int**)malloc(sizeof(int*)*nb_part);
   nBlock[i] = (int*)malloc(sizeof(int)*nb_part);
   typeBlock[i] = (int**)malloc(sizeof(int*)*nb_part);   
   printf("%i rank %i localCommSize %i code_name %s\n",i,currentRank,localComm_size,code_name);
  
   if( code_name == "code1") {
     char* codeCpl="code2";
     CWP_Cpl_create (code_name, cpl_id1, codeCpl, CWP_COMM_PAR_WITH_PART,
                    CWP_GEOM_LOCATION, nb_part,
                    CWP_DISPLACEMENT_STATIC, CWP_FREQ_CPL_TIME_STEP);
     
   }

   if( code_name == "code2") {
      char* codeCpl="code1";
      CWP_Cpl_create (code_name, cpl_id1, codeCpl, CWP_COMM_PAR_WITH_PART,
                    CWP_GEOM_LOCATION, nb_part,
                    CWP_DISPLACEMENT_STATIC, CWP_FREQ_CPL_TIME_STEP);
     
   }

  }//Loop on codes

  MPI_Barrier(MPI_COMM_WORLD);
  
/*******************************
          Loop on Code Names
*********************************/

  for (int i = 0; i < n_code_name; i++ ) {

   int localComm_size = localCommSize[i];
   code_name = codeNames[i];
   int currentRank = currentRankA[i];
       
    CWP_Visu_set(code_name, cpl_id1,1,Ensight,"text");


    char* geoModelfile;
    if(code_name == "code1") {
      char* geoModelfile2 = "meshes/sphere01";
      geoModelfile=geoModelfile2;
    }

    if(code_name == "code2") {
      char* geoModelfile2 = "meshes/sphere02";
      geoModelfile=geoModelfile2;
    }

    if(currentRank==0) {  
      char* s =(char*)malloc(sizeof(char)*500);
      sprintf(s,"python meshes/gmsh_mesh.py %i %s",nb_part,geoModelfile);
      system(s);
      
       system("cp meshes/sphere01Model.geo meshes/sphere02Model.geo");
    }

    geofile[i] =(char**)malloc(sizeof(char*)*nb_part);
    
    for(int i_part =0;i_part<nb_part;i_part++){

      geofile[i][i_part] =(char*)malloc(sizeof(char)*550); 
      sprintf(geofile[i][i_part],"%s_part%i.geo",geoModelfile,i_part);
      printf("geofile %s %i\n",geofile[i][i_part],localCommSize[i]);
      if(currentRank==0 /*&& code_name=="code1"*/) {
       _generate_gmsh_mesh(geofile[i][i_part],localComm_size,1);
      }
  /*    MPI_Barrier(MPI_COMM_WORLD);
      if(currentRank==0 && code_name=="code2") {
       _generate_gmsh_mesh(geofile[i][i_part],localComm_size,1);
      }
      MPI_Barrier(MPI_COMM_WORLD);*/
    }//Loop on i_part
    printf("codename %s currentRank %i i %i rank %i\n",code_name,currentRank,i,rank);
  }//Loop on codes

  MPI_Barrier(MPI_COMM_WORLD);

/*******************************
          Loop on Code Names
*********************************/
  for (int i = 0; i < n_code_name; i++ ) {


    int localComm_size = localCommSize[i];
    code_name = codeNames[i];
    int currentRank = currentRankA[i];

    for(int i_part =0;i_part<nb_part;i_part++){
   
    int len=0;
    for(len = 0; geofile[i][i_part][len] != '\0'; ++len);

    int lenfile;

    if((1+currentRankA[i])>=1000 && (1+currentRankA[i])<10000) lenfile=len+5;
    else if((1+currentRankA[i])>=100) lenfile=len+4;
    else if((1+currentRankA[i])>=10)  lenfile=len+3;  
    else lenfile=len+2;

    if(localCommSize[i]==1) {
     lenfile ==len;
    }


    char* root=(char*)malloc(sizeof(char)*(len-4));
  
    printf("filename %s %i\n",geofile[i][i_part],len);

    for (int j=0;j<len-4;j++) root[j]=geofile[i][i_part][j];   

    char* filename=(char*)malloc(sizeof(char)*(lenfile));  
    char* filename2=(char*)malloc(sizeof(char)*(10+2*lenfile));  

    if(localCommSize[i]!=1) sprintf(filename,"%s_%i.msh",root,1+currentRank);
    else sprintf(filename,"%s.msh",root);

    printf("filename %s %i\n",filename,currentRankA[i]);

    FILE* meshFile;
    meshFile = fopen(filename, "r");

    eltsConnec[i][i_part] = NULL;        // Connectivity
    coords[i][i_part]  = NULL; 
    nElBlock[i][i_part] = NULL;
    typeBlock[i][i_part] = NULL;

    _read_mesh_dim(meshFile,
              &(nVtx[i][i_part]), 
              &(nElts[i][i_part]),
              &(nBlock[i][i_part]),
              &(nElBlock[i][i_part]),
              &(typeBlock[i][i_part]));


    coords[i][i_part] = (double*) malloc(3 * nVtx[i][i_part] * sizeof(double));
    gnum_coord[i_part] = (int*) malloc(nVtx[i][i_part] * sizeof(int));
    eltsConnec[i][i_part] = (int**)malloc(nBlock[i][i_part]*sizeof(int*));
  
    printf("currentRank %i nVtx[i][i_part] %i nElts[i][i_part] %i nBlock[i][i_part] %i nElBlock[i][i_part][0] %i typeBlock[i][i_part][0] %i\n",
    currentRankA[i], nVtx[i][i_part],nElts[i][i_part],nBlock[i][i_part],nElBlock[i][i_part][0],typeBlock[i][i_part][0]);
    nBlockOld[i_part] = nBlock[i][i_part];
    for(int b=0;b<nBlock[i][i_part];b++) {
      //printf("Partition %i Block %i size = %i rank %i\n",i_part,b,nElBlock[i][i_part][b],rank);
      eltsConnec[i][i_part][b]=(int*)malloc(sizeForType(typeBlock[i][i_part][b])*nElBlock[i][i_part][b]*sizeof(int));
    }

    fclose(meshFile);    
    meshFile = fopen(filename/*"carre01.msh"*/, "r");


    int** connec_p=eltsConnec[i][i_part];    

    printf("filenamePP %s %i %i %i %i %i\n",filename,rank,nElts[i][i_part],nBlock[i][i_part],typeBlock[i][i_part][0],nVtx[i][i_part]);

      _read_mesh(meshFile,
                &nVtx[i][i_part], 
                &nElts[i][i_part],
                &nBlock[i][i_part],
                &(nElBlock[i][i_part]),
                &(typeBlock[i][i_part]),
                eltsConnec[i][i_part],
                coords[i][i_part],
                gnum_coord[i_part]);


    fclose(meshFile);
    free(filename);
    free(filename2);
    free(root);

    for(int b=0;b<nBlock[i][i_part];b++) {
      for(int j=0;j<nElBlock[i][i_part][b];j++) {
        printf("eltsConnec[%i][%i][%i][%i] rank %i %i\n",i,i_part,b,j,currentRank,eltsConnec[i][i_part][b][j]);
      }
    }

    printf("CWP_Mesh_interf_vtx_set %i\n",nVtx[i][i_part]);           
    CWP_Mesh_interf_vtx_set(codeNames[i], cpl_id1,i_part,nVtx[i][i_part],coords[i][i_part],NULL);                
    printf("After CWP_Mesh_interf_vtx_set %i %i\n",i_part,rank);  
    

   }//Loop on partition     


  int nb_elt_type =100;

  int** nb_block_type = (int**) malloc(nb_part*sizeof(int*));
  int** nb_block_type_max = (int**) malloc(nb_part*sizeof(int*));
   
  for(int i_part =0;i_part <nb_part;i_part++) {
    nb_block_type[i_part] = (int*) malloc(nb_elt_type*sizeof(int));
    nb_block_type_max[i_part] = (int*) malloc(nb_elt_type*sizeof(int));
    for(int l=0; l<nb_elt_type;l++) nb_block_type[i_part][l] = 0;
    
    for(int i_block =0;i_block <nBlock[i][i_part];i_block++) { 
        nb_block_type[i_part][typeBlock[i][i_part][i_block]]++;
    }
    MPI_Allreduce(nb_block_type[i_part], nb_block_type_max[i_part], nb_elt_type, MPI_INT, MPI_MAX,
                  localComm[i]);

    int nBlock_correct =0;
    for(int l=0; l<nb_elt_type;l++) nBlock_correct+=nb_block_type_max[i_part][l];

    int nBlock_old = nBlock[i][i_part];
    nBlock[i][i_part]=nBlock_correct;

    printf("nBlock_old %i %i\n",nBlock_old,nBlock_correct);
    
    int* nElBlockNew   =  (int*) malloc  ( nBlock[i][i_part]*sizeof(int)  );
    int* typeBlockNew  =  (int*) malloc  ( nBlock[i][i_part]*sizeof(int)  );
    int**eltsConnecNew =  (int**) malloc ( nBlock[i][i_part]*sizeof(int*) );

    int ind_block = 0;
    
    int* ibblock = (int*) malloc(nb_elt_type*sizeof(int));

    for(int bb=0; bb < nb_elt_type; bb++) ibblock[bb]=0;
    
    for(int l=0; l<nb_elt_type;l++) {
      for(int l2=0;l2<nb_block_type_max[i_part][l];l2++) {
        printf("ind_block %i %i \n",ind_block,i_part);
        if(l2<nb_block_type[i_part][l]){
        
          
          while(typeBlock[i][i_part][ ibblock[l] ] != l) {
            ibblock[l]++;
          }
          
          nElBlockNew [ind_block] = nElBlock[i][i_part][ ibblock[l] ];
          typeBlockNew[ind_block] = typeBlock[i][i_part][ ibblock[l] ];
          eltsConnecNew[ind_block]= eltsConnec[i][i_part][ ibblock[l] ];
          ibblock[l]++;
        }
        else {
          nElBlockNew [ind_block] = 0;
          typeBlockNew[ind_block] = l;  
          eltsConnecNew[ind_block] = (int*)malloc(sizeForType(typeBlockNew[ind_block])*nElBlockNew [ind_block]*sizeof(int)); 
        }
        
        ind_block++;
      }//end l2 loop
    }//end l loop
 
    realloc(nElBlock[i][i_part],0);

    nElBlock[i][i_part] = nElBlockNew;

    realloc(typeBlock[i][i_part],0);     
    typeBlock[i][i_part] = typeBlockNew;

    realloc(eltsConnec[i] [i_part],0);     
    eltsConnec[i][i_part] =eltsConnecNew;

 } // LOOP on i_part


  CWP_g_num_t*** GNUM=(CWP_g_num_t***) malloc(sizeof(CWP_g_num_t**)*nb_part);
            
  
  for(int i_part =0;i_part <nb_part;i_part++) {
 
    nElts[i][i_part]=0;
    GNUM[i_part] = (CWP_g_num_t**) malloc(sizeof(CWP_g_num_t*)*nBlock[i][i_part]);

    for(int i_block =0;i_block <nBlock[i][i_part];i_block++) { 
      GNUM[i_part][i_block] = (CWP_g_num_t*) malloc(sizeof(CWP_g_num_t)*nElBlock[i][i_part][i_block]);
      printf("Standard Block Add\n");
      CWP_Block_t cwp_block_t;
      if(typeBlock[i][i_part][i_block] == 1) cwp_block_t = CWP_BLOCK_EDGE2;
      if(typeBlock[i][i_part][i_block] == 2) cwp_block_t = CWP_BLOCK_FACE_TRIA3;
      if(typeBlock[i][i_part][i_block] == 3) cwp_block_t = CWP_BLOCK_FACE_QUAD4;
      if(typeBlock[i][i_part][i_block] == 15) cwp_block_t = CWP_BLOCK_NODE;
  
      int block_id = CWP_Mesh_interf_block_add(code_name, cpl_id1,cwp_block_t);

      printf("Standard block set %i\n",i_block);
      printf("\n");
      for(int g=0;g<nElBlock[i][i_part][i_block];g++) {
        if(eltsConnec[i][i_part][i_block][g]<0 || eltsConnec[i][i_part][i_block][g]>1000) printf("YYYY %i %i %i %i\n",i_part,i_block,g,eltsConnec[i][i_part][i_block][g]);
        GNUM[i_part][i_block][g]=-456;
      }
      CWP_Mesh_interf_block_std_set(code_name,cpl_id1,i_part,block_id,nElBlock[i][i_part][i_block],eltsConnec[i][i_part][i_block],NULL/*GNUM[i_part][i_block]*/);
      nElts[i][i_part]+=nElBlock[i][i_part][i_block];

    }//end i_block loop
    
               
  }//Loop on i_part

  }//Loop on codes
 
  MPI_Barrier(MPI_COMM_WORLD);
       
/*******************************
          Loop on Code Names
*********************************/

  double*** rank_data     = (double***)malloc(sizeof(double**)*n_code_name);
  double*** rank_data_vtx = (double***)malloc(sizeof(double**)*n_code_name);


  for (int i = 0; i < n_code_name; i++ ) {

    int localComm_size = localCommSize[i];
    code_name = codeNames[i];
    int currentRank = currentRankA[i];


    CWP_Mesh_interf_finalize(code_name,cpl_id1);

    CWP_Status_t visu_status = CWP_STATUS_OFF;


    
    if(code_name == "code1"){

      /* Rank field */
      CWP_Field_create (code_name,cpl_id1,"rank",CWP_DOUBLE,CWP_FIELD_STORAGE_BLOCK,1,
                        CWP_FIELD_VALUE_CELL_POINT,
                        CWP_FIELD_EXCH_SEND,
                        visu_status);
                        
      CWP_Field_create (code_name,cpl_id1,"rank_vtx",CWP_DOUBLE,CWP_FIELD_STORAGE_BLOCK,1,
                        CWP_FIELD_VALUE_NODE,
                        CWP_FIELD_EXCH_SEND,
                        visu_status);                        
                        
   } //end  if(code_name == "code1")
   
  
   if(code_name == "code2"){
              
    CWP_Field_create (code_name,cpl_id1,"rank",CWP_DOUBLE,CWP_FIELD_STORAGE_BLOCK,1,
                      CWP_FIELD_VALUE_CELL_POINT,
                      CWP_FIELD_EXCH_RECV,
                      visu_status);       
                      
    CWP_Field_create (code_name,cpl_id1,"rank_vtx",CWP_DOUBLE,CWP_FIELD_STORAGE_BLOCK,1,
                      CWP_FIELD_VALUE_NODE,
                      CWP_FIELD_EXCH_RECV,
                      visu_status);      
                                         
   }

   rank_data_vtx[i] = (double**)malloc(sizeof(double*)*nb_part);
   rank_data    [i] = (double**)malloc(sizeof(double*)*nb_part);
   
   for(int i_part =0;i_part <nb_part;i_part++) {    
      rank_data    [i][i_part] = (double*)malloc(sizeof(double)*nElts[i][i_part]);
      rank_data_vtx[i][i_part] = (double*)malloc(sizeof(double)*nVtx [i][i_part]);
      CWP_Field_data_set(code_name,cpl_id1,"rank"    ,i_part,rank_data    [i][i_part]);
      CWP_Field_data_set(code_name,cpl_id1,"rank_vtx",i_part,rank_data_vtx[i][i_part]);
   }   

   printf("After data set\n");

  }//Loop on codes


  MPI_Barrier(MPI_COMM_WORLD);
  
  printf("Before CWP_Geom_compute\n");
  
/*******************************
          Loop on Code Names
*********************************/
for (int i = 0; i < n_code_name; i++ ) {

  int* n_uncomputed_tgt;

  //TODO: Calcul géom piloté par la nature des champs
//Erreur si on crée un champ après le calcul
  code_name = codeNames[i];
  CWP_Geom_compute(code_name,cpl_id1, CWP_FIELD_VALUE_NODE, n_uncomputed_tgt);
  CWP_Geom_compute(code_name,cpl_id1, CWP_FIELD_VALUE_CELL_POINT, n_uncomputed_tgt);
  //Argument tag points localisé oui/non + print 

 // printf("After CWP_Geom_compute code %i\n",i);
 
  }//Loop on codes
  //while(1==1){}
  MPI_Barrier(MPI_COMM_WORLD);

  double recv_time = 0.150;

  for (int i_time=0;i_time<1;i_time++) {

/*******************************
          Loop on Code Names
*********************************/
  for (int i = 0; i < n_code_name; i++ ) {

    int localComm_size = localCommSize[i];
    code_name = codeNames[i];
    int currentRank = currentRankA[i];  
    
    printf("Before CWP_next_recv_time_set\n");
    CWP_next_recv_time_set(code_name,cpl_id1,recv_time);
    printf("After CWP_next_recv_time_set\n");

    if(code_name == "code1"){
      printf("TEST %i %s %s\n",rank,code_name,codeNames[i]);
      for(int i_part =0;i_part <nb_part;i_part++) {     
        for(int j=0;j<nElts[i][i_part];j++)
          rank_data[i][i_part][j]= rank;

        for(int j=0;j<nVtx[i][i_part];j++)
          rank_data_vtx[i][i_part][j]= coords[i][i_part][3*j];          
      }

      printf("CWP_Issend at %f\n",recv_time);
      CWP_Issend (code_name,cpl_id1,"rank");    
      CWP_Issend (code_name,cpl_id1,"rank_vtx");    
   
    }
    else { 
      CWP_Irecv (code_name,cpl_id1,"rank");
      CWP_Irecv (code_name,cpl_id1,"rank_vtx");
    }

  }// end code loop  

/*******************************
          Loop on Code Names
*********************************/
  for (int i = 0; i < n_code_name; i++ ) {
    int localComm_size = localCommSize[i];
    code_name = codeNames[i];
    int currentRank = currentRankA[i];  
    
    if(code_name == "code1"){
      CWP_Wait_issend (code_name,cpl_id1,"rank");
      CWP_Wait_issend (code_name,cpl_id1,"rank_vtx");
    }
    else { 
      CWP_Wait_irecv (code_name,cpl_id1,"rank"    );
      CWP_Wait_irecv (code_name,cpl_id1,"rank_vtx");
    }

   }//end codename loop
     
   recv_time+=0.1;  
  }// end i_time loop  

/*******************************
          Loop on Code Names
*********************************/
  for (int i = 0; i < n_code_name; i++ ) {

    int localComm_size = localCommSize[i];
    code_name = codeNames[i];
    int currentRank = currentRankA[i];  
 
  CWP_Mesh_interf_del(code_name, cpl_id1); 
  CWP_Cpl_del (code_name, cpl_id1);
           
 }//Loop on codenames

  MPI_Barrier(MPI_COMM_WORLD);


/******************************************/ 

  fflush(stdout);

  CWP_Finalize();
  MPI_Finalize();

  free (localComm);
  free (codeNames);
  free (is_coupled_rank);
  free (times_init);

  return 0;
}





