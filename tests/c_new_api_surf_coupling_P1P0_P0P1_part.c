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
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cwp.h"
#include "grid_mesh.h"


/*----------------------------------------------------------------------
 *
 * Dump status exchange
 *
 * parameters:
 *   status              <-- Exchange status
 *---------------------------------------------------------------------*/

static void _dumpStatus(FILE* outputFile, CWP_Status_t status)
{
  switch(status) {
  case CWP_STATUS_ON :
    fprintf(outputFile, "Exchange Ok\n");
    break;
  case CWP_STATUS_OFF :
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

/* static void _dumpNotLocatedPoints(FILE* outputFile, */
/*                                   const char *coupling_id, */
/*                                   const int nNotLocatedPoints) */
/* { */
/*   if ( nNotLocatedPoints > 0) { */
/*     fprintf(outputFile, "Not located points :\n"); */
/*     const int* notLocatedPoints = cwipi_get_not_located_points(coupling_id); */
/*     for(int i = 0; i < nNotLocatedPoints; i++) */
/*      fprintf(outputFile, "%i ", notLocatedPoints[i]); */
/*     fprintf(outputFile, "\n"); */
/*   } */
/* } */

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
           double       *randLevel)
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
    i++;
  }
}

//double *shapef = NULL;

static void _userInterpolation(const int                   interface_type,
                               const int                   n_src_vtcs,
                               const int                   n_src_std_elts,
                          //   const int                   n_src_poly,
                               const int                   n_tgt_pts,
                               const double                src_vtcs_coords[],
                               const int                   src_connec_idx[],
                               const int                   src_connec[],
                               const double                tgt_pts_coords[],
                               const int                   tgt_pts_target_location[],
                               const float                 tgt_pts_dist[],
                               const int                   tgt_pts_bary_coords_idx[],
                               const double                tgt_pts_bary_coords[],
                               const int                   stride,
                               const CWP_Dof_location_t     src_field_location,
                               const void                 *src_field,
                               const CWP_Dof_location_t     tgt_field_location,
                               void                       *tgt_field
                              )
{

  // Compute shapef

  int compute_shape_f = 0;

  if ( src_field_location == CWP_DOF_LOCATION_NODE ) {

    //printf("n_tgt_pts USER_INTERP %i\n",n_tgt_pts);

    for (int i = 0; i < n_tgt_pts; i++) {

      int ielt = tgt_pts_target_location[i];
      int ivertex[4];

      ivertex[0] = src_connec[src_connec_idx[ielt]  ] ;
      ivertex[1] = src_connec[src_connec_idx[ielt]+1] ;
      ivertex[2] = src_connec[src_connec_idx[ielt]+2] ;
      ivertex[3] = src_connec[src_connec_idx[ielt]+3] ;

     // if (shapef == NULL) {
        double *shapef = (double *) malloc(4 * n_src_std_elts * sizeof(double));
        compute_shape_f = 1;
     // }
      double *shapef_elt = shapef + 4 * ielt;

      //
    //  printf(" Compute shape function %i %i\n",i,n_tgt_pts);
      //

  /*    if (compute_shape_f == 1) {
        double deriv[4][2];
        double uv[2];
        double a[2][2];
        double b[2];
        double det_a;
        double x[2];
        double inv_a[2][2];

        for (int k = 0; k < 2; k++)
          uv[k] = 0.5;


        const int it_max = 100;
        for (int it = 0; it < it_max; it++) {

          shapef_elt[0] = (1 - uv[0]) * (1 - uv[1]);
          shapef_elt[1] = uv[0] * (1 - uv[1]);
          shapef_elt[2] = uv[0] * uv[1];
          shapef_elt[3] = (1 - uv[0]) * uv[1];

          deriv[0][0] = - (1 - uv[1]);
          deriv[0][1] = - (1 - uv[0]);
          deriv[1][0] =   (1 - uv[1]);
          deriv[1][1] = - uv[0];
          deriv[2][0] =   uv[1];
          deriv[2][1] =   uv[0];
          deriv[3][0] = - uv[1];
          deriv[3][1] =   (1 - uv[0]);

          for (int k = 0; k < 2; k++) {
            for (int l = 0; l < 2; l++)
              a[k][l] = 0.0;
          }

          b[0] = - tgt_pts_coords[3 * i    ];
          b[1] = - tgt_pts_coords[3 * i + 1];

          for (int k = 0; k < 4; k++) {

            b[0] += (shapef_elt[k] * src_vtcs_coords[3 * ivertex[k]   ]);
            b[1] += (shapef_elt[k] * src_vtcs_coords[3 * ivertex[k] +1]);

            for (int l = 0; l < 2; l++) {
              a[0][l]  -=  (deriv[k][l] * src_vtcs_coords[3 * ivertex[k]    ]);
              a[1][l]  -=  (deriv[k][l] * src_vtcs_coords[3 * ivertex[k] + 1]);
            }
          }

          det_a = a[0][0] * a[1][1] - a[0][1] * a[1][0];
          if (fabs(det_a) < 1e-12) {
            printf("matrice non inversible\n");
            exit(1);
          }

          double det_inv = 1./det_a;

          inv_a[0][0] =   det_inv * a[1][1];
          inv_a[0][1] = - det_inv * a[0][1];
          inv_a[1][0] =   det_inv * a[1][0];
          inv_a[1][1] =   det_inv * a[0][0];

          x[0] = inv_a[0][0] * b[0] + inv_a[0][1] * b[1];
          x[1] = inv_a[1][0] * b[0] + inv_a[1][1] * b[1];

          double dist = 0.0;

          for (int k = 0; k < 2; k++) {
            dist += x[k] * x[k];
            uv[k] += x[k];
          }

          if (dist <= 1e-5)
            break;
        }

        shapef_elt[0] = (1 - uv[0]) * (1 - uv[1]);
        shapef_elt[1] = uv[0] * (1 - uv[1]);
        shapef_elt[2] = uv[0] * uv[1];
        shapef_elt[3] = (1 - uv[0]) * uv[1];

      }
*/
      //
      // Insterpolation
      //
    //  printf("n_tgt_pts USER_INTERP %i\n",n_tgt_pts);
      ((double *) tgt_field)[i] = 0.0;
     //  while(1==1){}
      for (int k = 0; k < 4; k++) {
       ((double *) tgt_field)[i] = ((double *) src_field)[0];//ivertex[k]];
      }

      free(shapef);

    }
  }
  else {
    printf("Error in _userInterpolation : bad solver_type\n");
    exit(EXIT_FAILURE);
  }
}



/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
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

  srand(rank+time(0));

  int n_partition = 0;
  while(1.5 * (double)pow(n_partition, 2) < commWorldSize) n_partition++;

  if(pow(n_partition, 2) > commWorldSize) n_partition--;

  int size_code1 = (int)(pow(n_partition, 2));
  int size_code2 = (int)(pow(n_partition, 2));

  /* Read args from command line
   * --------------------------- */

  int nVertexSeg = 100;
  double randLevel = 0.4;

  _read_args(argc, argv, &nVertexSeg, &randLevel);

  /* Initialization
   * -------------- */

  int n_code_name = 0;
  char **codeName = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;

  char **codeCoupledName = NULL;
  int codeId;

  int rankCode1[ commWorldSize];
  int rankCode2[ commWorldSize];


  //Random distribution of codes on processes.
  int nsize = commWorldSize;
  int nranks[nsize];

  if(rank==0){
    time_t t;
    /* Intializes random number generator */
    srand((unsigned) time(&t));
    int size_code1_domain = size_code1;
    int size_code2_domain = size_code2;
    printf("size_code1 %i\n",size_code1);
    for(int i=0; i < commWorldSize; i++) {
      rankCode1[i] = 0;
      rankCode2[i] = 0;
    }
    for(int i=0; i < size_code1_domain; i++){
      int tmp = rand() % commWorldSize;
      while(rankCode1[tmp] == 1){
        tmp = rand() % commWorldSize;
      }
      rankCode1[tmp]=1;
    }

    for(int i=0; i < size_code2_domain; i++){
      int tmp = rand() % commWorldSize;
      while(rankCode2[tmp] == 1){
        tmp = rand() % commWorldSize;
      }
      rankCode2[tmp]=1;
    }

    for(int i=0; i < commWorldSize; i++) {
      if(rankCode1[i] == 0 && rankCode2[i] == 0){
        int t = rand()%commWorldSize;
        while(!(rankCode1[t] == 1 && rankCode2[t] == 1)){
          t = rand()%commWorldSize;
        }
        int h12 = rand()%2;
        if(h12 ==0) {
          rankCode1[t] = 0;
          rankCode1[i] = 1;
        }
        else {
          rankCode2[t] = 0;
          rankCode2[i] = 1;
        }
      }
    }

 /*  rankCode1[0] = 1;
   rankCode1[1] = 1;
   rankCode1[2] = 1;
   rankCode1[3] = 1;
   rankCode1[4] = 0;
   rankCode1[5] = 0;

   rankCode2[0] = 0;
   rankCode2[1] = 0;
   rankCode2[2] = 1;
   rankCode2[3] = 1;
   rankCode2[4] = 1;
   rankCode2[5] = 1;

 */
    int ind=0;
    for(int i=0; i < commWorldSize; i++) {
      if(rankCode2[i] || rankCode1[i]){
       nranks[ind]=i;
       ind++;
      }
    }



  } // if rank==0


  int nsizercv=nsize;

  int rankCode1Rcv = -1;
  int rankCode2Rcv = -1;

  MPI_Scatter(&rankCode1    ,1 ,MPI_INT,
              &rankCode1Rcv,1 ,MPI_INT,
              0,
              MPI_COMM_WORLD
             );

  MPI_Scatter(&rankCode2    ,1 ,MPI_INT,
              &rankCode2Rcv,1 ,MPI_INT,
              0,
              MPI_COMM_WORLD
             );

  MPI_Scatter(&rankCode2    ,1 ,MPI_INT,
              &rankCode2Rcv,1 ,MPI_INT,
              0,
              MPI_COMM_WORLD
             );



   if(rankCode1Rcv && rankCode2Rcv) {
    printf("I am process %i and I own the code1 and the code2.\n",rank);
  }
  else {
    if(rankCode1Rcv){
      printf("maxI am process %i and I own the code1.\n",rank);
    }
    if(rankCode2Rcv){
      printf("maxI am process %i and I own the code2.\n",rank);
    }
  }


  int* nbPartSeg = NULL;
  int* nbPart = NULL;

  if (rankCode1Rcv) {
    n_code_name = 1;
    codeName = malloc(sizeof(char*)*n_code_name);
    codeCoupledName = malloc(sizeof(char*)*n_code_name);
    codeName[0] = "code1";
    codeCoupledName[0] = "code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    nbPartSeg = (int*)malloc(sizeof(int)*n_code_name);
    nbPartSeg[0] = 1;
    nbPart = (int*)malloc(sizeof(int)*n_code_name);
    nbPart[0] = nbPartSeg[0]*nbPartSeg[0];
  }

  if (rankCode2Rcv) {
    n_code_name = 1;
    codeName = malloc(sizeof(char*)*n_code_name);
    codeCoupledName = malloc(sizeof(char*)*n_code_name);
    codeName[0] = "code2";
    codeCoupledName[0] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    nbPartSeg = (int*)malloc(sizeof(int)*n_code_name);
    nbPartSeg[0] = 2;
    nbPart = (int*)malloc(sizeof(int)*n_code_name);
    nbPart[0] = nbPartSeg[0]*nbPartSeg[0];
  }

  if (rankCode1Rcv && rankCode2Rcv) {
    n_code_name = 2;
    codeName        = malloc(sizeof(char*)*n_code_name);
    codeCoupledName = malloc(sizeof(char*)*n_code_name);
    codeName[0] = "code1";
    codeName[1] = "code2";
    codeCoupledName[0] = "code2";
    codeCoupledName[1] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
    is_coupled_rank[1] = CWP_STATUS_ON;
    nbPartSeg = (int*)malloc(sizeof(int)*n_code_name);
    nbPartSeg[0] = 1;
    nbPartSeg[1] = 2;
    nbPart = (int*)malloc(sizeof(int)*n_code_name);
    nbPart[0] = nbPartSeg[0]*nbPartSeg[0];
    nbPart[1] = nbPartSeg[1]*nbPartSeg[1];
  }

  char* fileName = (char *) malloc(sizeof(char) * 35);
  sprintf(fileName,"c_surf_coupling_P1P0_P0P1_%4.4d.txt",rank);

  outputFile = fopen(fileName,"w");

  free(fileName);
  //CWP_Output_file_set (outputFile);


  times_init = malloc(sizeof(double) * n_code_name);

  for (int i = 0; i < n_code_name; i++) {
    times_init[i] = 0;
  }


  MPI_Comm *localComm = malloc(sizeof(MPI_Comm)*n_code_name);
  CWP_Init(MPI_COMM_WORLD,
           n_code_name,
           (const char **) codeName,
           is_coupled_rank,
           times_init,
           localComm);

  /* Output redirection
   * ------------------ */

  int currentRank[n_code_name];
  int localCommSize[n_code_name];

  for(int i_code = 0; i_code < n_code_name; i_code++) {
    if(localComm[i_code]!=MPI_COMM_NULL) MPI_Comm_rank(localComm[i_code], &currentRank[i_code]);
    if(localComm[i_code]!=MPI_COMM_NULL) MPI_Comm_size(localComm[i_code], &localCommSize[i_code]);
  }

  fprintf(outputFile, "  Surface coupling test : P1P0_P0P1 with polygon\n");
  fprintf(outputFile, "\n");

  fprintf(outputFile, "\nDump after initialization\n");
  fprintf(outputFile, "-------------------------\n");


  /* Coupling creation
   * ----------------- */
//  if (rank == 0)
  //  printf("        Create coupling %i\n",rank);






  for(int i_code = 0; i_code < n_code_name; i_code++) {
    printf("        Create coupling %i code %i\n",rank,i_code);
    CWP_Cpl_create (codeName[i_code],             // Code name
                    "c_new_api_surf_cpl_P1P0_P0P1_part",  // Coupling id
                    codeCoupledName[i_code],         // Coupled application id
                    CWP_COMM_PAR_WITH_PART,  // Coupling type
                    CWP_SPATIAL_INTERP_FROM_LOCATION,       // Solver type
                    nbPart[i_code],          // Partition number
                    CWP_DYNAMIC_MESH_STATIC, // Mesh type
                    CWP_TIME_EXCH_ASYNCHRONOUS); // frequency
    printf("   After     Create coupling %i code %i\n",rank,i_code);
  }

  if (rank == 0)
    printf("        Set visu\n");

  for(int i_code = 0; i_code < n_code_name; i_code++) {
    CWP_Visu_set(codeName[i_code],                        // Code name
                 "c_new_api_surf_cpl_P1P0_P0P1_part",     // Coupling id
                 1,           // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT,     // Postprocessing format
                 "text");     // Postprocessing option
  }

  /* Mesh definition
   * --------------- */

  if (rank == 0)
    printf("        Create mesh\n");

  int** nVertex = NULL;                   // Number of vertex
  double ***coords = NULL;         // Vertex coordinates
  int** nElts = NULL;                     // Number of elements
  int ***eltsConnecPointer = NULL; // Connectivity index
  int ***eltsConnec = NULL;        // Connectivity

  /* Domain bounds */

  const double xmin = -10;
  const double xmax =  10;
  const double ymin = -10;
  const double ymax =  10;

  coords = (double ***) malloc(sizeof(double**) * n_code_name );
  eltsConnecPointer = (int ***) malloc(sizeof(int**) * n_code_name);
  eltsConnec = (int ***) malloc(sizeof(int**) * n_code_name);
  nVertex = (int **) malloc(sizeof(int*) * n_code_name);
  nElts = (int **) malloc(sizeof(int*) * n_code_name);
  srand (time(NULL));


  for(int i_code = 0; i_code < n_code_name; i_code++) {
    coords           [i_code] = (double **) malloc(sizeof(double*) * nbPart[i_code] );
    eltsConnecPointer[i_code] = (int **)    malloc(sizeof(int*) * nbPart[i_code]);
    eltsConnec       [i_code] = (int **)    malloc(sizeof(int*) * nbPart[i_code]);
    nVertex[i_code] = (int *)    malloc(sizeof(int) * nbPart[i_code]);
    nElts  [i_code] = (int *)    malloc(sizeof(int) * nbPart[i_code]);
  }


  for(int i_code = 0; i_code < n_code_name; i_code++) {
    randLevel =0.4;
    double xSegPart = (xmax-xmin)/(double)nbPartSeg[i_code];
    double ySegPart = (ymax-ymin)/(double)nbPartSeg[i_code];
    for(int u=0;u<nbPartSeg[i_code];u++) {
      for(int v=0;v<nbPartSeg[i_code];v++) {
        int i_part = nbPartSeg[i_code] * u + v;
        int nVertexSegPart = nVertexSeg/nbPartSeg[i_code];
        nVertex[i_code][i_part] = nVertexSegPart * nVertexSegPart;
        nElts[i_code][i_part] = (nVertexSegPart-1) * (nVertexSegPart-1);

        coords           [i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nVertex[i_code][i_part] );
        eltsConnecPointer[i_code][i_part] = (int *)    malloc(sizeof(int) * (nElts[i_code][i_part] + 1));
        eltsConnec       [i_code][i_part] = (int *)    malloc(sizeof(int) * 4 * nElts[i_code][i_part]);
        if(codeName[i_code]=="code1") {
          double xminPart = xmin + xSegPart * v;
          double yminPart = ymin + ySegPart * u;
          double xmaxPart = xminPart + xSegPart;
          double ymaxPart = yminPart + ySegPart;
          printf("INFO rank %i i_part %i x %f %f y %f %f nbPart %i test %i nVertex %i \n",
          rank,i_part,xminPart,yminPart,xmaxPart,ymaxPart,nbPart[i_code],nVertexSeg/nbPartSeg[i_code],nVertex[i_code][i_part]);

          grid_mesh(xminPart,
                    xmaxPart,
                    yminPart,
                    ymaxPart,
                    randLevel,
                    nVertexSegPart,
                    sqrt(localCommSize[i_code]),
                    coords           [i_code][i_part],
                    eltsConnecPointer[i_code][i_part],
                    eltsConnec       [i_code][i_part],
                    localComm        [i_code]);
        } //if code1
      }//loop on v
    }//loop on u
  } // loop on codes



  MPI_Barrier(MPI_COMM_WORLD);
   for(int i_code = 0; i_code < n_code_name; i_code++) {
    randLevel =0.2;
    double xSegPart = (xmax-xmin)/(double)nbPartSeg[i_code];
    double ySegPart = (ymax-ymin)/(double)nbPartSeg[i_code];
    for(int u=0;u<nbPartSeg[i_code];u++) {
      for(int v=0;v<nbPartSeg[i_code];v++) {
        int i_part = nbPartSeg[i_code] * u + v;
        int nVertexSegPart = nVertexSeg/nbPartSeg[i_code];
        if(codeName[i_code]=="code2") {
          double xminPart = xmin + xSegPart * v;
          double yminPart = ymin + ySegPart * u;
          double xmaxPart = xminPart + xSegPart;
          double ymaxPart = yminPart + ySegPart;
          printf("INFO rank %i i_part %i xmin %f ymin %f xmax %f ymax %f nbPart %i test %i nVertex %i \n",
          rank,i_part,xminPart,yminPart,xmaxPart,ymaxPart,nbPart[i_code],nVertexSeg/nbPartSeg[i_code],nVertex[i_code][i_part]);

          grid_mesh(xminPart,
                      xmaxPart,
                      yminPart,
                      ymaxPart,
                      randLevel,
                      nVertexSegPart,
                      sqrt(localCommSize[i_code]),
                      coords           [i_code][i_part],
                      eltsConnecPointer[i_code][i_part],
                      eltsConnec       [i_code][i_part],
                      localComm        [i_code]);
        } //if code1
      }//loop on v
    }//loop on u
  } // loop on codes


  fprintf(outputFile, "   Number of vertex   : %i\n", nVertex);
  fprintf(outputFile, "   Number of elements : %i\n", nElts);

  for(int i_code = 0; i_code < n_code_name; i_code++) {
    for(int i_part=0;i_part<nbPart[i_code];i_part++)  {
      CWP_Mesh_interf_vtx_set(codeName[i_code],             //Code name
                            "c_new_api_surf_cpl_P1P0_P0P1_part",  // Coupling id
                            i_part,
                            nVertex[i_code][i_part],
                            coords[i_code][i_part],
                            NULL);
    }
  }

  for(int i_code = 0; i_code < n_code_name; i_code++) {

    int block_id = CWP_Mesh_interf_block_add(codeName[i_code],             // Code name
                                           "c_new_api_surf_cpl_P1P0_P0P1_part",  // Coupling id
                                           CWP_BLOCK_FACE_POLY);

    for(int i_part=0;i_part<nbPart[i_code];i_part++)  {
      CWP_Mesh_interf_f_poly_block_set(codeName[i_code],             //Code name
                                     "c_new_api_surf_cpl_P1P0_P0P1_part",  // Coupling id
                                     i_part,
                                     block_id,
                                     nElts[i_code][i_part],
                                     eltsConnecPointer[i_code][i_part],
                                     eltsConnec       [i_code][i_part],
                                     NULL);

     }
  }

  for(int i_code = 0; i_code < n_code_name; i_code++) {
     CWP_Mesh_interf_finalize (codeName[i_code],
                              "c_new_api_surf_cpl_P1P0_P0P1_part"  // Coupling id
                            );
  }

  /* Fields exchange
   *     - Proc 0 : Send X coordinates
   *                Recv Y coordinates
   *     - Proc 1 : Send Y coordinates
   *                Recv X coordinates
   * --------------------------------- */

 // if (rank == 0)
    printf("        Exchange Code1 <-> Code2 %i\n",rank);

  double ***sendValues = (double ***) malloc(sizeof(double*) * n_code_name);
  double ***sendValues2 = (double ***) malloc(sizeof(double*) * n_code_name);
  double ***sendValues3 = (double ***) malloc(sizeof(double*) * n_code_name);
  double ***recvValues = (double ***) malloc(sizeof(double*) * n_code_name);
  double ***recvValues2 = (double ***) malloc(sizeof(double*) * n_code_name);
  double ***recvValues3 = (double ***) malloc(sizeof(double*) * n_code_name);
  double ***Values2Vertex = (double ***) malloc(sizeof(double*) * n_code_name);
  double ***recvValuesUser = (double ***) malloc(sizeof(double*) * n_code_name);

  int nbPoints = 3;
  int** nbPointsUser = (int **) malloc(sizeof(int*) * n_code_name);

  double*** coordsPointsUser = (int ***) malloc(sizeof(int**) * n_code_name);

  for(int i_code = 0; i_code < n_code_name; i_code++) {
    sendValues[i_code] = (double **) malloc(sizeof(double*) * nbPart[i_code]);
    sendValues2[i_code] = (double **) malloc(sizeof(double*) * nbPart[i_code]);
    sendValues3[i_code] = (double **) malloc(sizeof(double*) * nbPart[i_code]);

    nbPointsUser[i_code] = (int *) malloc(sizeof(int) * nbPart[i_code]);
    recvValues[i_code] = (double **) malloc(sizeof(double*) * nbPart[i_code]);
    recvValues2[i_code] = (double **) malloc(sizeof(double*) * nbPart[i_code]);
    recvValues3[i_code] = (double **) malloc(sizeof(double*) * nbPart[i_code]);
    Values2Vertex[i_code] = (double **) malloc(sizeof(double*) * nbPart[i_code]);
    recvValuesUser[i_code] = (double **) malloc(sizeof(double*) * nbPart[i_code]);
    coordsPointsUser[i_code] = (double **) malloc(sizeof(double*) * nbPart[i_code]);
    for(int i_part=0;i_part<nbPart[i_code];i_part++)  {

      if (codeName[i_code] == "code1") {
        sendValues[i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);
        sendValues2[i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);
        sendValues3[i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);

        recvValues[i_code][i_part] = (double *) malloc(sizeof(double) * nElts[i_code][i_part]);
        Values2Vertex[i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);
        nbPointsUser[i_code][i_part] = 0;
        recvValuesUser[i_code][i_part] = (double *) malloc(sizeof(double) * nbPointsUser[i_code][i_part]);
        coordsPointsUser[i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nbPointsUser[i_code][i_part]  );
        for (int i = 0; i <nVertex[i_code][i_part]; i++) {
          sendValues[i_code][i_part][i] = coords[i_code][i_part][3 * i];
        }
      }
      else {
        sendValues[i_code][i_part] = (double *) malloc(sizeof(double) * 3*nElts[i_code][i_part]     );
        sendValues2[i_code][i_part] = (double *) malloc(sizeof(double) * 3*nElts[i_code][i_part]     );
        sendValues3[i_code][i_part] = (double *) malloc(sizeof(double) * 3*nElts[i_code][i_part]     );

        recvValues[i_code][i_part] = (double *) malloc(sizeof(double) * 3* nVertex[i_code][i_part]);
        recvValues2[i_code][i_part] = (double *) malloc(sizeof(double) * 3* nVertex[i_code][i_part]);
        recvValues3[i_code][i_part] = (double *) malloc(sizeof(double) * 3* nVertex[i_code][i_part]);
        Values2Vertex[i_code][i_part] = (double *) malloc(sizeof(double) * nVertex[i_code][i_part]);
        nbPointsUser[i_code][i_part] = nbPoints;
        recvValuesUser[i_code][i_part] = (double *) malloc(sizeof(double) * nbPointsUser[i_code][i_part]);
        coordsPointsUser[i_code][i_part] = (double *) malloc(sizeof(double) * 3 * nbPointsUser[i_code][i_part]  );
        for (int i = 0; i <nElts[i_code][i_part]; i++) {
          sendValues[i_code][i_part][3*i] = rank;
          sendValues2[i_code][i_part][3*i] = i_part;
          sendValues3[i_code][i_part][3*i] = i;
        }
        for (int i = 0; i <nVertex[i_code][i_part]; i++) {
          Values2Vertex[i_code][i_part][i] = i;//coords[3 * i];
        }


        for (int i = 0; i <nbPointsUser[i_code][i_part]; i++) {
          coordsPointsUser[i_code][i_part][3*i  ]= 0.0+0.01*i_part+0.001*i;
          coordsPointsUser[i_code][i_part][3*i+1]= rank*0.1+0.01*i_part+0.001*i;
          coordsPointsUser[i_code][i_part][3*i+2]= 0.0;
        }
      }
    }//loop on part
  }
  /* Define fields to send (X coordinate or Y coordinate) */

  /* Exchange */

  int code1I = 1;
  int code2I = 0;
  int code3I = 0;
  int code4I = 0;
  int code5I = 0;
  int code6I = 1;
  int code7I = 1;


  int nNotLocatedPoints = 0;
  char *fieldName1;
  char *fieldName2;
  char *fieldName3;
  fieldName1 = "cooX";
  fieldName2 = "rank";
  fieldName3 = "userField";
  char *fieldName4 = "userInterpolation";
  char *fieldName5 = "sendRecvField";
  char *fieldName6 = "part";
  char *fieldName7 = "num";

  CWP_Status_t visu_status = CWP_STATUS_OFF;

  printf("        Field %i\n",rank);
  for(int i_code = 0; i_code < n_code_name; i_code++) {
    if (codeName[i_code] == "code1") {

     if(code1I==1)  CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       CWP_STATUS_ON);


     if(code6I==1)  CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName6,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       CWP_STATUS_ON);

     if(code7I==1)  CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName7,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       CWP_STATUS_ON);

     if(code2I==1) CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       CWP_STATUS_ON);
     if(code3I==1)
     CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName3,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       CWP_STATUS_OFF);


    if(code4I==1)
     CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName4,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       CWP_STATUS_ON);

    if(code5I==1)
     CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName5,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       CWP_STATUS_OFF);


     for(int i_part=0;i_part<nbPart[i_code];i_part++)  {
       if(code1I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName1,i_part, sendValues[i_code][i_part]);
       if(code6I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName6,i_part, sendValues2[i_code][i_part]);
       if(code7I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName7,i_part, sendValues3[i_code][i_part]);
       if(code2I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName2,i_part, recvValues[i_code][i_part]);
       if(code3I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName3,i_part, sendValues[i_code][i_part]);
       if(code4I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName4,i_part, sendValues[i_code][i_part]);
       if(code5I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName5,i_part, Values2Vertex[i_code][i_part]);
     }

     if(code4I==1) CWP_Interp_from_location_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName4,_userInterpolation);

    }
    else {
       if(code1I==1)
       CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       CWP_STATUS_ON);


       if(code6I==1)
       CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName6,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       CWP_STATUS_ON);

       if(code7I==1)
       CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName7,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       CWP_STATUS_ON);


       if(code2I==1)
       CWP_Field_create (codeName[i_code],
                     "c_new_api_surf_cpl_P1P0_P0P1_part",
                     fieldName2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_BLOCK,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     CWP_STATUS_ON);

       if(code3I==1)
       CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName3,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_USER,
                       CWP_FIELD_EXCH_RECV,
                       CWP_STATUS_OFF);

      if(code4I==1)
      CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName4,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       CWP_STATUS_ON);

      if(code5I==1)
      CWP_Field_create (codeName[i_code],
                       "c_new_api_surf_cpl_P1P0_P0P1_part",
                       fieldName5,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       CWP_STATUS_ON);

     for(int i_part=0;i_part<nbPart[i_code];i_part++)  {
       if(code1I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName1,i_part, recvValues[i_code][i_part]);
       if(code6I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName6,i_part, recvValues2[i_code][i_part]);
       if(code7I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName7,i_part, recvValues3[i_code][i_part]);
       if(code2I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName2,i_part, sendValues[i_code][i_part]);
       if(code3I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName3,i_part, recvValuesUser[i_code][i_part]);
       if(code4I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName4,i_part, recvValues[i_code][i_part]);
       if(code5I==1) CWP_Field_data_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName5,i_part, Values2Vertex[i_code][i_part]);

       if(code3I==1)
       CWP_User_tgt_pts_set(codeName[i_code],
                            "c_new_api_surf_cpl_P1P0_P0P1_part",
                            i_part,
                            nbPointsUser[i_code][i_part],
                            coordsPointsUser[i_code][i_part]);


     }//loop on part
    }// if

  }

  MPI_Barrier(MPI_COMM_WORLD);
  printf("Before compute\n");

  for(int i_code = 0; i_code < n_code_name; i_code++) {
    CWP_Spatial_interp_weights_compute(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part");
  }

  double recv_time = 0.150;

  printf("After compute\n");
  MPI_Barrier(MPI_COMM_WORLD);


  for(int i_code = 0; i_code < n_code_name; i_code++) {

    CWP_Next_recv_time_set(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",recv_time);

    if (codeName[i_code] == "code1") {
      if(code1I==1) CWP_Field_Issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName1);
      if(code1I==1) CWP_Wait_issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName1);

      if(code6I==1) CWP_Field_Issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName6);
      if(code6I==1) CWP_Wait_issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName6);

      if(code7I==1) CWP_Field_Issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName7);
      if(code7I==1) CWP_Wait_issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName7);

    }
    else {
     if(code1I==1) CWP_Field_Irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName1);
     if(code1I==1) CWP_Wait_irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName1);

     if(code6I==1) CWP_Field_Irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part"     ,fieldName6);
     if(code6I==1) CWP_Wait_irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName6);

     if(code7I==1) CWP_Field_Irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part"     ,fieldName7);
     if(code7I==1) CWP_Wait_irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName7);

    }
  }
   MPI_Barrier(MPI_COMM_WORLD);

  for(int i_code = 0; i_code < n_code_name; i_code++) {

    if (codeName[i_code] == "code1") {
      if(code2I==1) CWP_Field_Irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName2);
      if(code2I==1) CWP_Wait_irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName2);

    }
    else {
     if(code2I==1) CWP_Field_Issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName2);
     if(code2I==1) CWP_Wait_issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName2);
    }
  }

   MPI_Barrier(MPI_COMM_WORLD);

  for(int i_code = 0; i_code < n_code_name; i_code++) {

    if (codeName[i_code] == "code1") {

      if(code3I==1) CWP_Field_Issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName3);
      if(code3I==1) CWP_Wait_issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName3);

    }
    else {
     if(code3I==1) CWP_Field_Irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName3);
     if(code3I==1) CWP_Wait_irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName3);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for(int i_code = 0; i_code < n_code_name; i_code++) {
    if (codeName[i_code] == "code1") {
      if(code4I==1) CWP_Field_Issend      (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName4);
      if(code4I==1) CWP_Wait_issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName4);
    }
    else {
     if(code4I==1) CWP_Field_Irecv       (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName4);
     if(code4I==1) CWP_Wait_irecv  (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName4);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for(int i_code = 0; i_code < n_code_name; i_code++) {
    if (codeName[i_code] == "code1") {
      if(code5I==1) CWP_Field_Irecv      (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName5);
      if(code5I==1) CWP_Wait_irecv (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName5);
    }
    else {
     if(code5I==1) CWP_Field_Issend      (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName5);
     if(code5I==1) CWP_Wait_issend (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName5);
    }
  }


 /* for(int i_code = 0; i_code < n_code_name; i_code++) {
    if(code5I==1) CWP_Sendrecv (codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part",fieldName5);
  }
*/

  /* Coupling deletion
   * ----------------- */

  //if (rank == 0)
   printf("        Delete mesh\n",rank);
   MPI_Barrier(MPI_COMM_WORLD);
  for(int i_code = 0; i_code < n_code_name; i_code++)
     CWP_Mesh_interf_del(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part");

  if (rank == 0)
    printf("        Delete coupling\n");

  for(int i_code = 0; i_code < n_code_name; i_code++)
    CWP_Cpl_del(codeName[i_code],"c_new_api_surf_cpl_P1P0_P0P1_part");

  /* Freeing memory
   * -------------- */
  for(int i_code = 0; i_code < n_code_name; i_code++) {
    for(int i_part=0;i_part<nbPart[i_code];i_part++)  {
      free(coords[i_code][i_part]);
      free(sendValues[i_code][i_part]);
      free(recvValues[i_code][i_part]);
      free(Values2Vertex[i_code][i_part]);
    }
    free(coords    [i_code]);
    free(sendValues[i_code]);
    free(recvValues[i_code]);
    free(Values2Vertex[i_code]);
  }

  free(coords);
  free(sendValues);
  free(Values2Vertex);
  free(recvValues);
  free(srcName);

  /* Finalize
   * -------- */

  CWP_Finalize();

  MPI_Finalize();

  fclose(outputFile);

  return EXIT_SUCCESS;
}
