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
#include <assert.h>

#include "cwp.h"
#include "cwipi.h"
#include "cwipi_config.h"
#include "pdm_part.h"
#include "pdm_timer.h"
#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mesh_nodal.h"
#include "pdm_poly_surf_gen.h"
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
     "  -nx     <val>    Global number of vertices in the side of domaine.\n\n"
     "  -part   <val>    Part of active ranks in the coupling.\n\n"
     "  -s      <val>    Size of domain.\n\n"
     "  -h               this message.\n\n");

  exit(exit_code);
}

/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nx             --> Global number of vertices in the side of domaine
 *   part           --> Part of active ranks in the coupling
 *   s              --> Size of domain
 *
 *---------------------------------------------------------------------*/

static void
_read_args
(
 int            argc,
 char         **argv,
 CWP_g_num_t   *nx,
 double        *part,
 double        *s,
  int           *n_compute
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc+1)
        _usage(EXIT_FAILURE);
      else
        *nx = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *part = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *s = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-nc") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_compute = atoi(argv[i]);
    } 
    i++;
  }
}



int _randomGlobalInt(){
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int resultInt = -1;
  int* randInt = NULL;
  if(rank == 0){
    srand(time(NULL));
    randInt = (int*)malloc(sizeof(int)*size);
    for(int i=0;i<size;i++){
      int test = -1;
      int ok = 0;
      do
      {
        test = rand()%size;
        ok = 1;
        for(int j=0;j<i;j++){
          if(test == randInt[j] ){
            ok = 0;
            break;
          }
        }
      }
      while(!ok);
      randInt[i] = test;
    }
  }
  MPI_Scatter(randInt,1,MPI_INT,&resultInt,1,MPI_INT,0,MPI_COMM_WORLD);
  return resultInt;
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

  /* Read args from command line
   * --------------------------- */

  CWP_g_num_t nx = 10;
  double part = 1.;
  double s = 10.;
  int n_compute = 1;
  const double dev_limit = 0.05;

  _read_args (argc, argv, &nx, &part, &s, &n_compute);
  

  /* Coupled MPI process fraction for each code */
  double* cpl_frac = (double*) malloc(sizeof(double)*2);
  cpl_frac[0] = part;
  cpl_frac[1] = part;


  /* Non null interface mesh fraction for each code */
  double* non_null_mesh_frac = (double*) malloc(sizeof(double)*2);
  non_null_mesh_frac[0] = 1.0;
  non_null_mesh_frac[1] = 1.0;


  /* Init + create coupling
   * ---------------------- */


  assert (commWorldSize >=2);

  char* fileName = (char *) malloc(sizeof(char) * 44);
  sprintf(fileName,"c_new_surf_coupling_poly_P1P0_P0P1_%4.4d.txt",rank);

 // outputFile = fopen(fileName,"w");

  free(fileName);

  // CWP_Output_file_set(outputFile);

  int           n_code;
  double prop1 = 1./2.;
  double prop2 = 1./2.;
  double prop12 = 0./3.;

  int partial_covering_mpi_domains = 1;
  if( partial_covering_mpi_domains == 1){
    prop1 = 1./3.;
    prop2 = 1./3.;
    prop12 = 1./3.;
  }

  int randomGlobalInt = _randomGlobalInt();
  
  //printf("rankBis %i %i\n",rank,randomGlobalInt);
  
  if (randomGlobalInt < (double)commWorldSize * prop1) {
    n_code = 1;
  }
  else if (randomGlobalInt < (double)commWorldSize * (prop1+prop12)) {
    n_code = 2;
  }
  else{
    n_code = 1;
  }
  
  

  CWP_g_num_t* nxCode = (CWP_g_num_t*)malloc( n_code*sizeof(CWP_g_num_t)); 

  char** codeName = (char**)malloc(sizeof(char*)*n_code); 
  char** codeCoupledName = (char**)malloc(sizeof(char*)*n_code); 
  int* codeId = (int*)malloc(sizeof(int)*n_code); 
  MPI_Comm* localComm = malloc(sizeof(MPI_Comm)*n_code);
  int* localRank = (int*)malloc(sizeof(int)*n_code); 
  CWP_Status_t*  is_coupled_rank = (CWP_Status_t*)malloc(sizeof(CWP_Status_t)*n_code); 
  double*       time_init = (double*)malloc(sizeof(double)*n_code); 
  int* nb_part = (int*) malloc(sizeof(int)*n_code);

    
  if (randomGlobalInt < (double)commWorldSize * prop1) {
    codeName[0] = "code1";
    codeId[0] = 1;
    codeCoupledName[0] = "code2";
    
    if(randomGlobalInt< (double)commWorldSize * prop1 * cpl_frac[0] )
      is_coupled_rank[0] = CWP_STATUS_ON;
    else
      is_coupled_rank[0] = CWP_STATUS_OFF;
      
    time_init[0] = 0.0;
    nxCode[0] = nx;
    nb_part[0] = 2;
  }
  else if (randomGlobalInt < (double)commWorldSize * (prop1+prop12)) {
    codeName[0] = "code1";
    codeId[0] = 1;
    codeCoupledName[0] = "code2"; 

    if(randomGlobalInt< (double)commWorldSize * ( prop1 + cpl_frac[0]*prop12 ) )
      is_coupled_rank[0] = CWP_STATUS_ON;
    else
      is_coupled_rank[0] = CWP_STATUS_OFF;

    time_init[0] = 0.0;
    nxCode[0] = nx;
    nb_part[0] = 2;   
     
        
    codeName[1] = "code2";
    codeId[1] = 2;
    codeCoupledName[1] = "code1";
    if(randomGlobalInt< (double)commWorldSize * ( prop1 + cpl_frac[1]*prop12 )  )
      is_coupled_rank[1] = CWP_STATUS_ON;
    else
      is_coupled_rank[1] = CWP_STATUS_OFF;
      
    time_init[1] = 0.0;
    nxCode[1] = (CWP_g_num_t)(2.0*(double)nx);
    nb_part[1] = 4;
  }
  else if (randomGlobalInt < (double)commWorldSize * (prop1+prop12+prop2)) {
    codeName[0] = "code2";
    codeId[0] = 2;
    codeCoupledName[0] = "code1";
    
    if(randomGlobalInt< (double)commWorldSize * ( prop1+prop12 + cpl_frac[1]* prop2 )  )
      is_coupled_rank[0] = CWP_STATUS_ON;
    else
      is_coupled_rank[0] = CWP_STATUS_OFF;   
      
    time_init[0] = 0.0;  
    nxCode[0] =(CWP_g_num_t)(2.0*(double)nx); 
    nb_part[0] = 4;
  }
  

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char**)codeName,
           is_coupled_rank,
           time_init,
           localComm);

  for(int i_code=0; i_code<n_code; i_code++){
    printf("Rank %i: %i codes, codeName[%i] %s coupled %i prop1 %3.2f %3.2f\n",rank,n_code,i_code,codeName[i_code],is_coupled_rank[i_code],prop1,prop12);


    if(is_coupled_rank[i_code] == CWP_STATUS_ON){
      MPI_Comm_rank(localComm[i_code], &(localRank[i_code]));
      assert(localComm[i_code] != MPI_COMM_NULL);
    
      CWP_Cpl_create (codeName[i_code],                // Code name
                      "multipart_testcase",            // Coupling id
                      codeCoupledName[i_code],         // Coupled application id
                      CWP_COMM_PAR_WITH_PART,          // Coupling type
                      CWP_SPATIAL_INTERP_FROM_LOCATION,            // Solver type
                      nb_part[i_code],                 // Partition number
                      CWP_DISPLACEMENT_STATIC,         // Mesh type
                      CWP_FREQ_CPL_TIME_STEP);         // Postprocessing frequency

    }  
  }//end i_code loop

  if(rank==0)
    printf("After create_cpl\n");

  PDM_timer_t* timer = PDM_timer_create();
  PDM_timer_t* timer2 = PDM_timer_create();
  PDM_timer_init(timer);
  
  double time[3];
  
  int n_int = 1;
  double compute_time[n_compute];
  double compute_exch_time[n_int];

  double ***sendValues = (double***)malloc(sizeof(double**)*n_code);
  double ***recvValues = (double***)malloc(sizeof(double**)*n_code);
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  
  for(int i_code = 0;i_code<n_code;i_code++){
    if(is_coupled_rank[i_code] == CWP_STATUS_ON){
      CWP_Visu_set(codeName[i_code], // Code name
                   "multipart_testcase",     // Coupling id
                   1,           // Postprocessing frequency
                   Ensight,     // Postprocessing format
                   "text");     // Postprocessing option
      MPI_Comm connectableComm = CWP_Connectable_comm_get( codeName[i_code] );
      
      CWP_surf_gen_init(codeName[i_code], nxCode[i_code], nxCode[i_code], nb_part[i_code], &connectableComm, non_null_mesh_frac[i_code], s, (double)codeId[i_code]);
    }//end if(is_coupled_rank[i_code]
  }
  
  if(rank==0)
    printf("After gen_init\n");

  for(int i_code = 0;i_code<n_code;i_code++)
    if(is_coupled_rank[i_code] == CWP_STATUS_ON)    
      CWP_surf_gen_compute(codeName[i_code]);
  
  if(rank==0)
    printf("After gen_compute\n");

  int **eltsConnecPolyIndex = NULL;
  for(int i_code = 0;i_code<n_code;i_code++){
     if(is_coupled_rank[i_code] == CWP_STATUS_ON){    

        sendValues[i_code] = (double**)malloc(sizeof(double*)*nb_part[i_code]);
        recvValues[i_code] = (double**)malloc(sizeof(double*)*nb_part[i_code]);
      
        double **coords = (double**)malloc(sizeof(double*)*nb_part[i_code]);
      
        CWP_g_num_t **vtxGnum = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*)*nb_part[i_code]);
        int *nVtx = (int*)malloc(sizeof(int)*nb_part[i_code]);
    
        int *n_tri = (int*)malloc(sizeof(int)*nb_part[i_code]);
        int **eltsConnecTri = (int**)malloc(sizeof(int*)*nb_part[i_code]);
        CWP_g_num_t **eltsGnumTri = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*)*nb_part[i_code]);  
      
        int *n_quad = (int*)malloc(sizeof(int)*nb_part[i_code]);
        int **eltsConnecQuad = (int**)malloc(sizeof(int*)*nb_part[i_code]);  
        CWP_g_num_t **eltsGnumQuad = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*)*nb_part[i_code]);   
    
        int *n_poly2d = (int*)malloc(sizeof(int)*nb_part[i_code]);
        eltsConnecPolyIndex = (int**)malloc(sizeof(int*)*nb_part[i_code]);  
        int **eltsConnecPoly = (int**)malloc(sizeof(int*)*nb_part[i_code]);    


        int *n_faces = (int*)malloc(sizeof(int)*nb_part[i_code]);
        int **faceEdgeIdx = (int**)malloc(sizeof(int*)*nb_part[i_code]);  
        int **faceEdge = (int**)malloc(sizeof(int*)*nb_part[i_code]);    

        int *n_edges = (int*)malloc(sizeof(int)*nb_part[i_code]);
        int **edgeVtxIdx = (int**)malloc(sizeof(int*)*nb_part[i_code]);  
        int **edgeVtx = (int**)malloc(sizeof(int*)*nb_part[i_code]);    

        CWP_g_num_t **faceLNToGN = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*)*nb_part[i_code]);  
      

        CWP_g_num_t **eltsGnumPoly = (CWP_g_num_t**)malloc(sizeof(CWP_g_num_t*)*nb_part[i_code]);  
      
        int gnum_compute = 0;
        int face_edge = 1;

        for(int i_part=0;i_part<nb_part[i_code];i_part++){
    
          nVtx[i_part] = 0;
          coords[i_part] = NULL;
          vtxGnum[i_part] = NULL;
    
          int nElts = 0;
    
          n_tri[i_part] = 0;  
          eltsConnecTri[i_part] = NULL;
          eltsGnumTri[i_part] = NULL;
    
          n_quad[i_part] = 0;  
          eltsConnecQuad[i_part] = NULL;
          eltsGnumQuad[i_part] = NULL;
    
          n_poly2d[i_part] = 0;  
          eltsConnecPolyIndex[i_part] = NULL;
          eltsConnecPoly[i_part] = NULL;
          eltsGnumPoly[i_part] = NULL;  
    
          CWP_surf_gen_by_block_get( codeName[i_code], i_part,
                                     &nVtx[i_part]  , &coords[i_part] , &vtxGnum[i_part], &nElts,
                                     &n_tri[i_part] , &eltsConnecTri[i_part] , &eltsGnumTri[i_part],
                                     &n_quad[i_part], &eltsConnecQuad[i_part], &eltsGnumQuad[i_part],
                                     &n_poly2d[i_part], &eltsConnecPolyIndex[i_part], &eltsConnecPoly[i_part], &eltsGnumPoly[i_part]);

          if(face_edge == 1) {
             CWP_surf_face_edge_get( codeName[i_code], i_part,
                                     &nVtx[i_part]  , &coords[i_part] , &vtxGnum[i_part], 
                                     &n_faces[i_part],
                                     &faceEdgeIdx[i_part], &faceEdge[i_part], 
                                     &n_edges[i_part], &edgeVtxIdx[i_part], &edgeVtx[i_part], 
                                     & faceLNToGN[i_part]
                                   );

          }

          if(gnum_compute==0){
            vtxGnum     [i_part] = NULL;
            eltsGnumTri [i_part] = NULL;
            eltsGnumQuad[i_part] = NULL;
            eltsGnumPoly[i_part] = NULL;
          }


          CWP_Mesh_interf_vtx_set (codeName[i_code],             //Code name
                                   "multipart_testcase",             // Coupling id
                                   i_part,
                                   nVtx[i_part],
                                   coords[i_part],
                                   vtxGnum[i_part]);


          if(face_edge == 1) {

            CWP_Mesh_interf_from_faceedge_set (codeName[i_code],
                                               "multipart_testcase",  // Coupling id
                                               i_part,
                                               n_faces[i_part],
                                               faceEdgeIdx[i_part],faceEdge[i_part],
                                               n_edges[i_part],
                                               edgeVtxIdx[i_part], edgeVtx[i_part],
                                               faceLNToGN[i_part]
                                             );
          }
          else{

             int block_id_tri = CWP_Mesh_interf_block_add (codeName[i_code],             // Code name
                                                          "multipart_testcase",             // Coupling id
                                                          CWP_BLOCK_FACE_TRIA3);
    
            CWP_Mesh_interf_block_std_set (codeName[i_code],
                                         "multipart_testcase",  // Coupling id
                                         i_part,
                                         block_id_tri,
                                         n_tri[i_part],
                                         eltsConnecTri[i_part],
                                         eltsGnumTri[i_part]);

            int block_id_quad = CWP_Mesh_interf_block_add (codeName[i_code],             // Code name
                                                           "multipart_testcase",             // Coupling id
                                                           CWP_BLOCK_FACE_QUAD4);
      
            CWP_Mesh_interf_block_std_set (codeName[i_code],
                                         "multipart_testcase",  // Coupling id
                                         i_part,
                                         block_id_quad,
                                         n_quad[i_part],
                                         eltsConnecQuad[i_part],
                                         eltsGnumQuad[i_part]
                                         );

            int block_id_poly = CWP_Mesh_interf_block_add (codeName[i_code],             // Code name
                                                           "multipart_testcase",             // Coupling id
                                                           CWP_BLOCK_FACE_POLY);
        
             CWP_Mesh_interf_f_poly_block_set (codeName[i_code],
                                            "multipart_testcase",  // Coupling id
                                            i_part              ,
                                            block_id_poly       ,
                                            n_poly2d[i_part]            ,
                                            eltsConnecPolyIndex[i_part] ,
                                            eltsConnecPoly[i_part]      ,
                                            eltsGnumPoly[i_part]);
          }
                        
          if (codeId[i_code] == 1 ) {
            sendValues[i_code][i_part] = (double *) malloc(sizeof(double) * nVtx[i_part]);
            recvValues[i_code][i_part] = (double *) malloc(sizeof(double) * nElts);
            for (int i = 0; i < nVtx[i_part]; i++) {
              sendValues[i_code][i_part][i] = coords[i_part][3 * i];
            }
          }
          else {
          
            double* field_tri   ;
            double* field_quad  ;
            double* field_poly2d;            
            CWP_surf_gen_tri_field_get( codeName[i_code], i_part,  &field_tri);
            CWP_surf_gen_quad_field_get( codeName[i_code], i_part,  &field_quad);
            CWP_surf_gen_poly_field_get( codeName[i_code], i_part,  &field_poly2d);         
            sendValues[i_code][i_part] = (double *) malloc(sizeof(double) * nElts);
            recvValues[i_code][i_part] = (double *) malloc(sizeof(double) * nVtx[i_part]);
            
            int ind=0;

            for (int i = 0; i < n_tri[i_part]; i++) {
              sendValues[i_code][i_part][ind] = field_tri[i];
              ind++;
            }

            for (int i = 0; i < n_quad[i_part]; i++) {
              sendValues[i_code][i_part][ind] = field_quad[i];
              ind++;
            }   

            for (int i = 0; i < n_poly2d[i_part]; i++) {
              sendValues[i_code][i_part][ind] = field_poly2d[i];
              ind++;
            }    

          }                                       
                                           
        }//end i_part loop
     }//end if(is_coupled      
  }//end i_code loop    

  if(rank==0)
    printf("After Mesh init\n");
  
  for(int i_code = 0;i_code<n_code;i_code++){      
      
      if(is_coupled_rank[i_code] == CWP_STATUS_ON){    
        CWP_Mesh_interf_finalize (codeName[i_code],
                                  "multipart_testcase");  // Coupling id
      }
  }//end i_code loop

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("After interf_finalize\n");

  /*************************************
   * Fields exchange
   *     - Proc 0 : Send X coordinates
   *                Recv Y coordinates
   *     - Proc 1 : Send Y coordinates
   *                Recv X coordinates
   * --------------------------------- */

  if(rank==0)
    printf("        Exchange Code1 <-> Code2 %i\n",rank);

  /* Exchange */

  CWP_Status_t visu_status = CWP_STATUS_ON;
  char *fieldName1 = "cooX";
  char *fieldName2 = "rank";

  for(int i_code = 0;i_code<n_code;i_code++){  
    if (codeName[i_code] == "code1" && is_coupled_rank[i_code] == CWP_STATUS_ON) {
  /*    CWP_Field_create (codeName[i_code],
                       "multipart_testcase",
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_FIELD_VALUE_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
*/
      CWP_Field_create (codeName[i_code],
                        "multipart_testcase",
                       fieldName2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_FIELD_VALUE_CELL_POINT,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      
        
      for(int i_part=0;i_part<nb_part[i_code];i_part++){
    /*      CWP_Field_data_set(codeName[i_code],
                             "multipart_testcase",
                             fieldName1,
                             i_part,
                             sendValues[i_code][i_part]);
*/
          CWP_Field_data_set(codeName[i_code],
                            "multipart_testcase",
                            fieldName2,
                            i_part,
                            recvValues[i_code][i_part]);
      }
    }
    else if (codeName[i_code] == "code2" && is_coupled_rank[i_code] == CWP_STATUS_ON) {

   /*   CWP_Field_create (codeName[i_code],
                       "multipart_testcase",
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_FIELD_VALUE_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
 */
      CWP_Field_create (codeName[i_code],
                        "multipart_testcase",
                        fieldName2,
                        CWP_DOUBLE,
                        CWP_FIELD_STORAGE_BLOCK,
                        1,
                        CWP_FIELD_VALUE_CELL_POINT,
                        CWP_FIELD_EXCH_SEND,
                        visu_status);


        for(int i_part=0;i_part<nb_part[i_code];i_part++){
      /*    CWP_Field_data_set(codeName[i_code],
                             "multipart_testcase",
                             fieldName1,
                             i_part,
                             recvValues[i_code][i_part]);
        */                     
          CWP_Field_data_set(codeName[i_code],
                             "multipart_testcase",
                             fieldName2,
                             i_part,
                             sendValues[i_code][i_part]);

        }

     }
     
    }//end i_code loop     
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
       printf("After Fields");
    PDM_timer_init(timer);

    PDM_timer_resume(timer);
    
    PDM_timer_init(timer);
    PDM_timer_init(timer2);  
    
    PDM_timer_resume(timer);
    double mean = 0.0;
    double mean2 = 0.0;    
    double std_dev=0.0;
    int n_it=0;

    for(int i=0;i<n_compute;i++){
      PDM_timer_init(timer2);
      PDM_timer_resume(timer2);
      for(int i_code = 0;i_code<n_code;i_code++){  
        if (is_coupled_rank[i_code] == CWP_STATUS_ON ) 
          CWP_Spatial_interp_weights_compute(codeName[i_code],"multipart_testcase");
      }

      PDM_timer_hang_on(timer2);
      compute_time[i] = PDM_timer_elapsed(timer2);

      mean +=  compute_time[i];
      MPI_Allreduce(&mean,&mean2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      mean2 = mean2/((double)(i+1)*(double)commWorldSize);

      std_dev = 0.0;
      for(int h=0;h<=i;h++)
        std_dev += pow((compute_time[h] - mean2)/mean2 ,2);
      std_dev = sqrt(std_dev)/(double)(i+1);
      
      double std2;
      MPI_Allreduce(&std_dev,&std2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      std_dev = std2/(double)commWorldSize;  
      
      
      n_it = i;
      if(i>3 && std_dev<dev_limit) 
        i=n_compute+1;
      
      if(localRank[0]==0 ) printf("Survey exchange %i %5.4e\n",i,std_dev);
    }
    
    PDM_timer_hang_on(timer);

    MPI_Barrier(MPI_COMM_WORLD);

    if(localRank[0]==0)
      printf("New localization time %5.4e codeName[0] %s deviation %5.4e nb_it %i \n",  mean2,codeName[0],std_dev,n_it);
    double recv_time = 0.150;

    
    PDM_timer_init(timer);
    PDM_timer_resume(timer);
    
    mean = 0.0;
    std_dev=0.0;    

    MPI_Barrier(MPI_COMM_WORLD);



    for(int i =0; i<n_int;i++){
      PDM_timer_init(timer2);
      PDM_timer_resume(timer2);
      
      for(int i_code = 0;i_code<n_code;i_code++){  
        if(is_coupled_rank[i_code] == CWP_STATUS_ON){
          CWP_next_recv_time_set(codeName[i_code],"multipart_testcase",recv_time);
          recv_time *= 0.01;
      
          if (codeName[i_code] == "code1" ) {
           /*  CWP_Issend (codeName[i_code],"multipart_testcase",fieldName1);
             CWP_Wait_issend (codeName[i_code],"multipart_testcase",fieldName1);
           */  
             CWP_Irecv  (codeName[i_code],"multipart_testcase",fieldName2);
             CWP_Wait_irecv  (codeName[i_code],"multipart_testcase",fieldName2);
          }
          else if (codeName[i_code] == "code2") {
        /*    CWP_Irecv  (codeName[i_code],"multipart_testcase",fieldName1);
            CWP_Wait_irecv  (codeName[i_code],"multipart_testcase",fieldName1);      
          */  
            CWP_Issend (codeName[i_code],"multipart_testcase",fieldName2);
            CWP_Wait_issend (codeName[i_code],"multipart_testcase",fieldName2);
          }
        }
      }//end i_code loop
      PDM_timer_hang_on(timer2);
      compute_exch_time[i] = PDM_timer_elapsed(timer2);
      
      mean +=  compute_exch_time[i];
      MPI_Allreduce(&mean,&mean2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      mean2 = mean2/((double)(i+1)*(double)commWorldSize);

      std_dev = 0.0;
      for(int h=0;h<=i;h++)
        std_dev += pow((compute_exch_time[h] - mean2)/mean2 ,2);
      std_dev = sqrt(std_dev)/(double)(i+1);
      
      double std2;
      MPI_Allreduce(&std_dev,&std2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      std_dev = std2/(double)commWorldSize;        
      
      if(i>3 && std_dev<dev_limit) {
        i=n_int+1;
      }   
      if(localRank[0]==0) printf("Survey exchange %i %5.4e\n",i,std_dev);
      MPI_Barrier(MPI_COMM_WORLD);
    }  
  
  
    PDM_timer_hang_on(timer);
    
    if(localRank[0]==0)
      printf("New exchange time for %i iterations %5.4e s codeName[0] %s deviation %5.4e\n", n_int, mean2,codeName[0],std_dev); 


    MPI_Barrier(MPI_COMM_WORLD);

    for(int i_code = 0;i_code<n_code;i_code++){  
      if( is_coupled_rank[i_code] == CWP_STATUS_ON ) {
        CWP_Mesh_interf_del(codeName[i_code],"multipart_testcase");
        CWP_Cpl_del(codeName[i_code],"multipart_testcase");
        free(sendValues[i_code]);
        free(recvValues[i_code]);
      }
    }

  free(sendValues);
  free(recvValues);
  free(srcName);
  
  /* Finalize
   * -------- */

  CWP_Finalize();
   
  MPI_Finalize();
   
  //fclose(outputFile);
  return EXIT_SUCCESS;
}
