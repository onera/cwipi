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
#include "cwp_priv.h"
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
     "  -new             Use new CWIPI API.\n\n"
     "  -old             Use old CWIPI API.\n\n"
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
 *   new            --> New CWIPI API is used
 *   old            --> Old CWIPI API is used
 *
 *---------------------------------------------------------------------*/

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *nx,
 double        *part,
 double        *s,
 int           *new,
 int           *old,
  int           *n_compute
)
{
  int i = 1;
  int isNew = 0;
  int isOld = 0;

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
    else if (strcmp(argv[i], "-old") == 0) {
      if (isNew != 0) {
        printf("Error : new CWIPI is already selected\n");
        exit(1);
      }
      isOld = 1;
      if (i >= argc+1)
        _usage(EXIT_FAILURE);
      else {
        *old = 1;
        *new = 0;
      }
    }
    else if (strcmp(argv[i], "-nc") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_compute = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-new") == 0) {
      if (isOld != 0) {
        printf("Error : old CWIPI is already selected\n");
        exit(1);
      }
      isNew = 1;
      if (i >= argc+1)
        _usage(EXIT_FAILURE);
      else {
        *new = 1;
        *old = 0;
      }
    }
    i++;
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

  /* Read args from command line
   * --------------------------- */

  PDM_g_num_t nx = 10;
  double part = 1.;
  double s = 10.;
  int new = 0;
  int old = 1;
  const int nPart = 1;
  int n_compute = 10;
  const double dev_limit = 0.05;


  _read_args (argc, argv, &nx, &part, &s, &new, &old,&n_compute);


  /* Init + create coupling
   * ---------------------- */

  char *codeName;
  int codeId;
  char *codeCoupledName;

  assert (commWorldSize >=2);

  if (rank < commWorldSize / 2) {
    codeName = "code1";
    codeId = 1;
    codeCoupledName = "code2";
  }
  else {
    codeName = "code2";
    codeId = 2;
    codeCoupledName = "code1";
  }

  char* fileName = (char *) malloc(sizeof(char) * 44);
  sprintf(fileName,"c_new_surf_coupling_poly_P1P0_P0P1_%4.4d.txt",rank);

 // outputFile = fopen(fileName,"w");

  free(fileName);

  assert ((old == 1) || (new == 1));

  MPI_Comm localComm = MPI_COMM_NULL ;

  const int nb_part = 1;



  if (old) {
    //cwipi_set_output_listing(outputFile);
    cwipi_init(MPI_COMM_WORLD,
               codeName ,
               &localComm);

    cwipi_solver_type_t solver_type;

    if (codeId == 1) {
      solver_type = CWIPI_SOLVER_CELL_VERTEX;
    }
    else {
      solver_type = CWIPI_SOLVER_CELL_CENTER;
    }

    cwipi_create_coupling("old_cpl",                                 // Coupling id
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                          codeCoupledName,                           // Coupled application id
                          2,                                         // Geometric entities dimension
                          0.001,                                     // Geometric tolerance
                          CWIPI_STATIC_MESH,                         // Mesh type
                          solver_type,                               // Solver type
                          1,                                        // Postprocessing frequency
                          "EnSight Gold",                            // Postprocessing format
                          "text");                                   // Postprocessing option
  }
  else {
   // CWP_Output_file_set(outputFile);

    //const MPI_Comm      global_comm;
    const int           n_code = 1;
    const CWP_Status_t  is_coupled_rank = CWP_STATUS_ON;
    const double        time_init = 0.;

    CWP_Init(MPI_COMM_WORLD,
             1,
             (const char **) &(codeName),
             &is_coupled_rank,
             &time_init,
             &localComm);

    CWP_Cpl_create (codeName,                        // Code name
                    "new_cpl",                       // Coupling id
                    codeCoupledName,                 // Coupled application id
                    CWP_COMM_PAR_WITH_PART,          // Coupling type
                    CWP_SPATIAL_INTERP_FROM_LOCATION,               // Solver type
                    nb_part,                         // Partition number
                    CWP_DYNAMIC_MESH_STATIC,         // Mesh type
                    CWP_TIME_EXCH_CPL_TIME_STEP);         // Postprocessing frequency

    CWP_Visu_set(codeName, // Code name
                 "new_cpl",     // Coupling id
                 1,           // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT,     // Postprocessing format
                 "text");     // Postprocessing option
  }



  CWP_surf_gen_init("generator",nx, nx, nPart, &localComm, part, s, (double)codeId);
  CWP_surf_gen_compute("generator");


  int nVtx = 0;
  double * coords = NULL;
  CWP_g_num_t *vtxGnum = NULL;

  int nElts = 0;

  int localRank;
  MPI_Comm_rank(localComm, &localRank);


  /*
  double Rcurv = 10.0 * s ;
  for(int i =0; i<nVtx;i++){
     coords[3*i+2] = sqrt( pow(Rcurv,2) - pow(coords[3*i],2) - pow(coords[3*i+1],2) );
  */


  PDM_timer_t* timer = PDM_timer_create();
  PDM_timer_t* timer2 = PDM_timer_create();

  PDM_timer_init(timer);

  double time[3];

  int n_int = 1;
  double compute_time[n_compute];
  double compute_exch_time[n_int];

  if (old) {

    int *eltsConnecIndex = NULL;
    int *eltsConnec = NULL;
    CWP_g_num_t *eltsGnum = NULL;

    CWP_surf_gen_one_connectivity_get ("generator",0,
                                       &nVtx  , &coords , &vtxGnum,
                                       &nElts  , &eltsConnecIndex, &eltsConnec, &eltsGnum);

    cwipi_define_mesh("old_cpl",
                      nVtx,
                      nElts,
                      coords,
                      eltsConnecIndex,
                      eltsConnec);

    MPI_Barrier(MPI_COMM_WORLD);

    PDM_timer_init(timer);
    PDM_timer_init(timer2);

    double mean =  0.0;
    double std_dev=0.0;
    double mean2;

    int n_it=0;
    PDM_timer_resume(timer);
    for(int i=0;i<n_compute;i++){
      cwipi_update_location("old_cpl");
      PDM_timer_init(timer2);
      PDM_timer_resume(timer2);
      cwipi_locate("old_cpl");
      PDM_timer_hang_on(timer2);
      compute_time[i] = PDM_timer_elapsed(timer2);

      mean +=  compute_time[i];
      MPI_Allreduce(&mean,&mean2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      mean2 = mean2/((double)(i+1)*(double)commWorldSize);

      std_dev = 0.0;
      //printf("compute_time[i] %10.5e mmean %10.5e\n",compute_time[i],mean2);
      for(int h=0;h<=i;h++)
        std_dev += pow((compute_time[h] - mean2)/mean2 ,2);
      std_dev = sqrt(std_dev)/(double)(i+1);

      double std2;
      MPI_Allreduce(&std_dev,&std2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      std_dev = std2/(double)commWorldSize;
      n_it = i;
      if(i>3 && std_dev<dev_limit) {
        i=n_compute+1;
      }

      if(mean2*n_compute > 600)
         n_compute=n_compute/2;

      if(localRank==0 ) printf("Survey localization %i %5.4e\n",i,std_dev);
    }
    PDM_timer_hang_on(timer);
    if(localRank==0)
      printf("Old localization time %5.4e codeName %s deviation %5.4e nb_it %i\n",  mean2,codeName,std_dev,n_it);
  }

  else {

    int n_tri = 0;
    int n_quad = 0;
    int n_poly2d = 0;

    int *eltsConnecQuad = NULL;
    CWP_g_num_t *eltsGnumQuad = NULL;

    int *eltsConnecTri = NULL;
    CWP_g_num_t *eltsGnumTri = NULL;

    int *eltsConnecPolyIndex = NULL;
    int *eltsConnecPoly = NULL;
    CWP_g_num_t *eltsGnumPoly = NULL;


    CWP_surf_gen_by_block_get( "generator", 0,
                               &nVtx  , &coords , &vtxGnum, &nElts,
                               &n_tri , &eltsConnecTri , &eltsGnumTri,
                               &n_quad, &eltsConnecQuad, &eltsGnumQuad,
                               &n_poly2d, &eltsConnecPolyIndex, &eltsConnecPoly, &eltsGnumPoly);

    CWP_Mesh_interf_vtx_set (codeName,             //Code name
                             "new_cpl",             // Coupling id
                             0,
                             nVtx,
                             coords,
                             vtxGnum);

    int block_id = CWP_Mesh_interf_block_add (codeName,             // Code name
                                              "new_cpl",             // Coupling id
                                              CWP_BLOCK_FACE_TRIA3);

    CWP_Mesh_interf_block_std_set (codeName,
                                   "new_cpl",  // Coupling id
                                   0,
                                   block_id,
                                   n_tri,
                                   eltsConnecTri,
                                   eltsGnumTri);

    block_id = CWP_Mesh_interf_block_add (codeName,             // Code name
                                          "new_cpl",             // Coupling id
                                          CWP_BLOCK_FACE_QUAD4);

    CWP_Mesh_interf_block_std_set (codeName,
                                     "new_cpl",  // Coupling id
                                     0,
                                     block_id,
                                     n_quad,
                                     eltsConnecQuad,
                                     eltsGnumQuad
                                     );

    block_id = CWP_Mesh_interf_block_add (codeName,             // Code name
                                          "new_cpl",             // Coupling id
                                          CWP_BLOCK_FACE_POLY);

    CWP_Mesh_interf_f_poly_block_set (codeName,
                                      "new_cpl",  // Coupling id
                                      0,
                                      block_id,
                                      n_poly2d,
                                      eltsConnecPolyIndex,
                                      eltsConnecPoly       ,
                                      eltsGnumPoly);


    CWP_Mesh_interf_finalize (codeName,
                              "new_cpl");  // Coupling id
  }


  /*************************************
   * Fields exchange
   *     - Proc 0 : Send X coordinates
   *                Recv Y coordinates
   *     - Proc 1 : Send Y coordinates
   *                Recv X coordinates
   * --------------------------------- */

  printf("        Exchange Code1 <-> Code2 %i\n",rank);

  double *sendValues;
  double *recvValues;

  if (codeId == 1 ) {
    sendValues = (double *) malloc(sizeof(double) * nVtx);
    recvValues = (double *) malloc(sizeof(double) * nElts);
    for (int i = 0; i < nVtx; i++) {
      sendValues[i] = coords[3 * i];
    }
  }

  else {
    sendValues = (double *) malloc(sizeof(double) * nElts);
    recvValues = (double *) malloc(sizeof(double) * nVtx);
    for (int i = 0; i < nElts; i++) {
      sendValues[i] = rank;//i;//coords[3 * i];
    }
  }

  /* Exchange */

  char *fieldName1 = "cooX";
  char *fieldName2 = "rank";


  if (old) {

    int sRequest, rRequest;
    int tag = 1;

    double mean = 0.0;
    double std_dev=0.0;
    double mean2 = 0.0;

    PDM_timer_init(timer);
    PDM_timer_resume(timer);

    for(int i =0; i<n_int;i++){
      PDM_timer_init(timer2);
      PDM_timer_resume(timer2);

      if (codeId == 1) {
      cwipi_issend("old_cpl",
                   "ech1",
                   tag,
                   1,
                   1,
                   0.1,
                   fieldName1,
                   sendValues,
                   &sRequest);

      cwipi_wait_issend("old_cpl", sRequest);
/*
      cwipi_irecv("old_cpl",
                  "ech2",
                  tag,
                  1,
                  1,
                  0.1,
                  fieldName2,
                  recvValues,
                  &rRequest);
g

      cwipi_wait_irecv("old_cpl", rRequest);
  */  }

    else {

      cwipi_irecv("old_cpl",
                  "ech1",
                  tag,
                  1,
                  1,
                  0.1,
                  fieldName1,
                  recvValues,
                  &rRequest);

      cwipi_wait_irecv("old_cpl", rRequest);
/*
      cwipi_issend("old_cpl",
                   "ech2",
                   tag,
                   1,
                   1,
                   0.1,
                   fieldName2,
                   sendValues,
                   &sRequest);

      cwipi_wait_issend("old_cpl", sRequest);
 */   }

      PDM_timer_hang_on(timer2);
      compute_exch_time[i] = PDM_timer_elapsed(timer2);

      mean +=  compute_exch_time[i];
      MPI_Allreduce(&mean,&mean2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      mean2 = mean2/((double)(i+1)*(double)commWorldSize);

      std_dev = 0.0;
      for(int h=0;h<=i;h++)
        std_dev += pow((compute_exch_time[h] - mean2)/mean2 ,2);
      std_dev = sqrt(std_dev)/(double)(i+1);
      /*printf("compute_exch_time[i] %10.5e mmean %10.5e std_dev  %10.5e compute_time[i] - mean2 %10.5e\n",
              compute_exch_time[i],mean2,std_dev,compute_time[i] - mean2);
        */
      double std2;
      MPI_Allreduce(&std_dev,&std2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      std_dev = std2/(double)commWorldSize;

      if(i>3 && std_dev<dev_limit) {
        i=n_int+1;
      }

      if(mean2*n_compute > 600)
         n_compute=n_compute/2;


      if(localRank==0 ) printf("Survey localization %i %5.4e\n",i,std_dev);
    }

     PDM_timer_hang_on(timer);

     if(localRank==0)
       printf("Old exchange time for %i iterations %5.4e s codeName %s deviation %5.4e\n", 0, mean2,codeName,std_dev);

  }

  else {

    int nNotLocatedPoints = 0;

    CWP_Status_t visu_status = CWP_STATUS_ON;

    if (codeName == "code1") {
      CWP_Field_create (codeName,
                       "new_cpl",
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
/*
      CWP_Field_create (codeName,
                        "new_cpl",
                       fieldName2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
*/
      CWP_Field_data_set(codeName,
                         "new_cpl",
                         fieldName1,
                         0,
                         sendValues);
/*
      CWP_Field_data_set(codeName,
                         "new_cpl",
                         fieldName2,
                         0,
                         recvValues);
        */
    }
    else {

      CWP_Field_create (codeName,
                       "new_cpl",
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_BLOCK,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
/*
      CWP_Field_create (codeName,
                        "new_cpl",
                        fieldName2,
                        CWP_DOUBLE,
                        CWP_FIELD_STORAGE_BLOCK,
                        1,
                        CWP_DOF_LOCATION_CELL_CENTER,
                        CWP_FIELD_EXCH_SEND,
                        visu_status);

      CWP_Field_data_set(codeName,
                         "new_cpl",
                         fieldName2,
                         0,
                         sendValues);
*/
      CWP_Field_data_set(codeName,
                         "new_cpl",
                         fieldName1,
                         0,
                         recvValues);



    }

    MPI_Barrier(MPI_COMM_WORLD);

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
      CWP_Spatial_interp_weights_compute(codeName,"new_cpl");;
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

      if(localRank==0 ) printf("Survey exchange %i %5.4e\n",i,std_dev);
    }
    PDM_timer_hang_on(timer);

    if(localRank==0 && new==1)
      printf("New localization time %5.4e codeName %s deviation %5.4e nb_it %i \n",  mean2,codeName,std_dev,n_it);

    double recv_time = 0.150;


    PDM_timer_init(timer);
    PDM_timer_resume(timer);

    mean = 0.0;
    std_dev=0.0;

    MPI_Barrier(MPI_COMM_WORLD);

    for(int i =0; i<n_int;i++){
      PDM_timer_init(timer2);
      PDM_timer_resume(timer2);

      CWP_Next_recv_time_set(codeName,"new_cpl",recv_time);
      recv_time *= 0.01;
      if (codeName == "code1") {
        CWP_Field_Issend (codeName,"new_cpl",fieldName1);
        CWP_Field_wait_issend (codeName,"new_cpl",fieldName1);
 /*     CWP_Field_Irecv  (codeName,"new_cpl",fieldName2);
      CWP_Field_wait_irecv  (codeName,"new_cpl",fieldName2);
 */     }
      else {
        CWP_Field_Irecv  (codeName,"new_cpl",fieldName1);
        CWP_Field_wait_irecv  (codeName,"new_cpl",fieldName1);
  /*    CWP_Field_Issend (codeName,"new_cpl",fieldName2);
      CWP_Field_wait_issend (codeName,"new_cpl",fieldName2);
  */  }

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
      if(localRank==0) printf("Survey exchange %i %5.4e\n",i,std_dev);
      MPI_Barrier(MPI_COMM_WORLD);
    }


     PDM_timer_hang_on(timer);

     if(localRank==0)
       printf("New exchange time for %i iterations %5.4e s codeName %s deviation %5.4e\n", n_int, mean2,codeName,std_dev);

  }
  int nNotLocatedPoints = 0;

  CWP_Status_t visu_status = CWP_STATUS_OFF;

  MPI_Barrier(MPI_COMM_WORLD);

  if (old) {
    cwipi_delete_coupling("old_cpl");
  }

  else {
    CWP_Mesh_interf_del(codeName,"new_cpl");
    CWP_Cpl_del(codeName,"new_cpl");
  }

  free(sendValues);
  free(recvValues);

  free(srcName);

  /* Finalize
   * -------- */

  if (old) {
    cwipi_finalize();
  }
  else {
    CWP_Finalize();
  }

  MPI_Finalize();

  //fclose(outputFile);
  return EXIT_SUCCESS;
}
