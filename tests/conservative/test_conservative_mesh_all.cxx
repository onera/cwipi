/*
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


#include "cwipi.h"
#include "oldMesh.hxx"
#include "conservativeMesh.hxx"
#include "fvm_writer.h"



/*----------------------------------------------------------------------
 *                                                                     
 * Dump status exchange                                                
 *                                                                     
 * parameters:
 *   status              <-- Exchange status           
 *---------------------------------------------------------------------*/

static void _dumpStatus(FILE* outputFile, cwipi_exchange_status_t status)
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



int main(
         int argc,
         char *argv[]

         ){
  
    FILE *outputFile;

  MPI_Init(&argc, &argv);

  int rank;
  int commWorldSize;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

  char *srcName = (char *) malloc (sizeof(char) * (strlen(__FILE__) + 1));
  strcpy(srcName, __FILE__);
  char *srcBaseName = strrchr(srcName, '.');
  *srcBaseName = '\0';
  srcBaseName = strrchr(srcName, '/') + 1;

  if (rank == 0)
    printf("\nSTART: %s\n", srcBaseName);

  srand(rank+time(0));

  int n_partition = 0;
  while(2 * pow(n_partition, 2) < commWorldSize) n_partition++;

  int n2 = 2 * pow(n_partition, 2);

//   if (n2 != commWorldSize) {
//     if (rank == 0)
//       printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
//     MPI_Finalize();
//     return EXIT_SUCCESS;
//   }


    /* Read args from command line
   * --------------------------- */

  int nVertexSeg = 10;
  double randLevel = 0.1;

  _read_args(argc, argv, &nVertexSeg, &randLevel);

  /* Initialization
   * -------------- */

  char *codeName;
  int codeId;
  char *codeCoupledName;

  int third_size; 

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

  char* fileName = (char *) malloc(sizeof(char) * 50);
  sprintf(fileName,"c_conservative_mesh_all%4.4d.txt",rank);

  outputFile = fopen(fileName,"w");

  free(fileName);

  cwipi_set_output_listing(outputFile);

  MPI_Comm localComm;
  cwipi_init(MPI_COMM_WORLD,
             codeName ,
             &localComm);

  /* Output redirection
   * ------------------ */

  int currentRank;
  int localCommSize;

  MPI_Comm_rank(localComm, &currentRank);
  MPI_Comm_size(localComm, &localCommSize);

  fprintf(outputFile, "  Conservative Mesh Test\n");
  fprintf(outputFile, "\n");
  
  fprintf(outputFile, "\nDump after initialization\n");
  fprintf(outputFile, "-------------------------\n");
  cwipi_dump_application_properties();
  
  if (rank == 0)
    printf("        Create coupling\n");
  
  cwipi_solver_type_t solver_type;
  
  /**********************  Maillage test ***********************************/
  

  int nVertexSM = 100;
  
  double *coordsSM = (double *) malloc(nVertexSM * 3 * sizeof(double));

  for(int i = 0; i < 10; i++){
    for(int j = 0; j < 10 ; j++){
      coordsSM[3*(10*i+j)] = j;
      coordsSM[3*(10*i+j) + 1] = i;
      coordsSM[3*(10*i+j) + 2] = 0;
    }
  }

  int nEltsSM = 81;
  int *eltsConnecPointerSM = (int *) malloc((nEltsSM+1) * sizeof(int));
  int *eltsConnecSM = (int *) malloc( 4*nEltsSM* sizeof(int));

  for(int i = 0; i < nEltsSM + 1; i++)
    eltsConnecPointerSM[i] = 4*i;

  for(int i = 0 ; i < 9 ; i++){
    for(int j = 0; j < 9 ; j++){
      eltsConnecSM[4*(9*i + j)] = 10*i + j + 1;
      eltsConnecSM[4*(9*i + j) + 1] = 10*i + j + 2;
      eltsConnecSM[4*(9*i + j) + 2] = 10*(i + 1) + j + 2;
      eltsConnecSM[4*(9*i + j) + 3] = 10*(i + 1) + j + 1;
    }
  }


  MPI_Comm fvmComm = MPI_COMM_WORLD;

    
  cwipi::oldMesh* sourceMesh = new cwipi::oldMesh(MPI_COMM_WORLD,
                                             2 ,
                                             nVertexSM,
                                             nEltsSM,
                                             coordsSM,
                                             eltsConnecPointerSM,
                                             eltsConnecSM);

  int nVertexTM = 9;
  
  double *coordsTM = (double *) malloc(nVertexTM * 3 * sizeof(double));
  
  coordsTM[0] = -4.;
  coordsTM[1] = 7;
  coordsTM[2] = 0.;

  coordsTM[3] = 1.5;
  coordsTM[4] = 3.5;
  coordsTM[5] = 0.;

  coordsTM[6] = 6.;
  coordsTM[7] = -1.;
  coordsTM[8] = 0.;

  coordsTM[9]  = 2.;
  coordsTM[10] = 11.;
  coordsTM[11] = 0.;

  coordsTM[12] = 5.5;
  coordsTM[13] = 7.5;
  coordsTM[14] = 0;

  coordsTM[15] = 10.;
  coordsTM[16] = 3.;
  coordsTM[17] = 0.;

  coordsTM[18] = 6.;
  coordsTM[19] = 15.;
  coordsTM[20] = 0;

  coordsTM[21]  = 9.5;
  coordsTM[22] = 10.5;
  coordsTM[23] = 0;

  coordsTM[24] = 14.;
  coordsTM[25] = 7.;
  coordsTM[26] = 0;


  int nEltsTM = 4;
  int *eltsConnecPointerTM = (int *) malloc((nEltsTM+1) * sizeof(int));
  int *eltsConnecTM = (int *) malloc( 16* sizeof(int));

  eltsConnecPointerTM[0] = 0;
  eltsConnecPointerTM[1] = 4;
  eltsConnecPointerTM[2] = 8;
  eltsConnecPointerTM[3] = 12;
  eltsConnecPointerTM[4] = 16;


  eltsConnecTM[0] = 1;
  eltsConnecTM[1] = 2;
  eltsConnecTM[2] = 5;
  eltsConnecTM[3] = 4;

  eltsConnecTM[4] = 2;
  eltsConnecTM[5] = 3;
  eltsConnecTM[6] = 6;
  eltsConnecTM[7] = 5;

  eltsConnecTM[8] = 4;
  eltsConnecTM[9] = 5;
  eltsConnecTM[10] = 8;
  eltsConnecTM[11] = 7;

  eltsConnecTM[12] = 5;
  eltsConnecTM[13] = 6;
  eltsConnecTM[14] = 9;
  eltsConnecTM[15] = 8;

  cwipi::oldMesh* targetMesh = new cwipi::oldMesh( fvmComm,
                                             2 ,
                                             nVertexTM,
                                             nEltsTM,
                                             coordsTM,
                                             eltsConnecPointerTM,
                                             eltsConnecTM);

  cwipi::ConservativeMesh* intersMesh = new cwipi::ConservativeMesh(fvmComm,
                                                                    *sourceMesh,
                                                                    *targetMesh,
                                                                    0.01);


  cwipi::ConservativeMesh* intersMesh2 = new cwipi::ConservativeMesh(fvmComm,
                                                                    *targetMesh,
                                                                    *sourceMesh,
                                                                    0.01);
fvm::fvm_writer_t *fvmWriterIM = fvm::fvm_writer_init("MeshIMallInters",
                                                        "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                                        "Ensight Gold",
                                                        "text",
                                                        fvm::FVM_WRITER_FIXED_MESH);

  fvm::fvm_writer_export_nodal(fvmWriterIM, 
                               &(intersMesh->getIntersectionMesh()->getFvmNodal()));

  fvm::fvm_writer_finalize(fvmWriterIM);


  fvmWriterIM = fvm::fvm_writer_init("MeshIMallInters2",
                                     "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                     "Ensight Gold",
                                     "text",
                                     fvm::FVM_WRITER_FIXED_MESH);

  fvm::fvm_writer_export_nodal(fvmWriterIM, 
                               &(intersMesh2->getIntersectionMesh()->getFvmNodal()));

  fvm::fvm_writer_finalize(fvmWriterIM);


  free(srcName);
  free(eltsConnecPointerSM);
  free(eltsConnecPointerTM);
  free(eltsConnecSM);
  free(eltsConnecTM);
  free(coordsTM);
  free(coordsSM);
  delete sourceMesh;
  delete targetMesh;
  delete intersMesh2;
  delete intersMesh;

  return EXIT_SUCCESS;
}
