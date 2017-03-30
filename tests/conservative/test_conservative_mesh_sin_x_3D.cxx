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
#include "fvm_writer.h"


#include "cwipi.h"
#include "mesh.hxx"
#include "conservativeMesh.hxx"



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
  sprintf(fileName,"c_conservative_mesh_sin_x3D%4.4d.txt",rank);

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
    printf("Create coupling\n");
  
  cwipi_solver_type_t solver_type;

  /**************** Maillage intersection sin(x) et x  **************/

  int nSM = 30;
  int mSM = 5;
  int nVertexSM = nSM*(2*mSM - 1);

  double *coordsSM = (double *) malloc(nVertexSM * 3 * sizeof(double));

  for(int i = 0; i < nSM; i++){
    for(int j = -mSM + 1; j < mSM ; j++){
      coordsSM[3*((2*mSM - 1)*i+j + mSM - 1)] = 1.*j/mSM;
      coordsSM[3*((2*mSM - 1)*i+j + mSM - 1) + 1] = 1.*i/nSM;
      coordsSM[3*((2*mSM - 1)*i+j + mSM - 1) + 2] = 1.*sin(1.*j/(mSM));
    }
  }

  int nEltsSM = (nSM-1)*(2*mSM-2);
  int *eltsConnecPointerSM = (int *) malloc((nEltsSM+1) * sizeof(int));
  int *eltsConnecSM = (int *) malloc( 4*nEltsSM* sizeof(int));

  for(int i = 0; i < nEltsSM + 1; i++)
    eltsConnecPointerSM[i] = 4*i;

  for(int i = 0 ; i < nSM - 1 ; i++){
    for(int j = - mSM + 1; j < mSM-1 ; j++){
      eltsConnecSM[4*((2*mSM-2)*i + j + mSM - 1)] = (2*mSM - 1)*i + j + mSM - 1 + 1;
      eltsConnecSM[4*((2*mSM-2)*i + j + mSM - 1) + 1] = (2*mSM - 1)*i + j + mSM - 1 + 2;
      eltsConnecSM[4*((2*mSM-2)*i + j + mSM - 1) + 2] = (2*mSM - 1)*(i + 1) + j + mSM - 1 + 2;
      eltsConnecSM[4*((2*mSM-2)*i + j + mSM - 1) + 3] = (2*mSM - 1)*(i + 1) + j + mSM - 1 + 1;
    }
  }

  MPI_Comm fvmComm = MPI_COMM_WORLD;

    
  cwipi::Mesh* sourceMesh = new cwipi::Mesh(MPI_COMM_WORLD,
                                             2 ,
                                             nVertexSM,
                                             nEltsSM,
                                             coordsSM,
                                             eltsConnecPointerSM,
                                             eltsConnecSM);




  fvm::fvm_writer_t *fvmWriterSM = fvm::fvm_writer_init("MeshSMSinX3D",
                                                        "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                                        "Ensight Gold",
                                                        "text",
                                                        fvm::FVM_WRITER_FIXED_MESH);
  
 
  fvm::fvm_writer_export_nodal(fvmWriterSM, 
  &(sourceMesh->getFvmNodal()));
  

  fvm::fvm_writer_finalize(fvmWriterSM);

  int nTM = 30;
  int mTM = 5;

  int nVertexTM = nTM*(2*mTM - 1);

  double *coordsTM = (double *) malloc(nVertexTM * 3 * sizeof(double));

  for(int i = 0; i < nTM; i++){
    for(int j = -mTM + 1; j < mTM ; j++){
      coordsTM[3*((2*mTM - 1)*i+j + mTM - 1)] = 1.*j/mTM;
      coordsTM[3*((2*mTM - 1)*i+j + mTM - 1) + 1] = 1.*i/nTM;
      coordsTM[3*((2*mTM - 1)*i+j + mTM - 1) + 2] = 1.*j/mSM;

    }
  }

  int nEltsTM = (nTM-1)*(2*mTM-2);
  int *eltsConnecPointerTM = (int *) malloc((nEltsTM+1) * sizeof(int));
  int *eltsConnecTM = (int *) malloc( 4*nEltsTM* sizeof(int));

  for(int i = 0; i < nEltsTM + 1; i++)
    eltsConnecPointerTM[i] = 4*i;

  for(int i = 0 ; i < nTM - 1 ; i++){
    for(int j = - mTM + 1; j < mTM - 1 ; j++){
      eltsConnecTM[4*((2*mTM - 2)*i + j + mTM - 1)] = (2*mTM - 1)*i + j + mTM - 1 + 1;
      eltsConnecTM[4*((2*mTM - 2)*i + j + mTM - 1) + 1] = (2*mTM - 1)*i + j + mTM - 1 + 2;
      eltsConnecTM[4*((2*mTM - 2)*i + j + mTM - 1) + 2] = (2*mTM - 1)*(i + 1) + j + mTM - 1 + 2;
      eltsConnecTM[4*((2*mTM - 2)*i + j + mTM - 1) + 3] = (2*mTM - 1)*(i + 1) + j + mTM - 1 + 1;
    }
  }

    
  cwipi::Mesh* targetMesh = new cwipi::Mesh(MPI_COMM_WORLD,
                                             2 ,
                                             nVertexTM,
                                             nEltsTM,
                                             coordsTM,
                                             eltsConnecPointerTM,
                                             eltsConnecTM);


  fvm::fvm_writer_t *fvmWriterTM = fvm::fvm_writer_init("MeshTMSinX3D",
                                                        "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                                        "Ensight Gold",
                                                        "text",
                                                        fvm::FVM_WRITER_FIXED_MESH);

  fvm::fvm_writer_export_nodal(fvmWriterTM, 
                               &(targetMesh->getFvmNodal()));
  
                               fvm::fvm_writer_finalize(fvmWriterTM);



  cwipi::ConservativeMesh* intersMesh = new cwipi::ConservativeMesh(fvmComm,
                                                                    *sourceMesh,
                                                                    *targetMesh,
                                                                    0.5);

  fvm::fvm_writer_t *fvmWriterIM = fvm::fvm_writer_init("MeshIMSinX3D",
                                                       "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                                        "Ensight Gold",
                                                        "text",
                                                        fvm::FVM_WRITER_FIXED_MESH);

  fvm::fvm_writer_export_nodal(fvmWriterIM, 
                               &(intersMesh->getIntersectionMesh()->getFvmNodal()));

                               fvm::fvm_writer_finalize(fvmWriterIM);


  
  /*  for(int i = 0; i < n ; i++){
    for(int j = 0; j < n ; j++){
      coordsTM[3*(m*i+j)] = 1 - j;
      coordsTM[3*(m*i+j) + 1] = -i;
      coordsTM[3*(m*i+j) + 2] = 0.001;
    }
  }

  delete targetMesh;

  targetMesh = new cwipi::Mesh(MPI_COMM_WORLD,
                               3 ,
                               nVertexTM,
                               nEltsTM,
                               coordsTM,
                               eltsConnecPointerTM,
                               eltsConnecTM);


  cwipi::ConservativeMesh* intersMesh2 = new cwipi::ConservativeMesh(fvmComm,
                                                                    *sourceMesh,
                                                                    *targetMesh,
                                                                    0.01);


  for(int i = 0; i < n ; i++){
    for(int j = 0; j < n ; j++){
      coordsSM[3*(m*i+j)] = 1 - j;
      coordsSM[3*(m*i+j) + 1] = -i;
      coordsSM[3*(m*i+j) + 2] = 0.001;
    }
  }
  
  delete sourceMesh;

  sourceMesh = new cwipi::Mesh(MPI_COMM_WORLD,
                               3 ,
                               nVertexSM,
                               nEltsSM,
                               coordsSM,
                               eltsConnecPointerSM,
                               eltsConnecSM);

  


  cwipi::ConservativeMesh* intersMesh3 = new cwipi::ConservativeMesh(fvmComm,
                                                                    *sourceMesh,
                                                                    *targetMesh,
                                                                    0.01);




  fvmWriterIM = fvm::fvm_writer_init("MeshIMSinX3D2",
                                     "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                     "Ensight Gold",
                                     "text",
                                     fvm::FVM_WRITER_FIXED_MESH);

  fvm::fvm_writer_export_nodal(fvmWriterIM, 
                               &(intersMesh2->getIntersectionMesh()->getFvmNodal()));

  fvm::fvm_writer_finalize(fvmWriterIM);


  fvmWriterIM = fvm::fvm_writer_init("MeshIMSinX3D3",
                                     "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                     "Ensight Gold",
                                     "text",
                                     fvm::FVM_WRITER_FIXED_MESH);

  fvm::fvm_writer_export_nodal(fvmWriterIM, 
                               &(intersMesh3->getIntersectionMesh()->getFvmNodal()));

  fvm::fvm_writer_finalize(fvmWriterIM);


  */





  free(srcName);
  free(eltsConnecPointerSM);
  free(eltsConnecPointerTM);
  free(eltsConnecSM);
  free(eltsConnecTM);
  free(coordsTM);
  free(coordsSM);
  delete sourceMesh;
  delete targetMesh;
  delete intersMesh;
  /* delete intersMesh2;
     delete intersMesh3;*/

  return EXIT_SUCCESS;
}