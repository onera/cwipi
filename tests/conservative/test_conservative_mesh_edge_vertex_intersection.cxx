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
#include "mesh.hxx"
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
  sprintf(fileName,"c_conservative_mesh_intersection%4.4d.txt",rank);

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

  /**************  Maillages aretes confondus  ****************/
  /*
     11 : premier sommet du maillage source
     21 : second sommet du maillage cible

    


              /\
        _____/__\__
       |    /|   \/|    
       |   / |   /\|   
       |  /  |  /  |\  
       | /   | /   |/   
       |/____|/____/
       |     |    /|
      /|    /|   / |
     / |   / |  /  |
    /  |  /  | /   |
   /   |_/___|/____|
   \    /    /
    \  /    /
     \/    /
      \   /
       \_/
        
   */


  int nVertexSM = 9;
  
  double *coordsSM = (double *) malloc(nVertexSM * 3 * sizeof(double));
  
  coordsSM[0] = 0;
  coordsSM[1] = 0;
  coordsSM[2] = 0;

  coordsSM[3] = 1;
  coordsSM[4] = 0;
  coordsSM[5] = 0;

  coordsSM[6] = 2;
  coordsSM[7] = 0;
  coordsSM[8] = 0;

  coordsSM[9]  = 0;
  coordsSM[10] = 1;
  coordsSM[11] = 0;

  coordsSM[12] = 1;
  coordsSM[13] = 1;
  coordsSM[14] = 0;

  coordsSM[15] = 2;
  coordsSM[16] = 1;
  coordsSM[17] = 0;

  coordsSM[18] = 0;
  coordsSM[19] = 2;
  coordsSM[20] = 0;

  coordsSM[21] = 1;
  coordsSM[22] = 2;
  coordsSM[23] = 0;

  coordsSM[24] = 2;
  coordsSM[25] = 2;
  coordsSM[26] = 0;

  int nEltsSM = 4;
  int *eltsConnecPointerSM = (int *) malloc((nEltsSM+1) * sizeof(int));
  int *eltsConnecSM = (int *) malloc( 16* sizeof(int));

  eltsConnecPointerSM[0] = 0;
  eltsConnecPointerSM[1] = 4;
  eltsConnecPointerSM[2] = 8;
  eltsConnecPointerSM[3] = 12;
  eltsConnecPointerSM[4] = 16;

  eltsConnecSM[0] = 1;
  eltsConnecSM[1] = 2;
  eltsConnecSM[2] = 5;
  eltsConnecSM[3] = 4;

  eltsConnecSM[4] = 2;
  eltsConnecSM[5] = 3;
  eltsConnecSM[6] = 6;
  eltsConnecSM[7] = 5;

  eltsConnecSM[8] = 4;
  eltsConnecSM[9] = 5;
  eltsConnecSM[10] = 8;
  eltsConnecSM[11] = 7;

  eltsConnecSM[12] = 5;
  eltsConnecSM[13] = 6;
  eltsConnecSM[14] = 9;
  eltsConnecSM[15] = 8;

  MPI_Comm fvmComm = MPI_COMM_WORLD;

    
  cwipi::Mesh* sourceMesh = new cwipi::Mesh(MPI_COMM_WORLD,
                                             2 ,
                                             nVertexSM,
                                             nEltsSM,
                                             coordsSM,
                                             eltsConnecPointerSM,
                                             eltsConnecSM);


  int nVertexTM = 16;
  
  double *coordsTM = (double *) malloc(nVertexTM * 3 * sizeof(double));
  
  coordsTM[0] = -0.7;
  coordsTM[1] = 1.3;
  coordsTM[2] = 0;

  coordsTM[3] = -0.2;
  coordsTM[4] = 0.8;
  coordsTM[5] = 0;

  coordsTM[6] = 0.3;
  coordsTM[7] = 0.3;
  coordsTM[8] = 0;

  coordsTM[9]  = 0.8;
  coordsTM[10] = -0.2;
  coordsTM[11] = 0;

  coordsTM[12] = -0.2;
  coordsTM[13] = 1.8;
  coordsTM[14] = 0;

  coordsTM[15] = 0.3;
  coordsTM[16] = 1.3;
  coordsTM[17] = 0;

  coordsTM[18] = 0.8;
  coordsTM[19] = 0.8;
  coordsTM[20] = 0;

  coordsTM[21]  = 1.3;
  coordsTM[22] = 0.3;
  coordsTM[23] = 0;

  coordsTM[24] = 0.3;
  coordsTM[25] = 2.3;
  coordsTM[26] = 0;

  coordsTM[27] = 0.8;
  coordsTM[28] = 1.8;
  coordsTM[29] = 0;

  coordsTM[30] = 1.3;
  coordsTM[31] = 1.3;
  coordsTM[32] = 0;

  coordsTM[33]  = 1.8;
  coordsTM[34] = 0.8;
  coordsTM[35] = 0;

  coordsTM[36] = 0.8;
  coordsTM[37] = 2.8;
  coordsTM[38] = 0;

  coordsTM[39] = 1.3;
  coordsTM[40] = 2.3;
  coordsTM[41] = 0;

  coordsTM[42] = 1.8;
  coordsTM[43] = 1.8;
  coordsTM[44] = 0;

  coordsTM[45]  = 2.3;
  coordsTM[46] = 1.3;
  coordsTM[47] = 0;


  int nEltsTM = 9;
  int *eltsConnecPointerTM = (int *) malloc((nEltsTM+1) * sizeof(int));
  int *eltsConnecTM = (int *) malloc( 36* sizeof(int));

  eltsConnecPointerTM[0] = 0;
  eltsConnecPointerTM[1] = 4;
  eltsConnecPointerTM[2] = 8;
  eltsConnecPointerTM[3] = 12;
  eltsConnecPointerTM[4] = 16;
  eltsConnecPointerTM[5] = 20;
  eltsConnecPointerTM[6] = 24;
  eltsConnecPointerTM[7] = 28;
  eltsConnecPointerTM[8] = 32;
  eltsConnecPointerTM[9] = 36;

  eltsConnecTM[0] = 1;
  eltsConnecTM[1] = 2;
  eltsConnecTM[2] = 6;
  eltsConnecTM[3] = 5;

  eltsConnecTM[4] = 2;
  eltsConnecTM[5] = 3;
  eltsConnecTM[6] = 7;
  eltsConnecTM[7] = 6;

  eltsConnecTM[8] = 3;
  eltsConnecTM[9] = 4;
  eltsConnecTM[10] = 8;
  eltsConnecTM[11] = 7;

  eltsConnecTM[12] = 5;
  eltsConnecTM[13] = 6;
  eltsConnecTM[14] = 10;
  eltsConnecTM[15] = 9;

  eltsConnecTM[16] = 6;
  eltsConnecTM[17] = 7;
  eltsConnecTM[18] = 11;
  eltsConnecTM[19] = 10;

  eltsConnecTM[20] = 7;
  eltsConnecTM[21] = 8;
  eltsConnecTM[22] = 12;
  eltsConnecTM[23] = 11;

  eltsConnecTM[24] = 9;
  eltsConnecTM[25] = 10;
  eltsConnecTM[26] = 14;
  eltsConnecTM[27] = 13;

  eltsConnecTM[28] = 10;
  eltsConnecTM[29] = 11;
  eltsConnecTM[30] = 15;
  eltsConnecTM[31] = 14;

  eltsConnecTM[32] = 11;
  eltsConnecTM[33] = 12;
  eltsConnecTM[34] = 16;
  eltsConnecTM[35] = 15;


  cwipi::Mesh* targetMesh = new cwipi::Mesh( fvmComm,
                                             2 ,
                                             nVertexTM,
                                             nEltsTM,
                                             coordsTM,
                                             eltsConnecPointerTM,
                                             eltsConnecTM);



  cwipi::ConservativeMesh* intersMesh = new cwipi::ConservativeMesh(fvmComm,
                                                                    *targetMesh,
                                                                    *sourceMesh,
                                                                    0.01);

  cwipi::ConservativeMesh* intersMesh2 = new cwipi::ConservativeMesh(fvmComm,
                                                                    *targetMesh,
                                                                    *sourceMesh,
                                                                    0.01);


fvm::fvm_writer_t *fvmWriterIM = fvm::fvm_writer_init("MeshIMEdgeVertInters",
                                                        "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                                        "Ensight Gold",
                                                        "text",
                                                        fvm::FVM_WRITER_FIXED_MESH);

  fvm::fvm_writer_export_nodal(fvmWriterIM, 
                               &(intersMesh->getIntersectionMesh()->getFvmNodal()));

  fvm::fvm_writer_finalize(fvmWriterIM);


  fvmWriterIM = fvm::fvm_writer_init("MeshIMEdgeVertInters2",
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
  delete intersMesh;
  delete intersMesh2;

  return EXIT_SUCCESS;
}
