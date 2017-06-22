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
#include "creeMaillagePolygone2D.h"


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
  

  cwipi_solver_type_t solver_type;


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
  sprintf(fileName,"c_conservative_mesh_polygons%4.4d.txt",rank);

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
  




  int nVertexSM = 0;
  int nVertexTM = 0;
  int n = 0;
  int m = 0;
  int nEltsSM = 0;
  int nEltsTM = 0;
  double xmin = 0;
  double xmax = 0;
  double ymin = 0;
  double ymax = 0;
  
  bool begin = false;
                  
  double *coordsSM = NULL;
  int *eltsConnecPointerSM = NULL;
  int *eltsConnecSM = NULL;

  double *coordsTM = NULL;
  int *eltsConnecPointerTM = NULL;
  int *eltsConnecTM = NULL;

  const int order = 1;
  
  MPI_Comm fvmComm = MPI_COMM_WORLD;

  cwipi::oldMesh* sourceMesh;
  cwipi::oldMesh* targetMesh;
  cwipi::ConservativeMesh* intersMesh;

  fvm::fvm_writer_t *fvmWriterIM;

  char* name;
  
  /************* Tests Maillages aléatoires  ********************/  

    
  /*for (int i = 1 ; i < 6 ; i ++) {
    for (int j = 1 ; j < 6 ; j ++){ 
      for (int k = 1 ; k < 6 ; k ++){ 
        for (int l = 1 ; l < 6 ; l ++){
          
          for(double x = 1; x < 5 ; x++){
            for(double y = 1 ; y < 5 ; y ++){
              for(double s = 1 ; s < 5 ; s ++){  
              for(double t = 1 ; t < 5 ; t ++){*/
                  
  double x = 1.;
  double y = 1.;
  double s = 1.;
  double t = 1.;
  
  int i = 3;
  int j = 3;
  int k = 2;
  int l = 2; 
  
                /*if( x== 1. &&
                    y == 1. &&
                    s == 1. &&
                    t == 1. &&
                    
                    i == 1 &&
                    j == 1 &&
                    k == 1 &&
                    l == 1)
                  begin = true;
                  
                  
                  if (begin){*/
                    printf("x %f y %f s %f t %f i %d j %d k %d l %d \n",
                           x,y,s,t,i,j,k,l);
                    
                    nVertexSM = 0;
                    n = pow(2,i);
                    m = pow(2,j);
                    nEltsSM = 0;
                    xmin = -x;
                    xmax = x;
                    ymin = -y;
                    ymax = y;
                    
                  
                    
                    coordsSM = NULL;
                    eltsConnecPointerSM = NULL;
                    eltsConnecSM = NULL;
                    
                    creeMaillagePolygone2D(order,
                                           fvmComm,
                                           -1,
                                           1,
                                           -2,
                                           2,
                                           rank,
                                           n,
                                           m,
                                           &nVertexSM,
                                           &coordsSM,
                                           &nEltsSM,
                                           &eltsConnecPointerSM,
                                           &eltsConnecSM);
                    
                    sourceMesh = new cwipi::oldMesh(fvmComm,
                                                 2 ,
                                                 nVertexSM,
                                                 nEltsSM,
                                                 coordsSM,
                                                 eltsConnecPointerSM,
                                                 eltsConnecSM);
                    
                    
                    
                    nVertexTM = 0;
                    n = pow(2,k);
                    m = pow(2,l);
                    nEltsTM = 0;
                    xmin = -s;
                    xmax = s;
                    ymin = -t;
                    ymax = t;
                    
                    coordsTM = NULL;
                    eltsConnecPointerTM = NULL;
                    eltsConnecTM = NULL;
                    
                    creeMaillagePolygone2D(1,
                                           fvmComm,
                                           -2,
                                           2,
                                           -1,
                                           1,
                                           rank,
                                           n,
                                           m,
                                           &nVertexTM,
                                           &coordsTM,
                                           &nEltsTM,
                                           &eltsConnecPointerTM,
                                           &eltsConnecTM);
                    
                    
                    targetMesh = new cwipi::oldMesh(fvmComm,
                                                 2 ,
                                                 nVertexTM,
                                                 nEltsTM,
                                               coordsTM,
                                                 eltsConnecPointerTM,
                                                 eltsConnecTM);
                    
                    
                    
                    fvmWriterIM = fvm::fvm_writer_init("MeshSMPolygon",
                                                     "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                                       "Ensight Gold",
                                                       "text",
                                                       fvm::FVM_WRITER_FIXED_MESH);
                    
                    fvm::fvm_writer_export_nodal(fvmWriterIM, 
                                                 &(sourceMesh->getFvmNodal()));
                    
                    fvm::fvm_writer_finalize(fvmWriterIM);
                    
                    fvmWriterIM = fvm::fvm_writer_init("MeshTMPolygon",
                                                       "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                                       "Ensight Gold",
                                                       "text",
                                                       fvm::FVM_WRITER_FIXED_MESH);
                    
                    fvm::fvm_writer_export_nodal(fvmWriterIM, 
                                                 &(targetMesh->getFvmNodal()));
                    
                    fvm::fvm_writer_finalize(fvmWriterIM);
                    
                    if(sourceMesh->getNVertex() > targetMesh->getNVertex())                                         
                      intersMesh = new cwipi::ConservativeMesh(fvmComm,
                                                             *sourceMesh,
                                                             *targetMesh,
                                                             0.1);
                    else
                      intersMesh = new cwipi::ConservativeMesh(fvmComm,
                                                             *targetMesh,
                                                             *sourceMesh,
                                                             0.1);
                      
                    
                    
                  /*sprintf(name,
                    "MeshIMPolygon_x%f_y%f_s%f_t%f_i%d_j%d_k%d_l%d",
                    x,y,s,y,i,j,k,l);*/
                    
                    fvmWriterIM = fvm::fvm_writer_init("MeshIMPolygon",
                                                     "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                                     "Ensight Gold",
                                                     "text",
                                                     fvm::FVM_WRITER_FIXED_MESH);
                  
                  fvm::fvm_writer_export_nodal(fvmWriterIM, 
                                               &(intersMesh->getIntersectionMesh()->getFvmNodal()));
                  
                  fvm::fvm_writer_finalize(fvmWriterIM);
                  /*}
              }
              }
            }
          }
        }
      }
    }
    }*/
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

  return EXIT_SUCCESS;
}
