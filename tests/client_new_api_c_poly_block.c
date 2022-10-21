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
#include <string.h>
#include <math.h>

#include "cwp.h"
#include "pdm_io.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "client.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CWP_HEADER_SIZE    32

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

static int read_mesh_dim(FILE *f,
                         int *dimension,
                         int *nVertex,
                         int *nFace,
                         int *nElt,
                         int *lFaceConnec,
                         int *lCellConnec) {
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

static int read_mesh(FILE *f,
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
                     int *cellFace) {
  int i, j, r;

  // Read coordinates
  for (i = 0 ; i < nVertex ; i++) {
    for (j = 0 ; j < dimension ; j++) {
      r = fscanf(f, "%lf", coords + i * dimension + j);
      if (r == EOF)
        return EXIT_FAILURE;
    }
  }

  // Read face -> vertex connectivity index
  for (i = 0 ; i < nFace + 1 ; i++) {
    r = fscanf(f, "%d", faceVertexIdx + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read face -> vertex connectivity
  for (i = 0 ; i < lFaceConnec ; i++) {
    r = fscanf(f, "%d", faceVertex + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read cell -> face connectivity index
  for (i = 0 ; i < nElt + 1 ; i++) {
    r = fscanf(f, "%d", cellFaceIdx + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read cell -> face connectivity
  for (i = 0 ; i < lCellConnec ; i++) {
    r = fscanf(f, "%d", cellFace + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


/*=============================================================================
 * Util functions
 *============================================================================*/

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -c     Filename of the server configuration file.\n\n"
     "  -h     This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
 int            argc,
 char         **argv,
 char         **config  // filename for server ip adresses + ports
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *config = argv[i];
      }
    }

    else
      _usage(EXIT_FAILURE);
    i++;

  }

}

/*=============================================================================
 * Main
 *============================================================================*/

int
main
(
 int argc,
 char *argv[]
)
{
  // default
  char *config     = NULL;

  _read_args(argc,
             argv,
             &config);

  if (config == NULL) {
    config = (char *) "cwp_config_srv.txt";
  }

  // mpi
  int rank;
  int comm_world_size;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &rank);
  PDM_MPI_Comm_size(comm, &comm_world_size);

  // read host_name and port from config file
  // --> open
  PDM_io_file_t *read = NULL;
  PDM_l_num_t    ierr;

  PDM_io_open(config,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPI_SIMPLE,
              PDM_IO_MOD_READ,
              PDM_IO_NATIVE,
              comm,
              -1.,
              &read,
              &ierr);

  // --> global read offset in header
  char *buffer = malloc(CWP_HEADER_SIZE+1);
  for (int i = 0; i < CWP_HEADER_SIZE; i++) {
    buffer[i] = '\0';
  }

  PDM_io_global_read(read,
                     CWP_HEADER_SIZE * sizeof(char),
                     1,
                     buffer);

  char div_line[] = "\n";
  char *line = strtok(buffer, div_line);
  line = strtok(NULL, div_line);

  char div_word[] = " ";
  char *word = strtok(line, div_word);
  word = strtok(NULL, div_word);

  int offset = atoi(word);

  // --> block read hostname/port
  char *data = malloc(offset+1);

  for (int i = 0; i < offset+1; i++) {
    data[i] = '\0';
  }

  int one = 1;
  PDM_g_num_t rank_gnum = (PDM_g_num_t) (rank+1);

  PDM_io_par_interlaced_read(read,
                             PDM_STRIDE_CST_INTERLACED,
             (PDM_l_num_t *) &one,
               (PDM_l_num_t) offset,
               (PDM_l_num_t) one,
                             &rank_gnum,
                             data);

  char div[] = "/";
  char *str = strtok(data, div);
  char *server_name = malloc(strlen(str)+1);
  memcpy(server_name, str, strlen(str)+1);
  str = strtok(NULL, div);
  int server_port = atoi(str);

  // --> close
  PDM_io_close(read);
  PDM_io_free(read);

  // connect
  if (CWP_client_connect(server_name, server_port, CWP_CLIENTFLAG_VERBOSE) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "Client connexion failed\n");
    return -1;
  }

  FILE *meshFile;

  meshFile = fopen("../meshes/mesh_poly_d1", "r");

  int n_partition = 0;
  const int two = 2;
  while (two * pow(n_partition, two) < comm_world_size) n_partition++;

  int n2 = (int) (two * pow(n_partition, two));

  if (n2 != comm_world_size) {
    if (rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
    PDM_MPI_Finalize();
    return EXIT_SUCCESS;
  }

  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &comm_world_size);

  char *srcName = (char *) malloc(sizeof(char) * (strlen(__FILE__) + 1));
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
  const char **codeNames = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;

  if (rank == 0) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "cpoly";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }

  if (rank == 1) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code_name);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }

  times_init = malloc(sizeof(double) * n_code_name);

  for (int i = 0 ; i < n_code_name ; i++) {
    times_init[i] = 0;
  }

  CWP_client_Init(n_code_name,
                  codeNames,
                  is_coupled_rank,
                  times_init);

  char cpl_id1[] = "cpl_code1_code2";

  printf("Coupling creation\n");
  if (rank == 0) {
    CWP_client_Cpl_create("cpoly", cpl_id1, "code2", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                          CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, 1,
                          CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }

  if (rank == 1) {
    CWP_client_Cpl_create("code2", cpl_id1, "cpoly", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                          CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, 1,
                          CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }
  printf("Coupling created\n");

  /* Building of the local mesh */

  int dimension = 0;             // Dimension of the space
  int nVertex = 0;               // Number of points in the mesh
  int nFace = 0;                 // Number of face
  int nElements = 0;             // Number of cells
  int lFaceConnec = 0;
  int lCellConnec = 0;

  double *coords        = NULL;         // Coordinates of the points
  int    *faceVertexIdx = NULL;
  int    *faceVertex    = NULL;
  int    *cellFaceIdx   = NULL;
  int    *cellFace      = NULL;

  if (rank == 0)
    printf("        Read mesh\n");

  read_mesh_dim(meshFile, &dimension, &nVertex, &nFace, &nElements, &lFaceConnec, &lCellConnec);

  coords        = (double *) malloc(dimension * nVertex * sizeof(double));
  faceVertexIdx = (int    *) malloc((nFace + 1)         * sizeof(int   ));
  faceVertex    = (int    *) malloc(lFaceConnec         * sizeof(int   ));
  cellFaceIdx   = (int    *) malloc((nElements + 1)     * sizeof(int   ));
  cellFace      = (int    *) malloc(lCellConnec         * sizeof(int   ));

  read_mesh(meshFile,
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

  if (rank == 0) {
    printf("Visu Setting\n");
    CWP_client_Visu_set("cpoly", cpl_id1, 1, CWP_VISU_FORMAT_ENSIGHT, "binary");
    printf("Visu Set\n");

    printf("vtx_set\n");
    CWP_client_Mesh_interf_vtx_set("cpoly", cpl_id1, 0, nVertex, coords, NULL);

    printf("3D Cell Polyhedra Block Add\n");
    int block_id = CWP_Mesh_interf_block_add("cpoly", cpl_id1, CWP_BLOCK_CELL_POLY);

    printf("3D Cell Polyhedra Block Set\n");
    CWP_client_Mesh_interf_c_poly_block_set("cpoly", cpl_id1, 0, block_id,
                                            nElements,
                                            nFace,
                                            faceVertexIdx,
                                            faceVertex,
                                            cellFaceIdx,
                                            cellFace,
                                            NULL);

    CWP_g_num_t *cellGnum = malloc(sizeof(CWP_g_num_t) * nElements);
    int getNElements = -1;
    int getNFace = -1;
    int *getFaceVertexIdx = malloc(sizeof(int) * (nFace + 1));
    int *getFaceVertex    = malloc(sizeof(int) * getFaceVertexIdx[nFace]);
    int *getCellFaceIdx   = malloc(sizeof(int) * (nElements + 1));
    int *getCellFace      = malloc(sizeof(int) * getCellFaceIdx[nElements]);
    CWP_client_Mesh_interf_c_poly_block_get("cpoly", cpl_id1, 0, block_id,
                                            &getNElements,
                                            &getNFace,
                                            &getFaceVertexIdx,
                                            &getFaceVertex,
                                            &getCellFaceIdx,
                                            &getCellFace,
                                            &cellGnum);

    // Check output
    printf("nElements same ? %d\n", getNElements == nElements);
    printf("nFaces same ? %d\n", getNFace == nFace);
    int equal = -1;
    for (int i = 0; i < nElements + 1; i++) {
      equal = (getFaceVertexIdx[i] == faceVertexIdx[i]);
      if (equal == 0) {
        break;
      }
    }
    printf("FaceVertexIdx same ? %d\n", equal);
    equal = -1;
    for (int i = 0; i < faceVertexIdx[nElements]; i++) {
      equal = (getFaceVertex[i] == faceVertex[i]);
      if (equal == 0) {
        break;
      }
    }
    printf("FaceVertex same ? %d\n", equal);
    equal = -1;
    for (int i = 0; i < nFace + 1; i++) {
      equal = (getCellFaceIdx[i] == cellFaceIdx[i]);
      if (equal == 0) {
        break;
      }
    }
    printf("CellFaceIdx same ? %d\n", equal);
    equal = -1;
    for (int i = 0; i < cellFaceIdx[nFace]; i++) {
      equal = (getCellFace[i] == cellFace[i]);
      if (equal == 0) {
        break;
      }
    }
    printf("CellFace same ? %d\n", equal);

    printf("Interface Mesh deletion\n");
    CWP_client_Mesh_interf_del("cpoly", cpl_id1);
    printf("Interface Mesh deleted\n");
    fflush(stdout);
  }

  if (rank == 0) {
    CWP_client_Cpl_del("cpoly", cpl_id1);
  }

  if (rank == 1) {
    CWP_client_Cpl_del("code2", cpl_id1);
  }

  free(coords       );
  free(faceVertexIdx);
  free(faceVertex   );
  free(cellFaceIdx  );
  free(cellFace     );

  CWP_client_Finalize();
  PDM_MPI_Finalize();

  free(srcName);
  free(codeNames);
  free(is_coupled_rank);
  free(times_init);

  return 0;
}
