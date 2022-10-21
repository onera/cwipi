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
 * Read mesh
 *
 * parameters:
 *   f                   <-- Mesh file
 *   dimension           --> Dimension
 *   n_vtx               <-- number of vertices
 *   n_elmt              <-- number of elements
 *   nConnecVertex       <-- size of connectivity
 *   coords              --> vertices coordinates
 *   connecPointer       --> connectivity index
 *   connec              --> connectivity
 *---------------------------------------------------------------------*/

static int read_mesh(FILE    *f,
                     int     dimension,
                     int     n_vtx,
                     int     n_face,
                     int     n_elt,
                     int     lface_connec,
                     int     lcell_connec,
                     double *coords,
                     int    *face_vtx_idx,
                     int    *face_vtx,
                     int    *cell_face_idx,
                     int    *cell_face) {
  int i, j, r;

  // Read coordinates
  for (i = 0 ; i < n_vtx ; i++) {
    for (j = 0 ; j < dimension ; j++) {
      r = fscanf(f, "%lf", coords + i * dimension + j);
      if (r == EOF)
        return EXIT_FAILURE;
    }
  }

  // Read face -> vertex connectivity index
  for (i = 0 ; i < n_face + 1 ; i++) {
    r = fscanf(f, "%d", face_vtx_idx + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read face -> vertex connectivity
  for (i = 0 ; i < lface_connec ; i++) {
    r = fscanf(f, "%d", face_vtx + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read cell -> face connectivity index
  for (i = 0 ; i < n_elt + 1 ; i++) {
    r = fscanf(f, "%d", cell_face_idx + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read cell -> face connectivity
  for (i = 0 ; i < lcell_connec ; i++) {
    r = fscanf(f, "%d", cell_face + i);
    //if(cell_face[i]<0) printf("cell_face[%i] %i\n",i,cell_face[i]);
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


  /* Initialization
   * -------------- */

  int n_code = 0;
  const char **code_names = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;

  if (rank == 0) {
    n_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code1_cell_faces";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }

  if (rank == 1) {
    n_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }


  times_init = malloc(sizeof(double) * n_code);

  for (int i = 0 ; i < n_code ; i++) {
    times_init[i] = 0;
  }

  printf("CWIPI Initialization rank %i\n", rank);
  CWP_client_Init(n_code,
                  code_names,
                  is_coupled_rank,
                  times_init);

  /* Finalize
   * -------- */

  char cpl_id1[] = "cpl_code1_code2";

  printf("Coupling creation\n");

  if (rank == 0) {
    CWP_client_Cpl_create("code1_cell_faces", cpl_id1, "code2", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                          CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_DBBTREE, 1,
                          CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }

  if (rank == 1) {
    CWP_client_Cpl_create("code2", cpl_id1, "code1_cell_faces", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                          CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_DBBTREE, 1,
                          CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }
  printf("Coupling created \n");

  /* Building of the local mesh */

  int dimension    = 3;             // Dimension of the space
  int n_vtx        = 0;             // Number of points in the mesh
  int n_face       = 0;             // Number of face
  int n_elmt       = 0;             // Number of cells
  int lface_connec = 0;
  int lcell_connec = 0;

  double *coords        = NULL;         // Coordinates of the points
  int    *face_vtx_idx  = NULL;
  int    *face_vtx      = NULL;
  int    *cell_face_idx = NULL;
  int    *cell_face     = NULL;

  if (rank == 0) printf("        Read mesh\n");

  FILE *mesh_file;
  mesh_file = fopen("../meshes/mesh_poly_d1", "r");

  fscanf(mesh_file, "%d %d %d %d %d %d",
         &dimension,
         &n_vtx,
         &n_face,
         &n_elmt,
         &lface_connec,
         &lcell_connec);

  coords        = (double *) malloc(dimension * n_vtx * sizeof(double));
  face_vtx_idx  = (int    *) malloc((n_face + 1)      * sizeof(int   ));
  face_vtx      = (int    *) malloc(lface_connec      * sizeof(int   ));
  cell_face_idx = (int    *) malloc((n_elmt + 1)      * sizeof(int   ));
  cell_face     = (int    *) malloc(lcell_connec      * sizeof(int   ));

  read_mesh(mesh_file,
            dimension,
            n_vtx,
            n_face,
            n_elmt,
            lface_connec,
            lcell_connec,
            coords,
            face_vtx_idx,
            face_vtx,
            cell_face_idx,
            cell_face);

  fclose(mesh_file);

  if (rank == 0) {

    printf("Visu Setting\n");
    CWP_client_Visu_set("code1_cell_faces", cpl_id1, 1.0, CWP_VISU_FORMAT_ENSIGHT, "binary");
    printf("Visu Set\n");

    CWP_g_num_t *global_num_vtx = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_vtx);
    for (int i = 0; i < n_vtx; i++) {
      global_num_vtx[i] = i + 1;
    }

    printf("vtx_set\n");
    CWP_client_Mesh_interf_vtx_set("code1_cell_faces", cpl_id1, 0, n_vtx, coords, global_num_vtx);

    CWP_g_num_t *global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elmt);
    for (int i = 0; i < n_elmt; i++) {
      global_num[i] = i + 1;
    }

    printf("Cell_face Add and Setting\n");
    CWP_client_Mesh_interf_from_cellface_set("code1_cell_faces",
                                             cpl_id1,
                                             0,
                                             n_elmt,
                                             cell_face_idx,
                                             cell_face,
                                             n_face,
                                             face_vtx_idx,
                                             face_vtx,
                                             global_num); // TO DO: deal with NULL tab passed !


    printf("Interface Mesh deletion\n");
    CWP_client_Mesh_interf_del("code1_cell_faces", cpl_id1);
    printf("Interface Mesh deleted\n");
  }


  fflush(stdout);

  CWP_client_Finalize();

  PDM_MPI_Finalize();

  free(code_names);
  free(is_coupled_rank);
  free(times_init);
  free(coords       );
  free(face_vtx_idx );
  free(face_vtx     );
  free(cell_face_idx);
  free(cell_face    );

  return 0;
}
