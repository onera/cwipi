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

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include "cwp.h"
#include "pdm.h"
#include "pdm_part.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_logging.h"
#include "client.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_io.h"
#include "pdm_mpi.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CWP_HEADER_SIZE    32

/*=============================================================================
 * Private functions definitions
 *============================================================================*/

static void
create_dcube_from_nodal
        (
                PDM_MPI_Comm pdm_comm,
                PDM_Mesh_nodal_elt_t element_type,
                PDM_g_num_t n_vtx_seg,
                double length,
                double xmin,
                double ymin,
                double zmin,
                int **n_vtx,
                int **n_faces,
                int **n_cells,
                double ***coord,
                int ***face_vtx_idx,
                int ***face_vtx,
                int ***cell_face_idx,
                int ***cell_face,
                PDM_g_num_t ***vtx_ln_to_gn,
                PDM_g_num_t ***cell_ln_to_gn
        ) {
  int n_part = 1;

  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        xmin,
                                                        ymin,
                                                        zmin,
                                                        element_type,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_ordering_set(dcube, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build(dcube);

  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_g_num_t *vertex_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double *d_vertex_coord = PDM_DMesh_nodal_vtx_get(dmn);
  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  int d_n_vertices = vertex_distrib[i_rank + 1] - vertex_distrib[i_rank];

  PDM_dmesh_nodal_to_dmesh_t *dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, pdm_comm, PDM_OWNERSHIP_KEEP);
  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);
  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);
  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  PDM_dmesh_t *dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

  // Récupérer les connectivités : PDM_dmesh_connectivity_get
  int *d_face_vertex_idx, *d_cell_face_idx;
  PDM_g_num_t *d_face_vertex, *d_cell_face;
  int n_face_group = 0;
  int d_face_group_idx[1] = {0};
  PDM_g_num_t *d_face_group = NULL;
  int d_n_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                            &d_face_vertex,
                                            &d_face_vertex_idx,
                                            PDM_OWNERSHIP_KEEP);

  int d_n_cell = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                            &d_cell_face,
                                            &d_cell_face_idx,
                                            PDM_OWNERSHIP_KEEP);

  // Partitionnement
  int *d_cell_part = (int *) malloc(sizeof(int) * d_n_cell);

  PDM_part_split_t method = PDM_PART_SPLIT_PTSCOTCH;
  PDM_part_t *ppart_id = PDM_part_create(pdm_comm,
                                         method,
                                         "PDM_PART_RENUM_CELL_NONE",
                                         "PDM_PART_RENUM_FACE_NONE",
                                         0,
                                         NULL,
                                         0,
                                         NULL,
                                         n_part,
                                         d_n_cell,
                                         d_n_face,
                                         d_n_vertices,
                                         n_face_group,
                                         d_cell_face_idx,
                                         d_cell_face,
                                         NULL,
                                         NULL,
                                         0,
                                         d_cell_part,
                                         NULL,
                                         d_face_vertex_idx,
                                         d_face_vertex,
                                         NULL,
                                         d_vertex_coord,
                                         NULL,
                                         d_face_group_idx,
                                         d_face_group);

  // PDM_log_trace_connectivity_long (d_cell_face_idx, d_cell_face, d_n_cell, "d_cell_face : ");
  // PDM_log_trace_connectivity_long (d_face_vertex_idx, d_face_vertex, d_n_face, "d_face_vertex : ");
  // PDM_log_trace_array_double (d_vertex_coord, 3 * d_n_vertices, "d_coords : ");

  // free(d_face_vertex_idx);
  // free(d_face_vertex);
  // free(d_face_group);
  // free(d_vertex_coord);

  free(d_cell_part);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);

  // Get connectivity
  // Cells
  *n_cells = (int *) malloc(sizeof(int) * n_part);
  *cell_face_idx = (int **) malloc(sizeof(int *) * n_part);
  *cell_face = (int **) malloc(sizeof(int *) * n_part);
  *cell_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  // Faces
  *n_faces = (int *) malloc(sizeof(int) * n_part);
  *face_vtx_idx = (int **) malloc(sizeof(int *) * n_part);
  *face_vtx = (int **) malloc(sizeof(int *) * n_part);

  // Vertices
  *n_vtx = (int *) malloc(sizeof(int) * n_part);
  *coord = (double **) malloc(sizeof(double *) * n_part);
  *vtx_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int i_part = 0 ; i_part < n_part ; i_part++) {

    int _n_cells;
    int _n_faces;
    int _n_face_part_bound;
    int _n_vtx;
    int _n_proc;
    int _n_total_part;
    int _s_cell_face;
    int _s_face_vtx;
    int _s_face_group;
    int _n_edge_group2;

    PDM_part_part_dim_get(ppart_id,
                          i_part,
                          &_n_cells,
                          &_n_faces,
                          &_n_face_part_bound,
                          &_n_vtx,
                          &_n_proc,
                          &_n_total_part,
                          &_s_cell_face,
                          &_s_face_vtx,
                          &_s_face_group,
                          &_n_edge_group2);

    int *_cell_tag;
    int *_cell_face_idx;
    int *_cell_face;
    int *_face_tag;
    int *_face_cell;
    int *_face_vtx_idx;
    int *_face_vtx;
    int *_face_part_bound_proc_idx;
    int *_face_part_bound_part_idx;
    int *_face_part_bound;
    int *_vtx_tag;
    int *_face_group_idx;
    int *_face_group;
    double *_vtx_coords;

    PDM_g_num_t *_cell_ln_to_gn;
    PDM_g_num_t *_face_ln_to_gn;
    PDM_g_num_t *_vtx_ln_to_gn;
    PDM_g_num_t *_face_group_ln_to_gn;

    PDM_part_part_val_get(ppart_id,
                          i_part,
                          &_cell_tag,
                          &_cell_face_idx,
                          &_cell_face,
                          &_cell_ln_to_gn,
                          &_face_tag,
                          &_face_cell,
                          &_face_vtx_idx,
                          &_face_vtx,
                          &_face_ln_to_gn,
                          &_face_part_bound_proc_idx,
                          &_face_part_bound_part_idx,
                          &_face_part_bound,
                          &_vtx_tag,
                          &_vtx_coords,
                          &_vtx_ln_to_gn,
                          &_face_group_idx,
                          &_face_group,
                          &_face_group_ln_to_gn);

    // Cells
    *n_cells[i_part] = _n_cells;
    *cell_face_idx[i_part] = (int *) malloc(sizeof(int) * (_n_cells + 1));
    *cell_face[i_part] = (int *) malloc(sizeof(int) * _s_cell_face);
    *cell_ln_to_gn[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_cells);

    memcpy(*cell_face_idx[i_part], _cell_face_idx, (_n_cells + 1) * sizeof(int));
    memcpy(*cell_face[i_part], _cell_face, _s_cell_face * sizeof(int));
    memcpy(*cell_ln_to_gn[i_part], _cell_ln_to_gn, _n_cells * sizeof(PDM_g_num_t));

    PDM_log_trace_connectivity_int(*cell_face_idx[i_part], *cell_face[i_part], *n_cells[i_part], "final cell_face");

    // Faces
    *n_faces[i_part] = _n_faces;
    *face_vtx_idx[i_part] = (int *) malloc(sizeof(int) * (_n_faces + 1));
    *face_vtx[i_part] = (int *) malloc(sizeof(int) * _s_face_vtx);

    memcpy(*face_vtx_idx[i_part], _face_vtx_idx, (_n_faces + 1) * sizeof(int));
    memcpy(*face_vtx[i_part], _face_vtx, _s_face_vtx * sizeof(int));

    PDM_log_trace_connectivity_int(*face_vtx_idx[i_part], *face_vtx[i_part], *n_faces[i_part], "final face_vertex");

    // Vertices
    *n_vtx[i_part] = _n_vtx;
    *coord[i_part] = (double *) malloc(sizeof(double) * (3 * _n_vtx));
    *vtx_ln_to_gn[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_vtx);

    memcpy(*coord[i_part], _vtx_coords, 3 * _n_vtx * sizeof(double));
    memcpy(*vtx_ln_to_gn[i_part], _vtx_ln_to_gn, _n_vtx * sizeof(PDM_g_num_t));
  }

  PDM_part_free(ppart_id);
  PDM_dcube_nodal_gen_free(dcube);
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
  memcpy(server_name, str, strlen(str));
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

  // Input
  int n_part = 1;

  CWP_Field_exch_t exchDirection[2];

  bool cond_code1 = rank == 1;
  bool cond_code2 = rank == 0 || rank == 2;
  exchDirection[0] = CWP_FIELD_EXCH_SEND;
  exchDirection[1] = CWP_FIELD_EXCH_RECV;

  bool cond_both = cond_code1 && cond_code2;

  PDM_Mesh_nodal_elt_t element_type_code1 = PDM_MESH_NODAL_HEXA8;
  PDM_Mesh_nodal_elt_t element_type_code2 = PDM_MESH_NODAL_HEXA8;

  CWP_Block_t element_type_code1_cwp = CWP_BLOCK_CELL_HEXA8;
  CWP_Block_t element_type_code2_cwp = CWP_BLOCK_CELL_HEXA8;

  int n_vtx_seg_code1 = 4, n_vtx_seg_code2 = 3;
  double x_min_code1 = 0., x_min_code2 = 0.;
  double y_min_code1 = 0., y_min_code2 = 0.;
  double z_min_code1 = 0., z_min_code2 = 0.;


  // Define the number of codes per rank
  int n_code;
  if (cond_both) n_code = 2;
  else if (cond_code1 || cond_code2) n_code = 1;
  else n_code = 0;

  int *code_id = (int *) malloc(n_code * sizeof(int));
  const char **code_names = (const char **) malloc(n_code * sizeof(char *));
  const char **coupled_code_names = (const char **) malloc(n_code * sizeof(char *));
  double *time_init = (double *) malloc(n_code * sizeof(double));

  PDM_Mesh_nodal_elt_t *element_type = (PDM_Mesh_nodal_elt_t *) malloc(n_code * sizeof(PDM_Mesh_nodal_elt_t));

  CWP_Block_t *element_type_cwp = (CWP_Block_t *) malloc(n_code * sizeof(CWP_Block_t));

  PDM_g_num_t *n_vtx_seg = (PDM_g_num_t *) malloc(n_code * sizeof(PDM_g_num_t));

  double *x_min = (double *) malloc(n_code * sizeof(double));
  double *y_min = (double *) malloc(n_code * sizeof(double));
  double *z_min = (double *) malloc(n_code * sizeof(double));

  CWP_Status_t *is_active_rank = (CWP_Status_t *) malloc(n_code * sizeof(CWP_Status_t));

  // Define which rank works for which code
  if (cond_both) {
    code_id[0] = 1;
    code_id[1] = 2;
    code_names[0] = "code1";
    code_names[1] = "code2";
    coupled_code_names[0] = "code2";
    coupled_code_names[1] = "code1";
    time_init[0] = 0.;
    time_init[1] = 0.;
    is_active_rank[0] = CWP_STATUS_ON;
    is_active_rank[1] = CWP_STATUS_ON;
    element_type[0] = element_type_code1;
    element_type[1] = element_type_code2;
    element_type_cwp[0] = element_type_code1_cwp;
    element_type_cwp[1] = element_type_code2_cwp;
    n_vtx_seg[0] = n_vtx_seg_code1;
    n_vtx_seg[1] = n_vtx_seg_code2;
    x_min[0] = x_min_code1;
    x_min[1] = x_min_code2;
    y_min[0] = y_min_code1;
    y_min[1] = y_min_code2;
    z_min[0] = z_min_code1;
    z_min[1] = z_min_code2;
  }

  else if (cond_code1) {
    code_id[0] = 1;
    code_names[0] = "code1";
    coupled_code_names[0] = "code2";
    time_init[0] = 0.;
    is_active_rank[0] = CWP_STATUS_ON;
    element_type[0] = element_type_code1;
    element_type_cwp[0] = element_type_code1_cwp;
    n_vtx_seg[0] = n_vtx_seg_code1;
    x_min[0] = x_min_code1;
    y_min[0] = y_min_code1;
    z_min[0] = z_min_code1;
  }

  else if (cond_code2) {
    code_id[0] = 2;
    code_names[0] = "code2";
    coupled_code_names[0] = "code1";
    time_init[0] = 0.;
    is_active_rank[0] = CWP_STATUS_ON;
    element_type[0] = element_type_code2;
    element_type_cwp[0] = element_type_code2_cwp;
    n_vtx_seg[0] = n_vtx_seg_code2;
    x_min[0] = x_min_code2;
    y_min[0] = y_min_code2;
    z_min[0] = z_min_code2;
  }

  // Init cwipi
  CWP_client_Init(n_code,
                  (const char **) code_names,
                  is_active_rank,
                  time_init);

  printf("%d --- CWIPI initialized\n", rank);

  // Create coupling and visu
  const char *cpl_name = "c_new_api_disjoint_comms";
  CWP_Spatial_interp_t interp_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_DBBTREE;

  for (int i_code = 0 ; i_code < n_code ; ++i_code) {

    CWP_client_Cpl_create(code_names[i_code],
                          cpl_name,
                          coupled_code_names[i_code],
                          CWP_INTERFACE_VOLUME,
                          CWP_COMM_PAR_WITH_PART,
                          interp_method,
                          n_part,
                          CWP_DYNAMIC_MESH_STATIC,
                          CWP_TIME_EXCH_USER_CONTROLLED);

    printf("%d --- Coupling created between %s and %s\n", rank, code_names[i_code], coupled_code_names[i_code]);
  }

  for (int i_code = 0 ; i_code < n_code ; ++i_code) {

    CWP_client_Visu_set(code_names[i_code],
                        cpl_name,
                        1,
                        CWP_VISU_FORMAT_ENSIGHT,
                        "text");

    printf("%d : %s --- Visu set\n", rank, code_names[i_code]);
  }

  // Create geometry
  int **n_vtx = (int **) malloc(n_code * sizeof(int **));
  int **n_faces = (int **) malloc(n_code * sizeof(int **));
  int **n_cells = (int **) malloc(n_code * sizeof(int **));
  int *n_blocks = (int *) malloc(n_code * sizeof(int *));

  double ***coord = (double ***) malloc(n_code * sizeof(double ***));

  int ***face_vtx_idx = (int ***) malloc(n_code * sizeof(int ***));
  int ***face_vtx = (int ***) malloc(n_code * sizeof(int ***));
  int ***cell_face_idx = (int ***) malloc(n_code * sizeof(int ***));
  int ***cell_face = (int ***) malloc(n_code * sizeof(int ***));

  PDM_l_num_t ***connec = (PDM_l_num_t ***) malloc(n_code * sizeof(PDM_l_num_t ***));

  PDM_l_num_t ***face_vtx_nb = malloc(n_code * sizeof(PDM_l_num_t ***));
  PDM_l_num_t ***cell_face_nb = malloc(n_code * sizeof(PDM_l_num_t ***));

  PDM_g_num_t ***vtx_ln_to_gn = (PDM_g_num_t ***) malloc(n_code * sizeof(PDM_g_num_t ***));
  PDM_g_num_t ***cell_ln_to_gn = (PDM_g_num_t ***) malloc(n_code * sizeof(PDM_g_num_t ***));

  PDM_Mesh_nodal_t **mesh_nodal = (PDM_Mesh_nodal_t **) malloc(n_code * sizeof(PDM_Mesh_nodal_t **));

  // Create local comm to code using comm_split

  PDM_MPI_Comm LocalComm;
  int color;
  if (cond_code2) {
    color = 1;
  } else {
    color = 0;
  }
  PDM_MPI_Comm_split(comm, color, 0, &LocalComm);

  for (int i_code = 0 ; i_code < n_code ; ++i_code) {

    face_vtx_nb[i_code] = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t **) * n_part);
    cell_face_nb[i_code] = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t **) * n_part);
    connec[i_code] = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t **) * n_part);

    create_dcube_from_nodal(LocalComm,
                            element_type[i_code],
                            n_vtx_seg[i_code],
                            1.,
                            x_min[i_code],
                            y_min[i_code],
                            z_min[i_code],
                            &(n_vtx[i_code]), &n_faces[i_code],
                            &n_cells[i_code],
                            &coord[i_code],
                            &face_vtx_idx[i_code],
                            &face_vtx[i_code],
                            &cell_face_idx[i_code],
                            &cell_face[i_code],
                            &vtx_ln_to_gn[i_code],
                            &cell_ln_to_gn[i_code]);

    printf("%d : %s --- dcube created\n", rank, code_names[i_code]);

    mesh_nodal[i_code] = PDM_Mesh_nodal_create(n_part, LocalComm);

    printf("pointer value of my mesh is %p\n", mesh_nodal[i_code]);

    for (int i_part = 0 ; i_part < n_part ; ++i_part) {
      face_vtx_nb[i_code][i_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t *) * n_faces[i_code][i_part]);
      cell_face_nb[i_code][i_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t *) * n_cells[i_code][i_part]);

      for (int i = 0 ; i < n_faces[i_code][i_part] ; i++) {
        face_vtx_nb[i_code][i_part][i] = face_vtx_idx[i_code][i_part][i + 1] - face_vtx_idx[i_code][i_part][i];
      }

      for (int i = 0 ; i < n_cells[i_code][i_part] ; i++) {
        cell_face_nb[i_code][i_part][i] = cell_face_idx[i_code][i_part][i + 1] - cell_face_idx[i_code][i_part][i];
      }

      // Set coords
      CWP_client_Mesh_interf_vtx_set(code_names[i_code],
                                     cpl_name,
                                     i_part,
                                     n_vtx[i_code][i_part],
                                     coord[i_code][i_part],
                                     vtx_ln_to_gn[i_code][i_part]);

      printf("%d : %s --- Points set for part %d %d\n", rank, code_names[i_code], i_part, n_vtx[i_code][i_part]);

      // 1 - Set connectivities from nodal
//        CWP_Mesh_interf_from_cellface_set(code_names[i_code], cpl_name, i_part, n_cells[i_code][i_part], cell_face_idx[i_code][i_part], cell_face[i_code][i_part],
//                                          n_faces[i_code][i_part], face_vtx_idx[i_code][i_part], face_vtx[i_code][i_part], cell_ln_to_gn[i_code][i_part]);

      // 2 - Set connectivities by reverting to a standard block from nodal
      PDM_Mesh_nodal_coord_set(mesh_nodal[i_code],
                               i_part,
                               n_vtx[i_code][i_part],
                               coord[i_code][i_part],
                               vtx_ln_to_gn[i_code][i_part],
                               PDM_OWNERSHIP_USER);

      PDM_Mesh_nodal_cell3d_cellface_add(mesh_nodal[i_code],
                                         i_part,
                                         n_cells[i_code][i_part],
                                         n_faces[i_code][i_part],
                                         face_vtx_idx[i_code][i_part],
                                         face_vtx_nb[i_code][i_part],
                                         face_vtx[i_code][i_part],
                                         cell_face_idx[i_code][i_part],
                                         cell_face_nb[i_code][i_part],
                                         cell_face[i_code][i_part],
                                         cell_ln_to_gn[i_code][i_part],
                                         PDM_OWNERSHIP_USER);
    }

    n_blocks[i_code] = PDM_Mesh_nodal_n_blocks_get(mesh_nodal[i_code]);

    for (int i_block = 0 ; i_block < n_blocks[i_code] ; ++i_block) {
      for (int i_part = 0 ; i_part < n_part ; ++i_part) {
        PDM_Mesh_nodal_block_std_get(mesh_nodal[i_code], i_block, i_part, &connec[i_code][i_part]);
        PDM_Mesh_nodal_g_num_in_block_compute(mesh_nodal[i_code], i_block, PDM_OWNERSHIP_USER);
        PDM_g_num_t *g_num = PDM_Mesh_nodal_block_g_num_get(mesh_nodal[i_code], i_block, i_part);

        int block_id = CWP_client_Mesh_interf_block_add(code_names[i_code],
                                                        cpl_name,
                                                        element_type_cwp[i_code]);

        CWP_client_Mesh_interf_block_std_set(code_names[i_code],
                                             cpl_name,
                                             i_part,
                                             block_id,
                                             element_type_cwp[i_code],
                                             n_cells[i_code][i_part],
                                             connec[i_code][i_part],
                                             g_num);
      }
    }

    CWP_client_Mesh_interf_finalize(code_names[i_code], cpl_name);

    printf("%d : %s --- Geometry set\n", rank, code_names[i_code]);
  }

  // Create and initialise Fields: code1 -> code2
  const char *field_name = "coo";

  double **send_values = (double **) malloc(n_code * sizeof(double *));
  double **recv_values = (double **) malloc(n_code * sizeof(double *));
  for (int i_code = 0 ; i_code < n_code ; ++i_code) {
    send_values[i_code] = (double *) malloc(3 * n_vtx[i_code][0] * sizeof(double));
    recv_values[i_code] = (double *) malloc(3 * n_vtx[i_code][0] * sizeof(double));
    if (code_id[i_code] == 1 && (exchDirection[0] == CWP_FIELD_EXCH_SEND)) {
      for (int i = 0 ; i < 3 * n_vtx[i_code][0] ; i++) {
        send_values[i_code][i] = coord[i_code][0][i];
      }
    }
    if (code_id[i_code] == 2 && (exchDirection[1] == CWP_FIELD_EXCH_RECV)) {
      for (int i = 0 ; i < 3 * n_vtx[i_code][0] ; i++) {
        recv_values[i_code][i] = -code_id[i_code];
      }
    }

    if (code_id[i_code] == 1 && (exchDirection[0] == CWP_FIELD_EXCH_RECV)) {
      for (int i = 0 ; i < 3 * n_vtx[i_code][0] ; i++) {
        recv_values[i_code][i] = -code_id[i_code];
      }
    }
    if (code_id[i_code] == 2 && (exchDirection[1] == CWP_FIELD_EXCH_SEND)) {
      for (int i = 0 ; i < 3 * n_vtx[i_code][0] ; i++) {
        send_values[i_code][i] = coord[i_code][0][i];
      }
    }
  }

  for (int i_code = 0 ; i_code < n_code ; ++i_code) {

    if (code_id[i_code] == 1) {

      CWP_client_Field_create(code_names[i_code],
                              cpl_name,
                              field_name,
                              CWP_DOUBLE,
                              CWP_FIELD_STORAGE_INTERLEAVED,
                              3,
                              CWP_DOF_LOCATION_NODE,
                              exchDirection[0],
                              CWP_STATUS_ON);

      if (exchDirection[0] == CWP_FIELD_EXCH_SEND) {

        CWP_client_Field_data_set(code_names[i_code],
                                  cpl_name,
                                  field_name,
                                  0,
                                  CWP_FIELD_MAP_SOURCE,
                                  send_values[i_code],
                                  3,
                                  n_vtx[i_code][0]);

      }

      else {

        CWP_client_Field_data_set(code_names[i_code],
                                  cpl_name,
                                  field_name,
                                  0,
                                  CWP_FIELD_MAP_TARGET,
                                  recv_values[i_code],
                                  3,
                                  n_vtx[i_code][0]);

      }
    }

    if (code_id[i_code] == 2) {

      CWP_client_Field_create(code_names[i_code],
                              cpl_name,
                              field_name,
                              CWP_DOUBLE,
                              CWP_FIELD_STORAGE_INTERLEAVED,
                              3,
                              CWP_DOF_LOCATION_NODE,
                              exchDirection[1],
                              CWP_STATUS_ON);

      if (exchDirection[1] == CWP_FIELD_EXCH_RECV) {

        CWP_client_Field_data_set(code_names[i_code],
                                  cpl_name,
                                  field_name,
                                  0,
                                  CWP_FIELD_MAP_TARGET,
                                  recv_values[i_code],
                                  3,
                                  n_vtx[i_code][0]);
      }

      else {

        CWP_client_Field_data_set(code_names[i_code],
                                  cpl_name,
                                  field_name,
                                  0,
                                  CWP_FIELD_MAP_SOURCE,
                                  send_values[i_code],
                                  3,
                                  n_vtx[i_code][0]);
      }
    }
    printf("%d : %s --- Field created and data set\n", rank, code_names[i_code]);
  }

  PDM_MPI_Barrier(comm);

  // Compute weights
  for (int i_code = 0 ; i_code < n_code ; ++i_code) {
    printf("%d : %s --- Weights computed\n", rank, code_names[i_code]);
    fflush(stdout);

    CWP_client_Spatial_interp_weights_compute(code_names[i_code], cpl_name);

  }

  int n_computed_tgts = 0, n_uncomputed_tgts = 0, n_involved_srcs = 0;
  const int *computed_tgts = NULL, *uncomputed_tgts = NULL, *involved_srcs = NULL;

  if (exchDirection[0] == CWP_FIELD_EXCH_SEND) {
    if (cond_code2) {

      n_computed_tgts = CWP_client_N_computed_tgts_get("code2", cpl_name, field_name, 0);
      n_uncomputed_tgts = CWP_client_N_uncomputed_tgts_get("code2", cpl_name, field_name, 0);
      computed_tgts = CWP_client_Computed_tgts_get("code2", cpl_name, field_name, 0);
      uncomputed_tgts = CWP_client_Uncomputed_tgts_get("code2", cpl_name, field_name, 0);

      int i_code = n_code == 2 ? 1 : 0;

      printf("%d : %s --- n computed targets: %d\n", rank, code_names[i_code], n_computed_tgts);
      printf("%d : %s --- n uncomputed targets: %d\n", rank, code_names[i_code], n_uncomputed_tgts);
      if (n_computed_tgts != 0) {
        printf("%d : %s --- computed targets: ", rank, code_names[i_code]);
        for (int i = 0 ; i < n_computed_tgts ; ++i) printf("%d ", computed_tgts[i]);
        printf("\n");
      }
      if (n_uncomputed_tgts != 0) {
        printf("%d : %s --- uncomputed targets: ", rank, code_names[i_code]);
        for (int i = 0 ; i < n_uncomputed_tgts ; ++i) printf("%d ", uncomputed_tgts[i]);
        printf("\n");
      }
    }

    if (cond_code1) {

      n_involved_srcs = CWP_client_N_involved_srcs_get("code1", cpl_name, field_name, 0);
      involved_srcs = CWP_client_Involved_srcs_get("code1", cpl_name, field_name, 0);

      int i_code = n_code == 2 ? 1 : 0;

      printf("%d : %s --- n involved sources: %d\n", rank, code_names[i_code], n_involved_srcs);
      if (n_involved_srcs != 0) {
        printf("%d : %s --- involved sources: ", rank, code_names[i_code]);
        for (int i = 0 ; i < n_involved_srcs ; ++i) printf("%d ", involved_srcs[i]);
        printf("\n");
      }
    }
  }
  else {
    if (cond_code1) {

      n_computed_tgts = CWP_client_N_computed_tgts_get("code1", cpl_name, field_name, 0);
      n_uncomputed_tgts = CWP_client_N_uncomputed_tgts_get("code1", cpl_name, field_name, 0);
      computed_tgts = CWP_client_Computed_tgts_get("code1", cpl_name, field_name, 0);
      uncomputed_tgts = CWP_client_Uncomputed_tgts_get("code1", cpl_name, field_name, 0);

      printf("%d : %s --- n computed targets: %d\n", rank, code_names[0], n_computed_tgts);
      printf("%d : %s --- n uncomputed targets: %d\n", rank, code_names[0], n_uncomputed_tgts);
      if (n_computed_tgts != 0) {
        printf("%d : %s --- computed targets: ", rank, code_names[0]);
        for (int i = 0 ; i < n_computed_tgts ; ++i) printf("%d ", computed_tgts[i]);
        printf("\n");
      }
      if (n_uncomputed_tgts != 0) {
        printf("%d : %s --- uncomputed targets: ", rank, code_names[0]);
        for (int i = 0 ; i < n_uncomputed_tgts ; ++i) printf("%d ", uncomputed_tgts[i]);
        printf("\n");
      }
    }

    if (cond_code2) {

      n_involved_srcs = CWP_N_involved_srcs_get("code2", cpl_name, field_name, 0);
      involved_srcs = CWP_Involved_srcs_get("code2", cpl_name, field_name, 0);

      int i_code = n_code == 2 ? 1 : 0;

      printf("%d : %s --- n involved sources: %d\n", rank, code_names[i_code], n_involved_srcs);
      if (n_involved_srcs != 0) {
        printf("%d : %s --- involved sources: ", rank, code_names[i_code]);
        for (int i = 0 ; i < n_involved_srcs ; ++i) printf("%d ", involved_srcs[i]);
        printf("\n");
      }
    }
  }


  // // Send and receive field
  // for (int i_code = 0 ; i_code < n_code ; i_code++) {
  //   if (code_id[i_code] == 2) {
  //     if (exchDirection[1] == CWP_FIELD_EXCH_SEND) {

  //       CWP_client_Field_issend(code_names[i_code], cpl_name, field_name);

  //       printf("%d : %s --- Sent field\n", rank, code_names[i_code]);
  //     }
  //     else {

  //       CWP_client_Field_irecv(code_names[i_code], cpl_name, field_name);

  //       printf("%d : %s --- Received field\n", rank, code_names[i_code]);
  //     }
  //   }

  //   if (code_id[i_code] == 1) {
  //     if (exchDirection[0] == CWP_FIELD_EXCH_SEND) {

  //       CWP_client_Field_issend(code_names[i_code], cpl_name, field_name);

  //       printf("%d : %s --- Sent field\n", rank, code_names[i_code]);
  //     }
  //     else {

  //       CWP_client_Field_irecv(code_names[i_code], cpl_name, field_name);

  //       printf("%d : %s --- Received field\n", rank, code_names[i_code]);
  //     }
  //   }
  // }

  // for (int i_code = 0 ; i_code < n_code ; i_code++) {
  //   if (code_id[i_code] == 2) {
  //     if (exchDirection[1] == CWP_FIELD_EXCH_RECV) {

  //       CWP_client_Field_wait_irecv(code_names[i_code], cpl_name, field_name);

  //       printf("%d : %s --- wait Received field\n", rank, code_names[i_code]);
  //     }
  //     else {

  //       CWP_client_Field_wait_issend(code_names[i_code], cpl_name, field_name);

  //       printf("%d : %s --- wait Sent field\n", rank, code_names[i_code]);
  //     }
  //   }
  //   if (code_id[i_code] == 1) {
  //     if (exchDirection[0] == CWP_FIELD_EXCH_RECV) {

  //       CWP_client_Field_wait_irecv(code_names[i_code], cpl_name, field_name);

  //       printf("%d : %s --- wait Received field\n", rank, code_names[i_code]);
  //     }
  //     else {

  //       CWP_client_Field_wait_issend(code_names[i_code], cpl_name, field_name);

  //       printf("%d : %s --- wait Sent field\n", rank, code_names[i_code]);
  //     }
  //   }
  // }

  // for (int i_code = 0 ; i_code < n_code ; ++i_code) {
  //   if (code_id[i_code] == 2 && (exchDirection[1] == CWP_FIELD_EXCH_RECV)) {
  //     for (int i = 0 ; i < n_computed_tgts ; i++) {
  //       printf("%12.5e %12.5e\n", recv_values[i_code][3 * i], coord[i_code][0][3 * (computed_tgts[i] - 1)]);
  //     }
  //   }
  //   if (code_id[i_code] == 1 && (exchDirection[0] == CWP_FIELD_EXCH_RECV)) {
  //     for (int i = 0 ; i < n_computed_tgts ; i++) {
  //       printf("%12.5e %12.5e\n", recv_values[i_code][3 * i], coord[i_code][0][3 * (computed_tgts[i] - 1)]);
  //     }
  //   }
  // }

  // Delete interf
  for (int i_code = 0 ; i_code < n_code ; i_code++) {

    CWP_client_Mesh_interf_del(code_names[i_code], cpl_name);

    printf("%d : %s --- Interface deleted\n", rank, code_names[i_code]);
  }

  // Delete coupling
  for (int i_code = 0 ; i_code < n_code ; i_code++) {

    CWP_client_Cpl_del(code_names[i_code], cpl_name);

    printf("%d : %s --- Coupling deleted\n", rank, code_names[i_code]);
  }

  free(element_type);
  free(element_type_cwp);

  //Finalize cwipi

  CWP_client_Finalize();
  printf("%d --- CWIPI finalized\n", rank);

  PDM_MPI_Finalize();

  return 0;
}
