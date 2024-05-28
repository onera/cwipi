/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

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
#include <assert.h>
#include <string.h>
#include <time.h>

#include "cwipi.h"
#include "cwp.h"
#include "cwp_priv.h"

#define ABS(a)    ((a) <  0  ? -(a) : (a))
#define MAX(a, b) ((a) > (b) ?  (a) : (b))

typedef enum {
  CWP_VERSION_OLD,
  CWP_VERSION_NEW
} CWP_Version_t;


/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -n           <> Number of vertices in band width.\n\n"
         "  -no_random      Disable mesh randomization\n\n"
         "  -n_proc_data <> Number of processes where there are data \n\n"
         "  -h              this message.\n\n");

  exit(exit_code);
}


/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 *---------------------------------------------------------------------*/

static void
_read_args
(
 int                    argc,
 char                 **argv,
 CWP_Version_t         *version,
 int                   *n_vtx_seg1,
 int                   *n_vtx_seg2,
 double                *width,
 double                *depth,
 int                   *rotation,
 int                   *randomize,
 int                   *nProcData,
 CWP_Spatial_interp_t  *loc_method,
 char                 **output_filename,
 int                   *verbose
 )
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-new") == 0) {
      *version = CWP_VERSION_NEW;
    }
    else if (strcmp(argv[i], "-old") == 0) {
      *version = CWP_VERSION_OLD;
    }
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg1 = atoi(argv[i]);
        *n_vtx_seg2 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg1 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg2 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-width") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *width = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-depth") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *depth = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-output") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *output_filename = argv[i];
      }
    }
    else if (strcmp(argv[i], "-rot") == 0) {
      *rotation = 1;
    }
    else if (strcmp(argv[i], "-no_random") == 0) {
      *randomize = 0;
    }
    else if (strcmp(argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *nProcData = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }

    i++;
  }
}


static int
_set_rank_has_mesh(const MPI_Comm comm, const int nProcData, MPI_Comm *meshComm) {
  int current_rank_has_mesh = 1;
  int rank;
  int commSize;

  assert (comm != MPI_COMM_NULL);

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &commSize);

  if (nProcData > 0 && nProcData < commSize) {

    MPI_Comm comm_shared;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &comm_shared);

    int rankInNode = 0;
    MPI_Comm_rank(comm_shared, &rankInNode);

    int nNode = 0;
    int iNode = -1;
    int masterRank = (rankInNode == 0);

    int *rankInNodes = malloc(sizeof(int) * commSize);

    MPI_Allreduce(&masterRank, &nNode, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allgather(&rankInNode, 1, MPI_INT, rankInNodes, 1, MPI_INT, comm);

    current_rank_has_mesh = 0;

    for (int i = 0 ; i < rank ; i++) {
      if (rankInNodes[i] == 0) {
        iNode += 1;
      }
    }

    if (nProcData <= nNode) {
      if (iNode < nProcData && masterRank) {
        current_rank_has_mesh = 1;
      }
    }
    else {
      if (rankInNode < (nProcData / nNode)) {
        current_rank_has_mesh = 1;
      }
      if ((rankInNode == (nProcData / nNode)) && (iNode < (nProcData % nNode))) {
        current_rank_has_mesh = 1;
      }
    }

    MPI_Comm_split(comm, current_rank_has_mesh, rank, meshComm);
    free(rankInNodes);
  }
  else {
    MPI_Comm_dup(comm, meshComm);
  }

  return current_rank_has_mesh;
}


static void
_add_depth(const int n_pts, const double width, const double depth, double *coord) {
  for (int i = 0 ; i < n_pts ; i++) {
    double x = 2. * coord[3 * i] / width;
    double y = 2. * coord[3 * i + 1] / width;
    coord[3 * i + 2] = 0.5 * depth * (1. - (x * x + y * y));
  }
}


static void
_rotate(const int n_pts, double *coord) {
  double R[3][3] = {{0.9362934,  -0.2896295, 0.1986693},
                    {0.3129918,  0.9447025,  -0.0978434},
                    {-0.1593451, 0.1537920,  0.9751703}};

  for (int i = 0 ; i < n_pts ; i++) {
    double x = coord[3 * i];
    double y = coord[3 * i + 1];
    double z = coord[3 * i + 2];

    for (int j = 0 ; j < 3 ; j++) {
      coord[3 * i + j] = R[j][0] * x + R[j][1] * y + R[j][2] * z;
    }
  }
}


static void
_gen_mesh
(
 const int            activeRank,
 const MPI_Comm       comm,
 const CWP_g_num_t    n_vtx_seg,
 const double         length,
 const double         depth,
 const int            rotation,
 int                 *n_vtx,
 int                 *n_elt,
 double            **vtx_coord,
 int               **elt_vtx_idx,
 int               **elt_vtx
)
{
  if (activeRank) {
    CWPT_generate_mesh_rectangle_simplified(comm,
                                            n_vtx_seg,
                                            n_vtx,
                                            n_elt,
                                            vtx_coord,
                                            elt_vtx_idx,
                                            elt_vtx);
  }
  else {
    *n_vtx = 0;
    *n_elt = 0;
    *vtx_coord   = malloc(sizeof(double) * 0);
    *elt_vtx_idx = malloc(sizeof(int   ) * 1);
    *elt_vtx     = malloc(sizeof(int   ) * 0);
  }


  for (int i = 0; i < *n_vtx; i++) {
    (*vtx_coord)[3*i  ] = length * ((*vtx_coord)[3*i  ]/10. - 0.5);
    (*vtx_coord)[3*i+1] = length * ((*vtx_coord)[3*i+1]/5.  - 0.5);
  }

  _add_depth(*n_vtx, length, depth, *vtx_coord);

  if (rotation) {
    _rotate(*n_vtx, *vtx_coord);
  }
}

/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1_P1
 *
 *---------------------------------------------------------------------*/

int
main
(
 int   argc,
 char *argv[]
 )
{
  // Read args from command line
  CWP_Version_t         version         = CWP_VERSION_NEW;
  int                   n_vtx_seg1      = 10;
  int                   n_vtx_seg2      = 11;
  int                   randomize       = 1;
  int                   n_proc_data     = -1;
  CWP_Spatial_interp_t  loc_method      = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  int                   verbose         = 0;
  double                width           = 20.;
  double                depth           = 1.;
  int                   rotation        = 0;
  char                 *output_filename = NULL;

  _read_args(argc,
             argv,
             &version,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &width,
             &depth,
             &rotation,
             &randomize,
             &n_proc_data,
             &loc_method,
             &output_filename,
             &verbose);

  int filedump = 0;
  if (output_filename != NULL) {
    filedump = 1;
  }

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  assert(n_rank > 1);

  if (n_proc_data == 1) {
    n_proc_data = 2;
  }

  // Initialize CWIPI
  int n_part = 1;
  int n_code = 1;
  int code_id;
  const char **code_name         = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  int n_vtx_seg;
  if (i_rank < n_rank / 2) {
    code_id = 1;
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
    n_vtx_seg = n_vtx_seg1;
  }
  else {
    code_id = 2;
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
    n_vtx_seg = n_vtx_seg2;
  }

  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
  if (version == CWP_VERSION_OLD) {
    cwipi_init(MPI_COMM_WORLD, code_name[0], intra_comm);
  }

  else {
    CWP_Init(MPI_COMM_WORLD,
             n_code,
             (const char **) code_name,
             is_active_rank,
             intra_comm);
  }

  if (verbose && i_rank == 0) {
    printf("CWIPI Init OK\n");
  }


  // Create coupling
  const char *coupling_name = "c_new_to_old_api_surf_P1P1";

  if (version == CWP_VERSION_OLD) {
    cwipi_create_coupling(coupling_name,                             // Coupling id
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                          coupled_code_name[0],                      // Coupled application id
                          2,                                         // Geometric entities dimension
                          0.01,                                      // Geometric tolerance
                          CWIPI_STATIC_MESH,                         // Mesh type
                          CWIPI_SOLVER_CELL_VERTEX,                  // Solver type
                          -1,                                        // Postprocessing frequency
                          "EnSight Gold",                            // Postprocessing format
                          "text");
  }
  else {
    CWP_Cpl_create(code_name[0],
                   coupling_name,
                   coupled_code_name[0],
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   loc_method,
                   n_part,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

    // CWP_Visu_set (code_name[0],
    //               coupling_name,
    //               -1,
    //               CWP_VISU_FORMAT_ENSIGHT,
    //               "text");

  }

  if (verbose && i_rank == 0) {
    printf("Create coupling OK\n");
  }

  // Define mesh
  MPI_Comm mesh_comm;
  int _n_proc_data = n_proc_data;
  if (n_proc_data > 0) {
    if (code_id == 1) {
      _n_proc_data /= 2;
    }
    else {
      _n_proc_data -= n_proc_data / 2;
    }
  }
  int current_rank_has_mesh = _set_rank_has_mesh(intra_comm[0], _n_proc_data, &mesh_comm);

  int true_n_proc_data;
  MPI_Reduce(&current_rank_has_mesh, &true_n_proc_data, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("nb procs with mesh data = %d\n", true_n_proc_data);
  }

  int     n_vtx;
  int     n_elt;
  double *vtx_coord;
  int    *elt_vtx_idx;
  int    *elt_vtx;
  _gen_mesh(current_rank_has_mesh,
            mesh_comm,
            n_vtx_seg,
            width,
            depth,
            rotation,
            &n_vtx,
            &n_elt,
            &vtx_coord,
            &elt_vtx_idx,
            &elt_vtx);



  // Set interface mesh
  if (version == CWP_VERSION_OLD) {
    cwipi_define_mesh(coupling_name,
                      n_vtx,
                      n_elt,
                      vtx_coord,
                      elt_vtx_idx,
                      elt_vtx);
  }
  else {
    CWP_Mesh_interf_vtx_set(code_name[0], coupling_name, 0, n_vtx, vtx_coord, NULL);

    // CWP_Mesh_interf_from_faceedge_set(code_name[0],
    //                                   coupling_name,
    //                                   0,
    //                                   nFace[0],
    //                                   faceEdgeIdx[0],
    //                                   faceEdge[0],
    //                                   nEdge[0],
    //                                   edgeVtx[0],
    //                                   faceLNToGN[0]);
    CWP_Mesh_interf_from_facevtx_set(code_name[0],
                                     coupling_name,
                                     0,
                                     n_elt,
                                     elt_vtx_idx,
                                     elt_vtx,
                                     NULL);

    CWP_Mesh_interf_finalize(code_name[0], coupling_name);
  }

  if (verbose && i_rank == 0) {
    printf("Set mesh OK\n");
  }

  // Create and set fields
  double *send_val = NULL;
  double *recv_val = NULL;

  const char *field_name  = "cooX";
  const char *field_name2 = "coocooY";

  if (code_id == 1) {
    send_val = (double *) malloc(sizeof(double) * n_vtx);

    for (int i = 0 ; i < n_vtx ; i++) {
      send_val[i] = vtx_coord[3 * i];
    }
  }
  else {
    recv_val = (double *) malloc(sizeof(double) * n_vtx);
  }

  if (version == CWP_VERSION_NEW) {
    CWP_Status_t visu_status = CWP_STATUS_OFF;
    MPI_Barrier(MPI_COMM_WORLD);

    if (code_id == 1) {
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val);
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val);
    }
    else {
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val);
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val);
    }
  }

  if (verbose && i_rank == 0) {
    printf("Fields OK\n");
  }

  // Perform geometric algorithm
  CWP_timer_t *timer = CWP_timer_create();
  double t_start, t_end;

  MPI_Barrier(MPI_COMM_WORLD);
  CWP_timer_init(timer);

  t_start = CWP_timer_elapsed(timer);
  CWP_timer_resume(timer);

  // int n_unlocated = 0;
  int n_located = 0;
  const int *located = NULL;
  if (version == CWP_VERSION_OLD) {
    cwipi_locate(coupling_name);

    if (code_id != 1) {
      // n_unlocated = cwipi_get_n_not_located_points(coupling_name);
      n_located = cwipi_get_n_located_points(coupling_name);
      located   = cwipi_get_located_points  (coupling_name);
    }
  }
  else {
    CWP_Spatial_interp_property_set(code_name[0], coupling_name, "tolerance", CWP_DOUBLE, "1e-2");
    CWP_Spatial_interp_weights_compute(code_name[0], coupling_name);

    if (code_id != 1) {
      // n_unlocated = CWP_N_uncomputed_tgts_get(code_name[0], coupling_name, field_name2, 0);
      n_located = CWP_N_computed_tgts_get(code_name[0], coupling_name, field_name2, 0);
      located   = CWP_Computed_tgts_get  (code_name[0], coupling_name, field_name2, 0);
    }
  }

  CWP_timer_hang_on(timer);
  t_end = CWP_timer_elapsed(timer);
  CWP_timer_resume(timer);
  double geom_time = t_end - t_start;
  double max_geom_time;
  MPI_Reduce(&geom_time, &max_geom_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  FILE *output;
  if (filedump) {
    output = fopen(output_filename, "w");
  }

  if (i_rank == 0) {
    printf("\n\nGeometric algorithm :%12.5es\n", max_geom_time);

    if (filedump) {
      fprintf(output, "%d %12.5e ", n_rank, max_geom_time);
    }
  }

  if (verbose && i_rank == 0) {
    printf("Geometric algorithm OK\n");
  }

  //  Exchange interpolated fields 1
  MPI_Barrier(MPI_COMM_WORLD);

  CWP_timer_hang_on(timer);
  t_start = CWP_timer_elapsed(timer);
  CWP_timer_resume(timer);

  int request;
  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_issend(coupling_name, "ech", 0, 1, 1, 0.1, field_name, send_val, &request);
    }
    else {
      cwipi_irecv(coupling_name, "ech", 0, 1, 1, 0.1, field_name, recv_val, &request);
    }
  }

  else {
    if (code_id == 1) {
      CWP_Field_issend(code_name[0], coupling_name, field_name);
    }
    else {
      CWP_Field_irecv(code_name[0], coupling_name, field_name);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  CWP_timer_hang_on(timer);
  t_end = CWP_timer_elapsed(timer);
  double exch_time1 = t_end - t_start;
  double max_exch_time1;
  t_start = t_end;
  MPI_Reduce(&exch_time1, &max_exch_time1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  CWP_timer_resume(timer);

  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_wait_issend(coupling_name, request);
    }
    else {
      cwipi_wait_irecv(coupling_name, request);
    }
  }
  else {
    if (code_id == 1) {
      CWP_Field_wait_issend(code_name[0], coupling_name, field_name);
    }
    else {
      CWP_Field_wait_irecv(code_name[0], coupling_name, field_name);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  CWP_timer_hang_on(timer);
  t_end = CWP_timer_elapsed(timer);
  double exch_time = t_end - t_start;
  double max_exch_time;
  MPI_Reduce(&exch_time, &max_exch_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (i_rank == 0) {
    printf("Exchange 1 fields issend/irecv: %12.5es\n", max_exch_time1);
    printf("Exchange 1 fields wait        : %12.5es\n", max_exch_time);
    printf("Total Exchange 1              : %12.5es\n", max_exch_time1 + max_exch_time);

    if (filedump) {
      fprintf(output, "%12.5e ", max_exch_time1 + max_exch_time);
    }
  }

  double redondance_geom = max_exch_time1;
  max_geom_time += max_exch_time1;

  //  Exchange interpolated fields 2
  CWP_timer_resume(timer);
  MPI_Barrier(MPI_COMM_WORLD);

  CWP_timer_hang_on(timer);
  t_start = CWP_timer_elapsed(timer);
  CWP_timer_resume(timer);

  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_issend(coupling_name, "ech", 0, 1, 1, 0.1, field_name2, send_val, &request);
    }
    else {
      cwipi_irecv(coupling_name, "ech", 0, 1, 1, 0.1, field_name2, recv_val, &request);
    }
  }

  else {
    if (code_id == 1) {
      // CWP_Field_issend(code_name[0], coupling_name, field_name2);
    }
    else {
      // CWP_Field_irecv(code_name[0], coupling_name, field_name2);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  CWP_timer_hang_on(timer);
  t_end = CWP_timer_elapsed(timer);
  exch_time1 = t_end - t_start;
  CWP_UNUSED (max_exch_time1);
  t_start = t_end;
  MPI_Reduce(&exch_time1, &max_exch_time1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  CWP_timer_resume(timer);

  redondance_geom += -max_exch_time1;

  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_wait_issend(coupling_name, request);
    }
    else {
      cwipi_wait_irecv(coupling_name, request);
    }
  }
  else {
    if (code_id == 1) {
      // CWP_Field_wait_issend(code_name[0], coupling_name, field_name2);
    }
    else {
      // CWP_Field_wait_irecv(code_name[0], coupling_name, field_name2);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  CWP_timer_hang_on(timer);
  t_end = CWP_timer_elapsed(timer);
  exch_time = t_end - t_start;
  CWP_UNUSED (max_exch_time1);
  MPI_Reduce(&exch_time, &max_exch_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (i_rank == 0) {
    printf("Exchange 2 fields issend/irecv :%12.5es\n", max_exch_time1);
    printf("Exchange 2  fields wait        :%12.5es\n", max_exch_time);
    printf("Total exchange 2               :%12.5es\n", max_exch_time1 + max_exch_time);

    printf("\n\nTemps geometrie                            : %12.5es\n", max_geom_time);
    printf("Temps geometrie escompte (sans redondance) : %12.5es\n",
           max_geom_time - redondance_geom);
    printf("Temps un Echange aux noeuds                : %12.5es\n",
           max_exch_time1 + max_exch_time);

    if (filedump) {
      fprintf(output, "%12.5e\n", max_exch_time1 + max_exch_time);
      fclose(output);
    }
  }


  //  Check
  int check = EXIT_SUCCESS;
  if (1) {
    double max_err = 0.;
    if (code_id == 2) {
      for (int i = 0 ; i < n_located ; i++) {
        double err = ABS (recv_val[i] - vtx_coord[3 * (located[i] - 1)]);
        // if (err > 1.e-5) {
        //   printf("[%d] !! vtx %ld %d err = %g (x = %f, recv = %f)\n",
        //   i_rank, vtxLNToGN[0][(located[i] - 1)], located[i], err, vtxCoord[0][3*(located[i] - 1)], recv_val[i]);
        // }
        if (err > max_err) {
          max_err = err;
        }
      }
    }

    double global_max_err = 0.;
    MPI_Reduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (i_rank == 0) {
      printf("Max error = %g\n", global_max_err);
    }

    if (global_max_err > MAX(1e-9, 1e-3*depth)) {
      check = EXIT_FAILURE;
    }
  }

  //  Delete interface mesh
  if (version == CWP_VERSION_NEW) {
    CWP_Mesh_interf_del(code_name[0], coupling_name);
  }
  //  Delete coupling
  if (version == CWP_VERSION_OLD) {
    cwipi_delete_coupling(coupling_name);
  }
  else {
    CWP_Cpl_del(code_name[0], coupling_name);
  }

  // Free memory
  free(code_name);
  free(coupled_code_name);
  free(intra_comm);

  free(vtx_coord);
  free(elt_vtx_idx);
  free(elt_vtx);

  if (code_id == 1) {
    free(send_val);
  }
  else {
    free(recv_val);
  }
  CWP_timer_free(timer);

  //  Finalize CWIPI
  if (version == CWP_VERSION_OLD) {
    cwipi_finalize();
  }
  else {
    CWP_Finalize();
  }

  // Finalize MPI
  MPI_Finalize();

  return check;
}
