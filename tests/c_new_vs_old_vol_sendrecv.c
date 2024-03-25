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

#define ABS(a)   ((a) <  0  ? -(a) : (a))
#define MIN(a,b) ((a) < (b) ?  (a) : (b))
#define MAX(a,b) ((a) > (b) ?  (a) : (b))

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
         "  -n           <> Number of vertices in band length.\n\n"
         "  -no_random      Disable mesh randomization\n\n"
         "  -n_proc_data <> Number of processes where there are data \n\n"
         "  -h              this message.\n\n");

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
_read_args
(
  int                   argc, 
  char                **argv, 
  CWP_Version_t        *version, 
  int                  *n_vtx_seg1, 
  int                  *n_vtx_seg2,
  double               *length,
  double               *separation_x,
  double               *separation_y,
  double               *separation_z,
  int                  *deform, 
  double               *tolerance,
  int                  *randomize, 
  int                  *nProcData,
  CWP_Spatial_interp_t *loc_method, 
  char                **output_filename, 
  int                  *verbose,
  int                  *extension_depth_tgt,
  int                  *extension_depth_src,
  CWPT_Mesh_nodal_elt_t *elt_type
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
    else if (strcmp(argv[i], "-length") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *length = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-sep") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepy") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_y = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_z = atof(argv[i]);
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
    else if (strcmp(argv[i], "-def") == 0) {
      *deform = 1;
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
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
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
    else if (strcmp(argv[i], "-ext_depth") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_tgt = atoi(argv[i]);
        *extension_depth_src = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ext_depth_tgt") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_tgt = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ext_depth_src") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_src = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-elt_type") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = (CWPT_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static double R[3][3] =
{
  {-0.14547275709949994,  0.8415293589391187 , -0.5202557207618055 },
  { 0.9893622576902102 ,  0.12373586628506748, -0.07649678720582984},
  { 0.                 , -0.5258495730132333 , -0.8505775840931856 }
};

static void
_rotate(const int n_pts, double *coord) {
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
_unrotate(const int n_pts, double *coord) {
  for (int i = 0 ; i < n_pts ; i++) {
    double x = coord[3 * i];
    double y = coord[3 * i + 1];
    double z = coord[3 * i + 2];

    for (int j = 0 ; j < 3 ; j++) {
      coord[3 * i + j] = R[0][j] * x + R[1][j] * y + R[2][j] * z;
    }
  }
}



static void
_cube_mesh
(
 const int                     activeRank,  
 const MPI_Comm                comm,
 const int                     n_part,
 const CWPT_split_dual_t       part_method,
 const CWP_g_num_t             n_vtx_seg,
 const double                  xmin,
 const double                  ymin,
 const double                  zmin,
 const double                  length,
 const int                     deform,
 const int                     part_extension_depth,
 const CWPT_Mesh_nodal_elt_t   elt_type,
 int                         **pn_cell,
 int                         **pn_face,
 int                         **pn_vtx,
 int                        ***pcell_face_idx,
 int                        ***pcell_face,
 int                        ***pcell_vtx,
 int                        ***pface_vtx_idx,
 int                        ***pface_vtx,
 double                     ***pvtx_coord,
 CWP_g_num_t                ***pcell_ln_to_gn,
 CWP_g_num_t                ***pface_ln_to_gn,
 CWP_g_num_t                ***pvtx_ln_to_gn
)
{
  if (activeRank) {

    // Unused variables
    int          *unused_n_edge                = NULL;
    int         **unused_edge_vtx              = NULL;
    int         **unused_face_edge             = NULL;
    CWP_g_num_t **unused_edge_ln_to_gn         = NULL;
    int          *unused_n_surface             = NULL;
    int         **unused_surface_face_idx      = NULL;
    int         **unused_surface_face          = NULL;
    CWP_g_num_t **unused_surface_face_ln_to_gn = NULL;
    int          *unused_n_ridge               = NULL;
    int         **unused_ridge_edge_idx        = NULL;
    int         **unused_ridge_edge            = NULL;
    CWP_g_num_t **unused_ridge_edge_ln_to_gn   = NULL;

    CWPT_generate_mesh_parallelepiped_ngon(comm,
                                           elt_type,
                                           1,
                                           NULL,
                                           xmin,
                                           ymin,
                                           zmin,
                                           length,
                                           length,
                                           length,
                                           n_vtx_seg,
                                           n_vtx_seg,
                                           n_vtx_seg,
                                           n_part,
                                           part_method,
                                           pn_vtx,
                                           &unused_n_edge,
                                           pn_face,
                                           pn_cell,
                                           pvtx_coord,
                                           &unused_edge_vtx,
                                           pface_vtx_idx,
                                           &unused_face_edge,
                                           pface_vtx,
                                           pcell_face_idx,
                                           pcell_face,
                                           pvtx_ln_to_gn,
                                           &unused_edge_ln_to_gn,
                                           pface_ln_to_gn,
                                           pcell_ln_to_gn,
                                           &unused_n_surface,
                                           &unused_surface_face_idx,
                                           &unused_surface_face,
                                           &unused_surface_face_ln_to_gn,
                                           &unused_n_ridge,
                                           &unused_ridge_edge_idx,
                                           &unused_ridge_edge,
                                           &unused_ridge_edge_ln_to_gn);

    // Free unused variables
    for (int i_part = 0; i_part < n_part; i_part++) {
      free(unused_edge_vtx             [i_part]);
      free(unused_face_edge            [i_part]);
      free(unused_edge_ln_to_gn        [i_part]);
      free(unused_surface_face_idx     [i_part]);
      free(unused_surface_face         [i_part]);
      free(unused_surface_face_ln_to_gn[i_part]);
      free(unused_ridge_edge_idx       [i_part]);
      free(unused_ridge_edge           [i_part]);
      free(unused_ridge_edge_ln_to_gn  [i_part]);
    }
    free(unused_n_edge               );
    free(unused_edge_vtx             );
    free(unused_face_edge            );
    free(unused_edge_ln_to_gn        );
    free(unused_n_surface            );
    free(unused_surface_face_idx     );
    free(unused_surface_face         );
    free(unused_surface_face_ln_to_gn);
    free(unused_n_ridge              );
    free(unused_ridge_edge_idx       );
    free(unused_ridge_edge           );
    free(unused_ridge_edge_ln_to_gn  );

    // TO DO : part_extension_depth > 0 => create extension
  }
  else {
    *pn_cell = (int *) malloc (sizeof(int *) * n_part);
    *pn_face = (int *) malloc (sizeof(int *) * n_part);
    *pn_vtx = (int *) malloc (sizeof(int *) * n_part);
    *pcell_face_idx = (int **) malloc (sizeof(int *) * n_part);
    *pcell_face = (int **) malloc (sizeof(int *) * n_part);
    *pcell_vtx = (int **) malloc (sizeof(int *) * n_part);
    *pface_vtx_idx = (int **) malloc (sizeof(int *) * n_part);
    *pface_vtx = (int **) malloc (sizeof(int *) * n_part);
    *pvtx_coord = (double **) malloc (sizeof(int *) * n_part);
    *pcell_ln_to_gn = (CWP_g_num_t **) malloc (sizeof(int *) * n_part);
    *pface_ln_to_gn = (CWP_g_num_t **) malloc (sizeof(int *) * n_part);
    *pvtx_ln_to_gn = (CWP_g_num_t **) malloc (sizeof(int *) * n_part);
    for (int ipart = 0 ; ipart < n_part ; ipart++) {
      int _nFace = 0;
      int _nCell = 0;
      int _nVtx = 0;
      (*pn_cell)[ipart] = _nCell;
      (*pn_face)[ipart] = _nFace;
      (*pn_vtx)[ipart] = _nVtx;
      (*pcell_face_idx)[ipart] = (int *) malloc(sizeof(int) * (_nCell + 1));
      (*pcell_face_idx)[ipart][0] = 0;
      (*pcell_face)[ipart] = (int *) malloc(sizeof(int) * 0); 
      (*pcell_vtx)[ipart] = (int *) malloc(sizeof(int) * 0);
      (*pface_vtx_idx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
      (*pface_vtx_idx)[ipart][0] = 0;
      (*pface_vtx)[ipart] = (int *) malloc(sizeof(int) * 0);
      (*pvtx_coord)[ipart] = (double *) malloc(sizeof(double) * 0);
      (*pcell_ln_to_gn)[ipart] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * _nCell);
      (*pface_ln_to_gn)[ipart] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * _nFace);
      (*pvtx_ln_to_gn)[ipart] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * _nVtx);
    }
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

  return current_rank_has_mesh;
}


/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {
  // Read args from command line
  CWP_Version_t version           = CWP_VERSION_NEW;
  int n_vtx_seg1                  = 10;
  int n_vtx_seg2                  = 10;
  int randomize                   = 1;
  int n_proc_data                 = -1;

  CWP_Spatial_interp_t loc_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;
  // CWP_Spatial_interp_t loc_method = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  int verbose                     = 0;

  double length                   = 20.;
  int deform                      = 0;

  double      separation_x        = 2.;
  double      separation_y        = 0.;
  double      separation_z        = 0.;

  double      tolerance           = 1e-2;

  char* output_filename           = NULL;

  int         extension_depth_tgt = 0;
  int         extension_depth_src = 0;

  CWPT_Mesh_nodal_elt_t elt_type = CWPT_MESH_NODAL_TETRA4;

  _read_args (argc, 
              argv, 
             &version, 
             &n_vtx_seg1, 
             &n_vtx_seg2,
             &length,
             &separation_x,
             &separation_y,
             &separation_z,
             &deform, 
             &tolerance,
             &randomize, 
             &n_proc_data,
             &loc_method, 
             &output_filename, 
             &verbose,
             &extension_depth_tgt,
             &extension_depth_src,
             &elt_type);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  assert (comm_world_size > 1);

  if (n_proc_data == 1) {
    n_proc_data = 2;
  }

  // Initialize CWIPI
  int n_part = 1;
  int n_code = 1;
  int code_id;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  int n_vtx_seg;
  if (rank < comm_world_size / 2) {
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

  if (verbose && rank == 0) {
    printf("CWIPI Init OK\n");
  }

  // Create coupling
  const char *coupling_name = "c_new_vs_old_sendrecv";

  if (version == CWP_VERSION_OLD) {
    cwipi_create_coupling(coupling_name,                             // Coupling id
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                          coupled_code_name[0],                      // Coupled application id
                          3,                                         // Geometric entities dimension
                          tolerance,                                 // Geometric tolerance
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
                   CWP_INTERFACE_VOLUME,
                   CWP_COMM_PAR_WITH_PART,
                   loc_method,
                   n_part,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

    CWP_Visu_set (code_name[0],
                  coupling_name,
                  1,
                  CWP_VISU_FORMAT_ENSIGHT,
                  "text");

  }

  if (verbose && rank == 0) {
    printf("Create coupling OK\n");
  }

  // Define mesh

  int _n_proc_data = n_proc_data;
  if (n_proc_data > 0) {
    if (code_id == 1) {
      _n_proc_data /= 2;
    }
    else {
      _n_proc_data -= n_proc_data / 2;
    }
  }
  MPI_Comm mesh_comm;
  int current_rank_has_mesh = _set_rank_has_mesh(intra_comm[0], _n_proc_data, &mesh_comm);

  int true_n_proc_data;
  MPI_Reduce(&current_rank_has_mesh, &true_n_proc_data, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    printf("nb procs with mesh data = %d\n", true_n_proc_data);
  }

  double xmin = -0.5 * length;
  double ymin = -0.5 * length;
  double zmin = -0.5 * length;

  double xyz_min[3] = {xmin, ymin, zmin};
  double xyz_max[3] = {xmin + length, ymin + length, zmin + length};

//  int init_random = (int) time(NULL);
  int init_random = 5;

  MPI_Comm code_mesh_comm;
  MPI_Comm_split(mesh_comm, code_id, rank, &code_mesh_comm);

  if (code_id == 2) {
    init_random++;
    xmin += separation_x;
    ymin += separation_y;
    zmin += separation_z;
  }

  int                         *pn_cell;
  int                         *pn_face;
  int                         *pn_vtx;
  int                        **pcell_face_idx;
  int                        **pcell_face;
  int                        **pcell_vtx;
  int                        **pface_vtx_idx;
  int                        **pface_vtx;
  double                     **pvtx_coord;
  CWP_g_num_t                **pcell_ln_to_gn;
  CWP_g_num_t                **pface_ln_to_gn;
  CWP_g_num_t                **pvtx_ln_to_gn;

  _cube_mesh (current_rank_has_mesh,
              code_mesh_comm,
              1,
              CWPT_SPLIT_DUAL_WITH_HILBERT,
              n_vtx_seg,
              xmin,
              ymin,
              zmin,
              length,
              deform,
              extension_depth_src,
              elt_type,
             &pn_cell,
             &pn_face,
             &pn_vtx,
             &pcell_face_idx,
             &pcell_face,
             &pcell_vtx,
             &pface_vtx_idx,
             &pface_vtx,
             &pvtx_coord,
             &pcell_ln_to_gn,
             &pface_ln_to_gn,
             &pvtx_ln_to_gn);

  int *cellVtxIdx = malloc(sizeof(int) * (pn_cell[0] + 1));
  cellVtxIdx[0] = 0;

  CWP_Block_t block_type;

  if (elt_type == CWPT_MESH_NODAL_TETRA4) {
    block_type = CWP_BLOCK_CELL_TETRA4;
    for (int i = 0; i < pn_cell[0]; i++) {
      cellVtxIdx[i+1] = 4 + cellVtxIdx[i];  
    }
  }

  else if (elt_type == CWPT_MESH_NODAL_HEXA8) {
    block_type = CWP_BLOCK_CELL_HEXA8;
    for (int i = 0; i < pn_cell[0]; i++) {
      cellVtxIdx[i+1] = 8 + cellVtxIdx[i];  
    }
  }

  else if (elt_type == CWPT_MESH_NODAL_PYRAMID5) {
    block_type = CWP_BLOCK_CELL_PYRAM5;
    for (int i = 0; i < pn_cell[0]; i++) {
      cellVtxIdx[i+1] = 5 + cellVtxIdx[i];  
    }
  }

  else if (elt_type == CWPT_MESH_NODAL_PRISM6) {
    block_type = CWP_BLOCK_CELL_PRISM6;
    for (int i = 0; i < pn_cell[0]; i++) {
      cellVtxIdx[i+1] = 6 + cellVtxIdx[i];  
    }
  }

  // Set interface mesh
  if (version == CWP_VERSION_OLD) {
    cwipi_define_mesh(coupling_name, pn_vtx[0], pn_cell[0], pvtx_coord[0], cellVtxIdx, pcell_vtx[0]);
  }
  else {
    CWP_Mesh_interf_vtx_set(code_name[0], coupling_name, 0, pn_vtx[0], pvtx_coord[0], pvtx_ln_to_gn[0]);

    int block_id = CWP_Mesh_interf_block_add (code_name[0],
                                              coupling_name,
                                              block_type);  


    CWP_Mesh_interf_block_std_set (code_name[0],
                                   coupling_name,
                                   0,
                                   block_id,
                                   pn_cell[0],
                                   pcell_vtx[0],
                                   pcell_ln_to_gn[0]); 


    CWP_Mesh_interf_finalize(code_name[0], coupling_name);
  }

  if (verbose && rank == 0) {
    printf("Set mesh OK\n");
  }

  // Create and set fields
  double *send_val = (double *) malloc(sizeof(double) * pn_vtx[0]);
  double *recv_val = (double *) malloc(sizeof(double) * pn_vtx[0]);

  const char *field_name  = "cooX";

  for (int i = 0 ; i < pn_vtx[0] ; i++) {
    send_val[i] = pvtx_coord[0][3 * i];
  }

  if (version == CWP_VERSION_NEW) {
    CWP_Status_t visu_status = CWP_STATUS_ON;
    MPI_Barrier(MPI_COMM_WORLD);

    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SENDRECV,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);
  }

  if (verbose && rank == 0) {
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
      located = cwipi_get_located_points(coupling_name);
    }

  }
  else {
    // PDM_part_to_block_global_statistic_reset();
    // PDM_block_to_part_global_statistic_reset();

    char char_tol[99];
    sprintf(char_tol, "%e", tolerance);
    CWP_Spatial_interp_property_set(code_name[0], coupling_name, "tolerance", CWP_DOUBLE, char_tol);
    CWP_Spatial_interp_weights_compute(code_name[0], coupling_name);
  
    // n_unlocated = CWP_N_uncomputed_tgts_get(code_name[0], coupling_name, field_name, 0);
    n_located = CWP_N_computed_tgts_get(code_name[0], coupling_name, field_name, 0);
    located =  CWP_Computed_tgts_get(code_name[0], coupling_name, field_name, 0);
  }

  CWP_timer_hang_on(timer);
  t_end = CWP_timer_elapsed(timer);
  CWP_timer_resume(timer);
  double geom_time = t_end - t_start;
  double max_geom_time;
  MPI_Reduce(&geom_time, &max_geom_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("\n\nGeometric algorithm :%12.5es\n", max_geom_time);
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
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val);
      CWP_Field_issend(code_name[0], coupling_name, field_name);

      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val);
      CWP_Field_irecv(code_name[0], coupling_name, field_name);
    }
    else {
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val);
      CWP_Field_irecv(code_name[0], coupling_name, field_name);

      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val);
      CWP_Field_issend(code_name[0], coupling_name, field_name);
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
      CWP_Field_wait_irecv(code_name[0], coupling_name, field_name);
    }
    else {
      CWP_Field_wait_irecv(code_name[0], coupling_name, field_name);
      CWP_Field_wait_issend(code_name[0], coupling_name, field_name);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  CWP_timer_hang_on(timer);
  t_end = CWP_timer_elapsed(timer);
  double exch_time = t_end - t_start;
  double max_exch_time;
  MPI_Reduce(&exch_time, &max_exch_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("Exchange 1 fields issend/irecv: %12.5es\n", max_exch_time1);
    printf("Exchange 1 fields wait        : %12.5es\n", max_exch_time);
    printf("Total Exchange 1              : %12.5es\n", max_exch_time1 + max_exch_time);
  }

  // double redondance_geom = max_exch_time1;
  max_geom_time += max_exch_time1;

  //  Check
  if (1) {
    double max_err = 0.;
    CWP_g_num_t n_wrong = 0;
    if (code_id == 2) {
      for (int i = 0 ; i < n_located ; i++) {

        int pt_id = located[i] -1;
        double coord[3] = {pvtx_coord[0][3*pt_id], pvtx_coord[0][3*pt_id+1], pvtx_coord[0][3*pt_id+2]};
        if (loc_method == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE ||
            loc_method == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE) {
          if (deform) {
            _unrotate(1, coord);
          }
          for (int j = 0; j < 3; j++) {
            coord[j] = MIN(MAX(coord[j], xyz_min[j]), xyz_max[j]);
          }
          if (deform) {
            _rotate(1, coord);
          }
        }

        double err = ABS (recv_val[i] - coord[0]);
        if (err > 1.e-4) {
          n_wrong++;
          printf("[%d] !! vtx %ld %d err = %g (x = %f, recv = %f)\n",
          rank, pvtx_ln_to_gn[0][(located[i] - 1)], located[i], err, pvtx_coord[0][3*(located[i]-1)], recv_val[i]);
        }
        if (err > max_err) {
          max_err = err;
        }
      }
    }

    double global_max_err = 0.;
    MPI_Reduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    CWP_g_num_t global_n_wrong = 0;
    MPI_Reduce(&n_wrong, &global_n_wrong, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      printf("Max error = %g (%ld wrong points)\n", global_max_err, global_n_wrong);
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
    CWP_Time_step_end(code_name[0]);
    CWP_Cpl_del(code_name[0], coupling_name);
  }
  // Free memory
  free(code_name);
  free(coupled_code_name);
  free(intra_comm);
  
  if (current_rank_has_mesh) {
    for (int ipart = 0; ipart < n_part; ipart++) {
      free (pcell_face_idx[ipart]);
      free (pcell_face[ipart]);
      free (pcell_vtx[ipart]);
      free (pface_vtx_idx[ipart]);
      free (pface_vtx[ipart]);
      free (pvtx_coord[ipart]);
      free (pcell_ln_to_gn[ipart]);
      free (pface_ln_to_gn[ipart]);
      free (pvtx_ln_to_gn[ipart]);
    }
  }

  free(pn_vtx);
  free(pn_cell);
  free(pn_face);
  free(pvtx_coord);
  free(pvtx_ln_to_gn);
  free(pcell_face_idx);
  free(pcell_face);
  free(pcell_vtx);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pcell_ln_to_gn);
  free(pface_ln_to_gn);
  free(cellVtxIdx);
  free(send_val);
  free(recv_val);

  CWP_timer_free(timer);

  //  Finalize CWIPI
  if (version == CWP_VERSION_OLD) {
    cwipi_finalize();
  }
  else {
    CWP_Finalize();
  }

  // if (version == CWP_VERSION_NEW) {
  //   double min_elaps_create_ptb;
  //   double max_elaps_create_ptb;
  //   double min_cpu_create_ptb;
  //   double max_cpu_create_ptb;
  //   double min_elaps_create2_ptb;
  //   double max_elaps_create2_ptb;
  //   double min_cpu_create2_ptb;
  //   double max_cpu_create2_ptb;
  //   double min_elaps_exch_ptb;
  //   double max_elaps_exch_ptb;
  //   double min_cpu_exch_ptb;
  //   double max_cpu_exch_ptb;

  //   PDM_part_to_block_global_timer_get(PDM_MPI_COMM_WORLD,
  //                                      &min_elaps_create_ptb,
  //                                      &max_elaps_create_ptb,
  //                                      &min_cpu_create_ptb,
  //                                      &max_cpu_create_ptb,
  //                                      &min_elaps_create2_ptb,
  //                                      &max_elaps_create2_ptb,
  //                                      &min_cpu_create2_ptb,
  //                                      &max_cpu_create2_ptb,
  //                                      &min_elaps_exch_ptb,
  //                                      &max_elaps_exch_ptb,
  //                                      &min_cpu_exch_ptb,
  //                                      &max_cpu_exch_ptb);

  //   double min_elaps_create_btp;
  //   double max_elaps_create_btp;
  //   double min_cpu_create_btp;
  //   double max_cpu_create_btp;
  //   double min_elaps_exch_btp;
  //   double max_elaps_exch_btp;
  //   double min_cpu_exch_btp;
  //   double max_cpu_exch_btp;

  //   PDM_block_to_part_global_timer_get(PDM_MPI_COMM_WORLD,
  //                                      &min_elaps_create_btp,
  //                                      &max_elaps_create_btp,
  //                                      &min_cpu_create_btp,
  //                                      &max_cpu_create_btp,
  //                                      &min_elaps_exch_btp,
  //                                      &max_elaps_exch_btp,
  //                                      &min_cpu_exch_btp,
  //                                      &max_cpu_exch_btp);

  //   if (rank == 0) {
  //     printf("Global time in PDM_part_to_block : \n");
  //     printf("   - ptb min max elaps create  : %12.5e %12.5e\n",
  //            min_elaps_create_ptb,
  //            max_elaps_create_ptb);
  //     printf("   - ptb min max elaps create2 : %12.5e %12.5e\n",
  //            min_elaps_create2_ptb,
  //            max_elaps_create2_ptb);
  //     printf("   - ptb min max elaps exch    : %12.5e %12.5e\n",
  //            min_elaps_exch_ptb,
  //            max_elaps_exch_ptb);
  //     fflush(stdout);

  //     printf("Global time in PDM_block_to_part : \n");
  //     printf("   - btp min max elaps create  : %12.5e %12.5e\n",
  //            min_elaps_create_btp,
  //            max_elaps_create_btp);
  //     printf("   - btp min max elaps exch    : %12.5e %12.5e\n",
  //            min_elaps_exch_btp,
  //            max_elaps_exch_btp);
  //     fflush(stdout);
  //   }
  // }

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
