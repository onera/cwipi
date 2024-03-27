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

#include "grid_mesh.h"



#define ABS(a)   ((a) <  0  ? -(a) : (a))
#define MAX(a,b) ((a) > (b) ?  (a) : (b))

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
  int                    argc,
  char                 **argv,
  int                   *n_vtx_seg1,
  int                   *n_vtx_seg2,
  double                *tolerance,
  int                   *randomize
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
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
    else if (strcmp(argv[i], "-no_random") == 0) {
      *randomize = 0;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_gen_mesh
(
 const MPI_Comm            comm,
 const int                 n_part,
 const CWPT_split_dual_t   part_method,
 const CWP_g_num_t         n_vtx_seg,
 const int                 randomize,
 const int                 random_seed,
 int                     **pn_face,
 int                     **pn_vtx,
 int                    ***pface_vtx_idx,
 int                    ***pface_vtx,
 double                 ***pvtx_coord,
 CWP_g_num_t            ***pface_ln_to_gn,
 CWP_g_num_t            ***pvtx_ln_to_gn
 )
{
  CWP_UNUSED(randomize);
  CWP_UNUSED(random_seed);

  int          *pn_edge        = NULL;
  int         **pface_edge     = NULL;
  int         **pedge_vtx      = NULL;
  CWP_g_num_t **pedge_ln_to_gn = NULL;
  CWPT_generate_mesh_rectangle_ngon(comm,
                                    (CWPT_Mesh_nodal_elt_t) CWP_BLOCK_FACE_QUAD4,
                                    0.,
                                    0.,
                                    0.,
                                    1.,
                                    1.,
                                    n_vtx_seg,
                                    n_vtx_seg,
                                    n_part,
                                    part_method,
                                    0,
                                    pn_vtx,
                                    &pn_edge,
                                    pn_face,
                                    pvtx_coord,
                                    &pedge_vtx,
                                    pface_vtx_idx,
                                    &pface_edge,
                                    pface_vtx,
                                    pvtx_ln_to_gn,
                                    &pedge_ln_to_gn,
                                    pface_ln_to_gn);

  for (int i = 0; i < n_part; i++) {
    free(pface_edge    [i]);
    free(pedge_vtx     [i]);
    free(pedge_ln_to_gn[i]);
  }
  free(pn_edge);
  free(pface_edge    );
  free(pedge_vtx     );
  free(pedge_ln_to_gn);
}



/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  // Read args from command line
  int    n_vtx_seg1            = 4;
  int    n_vtx_seg2            = 4;
  int    randomize             = 1;
  double tolerance             = 1e-2;

  _read_args(argc,
             argv,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &tolerance,
             &randomize);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &comm_world_size);

  assert (comm_world_size > 1);


  // Initialize CWIPI
  int n_part = 1;
  int n_code = 1;
  int code_id[2];
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  int n_vtx_seg;
  CWP_Comm_t comm_type;
  if (rank < comm_world_size / 2) {
    code_id[0] = 1;
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
    n_vtx_seg = n_vtx_seg1;
    // comm_type = CWP_COMM_PAR_WITH_PART;
    comm_type = CWP_COMM_PAR_WITHOUT_PART;
  }
  else {
    code_id[0] = 2;
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
    n_vtx_seg = n_vtx_seg2;
    // comm_type = CWP_COMM_PAR_WITH_PART;
    comm_type = CWP_COMM_PAR_WITHOUT_PART;
  }


  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);

  CWP_Init(comm,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);


  if (rank == 0) {
    printf("CWIPI Init OK\n");
    fflush(stdout);
  }


  // Create coupling
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1_without_part";
  CWP_Cpl_create(code_name[0],                                          // Code name
                 cpl_name,                                              // Coupling id
                 coupled_code_name[0],                                  // Coupled application id
                 CWP_INTERFACE_SURFACE,
                 comm_type,                                             // Coupling type
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                 n_part,                                                // Partition number
                 CWP_DYNAMIC_MESH_STATIC,                               // Mesh displacement type
                 CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency

  MPI_Barrier(comm);

  if (rank == 0) {
    printf("Create coupling OK\n");
    fflush(stdout);
  }

  CWP_Visu_set(code_name[0],            // Code name
               cpl_name,                // Coupling id
               1,                       // Postprocessing frequency
               CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
               "text");                 // Postprocessing option

  MPI_Barrier(comm);

  if (rank == 0) {
    printf("Set visu OK\n");
    fflush(stdout);
  }


  // Define mesh

  int          *pn_face        = NULL;
  int          *pn_vtx         = NULL;
  int         **pface_vtx_idx  = NULL;
  int         **pface_vtx      = NULL;
  double      **pvtx_coord     = NULL;
  CWP_g_num_t **pface_ln_to_gn = NULL;
  CWP_g_num_t **pvtx_ln_to_gn  = NULL;
  _gen_mesh(intra_comm[0],
            n_part,
            CWPT_SPLIT_DUAL_WITH_HILBERT,
            n_vtx_seg,
            randomize,
            code_id[0],
            &pn_face,
            &pn_vtx,
            &pface_vtx_idx,
            &pface_vtx,
            &pvtx_coord,
            &pface_ln_to_gn,
            &pvtx_ln_to_gn);

  CWP_Mesh_interf_vtx_set(code_name[0],
                          cpl_name,
                          0,
                          pn_vtx[0],
                          pvtx_coord[0],
                          pvtx_ln_to_gn[0]);

  int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                           cpl_name,
                                           CWP_BLOCK_FACE_POLY);

  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   cpl_name,
                                   0,
                                   block_id,
                                   pn_face[0],
                                   pface_vtx_idx[0],
                                   pface_vtx[0],
                                   pface_ln_to_gn[0]);

  CWP_Mesh_interf_finalize(code_name[0], cpl_name);

  if (rank == 0) {
    printf("Set mesh OK\n");
    fflush(stdout);
  }


  // Create and set fields
  // field1: code1 -> code2
  // field2: code2 -> code1
  const char *field_name1 = "cooX_t0";
  const char *field_name2 = "code2_elt_gnum";
  double *send_val = NULL;
  double *recv_val = NULL;
  if (code_id[0] == 1) {
    send_val = (double *) malloc(sizeof(double *) * pn_vtx[0]);
    recv_val = (double *) malloc(sizeof(double *) * pn_face[0]);
    for (int i = 0; i < pn_vtx[0]; i++) {
      send_val[i] = pvtx_coord[0][3*i];
    }
  }
  else {
    send_val = (double *) malloc(sizeof(double) * pn_face[0]);
    recv_val = (double *) malloc(sizeof(double) * pn_vtx[0]);
    for (int i = 0; i < pn_face[0]; i++) {
      send_val[i] = (double) pface_ln_to_gn[0][i];
    }
  }

  CWP_Status_t visu_status = CWP_STATUS_ON;


  if (code_id[0] == 1) {
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    CWP_Field_data_set(code_name[0],
                       cpl_name,
                       field_name1,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_val);

    CWP_Involved_srcs_bcast_enable(code_name[0],
                                   cpl_name,
                                   field_name1);

    CWP_Field_data_set(code_name[0],
                       cpl_name,
                       field_name2,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_val);

    CWP_Computed_tgts_bcast_enable(code_name[0],
                                   cpl_name,
                                   field_name2);
  }
  else {
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);

    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    CWP_Field_data_set(code_name[0],
                       cpl_name,
                       field_name1,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_val);

    CWP_Computed_tgts_bcast_enable(code_name[0],
                                   cpl_name,
                                   field_name1);

    CWP_Field_data_set(code_name[0],
                       cpl_name,
                       field_name2,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_val);

    CWP_Involved_srcs_bcast_enable(code_name[0],
                                   cpl_name,
                                   field_name2);
  }

  MPI_Barrier(comm);

  if (rank == 0) {
    printf("Fields OK\n");
    fflush(stdout);
  }

  MPI_Barrier(comm);


  char char_tol[99];
  sprintf(char_tol, "%e", tolerance);
  CWP_Spatial_interp_property_set(code_name[0], cpl_name, "tolerance", CWP_DOUBLE, char_tol);
  CWP_Spatial_interp_weights_compute(code_name[0], cpl_name);



  if (code_id[0] == 1) {
    CWP_Field_issend(code_name[0], cpl_name, field_name1);
    CWP_Field_irecv (code_name[0], cpl_name, field_name2);
  }
  else {
    CWP_Field_irecv (code_name[0], cpl_name, field_name1);
    CWP_Field_issend(code_name[0], cpl_name, field_name2);
  }


  if (code_id[0] == 1) {
    CWP_Field_wait_issend(code_name[0], cpl_name, field_name1);
    CWP_Field_wait_irecv (code_name[0], cpl_name, field_name2);
  }
  else {
    CWP_Field_wait_irecv (code_name[0], cpl_name, field_name1);
    CWP_Field_wait_issend(code_name[0], cpl_name, field_name2);
  }


  double max_err = 0;
  CWP_g_num_t n_wrong = 0;
  if (code_id[0] != 1) {
    int n_located = CWP_N_computed_tgts_get(code_name[0], cpl_name, field_name1, 0);

    const int *located = CWP_Computed_tgts_get(code_name[0], cpl_name, field_name1, 0);

    for (int i = 0; i < n_located; i++) {
      int vtx_id = located[i] - 1;
      double err = ABS(pvtx_coord[0][3*vtx_id] - recv_val[vtx_id]);
      max_err = MAX(max_err, err);
      if (err > 1e-6) {
        n_wrong++;
        printf("!!! vtx %ld : err = %e\n", pvtx_ln_to_gn[0][vtx_id], err);
      }
    }

  }
  double global_max_err = 0.;
  MPI_Reduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  CWP_g_num_t global_n_wrong = 0;
  MPI_Reduce(&n_wrong, &global_n_wrong, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("Max error = %g (%ld wrong points)\n", global_max_err, global_n_wrong);
    fflush(stdout);
  }

  CWP_Time_step_end(code_name[0]);
  CWP_Mesh_interf_del(code_name[0], cpl_name);
  CWP_Cpl_del(code_name[0], cpl_name);


  for (int i_part = 0; i_part < n_part; i_part++) {
    free(pface_vtx_idx[i_part]);
    free(pface_vtx[i_part]);
    free(pvtx_coord[i_part]);
    free(pface_ln_to_gn[i_part]);
    free(pvtx_ln_to_gn[i_part]);
  }
  free(pn_face);
  free(pn_vtx);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pvtx_coord);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn);

  free(send_val);
  free(recv_val);


  free(coupled_code_name);
  free(code_name);
  free(intra_comm);

  //  Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
