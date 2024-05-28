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
  int                   *n_part1,
  int                   *n_part2,
  double                *tolerance,
  int                   *randomize,
  int                   *verbose
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
     else if (strcmp(argv[i], "-n_part1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part1 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part2 = atoi(argv[i]);
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
    else if (strcmp(argv[i], "-verbose") == 0) {
      *verbose = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_gen_mesh
(
 const MPI_Comm             comm,
 const int                  n_part,
 const CWPT_split_dual_t    part_method,
 const CWP_g_num_t          n_vtx_seg,
 const int                  randomize,
 const int                  is_partitioned,
 int                      **pn_face,
 int                      **pn_vtx,
 int                     ***pface_vtx_idx,
 int                     ***pface_vtx,
 double                  ***pvtx_coord,
 CWP_g_num_t             ***pface_ln_to_gn,
 CWP_g_num_t             ***pvtx_ln_to_gn
 )
{
  int i_rank;
  MPI_Comm_rank(comm, &i_rank);

  MPI_Comm mesh_comm;
  if (is_partitioned) {
    MPI_Comm_dup(comm, &mesh_comm);
  }
  else {
    MPI_Comm_split(comm, i_rank, i_rank, &mesh_comm);
  }


  int          *pn_edge        = NULL;
  int         **pedge_vtx      = NULL;
  int         **pface_edge     = NULL;
  CWP_g_num_t **pedge_ln_to_gn = NULL;

  double random_factor = 0.;
  if (randomize) {
    random_factor = 1.;
  }

  CWPT_generate_mesh_rectangle_ngon(mesh_comm,
                                    CWPT_MESH_NODAL_POLY_2D,
                                    0.,
                                    0.,
                                    0.,
                                    1.,
                                    1.,
                                    n_vtx_seg,
                                    n_vtx_seg,
                                    n_part,
                                    part_method,
                                    random_factor,
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

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(pedge_vtx     [i_part]);
    free(pface_edge    [i_part]);
    free(pedge_ln_to_gn[i_part]);
  }
  free(pn_edge       );
  free(pedge_vtx     );
  free(pface_edge    );
  free(pedge_ln_to_gn);

  MPI_Comm_free(&mesh_comm);
}

/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P0P0
 *
 *---------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  // Read args from command line
  int    n_vtx_seg1            = 4;
  int    n_vtx_seg2            = 4;
  int    n_part1               = 1;
  int    n_part2               = 1;
  int    randomize             = 1;
  double tolerance             = 1e-2;
  int    verbose               = 0;
  CWPT_split_dual_t part_method = CWPT_SPLIT_DUAL_WITH_HILBERT;


  _read_args(argc,
             argv,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &n_part1,
             &n_part2,
             &tolerance,
             &randomize,
             &verbose);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  assert(n_rank > 1);

  FILE *file_log = NULL;
  if (verbose) {
    char file_log_name[999];
    sprintf(file_log_name, "c_new_api_surf_P0P0_without_part_%d.log", i_rank);
    file_log = fopen(file_log_name, "w");
    if (file_log == NULL) {
      printf("Warning : failed to open %s\n", file_log_name);
      verbose = 0;
    }
  }

  // Initialize CWIPI

  const char **code_name         = malloc(sizeof(char *) * 2);
  const char **coupled_code_name = malloc(sizeof(char *) * 2);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;


  int has_code[2] = {0, 0};
  if (i_rank < (2*n_rank) / 3) {
    has_code[0] = 1;
  }
  if (i_rank >= n_rank / 3) {
    has_code[1] = 1;
  }

  const char *all_code_names[2] = {"code1", "code2"};
  int all_n_vtx_seg[2] = {n_vtx_seg1, n_vtx_seg2};
  int all_n_part   [2] = {n_part1,    n_part2};
  CWP_Comm_t all_comm_type[2] = {CWP_COMM_PAR_WITHOUT_PART, CWP_COMM_PAR_WITH_PART};

  int n_code = 0;
  int n_vtx_seg[2];
  int n_part   [2];
  int code_id  [2];
  CWP_Comm_t comm_type[2];
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      n_vtx_seg        [n_code] = all_n_vtx_seg [icode];
      n_part           [n_code] = all_n_part    [icode];
      comm_type        [n_code] = all_comm_type[icode];
      if (all_comm_type[icode] == CWP_COMM_PAR_WITHOUT_PART) {
        n_part[n_code] = 1;
      }
      // fprintf(file_log, "%s\n", code_name[n_code]);
      n_code++;
    }
  }

  // fprintf(file_log, "n_code = %d : ", n_code);
  // for (int i_code = 0; i_code < n_code; i_code++) {
  //   fprintf(file_log, "%s ", code_name[i_code]);
  // }
  // fprintf(file_log, "\n");

  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);

  CWP_Init(comm,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  if (i_rank == 0) {
    printf("CWIPI Init OK\n");
  }


  // Create coupling
  const char *cpl_name = "c_new_api_surf_P0P0_without_part";
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;

  for (int i_code = 0; i_code < n_code; i_code++) {
    // fprintf(file_log, ">> CWP_Cpl_create %s\n", code_name[i_code]);
    CWP_Cpl_create(code_name[i_code],                                     // Code name
                   cpl_name,                                              // Coupling id
                   coupled_code_name[i_code],                             // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   comm_type[i_code],                                     // Coupling type
                   spatial_interp,
                   n_part[i_code],                                        // Number of partitions
                   CWP_DYNAMIC_MESH_STATIC,                               // Mesh displacement type
                   CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency
  }



  for (int i_code = 0; i_code < n_code; i_code++) {
    CWP_Visu_set(code_name[i_code],       // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }
  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Create coupling OK\n");
  }

  // Mesh definition
  int          **pn_face        = malloc(sizeof(int          *) * n_code);
  int          **pn_vtx         = malloc(sizeof(int          *) * n_code);
  int         ***pface_vtx_idx  = malloc(sizeof(int         **) * n_code);
  int         ***pface_vtx      = malloc(sizeof(int         **) * n_code);
  double      ***pvtx_coord     = malloc(sizeof(double      **) * n_code);
  CWP_g_num_t ***pface_ln_to_gn = malloc(sizeof(CWP_g_num_t **) * n_code);
  CWP_g_num_t ***pvtx_ln_to_gn  = malloc(sizeof(CWP_g_num_t **) * n_code);


  for (int i_code = 0; i_code < n_code; i_code++) {
    _gen_mesh(intra_comm[i_code],
              n_part[i_code],
              part_method,
              n_vtx_seg[i_code],
              randomize,
              (comm_type[i_code] == CWP_COMM_PAR_WITH_PART),
              &pn_face       [i_code],
              &pn_vtx        [i_code],
              &pface_vtx_idx [i_code],
              &pface_vtx     [i_code],
              &pvtx_coord    [i_code],
              &pface_ln_to_gn[i_code],
              &pvtx_ln_to_gn [i_code]);

    int block_id = CWP_Mesh_interf_block_add(code_name[i_code],
                                             cpl_name,
                                             CWP_BLOCK_FACE_POLY);
    for (int i = 0; i < n_part[i_code]; i++) {
      CWP_Mesh_interf_vtx_set(code_name[i_code],
                              cpl_name,
                              i,
                              pn_vtx       [i_code][i],
                              pvtx_coord   [i_code][i],
                              pvtx_ln_to_gn[i_code][i]);


      CWP_Mesh_interf_f_poly_block_set(code_name[i_code],
                                       cpl_name,
                                       i,
                                       block_id,
                                       pn_face       [i_code][i],
                                       pface_vtx_idx [i_code][i],
                                       pface_vtx     [i_code][i],
                                       pface_ln_to_gn[i_code][i]);
    }

    CWP_Mesh_interf_finalize(code_name[i_code], cpl_name);
  }
  MPI_Barrier(comm);

  if (i_rank == 0) {
    printf("Set mesh OK\n");
  }


  // Create and set fields
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name1 = "field1";
  const char *field_name2 = "field2";

  double ***send_val = malloc(sizeof(double **) * n_code);
  double ***recv_val = malloc(sizeof(double **) * n_code);
  for (int i_code = 0; i_code < n_code; i_code++) {
    send_val[i_code] = malloc(sizeof(double *) * n_part[i_code]);
    recv_val[i_code] = malloc(sizeof(double *) * n_part[i_code]);
    for (int i = 0; i < n_part[i_code]; i++) {
      send_val[i_code][i] = malloc(sizeof(double) * pn_face[i_code][i]);
      recv_val[i_code][i] = malloc(sizeof(double) * pn_face[i_code][i]);
    }

    if (code_id[i_code] == 1) {
      for (int ipart = 0; ipart < n_part[i_code]; ipart++) {
        for (int i = 0 ; i < pn_face[i_code][ipart]; i++) {
          send_val[i_code][ipart][i] = (double) rand() / (double) RAND_MAX;
        }
      }

      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Time_step_beg(code_name[i_code],
                        0.0);

      for (int i = 0; i < n_part[i_code]; i++) {
        CWP_Field_data_set(code_name[i_code],
                           cpl_name,
                           field_name1,
                           i,
                           CWP_FIELD_MAP_SOURCE,
                           send_val[i_code][i]);
      }

      CWP_Involved_srcs_bcast_enable(code_name[i_code],
                                     cpl_name,
                                     field_name1);

      for (int i = 0; i < n_part[i_code]; i++) {
        CWP_Field_data_set(code_name[i_code],
                           cpl_name,
                           field_name2,
                           i,
                           CWP_FIELD_MAP_TARGET,
                           recv_val[i_code][i]);
      }

      CWP_Computed_tgts_bcast_enable(code_name[i_code],
                                     cpl_name,
                                     field_name2);
    }

    if (code_id[i_code] == 2) {
      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Time_step_beg(code_name[i_code],
                        0.0);

      for (int i = 0; i < n_part[i_code]; i++) {
        CWP_Field_data_set(code_name[i_code],
                           cpl_name,
                           field_name1,
                           i,
                           CWP_FIELD_MAP_TARGET,
                           recv_val[i_code][i]);
      }

      for (int i = 0; i < n_part[i_code]; i++) {
        CWP_Field_data_set(code_name[i_code],
                           cpl_name,
                           field_name2,
                           i,
                           CWP_FIELD_MAP_SOURCE,
                           send_val[i_code][i]);
      }
    }
  }

  MPI_Barrier(comm);

  if (i_rank == 0) {
    printf("Create fields OK\n");
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    CWP_Spatial_interp_weights_compute(code_name[i_code], cpl_name);
  }
  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Interpolation weights computation OK\n");
  }

  MPI_Barrier(comm);

  for (int i_code = 0; i_code < n_code; i_code++) {
    if (code_id[i_code] == 1) {
      CWP_Field_issend(code_name[i_code], cpl_name, field_name1);
      CWP_Field_irecv (code_name[i_code], cpl_name, field_name2);
    }
    else {
      CWP_Field_irecv (code_name[i_code], cpl_name, field_name1);
      CWP_Field_issend(code_name[i_code], cpl_name, field_name2);
    }

    if (verbose) {
      fprintf(file_log, "\n\n--- %s ---\n", code_name[i_code]);
    }


    if (code_id[i_code] == 1) {
      CWP_Field_wait_issend(code_name[i_code], cpl_name, field_name1);
      CWP_Field_wait_irecv (code_name[i_code], cpl_name, field_name2);
      for (int ipart = 0; ipart < n_part[i_code]; ipart++) {

        int n_involved_src = CWP_N_involved_srcs_get(code_name[i_code],
                                                     cpl_name,
                                                     field_name1,
                                                     ipart);
        const int *involved_src = CWP_Involved_srcs_get(code_name[i_code],
                                                        cpl_name,
                                                        field_name1,
                                                        ipart);
        if (verbose) {
          fprintf(file_log, "Field 1 :\n");
          fprintf(file_log, "  part %d, involved_src : ", ipart);
          for (int i = 0; i < n_involved_src; i++) {
            fprintf(file_log, "%d ", involved_src[i]);
          }
          fprintf(file_log, "\n");
        }

        int n_computed_tgt = CWP_N_computed_tgts_get(code_name[i_code],
                                                     cpl_name,
                                                     field_name2,
                                                     ipart);
        const int *computed_tgt = CWP_Computed_tgts_get(code_name[i_code],
                                                        cpl_name,
                                                        field_name2,
                                                        ipart);
        if (verbose) {
          fprintf(file_log, "Field 2 :\n");
          fprintf(file_log, "  part %d, computed_tgt : ", ipart);
          for (int i = 0; i < n_computed_tgt; i++) {
            fprintf(file_log, "%d ", computed_tgt[i]);
          }
          fprintf(file_log, "\n");
        }
      }
    }
    else {
      CWP_Field_wait_irecv (code_name[i_code], cpl_name, field_name1);
      CWP_Field_wait_issend(code_name[i_code], cpl_name, field_name2);
      for (int ipart = 0; ipart < n_part[i_code]; ipart++) {
        int n_computed_tgt = CWP_N_computed_tgts_get(code_name[i_code],
                                                     cpl_name,
                                                     field_name1,
                                                     ipart);
        const int *computed_tgt = CWP_Computed_tgts_get(code_name[i_code],
                                                        cpl_name,
                                                        field_name1,
                                                        ipart);
        if (verbose) {
          fprintf(file_log, "Field 1 :\n");
          fprintf(file_log, "  part %d, computed_tgt : ", ipart);
          for (int i = 0; i < n_computed_tgt; i++) {
            fprintf(file_log, "%d ", computed_tgt[i]);
          }
          fprintf(file_log, "\n");
        }
        int n_involved_src = CWP_N_involved_srcs_get(code_name[i_code],
                                                     cpl_name,
                                                     field_name2,
                                                     ipart);
        const int *involved_src = CWP_Involved_srcs_get(code_name[i_code],
                                                        cpl_name,
                                                        field_name2,
                                                        ipart);
        if (verbose) {
          fprintf(file_log, "Field 2 :\n");
          fprintf(file_log, "  part %d, involved_src : ", ipart);
          for (int i = 0; i < n_involved_src; i++) {
            fprintf(file_log, "%d ", involved_src[i]);
          }
          fprintf(file_log, "\n");
        }
      }
    }
  }

  if (i_rank == 0) {
    printf("Exchange fields OK\n");
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    CWP_Time_step_end(code_name[i_code]);

    CWP_Mesh_interf_del(code_name[i_code], cpl_name);

    CWP_Cpl_del(code_name[i_code], cpl_name);
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    for (int i_part = 0 ; i_part < n_part[i_code]; i_part++) {
      free(pface_vtx_idx [i_code][i_part]);
      free(pface_vtx     [i_code][i_part]);
      free(pvtx_coord    [i_code][i_part]);
      free(pface_ln_to_gn[i_code][i_part]);
      free(pvtx_ln_to_gn [i_code][i_part]);
      free(send_val      [i_code][i_part]);
      free(recv_val      [i_code][i_part]);
    }
    free(pn_face       [i_code]);
    free(pn_vtx        [i_code]);
    free(pface_vtx_idx [i_code]);
    free(pface_vtx     [i_code]);
    free(pvtx_coord    [i_code]);
    free(pface_ln_to_gn[i_code]);
    free(pvtx_ln_to_gn [i_code]);
    free(send_val      [i_code]);
    free(recv_val      [i_code]);
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

  if (verbose) {
    fclose(file_log);
  }

  //  Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}

