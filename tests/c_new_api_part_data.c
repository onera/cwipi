/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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

#include "pdm_array.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"

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
  int                  *verbose
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  int verbose                     = 0;

  _read_args (argc,
              argv,
             &verbose);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  assert (comm_world_size % 2 == 0);

  // Initialize CWIPI
  int n_part = 1;
  int n_code = 1;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t *is_active_rank = malloc(sizeof(CWP_Status_t) * n_code);
  double *time_init = malloc(sizeof(double) * n_code);

  if (rank % 2 == 0) {
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
  }
  else {
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
  }

  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           time_init,
           intra_comm);

  if (verbose && rank == 0) {
    printf("CWIPI Init OK\n");
  }

  // Create coupling
  const char *coupling_name = "c_surf_cpl_P1P1";

  CWP_Spatial_interp_t loc_method = CWP_SPATIAL_INTERP_FROM_CLOSEST_POINT_LEAST_SQUARES;
  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_VOLUME,
                 CWP_COMM_PAR_WITH_PART,
                 loc_method,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  // Partitionned data exchange
  const char *part_data_name = "schtroumpf";

  // --> create
  CWP_PartData_exch_t side;
  if (rank % 2 == 0) {
    side = CWP_PARTDATA_SEND;
  }
  else {
    side = CWP_PARTDATA_RECV;
  }

  int n_elt;
  if (rank % 2 == 0) {
     n_elt = 8;
  } else {
     n_elt = 3;
  }
  int *n_elts = malloc(sizeof(int *) * n_part);
  for (int i_part; i_part < n_part; i_part++) {
    n_elts[i_part] = n_elt;
  }

  CWP_g_num_t **gnum_elt = malloc(sizeof(CWP_g_num_t *) * n_part);
  for (int i_part; i_part < n_part; i_part++) {
    gnum_elt[i_part] = malloc(sizeof(CWP_g_num_t) * n_elt);
  }

  for (int i_part; i_part < n_part; i_part++) {
    for (int i = 0; i < n_elt; i++) {
      if (rank % 2 == 0) {
        gnum_elt[i_part][i] = n_elt * rank + i + 1;
      } else {
        gnum_elt[i_part][i] = n_elt * (rank - 1) +i + 1;
      }
    }
  }

  CWP_Part_data_create(code_name[0],
                       coupling_name,
                       part_data_name,
                       side,
                       gnum_elt,
                       n_elts,
                       n_part);

  // --> exchange
  int **part1_to_part2_data = NULL;
  int send_request = -1;
  int **part2_data = NULL;
  int recv_request = -1;
  int n_comp = 3;

  if (rank % 2 == 0) {

    part1_to_part2_data = malloc(sizeof(int *) * 3 * n_elt);
    for (int i_part; i_part < n_part; i_part++) {
      part1_to_part2_data[i_part] = malloc(sizeof(int) * n_elt * n_comp);
    }

    for (int i_part; i_part < n_part; i_part++) {
      for (int i = 0; i < n_elt; i++) {
        for (int i_comp = 0; i_comp < n_comp; i_comp++) {
          part1_to_part2_data[i_part][3*i + i_comp] = 10*i + i_comp;
        }
      }
    }

    CWP_Part_data_issend(code_name[0],
                         coupling_name,
                         part_data_name,
                         sizeof(int),
                         n_comp,
                         (void **) part1_to_part2_data,
                         &send_request);
  }
  else {

    CWP_Part_data_irecv(code_name[0],
                        coupling_name,
                        part_data_name,
                        sizeof(int),
                        n_comp,
                        (void ***) &part2_data,
                        &recv_request);

  }

  MPI_Barrier(MPI_COMM_WORLD);

  // --> wait
  if (rank % 2 == 0) {

    CWP_Part_data_wait_issend(code_name[0],
                              coupling_name,
                              part_data_name,
                              &send_request);
  }
  else {

    CWP_Part_data_wait_irecv(code_name[0],
                             coupling_name,
                             part_data_name,
                             &recv_request);
  }

  // --> check
  for (int i_part; i_part < n_part; i_part++) {
      for (int i = 0; i < n_elt; i++) {
        for (int i_comp = 0; i_comp < n_comp; i_comp++) {
          if (rank % 2 == 0) {
            printf("%d - s[%d][%d][%d] : %d\n", rank, i_part, i, i_comp, part1_to_part2_data[i_part][3*i + i_comp]);
            fflush(stdout);
          }
          else {
            printf("%d - r[%d][%d][%d] : %d\n", rank, i_part, i, i_comp, part2_data[i_part][3*i + i_comp]);
            fflush(stdout);
          }
        }
      }
    }

  MPI_Barrier(MPI_COMM_WORLD);

  // Delete part_data object
  CWP_Part_data_del(code_name[0],
                    coupling_name,
                    part_data_name);

  // Delete coupling
  CWP_Cpl_del(code_name[0], coupling_name);

  // Free memory
  free(code_name);
  free(coupled_code_name);
  free(is_active_rank);
  free(time_init);
  free(intra_comm);

  if (part1_to_part2_data != NULL) {
    for (int i_part; i_part < n_part; i_part++) {
      free(part1_to_part2_data[i_part]);
    }
    free(part1_to_part2_data);
  }
  if (part2_data != NULL) {
    for (int i_part; i_part < n_part; i_part++) {
      if (part2_data[i_part] != NULL) free(part2_data[i_part]);
    }
    free(part2_data);
  }
  if (gnum_elt != NULL) {
    for (int i_part; i_part < n_part; i_part++) {
      if (gnum_elt[i_part] != NULL) free(gnum_elt[i_part]);
    }
    free(gnum_elt);
  }
  if (n_elts != NULL) free(n_elts);

  // Finalize cwipi
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}

