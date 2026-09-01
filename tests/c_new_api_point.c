/*
  This file is part of the CWIPI library.

  Copyright (C) 2023  ONERA

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
#include <math.h>

#include "cwp.h"
#include "cwipi_config.h"
#include "cwp_priv.h"

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
         "  -v              verbose.\n\n"
         "  -n_part1        number of parts for code1.\n\n"
         "  -n_part2        number of parts for code2.\n\n"
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
  int                   *verbose,
  int                    all_n_part[],
  int                   *visu
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
    else if (strcmp(argv[i], "-all_n_part1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_part[0] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-all_n_part2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_part[1] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/*----------------------------------------------------------------------
 *
 * Main : Point coupling interface
 *
 *---------------------------------------------------------------------*/

int
main
(
 int   argc,
 char *argv[]
 )
{
  int                  verbose        = 0;
  int                  all_n_part[2]  = {1, 1};
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_IDENTITY;
  int                  visu           = 0;

  _read_args(argc,
             argv,
             &verbose,
             all_n_part,
             &visu);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  const char *all_code_names[2] = {"code1", "code2"};
  int n_code = 2;

  int           *code_id           = malloc(sizeof(int         ) * n_code);
  const char   **code_name         = malloc(sizeof(char       *) * n_code);
  const char   **coupled_code_name = malloc(sizeof(char       *) * n_code);
  CWP_Status_t   is_active_rank    = CWP_STATUS_ON;
  MPI_Comm      *intra_comm        = malloc(sizeof(MPI_Comm    ) * n_code);

  for (int icode = 0; icode < n_code; icode++) {
    code_id          [icode] = icode+1;
    code_name        [icode] = all_code_names[icode];
    coupled_code_name[icode] = all_code_names[(icode+1)%n_code];
  }

  // Set up
  CWP_Init(comm,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("CWIPI Init OK\n");
    fflush(stdout);
  }

  /* Create coupling */
  const char *cpl_name = "c_new_api_point";

  for (int icode = 0; icode < n_code; icode++) {
    CWP_Cpl_create(code_name[icode],
                   cpl_name,
                   coupled_code_name[icode],
                   CWP_INTERFACE_POINT,
                   CWP_COMM_PAR_WITH_PART,
                   spatial_interp,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }


  for (int icode = 0; icode < n_code; icode++) {
    CWP_Visu_set(code_name[icode],        // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Create coupling OK\n");
    fflush(stdout);
  }

  /* Define interface mesh */
  int     **pn_vtx       = malloc(sizeof(int     *) * n_code);
  double ***pvtx_coord   = malloc(sizeof(double **) * n_code);
  int    ***point_connec = malloc(sizeof(int    **) * n_code);

  for (int icode = 0; icode < n_code; icode++) {
    pn_vtx    [icode]   = malloc(sizeof(int     ) * all_n_part[icode]);
    pvtx_coord[icode]   = malloc(sizeof(double *) * all_n_part[icode]);
    point_connec[icode] = malloc(sizeof(int    *) * all_n_part[icode]);

    int block_id = CWP_Mesh_interf_block_add(code_name[icode],
                                             cpl_name,
                                             CWP_BLOCK_NODE);

    for (int ipart = 0; ipart < all_n_part[icode]; ipart++) {

      if (ipart == 0 && i_rank == icode) {
        pn_vtx    [icode][ipart]    = 1;
        pvtx_coord[icode][ipart]    = malloc(sizeof(double) * 3);
        pvtx_coord[icode][ipart][0] = 0.0;
        pvtx_coord[icode][ipart][1] = 0.0;
        pvtx_coord[icode][ipart][2] = 0.0;
      }
      else {
        pn_vtx[icode][ipart]     = 0;
        pvtx_coord[icode][ipart] = malloc(sizeof(double) * 0);
      }

      CWP_Mesh_interf_vtx_set(code_name[icode],
                              cpl_name,
                              ipart,
                              pn_vtx    [icode][ipart],
                              pvtx_coord[icode][ipart],
                              NULL);

    point_connec[icode][ipart] = malloc(sizeof(int) * pn_vtx[icode][ipart]);
    for (int i = 0; i < pn_vtx[icode][ipart]; i++) {
      point_connec[icode][ipart][i] = i + 1;
    }

    CWP_Mesh_interf_block_std_set(code_name[icode],
                                cpl_name,
                                ipart,
                                block_id,
                                pn_vtx      [icode][ipart],
                                point_connec[icode][ipart],
                                NULL);
    }

    CWP_Mesh_interf_finalize(code_name[icode], cpl_name);
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Set mesh OK\n");
    fflush(stdout);
  }

  /* Define fields */
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name1 = "all_coords";
  const char *field_name2 = "coordX";

  double ***field1_val = malloc(sizeof(double **) * n_code);
  double ***field2_val = malloc(sizeof(double **) * n_code);

  for (int icode = 0; icode < n_code; icode++) {

    field1_val[icode] = malloc(sizeof(double *) * all_n_part[icode]);
    field2_val[icode] = malloc(sizeof(double *) * all_n_part[icode]);

    if (code_id[icode] == 1) {
      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       3,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Time_step_beg(code_name[icode],
                        0.0);

      for (int ipart = 0; ipart < all_n_part[icode]; ipart++) {
        field1_val[icode][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart] * 3);
        field2_val[icode][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart]);
        for (int i = 0; i < 3*pn_vtx[icode][ipart]; i++) {
          field1_val[icode][ipart][i] = pvtx_coord[icode][ipart][i];
        }

        CWP_Field_data_set(code_name[icode],
                           cpl_name,
                           field_name1,
                           ipart,
                           CWP_FIELD_MAP_SOURCE,
                           field1_val[icode][ipart]);

        CWP_Field_data_set(code_name[icode],
                           cpl_name,
                           field_name2,
                           ipart,
                           CWP_FIELD_MAP_TARGET,
                           field2_val[icode][ipart]);
      }
    }
    else {
      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       3,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Time_step_beg(code_name[icode],
                        0.0);

      for (int ipart = 0; ipart < all_n_part[icode]; ipart++) {
        field1_val[icode][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart] * 3);
        field2_val[icode][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart]);
        for (int i = 0; i < pn_vtx[icode][ipart]; i++) {
          field2_val[icode][ipart][i] = pvtx_coord[icode][ipart][3*i];
        }

        CWP_Field_data_set(code_name[icode],
                           cpl_name,
                           field_name1,
                           ipart,
                           CWP_FIELD_MAP_TARGET,
                           field1_val[icode][ipart]);

        CWP_Field_data_set(code_name[icode],
                           cpl_name,
                           field_name2,
                           ipart,
                           CWP_FIELD_MAP_SOURCE,
                           field2_val[icode][ipart]);
      }
    }
  }

  /* Exchange fields */
  for (int icode = 0; icode < n_code; icode++) {
    CWP_Spatial_interp_weights_compute(code_name[icode], cpl_name);
  }

  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 1) {
      CWP_Field_issend(code_name[icode], cpl_name, field_name1);
      CWP_Field_irecv (code_name[icode], cpl_name, field_name2);
    }
    else {
      CWP_Field_irecv (code_name[icode], cpl_name, field_name1);
      CWP_Field_issend(code_name[icode], cpl_name, field_name2);
    }
  }

  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 1) {
      CWP_Field_wait_issend(code_name[icode], cpl_name, field_name1);
      CWP_Field_wait_irecv (code_name[icode], cpl_name, field_name2);
    }
    else {
      CWP_Field_wait_irecv (code_name[icode], cpl_name, field_name1);
      CWP_Field_wait_issend(code_name[icode], cpl_name, field_name2);
    }
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Exchange fields OK\n");
    fflush(stdout);
  }

  /* Finalize */
  for (int icode = 0; icode < n_code; icode++) {
    for (int ipart = 0; ipart < all_n_part[icode]; ipart++) {
      free(pvtx_coord  [icode][ipart]);
      free(point_connec[icode][ipart]);
      free(field1_val  [icode][ipart]);
      free(field2_val  [icode][ipart]);
    }
    free(pn_vtx      [icode]);
    free(point_connec[icode]);
    free(pvtx_coord  [icode]);
    free(field1_val  [icode]);
    free(field2_val  [icode]);
  }
  free(pn_vtx      );
  free(pvtx_coord  );
  free(point_connec);
  free(field1_val  );
  free(field2_val  );

  for (int icode = 0; icode < n_code; icode++) {
    CWP_Time_step_end(code_name[icode]);
    CWP_Mesh_interf_del(code_name[icode], cpl_name);
    CWP_Cpl_del        (code_name[icode], cpl_name);
  }
  free(code_id);
  free(coupled_code_name);
  free(code_name);
  free(intra_comm);

  CWP_Finalize();

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("End\n");
    fflush(stdout);
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}
