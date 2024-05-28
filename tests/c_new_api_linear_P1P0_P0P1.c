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

#define ABS(a)    ((a) < 0   ? -(a) : (a))
#define MAX(a, b) ((a) > (b) ?  (a) : (b))

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
         "  -n_rank1        number of MPI ranks for code1.\n\n"
         "  -n_rank2        number of MPI ranks for code2.\n\n"
         "  -v              verbose.\n\n"
         "  -n1             square root of number of vertices for code1.\n\n"
         "  -n2             square root of number of vertices for code2.\n\n"
         "  -swap_codes     swap rank order of code1 and 2.\n\n"
         "  -algo           spatial interpolation algorithm.\n\n"
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
  int                   *swap_codes,
  CWP_g_num_t            all_gn_vtx[],
  int                    all_n_rank[],
  int                    all_n_part[],
  CWP_Spatial_interp_t  *spatial_interp_algo,
  double                *tolerance,
  int                   *n_neighbors,
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
    else if (strcmp(argv[i], "-swap_codes") == 0) {
      *swap_codes = 1;
    }
    else if (strcmp(argv[i], "-n1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_gn_vtx[0] = atol(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_gn_vtx[1] = atol(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_rank1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_rank[0] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_rank2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_rank[1] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_part[0] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_part[1] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-algo") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *spatial_interp_algo = (CWP_Spatial_interp_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-tol") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *tolerance = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_cls") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_neighbors = atoi(argv[i]);
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


static void
_gen_mesh
(
 const MPI_Comm             comm,
 const CWP_g_num_t          gn_vtx,
 const int                  n_part,
 int                      **pn_elt,
 int                      **pn_vtx,
 int                     ***pelt_vtx,
 double                  ***pvtx_coord,
 CWP_g_num_t             ***pelt_ln_to_gn,
 CWP_g_num_t             ***pvtx_ln_to_gn
 )
{
  int i_rank;
  MPI_Comm_rank(comm, &i_rank);

  int n_rank;
  MPI_Comm_size(comm, &n_rank);


  double step = 6.28 / (double) (gn_vtx - 1);


  *pn_elt        = malloc(sizeof(int          ) * n_part);
  *pn_vtx        = malloc(sizeof(int          ) * n_part);
  *pelt_vtx      = malloc(sizeof(int         *) * n_part);
  *pvtx_coord    = malloc(sizeof(double      *) * n_part);
  *pelt_ln_to_gn = malloc(sizeof(CWP_g_num_t *) * n_part);
  *pvtx_ln_to_gn = malloc(sizeof(CWP_g_num_t *) * n_part);

  int tn_part = n_part * n_rank;

  if (tn_part == 1) {
    (*pn_elt)[0] = (int) gn_vtx;
    (*pn_vtx)[0] = (int) gn_vtx;

    (*pelt_vtx     )[0] = malloc(sizeof(int        ) * (*pn_elt)[0] * 2);
    (*pvtx_coord   )[0] = malloc(sizeof(double     ) * (*pn_vtx)[0] * 3);
    (*pelt_ln_to_gn)[0] = malloc(sizeof(CWP_g_num_t) * (*pn_elt)[0]);
    (*pvtx_ln_to_gn)[0] = malloc(sizeof(CWP_g_num_t) * (*pn_vtx)[0]);

    for (int i = 0; i < (*pn_elt)[0]; i++) {
      (*pelt_ln_to_gn)[0][i] = i + 1;
      (*pelt_vtx     )[0][2*i  ] = i+1;
      (*pelt_vtx     )[0][2*i+1] = (i+1)%gn_vtx + 1;
    }

    for (int i = 0; i < (*pn_vtx)[0]; i++) {
      (*pvtx_ln_to_gn)[0][i] = i + 1;

      double r = 1 + 0.3*cos(5*i*step);
      double t = i*step + 0.1*sin(5*i*step);

      (*pvtx_coord)[0][3*i  ] = r*cos(t);
      (*pvtx_coord)[0][3*i+1] = r*sin(t);
      (*pvtx_coord)[0][3*i+2] = 0.1*cos(5*i*step);
    }

  }
  else {
    int div = (int) (gn_vtx / tn_part);
    int rem = (int) (gn_vtx % tn_part);
    CWP_g_num_t *offset_elt = malloc(sizeof(CWP_g_num_t) * (tn_part+1));
    offset_elt[0] = 0;
    for (int i_part = 0; i_part < tn_part; i_part++) {
      offset_elt[i_part+1] = offset_elt[i_part] + div;
      if (i_part < rem) {
        offset_elt[i_part+1]++;
      }
    }

    for (int i_part = 0; i_part < n_part; i_part++) {
      int j_part = n_part*i_rank + i_part;
      (*pn_elt)[i_part] = offset_elt[j_part+1] - offset_elt[j_part];
      (*pn_vtx)[i_part] = (*pn_elt)[i_part] + 1;

      (*pelt_vtx     )[i_part] = malloc(sizeof(int        ) * (*pn_elt)[i_part] * 2);
      (*pvtx_coord   )[i_part] = malloc(sizeof(double     ) * (*pn_vtx)[i_part] * 3);
      (*pelt_ln_to_gn)[i_part] = malloc(sizeof(CWP_g_num_t) * (*pn_elt)[i_part]);
      (*pvtx_ln_to_gn)[i_part] = malloc(sizeof(CWP_g_num_t) * (*pn_vtx)[i_part]);


      for (int i = 0; i < (*pn_elt)[i_part]; i++) {
        (*pelt_ln_to_gn)[i_part][i] = offset_elt[j_part] + i + 1;
        (*pelt_vtx     )[i_part][2*i  ] = i+1;
        (*pelt_vtx     )[i_part][2*i+1] = i+2;
      }

      for (int i = 0; i <= (*pn_elt)[i_part]; i++) {
        CWP_g_num_t g = (offset_elt[j_part] + i) % gn_vtx;

        (*pvtx_ln_to_gn)[i_part][i] = g + 1;

        double r = 1 + 0.3*cos(5*g*step);
        double t = g*step + 0.1*sin(5*g*step);

        (*pvtx_coord)[i_part][3*i  ] = r*cos(t);
        (*pvtx_coord)[i_part][3*i+1] = r*sin(t);
        (*pvtx_coord)[i_part][3*i+2] = 0.1*cos(5*g*step);
      }
    }

    free(offset_elt);
  }




}


/*----------------------------------------------------------------------
 *
 * Main : Linear coupling interface
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
  int                  swap_codes     = 0;
  CWP_g_num_t          all_gn_vtx[2]  = {100, 50};
  int                  all_n_rank[2]  = {-1, -1};
  int                  all_n_part[2]  = {1, 1};
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  double               tolerance      = 1e-2;
  int                  n_neighbors    = 5;
  int                  visu           = 0;

  _read_args(argc,
             argv,
             &verbose,
             &swap_codes,
             all_gn_vtx,
             all_n_rank,
             all_n_part,
             &spatial_interp,
             &tolerance,
             &n_neighbors,
             &visu);


  // Initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  for (int i = 0; i < 2; i++) {
    if (all_n_rank[i] <= 0) {
      all_n_rank[i] = n_rank;
    }
  }


  const char *all_code_names[2] = {"code1", "code2"};
  int has_code[2] = {0, 0};


  has_code[0] = i_rank <  all_n_rank[0];
  has_code[1] = i_rank >= n_rank - all_n_rank[1];

  int n_code = has_code[0] + has_code[1];

  int           *code_id           = malloc(sizeof(int         ) * n_code);
  const char   **code_name         = malloc(sizeof(char       *) * n_code);
  const char   **coupled_code_name = malloc(sizeof(char       *) * n_code);
  CWP_Status_t   is_active_rank    = CWP_STATUS_ON;
  MPI_Comm      *intra_comm        = malloc(sizeof(MPI_Comm    ) * n_code);
  int            n_part[2];
  CWP_g_num_t    gn_vtx[2];

  n_code = 0;
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      gn_vtx           [n_code] = all_gn_vtx[icode];
      n_part           [n_code] = all_n_part[icode];

      //if (verbose) {
      //  log_trace("Running %s, coupled with %s\n",
      //            code_name[n_code], coupled_code_name[n_code]);
      //}
      n_code++;
    }
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
  const char *cpl_name = "c_new_api_linear";

  for (int icode = 0; icode < n_code; icode++) {
    CWP_Cpl_create(code_name[icode],
                   cpl_name,
                   coupled_code_name[icode],
                   CWP_INTERFACE_LINEAR,
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
  int          **pn_vtx        = malloc(sizeof(int          *) * n_code);
  int          **pn_elt        = malloc(sizeof(int          *) * n_code);
  double      ***pvtx_coord    = malloc(sizeof(double      **) * n_code);
  int         ***pelt_vtx      = malloc(sizeof(int         **) * n_code);
  CWP_g_num_t ***pelt_ln_to_gn = malloc(sizeof(CWP_g_num_t **) * n_code);
  CWP_g_num_t ***pvtx_ln_to_gn = malloc(sizeof(CWP_g_num_t **) * n_code);

  for (int icode = 0; icode < n_code; icode++) {
    _gen_mesh(intra_comm[icode],
              gn_vtx[icode],
              n_part[icode],
              &pn_elt       [icode],
              &pn_vtx       [icode],
              &pelt_vtx     [icode],
              &pvtx_coord   [icode],
              &pelt_ln_to_gn[icode],
              &pvtx_ln_to_gn[icode]);

    int block_id = CWP_Mesh_interf_block_add(code_name[icode],
                                             cpl_name,
                                             CWP_BLOCK_EDGE2);

    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      CWP_Mesh_interf_vtx_set(code_name[icode],
                              cpl_name,
                              ipart,
                              pn_vtx       [icode][ipart],
                              pvtx_coord   [icode][ipart],
                              pvtx_ln_to_gn[icode][ipart]);

      CWP_Mesh_interf_block_std_set(code_name[icode],
                                    cpl_name,
                                    ipart,
                                    block_id,
                                    pn_elt       [icode][ipart],
                                    pelt_vtx     [icode][ipart],
                                    pelt_ln_to_gn[icode][ipart]);
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

    field1_val[icode] = malloc(sizeof(double *) * n_part[icode]);
    field2_val[icode] = malloc(sizeof(double *) * n_part[icode]);

    if (code_id[icode] == 1) {
      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       3,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Time_step_beg(code_name[icode],
                        0.0);

      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
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

      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        field1_val[icode][ipart] = malloc(sizeof(double) * pn_elt[icode][ipart] * 3);
        field2_val[icode][ipart] = malloc(sizeof(double) * pn_elt[icode][ipart]);
        for (int i = 0; i < pn_elt[icode][ipart]; i++) {
          int vtx_id0 = pelt_vtx[icode][ipart][2*i  ] - 1;
          int vtx_id1 = pelt_vtx[icode][ipart][2*i+1] - 1;
          field2_val[icode][ipart][i] = 0.5*(pvtx_coord[icode][ipart][3*vtx_id0] + pvtx_coord[icode][ipart][3*vtx_id1]);
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
    char char_param[99];
    sprintf(char_param, "%e", tolerance);
    CWP_Spatial_interp_property_set(code_name[icode],
                                    cpl_name,
                                    "tolerance",
                                    CWP_DOUBLE,
                                    char_param);

    sprintf(char_param, "%d", n_neighbors);
    CWP_Spatial_interp_property_set(code_name[icode],
                                    cpl_name,
                                    "n_neighbors",
                                    CWP_INT,
                                    char_param);

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


  /* Check interpolated fields */
  // for (int icode = 0; icode < n_code; icode++) {
  //   if (code_id[icode] == 1) {

  //   }
  //   else {

  //   }
  // }



  /* Finalize */
  for (int icode = 0; icode < n_code; icode++) {
    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      free(pvtx_coord   [icode][ipart]);
      free(pvtx_ln_to_gn[icode][ipart]);
      free(pelt_vtx     [icode][ipart]);
      free(pelt_ln_to_gn[icode][ipart]);
      free(field1_val   [icode][ipart]);
      free(field2_val   [icode][ipart]);
    }
    free(pn_vtx       [icode]);
    free(pn_elt       [icode]);
    free(pvtx_coord   [icode]);
    free(pvtx_ln_to_gn[icode]);
    free(pelt_vtx     [icode]);
    free(pelt_ln_to_gn[icode]);
    free(field1_val   [icode]);
    free(field2_val   [icode]);
  }
  free(pn_vtx       );
  free(pn_elt       );
  free(pvtx_coord   );
  free(pvtx_ln_to_gn);
  free(pelt_vtx     );
  free(pelt_ln_to_gn);
  free(field1_val   );
  free(field2_val   );


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
