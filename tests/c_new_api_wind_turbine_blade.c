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

#include "pdm_multipart.h"
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
         "  -f1             first filename.\n\n"
         "  -f2             second filename.\n\n"
         "  -v              verbose.\n\n"
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
  int            argc,
  char         **argv,
  int           *verbose,
  char         **filename1,
  char         **filename2
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
    else if (strcmp(argv[i], "-f1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *filename1 = argv[i];
      }
    }
    else if (strcmp(argv[i], "-f2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *filename2 = argv[i];
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

/*----------------------------------------------------------------------
 *
 * Read VTK file
 *
 *---------------------------------------------------------------------*/

static void
_read_vtk
(
 MPI_Comm     comm,
 const char  *filename,
 int         *n_vtx,
 int         *n_face,
 int        **face_vtx,
 double     **vtx_coord,
 double     **vtx_field
)
{
  // MPI
  int i_rank;
  int n_rank;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  char line[999];

  if (i_rank == 0) {

    FILE *f = fopen(filename, "r");

    if (f == NULL) {
      PDM_error(__FILE__, __LINE__, 0, "Could not read file %s\n", filename);
    }

    // PASS 1
    while (1) {

      int stat = fscanf(f, "%s", line);

      if (stat == EOF) {
        // End of file
        break;
      }

      if (strstr(line, "POINTS") != NULL) {
        // Get dimension
        fscanf(f, "%d", &(*n_vtx));
        char *data_type = malloc(sizeof(char) * 10);
        fscanf(f, "%s", data_type);
        free(data_type);

        // Malloc
        (*vtx_coord) = malloc(sizeof(double) * 3 * ((*n_vtx)));

        // Get coordinates
        for (int i = 0; i < (*n_vtx); i++) {
          for (int j = 0; j < 3; j++) {
            fscanf(f, "%le", &(*vtx_coord)[3*i + j]);
          }
        }
      } // end if POINTS

      if (strstr(line, "CELLS") != NULL) {
        // Get dimension
        int n_cell = 0;
        fscanf(f, "%d", &n_cell);
        int size = 0;
        fscanf(f, "%d", &size);

        // Malloc
        (*face_vtx) = malloc(sizeof(int) * size);

        // Get connecitivity
        int n_vtx_cell = 0;
        for (int i = 0; i < n_cell; i++) {
          fscanf(f, "%d", &n_vtx_cell);
          if (n_vtx_cell == 3) {
            for (int j = 0; j < 3; j++) {
             fscanf(f, "%d", &(*face_vtx)[3*(*n_face)+j]);
            }
            (*n_face)++;
          } else {
            int other_cell_vtx = 0;
            for (int j = 0; j < n_vtx_cell; j++) {
              fscanf(f, "%d", &other_cell_vtx);
            }
          }
        }

        (*face_vtx) = realloc((*face_vtx), 3 * (*n_face) * sizeof(int));
      } // end if CELLS

      if (strstr(line, "LOOKUP_TABLE") != NULL) {
        char *table = malloc(sizeof(char) * 10);
        fscanf(f, "%s", table);
        free(table);

        // Malloc
        (*vtx_field) = malloc(sizeof(double) * (*n_vtx));

        // Get connectivity
        for (int i = 0; i < (*n_vtx); i++) {
          fscanf(f, "%lf", &(*vtx_field)[i]);
        }
      } // end if POINT_DATA

    } // while not EOF
  }

}

/*----------------------------------------------------------------------
 *
 * Main : preCICE wind turbine blade test case
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  int verbose = 0;
  char *filename1 = NULL;
  char *filename2 = NULL;

  _read_args (argc,
              argv,
              &verbose,
              &filename1,
              &filename2);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // i_rank == 0 read mesh
  int         n_vtx  = 0;
  int         n_face = 0;
  int        *face_vtx  = NULL;
  double     *vtx_coord = NULL;
  double     *vtx_field = NULL;
  _read_vtk(MPI_COMM_WORLD,
            filename1,
            &n_vtx,
            &n_face,
            &face_vtx,
            &vtx_coord,
            &vtx_field);

  if (verbose && i_rank == 0) {
    log_trace("n_vtx : %d\n", n_vtx);
    log_trace("n_face : %d\n", n_face);
    PDM_log_trace_array_double(vtx_coord, 3*n_vtx, "vtx_coord: ");
    PDM_log_trace_array_int(face_vtx, 3*n_face, "face_vtx: ");
    PDM_log_trace_array_double(vtx_field, n_vtx, "vtx_field: ");
  }

  // DMesh nodal
  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(PDM_MPI_COMM_WORLD,
                                                  3,
                                                  n_vtx,
                                                  0,
                                                  n_face,
                                                  0);

  PDM_DMesh_nodal_coord_set(dmn,
                            n_vtx,
                            vtx_coord,
                            PDM_OWNERSHIP_KEEP);

  int id_section = PDM_DMesh_nodal_section_add(dmn,
                                               PDM_GEOMETRY_KIND_SURFACIC,
                                               PDM_MESH_NODAL_TRIA3);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_MESH_NODAL_TRIA3,
                                  id_section,
                                  n_face,
                                  face_vtx, // should be PDM_g_num
                                  PDM_OWNERSHIP_KEEP);

  // Multipart
  int n_part = 1;
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  PDM_multipart_t *mpart = PDM_multipart_create(1,
                                                &n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                PDM_MPI_COMM_WORLD,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_register_dmesh_nodal(mpart, 0, dmn);
  PDM_multipart_run_ppart(mpart);

  // free
  PDM_DMesh_nodal_free(dmn);
  free(face_vtx );
  free(vtx_coord);

  // Partitionned mesh
  double *pvtx_coord = NULL;
  PDM_multipart_part_vtx_coord_get(mpart,
                                   0,
                                   0,
                                   &pvtx_coord,
                                   PDM_OWNERSHIP_KEEP);

  int  dn_section = 0;
  int *dn_elt = NULL;
  int  dn_cell = 0;
  int  dn_face = 0;
  int  dn_face_part_bound = 0;
  int  dn_vtx = 0;
  int  dn_proc = 0;
  int  dn_total_part = 0;
  int  ds_cell_face = 0;
  int  ds_face_vtx = 0;
  int  ds_face_bound = 0;
  int  dn_bound_groups = 0;
  int  ds_face_join = 0;
  int  dn_join_groups = 0;
  PDM_multipart_part_dim_get(mpart,
                             0,
                             0,
                             &dn_section,
                             &dn_elt,
                             &dn_cell,
                             &dn_face,
                             &dn_face_part_bound,
                             &dn_vtx,
                             &dn_proc,
                             &dn_total_part,
                             &ds_cell_face,
                             &ds_face_vtx,
                             &ds_face_bound,
                             &dn_bound_groups,
                             &ds_face_join,
                             &dn_join_groups);

  int         **delt_vtx_idx = NULL;
  int         **delt_vtx = NULL;
  PDM_g_num_t **delt_section_ln_to_gn = NULL;
  int          *dcell_tag = NULL;
  int          *dcell_face_idx = NULL;
  int          *dcell_face = NULL;
  PDM_g_num_t  *dcell_ln_to_gn = NULL;
  int          *dface_tag = NULL;
  int          *dface_cell = NULL;
  int          *dface_vtx_idx = NULL;
  int          *dface_vtx = NULL;
  PDM_g_num_t  *dface_ln_to_gn = NULL;
  int          *dface_part_bound_proc_idx = NULL;
  int          *dface_part_bound_part_idx = NULL;
  int          *dface_part_bound = NULL;
  int          *dvtx_tag = NULL;
  double       *dvtx = NULL;
  PDM_g_num_t  *dvtx_ln_to_gn = NULL;
  int          *dface_bound_idx = NULL;
  int          *dface_bound = NULL;
  PDM_g_num_t  *dface_bound_ln_to_gn = NULL;
  int          *dface_join_idx = NULL;
  int          *dface_join = NULL;
  PDM_g_num_t  *dface_join_ln_to_gn = NULL;
  PDM_multipart_part_val_get(mpart,
                             0,
                             0,
                             &delt_vtx_idx,
                             &delt_vtx,
                             &delt_section_ln_to_gn,
                             &dcell_tag,
                             &dcell_face_idx,
                             &dcell_face,
                             &dcell_ln_to_gn,
                             &dface_tag,
                             &dface_cell,
                             &dface_vtx_idx,
                             &dface_vtx,
                             &dface_ln_to_gn,
                             &dface_part_bound_proc_idx,
                             &dface_part_bound_part_idx,
                             &dface_part_bound,
                             &dvtx_tag,
                             &dvtx,
                             &dvtx_ln_to_gn,
                             &dface_bound_idx,
                             &dface_bound,
                             &dface_bound_ln_to_gn,
                             &dface_join_idx,
                             &dface_join,
                             &dface_join_ln_to_gn);

  // free
  PDM_multipart_free(mpart);

  // Partition field (part_to_part?)

  // free
  free(vtx_field);

  // Set up

  // Mesh

  // Field

  // Exchange

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
