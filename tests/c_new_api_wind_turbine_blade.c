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
         "  -f              filename.\n\n"
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
  char         **filename
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
    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *filename = argv[i];
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

  // data
  int         _n_vtx     = *n_vtx;
  int         _n_face    = * n_face;
  int        *_face_vtx  = *face_vtx;
  double     *_vtx_coord = *vtx_coord;
  double     *_vtx_field = *vtx_field;

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
        fscanf(f, "%d", &_n_vtx);

        // Malloc
        _vtx_coord = malloc(sizeof(double) * 3 * (_n_vtx));

        // Get coordinates
        for (int i = 0; i < _n_vtx; i++) {
          for (int j = 0; j < 3; j++) {
            fscanf(f, "%lf", &_vtx_coord[3*i + j]);
          }
        }

        log_trace("_n_vtx : %d\n", _n_vtx);
        PDM_log_trace_array_double(_vtx_coord, 3*_n_vtx, "_vtx_coord: ");
      } // end if POINTS

      if (strstr(line, "CELLS") != NULL) {
        // Get dimension
        int n_cell = 0;
        fscanf(f, "%d", &n_cell);
        int size = 0;
        fscanf(f, "%d", &size);

        // Malloc
        _face_vtx = malloc(sizeof(int) * size);

        // Get connecitivity
        int n_vtx_cell = 0;
        for (int i = 0; i < n_cell; i++) {
          fscanf(f, "%d", &n_vtx_cell);
          if (n_vtx_cell == 3) {
            _n_face++;
            for (int j = 0; j < 3; j++) {
              fscanf(f, "%d", &_face_vtx[3*i+j]);
            }
          } else {
            int other_cell_vtx = 0;
            for (int j = 0; j < n_vtx_cell; j++) {
              fscanf(f, "%d", &other_cell_vtx);
            }
          }
        }

        _face_vtx = realloc(_face_vtx, 3 * _n_face * sizeof(int));

        log_trace("_n_face : %d\n", _n_face);
        PDM_log_trace_array_int(_face_vtx, 3*_n_face, "_face_vtx: ");
      } // end if CELLS

      if (strstr(line, "LOOKUP_TABLE default") != NULL) {
        // Malloc
        _vtx_field = malloc(sizeof(double) * _n_vtx);

        // Get connectivity
        for (int i = 0; i < _n_vtx; i++) {
          fscanf(f, "%lf", &_vtx_field[i]);
        }
        PDM_log_trace_array_double(_vtx_field, _n_vtx, "_vtx_field: ");
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
  char *filename  = NULL;

  _read_args (argc,
              argv,
              &verbose,
              &filename);

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
            filename,
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

  // Partition data

  // Set up

  // Mesh

  // Field

  // Exchange

  // free
  free(face_vtx );
  free(vtx_coord);
  free(vtx_field);

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
