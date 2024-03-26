/*
  This file is part of the CWIPI library.

  Copyright (C) 2022-2023  ONERA

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

#include "cwipi.h"
#include "cwp.h"
#include "cwp_priv.h"

#define ABS(a)    ((a) < 0   ? -(a) : (a))
#define MIN(a, b) ((a) < (b) ?  (a) : (b))
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
                                    (CWPT_Mesh_nodal_elt_t) CWP_BLOCK_FACE_POLY,
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

  _read_args(argc,
             argv,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &n_part1,
             &n_part2,
             &tolerance,
             &randomize);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  assert(n_rank > 1);

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

  int n_code = 0;
  int n_vtx_seg[2];
  int n_part   [2];
  int code_id  [2];
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      n_vtx_seg        [n_code] = all_n_vtx_seg [icode];
      n_part           [n_code] = all_n_part    [icode];
      n_code++;
    }
  }

  // log_trace("n_code = %d : ", n_code);
  // for (int i_code = 0; i_code < n_code; i_code++) {
  //   log_trace("%s ", code_name[i_code]);
  // }
  // log_trace("\n");

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
  const char *cpl_name = "c_new_api_surf_P0P0";
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;

  for (int i_code = 0; i_code < n_code; i_code++) {
    CWP_Cpl_create(code_name[i_code],                                     // Code name
                   cpl_name,                                              // Coupling id
                   coupled_code_name[i_code],                             // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,                                // Coupling type
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
              CWPT_SPLIT_DUAL_WITH_HILBERT,
              n_vtx_seg[i_code],
              randomize,
              code_id        [i_code],
              &pn_face       [i_code],
              &pn_vtx        [i_code],
              &pface_vtx_idx [i_code],
              &pface_vtx     [i_code],
              &pvtx_coord    [i_code],
              &pface_ln_to_gn[i_code],
              &pvtx_ln_to_gn [i_code]);

    // int block_id = CWP_Mesh_interf_block_add(code_name[i_code],
    //                                          cpl_name,
    //                                          CWP_BLOCK_FACE_POLY);
    for (int i = 0; i < n_part[i_code]; i++) {
      CWP_Mesh_interf_vtx_set(code_name[i_code],
                              cpl_name,
                              i,
                              pn_vtx       [i_code][i],
                              pvtx_coord   [i_code][i],
                              pvtx_ln_to_gn[i_code][i]);


      // CWP_Mesh_interf_f_poly_block_set(code_name[i_code],
      //                                  cpl_name,
      //                                  i,
      //                                  block_id,
      //                                  pn_face       [i_code][i],
      //                                  pface_vtx_idx [i_code][i],
      //                                  pface_vtx     [i_code][i],
      //                                  pface_ln_to_gn[i_code][i]);
      CWP_Mesh_interf_from_facevtx_set(code_name[i_code],
                                       cpl_name,
                                       i,
                                       pn_face       [i_code][i],
                                       pface_vtx_idx [i_code][i],
                                       pface_vtx     [i_code][i],
                                       NULL);//pface_ln_to_gn[i_code][i]);
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
    }
    else {
      CWP_Field_irecv (code_name[i_code], cpl_name, field_name1);
    }


    if (code_id[i_code] == 1) {
      CWP_Field_wait_issend(code_name[i_code], cpl_name, field_name1);
    }
    else {
      CWP_Field_wait_irecv (code_name[i_code], cpl_name, field_name1);
    }
  }

  if (i_rank == 0) {
    printf("Exchange fields OK\n");
  }


  int check = EXIT_SUCCESS;
  if (spatial_interp == CWP_SPATIAL_INTERP_FROM_INTERSECTION) {

    double l_mass[2] = {0., 0.};
    double l_area[2] = {0., 0.};

    for (int icode = 0; icode < n_code; icode++) {
      if (is_active_rank == CWP_STATUS_ON) {
        for (int ipart = 0; ipart < n_part[icode]; ipart++) {

          double *field_val = NULL;
          if (code_id[icode] == 1) {
            field_val = send_val[icode][ipart];
          }
          else {
            field_val = recv_val[icode][ipart];
          }

          // Compute polygon area and mass
          int    n_vtx_face = 0;
          double area       = 0.;
          double face_center[3];
          double edge_center[3];
          double v1v2[3];
          double vectFECenter[3];
          for (int i_face = 0; i_face < pn_face[icode][ipart]; i_face++) {

            area = 0;
            // Compute face center
            for (int j = 0; j < 3; j++) {
              face_center[j]  = 0.;
            }

            n_vtx_face = pface_vtx_idx[icode][ipart][i_face+1]-pface_vtx_idx[icode][ipart][i_face];
            for (int i_vtx = 0; i_vtx < n_vtx_face; i_vtx++) {
              int vtx_id = pface_vtx[icode][ipart][pface_vtx_idx[icode][ipart][i_face] + i_vtx]-1;

              for (int j = 0; j < 3; j++) {
                face_center[j] += pvtx_coord[icode][ipart][3*vtx_id + j];
              }
            }

            for (int j = 0; j < 3; j++) {
              face_center[j] /= n_vtx_face;
            }

            // Compute area
            for (int i_vtx = 0; i_vtx < n_vtx_face; i_vtx++) {
              int vtx_id = pface_vtx[icode][ipart][pface_vtx_idx[icode][ipart][i_face] + i_vtx]-1;
              int vtx_jd = pface_vtx[icode][ipart][pface_vtx_idx[icode][ipart][i_face] + (i_vtx + 1) % n_vtx_face] - 1;

              for (int j = 0; j < 3; j++) {
                edge_center[j] = 0.5 * (pvtx_coord[icode][ipart][3*vtx_jd + j] + pvtx_coord[icode][ipart][3*vtx_id + j]);
                v1v2[j]        = pvtx_coord[icode][ipart][3*vtx_jd + j] - pvtx_coord[icode][ipart][3*vtx_id + j];
              }

              for (int j = 0; j < 3; j++) {
                vectFECenter[j] = edge_center[j] - face_center[j];
              }

              double surface_vectorTria[3];
              surface_vectorTria[0] = vectFECenter[1] * v1v2[2] - v1v2[1] * vectFECenter[2];
              surface_vectorTria[1] = v1v2[0] * vectFECenter[2] - vectFECenter[0] * v1v2[2];
              surface_vectorTria[2] = vectFECenter[0] * v1v2[1] - v1v2[0] * vectFECenter[1];

              for (int i = 0; i < 3; i++)
                surface_vectorTria[i] *= 0.5;

              const double areaTri = sqrt(surface_vectorTria[0] * surface_vectorTria[0] + surface_vectorTria[1] * surface_vectorTria[1] + surface_vectorTria[2] * surface_vectorTria[2]);

              area += areaTri;
            }

            l_area[code_id[icode]-1] += area;
            l_mass[code_id[icode]-1] += field_val[i_face] * area;
          }
        }
      }
    }


    double g_mass[2];
    MPI_Allreduce(l_mass, g_mass, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double g_area[2];
    MPI_Allreduce(l_area, g_area, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double rel_error_area = ABS(g_area[0] - g_area[1])/ABS(g_area[1]);
    double rel_error_mass = ABS(g_mass[0] - g_mass[1])/ABS(g_mass[1]);

    if (i_rank == 0) {
      printf("g_area = %20.16e / %20.16e, relative diff = %e\n",
             g_area[0], g_area[1], rel_error_area);

      printf("g_mass   = %20.16e / %20.16e, relative diff = %e\n",
             g_mass[0], g_mass[1], rel_error_mass);
    }

    if (rel_error_area > 1e-12 || rel_error_mass > 1e-12) {
      check = EXIT_FAILURE;
    }
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

  //  Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return check;
}

