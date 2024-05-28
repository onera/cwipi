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
  int                   *randomize,
  int                   *variable_mesh
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
    else if (strcmp(argv[i], "-variable_mesh") == 0) {
      *variable_mesh = 1;
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
  int    variable_mesh         = 0;

  _read_args(argc,
             argv,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &tolerance,
             &randomize,
             &variable_mesh);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  assert (comm_world_size > 1);


  // Initialize CWIPI
  int n_part = 1;
  int n_code = 1;
  int code_id[2];
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  int n_vtx_seg;
  if (rank < comm_world_size / 2) {
    code_id[0] = 1;
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
    n_vtx_seg = n_vtx_seg1;
  }
  else {
    code_id[0] = 2;
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
    n_vtx_seg = n_vtx_seg2;
  }


  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);


  if (rank == 0) {
    printf("CWIPI Init OK\n");
  }


  // Create coupling
  CWP_Dynamic_mesh_t displacement = CWP_DYNAMIC_MESH_DEFORMABLE;
  if (variable_mesh) {
    displacement = CWP_DYNAMIC_MESH_VARIABLE;
  }
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1_dynamic";
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Cpl_create(code_name[i_code],              // Code name
                   cpl_name,                       // Coupling id
                   coupled_code_name[i_code],      // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,         // Coupling type
                   spatial_interp,
                   n_part,                         // Partition number
                   displacement,                   // Mesh displacement type
                   CWP_TIME_EXCH_USER_CONTROLLED); // Postprocessing frequency
  }


  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Visu_set(code_name[i_code],       // Code name
                 cpl_name,                // Coupling id
                 2,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

  if (rank == 0) {
    printf("Create coupling OK\n");
  }



  // Mesh definition
  int          **pn_face        = (int          **) malloc(sizeof(int          *) * n_code);
  int          **pn_vtx         = (int          **) malloc(sizeof(int          *) * n_code);
  int         ***pface_vtx_idx  = (int         ***) malloc(sizeof(int         **) * n_code);
  int         ***pface_vtx      = (int         ***) malloc(sizeof(int         **) * n_code);
  double      ***pvtx_coord     = (double      ***) malloc(sizeof(double      **) * n_code);
  CWP_g_num_t ***pface_ln_to_gn = (CWP_g_num_t ***) malloc(sizeof(CWP_g_num_t **) * n_code);
  CWP_g_num_t ***pvtx_ln_to_gn  = (CWP_g_num_t ***) malloc(sizeof(CWP_g_num_t **) * n_code);


  for (int i_code = 0 ; i_code < n_code ; i_code++) {

    _gen_mesh(intra_comm[i_code],
              n_part,
              CWPT_SPLIT_DUAL_WITH_HILBERT,
              n_vtx_seg,
              randomize,
              code_id[i_code],
              &pn_face[i_code],
              &pn_vtx[i_code],
              &pface_vtx_idx[i_code],
              &pface_vtx[i_code],
              &pvtx_coord[i_code],
              &pface_ln_to_gn[i_code],
              &pvtx_ln_to_gn[i_code]);


    if (!variable_mesh) {
      CWP_Mesh_interf_del(code_name[i_code],
                          cpl_name);
      CWP_Mesh_interf_vtx_set(code_name[i_code],
                              cpl_name,
                              0,
                              pn_vtx[i_code][0],
                              pvtx_coord[i_code][0],
                              pvtx_ln_to_gn[i_code][0]);

      int block_id = CWP_Mesh_interf_block_add(code_name[i_code],
                                               cpl_name,
                                               CWP_BLOCK_FACE_POLY);

      CWP_Mesh_interf_f_poly_block_set(code_name[i_code],
                                       cpl_name,
                                       0,
                                       block_id,
                                       pn_face[i_code][0],
                                       pface_vtx_idx[i_code][0],
                                       pface_vtx[i_code][0],
                                       pface_ln_to_gn[i_code][0]);

      CWP_Mesh_interf_finalize(code_name[i_code], cpl_name);
    }
  }

  if (rank == 0) {
    printf("Set mesh OK\n");
  }

  // Create and set fields
  // field1: code1 -> code2
  // field2: code2 -> code1
  const char *field_name1 = "cooX_t0";
  const char *field_name2 = "code2_elt_gnum";
  double **send_val = (double **) malloc(sizeof(double *) * n_code);
  double **recv_val = (double **) malloc(sizeof(double *) * n_code);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (code_id[i_code] == 1) {
      send_val[i_code] = (double *) malloc(sizeof(double) * 3 * pn_vtx[i_code][0]);
      recv_val[i_code] = (double *) malloc(sizeof(double) * pn_face[i_code][0]);
      for (int i = 0 ; i < pn_vtx[i_code][0] ; i++) {
        send_val[i_code][3*i  ] = pvtx_coord[i_code][0][3*i  ];
        send_val[i_code][3*i+1] = pvtx_coord[i_code][0][3*i+1];
        send_val[i_code][3*i+2] = pvtx_coord[i_code][0][3*i+2];
      }
    }
    else {
      send_val[i_code] = (double *) malloc(sizeof(double) * pn_face[i_code][0]);
      recv_val[i_code] = (double *) malloc(sizeof(double) * 3 * pn_vtx[i_code][0]);
      for (int i = 0 ; i < pn_face[i_code][0] ; i++) {
        send_val[i_code][i] = (double) pface_ln_to_gn[i_code][0][i];
      }
    }
  }

  CWP_Status_t visu_status = CWP_STATUS_ON;

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (code_id[i_code] == 1) {
      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       3,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         field_name1,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val[i_code]);

      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val[i_code]);
    }
    else {
      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       3,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         field_name1,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val[i_code]);

      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val[i_code]);
    }
  }





  double recv_time = 0.;
  for (int step = 0; step < 10; step++) {

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
      printf("  Step %d\n", step);
    }

    // Mesh rotation and new localisation
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      if (code_id[i_code] == 1) {
        mesh_rotate(pvtx_coord[i_code][0], pn_vtx[i_code][0], recv_time);
      } else {
        mesh_rotate(pvtx_coord[i_code][0], pn_vtx[i_code][0], 3 * recv_time);
      }

      // Start time step
      CWP_Time_step_beg(code_name[i_code],
                        recv_time);

      if (variable_mesh && step%2 == 0) {
        CWP_Mesh_interf_del(code_name[i_code],
                            cpl_name);

        CWP_Mesh_interf_vtx_set(code_name[i_code],
                                cpl_name,
                                0,
                                pn_vtx[i_code][0],
                                pvtx_coord[i_code][0],
                                pvtx_ln_to_gn[i_code][0]);

        int block_id = CWP_Mesh_interf_block_add(code_name[i_code],
                                                 cpl_name,
                                                 CWP_BLOCK_FACE_POLY);

        CWP_Mesh_interf_f_poly_block_set(code_name[i_code],
                                         cpl_name,
                                         0,
                                         block_id,
                                         pn_face[i_code][0],
                                         pface_vtx_idx[i_code][0],
                                         pface_vtx[i_code][0],
                                         pface_ln_to_gn[i_code][0]);

        CWP_Mesh_interf_finalize(code_name[i_code], cpl_name);
      }
    }

    // Separate loops to avoid deadlock if multiple codes on same MPI rank
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      CWP_Spatial_interp_property_set(code_name[i_code],
                                      cpl_name,
                                      "n_neighbors",
                                      CWP_INT,
                                      "1");

      CWP_Spatial_interp_weights_compute(code_name[i_code], cpl_name);
    }

    MPI_Barrier(MPI_COMM_WORLD);


    if (step%3 == 0) {
      for (int i_code = 0 ; i_code < n_code ; i_code++) {
        if (code_id[i_code] == 1) {
          CWP_Field_issend(code_name[i_code], cpl_name, field_name1);
          CWP_Field_irecv (code_name[i_code], cpl_name, field_name2);
        }
        else {
          CWP_Field_irecv (code_name[i_code], cpl_name, field_name1);
          CWP_Field_issend(code_name[i_code], cpl_name, field_name2);
        }
      }


      for (int i_code = 0 ; i_code < n_code ; i_code++) {
        if (code_id[i_code] == 1) {
          CWP_Field_wait_issend(code_name[i_code], cpl_name, field_name1);
          CWP_Field_wait_irecv (code_name[i_code], cpl_name, field_name2);
        }
        else {
          CWP_Field_wait_irecv (code_name[i_code], cpl_name, field_name1);
          CWP_Field_wait_issend(code_name[i_code], cpl_name, field_name2);
        }
      }
    }

    // End time step
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      CWP_Time_step_end(code_name[i_code]);
    }

    // Increase
    recv_time += 1.;

  }


  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Mesh_interf_del(code_name[i_code], cpl_name);
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Cpl_del(code_name[i_code], cpl_name);
  }



  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    for (int i_part = 0 ; i_part < n_part ; i_part++) {
      free(pface_vtx_idx[i_code][i_part]);
      free(pface_vtx[i_code][i_part]);
      free(pvtx_coord[i_code][i_part]);
      free(pface_ln_to_gn[i_code][i_part]);
      free(pvtx_ln_to_gn[i_code][i_part]);
    }
    free(pn_face[i_code]);
    free(pn_vtx[i_code]);
    free(pface_vtx_idx[i_code]);
    free(pface_vtx[i_code]);
    free(pvtx_coord[i_code]);
    free(pface_ln_to_gn[i_code]);
    free(pvtx_ln_to_gn[i_code]);

    free(send_val[i_code]);
    free(recv_val[i_code]);
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
