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

#include "grid_mesh.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dmesh.h"
#include "pdm_array.h"
#include "pdm_logging.h"

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
  PDM_split_dual_t      *part_method,
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
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}












static void
_gen_mesh
(
 const PDM_MPI_Comm        comm,
 const int                 n_part,
 const PDM_split_dual_t    part_method,
 const PDM_g_num_t         n_vtx_seg,
 const int                 randomize,
 const int                 random_seed,
 int                     **pn_face,
 int                     **pn_vtx,
 int                    ***pface_vtx_idx,
 int                    ***pface_vtx,
 double                 ***pvtx_coord,
 PDM_g_num_t            ***pface_ln_to_gn,
 PDM_g_num_t            ***pvtx_ln_to_gn
 )
{
  /* Generate a distributed polygonal mesh */
  PDM_g_num_t  ng_face         = 0;
  PDM_g_num_t  ng_vtx          = 0;
  PDM_g_num_t  ng_edge         = 0;
  int          dn_vtx          = 0;
  double      *dvtx_coord      = NULL;
  int          dn_face         = 0;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  PDM_g_num_t *dface_edge      = NULL;
  int          dn_edge         = 0;
  PDM_g_num_t *dedge_vtx       = NULL;
  PDM_g_num_t *dedge_face      = NULL;
  int          n_edge_group    = 0;
  int         *dedge_group_idx = NULL;
  PDM_g_num_t *dedge_group     = NULL;

  PDM_poly_surf_gen(comm,
                    0.,
                    0.1,
                    0.,
                    0.1,
                    randomize,
                    random_seed,
                    n_vtx_seg,
                    n_vtx_seg,
                    &ng_face,
                    &ng_vtx,
                    &ng_edge,
                    &dn_vtx,
                    &dvtx_coord,
                    &dn_face,
                    &dface_vtx_idx,
                    &dface_vtx,
                    &dface_edge,
                    &dn_edge,
                    &dedge_vtx,
                    &dedge_face,
                    &n_edge_group,
                    &dedge_group_idx,
                    &dedge_group);

  /* Spit the mesh */
  int n_zone = 1;
  int *n_part_zones = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0] = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                n_part_zones,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);
  free(n_part_zones);

  int n_join = 0;
  PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                        dn_face,
                                        dn_edge,
                                        -1,
                                        dn_vtx,
                                        n_edge_group,
                                        n_join,
                                        comm);

  int *djoins_ids = malloc (sizeof(int) * n_join);
  int *dedge_join_idx = malloc (sizeof(int) * (n_join + 1));
  dedge_join_idx[0] = 0;
  PDM_g_num_t *dedge_join = malloc (sizeof(PDM_g_num_t) * dedge_join_idx[n_join]);

  int *dedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);

  PDM_dmesh_set(dmesh,
                dvtx_coord,
                dedge_vtx_idx,
                dedge_vtx,
                dedge_face,
                dedge_group_idx,
                dedge_group,
                djoins_ids,
                dedge_join_idx,
                dedge_join);

  PDM_multipart_register_block(mpart, 0, dmesh);

  /* Connection between zones */
  int n_total_joins = 0;
  int *join_to_opposite = malloc(sizeof(int) * n_total_joins);
  PDM_multipart_register_joins(mpart, n_total_joins, join_to_opposite);

  /* Run */
  PDM_multipart_run_ppart(mpart);

  free(djoins_ids);
  free(dedge_join_idx);
  free(dedge_join);
  free(join_to_opposite);



  /* Get partitioned mesh */
  *pn_face        = (int *)          malloc(sizeof(int *)          * n_part);
  *pn_vtx         = (int *)          malloc(sizeof(int *)          * n_part);
  *pface_vtx_idx  = (int **)         malloc(sizeof(int **)         * n_part);
  *pface_vtx      = (int **)         malloc(sizeof(int **)         * n_part);
  *pface_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_coord     = (double **)      malloc(sizeof(double **)      * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_face;
    int n_edge;
    int n_edge_part_bound;
    int n_vtx;
    int n_proc;
    int n_t_part;
    int s_face_edge;
    int s_edge_vtx;
    int s_edge_join;
    int s_edge_group;

    int n_groups, n_joins;
    int n_section;
    int *n_elt;

    int         *face_tag;
    int         *face_edge_idx;
    int         *face_edge;
    PDM_g_num_t *face_ln_to_gn;
    int         *edge_tag;
    int         *edge_face;
    int         *edge_vtx_idx;
    int         *edge_vtx;
    PDM_g_num_t *edge_ln_to_gn;
    int         *edge_part_bound_proc_idx;
    int         *edge_part_bound_part_idx;
    int         *edge_part_bound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *edge_group_idx;
    int         *edge_group;
    PDM_g_num_t *edge_group_ln_to_gn;
    PDM_g_num_t *edge_join_ln_to_gn;
    int         *edge_join_idx, *edge_join;
    int         **elt_vtx_idx;
    int         **elt_vtx;
    PDM_g_num_t **elt_section_ln_to_gn;

    PDM_multipart_part_dim_get(mpart,
                               0,
                               i_part,
                               &n_section,
                               &n_elt,
                               &n_face,
                               &n_edge,
                               &n_edge_part_bound,
                               &n_vtx,
                               &n_proc,
                               &n_t_part,
                               &s_face_edge,
                               &s_edge_vtx,
                               &s_edge_group,
                               &n_groups,
                               &s_edge_join,
                               &n_joins);

    PDM_multipart_part_val_get(mpart,
                               0,
                               i_part,
                               &elt_vtx_idx,
                               &elt_vtx,
                               &elt_section_ln_to_gn,
                               &face_tag,
                               &face_edge_idx,
                               &face_edge,
                               &face_ln_to_gn,
                               &edge_tag,
                               &edge_face,
                               &edge_vtx_idx,
                               &edge_vtx,
                               &edge_ln_to_gn,
                               &edge_part_bound_proc_idx,
                               &edge_part_bound_part_idx,
                               &edge_part_bound,
                               &vtx_tag,
                               &vtx,
                               &vtx_ln_to_gn,
                               &edge_group_idx,
                               &edge_group,
                               &edge_group_ln_to_gn,
                               &edge_join_idx,
                               &edge_join,
                               &edge_join_ln_to_gn);

    *(pn_face)[i_part] = n_face;
    *(pn_vtx)[i_part]  = n_vtx;

    /* Vertices */
    (*pvtx_coord)[i_part] = (double *) malloc(sizeof(double) * 3 * n_vtx);
    memcpy((*pvtx_coord)[i_part], vtx, sizeof(double) * 3 * n_vtx);

    (*pvtx_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_vtx);
    memcpy((*pvtx_ln_to_gn)[i_part], vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_vtx);


    /* Faces */
    (*pface_vtx_idx)[i_part] = (int *) malloc(sizeof(int) * (n_face + 1));
    memcpy((*pface_vtx_idx)[i_part], face_edge_idx, sizeof(int) * (n_face + 1));

    PDM_compute_face_vtx_from_face_and_edge(n_face,
                                            face_edge_idx,
                                            face_edge,
                                            edge_vtx,
                                            *pface_vtx + i_part);

    (*pface_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_face);
    memcpy((*pface_ln_to_gn)[i_part], face_ln_to_gn, sizeof(PDM_g_num_t) * n_face);
  }
  PDM_multipart_free(mpart);
  PDM_dmesh_free(dmesh);

  free(dvtx_coord);
  free(dface_vtx_idx);
  free(dface_vtx);
  free(dface_edge);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face);
  free(dedge_group_idx);
  free(dedge_group);
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

#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#else
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
#endif
#endif


  _read_args(argc,
             argv,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &part_method,
             &tolerance,
             &randomize);

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
  CWP_Status_t *is_active_rank = malloc(sizeof(CWP_Status_t) * n_code);
  double *time_init = malloc(sizeof(double) * n_code);

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
  is_active_rank[0] = CWP_STATUS_ON;
  time_init[0] = 0.;

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           time_init,
           intra_comm);


  if (rank == 0) {
    printf("CWIPI Init OK\n");
  }


  // Create coupling
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1_dynamic";
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Cpl_create(code_name[i_code],                                     // Code name
                   cpl_name,                                              // Coupling id
                   coupled_code_name[i_code],                             // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,                                // Coupling type
                   CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                   n_part,                                                // Partition number
                   CWP_DYNAMIC_MESH_DEFORMABLE,                           // Mesh displacement type
                   CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency
  }


  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Visu_set(code_name[i_code],       // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
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
  PDM_g_num_t ***pface_ln_to_gn = (PDM_g_num_t ***) malloc(sizeof(PDM_g_num_t **) * n_code);
  PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc(sizeof(PDM_g_num_t **) * n_code);


  for (int i_code = 0 ; i_code < n_code ; i_code++) {

    PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) intra_comm);

    _gen_mesh(mesh_comm,
              n_part,
              part_method,
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
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Mesh_interf_finalize(code_name[i_code], cpl_name);
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

  for (int step = 0; step < 2; step++) {

    recv_time += 1.;

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

      CWP_next_recv_time_set(code_name[i_code],
                             cpl_name,
                             recv_time);

      CWP_Spatial_interp_weights_compute(code_name[i_code], cpl_name);
    }


    MPI_Barrier(MPI_COMM_WORLD);


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
  free(is_active_rank);
  free(intra_comm);
  free(time_init);

  //  Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
