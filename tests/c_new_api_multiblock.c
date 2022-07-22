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

#include "pdm_poly_surf_gen.h"
#include "pdm_part.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

#include "pdm_multipart.h"
#include "pdm_dcube_gen.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include "pdm_array.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"

#include "pdm_part_extension.h"
#include "pdm_vtk.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_poly_vol_gen.h"


#define ABS(a) ((a) <  0  ? -(a) : (a))

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
  double                *length,
  double                *separation_x,
  double                *separation_y,
  double                *separation_z,
  int                   *deform,
  double                *tolerance,
  int                   *randomize,
  int                   *nProcData,
  PDM_split_dual_t      *part_method,
  char                 **output_filename,
  int                   *verbose,
  int                   *use_gnum
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
    else if (strcmp(argv[i], "-length") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *length = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-sep") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepy") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_y = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_z = atof(argv[i]);
    }

    else if (strcmp(argv[i], "-output") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *output_filename = argv[i];
      }
    }
    else if (strcmp(argv[i], "-def") == 0) {
      *deform = 1;
    }
    else if (strcmp(argv[i], "-no_random") == 0) {
      *randomize = 0;
    }
    else if (strcmp(argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *nProcData = atoi(argv[i]);
      }
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
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else if (strcmp(argv[i], "-no_gnum") == 0) {
      *use_gnum = 0;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_gen_mesh
(
 const int                 active_rank,
 const PDM_MPI_Comm        comm,
 const int                 n_part,
 const PDM_split_dual_t    part_method,
 const PDM_g_num_t         n,
 const double              xmin,
 const double              ymin,
 const double              zmin,
 const double              length,
 const int                 randomize,
 const int                 random_seed,
 int                     **pn_cell,
 int                     **pn_face,
 int                     **pn_vtx,
 int                    ***pcell_face_idx,
 int                    ***pcell_face,
 int                    ***pface_vtx_idx,
 int                    ***pface_vtx,
 double                 ***pvtx_coord,
 PDM_g_num_t            ***pcell_ln_to_gn,
 PDM_g_num_t            ***pface_ln_to_gn,
 PDM_g_num_t            ***pvtx_ln_to_gn
 )
{
  PDM_multipart_t *mpart = NULL;
  PDM_dmesh_t     *dmesh = NULL;

  PDM_g_num_t  ng_cell;
  PDM_g_num_t  ng_face;
  PDM_g_num_t  ng_vtx;
  int          n_face_group;
  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  int         *dcell_face_idx  = NULL;
  PDM_g_num_t *dcell_face      = NULL;
  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  double      *dvtx_coord      = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  if (active_rank) {

    PDM_poly_vol_gen(comm,
                     xmin,
                     ymin,
                     zmin,
                     length,
                     length,
                     length,
                     n,
                     n,
                     n,
                     randomize,
                     random_seed,
                     &ng_cell,
                     &ng_face,
                     &ng_vtx,
                     &n_face_group,
                     &dn_cell,
                     &dn_face,
                     &dn_vtx,
                     &dcell_face_idx,
                     &dcell_face,
                     &dface_cell,
                     &dface_vtx_idx,
                     &dface_vtx,
                     &dvtx_coord,
                     &dface_group_idx,
                     &dface_group);


    mpart = PDM_multipart_create(1,
                                 &n_part,
                                 PDM_FALSE,
                                 part_method,
                                 PDM_PART_SIZE_HOMOGENEOUS,
                                 NULL,
                                 comm,
                                 PDM_OWNERSHIP_KEEP);

    /* Generate dmesh */
    int n_join = 0;
    dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                             dn_cell,
                             dn_face,
                             -1, // dn_edge
                             dn_vtx,
                             n_face_group,
                             n_join,
                             comm);

    int *djoins_ids = malloc (sizeof(int) * n_join);
    int *dface_join_idx = malloc (sizeof(int) * (n_join + 1));
    dface_join_idx[0] = 0;
    PDM_g_num_t *dface_join = malloc (sizeof(PDM_g_num_t) * dface_join_idx[n_join]);

    PDM_dmesh_set(dmesh,
                  dvtx_coord,
                  dface_vtx_idx,
                  dface_vtx,
                  dface_cell,
                  dface_group_idx,
                  dface_group,
                  djoins_ids,
                  dface_join_idx,
                  dface_join);

    PDM_multipart_register_block(mpart, 0, dmesh);

    /* Connection between zones */
    int n_total_joins = 0;
    int *join_to_opposite = malloc(sizeof(int) * n_total_joins);
    PDM_multipart_register_joins(mpart, n_total_joins, join_to_opposite);

    /* Run */
    PDM_multipart_run_ppart(mpart);

    free(djoins_ids);
    free(dface_join_idx);
    free(dface_join);
    free(join_to_opposite);
  } // end if (active_rank)


  *pn_cell        = (int *)          malloc(sizeof(int *)          * n_part);
  *pn_face        = (int *)          malloc(sizeof(int *)          * n_part);
  *pn_vtx         = (int *)          malloc(sizeof(int *)          * n_part);
  *pcell_face_idx = (int **)         malloc(sizeof(int **)         * n_part);
  *pcell_face     = (int **)         malloc(sizeof(int **)         * n_part);
  *pface_vtx_idx  = (int **)         malloc(sizeof(int **)         * n_part);
  *pface_vtx      = (int **)         malloc(sizeof(int **)         * n_part);
  *pcell_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pface_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_coord     = (double **)      malloc(sizeof(double **)      * n_part);

  if (active_rank) {

    for (int i_part = 0; i_part < n_part; i_part++) {
      int n_cell;
      int n_face;
      int n_face_part_bound;
      int n_vtx;
      int n_proc;
      int n_t_part;
      int s_cell_face;
      int s_face_vtx;
      int s_face_join;
      int s_face_group;

      int n_groups, n_joins;
      int n_section;
      int *n_elt;

      int         *cell_tag;
      int         *cell_face_idx;
      int         *cell_face;
      PDM_g_num_t *cell_ln_to_gn;
      int         *face_tag;
      int         *face_cell;
      int         *face_vtx_idx;
      int         *face_vtx;
      PDM_g_num_t *face_ln_to_gn;
      int         *face_part_bound_proc_idx;
      int         *face_part_bound_part_idx;
      int         *face_part_bound;
      int         *vtx_tag;
      double      *vtx;
      PDM_g_num_t *vtx_ln_to_gn;
      int         *face_group_idx;
      int         *face_group;
      PDM_g_num_t *face_group_ln_to_gn;
      PDM_g_num_t *face_join_ln_to_gn;
      int         *face_join_idx, *face_join;
      int         **elt_vtx_idx;
      int         **elt_vtx;
      PDM_g_num_t **elt_section_ln_to_gn;

      PDM_multipart_part_dim_get(mpart,
                                 0,
                                 i_part,
                                 &n_section,
                                 &n_elt,
                                 &n_cell,
                                 &n_face,
                                 &n_face_part_bound,
                                 &n_vtx,
                                 &n_proc,
                                 &n_t_part,
                                 &s_cell_face,
                                 &s_face_vtx,
                                 &s_face_group,
                                 &n_groups,
                                 &s_face_join,
                                 &n_joins);

      PDM_multipart_part_val_get(mpart,
                                 0,
                                 i_part,
                                 &elt_vtx_idx,
                                 &elt_vtx,
                                 &elt_section_ln_to_gn,
                                 &cell_tag,
                                 &cell_face_idx,
                                 &cell_face,
                                 &cell_ln_to_gn,
                                 &face_tag,
                                 &face_cell,
                                 &face_vtx_idx,
                                 &face_vtx,
                                 &face_ln_to_gn,
                                 &face_part_bound_proc_idx,
                                 &face_part_bound_part_idx,
                                 &face_part_bound,
                                 &vtx_tag,
                                 &vtx,
                                 &vtx_ln_to_gn,
                                 &face_group_idx,
                                 &face_group,
                                 &face_group_ln_to_gn,
                                 &face_join_idx,
                                 &face_join,
                                 &face_join_ln_to_gn);

      *(pn_cell)[i_part] = n_cell;
      *(pn_face)[i_part] = n_face;
      *(pn_vtx)[i_part]  = n_vtx;

      /* Vertices */
      (*pvtx_coord)[i_part] = (double *) malloc(sizeof(double) * 3 * n_vtx);
      memcpy((*pvtx_coord)[i_part], vtx, sizeof(double) * 3 * n_vtx);

      (*pvtx_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_vtx);
      memcpy((*pvtx_ln_to_gn)[i_part], vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_vtx);


      /* Cells */
      (*pcell_face_idx)[i_part] = (int *) malloc(sizeof(int) * (n_cell + 1));
      memcpy((*pcell_face_idx)[i_part], cell_face_idx, sizeof(int) * (n_cell + 1));

      s_cell_face = cell_face_idx[n_cell];
      (*pcell_face)[i_part] = (int *) malloc(sizeof(int) * s_cell_face);
      memcpy((*pcell_face)[i_part], cell_face, sizeof(int) * cell_face_idx[n_cell]);

      (*pcell_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_cell);
      memcpy((*pcell_ln_to_gn)[i_part], cell_ln_to_gn, sizeof(PDM_g_num_t) * n_cell);


      /* Faces */
      (*pface_vtx_idx)[i_part] = (int *) malloc(sizeof(int) * (n_face + 1));
      memcpy((*pface_vtx_idx)[i_part], face_vtx_idx, sizeof(int) * (n_face + 1));

      s_face_vtx = face_vtx_idx[n_face];
      (*pface_vtx)[i_part] = (int *) malloc(sizeof(int) * s_face_vtx);
      memcpy((*pface_vtx)[i_part], face_vtx, sizeof(int) * face_vtx_idx[n_face]);

      (*pface_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_face);
      memcpy((*pface_ln_to_gn)[i_part], face_ln_to_gn, sizeof(PDM_g_num_t) * n_face);
    }

    PDM_multipart_free(mpart);
    PDM_dmesh_free(dmesh);

    free(dcell_face_idx);
    free(dcell_face);
    free(dface_cell);
    free(dface_vtx_idx);
    free(dface_vtx);
    free(dvtx_coord);
    free(dface_group_idx);
    free(dface_group);
  }

  else {

    for (int i_part = 0; i_part < n_part; i_part++) {
      int n_cell = 0;
      int n_face = 0;
      int n_vtx  = 0;

      *(pn_cell)[i_part] = n_cell;
      *(pn_face)[i_part] = n_face;
      *(pn_vtx)[i_part]  = n_vtx;

      /* Vertices */
      (*pvtx_coord)[i_part] = (double *) malloc(sizeof(double) * 3 * n_vtx);

      (*pvtx_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_vtx);


      /* Cells */
      (*pcell_face_idx)[i_part] = (int *) malloc(sizeof(int) * (n_cell + 1));
      (*pcell_face_idx)[i_part][0] = 0;

      int s_cell_face = (*pcell_face_idx)[i_part][n_cell];
      (*pcell_face)[i_part] = (int *) malloc(sizeof(int) * s_cell_face);

      (*pcell_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_cell);


      /* Faces */
      (*pface_vtx_idx)[i_part] = (int *) malloc(sizeof(int) * (n_face + 1));
      (*pface_vtx_idx)[i_part][0] = 0;

      int s_face_vtx = (*pface_vtx_idx)[i_part][n_face];
      (*pface_vtx)[i_part] = (int *) malloc(sizeof(int) * s_face_vtx);

      (*pface_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_face);
    }

  }


}



/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  // Read args from command line
  int n_vtx_seg1                  = 4;
  int n_vtx_seg2                  = 4;
  int randomize                   = 1;
  int n_proc_data                 = -1;

#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#else
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_HILBERT;
#endif
#endif
  int verbose                     = 0;

  double length                   = 20.;
  int    deform                   = 0;

  double      separation_x        = 2.;
  double      separation_y        = 0.;
  double      separation_z        = 0.;

  double      tolerance           = 1e-2;

  char* output_filename           = NULL;
  int filedump                    = 0;

  int         use_gnum            = 1;

  _read_args(argc,
             argv,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &length,
             &separation_x,
             &separation_y,
             &separation_z,
             &deform,
             &tolerance,
             &randomize,
             &n_proc_data,
             &part_method,
             &output_filename,
             &verbose,
             &use_gnum);

  if (output_filename !=NULL) {
    filedump = 1;
  }

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  assert (comm_world_size > 1);

  if (n_proc_data == 1) {
    n_proc_data = 2;
  }




  // Initialize CWIPI
  int n_part = 1;
  int n_code = 1;
  int code_id;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t *is_active_rank = malloc(sizeof(CWP_Status_t) * n_code);
  double *time_init = malloc(sizeof(double) * n_code);

  int n_vtx_seg;
  if (rank < comm_world_size / 2) {
    code_id = 1;
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
    n_vtx_seg = n_vtx_seg1;
  }
  else {
    code_id = 2;
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


  if (verbose && rank == 0) {
    printf("CWIPI Init OK\n");
  }


  // Create coupling
  const char *coupling_name = "c_surf_cpl_P1P1";

  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_VOLUME,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  CWP_Visu_set(code_name[0], coupling_name, 1, CWP_VISU_FORMAT_ENSIGHT, "text");

  if (verbose && rank == 0) {
    printf("Create coupling OK\n");
  }

  // Define mesh
  PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) intra_comm);

  int _n_proc_data = n_proc_data;
  if (n_proc_data > 0) {
    if (code_id == 1) {
      _n_proc_data /= 2;
    }
    else {
      _n_proc_data -= n_proc_data / 2;
    }
  }
  int current_rank_has_mesh = 1;//_set_rank_has_mesh(intra_comm[0], _n_proc_data, &mesh_comm);

  // int true_n_proc_data;
  // MPI_Reduce(&current_rank_has_mesh, &true_n_proc_data, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  // if (rank == 0) {
  //   printf("nb procs with mesh data = %d\n", true_n_proc_data);
  // }


  double xmin = -0.5 * length;
  double ymin = -0.5 * length;
  double zmin = -0.5 * length;
  //  int init_random = (int) time(NULL);
  int init_random = 5;

  PDM_MPI_Comm code_mesh_comm;
  PDM_MPI_Comm_split(mesh_comm, code_id, rank, &code_mesh_comm);

  if (code_id == 2) {
    init_random++;
    xmin += separation_x;
    ymin += separation_y;
    zmin += separation_z;
  }


  int          *pn_cell        = NULL;
  int          *pn_face        = NULL;
  int          *pn_vtx         = NULL;
  int         **pcell_face_idx = NULL;
  int         **pcell_face     = NULL;
  int         **pface_vtx_idx  = NULL;
  int         **pface_vtx      = NULL;
  double      **pvtx_coord     = NULL;
  PDM_g_num_t **pcell_ln_to_gn = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;
  PDM_g_num_t **pvtx_ln_to_gn  = NULL;

  _gen_mesh(current_rank_has_mesh,
            code_mesh_comm,
            n_part,
            part_method,
            n_vtx_seg,
            xmin,
            ymin,
            zmin,
            length,
            randomize,
            init_random,
            &pn_cell,
            &pn_face,
            &pn_vtx,
            &pcell_face_idx,
            &pcell_face,
            &pface_vtx_idx,
            &pface_vtx,
            &pvtx_coord,
            &pcell_ln_to_gn,
            &pface_ln_to_gn,
            &pvtx_ln_to_gn);


  // Set interface mesh
  PDM_g_num_t *_pvtx_ln_to_gn  = NULL;
  PDM_g_num_t *_pcell_ln_to_gn = NULL;
  if (use_gnum) {
    _pvtx_ln_to_gn  = pvtx_ln_to_gn[0];
    _pcell_ln_to_gn = pcell_ln_to_gn[0];
  }

  CWP_Mesh_interf_vtx_set(code_name[0],
                          coupling_name,
                          0,
                          pn_vtx[0],
                          pvtx_coord[0],
                          _pvtx_ln_to_gn);

  CWP_Mesh_interf_from_cellface_set(code_name[0],
                                    coupling_name,
                                    0,
                                    pn_cell[0],
                                    pcell_face_idx[0],
                                    pcell_face[0],
                                    pn_face[0],
                                    pface_vtx_idx[0],
                                    pface_vtx[0],
                                    _pcell_ln_to_gn);

  CWP_Mesh_interf_finalize(code_name[0],
                           coupling_name);

  if (verbose && rank == 0) {
    printf("Set mesh OK\n");
  }


  // Create and set fields
  double *send_val = NULL;
  double *recv_val = NULL;

  const char *field_name  = "coo";

  if (code_id == 1) {
    send_val = (double *) malloc(sizeof(double) * pn_vtx[0] * 3);
    memcpy(send_val, pvtx_coord[0], sizeof(double) * pn_vtx[0] * 3);
  } else {
    recv_val = (double *) malloc(sizeof(double) * pn_vtx[0] * 3);
  }

  CWP_Status_t visu_status = CWP_STATUS_ON;
  MPI_Barrier(MPI_COMM_WORLD);

  if (code_id == 1) {
    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     3,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);
    CWP_Field_data_set(code_name[0],
                       coupling_name,
                       field_name,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_val);
  } else {
    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     3,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);
    CWP_Field_data_set(code_name[0],
                       coupling_name,
                       field_name,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_val);
  }

  if (verbose && rank == 0) {
    printf("Fields OK\n");
  }

  // Perform geometric algorithm
  CWP_Spatial_interp_property_set(code_name[0], coupling_name, "tolerance", "double", "1e-2");
  CWP_Spatial_interp_weights_compute(code_name[0], coupling_name);

  int n_unlocated = 0;
  int n_located = 0;
  const int *located = NULL;
  if (code_id != 1) {
    n_unlocated = CWP_N_uncomputed_tgts_get(code_name[0], coupling_name, field_name, 0);
    n_located   = CWP_N_computed_tgts_get  (code_name[0], coupling_name, field_name, 0);
    located     = CWP_Computed_tgts_get    (code_name[0], coupling_name, field_name, 0);
  }


  //  Exchange interpolated fields
  if (code_id == 1) {
    CWP_Field_issend(code_name[0], coupling_name, field_name);
  }
  else {
    CWP_Field_irecv(code_name[0], coupling_name, field_name);
  }



  if (code_id == 1) {
    CWP_Field_wait_issend(code_name[0], coupling_name, field_name);
  }
  else {
    CWP_Field_wait_irecv(code_name[0], coupling_name, field_name);
  }

  //  Check
  double max_err = 0.;
  int    n_err   = 0;
  if (code_id == 2) {

    if (1) {
      log_trace("recv_val / coord = \n");
      for (int i = 0 ; i < n_located; i++) {
        int ivtx = located[i] - 1;
        log_trace("%d ("PDM_FMT_G_NUM"): %f %f %f / %f %f %f\n",
                  located[i],
                  pvtx_ln_to_gn[0][located[i]-1],
                  recv_val[3*i], recv_val[3*i+1], recv_val[3*i+2],
                  pvtx_coord[0][3*ivtx], pvtx_coord[0][3*ivtx+1], pvtx_coord[0][3*ivtx+2]);
      }
    }


    for (int i = 0 ; i < n_located ; i++) {
      int wrong = 0;

      for (int j = 0; j < 3; j++) {
        double err = ABS (recv_val[3*i + j] - pvtx_coord[0][3 * (located[i] -1) + j]);
        if (err > 1.e-4) {
          wrong = 1;
          printf("[%d] !! vtx "PDM_FMT_G_NUM" %d err = %g (coord#%d = %f, recv = %f)\n",
                 rank, pvtx_ln_to_gn[0][(located[i] - 1)], located[i], err, j, pvtx_coord[0][3*(located[i]-1) + j], recv_val[3*i + j]);
        }
        if (err > max_err) {
          max_err = err;
        }
      }

      if (wrong) {
        n_err++;
      }
    }

  }

  double global_max_err = 0.;
  MPI_Reduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    printf("Max error = %g\n", global_max_err);
  }

  int global_n_err = 0.;
  MPI_Reduce(&n_err, &global_n_err, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    printf("N error = %d\n", global_n_err);
  }



  //  Delete interface mesh
  CWP_Mesh_interf_del(code_name[0], coupling_name);

  //  Delete coupling
  CWP_Cpl_del(code_name[0], coupling_name);

  // Free memory
  free(code_name);
  free(coupled_code_name);
  free(is_active_rank);
  // free(time_init);
  free(intra_comm);

  if (current_rank_has_mesh) {
    for (int ipart = 0; ipart < n_part; ipart++) {
      free(pcell_face_idx[ipart]);
      free(pcell_face[ipart]);
      free(pface_vtx_idx[ipart]);
      free(pface_vtx[ipart]);
      free(pvtx_coord[ipart]);
      free(pcell_ln_to_gn[ipart]);
      free(pface_ln_to_gn[ipart]);
      free(pvtx_ln_to_gn[ipart]);
    }
  }

  free(pn_vtx);
  free(pn_cell);
  free(pn_face);
  free(pvtx_coord);
  free(pvtx_ln_to_gn);
  free(pcell_face_idx);
  free(pcell_face);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pcell_ln_to_gn);
  free(pface_ln_to_gn);

  if (code_id == 1) {
    free(send_val);
  }
  else {
    free(recv_val);
  }

  //  Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
