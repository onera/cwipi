#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_para_graph_dual.h"
#include "pdm_multipart.h"
#include "pdm_part.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_dmesh.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
#include "pdm_geom_elem.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int           argc,
           char        **argv,
           PDM_g_num_t  *nx,
           PDM_g_num_t  *ny,
           PDM_g_num_t  *nz,
           double       *lengthx,
           double       *lengthy,
           double       *lengthz,
           int          *n_part,
           int          *randomize,
           int          *random_seed,
           int          *post,
           int          *method,
           int          *use_multipart)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
        *ny = (PDM_g_num_t) _n;
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ny") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *ny = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthx = atof(argv[i]);
        *lengthy = atof(argv[i]);
        *lengthz = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-lx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthx = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ly") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthy = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-lz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthz = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      *randomize = 1;
    }
    else if (strcmp(argv[i], "-seed") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *random_seed = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else if (strcmp(argv[i], "-multipart") == 0) {
      *use_multipart = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  /*
   *  Set default values
   */
  PDM_g_num_t nx = 10;
  PDM_g_num_t ny = 10;
  PDM_g_num_t nz = 10;

  double lengthx = 1.;
  double lengthy = 1.;
  double lengthz = 1.;

  int n_part      = 1;
  int post        = 0;
  int randomize   = 0;
  int random_seed = 0;

#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  int use_multipart = 0;

  /*
   *  Read args
   */
  _read_args (argc,
              argv,
              &nx,
              &ny,
              &nz,
              &lengthx,
              &lengthy,
              &lengthz,
              &n_part,
              &randomize,
              &random_seed,
              &post,
              (int *) &method,
              &use_multipart);


  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  /*
   *  Create distributed mesh
   */
  double xmin = 0.;
  double ymin = 0.;
  double zmin = 0.;

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

  PDM_poly_vol_gen (comm,
                    xmin,
                    ymin,
                    zmin,
                    lengthx,
                    lengthy,
                    lengthz,
                    nx,
                    ny,
                    nz,
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

  if (i_rank == 0) printf("ng_cell = "PDM_FMT_G_NUM", ng_face = "PDM_FMT_G_NUM", ng_vtx = "PDM_FMT_G_NUM"\n", ng_cell, ng_face, ng_vtx);

  /*
   *  Create mesh partitions
   */
  PDM_multipart_t *mpart = NULL;
  PDM_part_t      *ppart = NULL;
  int n_join = 0;
  PDM_dmesh_t *dmesh;

  if (use_multipart) {
    /* Initialize multipart */
    mpart = PDM_multipart_create(1, &n_part, PDM_FALSE,
				 method, PDM_PART_SIZE_HOMOGENEOUS,
				 NULL, comm, PDM_OWNERSHIP_KEEP);

    /* Generate dmesh */
    dmesh = PDM_dmesh_create (PDM_OWNERSHIP_KEEP,
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

    PDM_dmesh_set (dmesh,
                   dvtx_coord,
                   dface_vtx_idx,
                   dface_vtx,
                   dface_cell,
                   dface_group_idx,
                   dface_group,
                   djoins_ids,
                   dface_join_idx,
                   dface_join);

    PDM_multipart_register_block (mpart, 0, dmesh);

    /* Connection between zones */
    int n_total_joins = 0;
    int *join_to_opposite = malloc(sizeof(int) * n_total_joins);
    PDM_multipart_register_joins (mpart, n_total_joins, join_to_opposite);

    /* Run */
    PDM_multipart_run_ppart (mpart);

    free (djoins_ids);
    free (dface_join_idx);
    free (dface_join);
    free (join_to_opposite);
  }

  else {
    int have_dcell_part = 0;

    int *dcell_part = (int *) malloc(dn_cell*sizeof(int));
    int *renum_properties_cell = NULL;
    int *renum_properties_face = NULL;
    int n_property_cell = 0;
    int n_property_face = 0;

    // id_part = 0;

    ppart = PDM_part_create(PDM_MPI_COMM_WORLD,
                            (PDM_part_split_t) method,
                            "PDM_PART_RENUM_CELL_NONE",
                            "PDM_PART_RENUM_FACE_NONE",
                            n_property_cell,
                            renum_properties_cell,
                            n_property_face,
                            renum_properties_face,
                            n_part,
                            dn_cell,
                            dn_face,
                            dn_vtx,
                            n_face_group,
                            NULL,
                            NULL,
                            NULL,
                            NULL,
                            have_dcell_part,
                            dcell_part,
                            dface_cell,
                            dface_vtx_idx,
                            dface_vtx,
                            NULL,
                            dvtx_coord,
                            NULL,
                            dface_group_idx,
                            dface_group);
    free(dcell_part);
  }


  if (post) {
    /* Prepare writer */
    PDM_writer_t *id_cs = PDM_writer_create ("Ensight",
                                   PDM_WRITER_FMT_ASCII,
                                   PDM_WRITER_TOPO_CONSTANTE,
                                   PDM_WRITER_OFF,
                                   "test_polyvol",
                                   "polyvol",
                                   PDM_MPI_COMM_WORLD,
                                   PDM_IO_ACCES_MPI_SIMPLE,
                                   1.,
                                   NULL);

    int id_geom = PDM_writer_geom_create (id_cs,
                                          "mesh",
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_OFF,
                                          n_part);

    // Cell local id
    int id_var_cell_g_num = PDM_writer_var_create (id_cs,
                                                   PDM_WRITER_OFF,
                                                   PDM_WRITER_VAR_SCALAIRE,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "cell_g_num");

    int id_var_num_part = PDM_writer_var_create (id_cs,
                                                 PDM_WRITER_OFF,
                                                 PDM_WRITER_VAR_SCALAIRE,
                                                 PDM_WRITER_VAR_ELEMENTS,
                                                 "num_part");

    int id_var_vtx_g_num = PDM_writer_var_create (id_cs,
                                                  PDM_WRITER_OFF,
                                                  PDM_WRITER_VAR_SCALAIRE,
                                                  PDM_WRITER_VAR_SOMMETS,
                                                  "vtx_g_num");

    int id_var_coo_x = PDM_writer_var_create (id_cs,
                                              PDM_WRITER_ON,
                                              PDM_WRITER_VAR_SCALAIRE,
                                              PDM_WRITER_VAR_SOMMETS,
                                              "coo_x");

    int id_var_coo_xyz = PDM_writer_var_create (id_cs,
                                                PDM_WRITER_ON,
                                                PDM_WRITER_VAR_VECTEUR,
                                                PDM_WRITER_VAR_SOMMETS,
                                                "coo_xyz");

    PDM_writer_step_beg (id_cs, 0.);

    /* Write geometry */
    int **face_vtx_n  = malloc (sizeof(int *) * n_part);
    int **cell_face_n = malloc (sizeof(int *) * n_part);

    PDM_real_t **val_cell_g_num = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_num_part   = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_vtx_g_num  = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_coo_x      = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_coo_xyz    = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {
      int n_cell;
      int n_face;
      int n_face_part_bound;
      int n_vtx;
      int n_proc;
      int n_t_part;
      int s_cell_face;
      int s_face_vtx;
      int s_face_group;
      int n_edge_group2;

      int n_bounds, n_joins, n_part_joins;
      int scell_face, sface_vtx, sface_bound, sface_join;
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
      PDM_g_num_t *face_bound_ln_to_gn, *face_join_ln_to_gn;
      int         *face_bound_idx, *face_bound, *face_join_idx, *face_join;
      int         **elt_vtx_idx;
      int         **elt_vtx;
      PDM_g_num_t **elt_section_ln_to_gn;

      if (use_multipart) {
        PDM_multipart_part_dim_get (mpart,
                                    0,
                                    i_part,
                                    &n_section,
                                    &n_elt,
                                    &n_cell,
                                    &n_face,
                                    &n_part_joins,
                                    &n_vtx,
                                    &n_proc,
                                    &n_t_part,
                                    &scell_face,
                                    &sface_vtx,
                                    &sface_bound,
                                    &n_bounds,
                                    &sface_join,
                                    &n_joins);

        PDM_multipart_part_val_get (mpart,
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
                                    &face_bound_idx,
                                    &face_bound,
                                    &face_bound_ln_to_gn,
                                    &face_join_idx,
                                    &face_join,
                                    &face_join_ln_to_gn);
      }

      else {
        PDM_part_part_dim_get (ppart,
                               i_part,
                               &n_cell,
                               &n_face,
                               &n_face_part_bound,
                               &n_vtx,
                               &n_proc,
                               &n_t_part,
                               &s_cell_face,
                               &s_face_vtx,
                               &s_face_group,
                               &n_edge_group2);

        PDM_part_part_val_get (ppart,
                               i_part,
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
                               &face_group_ln_to_gn);
      }

      PDM_writer_geom_coord_set (id_cs,
                                 id_geom,
                                 i_part,
                                 n_vtx,
                                 vtx,
                                 vtx_ln_to_gn);

      face_vtx_n[i_part] = malloc (sizeof(int) * n_face);
      for (int i = 0; i < n_face; i++) {
        face_vtx_n[i_part][i] = face_vtx_idx[i+1] - face_vtx_idx[i];
      }

      cell_face_n[i_part] = malloc (sizeof(int) * n_cell);
      for (int i = 0; i < n_cell; i++) {
        cell_face_n[i_part][i] = cell_face_idx[i+1] - cell_face_idx[i];
      }

      PDM_writer_geom_cell3d_cellface_add (id_cs,
                                           id_geom,
                                           i_part,
                                           n_cell,
                                           n_face,
                                           face_vtx_idx,
                                           face_vtx_n[i_part],
                                           face_vtx,
                                           cell_face_idx,
                                           cell_face_n[i_part],
                                           cell_face,
                                           cell_ln_to_gn);

      val_cell_g_num[i_part] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);
      val_num_part[i_part]   = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);
      for (int i = 0; i < n_cell; i++) {
        val_cell_g_num[i_part][i] = (PDM_real_t) cell_ln_to_gn[i];
        val_num_part[i_part][i] = (PDM_real_t) (i_rank*n_part + i_part);
      }

      val_vtx_g_num[i_part] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_vtx);
      val_coo_x[i_part]     = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_vtx);
      val_coo_xyz[i_part]   = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_vtx * 3);
      for (int i = 0; i < n_vtx; i++) {
        val_vtx_g_num[i_part][i]   = (PDM_real_t) vtx_ln_to_gn[i];
        val_coo_x[i_part][i]       = vtx[3*i];
        val_coo_xyz[i_part][3*i  ] = vtx[3*i  ];
        val_coo_xyz[i_part][3*i+1] = vtx[3*i+1];
        val_coo_xyz[i_part][3*i+2] = vtx[3*i+2];
      }
    }

    PDM_writer_geom_write (id_cs,
                           id_geom);

    // write variables
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_writer_var_set (id_cs,
                          id_var_cell_g_num,
                          id_geom,
                          i_part,
                          val_cell_g_num[i_part]);

      PDM_writer_var_set (id_cs,
                          id_var_num_part,
                          id_geom,
                          i_part,
                          val_num_part[i_part]);

      PDM_writer_var_set (id_cs,
                          id_var_vtx_g_num,
                          id_geom,
                          i_part,
                          val_vtx_g_num[i_part]);
    }

    PDM_writer_var_write (id_cs,
                          id_var_cell_g_num);
    PDM_writer_var_write (id_cs,
                          id_var_num_part);
    PDM_writer_var_write (id_cs,
                          id_var_vtx_g_num);

    // PDM_writer_var_free (id_cs,
    //                      id_var_cell_g_num);
    // PDM_writer_var_free (id_cs,
    //                      id_var_num_part);
    // PDM_writer_var_free (id_cs,
    //                      id_var_vtx_g_num);

    for (int nstep = 0; nstep < 10; nstep++) {

      double tstep = nstep * 0.01;

      if (nstep > 0) {
        PDM_writer_step_beg(id_cs, tstep);
      }

      for (int i_part = 0; i_part < n_part; i_part++) {
        PDM_writer_var_set (id_cs,
                            id_var_coo_x,
                            id_geom,
                            i_part,
                            val_coo_x[i_part]);
        PDM_writer_var_set (id_cs,
                            id_var_coo_xyz,
                            id_geom,
                            i_part,
                            val_coo_xyz[i_part]);
      }

      PDM_writer_var_write (id_cs,
                            id_var_coo_x);
      PDM_writer_var_write (id_cs,
                            id_var_coo_xyz);

      PDM_writer_var_data_free (id_cs,
                                id_var_coo_x);
      PDM_writer_var_data_free (id_cs,
                                id_var_coo_xyz);

      PDM_writer_step_end (id_cs);
    }


    for (int i_part = 0; i_part < n_part; i_part++) {
      free (val_cell_g_num[i_part]);
      free (val_num_part[i_part]);
      free (val_vtx_g_num[i_part]);
      free (val_coo_x[i_part]);
      free (val_coo_xyz[i_part]);
      free (cell_face_n[i_part]);
      free (face_vtx_n[i_part]);
    }
    free (val_cell_g_num);
    free (val_num_part);
    free (val_vtx_g_num);
    free (val_coo_x);
    free (val_coo_xyz);
    free (cell_face_n);
    free (face_vtx_n);

    // PDM_writer_var_free (id_cs,
    //                      id_var_coo_x);
    // PDM_writer_var_free (id_cs,
    //                      id_var_coo_xyz);

    // PDM_writer_geom_data_free (id_cs, id_geom);
    // PDM_writer_geom_free (id_cs, id_geom);
    PDM_writer_free (id_cs);
  }

  /*
   *  Finalize
   */
  free (dcell_face_idx);
  free (dcell_face);
  free (dface_cell);
  free (dface_vtx_idx);
  free (dface_vtx);
  free (dvtx_coord);
  free (dface_group_idx);
  free (dface_group);

  if (use_multipart) {
    PDM_multipart_free (mpart);
    PDM_dmesh_free (dmesh);
  }
  else {
    PDM_part_free (ppart);
  }

  PDM_MPI_Finalize();

  return 0;
}
