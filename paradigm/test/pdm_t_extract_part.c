#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dcube_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_gnum.h"
#include "pdm_part_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_extract_part.h"
#include "pdm_vtk.h"
#include "pdm_dmesh.h"
#include "pdm_unique.h"
#include "pdm_part_geom.h"
#include "pdm_logging.h"
#include "pdm_priv.h"

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
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
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
 * \param [inout]   part_method Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
           int           *post,
           int           *part_method)
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
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = 1;
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

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
     (int *) &part_method);

  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell = NULL;
  int          *dface_vtx_idx = NULL;
  PDM_g_num_t  *dface_vtx = NULL;
  double       *dvtx_coord = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  PDM_dcube_t* dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          0.,
                                          0.,
                                          0.,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                        &n_face_group,
                        &dn_cell,
                        &dn_face,
                        &dn_vtx,
                        &dface_vtxL,
                        &dFaceGroupL);

  PDM_dcube_gen_data_get(dcube,
                         &dface_cell,
                         &dface_vtx_idx,
                         &dface_vtx,
                         &dvtx_coord,
                         &dface_group_idx,
                         &dface_group);

  /*
   * Create dmesh
   */
  int n_jn = 0;
  PDM_dmesh_t* dm = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                     dn_cell,
                                     dn_face,
                                     -1, // dn_edge
                                     dn_vtx,
                                     n_face_group, // n_bnd
                                     n_jn,
                                     comm);

  int *dface_join_idx = (int *) malloc( (n_jn+1) * sizeof(int));
  dface_join_idx[0] = 0;
  PDM_dmesh_set(dm,
                dvtx_coord,
                dface_vtx_idx,
                dface_vtx,
                dface_cell,
                dface_group_idx,
                dface_group,
                NULL,
                dface_join_idx,
                NULL);

  /*
   * Partitionnement
   */
  int n_zone = 1;
  int n_part_zones = n_part;
  PDM_multipart_t *mpart_id = PDM_multipart_create(n_zone,
						   &n_part_zones,
						   PDM_FALSE,
						   part_method,
						   PDM_PART_SIZE_HOMOGENEOUS,
						   NULL,
						   comm,
						   PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart_id, -1, "PDM_PART_RENUM_CELL_NONE", NULL, "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_register_block(mpart_id, 0, dm);
  PDM_multipart_run_ppart(mpart_id);

  /*
   * Get the partition zone
   */
  int i_zone = 0;

  double      **cell_center             = (double      **) malloc( n_part_zones * sizeof(double      *));
  int         **selected_l_num          = (int         **) malloc( n_part_zones * sizeof(int         *));
  PDM_g_num_t **pcell_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pface_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pvtx_ln_to_gn           = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  int          *pn_cell                 = (int          *) malloc( n_part_zones * sizeof(int          ));
  int          *pn_face                 = (int          *) malloc( n_part_zones * sizeof(int          ));
  int          *pn_vtx                  = (int          *) malloc( n_part_zones * sizeof(int          ));
  int          *pn_select_cell          = (int          *) malloc( n_part_zones * sizeof(int          ));
  // double      **weight                  = (double      **) malloc( n_part_zones * sizeof(double      *));
  int         **pcell_face              = (int         **) malloc( n_part_zones * sizeof(int         *));
  int         **pcell_face_idx          = (int         **) malloc( n_part_zones * sizeof(int         *));
  int         **pface_vtx               = (int         **) malloc( n_part_zones * sizeof(int         *));
  int         **pface_vtx_idx           = (int         **) malloc( n_part_zones * sizeof(int         *));
  double      **pvtx_coord              = (double      **) malloc( n_part_zones * sizeof(double      *));
  // double      **tmp_extract_cell_center = (double      **) malloc( n_part_zones * sizeof(double      *));

  for (int i_part = 0; i_part < n_part_zones; i_part++){

    int n_proc, tn_part;
    int n_cell, n_face, n_vtx, n_bounds, n_joins, n_part_joins;
    int scell_face, sface_vtx, sface_bound, sface_join;
    int  n_section;
    int* n_elt;

    PDM_multipart_part_dim_get(mpart_id, i_zone, i_part, &n_section, &n_elt,
                               &n_cell, &n_face, &n_part_joins, &n_vtx, &n_proc, &tn_part,
                               &scell_face, &sface_vtx, &sface_bound, &n_bounds, &sface_join, &n_joins);

    double       *vtx;
    int          *cell_face_idx, *cell_face, *face_cell, *face_vtx_idx, *face_vtx;
    int          *face_bound_idx, *face_bound, *face_join_idx, *face_join;
    int          *face_part_bound_proc_idx, *face_part_bound_part_idx, *face_part_bound;
    PDM_g_num_t  *cell_ln_to_gn, *face_ln_to_gn, *vtx_ln_to_gn, *face_bound_ln_to_gn, *face_join_ln_to_gn;
    int          *cell_tag, *face_tag, *vtx_tag;
    int         **elt_vtx_idx;
    int         **elt_vtx;
    PDM_g_num_t **elt_section_ln_to_gn;


    PDM_multipart_part_val_get(mpart_id, i_zone, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                               &cell_tag, &cell_face_idx, &cell_face, &cell_ln_to_gn,
                               &face_tag, &face_cell, &face_vtx_idx, &face_vtx, &face_ln_to_gn,
                               &face_part_bound_proc_idx, &face_part_bound_part_idx, &face_part_bound,
                               &vtx_tag, &vtx, &vtx_ln_to_gn, &face_bound_idx, &face_bound,
                               &face_bound_ln_to_gn, &face_join_idx, &face_join, &face_join_ln_to_gn);

    pn_cell       [i_part] = n_cell;
    pcell_ln_to_gn[i_part] = cell_ln_to_gn;
    pface_ln_to_gn[i_part] = face_ln_to_gn;
    pvtx_ln_to_gn [i_part] = vtx_ln_to_gn;
    pcell_face    [i_part] = cell_face;
    pcell_face_idx[i_part] = cell_face_idx;
    pn_face       [i_part] = n_face;
    pn_vtx        [i_part] = n_vtx;

    pface_vtx    [i_part] = face_vtx;
    pface_vtx_idx[i_part] = face_vtx_idx;
    pvtx_coord   [i_part] = vtx;

    // char filename[999];
    // sprintf(filename, "mesh_before_extract_%3.3d_%3.3d.vtk", i_part, i_rank);
    // PDM_vtk_write_polydata(filename,
    //                        pn_vtx[i_part],
    //                        pvtx_coord[i_part],
    //                        pvtx_ln_to_gn[i_part],
    //                        pn_face[i_part],
    //                        pface_vtx_idx[i_part],
    //                        pface_vtx[i_part],
    //                        pface_ln_to_gn[i_part],
    //                        NULL);

    /*
     * Compute center-cell and extract cells corresponding to criteria
     */
    double *face_center         = (double *) malloc( 3 * n_face * sizeof(double));

    for(int i_face = 0; i_face < n_face; ++i_face) {
      face_center[3*i_face  ] = 0.;
      face_center[3*i_face+1] = 0.;
      face_center[3*i_face+2] = 0.;
      int n_vtx_on_face = face_vtx_idx[i_face+1] - face_vtx_idx[i_face];
      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = PDM_ABS(face_vtx[idx_vtx])-1;
        face_center[3*i_face  ] += vtx[3*i_vtx  ];
        face_center[3*i_face+1] += vtx[3*i_vtx+1];
        face_center[3*i_face+2] += vtx[3*i_vtx+2];
      }
      face_center[3*i_face  ] = face_center[3*i_face  ] / n_vtx_on_face;
      face_center[3*i_face+1] = face_center[3*i_face+1] / n_vtx_on_face;
      face_center[3*i_face+2] = face_center[3*i_face+2] / n_vtx_on_face;
    }

    cell_center[i_part] = (double *) malloc( 3 * n_cell * sizeof(double));

    for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

      cell_center[i_part][3*i_cell  ] = 0.;
      cell_center[i_part][3*i_cell+1] = 0.;
      cell_center[i_part][3*i_cell+2] = 0.;

      int n_face_on_cell = cell_face_idx[i_cell+1] - cell_face_idx[i_cell];

      for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {

        int i_face = PDM_ABS(cell_face[idx_face])-1;
        cell_center[i_part][3*i_cell  ] += face_center[3*i_face  ];
        cell_center[i_part][3*i_cell+1] += face_center[3*i_face+1];
        cell_center[i_part][3*i_cell+2] += face_center[3*i_face+2];
      }
      cell_center[i_part][3*i_cell  ] = cell_center[i_part][3*i_cell  ] / n_face_on_cell;
      cell_center[i_part][3*i_cell+1] = cell_center[i_part][3*i_cell+1] / n_face_on_cell;
      cell_center[i_part][3*i_cell+2] = cell_center[i_part][3*i_cell+2] / n_face_on_cell;
    }

    free(face_center);

    selected_l_num[i_part] = (int         *) malloc(  n_cell          * sizeof(int        ));

    /*
     * Sub-part
     */
    double bbox[6];
    // bbox[0] = 0.3;
    // bbox[1] = 0.3;
    // bbox[2] = 0.35;
    // bbox[3] = 0.7;
    // bbox[4] = 0.7;
    // bbox[5] = 0.65;
    bbox[0] = 0.5;
    bbox[1] = 0.65;
    bbox[2] = 0.75;
    bbox[3] = 1.25;
    bbox[4] = 1.25;
    bbox[5] = 1.25;

    int n_select_cell = 0;
    for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

      int inside = 1;
      for(int i = 0; i < 3; ++i) {
        if (cell_center[i_part][3*i_cell+i] > bbox[i+3] || cell_center[i_part][3*i_cell+i] < bbox[i]) {
          inside = 0;
        }
      }
      if(inside == 1) {
        selected_l_num[i_part][n_select_cell]     = i_cell;
        n_select_cell++;

      }
    }

    selected_l_num[i_part] = realloc(selected_l_num[i_part], n_select_cell * sizeof(int        ));
    pn_select_cell[i_part] = n_select_cell;

  }
  free(dface_join_idx);

  /*
   * Extract
   */
  int n_part_out = 1;
  PDM_bool_t equilibrate = PDM_FALSE;
  // PDM_bool_t equilibrate = PDM_TRUE;
  PDM_extract_part_t* extrp = PDM_extract_part_create(3,
                                                      n_part,
                                                      n_part_out,
                                                      equilibrate,
                                                      PDM_SPLIT_DUAL_WITH_PTSCOTCH,
                                                      PDM_OWNERSHIP_KEEP,
                                                      comm);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_extract_part_part_set(extrp,
                              i_part,
                              pn_cell[i_part],
                              pn_face[i_part],
                              -1, // pn_edge[i_part],
                              pn_vtx[i_part],
                              pcell_face_idx[i_part],
                              pcell_face[i_part],
                              NULL, //pface_edge_idx[i_part],
                              NULL, //pface_edge[i_part],
                              NULL, //pedge_vtx[i_part],
                              pface_vtx_idx[i_part],
                              pface_vtx[i_part],
                              pcell_ln_to_gn[i_part],
                              pface_ln_to_gn[i_part],
                              NULL, //pedge_ln_to_gn[i_part],
                              pvtx_ln_to_gn[i_part],
                              pvtx_coord[i_part]);

    PDM_extract_part_selected_lnum_set(extrp,
                                       i_part,
                                       pn_select_cell[i_part],
                                       selected_l_num[i_part]);

    // PDM_log_trace_array_int(selected_l_num[i_part], pn_select_cell[i_part], "selected_l_num ::");

  }


  PDM_extract_part_compute(extrp);

  int          *pn_extract_face        = malloc(n_part_out * sizeof(int          ));
  int          *pn_extract_vtx         = malloc(n_part_out * sizeof(int          ));
  int         **pextract_face_vtx      = malloc(n_part_out * sizeof(int         *));
  int         **pextract_face_vtx_idx  = malloc(n_part_out * sizeof(int         *));
  double      **pextract_vtx           = malloc(n_part_out * sizeof(double      *));
  PDM_g_num_t **pextract_face_ln_to_gn = malloc(n_part_out * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pextract_vtx_ln_to_gn  = malloc(n_part_out * sizeof(PDM_g_num_t *));


  for(int i_part = 0; i_part < n_part_out; ++i_part) {

    pn_extract_face[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_FACE);

    pn_extract_vtx[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                           i_part,
                                                           PDM_MESH_ENTITY_VERTEX);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                      &pextract_face_vtx[i_part],
                                      &pextract_face_vtx_idx[i_part],
                                      PDM_OWNERSHIP_KEEP);

    PDM_extract_part_vtx_coord_get(extrp,
                                   i_part,
                                   &pextract_vtx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &pextract_face_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_VERTEX,
                                  &pextract_vtx_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);
  }

  /*
   * Export vtk en légende
   */
  if(0 == 1) {
    for(int i_part = 0; i_part < n_part_out; ++i_part) {

      char filename[999];
      sprintf(filename, "extract_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                pn_extract_vtx[i_part],
                                pextract_vtx[i_part],
                                NULL, NULL);

      PDM_log_trace_connectivity_int(pextract_face_vtx_idx[i_part],
                                     pextract_face_vtx    [i_part],
                                     pn_extract_face[i_part], " pextract_face_vtx :: ");

      sprintf(filename, "extract_face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_polydata(filename,
                             pn_extract_vtx[i_part],
                             pextract_vtx[i_part],
                             pextract_vtx_ln_to_gn[i_part],
                             pn_extract_face[i_part],
                             pextract_face_vtx_idx[i_part],
                             pextract_face_vtx[i_part],
                             pextract_face_ln_to_gn[i_part],
                             NULL);
    }
  }

  free(pn_extract_face);
  free(pn_extract_vtx);
  free(pextract_face_vtx     );
  free(pextract_face_vtx_idx );
  free(pextract_vtx          );
  free(pextract_face_ln_to_gn);
  free(pextract_vtx_ln_to_gn );

  PDM_extract_part_free(extrp);


  for (int i_part = 0; i_part < n_part_zones; i_part++){
    free(cell_center       [i_part]);
    free(selected_l_num    [i_part]);
  }
  free(cell_center);
  free(selected_l_num);
  free(pn_cell);
  free(pn_face);
  free(pn_vtx);
  free(pn_select_cell);

  free(pcell_ln_to_gn  );
  free(pface_ln_to_gn  );
  free(pvtx_ln_to_gn  );
  free(pcell_face    );
  free(pcell_face_idx);
  free(pface_vtx     );
  free(pface_vtx_idx );
  free(pvtx_coord    );

  PDM_multipart_free(mpart_id);
  PDM_dcube_gen_free(dcube);
  PDM_dmesh_free(dm);


  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize();

  return 0;
}
