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

  PDM_g_num_t* distrib_cell = PDM_compute_entity_distribution(comm, dn_cell);

  /*
   * Get the partition zone
   */
  int i_zone = 0;

  PDM_g_num_t **pcell_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pface_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pvtx_ln_to_gn           = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  int          *pn_cell                 = (int          *) malloc( n_part_zones * sizeof(int          ));
  int          *pn_face                 = (int          *) malloc( n_part_zones * sizeof(int          ));
  int          *pn_vtx                  = (int          *) malloc( n_part_zones * sizeof(int          ));

  int         **pcell_face              = (int         **) malloc( n_part_zones * sizeof(int         *));
  int         **pcell_face_idx          = (int         **) malloc( n_part_zones * sizeof(int         *));
  int         **pface_vtx               = (int         **) malloc( n_part_zones * sizeof(int         *));
  int         **pface_vtx_idx           = (int         **) malloc( n_part_zones * sizeof(int         *));
  double      **pvtx_coord              = (double      **) malloc( n_part_zones * sizeof(double      *));

  PDM_g_num_t **target_g_num   = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  int          *pn_target_cell = (int          *) malloc( n_part_zones * sizeof(int          ));

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

    /* On prends tout les gnum impaire dans un ordre reverse */
    int n_target_cell = 0;
    target_g_num[i_part] = malloc(dn_cell * sizeof(PDM_g_num_t));

    for(int i = 0; i < dn_cell; ++i ) {
      PDM_g_num_t gnum = distrib_cell[i_rank+1] - i;
      if(gnum % 2 == 0) {
        target_g_num[i_part][n_target_cell++] = gnum;
      }
    }

    target_g_num  [i_part] = realloc(target_g_num[i_part], n_target_cell * sizeof(PDM_g_num_t));
    pn_target_cell[i_part] = n_target_cell;

    if(0 == 1) {
      PDM_log_trace_array_long(target_g_num  [i_part], pn_target_cell[i_part], "target_g_num :: ");
    }


  }
  free(dface_join_idx);
  free(distrib_cell);

  /*
   * Extract
   */
  int n_part_out = 1;
  // PDM_bool_t equilibrate = PDM_FALSE;
  PDM_bool_t equilibrate = PDM_TRUE;
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

    PDM_extract_part_target_gnum_set(extrp,
                                     i_part,
                                     pn_target_cell[i_part],
                                     target_g_num  [i_part]);

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

    // PDM_g_num_t* pextract_parent_cell_ln_to_gn = NULL;
    // int n_cell = PDM_extract_part_parent_ln_to_gn_get(extrp,
    //                                                   i_part,
    //                                                   PDM_MESH_ENTITY_CELL,
    //                                                   &pextract_parent_cell_ln_to_gn,
    //                                                   PDM_OWNERSHIP_KEEP);

    // for(int i = 0; i < n_cell; ++i) {
    //   assert(pextract_parent_cell_ln_to_gn[i] == target_g_num[i_part][i]);
    // }
    // PDM_log_trace_array_long(pextract_parent_cell_ln_to_gn, n_cell, "pextract_parent_cell_ln_to_gn ::");

  }

  /*
   * Export vtk en légende
   */
  if(post) {
    for(int i_part = 0; i_part < n_part_out; ++i_part) {

      char filename[999];
      sprintf(filename, "extract_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                pn_extract_vtx[i_part],
                                pextract_vtx[i_part],
                                NULL, NULL);

      // PDM_log_trace_connectivity_int(pextract_face_vtx_idx[i_part],
      //                                pextract_face_vtx    [i_part],
      //                                pn_extract_face[i_part], " pextract_face_vtx :: ");

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
    free(target_g_num    [i_part]);
  }
  free(pn_target_cell);
  free(target_g_num);
  free(pn_cell);
  free(pn_face);
  free(pn_vtx);

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

  PDM_MPI_Finalize();

  printf("-- Fin test\n");
  fflush(stdout);

  return 0;
}
