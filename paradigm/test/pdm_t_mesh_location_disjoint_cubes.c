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
#include "pdm_part.h"
#include "pdm_multipart.h"
#include "pdm_dcube_gen.h"
#include "pdm_mesh_location.h"
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
#include "pdm_dmesh_nodal_to_dmesh.h"
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
     "  -t      <level>  Bounding boxes tolerance.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -p      <level>  Number of points to locate.\n\n"
     "  -octree          Use octree-based method.\n\n"
     "  -dbbree          Use dbbtree-based method.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scotch       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}



/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   tolerance  Bounding boxes tolerance
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   post       Ensight outputs status
 * \param [inout]   method     Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                          argc,
           char                       **argv,
           PDM_g_num_t                 *n_vtx_seg,
           double                      *length,
           double                      *separation_x,
           double                      *separation_y,
           double                      *separation_z,
           int                         *deform,
           double                      *tolerance,
           double                      *marge,
           int                         *n_part,
           PDM_g_num_t                 *n_pts,
           int                         *post,
           int                         *part_method,
           PDM_mesh_location_method_t  *loc_method,
           int                         *disable_uvw,
           int                         *use_tgt_nodes,
           int                         *extension_depth_tgt,
           int                         *extension_depth_src)
{
  int i = 1;

  PDM_UNUSED (post);

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
    else if (strcmp(argv[i], "-def") == 0) {
      *deform = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-m") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *marge = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_pts = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_PART_SPLIT_HILBERT;
    }
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = PDM_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DBBTREE;
    }
    else if (strcmp(argv[i], "-no_uvw") == 0) {
      *disable_uvw = 1;
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-nodes") == 0) {
      *use_tgt_nodes = 1;
    }
    else if (strcmp(argv[i], "-ext_depth") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_tgt = atoi(argv[i]);
        *extension_depth_src = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ext_depth_tgt") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_tgt = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ext_depth_src") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_src = atoi(argv[i]);
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void _rotate (const int  n_pts,
                     double    *coord)
{
  double R[3][3] = {{0.9362934, -0.2896295, 0.1986693},
                    {0.3129918,  0.9447025, -0.0978434},
                    {-0.1593451,  0.1537920,  0.9751703}};

  for (int i = 0; i < n_pts; i++) {
    double x = coord[3*i];
    double y = coord[3*i+1];
    double z = coord[3*i+2];

    for (int j = 0; j < 3; j++) {
      coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
    }
  }
}

static void
_cube_mesh
(
 const PDM_MPI_Comm        comm,
 const int                 n_part,
 const PDM_split_dual_t    part_method,
 const PDM_g_num_t         n_vtx_seg,
 const double              xmin,
 const double              ymin,
 const double              zmin,
 const double              length,
 const int                 deform,
 const int                 part_extension_depth,
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
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         xmin,
                                                         ymin,
                                                         zmin,
                                                         PDM_MESH_NODAL_HEXA8,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];


  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  PDM_dmesh_t*dmesh2 = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh2);

  int         *dface_vtx_idx;
  PDM_g_num_t *dface_vtx;
  int dn_face = PDM_dmesh_connectivity_get(dmesh2, PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                           &dface_vtx,
                                           &dface_vtx_idx,
                                           PDM_OWNERSHIP_KEEP);

  int         *dface_cell_idx = NULL;
  PDM_g_num_t *dface_cell = NULL;
  dn_face = PDM_dmesh_connectivity_get(dmesh2, PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                       &dface_cell,
                                       &dface_cell_idx,
                                       PDM_OWNERSHIP_KEEP);
  assert(dface_cell_idx == NULL);

  int         *dcell_face_idx;
  PDM_g_num_t *dcell_face;
  int dn_cell = PDM_dmesh_connectivity_get(dmesh2, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                           &dcell_face,
                                           &dcell_face_idx,
                                           PDM_OWNERSHIP_KEEP);

  PDM_g_num_t *dface_group     = NULL;
  int         *dface_group_idx = NULL;
  int n_face_group = PDM_dmesh_bound_get(dmesh2,
                                         PDM_BOUND_TYPE_FACE,
                                         &dface_group,
                                         &dface_group_idx,
                                         PDM_OWNERSHIP_KEEP);




  if (deform) {
    /*for (int i = 0; i < dn_vtx; i++) {
      double x = dvtx_coord[3*i];
      double z = dvtx_coord[3*i + 2];

      dvtx_coord[3*i]     += 0.1 * z * z;
      dvtx_coord[3*i + 2] += 0.2 * cos(PDM_PI * x);
      }*/
    _rotate (dn_vtx,
             dvtx_coord);
  }

  /*
   *  Create mesh partitiions
   */

  /* Initialize multipart */
  PDM_multipart_t *mpart = PDM_multipart_create(1,
                                                &n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  int n_join = 0;

  /* Generate dmesh */
  PDM_dmesh_t *dmesh = PDM_dmesh_create (PDM_OWNERSHIP_KEEP,
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
  PDM_multipart_register_joins (mpart,
                                n_total_joins,
                                join_to_opposite);

  /* Run */
  PDM_multipart_run_ppart (mpart);

  free (djoins_ids);
  free (dface_join_idx);
  free (dface_join);
  free (join_to_opposite);


  // PDM_dcube_gen_free(dcube);


  PDM_part_extension_t *part_ext = NULL;

  if (part_extension_depth > 0) {
    part_ext = PDM_part_extension_create(1,
                                         &n_part,
                                         PDM_EXTEND_FROM_FACE,
                                         part_extension_depth,
                                         comm,
                                         PDM_OWNERSHIP_KEEP);
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

      PDM_multipart_part_dim_get (mpart,
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
                                  &face_group_idx,
                                  &face_group,
                                  &face_group_ln_to_gn,
                                  &face_join_idx,
                                  &face_join,
                                  &face_join_ln_to_gn);


      PDM_part_extension_set_part(part_ext,
                                  0,
                                  i_part,
                                  n_cell,
                                  n_face,
                                  n_face_part_bound,
                                  n_face_group,
                                  0,   // n_edge
                                  n_vtx,
                                  cell_face_idx,
                                  cell_face,
                                  face_cell,
                                  NULL, // face_edge_idx
                                  NULL, // face_edge
                                  face_vtx_idx,
                                  face_vtx,
                                  NULL, //edge_vtx
                                  face_group_idx,
                                  face_group,
                                  NULL, // face_join_idx
                                  NULL, // face_join
                                  face_part_bound_proc_idx,
                                  face_part_bound_part_idx,
                                  face_part_bound,
                                  NULL, // vtx_part_bound_proc_idx
                                  NULL, // vtx_part_bound_part_idx
                                  NULL, // vtx_part_bound
                                  cell_ln_to_gn,
                                  face_ln_to_gn,
                                  NULL, // edge_ln_to_gn
                                  vtx_ln_to_gn,
                                  face_group_ln_to_gn,
                                  vtx);
    }

    PDM_part_extension_compute(part_ext);
  }


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

    PDM_multipart_part_dim_get (mpart,
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
                                &face_group_idx,
                                &face_group,
                                &face_group_ln_to_gn,
                                &face_join_idx,
                                &face_join,
                                &face_join_ln_to_gn);

    int n_ext_vtx  = 0;
    int n_ext_face = 0;
    int n_ext_cell = 0;
    double      *ext_vtx_coord     = NULL;
    PDM_g_num_t *ext_vtx_ln_to_gn  = NULL;
    int         *ext_cell_face     = NULL;
    int         *ext_cell_face_idx = NULL;
    PDM_g_num_t *ext_cell_ln_to_gn = NULL;
    int         *ext_face_vtx      = NULL;
    int         *ext_face_vtx_idx  = NULL;
    PDM_g_num_t *ext_face_ln_to_gn = NULL;

    if (part_extension_depth > 0) {
      /* Vertices */
      n_ext_vtx = PDM_part_extension_coord_get(part_ext,
                                               0,
                                               i_part,
                                               &ext_vtx_coord);

      PDM_part_extension_ln_to_gn_get(part_ext,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_VERTEX,
                                      &ext_vtx_ln_to_gn);


      /* Cells */
      n_ext_cell = PDM_part_extension_connectivity_get(part_ext,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &ext_cell_face,
                                                       &ext_cell_face_idx);

      PDM_part_extension_ln_to_gn_get(part_ext,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &ext_cell_ln_to_gn);


      /* Faces */
      n_ext_face = PDM_part_extension_connectivity_get(part_ext,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                       &ext_face_vtx,
                                                       &ext_face_vtx_idx);
      PDM_part_extension_ln_to_gn_get(part_ext,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    &ext_face_ln_to_gn);
    }


    *(pn_cell)[i_part] = n_cell + n_ext_cell;
    *(pn_face)[i_part] = n_face + n_ext_face;
    *(pn_vtx)[i_part]  = n_vtx  + n_ext_vtx;

    /* Vertices */
    (*pvtx_coord)[i_part] = (double *) malloc(sizeof(double) * 3 * (n_vtx + n_ext_vtx));
    memcpy((*pvtx_coord)[i_part], vtx, sizeof(double) * 3 * n_vtx);

    (*pvtx_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_vtx + n_ext_vtx));
    memcpy((*pvtx_ln_to_gn)[i_part], vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_vtx);


    /* Cells */
    (*pcell_face_idx)[i_part] = (int *) malloc(sizeof(int) * (n_cell + n_ext_cell + 1));
    memcpy((*pcell_face_idx)[i_part], cell_face_idx, sizeof(int) * (n_cell + 1));

    s_cell_face = cell_face_idx[n_cell];
    if (part_extension_depth > 0) {
      s_cell_face += ext_cell_face_idx[n_ext_cell];
    }
    (*pcell_face)[i_part] = (int *) malloc(sizeof(int) * s_cell_face);
    memcpy((*pcell_face)[i_part], cell_face, sizeof(int) * cell_face_idx[n_cell]);

    (*pcell_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_cell + n_ext_cell));
    memcpy((*pcell_ln_to_gn)[i_part], cell_ln_to_gn, sizeof(PDM_g_num_t) * n_cell);


    /* Faces */
    (*pface_vtx_idx)[i_part] = (int *) malloc(sizeof(int) * (n_face + n_ext_face + 1));
    memcpy((*pface_vtx_idx)[i_part], face_vtx_idx, sizeof(int) * (n_face + 1));

    s_face_vtx = face_vtx_idx[n_face];
    if (part_extension_depth > 0) {
      s_face_vtx += ext_face_vtx_idx[n_ext_face];
    }
    (*pface_vtx)[i_part] = (int *) malloc(sizeof(int) * s_face_vtx);
    memcpy((*pface_vtx)[i_part], face_vtx, sizeof(int) * face_vtx_idx[n_face]);

    (*pface_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_face + n_ext_face));
    memcpy((*pface_ln_to_gn)[i_part], face_ln_to_gn, sizeof(PDM_g_num_t) * n_face);


    if (part_extension_depth > 0) {
      /* Vertices */
      memcpy((*pvtx_coord)[i_part] + 3*n_vtx, ext_vtx_coord, sizeof(double) * 3 * n_ext_vtx);
      memcpy((*pvtx_ln_to_gn)[i_part] + n_vtx, ext_vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_ext_vtx);

      /* Cells */
      for (int i = 1; i <= n_ext_cell; i++) {
        (*pcell_face_idx)[i_part][n_cell + i] = cell_face_idx[n_cell] + ext_cell_face_idx[i];
      }

      memcpy((*pcell_face)[i_part] + cell_face_idx[n_cell],
             ext_cell_face,
             sizeof(int) * ext_cell_face_idx[n_ext_cell]);

      memcpy((*pcell_ln_to_gn)[i_part] + n_cell, ext_cell_ln_to_gn, sizeof(PDM_g_num_t) * n_ext_cell);

      /* Faces */
      for (int i = 1; i <= n_ext_face; i++) {
        (*pface_vtx_idx)[i_part][n_face + i] = face_vtx_idx[n_face] + ext_face_vtx_idx[i];
      }

      memcpy((*pface_vtx)[i_part] + face_vtx_idx[n_face],
             ext_face_vtx,
             sizeof(int) * ext_face_vtx_idx[n_ext_face]);

      memcpy((*pface_ln_to_gn)[i_part] + n_face, ext_face_ln_to_gn, sizeof(PDM_g_num_t) * n_ext_face);
    }
  }

  PDM_multipart_free (mpart);
  PDM_dmesh_free (dmesh);
  PDM_part_extension_free (part_ext);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_dcube_nodal_gen_free(dcube);
}


static void _export_point_cloud
(
 char         *filename,
 int           n_part,
 int          *n_pts,
 double      **coord,
 PDM_g_num_t **g_num,
 PDM_g_num_t **location
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\npoints\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  int n_pts_t = 0;
  for (int ipart = 0; ipart < n_part; ipart++) {
    n_pts_t += n_pts[ipart];
  }

  fprintf(f, "POINTS %d double\n", n_pts_t);
  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int i = 0; i < n_pts[ipart]; i++) {
      for (int j = 0; j < 3; j++) {
        fprintf(f, "%f ", coord[ipart][3*i + j]);
      }
      fprintf(f, "\n");
    }
  }

  fprintf(f, "CELLS %d %d\n", n_pts_t, 2*n_pts_t);
  for (int i = 0; i < n_pts_t; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts_t);
  for (int i = 0; i < n_pts_t; i++) {
    fprintf(f, "1\n");
  }

  fprintf(f, "CELL_DATA %d\n", n_pts_t);
  fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int i = 0; i < n_pts[ipart]; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", g_num[ipart][i]);
    }
  }

  if (location != NULL) {
    fprintf(f, "FIELD FieldData 1\n");
    fprintf(f, "location 1 %d int\n", n_pts_t);
    for (int ipart = 0; ipart < n_part; ipart++) {
      for (int i = 0; i < n_pts[ipart]; i++) {
        fprintf(f, ""PDM_FMT_G_NUM"\n", location[ipart][i]);
      }
    }
  }

  fclose(f);
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

  PDM_g_num_t n_vtx_seg           = 10;
  double      length              = 1.;
  double      separation_x        = 2.;
  double      separation_y        = 0.;
  double      separation_z        = 0.;
  int         deform              = 0;
  double      tolerance           = 1e-6;
  double      marge               = 0.;
  int         n_part              = 1;
  int         post                = 0;
  int         extension_depth_tgt = 0;
  int         extension_depth_src = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  PDM_g_num_t n_pts = 10;
  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;
  int disable_uvw = 0;
  int use_tgt_nodes = 0;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_seg,
              &length,
              &separation_x,
              &separation_y,
              &separation_z,
              &deform,
              &tolerance,
              &marge,
              &n_part,
              &n_pts,
              &post,
              (int *) &part_method,
              &loc_method,
              &disable_uvw,
              &use_tgt_nodes,
              &extension_depth_tgt,
              &extension_depth_src);


  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  const double xmin = 0;
  const double ymin = 0;
  const double zmin = 0;

  /*
   *  Source cube
   */
  int          *src_n_cell        = NULL;
  int          *src_n_face        = NULL;
  int          *src_n_vtx         = NULL;
  int         **src_cell_face_idx = NULL;
  int         **src_cell_face     = NULL;
  int         **src_face_vtx_idx  = NULL;
  int         **src_face_vtx      = NULL;
  double      **src_vtx_coord     = NULL;
  PDM_g_num_t **src_cell_ln_to_gn = NULL;
  PDM_g_num_t **src_face_ln_to_gn = NULL;
  PDM_g_num_t **src_vtx_ln_to_gn  = NULL;

  _cube_mesh (comm,
              n_part,
              part_method,
              n_vtx_seg,
              xmin,
              ymin,
              zmin,
              length,
              0,//deform
              extension_depth_src,
              &src_n_cell,
              &src_n_face,
              &src_n_vtx,
              &src_cell_face_idx,
              &src_cell_face,
              &src_face_vtx_idx,
              &src_face_vtx,
              &src_vtx_coord,
              &src_cell_ln_to_gn,
              &src_face_ln_to_gn,
              &src_vtx_ln_to_gn);

  /*
   *  Target cube
   */
  int          *tgt_n_cell        = NULL;
  int          *tgt_n_face        = NULL;
  int          *tgt_n_vtx         = NULL;
  int         **tgt_cell_face_idx = NULL;
  int         **tgt_cell_face     = NULL;
  int         **tgt_face_vtx_idx  = NULL;
  int         **tgt_face_vtx      = NULL;
  double      **tgt_vtx_coord     = NULL;
  PDM_g_num_t **tgt_cell_ln_to_gn = NULL;
  PDM_g_num_t **tgt_face_ln_to_gn = NULL;
  PDM_g_num_t **tgt_vtx_ln_to_gn  = NULL;

  _cube_mesh (comm,
              n_part,
              part_method,
              n_vtx_seg,
              xmin + separation_x*length,
              ymin + separation_y*length,
              zmin + separation_z*length,
              length,
              deform,
              extension_depth_tgt,
              &tgt_n_cell,
              &tgt_n_face,
              &tgt_n_vtx,
              &tgt_cell_face_idx,
              &tgt_cell_face,
              &tgt_face_vtx_idx,
              &tgt_face_vtx,
              &tgt_vtx_coord,
              &tgt_cell_ln_to_gn,
              &tgt_face_ln_to_gn,
              &tgt_vtx_ln_to_gn);

  /*
   *  Mesh location structure initialization
   */
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED,
                                                            1,
                                                            comm,
                                                            PDM_OWNERSHIP_KEEP);

  PDM_mesh_location_reverse_results_enable(mesh_loc);

  /* Set target point cloud */
  PDM_mesh_location_n_part_cloud_set (mesh_loc,
                                      0,
                                      n_part);

  int          *n_tgt     = NULL;
  PDM_g_num_t **tgt_g_num = NULL;
  double      **tgt_coord = NULL;
  if (use_tgt_nodes) {
    n_tgt     = tgt_n_vtx;
    tgt_g_num = tgt_vtx_ln_to_gn;
    tgt_coord = tgt_vtx_coord;
  }
  else {
    n_tgt     = tgt_n_cell;
    tgt_g_num = tgt_cell_ln_to_gn;
    tgt_coord = (double **) malloc (sizeof(double *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      const int is_oriented = 1;
      double *cell_volume = (double *) malloc(sizeof(double) * n_tgt[ipart]);
      tgt_coord[ipart] = (double *) malloc(sizeof(double) * 3 * n_tgt[ipart]);

      PDM_geom_elem_polyhedra_properties (is_oriented,
                                          tgt_n_cell[ipart],
                                          tgt_n_face[ipart],
                                          tgt_face_vtx_idx[ipart],
                                          tgt_face_vtx[ipart],
                                          tgt_cell_face_idx[ipart],
                                          tgt_cell_face[ipart],
                                          tgt_n_vtx[ipart],
                                          tgt_vtx_coord[ipart],
                                          cell_volume,
                                          tgt_coord[ipart],
                                          NULL,
                                          NULL);

      free (cell_volume);
    }
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_mesh_location_cloud_set (mesh_loc,
                                 0,
                                 ipart,
                                 n_tgt[ipart],
                                 tgt_coord[ipart],
                                 tgt_g_num[ipart]);
  }


  /* Set source mesh */
  PDM_mesh_location_mesh_global_data_set (mesh_loc,
                                          n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_mesh_location_part_set (mesh_loc,
                                ipart,
                                src_n_cell[ipart],
                                src_cell_face_idx[ipart],
                                src_cell_face[ipart],
                                src_cell_ln_to_gn[ipart],
                                src_n_face[ipart],
                                src_face_vtx_idx[ipart],
                                src_face_vtx[ipart],
                                src_face_ln_to_gn[ipart],
                                src_n_vtx[ipart],
                                src_vtx_coord[ipart],
                                src_vtx_ln_to_gn[ipart]);
  }


  /* Set location parameters */
  PDM_mesh_location_tolerance_set (mesh_loc,
                                   tolerance);

  PDM_mesh_location_method_set (mesh_loc,
                                loc_method);

  /*
   *  Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  PDM_mesh_location_compute (mesh_loc);

  PDM_mesh_location_dump_times (mesh_loc);




  /*
   *  Check result from target PoV
   */
  PDM_g_num_t **tgt_location   = malloc (sizeof(PDM_g_num_t *) * n_part);
  double      **tgt_proj_coord = malloc (sizeof(double *)      * n_part);

  int n_wrong = 0;
  const PDM_g_num_t n_cell_seg = n_vtx_seg - 1;
  const double cell_side = length / ((double) n_cell_seg);

  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_located = PDM_mesh_location_n_located_get (mesh_loc,
                                                     0,//i_point_cloud,
                                                     ipart);

    int *located = PDM_mesh_location_located_get (mesh_loc,
                                                  0,//i_point_cloud,
                                                  ipart);

    int n_unlocated = PDM_mesh_location_n_unlocated_get (mesh_loc,
                                                         0,//i_point_cloud,
                                                         ipart);

    int *unlocated = PDM_mesh_location_unlocated_get (mesh_loc,
                                                      0,//i_point_cloud,
                                                      ipart);

    printf("[%d] part %d, n_located = %d, n_unlocated = %d\n", i_rank, ipart, n_located, n_unlocated);
    assert(n_located + n_unlocated == n_tgt[ipart]);

    tgt_location[ipart] = PDM_array_const_gnum (n_tgt[ipart], 0);
    tgt_proj_coord[ipart] = malloc (sizeof(double) * n_tgt[ipart] * 3);

    PDM_g_num_t *p_location    = NULL;
    double      *p_dist2  = NULL;
    double      *p_proj_coord  = NULL;
    PDM_mesh_location_point_location_get (mesh_loc,
                                          0,//i_point_cloud,
                                          ipart,
                                          &p_location,
                                          &p_dist2,
                                          &p_proj_coord);

    for (int j = 0; j < n_located; j++) {
      int i = located[j] - 1;
      tgt_location[ipart][i] = p_location[j];
      for (int k = 0; k < 3; k++) {
        tgt_proj_coord[ipart][3*i+k] = p_proj_coord[3*j+k];
      }
    }

    for (int j = 0; j < n_unlocated; j++) {
      int i = unlocated[j] - 1;
      tgt_location[ipart][i] = -1;
      for (int k = 0; k < 3; k++) {
        tgt_proj_coord[ipart][3*i+k] = tgt_coord[ipart][3*i+k];
      }
    }


    if (1) {//!deform) {

      for (int k1 = 0; k1 < n_located; k1++) {
        int ipt = located[k1] - 1;
        double *p = tgt_coord[ipart] + 3*ipt;

        int i = (int) floor (p[0] / cell_side);
        int j = (int) floor (p[1] / cell_side);
        int k = (int) floor (p[2] / cell_side);

        PDM_g_num_t box_gnum = 1 + i + n_cell_seg*(j + n_cell_seg*k);

        if (p[0] < -tolerance || p[0] > length + tolerance ||
            p[1] < -tolerance || p[1] > length + tolerance ||
            p[2] < -tolerance || p[2] > length + tolerance) {
          box_gnum = -1;
        }

        if (p_location[k1] != box_gnum) {
          double cell_min[3] = {cell_side * i,     cell_side * j,     cell_side * k};
          double cell_max[3] = {cell_side * (i+1), cell_side * (j+1), cell_side * (k+1)};

          double dist = HUGE_VAL;
          for (int idim = 0; idim < 3; idim++) {
            double _dist1 = PDM_ABS (p[idim] - cell_min[idim]);
            double _dist2 = PDM_ABS (p[idim] - cell_max[idim]);
            double _dist = PDM_MIN (_dist1, _dist2);
            dist = PDM_MIN (dist, _dist);
          }

          if (dist > tolerance) {
            n_wrong++;
          }
        }
      }



      for (int k1 = 0; k1 < n_unlocated; k1++) {
        int ipt = unlocated[k1] - 1;

        double x = tgt_coord[ipart][3*ipt];
        double y = tgt_coord[ipart][3*ipt+1];
        double z = tgt_coord[ipart][3*ipt+2];
        if (x >= xmin && x <= xmin + length &&
            y >= ymin && y <= ymin + length &&
            z >= zmin && z <= zmin + length) {
          n_wrong++;
        }
      }

    }
  }

  int g_n_wrong;
  PDM_MPI_Allreduce (&n_wrong, &g_n_wrong, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  if (i_rank == 0) {
    printf("Viewed from target: g_n_wrong = %d / "PDM_FMT_G_NUM"\n", g_n_wrong, (n_vtx_seg-1)*(n_vtx_seg-1)*(n_vtx_seg-1));
  }









  /*
   *  Check result from source PoV
   */
  if (1) {//!deform) {

    n_wrong = 0;

    for (int ipart = 0; ipart < n_part; ipart++) {
      int *cell_vtx_idx;
      int *cell_vtx;
      PDM_mesh_location_cell_vertex_get (mesh_loc,
                                         ipart,
                                         &cell_vtx_idx,
                                         &cell_vtx);

      int         *elt_pts_inside_idx;
      PDM_g_num_t *points_gnum;
      double      *points_coords;
      double      *points_uvw;
      int         *points_weights_idx;
      double      *points_weights;
      double      *points_dist2;
      double      *points_projected_coords;

      PDM_mesh_location_points_in_elt_get (mesh_loc,
                                           ipart,
                                           0,//i_point_cloud,
                                           &elt_pts_inside_idx,
                                           &points_gnum,
                                           &points_coords,
                                           &points_uvw,
                                           &points_weights_idx,
                                           &points_weights,
                                           &points_dist2,
                                           &points_projected_coords);

      for (int i = 0; i < src_n_cell[ipart]; i++) {

        PDM_g_num_t ck = (src_cell_ln_to_gn[ipart][i] - 1) / (n_cell_seg * n_cell_seg);
        PDM_g_num_t ci = (src_cell_ln_to_gn[ipart][i] - 1) % n_cell_seg;
        PDM_g_num_t cj = (src_cell_ln_to_gn[ipart][i] - 1 - ck*n_cell_seg*n_cell_seg) / n_cell_seg;

        for (int j = elt_pts_inside_idx[i]; j < elt_pts_inside_idx[i+1]; j++) {
          double *p = points_coords + 3*j;

          PDM_g_num_t pi = (PDM_g_num_t) floor (p[0] / cell_side);
          PDM_g_num_t pj = (PDM_g_num_t) floor (p[1] / cell_side);
          PDM_g_num_t pk = (PDM_g_num_t) floor (p[2] / cell_side);

          if (ci != pi || cj != pj || ck != pk) {

            double cell_min[3] = {cell_side * ci,     cell_side * cj,     cell_side * ck};
            double cell_max[3] = {cell_side * (ci+1), cell_side * (cj+1), cell_side * (ck+1)};

            double dist = HUGE_VAL;
            for (int idim = 0; idim < 3; idim++) {
              double _dist1 = PDM_ABS (p[idim] - cell_min[idim]);
              double _dist2 = PDM_ABS (p[idim] - cell_max[idim]);
              double _dist = PDM_MIN (_dist1, _dist2);
              dist = PDM_MIN (dist, _dist);
            }

            if (dist > tolerance) {
              //printf("!!! part %d, from source cell "PDM_FMT_G_NUM", point "PDM_FMT_G_NUM"\n", ipart, src_g_num[ipart][i], points_gnum[j]);
              n_wrong++;
            }
          }
        }
      }
    }


    PDM_MPI_Allreduce (&n_wrong, &g_n_wrong, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    if (i_rank == 0) {
      printf("Viewed from source: g_n_wrong = %d / "PDM_FMT_G_NUM"\n", g_n_wrong, (n_vtx_seg-1)*(n_vtx_seg-1)*(n_vtx_seg-1));
    }
  }




  if (post) {
    char filename[999];

    sprintf(filename, "tgt_location_%3.3d.vtk", i_rank);
    _export_point_cloud (filename,
                         n_part,
                         n_tgt,
                         tgt_coord,
                         tgt_g_num,
                         tgt_location);

    sprintf(filename, "tgt_proj_coord_%3.3d.vtk", i_rank);

    _export_point_cloud (filename,
                         n_part,
                         n_tgt,
                         tgt_proj_coord,
                         tgt_g_num,
                         tgt_location);
  }


  /*
   *  Free memory
   */
  for (int ipart = 0; ipart < n_part; ipart++) {
    free(src_cell_face_idx[ipart]);
    free(src_cell_face    [ipart]);
    free(src_face_vtx_idx [ipart]);
    free(src_face_vtx     [ipart]);
    free(src_vtx_coord    [ipart]);
    free(src_cell_ln_to_gn[ipart]);
    free(src_face_ln_to_gn[ipart]);
    free(src_vtx_ln_to_gn [ipart]);

    free(tgt_cell_face_idx[ipart]);
    free(tgt_cell_face    [ipart]);
    free(tgt_face_vtx_idx [ipart]);
    free(tgt_face_vtx     [ipart]);
    free(tgt_vtx_coord    [ipart]);
    free(tgt_cell_ln_to_gn[ipart]);
    free(tgt_face_ln_to_gn[ipart]);
    free(tgt_vtx_ln_to_gn [ipart]);

    free (tgt_location[ipart]);
    free (tgt_proj_coord[ipart]);

    if (!use_tgt_nodes) {
      free(tgt_coord[ipart]);
    }
  }

  free(src_n_cell       );
  free(src_n_face       );
  free(src_n_vtx        );
  free(src_cell_face_idx);
  free(src_cell_face    );
  free(src_face_vtx_idx );
  free(src_face_vtx     );
  free(src_vtx_coord    );
  free(src_cell_ln_to_gn);
  free(src_face_ln_to_gn);
  free(src_vtx_ln_to_gn );

  free(tgt_n_cell       );
  free(tgt_n_face       );
  free(tgt_n_vtx        );
  free(tgt_cell_face_idx);
  free(tgt_cell_face    );
  free(tgt_face_vtx_idx );
  free(tgt_face_vtx     );
  free(tgt_vtx_coord    );
  free(tgt_cell_ln_to_gn);
  free(tgt_face_ln_to_gn);
  free(tgt_vtx_ln_to_gn );

  free (tgt_location);
  free (tgt_proj_coord);

  if (!use_tgt_nodes) {
    free(tgt_coord);
  }

  PDM_mesh_location_free (mesh_loc);
                          

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return g_n_wrong;
}

