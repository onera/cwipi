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
#include "pdm_unique.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_to_part.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*
 *
 * On a mesh_a et mesh_b et mesh_c
 *  Il faut être bijectif
 * 1/ mesh_a =
 *
 */


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
           int                         *extension_depth_src,
           PDM_Mesh_nodal_elt_t        *elt_type)
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
    else if (strcmp(argv[i], "-elt_type") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = atoi(argv[i]);
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
 const PDM_MPI_Comm            comm,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 const PDM_g_num_t             n_vtx_seg,
 const double                  xmin,
 const double                  ymin,
 const double                  zmin,
 const double                  length,
 const int                     deform,
 const int                     part_extension_depth,
 const PDM_Mesh_nodal_elt_t    elt_type,
 int                         **pn_cell,
 int                         **pn_face,
 int                         **pn_vtx,
 int                         **pn_cell_ext,
 int                         **pn_face_ext,
 int                         **pn_vtx_ext,
 int                        ***pcell_face_idx,
 int                        ***pcell_face,
 int                        ***pface_vtx_idx,
 int                        ***pface_vtx,
 double                     ***pvtx_coord,
 PDM_g_num_t                ***pcell_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn,
 PDM_g_num_t                ***pvtx_ln_to_gn,
 int                        ***pface_group_idx,
 int                        ***pface_group,
 PDM_g_num_t                ***pface_group_ln_to_gn
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
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];


  PDM_dmesh_nodal_to_dmesh_t *dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  PDM_dmesh_t *dmesh2 = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh2);

  int dn_cell = 0;
  int dn_face = 0;
  int dn_edge = -1;
  int n_face_group = 0;

  int         *dcell_face_idx  = NULL;
  PDM_g_num_t *dcell_face      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  int         *dface_cell_idx  = NULL;
  PDM_g_num_t *dface_cell      = NULL;
  PDM_g_num_t *dface_group     = NULL;
  int         *dface_group_idx = NULL;

  dn_face = PDM_dmesh_connectivity_get(dmesh2, PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                       &dface_vtx,
                                       &dface_vtx_idx,
                                       PDM_OWNERSHIP_KEEP);

  PDM_dmesh_connectivity_get(dmesh2, PDM_CONNECTIVITY_TYPE_FACE_CELL,
                             &dface_cell,
                             &dface_cell_idx,
                             PDM_OWNERSHIP_KEEP);
  assert(dface_cell_idx == NULL);

  dn_cell = PDM_dmesh_connectivity_get(dmesh2, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                       &dcell_face,
                                       &dcell_face_idx,
                                       PDM_OWNERSHIP_KEEP);

  n_face_group = PDM_dmesh_bound_get(dmesh2,
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
                                         dn_edge,
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


  *pn_cell        = (int          *) malloc(sizeof(int         * ) * n_part);
  *pn_face        = (int          *) malloc(sizeof(int         * ) * n_part);
  *pn_vtx         = (int          *) malloc(sizeof(int         * ) * n_part);
  *pcell_face_idx = (int         **) malloc(sizeof(int         **) * n_part);
  *pcell_face     = (int         **) malloc(sizeof(int         **) * n_part);
  *pface_vtx_idx  = (int         **) malloc(sizeof(int         **) * n_part);
  *pface_vtx      = (int         **) malloc(sizeof(int         **) * n_part);
  *pcell_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pface_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_coord     = (double      **) malloc(sizeof(double      **) * n_part);


  *pn_cell_ext    = (int *)          malloc(sizeof(int *)          * n_part);
  *pn_face_ext    = (int *)          malloc(sizeof(int *)          * n_part);
  *pn_vtx_ext     = (int *)          malloc(sizeof(int *)          * n_part);

  *pface_group_idx      = (int         **) malloc(sizeof(int         **) * n_part );
  *pface_group          = (int         **) malloc(sizeof(int         **) * n_part );
  *pface_group_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part );

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
    double      *ext_vtx_coord           = NULL;
    PDM_g_num_t *ext_vtx_ln_to_gn        = NULL;
    int         *ext_cell_face           = NULL;
    int         *ext_cell_face_idx       = NULL;
    PDM_g_num_t *ext_cell_ln_to_gn       = NULL;
    int         *ext_face_vtx            = NULL;
    int         *ext_face_vtx_idx        = NULL;
    PDM_g_num_t *ext_face_ln_to_gn       = NULL;
    int         *ext_face_group          = NULL;
    int         *ext_face_group_idx      = NULL;
    PDM_g_num_t *ext_face_group_ln_to_gn = NULL;

    int n_ext_face_group = 0;
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

      int n_group = PDM_part_extension_group_get(part_ext,
                                                 0,
                                                 i_part,
                                                 PDM_MESH_ENTITY_FACE,
                                                 &ext_face_group,
                                                 &ext_face_group_idx,
                                                 &ext_face_group_ln_to_gn);

      n_ext_face_group = ext_face_group_idx[n_group];

    }


    *(pn_cell)[i_part] = n_cell + n_ext_cell;
    *(pn_face)[i_part] = n_face + n_ext_face;
    *(pn_vtx )[i_part] = n_vtx  + n_ext_vtx;

    *(pn_cell_ext)[i_part] = n_ext_cell;
    *(pn_face_ext)[i_part] = n_ext_face;
    *(pn_vtx_ext )[i_part] = n_ext_vtx;

    *(pface_group_idx     )[i_part] = malloc( (n_groups+1)                                * sizeof(int        ));
    *(pface_group         )[i_part] = malloc( (face_group_idx[n_groups]+n_ext_face_group) * sizeof(int        ));
    *(pface_group_ln_to_gn)[i_part] = malloc( (face_group_idx[n_groups]+n_ext_face_group) * sizeof(PDM_g_num_t));

    int         *_pface_group_idx      = *(pface_group_idx     )[i_part];
    int         *_pface_group          = *(pface_group         )[i_part];
    PDM_g_num_t *_pface_group_ln_to_gn = *(pface_group_ln_to_gn)[i_part];

    _pface_group_idx[0] = 0;
    for(int i_group = 0; i_group < n_groups+1; ++i_group) {
      _pface_group_idx[i_group] = 0;
    }

    for(int i_group = 0; i_group < n_groups; ++i_group) {
      _pface_group_idx[i_group+1] = _pface_group_idx[i_group] + face_group_idx[i_group+1] - face_group_idx[i_group];
      if(ext_face_group != NULL) {
        int dn = ext_face_group_idx[i_group+1] - ext_face_group_idx[i_group];
        _pface_group_idx[i_group+1] += dn;
      }
    }

    for(int i_group = 0; i_group < n_groups; ++i_group) {

      int idx_write = _pface_group_idx[i_group];
      for(int idx_face = face_group_idx[i_group]; idx_face < face_group_idx[i_group+1]; ++idx_face) {
        _pface_group         [idx_write] = face_group         [idx_face];
        _pface_group_ln_to_gn[idx_write] = face_group_ln_to_gn[idx_face];
        idx_write++;
      }

      if(ext_face_group != NULL) {
        for(int idx_face = ext_face_group_idx[i_group]; idx_face < ext_face_group_idx[i_group+1]; ++idx_face) {
          _pface_group         [idx_write] = ext_face_group         [idx_face];
          _pface_group_ln_to_gn[idx_write] = ext_face_group_ln_to_gn[idx_face];
          idx_write++;
        }
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_int(_pface_group_idx, n_groups+1                , "_pface_group_idx ::");
      PDM_log_trace_array_int(_pface_group    , _pface_group_idx[n_groups], "_pface_group     ::");
    }


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
  for (int i_part = 0; i_part < n_part; i_part++) {
    n_pts_t += n_pts[i_part];
  }

  fprintf(f, "POINTS %d double\n", n_pts_t);
  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i = 0; i < n_pts[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        fprintf(f, "%f ", coord[i_part][3*i + j]);
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
  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i = 0; i < n_pts[i_part]; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", g_num[i_part][i]);
    }
  }

  if (location != NULL) {
    fprintf(f, "FIELD FieldData 1\n");
    fprintf(f, "location 1 %d int\n", n_pts_t);
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < n_pts[i_part]; i++) {
        fprintf(f, ""PDM_FMT_G_NUM"\n", location[i_part][i]);
      }
    }
  }

  fclose(f);
}


// /**
//  *   Compute volume and center_cell
//  */
// static
// void
// _compute_centers
// (
//   int           n_part,
//   int          *n_cell,
//   int          *n_face,
//   int          *n_vtx,
//   int         **cell_face_idx,
//   int         **cell_face,
//   int         **face_vtx_idx,
//   int         **face_vtx,
//   double      **vtx_coord,
//   double     ***cell_center_out,
//   double     ***cell_volume_out,
//   double     ***face_center_out
// )
// {

//   double **cell_volume = (double ** ) malloc( n_part * sizeof(double));
//   double **cell_center = (double ** ) malloc( n_part * sizeof(double));
//   double **face_center = (double ** ) malloc( n_part * sizeof(double));
//   for (int i_part = 0; i_part < n_part; i_part++) {
//     const int is_oriented = 1;
//     cell_volume[i_part] = (double *) malloc(sizeof(double) *     n_cell[i_part]);
//     cell_center[i_part] = (double *) malloc(sizeof(double) * 3 * n_cell[i_part]);

//     face_center[i_part] = (double *) malloc(sizeof(double) * 3 * n_face[i_part]);

//     double* face_normal = (double *) malloc(sizeof(double) * 3 * n_face[i_part]);

//     PDM_geom_elem_polyhedra_properties(is_oriented,
//                                        n_cell       [i_part],
//                                        n_face       [i_part],
//                                        face_vtx_idx [i_part],
//                                        face_vtx     [i_part],
//                                        cell_face_idx[i_part],
//                                        cell_face    [i_part],
//                                        n_vtx        [i_part],
//                                        vtx_coord    [i_part],
//                                        cell_volume  [i_part],
//                                        cell_center  [i_part],
//                                        NULL,
//                                        NULL);

//     PDM_geom_elem_polygon_properties(n_face      [i_part],
//                                      face_vtx_idx[i_part],
//                                      face_vtx    [i_part],
//                                      vtx_coord   [i_part],
//                                      face_normal,
//                                      face_center [i_part],
//                                      NULL,
//                                      NULL);
//     free(face_normal);

//   }

//   *cell_center_out = cell_center;
//   *cell_volume_out = cell_volume;
//   *face_center_out = face_center;

// }



// /**
//  *
//  *
//  *
//  */
// static
// void
// _prepare_target_cloud
// (
//   int            n_part,
//   int            n_group,
//   int           *n_cell,
//   int           *n_face,
//   int           *n_vtx,
//   PDM_g_num_t  **cell_ln_to_gn,
//   PDM_g_num_t  **face_ln_to_gn,
//   PDM_g_num_t  **vtx_ln_to_gn,
//   PDM_g_num_t  **face_group_ln_to_gn,
//   int          **cell_face_idx,
//   int          **cell_face,
//   int          **face_vtx_idx,
//   int          **face_vtx,
//   int          **face_group_idx,
//   int          **face_group,
//   double       **vtx_coord,
//   double       **cell_center,
//   double       **face_center,
//   int            extract_center_depth,
//   int            extract_bnd_faces,
//   int          **n_extract_cell,
//   int          **n_extract_face,
//   PDM_g_num_t ***extract_center_ln_to_gn,
//   double      ***extract_center_coord,
//   PDM_g_num_t ***extract_face_bnd_ln_to_gn,
//   double      ***extract_face_bnd_coord
// )
// {
//   PDM_UNUSED(n_vtx);
//   PDM_UNUSED(face_ln_to_gn);
//   PDM_UNUSED(vtx_ln_to_gn);
//   PDM_UNUSED(face_vtx_idx);
//   PDM_UNUSED(face_vtx);
//   PDM_UNUSED(vtx_coord);


//   /*
//    * Extract cell with a given depth and also face center
//    */
//   *n_extract_cell          = malloc(n_part * sizeof(int          ));
//   *extract_center_coord    = malloc(n_part * sizeof(double      *));
//   *extract_center_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));

//   *n_extract_face            = malloc(n_part * sizeof(int          ));
//   *extract_face_bnd_coord    = malloc(n_part * sizeof(double      *));
//   *extract_face_bnd_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));

//   int          *_n_extract_cell          = *n_extract_cell;
//   double      **_extract_center_coord    = *extract_center_coord;
//   PDM_g_num_t **_extract_center_ln_to_gn = *extract_center_ln_to_gn;

//   int          *_n_extract_face            = *n_extract_face;
//   double      **_extract_face_bnd_coord    = *extract_face_bnd_coord;
//   PDM_g_num_t **_extract_face_bnd_ln_to_gn = *extract_face_bnd_ln_to_gn;

//   /* Management of void */
//   for (int i_part = 0; i_part < n_part; i_part++) {
//     _n_extract_cell[i_part] = 0;
//     _n_extract_face[i_part] = 0;
//     _extract_center_coord     [i_part] = NULL;
//     _extract_center_ln_to_gn  [i_part] = NULL;
//     _extract_face_bnd_coord   [i_part] = NULL;
//     _extract_face_bnd_ln_to_gn[i_part] = NULL;
//   }

//   if(extract_center_depth > 0) {
//     for (int i_part = 0; i_part < n_part; i_part++) {

//       int* face_cell_idx = NULL;
//       int* face_cell     = NULL;
//       PDM_connectivity_transpose(n_cell       [i_part],
//                                  n_face       [i_part],
//                                  cell_face_idx[i_part],
//                                  cell_face    [i_part],
//                                  &face_cell_idx,
//                                  &face_cell);

//       int* extract_cell = malloc(n_face[i_part] * sizeof(int)); // Surallocated buyt needed

//       _n_extract_cell[i_part] = 0;

//       for(int i_group = 0; i_group < n_group; ++i_group) {
//         for(int idx_face = face_group_idx[i_part][i_group]; idx_face < face_group_idx[i_part][i_group+1]; ++idx_face) {
//           int i_face = face_group[i_part][idx_face]-1;
//           for(int idx_cell = face_cell_idx[i_face]; idx_cell < face_cell_idx[i_face+1]; ++idx_cell) {
//             extract_cell[_n_extract_cell[i_part]++] = PDM_ABS(face_cell[idx_cell]);
//           }
//         }
//       }

//       // sort and extract coordinates
//       _n_extract_cell[i_part] = PDM_inplace_unique(extract_cell, 0, _n_extract_cell[i_part]-1);

//       extract_cell = realloc(extract_cell, _n_extract_cell[i_part] * sizeof(int));

//       _extract_center_coord   [i_part] = (double      *) malloc(3 * _n_extract_cell[i_part] * sizeof(double     ));
//       _extract_center_ln_to_gn[i_part] = (PDM_g_num_t *) malloc(    _n_extract_cell[i_part] * sizeof(PDM_g_num_t));

//       for(int i = 0; i < _n_extract_cell[i_part]; ++i) {
//         int i_cell = extract_cell[i]-1;
//         _extract_center_coord   [i_part][3*i  ] = cell_center[i_part][3*i_cell  ];
//         _extract_center_coord   [i_part][3*i+1] = cell_center[i_part][3*i_cell+1];
//         _extract_center_coord   [i_part][3*i+2] = cell_center[i_part][3*i_cell+2];
//         _extract_center_ln_to_gn[i_part][i] = cell_ln_to_gn[i_part][i_cell];
//       }

//       free(face_cell_idx);
//       free(face_cell);
//       free(extract_cell);

//     }
//   }


//   if(extract_bnd_faces == 1) {

//     for (int i_part = 0; i_part < n_part; i_part++) {
//       _n_extract_face[i_part] = 0;

//       _extract_face_bnd_coord   [i_part] = malloc(face_group_idx[i_part][n_group] * sizeof(double     ));
//       _extract_face_bnd_ln_to_gn[i_part] = malloc(face_group_idx[i_part][n_group] * sizeof(PDM_g_num_t));


//       /* Dans cette exemple on considère que toutes les frontières sont BCOVerlap */
//       for(int i_group = 0; i_group < n_group; ++i_group) {
//         for(int idx_face = face_group_idx[i_part][i_group]; idx_face < face_group_idx[i_part][i_group+1]; ++idx_face) {
//           int i_face = face_group[i_part][idx_face]-1;

//           int idx_write = _n_extract_face[i_part]++;

//           _extract_face_bnd_coord   [i_part][3*idx_write  ] = face_center[i_part][3*i_face  ];
//           _extract_face_bnd_coord   [i_part][3*idx_write+1] = face_center[i_part][3*i_face+1];
//           _extract_face_bnd_coord   [i_part][3*idx_write+2] = face_center[i_part][3*i_face+2];
//           // _extract_face_bnd_ln_to_gn[i_part][idx_write] = face_ln_to_gn[i_part][i_face]; // Change rien
//           _extract_face_bnd_ln_to_gn[i_part][idx_write] = face_group_ln_to_gn[i_part][idx_face];

//         }
//       }
//       assert(_n_extract_face[i_part] == face_group_idx[i_part][n_group]);
//     }
//   }


// }



/**
 *
 * \brief  Visu
 *
 */

static void
_visu
(
 const char      *name_chr,
 const int        n_part,
 int              n_cell[],
 int              n_face[],
 int              n_vtx[],
 int             *cell_face_idx[],
 int             *cell_face[],
 PDM_g_num_t     *cell_ln_to_gn[],
 int             *face_vtx_idx[],
 int             *face_vtx[],
 PDM_g_num_t     *face_ln_to_gn[],
 double          *vtx[],
 PDM_g_num_t     *vtx_ln_to_gn[],
 double          *field[],
 double          *block_interp[]
)
{
  PDM_UNUSED(face_ln_to_gn);

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);


  PDM_writer_t *id_cs = PDM_writer_create("Ensight",
                                          PDM_WRITER_FMT_ASCII,
                                          PDM_WRITER_TOPO_CST,
                                          PDM_WRITER_OFF,
                                          "test_chimera",
                                          name_chr,
                                          PDM_MPI_COMM_WORLD,
                                          PDM_IO_KIND_MPI_SIMPLE,
                                          1.,
                                          NULL);

  /* Creation de la geometrie */

  int id_geom = PDM_writer_geom_create(id_cs,
                                       "cube_geom",
                                       n_part);

  int *n_part_procs = (int *) malloc(sizeof(int) * n_rank);

  PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
                     (void *) n_part_procs, 1, PDM_MPI_INT,
                     PDM_MPI_COMM_WORLD);

  int *distrib_part = (int *) malloc(sizeof(int) * (n_rank + 1));

  distrib_part[0] = 0;
  for (int i = 0; i < n_rank; i++) {
    distrib_part[i+1] = distrib_part[i] + n_part_procs[i];
  }

  free(n_part_procs);

  /* Creation des variables */

  int id_var_num_part = PDM_writer_var_create(id_cs,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");
  int id_var_field = PDM_writer_var_create(id_cs,
                                           PDM_WRITER_OFF,
                                           PDM_WRITER_VAR_SCALAR,
                                           PDM_WRITER_VAR_ELEMENTS,
                                           "field");

  int id_var_blk_interp = PDM_writer_var_create(id_cs,
                                                PDM_WRITER_OFF,
                                                PDM_WRITER_VAR_SCALAR,
                                                PDM_WRITER_VAR_ELEMENTS,
                                                "blk_interp");

  PDM_real_t **val_num_part = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
  int **face_vtx_n = (int **) malloc(sizeof(int *) * n_part);
  int **cell_face_n = (int **) malloc(sizeof(int *) * n_part);
  PDM_writer_step_beg (id_cs, 0.);
  for (int i_part = 0; i_part < n_part; i_part++) {

    face_vtx_n[i_part] = (int *) malloc(sizeof(int) * n_face[i_part]);
    cell_face_n[i_part] = (int *) malloc(sizeof(int) * n_cell[i_part]);

    for (int i = 0; i < n_face[i_part]; i++) {
      face_vtx_n[i_part][i] = face_vtx_idx[i_part][i+1] - face_vtx_idx[i_part][i];
    }

    for (int i = 0; i < n_cell[i_part]; i++) {
      cell_face_n[i_part][i] = cell_face_idx[i_part][i+1] - cell_face_idx[i_part][i];
    }

    PDM_writer_geom_coord_set(id_cs,
                              id_geom,
                              i_part,
                              n_vtx[i_part],
                              vtx[i_part],
                              vtx_ln_to_gn[i_part],
                              PDM_OWNERSHIP_USER);

    /* Construction de la connectivite pour sortie graphique */

    PDM_writer_geom_cell3d_cellface_add (id_cs,
                                         id_geom,
                                         i_part,
                                         n_cell       [i_part],
                                         n_face       [i_part],
                                         face_vtx_idx [i_part],
                                         face_vtx_n   [i_part],
                                         face_vtx     [i_part],
                                         cell_face_idx[i_part],
                                         cell_face_n  [i_part],
                                         cell_face    [i_part],
                                         cell_ln_to_gn[i_part]);

  }

  PDM_writer_geom_write(id_cs,
                        id_geom);

  for (int i_part = 0; i_part < n_part; i_part++) {

    val_num_part[i_part] = (double *) malloc(sizeof(double) * n_cell[i_part]);
    for (int i = 0; i < n_cell[i_part]; i++) {
      val_num_part[i_part][i] = i_part + 1 + distrib_part[i_rank];
    }
    // PDM_log_trace_array_double(val_num_part[i_part], n_cell[i_part], "val_num_part :: ");

    PDM_writer_var_set(id_cs,
                       id_var_num_part,
                       id_geom,
                       i_part,
                       val_num_part[i_part]);

    PDM_writer_var_set(id_cs,
                       id_var_field,
                       id_geom,
                       i_part,
                       field[i_part]);

    PDM_writer_var_set(id_cs,
                       id_var_blk_interp,
                       id_geom,
                       i_part,
                       block_interp[i_part]);
  }

  PDM_writer_var_write(id_cs, id_var_field);
  PDM_writer_var_write(id_cs, id_var_num_part);
  PDM_writer_var_write(id_cs, id_var_blk_interp);

  // PDM_writer_var_free(id_cs,
  //                     id_var_num_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free (face_vtx_n  [i_part]);
    free (cell_face_n [i_part]);
    free (val_num_part[i_part]);
  }
  free (distrib_part);
  free (face_vtx_n);
  free (val_num_part);
  free (cell_face_n);

  // PDM_writer_geom_free(id_cs,
  //                      id_geom);

  PDM_writer_step_end(id_cs);
  PDM_writer_free(id_cs);

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
  double      separation_x        = 0.;
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

  PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_HEXA8;
  //  2 -> tria
  //  3 -> quad
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  // 11 -> tria_ho
  // 12 -> quad_ho
  // 13 -> tetra_ho
  // 14 -> pyramid_ho
  // 15 -> prism_ho
  // 16 -> hexa_ho

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
              &extension_depth_src,
              &elt_type);


  printf("separation_x : %12.5e\n", separation_x);

  assert(//elt_type == PDM_MESH_NODAL_BAR2     ||
         //elt_type == PDM_MESH_NODAL_TRIA3    ||
         //elt_type == PDM_MESH_NODAL_QUAD4    ||
         elt_type == PDM_MESH_NODAL_TETRA4   ||
         elt_type == PDM_MESH_NODAL_PYRAMID5 ||
         elt_type == PDM_MESH_NODAL_PRISM6   ||
         elt_type == PDM_MESH_NODAL_HEXA8);


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
  int          *src_n_cell         = NULL;
  int          *src_n_face         = NULL;
  int          *src_n_vtx          = NULL;
  int          *src_n_cell_ext     = NULL;
  int          *src_n_face_ext     = NULL;
  int          *src_n_vtx_ext      = NULL;
  int         **src_cell_face_idx  = NULL;
  int         **src_cell_face      = NULL;
  int         **src_face_vtx_idx   = NULL;
  int         **src_face_vtx       = NULL;
  double      **src_vtx_coord      = NULL;
  PDM_g_num_t **src_cell_ln_to_gn  = NULL;
  PDM_g_num_t **src_face_ln_to_gn  = NULL;
  PDM_g_num_t **src_vtx_ln_to_gn   = NULL;

  int         **src_face_group_idx = NULL;
  int         **src_face_group     = NULL;
  PDM_g_num_t **src_group_ln_to_gn = NULL;

  int n_group = 6;
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
              elt_type,
              &src_n_cell,
              &src_n_face,
              &src_n_vtx,
              &src_n_cell_ext,
              &src_n_face_ext,
              &src_n_vtx_ext,
              &src_cell_face_idx,
              &src_cell_face,
              &src_face_vtx_idx,
              &src_face_vtx,
              &src_vtx_coord,
              &src_cell_ln_to_gn,
              &src_face_ln_to_gn,
              &src_vtx_ln_to_gn,
              &src_face_group_idx,
              &src_face_group,
              &src_group_ln_to_gn);

  /*
   *  Target cube
   */
  int          *tgt_n_cell         = NULL;
  int          *tgt_n_face         = NULL;
  int          *tgt_n_vtx          = NULL;
  int          *tgt_n_cell_ext     = NULL;
  int          *tgt_n_face_ext     = NULL;
  int          *tgt_n_vtx_ext      = NULL;
  int         **tgt_cell_face_idx  = NULL;
  int         **tgt_cell_face      = NULL;
  int         **tgt_face_vtx_idx   = NULL;
  int         **tgt_face_vtx       = NULL;
  double      **tgt_vtx_coord      = NULL;
  PDM_g_num_t **tgt_cell_ln_to_gn  = NULL;
  PDM_g_num_t **tgt_face_ln_to_gn  = NULL;
  PDM_g_num_t **tgt_vtx_ln_to_gn   = NULL;
  int         **tgt_face_group_idx = NULL;
  int         **tgt_face_group     = NULL;
  PDM_g_num_t **tgt_group_ln_to_gn = NULL;

  _cube_mesh (comm,
              n_part,
              part_method,
              n_vtx_seg,
              xmin + separation_x*length,
              ymin + separation_y*length,
              zmin + separation_z*length,
              length/4,
              deform,
              extension_depth_tgt,
              PDM_MESH_NODAL_HEXA8,
              &tgt_n_cell,
              &tgt_n_face,
              &tgt_n_vtx,
              &tgt_n_cell_ext,
              &tgt_n_face_ext,
              &tgt_n_vtx_ext,
              &tgt_cell_face_idx,
              &tgt_cell_face,
              &tgt_face_vtx_idx,
              &tgt_face_vtx,
              &tgt_vtx_coord,
              &tgt_cell_ln_to_gn,
              &tgt_face_ln_to_gn,
              &tgt_vtx_ln_to_gn,
              &tgt_face_group_idx,
              &tgt_face_group,
              &tgt_group_ln_to_gn);

  /*
   *  Mesh location structure initialization
   */
  // int n_point_cloud = 1; // 1 : adjacent cell center | 2 : center_interface
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED,
                                                            1,
                                                            comm,
                                                            PDM_OWNERSHIP_KEEP);

  PDM_mesh_location_reverse_results_enable(mesh_loc);

  /* Set target point cloud */
  PDM_mesh_location_n_part_cloud_set (mesh_loc,
                                      0,
                                      n_part);


  /*
   * Extraction des surfaces + centre voisins
   */
  // int tgt_have_adjacent_cell_center = 1;
  // int tgt_have_face_center          = 1;

  int          *n_tgt     = malloc(n_part * sizeof(int          ));
  PDM_g_num_t **tgt_g_num = malloc(n_part * sizeof(PDM_g_num_t *));
  double      **tgt_coord = NULL;
  if (use_tgt_nodes) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      n_tgt[i_part]     = tgt_n_vtx[i_part];
    }
    tgt_g_num = tgt_vtx_ln_to_gn;
    tgt_coord = tgt_vtx_coord;
  }
  else {
    // n_tgt     = tgt_n_cell;
    tgt_coord = (double **) malloc (sizeof(double *) * n_part);


    for (int i_part = 0; i_part < n_part; i_part++) {
      const int is_oriented = 1;
      double *cell_volume = (double *) malloc(sizeof(double) *     tgt_n_cell[i_part]);
      double *cell_center = (double *) malloc(sizeof(double) * 3 * tgt_n_cell[i_part]);

      PDM_geom_elem_polyhedra_properties (is_oriented,
                                          tgt_n_cell[i_part],
                                          tgt_n_face[i_part],
                                          tgt_face_vtx_idx[i_part],
                                          tgt_face_vtx[i_part],
                                          tgt_cell_face_idx[i_part],
                                          tgt_cell_face[i_part],
                                          tgt_n_vtx[i_part],
                                          tgt_vtx_coord[i_part],
                                          cell_volume,
                                          cell_center,
                                          NULL,
                                          NULL);

      free (cell_volume);
      n_tgt[i_part] = 0;

      int* extract_cell = malloc(tgt_n_face[i_part] * sizeof(int));
      int* face_cell_idx = NULL;
      int* face_cell     = NULL;
      PDM_connectivity_transpose(tgt_n_cell       [i_part],
                                 tgt_n_face       [i_part],
                                 tgt_cell_face_idx[i_part],
                                 tgt_cell_face    [i_part],
                                 &face_cell_idx,
                                 &face_cell);

      for(int i_group = 0; i_group < n_group; ++i_group) {
        for(int idx_face = tgt_face_group_idx[i_part][i_group]; idx_face < tgt_face_group_idx[i_part][i_group+1]; ++idx_face) {
          int i_face = tgt_face_group[i_part][idx_face]-1;
          for(int idx_cell = face_cell_idx[i_face]; idx_cell < face_cell_idx[i_face+1]; ++idx_cell) {
            extract_cell[n_tgt[i_part]++] = PDM_ABS(face_cell[idx_cell]);
          }
        }
      }

      // sort and extract coordinates
      n_tgt[i_part] = PDM_inplace_unique(extract_cell, 0, n_tgt[i_part]-1);

      extract_cell = realloc(extract_cell, n_tgt[i_part] * sizeof(int));
      tgt_coord[i_part] = (double      *) malloc(3 * n_tgt[i_part] * sizeof(double));
      tgt_g_num[i_part] = (PDM_g_num_t *) malloc(    n_tgt[i_part] * sizeof(PDM_g_num_t));

      if(0 == 1) {
        PDM_log_trace_array_int(extract_cell, n_tgt[i_part], "extract_cell :: ");
      }

      for(int i = 0; i < n_tgt[i_part]; ++i) {
        int i_cell = extract_cell[i]-1;
        tgt_coord[i_part][3*i  ] = cell_center[3*i_cell  ];
        tgt_coord[i_part][3*i+1] = cell_center[3*i_cell+1];
        tgt_coord[i_part][3*i+2] = cell_center[3*i_cell+2];
        tgt_g_num[i_part][i] = tgt_cell_ln_to_gn[i_part][i_cell];
      }


      free(cell_center);
      free(extract_cell);
      free(face_cell_idx);
      free(face_cell);

    }
  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_mesh_location_cloud_set (mesh_loc,
                                 0,
                                 i_part,
                                 n_tgt[i_part],
                                 tgt_coord[i_part],
                                 tgt_g_num[i_part]);
  }


  /* Set source mesh */
  PDM_mesh_location_mesh_global_data_set (mesh_loc,
                                          n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    // Only include interior
    int n_cell_without_ext = src_n_cell[i_part]-src_n_cell_ext[i_part];
    int n_face_without_ext = src_n_face[i_part]-src_n_face_ext[i_part];
    int n_vtx_without_ext  = src_n_vtx [i_part]-src_n_vtx_ext [i_part];

    PDM_mesh_location_part_set (mesh_loc,
                                i_part,
                                n_cell_without_ext,
                                src_cell_face_idx[i_part],
                                src_cell_face    [i_part],
                                src_cell_ln_to_gn[i_part],
                                n_face_without_ext,
                                src_face_vtx_idx [i_part],
                                src_face_vtx     [i_part],
                                src_face_ln_to_gn[i_part],
                                n_vtx_without_ext,
                                src_vtx_coord    [i_part],
                                src_vtx_ln_to_gn [i_part]);
  }


  /* Set location parameters */
  PDM_mesh_location_tolerance_set (mesh_loc, tolerance);

  PDM_mesh_location_method_set (mesh_loc, loc_method);

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


  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_located = PDM_mesh_location_n_located_get (mesh_loc,
                                                     0,//i_point_cloud,
                                                     i_part);

    int *located = PDM_mesh_location_located_get (mesh_loc,
                                                  0,//i_point_cloud,
                                                  i_part);

    int n_unlocated = PDM_mesh_location_n_unlocated_get (mesh_loc,
                                                         0,//i_point_cloud,
                                                         i_part);

    int *unlocated = PDM_mesh_location_unlocated_get (mesh_loc,
                                                      0,//i_point_cloud,
                                                      i_part);

    printf("[%d] part %d, n_located = %d, n_unlocated = %d\n", i_rank, i_part, n_located, n_unlocated);
    assert(n_located + n_unlocated == n_tgt[i_part]);

    tgt_location  [i_part] = PDM_array_const_gnum (n_tgt[i_part], 0);
    tgt_proj_coord[i_part] = malloc (sizeof(double) * n_tgt[i_part] * 3);

    PDM_g_num_t *p_location    = NULL;
    double      *p_dist2  = NULL;
    double      *p_proj_coord  = NULL;
    PDM_mesh_location_point_location_get (mesh_loc,
                                          0,//i_point_cloud,
                                          i_part,
                                          &p_location,
                                          &p_dist2,
                                          &p_proj_coord);

    for (int j = 0; j < n_located; j++) {
      int i = located[j] - 1;
      tgt_location[i_part][i] = p_location[j];
      for (int k = 0; k < 3; k++) {
        tgt_proj_coord[i_part][3*i+k] = p_proj_coord[3*j+k];
      }
    }

    for (int j = 0; j < n_unlocated; j++) {
      int i = unlocated[j] - 1;
      tgt_location[i_part][i] = -1;
      for (int k = 0; k < 3; k++) {
        tgt_proj_coord[i_part][3*i+k] = tgt_coord[i_part][3*i+k];
      }
    }

  }


  /*
   * Create fields on src mesh
   */
  double **src_cell_volume = malloc( n_part * sizeof(double *));
  double **src_cell_center = malloc( n_part * sizeof(double *));
  double **src_field       = malloc( n_part * sizeof(double *));
  double **src_grad_field  = malloc( n_part * sizeof(double *));
  double **src_blk_interp  = malloc( n_part * sizeof(double *));
  for (int i_part = 0; i_part < n_part; i_part++) {
    const int is_oriented = 1;
    src_cell_volume[i_part] = (double *) malloc(sizeof(double) *     src_n_cell[i_part]);
    src_cell_center[i_part] = (double *) malloc(sizeof(double) * 3 * src_n_cell[i_part]);
    src_field      [i_part] = (double *) malloc(sizeof(double) *     src_n_cell[i_part]);
    src_grad_field [i_part] = (double *) malloc(sizeof(double) * 3 * src_n_cell[i_part]);
    src_blk_interp [i_part] = (double *) malloc(sizeof(double) *     src_n_cell[i_part]);

    PDM_geom_elem_polyhedra_properties (is_oriented,
                                        src_n_cell[i_part],
                                        src_n_face[i_part],
                                        src_face_vtx_idx[i_part],
                                        src_face_vtx[i_part],
                                        src_cell_face_idx[i_part],
                                        src_cell_face[i_part],
                                        src_n_vtx[i_part],
                                        src_vtx_coord[i_part],
                                        src_cell_volume[i_part],
                                        src_cell_center[i_part],
                                        NULL,
                                        NULL);


    // int n_cell_without_ext = src_n_cell[i_part]-src_n_cell_ext[i_part];

    for(int i = 0; i < src_n_cell[i_part]; ++i) {
      src_blk_interp[i_part][i] = -1;
      src_field     [i_part][i] = -100000;
      src_grad_field[i_part][3*i  ] = -100000;
      src_grad_field[i_part][3*i+1] = -100000;
      src_grad_field[i_part][3*i+2] = -100000;
    }

    for(int i = 0; i < src_n_cell[i_part]; ++i) {
      double xc = src_cell_center[i_part][3*i  ];
      double yc = src_cell_center[i_part][3*i+1];
      double zc = src_cell_center[i_part][3*i+2];

      // src_field     [i_part][i    ] = sin(xc + yc + zc);
      // src_grad_field[i_part][3*i  ] = cos(xc + yc + zc);
      // src_grad_field[i_part][3*i+1] = cos(xc + yc + zc);
      // src_grad_field[i_part][3*i+2] = cos(xc + yc + zc);

      src_field     [i_part][i    ] = sin(xc*xc + yc*yc + zc*zc);
      src_grad_field[i_part][3*i  ] = 2. * xc * cos(xc + yc + zc);
      src_grad_field[i_part][3*i+1] = 2. * yc * cos(xc + yc + zc);
      src_grad_field[i_part][3*i+2] = 2. * zc * cos(xc + yc + zc);

    }
  }


  /*
   *  Interpolate in src mesh
   */
  double      **src_tgt_field = malloc(n_part * sizeof(double      *));
  int         **src_tgt_idx   = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **src_tgt_gnum  = malloc(n_part * sizeof(PDM_g_num_t *));
  for (int i_part = 0; i_part < n_part; i_part++) {
    int *cell_vtx_idx;
    int *cell_vtx;
    PDM_mesh_location_cell_vertex_get (mesh_loc,
                                       i_part,
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
                                         i_part,
                                         0,//i_point_cloud,
                                         &elt_pts_inside_idx,
                                         &points_gnum,
                                         &points_coords,
                                         &points_uvw,
                                         &points_weights_idx,
                                         &points_weights,
                                         &points_dist2,
                                         &points_projected_coords);

    src_tgt_idx [i_part] = elt_pts_inside_idx;
    src_tgt_gnum[i_part] = points_gnum;

    // PDM_log_trace_connectivity_int(elt_pts_inside_idx, points_gnum, src_n_cell[i_part], "elt_pts_inside_idx ::");

    int n_cell_without_ext = src_n_cell[i_part]-src_n_cell_ext[i_part];
    int n_tgt_to_interp = elt_pts_inside_idx[n_cell_without_ext];

    src_tgt_field[i_part] = malloc(n_tgt_to_interp * sizeof(double));

    // int n_face_without_ext = src_n_face[i_part]-src_n_face_ext[i_part];
    // int n_vtx_without_ext  = src_n_vtx [i_part]-src_n_vtx_ext [i_part];
    for (int i = 0; i < n_cell_without_ext; i++) {
      for (int j = elt_pts_inside_idx[i]; j < elt_pts_inside_idx[i+1]; j++) {
        double tgt_x = points_coords[3*j  ];
        double tgt_y = points_coords[3*j+1];
        double tgt_z = points_coords[3*j+2];

        double xc = src_cell_center[i_part][3*i  ];
        double yc = src_cell_center[i_part][3*i+1];
        double zc = src_cell_center[i_part][3*i+2];

        double val   = src_field     [i_part][i    ];
        double gvalx = src_grad_field[i_part][3*i  ];
        double gvaly = src_grad_field[i_part][3*i+1];
        double gvalz = src_grad_field[i_part][3*i+2];

        // Do interpolation
        // src_tgt_field[i_part][j] = (double) j; // val + (xc - tgt_x) * gvalx + (yc - tgt_y) * gvaly + (zc - tgt_z) * gvalz;
        src_tgt_field[i_part][j] = val + (tgt_x - xc) * gvalx + (tgt_y - yc) * gvaly + (tgt_z - zc) * gvalz;

      }
    }

    // PDM_log_trace_array_double(src_tgt_field[i_part], n_tgt_to_interp, "src_tgt_field ::");
  }

  /*
   * Send to target
   */
  int* src_cell_without_ext = malloc(n_part * sizeof(int));
  int* tgt_cell_without_ext = malloc(n_part * sizeof(int));
  int* src_face_without_ext = malloc(n_part * sizeof(int));
  int* tgt_face_without_ext = malloc(n_part * sizeof(int));
  int* src_vtx_without_ext  = malloc(n_part * sizeof(int));
  int* tgt_vtx_without_ext  = malloc(n_part * sizeof(int));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    src_cell_without_ext[i_part] = src_n_cell[i_part] - src_n_cell_ext[i_part];
    tgt_cell_without_ext[i_part] = tgt_n_cell[i_part] - tgt_n_cell_ext[i_part];
    src_face_without_ext[i_part] = src_n_face[i_part] - src_n_face_ext[i_part];
    tgt_face_without_ext[i_part] = tgt_n_face[i_part] - tgt_n_face_ext[i_part];
    src_vtx_without_ext [i_part] = src_n_vtx [i_part] - src_n_vtx_ext [i_part];
    tgt_vtx_without_ext [i_part] = tgt_n_vtx [i_part] - tgt_n_vtx_ext [i_part];
  }

  PDM_part_to_part_t *ptp = PDM_part_to_part_create ((const PDM_g_num_t **) src_cell_ln_to_gn,
                                                                            src_cell_without_ext,
                                                                            n_part,
                                                     (const PDM_g_num_t **) tgt_cell_ln_to_gn,
                                                                            tgt_cell_without_ext,
                                                                            n_part,
                                                     (const int         **) src_tgt_idx,
                                                     (const PDM_g_num_t **) src_tgt_gnum,
                                                                            comm);

  /*
   * On peut donc echanger des src -> tgt
   */
  double **part_data = NULL;
  int request_exch;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
         (const int  **) NULL,
         (const void **) src_tgt_field,
                         NULL,
              (void ***) &part_data,
                         &request_exch);

  PDM_part_to_part_iexch_wait(ptp, request_exch);

  int  *n_ref_tgt;
  int **ref_tgt;
  PDM_part_to_part_ref_lnum2_get (ptp,
                                  &n_ref_tgt,
                                  &ref_tgt);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);


  double **tgt_field       = malloc( n_part * sizeof(double *));
  double **tgt_blk_interp  = malloc( n_part * sizeof(double *));

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_cell_without_ext = tgt_n_cell[i_part]-tgt_n_cell_ext[i_part];

    tgt_field      [i_part] = (double *) malloc(sizeof(double) * n_cell_without_ext);
    tgt_blk_interp [i_part] = (double *) malloc(sizeof(double) * n_cell_without_ext);

    for(int i = 0; i < n_cell_without_ext; ++i) {
      tgt_field      [i_part][i] = -1.;
      tgt_blk_interp [i_part][i] = -1.;
    }

    /* Enforce with interpolation */
    for(int i = 0; i < n_ref_tgt[i_part]; ++i) {
      int i_tgt = ref_tgt[i_part][i]-1;
      for(int j = gnum1_come_from_idx[i_part][i]; j < gnum1_come_from_idx[i_part][i+1]; ++j) {
        tgt_field     [i_part][i_tgt] = part_data[i_part][j];
        tgt_blk_interp[i_part][i_tgt] = 0;
      }
    }

    // PDM_log_trace_connectivity_int(gnum1_come_from_idx[i_part], ref_tgt[i_part], n_ref_tgt[i_part], "ref_tgt :: ");
    // PDM_log_trace_array_double(part_data[i_part], gnum1_come_from_idx[i_part][n_ref_tgt[i_part]], "recv_field :: ");

  }


  PDM_part_to_part_free(ptp);

  _visu ("source_mesh",
         n_part,
         src_cell_without_ext,
         src_face_without_ext,
         src_vtx_without_ext,
         src_cell_face_idx,
         src_cell_face,
         src_cell_ln_to_gn,
         src_face_vtx_idx,
         src_face_vtx,
         src_face_ln_to_gn,
         src_vtx_coord,
         src_vtx_ln_to_gn,
         src_field,
         src_blk_interp);

  _visu ("target_mesh",
         n_part,
         tgt_cell_without_ext,
         tgt_face_without_ext,
         tgt_vtx_without_ext,
         tgt_cell_face_idx,
         tgt_cell_face,
         tgt_cell_ln_to_gn,
         tgt_face_vtx_idx,
         tgt_face_vtx,
         tgt_face_ln_to_gn,
         tgt_vtx_coord,
         tgt_vtx_ln_to_gn,
         tgt_field,
         tgt_blk_interp);



  free(src_cell_without_ext);
  free(tgt_cell_without_ext);
  free(src_face_without_ext);
  free(tgt_face_without_ext);
  free(src_vtx_without_ext );
  free(tgt_vtx_without_ext );
  free(src_tgt_idx);
  free(src_tgt_gnum);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(src_tgt_field [i_part]);
    free(src_blk_interp[i_part]);
  }
  free(src_tgt_field);
  free(src_blk_interp);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(tgt_field     [i_part]);
    free(tgt_blk_interp[i_part]);
    free(part_data     [i_part]);
  }
  free(tgt_field);
  free(tgt_blk_interp);
  free(part_data);

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
  for (int i_part = 0; i_part < n_part; i_part++) {
    free(src_cell_volume[i_part]);
    free(src_cell_center[i_part]);
    free(src_field      [i_part]);
    free(src_grad_field [i_part]);
  }

  free(src_cell_volume);
  free(src_cell_center);
  free(src_field);
  free(src_grad_field);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(src_cell_face_idx[i_part]);
    free(src_cell_face    [i_part]);
    free(src_face_vtx_idx [i_part]);
    free(src_face_vtx     [i_part]);
    free(src_vtx_coord    [i_part]);
    free(src_cell_ln_to_gn[i_part]);
    free(src_face_ln_to_gn[i_part]);
    free(src_vtx_ln_to_gn [i_part]);

    free(src_face_group_idx[i_part]);
    free(src_face_group    [i_part]);
    free(src_group_ln_to_gn[i_part]);

    free(tgt_cell_face_idx [i_part]);
    free(tgt_cell_face     [i_part]);
    free(tgt_face_vtx_idx  [i_part]);
    free(tgt_face_vtx      [i_part]);
    free(tgt_vtx_coord     [i_part]);
    free(tgt_cell_ln_to_gn [i_part]);
    free(tgt_face_ln_to_gn [i_part]);
    free(tgt_vtx_ln_to_gn  [i_part]);
    free(tgt_face_group_idx[i_part]);
    free(tgt_face_group    [i_part]);
    free(tgt_group_ln_to_gn[i_part]);

    free (tgt_location[i_part]);
    free (tgt_proj_coord[i_part]);

    free(tgt_g_num[i_part]);

    if (!use_tgt_nodes) {
      free(tgt_coord[i_part]);
    }
  }
  free(tgt_g_num);
  free(n_tgt);
  free(src_n_cell       );
  free(src_n_face       );
  free(src_n_vtx        );
  free(src_n_cell_ext   );
  free(src_n_face_ext   );
  free(src_n_vtx_ext    );
  free(src_cell_face_idx);
  free(src_cell_face    );
  free(src_face_vtx_idx );
  free(src_face_vtx     );
  free(src_vtx_coord    );
  free(src_cell_ln_to_gn);
  free(src_face_ln_to_gn);
  free(src_vtx_ln_to_gn );

  free(src_face_group_idx);
  free(src_face_group    );
  free(src_group_ln_to_gn);

  free(tgt_n_cell       );
  free(tgt_n_face       );
  free(tgt_n_vtx        );
  free(tgt_n_cell_ext   );
  free(tgt_n_face_ext   );
  free(tgt_n_vtx_ext    );
  free(tgt_cell_face_idx);
  free(tgt_cell_face    );
  free(tgt_face_vtx_idx );
  free(tgt_face_vtx     );
  free(tgt_vtx_coord    );
  free(tgt_cell_ln_to_gn);
  free(tgt_face_ln_to_gn);
  free(tgt_vtx_ln_to_gn );

  free(tgt_face_group_idx);
  free(tgt_face_group    );
  free(tgt_group_ln_to_gn);

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

  return 0;
}

