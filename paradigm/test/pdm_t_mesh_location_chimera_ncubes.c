#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

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
#include "pdm_inside_cloud_surf.h"
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
     "  -doctree         Use doctree-based method.\n\n"
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
      *part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
    }
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = PDM_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DBBTREE;
    }
    else if (strcmp(argv[i], "-doctree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DOCTREE;
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

static void _rotate_by_angle (const int  n_pts,
                              double    *coord,
                              double     angle)
{

  // for (int i = 0; i < n_pts; i++) {
  //   double x = coord[3*i];
  //   double y = coord[3*i+1];
  //   double z = coord[3*i+2];
  //   double theta = z;
  //   coord[3*i+1] = y * cos(theta);
  //   coord[3*i+2] = z * sin(theta);
  // }
  double Rx[3][3] = {{1.,          0,           0},
                     {0., cos(angle), -sin(angle)},
                     {0., sin(angle),  cos(angle)}};

  for (int i = 0; i < n_pts; i++) {
    double x = coord[3*i];
    double y = coord[3*i+1];
    double z = coord[3*i+2];

    for (int j = 0; j < 3; j++) {
      coord[3*i+j] = Rx[j][0]*x + Rx[j][1]*y + Rx[j][2]*z;
    }
  }
}

static void _rotate_rand (const int  n_pts,
                          double    *coord)
{
  double i_rand_max = 1. / ((double) RAND_MAX);
  double pi = 4 * atan(1.);
  double theta[3] = {(double) rand() * i_rand_max * pi ,
                     (double) rand() * i_rand_max * pi ,
                     (double) rand() * i_rand_max * pi };

  double Rx[3][3] = {{1.,             0,             0},
                     {0., cos(theta[0]), -sin(theta[0])},
                     {0., sin(theta[0]),  cos(theta[0])}};

  double Ry[3][3] = {{  cos(theta[1]) , 0,  sin(theta[1])},
                     {  0             , 1,  0            },
                     {-sin(theta[1])  , 0,  cos(theta[1])}};

  double Rz[3][3] = {{cos(theta[2]), -sin(theta[2]), 0},
                     {sin(theta[2]),  cos(theta[2]), 0},
                     {0            ,  0            , 1}};

  for (int i = 0; i < n_pts; i++) {
    double x = coord[3*i  ];
    double y = coord[3*i+1];
    double z = coord[3*i+2];

    for (int j = 0; j < 3; j++) {
      coord[3*i+j] = Rx[j][0]*x + Rx[j][1]*y + Rx[j][2]*z;
    }
    x = coord[3*i  ];
    y = coord[3*i+1];
    z = coord[3*i+2];

    for (int j = 0; j < 3; j++) {
      coord[3*i+j] = Ry[j][0]*x + Ry[j][1]*y + Ry[j][2]*z;
    }
    x = coord[3*i  ];
    y = coord[3*i+1];
    z = coord[3*i+2];
    for (int j = 0; j < 3; j++) {
      coord[3*i+j] = Rz[j][0]*x + Rz[j][1]*y + Rz[j][2]*z;
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

  if(i_rank == 0)  {
    printf("n_vtx_seg = "PDM_FMT_G_NUM" \n", n_vtx_seg);
    printf("elt_type  = %i \n", elt_type);
  }
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



  if (deform == 1) {
    /*for (int i = 0; i < dn_vtx; i++) {
      double x = dvtx_coord[3*i];
      double z = dvtx_coord[3*i + 2];

      dvtx_coord[3*i]     += 0.1 * z * z;
      dvtx_coord[3*i + 2] += 0.2 * cos(PDM_PI * x);
      }*/
    if(1 == 0) {
      _rotate (dn_vtx, dvtx_coord);
    }
    _rotate_rand(dn_vtx, dvtx_coord);
  } else if( deform >= 3) {
    double angle = (deform-2) * PDM_PI/6;
    // printf("angle = %12.5e \n", angle);
    _rotate_by_angle (dn_vtx, dvtx_coord, angle);
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


// static void _export_point_cloud
// (
//  char         *filename,
//  int           n_part,
//  int          *n_pts,
//  double      **coord,
//  PDM_g_num_t **g_num,
//  PDM_g_num_t **location
//  )
// {
//   FILE *f = fopen(filename, "w");

//   fprintf(f, "# vtk DataFile Version 2.0\npoints\nASCII\nDATASET UNSTRUCTURED_GRID\n");

//   int n_pts_t = 0;
//   for (int i_part = 0; i_part < n_part; i_part++) {
//     n_pts_t += n_pts[i_part];
//   }

//   fprintf(f, "POINTS %d double\n", n_pts_t);
//   for (int i_part = 0; i_part < n_part; i_part++) {
//     for (int i = 0; i < n_pts[i_part]; i++) {
//       for (int j = 0; j < 3; j++) {
//         fprintf(f, "%f ", coord[i_part][3*i + j]);
//       }
//       fprintf(f, "\n");
//     }
//   }

//   fprintf(f, "CELLS %d %d\n", n_pts_t, 2*n_pts_t);
//   for (int i = 0; i < n_pts_t; i++) {
//     fprintf(f, "1 %d\n", i);
//   }

//   fprintf(f, "CELL_TYPES %d\n", n_pts_t);
//   for (int i = 0; i < n_pts_t; i++) {
//     fprintf(f, "1\n");
//   }

//   fprintf(f, "CELL_DATA %d\n", n_pts_t);
//   fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
//   for (int i_part = 0; i_part < n_part; i_part++) {
//     for (int i = 0; i < n_pts[i_part]; i++) {
//       fprintf(f, ""PDM_FMT_G_NUM"\n", g_num[i_part][i]);
//     }
//   }

//   if (location != NULL) {
//     fprintf(f, "FIELD FieldData 1\n");
//     fprintf(f, "location 1 %d int\n", n_pts_t);
//     for (int i_part = 0; i_part < n_part; i_part++) {
//       for (int i = 0; i < n_pts[i_part]; i++) {
//         fprintf(f, ""PDM_FMT_G_NUM"\n", location[i_part][i]);
//       }
//     }
//   }

//   fclose(f);
// }


/**
 *   Compute volume and center_cell
 */
static
void
_compute_centers
(
  int           n_part,
  int          *n_cell,
  int          *n_face,
  int          *n_vtx,
  int         **cell_face_idx,
  int         **cell_face,
  int         **face_vtx_idx,
  int         **face_vtx,
  double      **vtx_coord,
  double     ***cell_center_out,
  double     ***cell_volume_out,
  double     ***face_center_out
)
{

  double **cell_volume = (double ** ) malloc( n_part * sizeof(double));
  double **cell_center = (double ** ) malloc( n_part * sizeof(double));
  double **face_center = (double ** ) malloc( n_part * sizeof(double));
  for (int i_part = 0; i_part < n_part; i_part++) {
    const int is_oriented = 1;
    cell_volume[i_part] = (double *) malloc(sizeof(double) *     n_cell[i_part]);
    cell_center[i_part] = (double *) malloc(sizeof(double) * 3 * n_cell[i_part]);

    face_center[i_part] = (double *) malloc(sizeof(double) * 3 * n_face[i_part]);

    double* face_normal = (double *) malloc(sizeof(double) * 3 * n_face[i_part]);

    PDM_geom_elem_polyhedra_properties(is_oriented,
                                       n_cell       [i_part],
                                       n_face       [i_part],
                                       face_vtx_idx [i_part],
                                       face_vtx     [i_part],
                                       cell_face_idx[i_part],
                                       cell_face    [i_part],
                                       n_vtx        [i_part],
                                       vtx_coord    [i_part],
                                       cell_volume  [i_part],
                                       cell_center  [i_part],
                                       NULL,
                                       NULL);

    PDM_geom_elem_polygon_properties(n_face      [i_part],
                                     face_vtx_idx[i_part],
                                     face_vtx    [i_part],
                                     vtx_coord   [i_part],
                                     face_normal,
                                     face_center [i_part],
                                     NULL,
                                     NULL);
    free(face_normal);

  }

  *cell_center_out = cell_center;
  *cell_volume_out = cell_volume;
  *face_center_out = face_center;

}


/**
 *   Compute volume and center_cell
 */
static
void
_init_fields
(
  int           n_part,
  int          *n_cell,
  double      **cell_center,
  double     ***field,
  double     ***grad_field,
  double     ***blk_interp_from,
  double     ***blk_interp_vol,
  double     ***cell_nat
)
{
  *field           = (double ** ) malloc( n_part * sizeof(double *));
  *grad_field      = (double ** ) malloc( n_part * sizeof(double *));
  *blk_interp_from = (double ** ) malloc( n_part * sizeof(double *));
  *blk_interp_vol  = (double ** ) malloc( n_part * sizeof(double *));
  *cell_nat        = (double ** ) malloc( n_part * sizeof(double *));

  double **_field           = *field;
  double **_grad_field      = *grad_field;
  double **_blk_interp_from = *blk_interp_from;
  double **_blk_interp_vol  = *blk_interp_vol;
  double **_cell_nat        = *cell_nat;

  for (int i_part = 0; i_part < n_part; i_part++) {
    _field          [i_part] = (double *) malloc(    n_cell[i_part] * sizeof(double));
    _grad_field     [i_part] = (double *) malloc(3 * n_cell[i_part] * sizeof(double));
    _blk_interp_from[i_part] = (double *) malloc(    n_cell[i_part] * sizeof(double));
    _blk_interp_vol [i_part] = (double *) malloc(    n_cell[i_part] * sizeof(double));
    _cell_nat       [i_part] = (double *) malloc(    n_cell[i_part] * sizeof(double));

    for(int i = 0; i < n_cell[i_part]; ++i) {
      double xc = cell_center[i_part][3*i  ];
      double yc = cell_center[i_part][3*i+1];
      double zc = cell_center[i_part][3*i+2];

      // field     [i_part][i    ] = sin(xc + yc + zc);
      // grad_field[i_part][3*i  ] = cos(xc + yc + zc);
      // grad_field[i_part][3*i+1] = cos(xc + yc + zc);
      // grad_field[i_part][3*i+2] = cos(xc + yc + zc);

      _field          [i_part][i    ] = sin(xc*xc + yc*yc + zc*zc);
      _grad_field     [i_part][3*i  ] = 2. * xc * cos(xc + yc + zc);
      _grad_field     [i_part][3*i+1] = 2. * yc * cos(xc + yc + zc);
      _grad_field     [i_part][3*i+2] = 2. * zc * cos(xc + yc + zc);
      _blk_interp_from[i_part][i    ] = -1;
      _blk_interp_vol [i_part][i    ] = DBL_MAX ;
      _cell_nat       [i_part][i    ] = -4.; // Normal
    }
  }
}


/**
 *
 *
 *
 */
static
void
_prepare_target_cloud
(
  int            n_part,
  int            n_group,
  int           *n_cell,
  int           *n_face,
  int           *n_vtx,
  PDM_g_num_t  **cell_ln_to_gn,
  PDM_g_num_t  **face_ln_to_gn,
  PDM_g_num_t  **vtx_ln_to_gn,
  PDM_g_num_t  **face_group_ln_to_gn,
  int          **cell_face_idx,
  int          **cell_face,
  int          **face_vtx_idx,
  int          **face_vtx,
  int          **face_group_idx,
  int          **face_group,
  double       **vtx_coord,
  double       **cell_center,
  double       **face_center,
  double       **cell_nat,
  int            extract_center_depth,
  int            extract_bnd_faces,
  int          **n_extract_cell,
  int          **n_extract_face,
  PDM_g_num_t ***extract_center_ln_to_gn,
  double      ***extract_center_coord,
  PDM_g_num_t ***extract_face_bnd_ln_to_gn,
  double      ***extract_face_bnd_coord
)
{

  PDM_UNUSED(n_vtx);
  PDM_UNUSED(face_ln_to_gn);
  PDM_UNUSED(vtx_ln_to_gn);
  PDM_UNUSED(face_vtx_idx);
  PDM_UNUSED(face_vtx);
  PDM_UNUSED(vtx_coord);

  /*
   * Extract cell with a given depth and also face center
   */
  *n_extract_cell          = malloc(n_part * sizeof(int          ));
  *extract_center_coord    = malloc(n_part * sizeof(double      *));
  *extract_center_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));

  *n_extract_face            = malloc(n_part * sizeof(int          ));
  *extract_face_bnd_coord    = malloc(n_part * sizeof(double      *));
  *extract_face_bnd_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));

  int          *_n_extract_cell          = *n_extract_cell;
  double      **_extract_center_coord    = *extract_center_coord;
  PDM_g_num_t **_extract_center_ln_to_gn = *extract_center_ln_to_gn;

  int          *_n_extract_face            = *n_extract_face;
  double      **_extract_face_bnd_coord    = *extract_face_bnd_coord;
  PDM_g_num_t **_extract_face_bnd_ln_to_gn = *extract_face_bnd_ln_to_gn;

  /* Management of void */
  for (int i_part = 0; i_part < n_part; i_part++) {
    _n_extract_cell[i_part] = 0;
    _n_extract_face[i_part] = 0;
    _extract_center_coord     [i_part] = NULL;
    _extract_center_ln_to_gn  [i_part] = NULL;
    _extract_face_bnd_coord   [i_part] = NULL;
    _extract_face_bnd_ln_to_gn[i_part] = NULL;
  }

  if(extract_center_depth > 0) {
    for (int i_part = 0; i_part < n_part; i_part++) {

      int* face_cell_idx = NULL;
      int* face_cell     = NULL;
      PDM_connectivity_transpose(n_cell       [i_part],
                                 n_face       [i_part],
                                 cell_face_idx[i_part],
                                 cell_face    [i_part],
                                 &face_cell_idx,
                                 &face_cell);

      int* extract_cell = malloc(n_face[i_part] * sizeof(int)); // Surallocated buyt needed

      _n_extract_cell[i_part] = 0;

      for(int i_group = 0; i_group < n_group; ++i_group) {
        for(int idx_face = face_group_idx[i_part][i_group]; idx_face < face_group_idx[i_part][i_group+1]; ++idx_face) {
          int i_face = face_group[i_part][idx_face]-1;
          for(int idx_cell = face_cell_idx[i_face]; idx_cell < face_cell_idx[i_face+1]; ++idx_cell) {
            extract_cell[_n_extract_cell[i_part]++] = PDM_ABS(face_cell[idx_cell]);
          }
        }
      }

      // sort and extract coordinates
      _n_extract_cell[i_part] = PDM_inplace_unique(extract_cell, 0, _n_extract_cell[i_part]-1);

      extract_cell = realloc(extract_cell, _n_extract_cell[i_part] * sizeof(int));

      _extract_center_coord   [i_part] = (double      *) malloc(3 * _n_extract_cell[i_part] * sizeof(double     ));
      _extract_center_ln_to_gn[i_part] = (PDM_g_num_t *) malloc(    _n_extract_cell[i_part] * sizeof(PDM_g_num_t));

      for(int i = 0; i < _n_extract_cell[i_part]; ++i) {
        int i_cell = extract_cell[i]-1;
        _extract_center_coord   [i_part][3*i  ] = cell_center[i_part][3*i_cell  ];
        _extract_center_coord   [i_part][3*i+1] = cell_center[i_part][3*i_cell+1];
        _extract_center_coord   [i_part][3*i+2] = cell_center[i_part][3*i_cell+2];
        _extract_center_ln_to_gn[i_part][i] = cell_ln_to_gn[i_part][i_cell];
        cell_nat                [i_part][i_cell] = 1.; // Interpolated
      }

      free(face_cell_idx);
      free(face_cell);
      free(extract_cell);

    }
  }


  if(extract_bnd_faces == 1) {

    for (int i_part = 0; i_part < n_part; i_part++) {
      _n_extract_face[i_part] = 0;

      _extract_face_bnd_coord   [i_part] = malloc(face_group_idx[i_part][n_group] * sizeof(double     ));
      _extract_face_bnd_ln_to_gn[i_part] = malloc(face_group_idx[i_part][n_group] * sizeof(PDM_g_num_t));


      /* Dans cette exemple on considère que toutes les frontières sont BCOVerlap */
      for(int i_group = 0; i_group < n_group; ++i_group) {
        for(int idx_face = face_group_idx[i_part][i_group]; idx_face < face_group_idx[i_part][i_group+1]; ++idx_face) {
          int i_face = face_group[i_part][idx_face]-1;

          int idx_write = _n_extract_face[i_part]++;

          _extract_face_bnd_coord   [i_part][3*idx_write  ] = face_center[i_part][3*i_face  ];
          _extract_face_bnd_coord   [i_part][3*idx_write+1] = face_center[i_part][3*i_face+1];
          _extract_face_bnd_coord   [i_part][3*idx_write+2] = face_center[i_part][3*i_face+2];
          // _extract_face_bnd_ln_to_gn[i_part][idx_write] = face_ln_to_gn[i_part][i_face]; // Change rien
          _extract_face_bnd_ln_to_gn[i_part][idx_write] = face_group_ln_to_gn[i_part][idx_face];

        }
      }
      assert(_n_extract_face[i_part] == face_group_idx[i_part][n_group]);
    }
  }


}


/**
 *
 *
 *
 */
static
void
_prepare_external_faces
(
  PDM_MPI_Comm   comm,
  int            n_part,
  int            n_group,
  int           *n_face,
  int           *n_vtx,
  PDM_g_num_t  **face_ln_to_gn,
  PDM_g_num_t  **vtx_ln_to_gn,
  PDM_g_num_t  **face_group_ln_to_gn,
  int          **face_vtx_idx,
  int          **face_vtx,
  int          **face_group_idx,
  int          **face_group,
  double       **vtx_coord,
  int          **n_external_face,
  int          **n_external_vtx,
  PDM_g_num_t ***external_face_ln_to_gn,
  PDM_g_num_t ***external_vtx_ln_to_gn,
  double      ***external_vtx_coord,
  int         ***external_face_vtx_idx,
  int         ***external_face_vtx
)
{
  PDM_UNUSED(n_face);
  PDM_UNUSED(face_group_ln_to_gn);

  *n_external_face        = malloc(n_part * sizeof(int          ));
  *n_external_vtx         = malloc(n_part * sizeof(int          ));
  *external_face_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));
  *external_vtx_ln_to_gn  = malloc(n_part * sizeof(PDM_g_num_t *));

  *external_vtx_coord    = malloc(n_part * sizeof(double      *));
  *external_face_vtx_idx = malloc(n_part * sizeof(int         *));
  *external_face_vtx     = malloc(n_part * sizeof(int         *));

  int          *_n_external_face        = *n_external_face;
  int          *_n_external_vtx         = *n_external_vtx;
  PDM_g_num_t **_external_face_ln_to_gn = *external_face_ln_to_gn;
  PDM_g_num_t **_external_vtx_ln_to_gn  = *external_vtx_ln_to_gn;

  double      **_external_vtx_coord    = *external_vtx_coord;
  int         **_external_face_vtx_idx = *external_face_vtx_idx;
  int         **_external_face_vtx     = *external_face_vtx;

  /* Management of void */
  for (int i_part = 0; i_part < n_part; i_part++) {
    _n_external_face       [i_part] = 0;
    _n_external_vtx        [i_part] = 0;
    _external_face_ln_to_gn[i_part] = NULL;
    _external_vtx_ln_to_gn [i_part] = NULL;
    _external_vtx_coord    [i_part] = NULL;
    _external_face_vtx_idx [i_part] = NULL;
    _external_face_vtx     [i_part] = NULL;
  }

  for (int i_part = 0; i_part < n_part; i_part++) {

    int *vtx_tag = malloc(n_vtx[i_part] * sizeof(int));
    for(int i = 0; i < n_vtx[i_part]; ++i ) {
      vtx_tag[i] = -1;
    }

    _external_face_vtx_idx[i_part] = malloc( (face_group_idx[i_part][n_group]+1) * sizeof(int));
    _external_face_vtx_idx[i_part][0] = 0;

    /* Count */
    for(int i_group = 0; i_group < n_group; ++i_group) {
      for(int idx_face = face_group_idx[i_part][i_group]; idx_face < face_group_idx[i_part][i_group+1]; ++idx_face) {
        int i_face = face_group[i_part][idx_face]-1;

        int n_vtx_on_face = face_vtx_idx[i_part][i_face+1] - face_vtx_idx[i_part][i_face];
        _external_face_vtx_idx[i_part][_n_external_face[i_part]+1]  = _external_face_vtx_idx[i_part][_n_external_face[i_part]] + n_vtx_on_face;
        _n_external_face[i_part]++;
      }
    }

    _external_face_vtx[i_part] = malloc( (_external_face_vtx_idx[i_part][_n_external_face[i_part]]) * sizeof(int));

    _external_face_ln_to_gn[i_part] = malloc(_n_external_face[i_part] * sizeof(PDM_g_num_t));;

    int max_size = _external_face_vtx_idx[i_part][_n_external_face[i_part]];
    _external_vtx_ln_to_gn [i_part] = malloc(    max_size * sizeof(PDM_g_num_t));;
    _external_vtx_coord    [i_part] = malloc(3 * max_size * sizeof(double     ));;

    _n_external_face[i_part] = 0;
    _n_external_vtx [i_part] = 0;

    /*
     * Fill
     */
    int idx_write = 0;
    for(int i_group = 0; i_group < n_group; ++i_group) {
      for(int idx_face = face_group_idx[i_part][i_group]; idx_face < face_group_idx[i_part][i_group+1]; ++idx_face) {
        int i_face = face_group[i_part][idx_face]-1;

        for(int idx_vtx = face_vtx_idx[i_part][i_face]; idx_vtx < face_vtx_idx[i_part][i_face+1]; ++idx_vtx) {
          int i_vtx = face_vtx[i_part][idx_vtx]-1;

          if(vtx_tag[i_vtx] == -1) {
            vtx_tag[i_vtx] = _n_external_vtx [i_part]+1;
            _external_vtx_ln_to_gn[i_part][_n_external_vtx [i_part]] = vtx_ln_to_gn[i_part][i_vtx];
            _external_face_vtx[i_part][idx_write++] = vtx_tag[i_vtx];

            _external_vtx_coord[i_part][3*_n_external_vtx[i_part]  ] = vtx_coord[i_part][3*i_vtx  ];
            _external_vtx_coord[i_part][3*_n_external_vtx[i_part]+1] = vtx_coord[i_part][3*i_vtx+1];
            _external_vtx_coord[i_part][3*_n_external_vtx[i_part]+2] = vtx_coord[i_part][3*i_vtx+2];

            _n_external_vtx [i_part]++;
          } else {
            _external_face_vtx[i_part][idx_write++] = vtx_tag[i_vtx];
          }

        }

        _external_face_ln_to_gn[i_part][_n_external_face[i_part]] = face_ln_to_gn[i_part][i_face];
        _n_external_face[i_part]++;
      }
    }

    /* Realloc */
    _external_vtx_ln_to_gn [i_part] = realloc(_external_vtx_ln_to_gn [i_part],     _n_external_vtx[i_part] * sizeof(PDM_g_num_t));
    _external_vtx_coord    [i_part] = realloc(_external_vtx_coord    [i_part], 3 * _n_external_vtx[i_part] * sizeof(double     ));

    free(vtx_tag);
  }


  /*
   * Generate new global numbering
   */
  PDM_gen_gnum_t* gen_gnum_face = PDM_gnum_create(3, n_part, PDM_TRUE, 1e-6, comm, PDM_OWNERSHIP_USER);
  PDM_gen_gnum_t* gen_gnum_vtx  = PDM_gnum_create(3, n_part, PDM_TRUE, 1e-6, comm, PDM_OWNERSHIP_USER);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_gnum_set_from_parents(gen_gnum_face, i_part, _n_external_face[i_part], _external_face_ln_to_gn[i_part]);
    PDM_gnum_set_from_parents(gen_gnum_vtx , i_part, _n_external_vtx [i_part], _external_vtx_ln_to_gn [i_part]);
  }
  PDM_gnum_compute(gen_gnum_face);
  PDM_gnum_compute(gen_gnum_vtx );

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(_external_face_ln_to_gn[i_part]);
    free(_external_vtx_ln_to_gn [i_part]);
    _external_face_ln_to_gn[i_part] = PDM_gnum_get(gen_gnum_face, i_part);
    _external_vtx_ln_to_gn [i_part] = PDM_gnum_get(gen_gnum_vtx , i_part);
  }

  PDM_gnum_free(gen_gnum_face);
  PDM_gnum_free(gen_gnum_vtx);

}



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
 double          *block_interp[],
 double          *cell_nat[],
 double          *mask[]
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
  int id_var_cell_nat = PDM_writer_var_create(id_cs,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "cell_nat");
  int id_var_mask = PDM_writer_var_create(id_cs,
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          "mask");

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
    PDM_writer_var_set(id_cs,
                       id_var_cell_nat,
                       id_geom,
                       i_part,
                       cell_nat[i_part]);
    PDM_writer_var_set(id_cs,
                       id_var_mask,
                       id_geom,
                       i_part,
                       mask[i_part]);
  }

  PDM_writer_var_write(id_cs, id_var_field);
  PDM_writer_var_write(id_cs, id_var_num_part);
  PDM_writer_var_write(id_cs, id_var_blk_interp);
  PDM_writer_var_write(id_cs, id_var_cell_nat);
  PDM_writer_var_write(id_cs, id_var_mask);

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
// mpirun -np 24 ./test/pdm_t_mesh_location_chimera_ncubes -post -sepx 0.33 -sepy 0.33 -sepz 0.33 -def -n 100 -ext_depth 1
// case_type = 1 --> mpirun -np 24 ./test/pdm_t_mesh_location_chimera_ncubes -post -sepx 0.33 -sepy 0.33 -sepz 0.33 -def -l 2 -n 100
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
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#else
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;
#endif
#endif

  PDM_g_num_t n_pts = 10;
  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;
  // PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_DBBTREE;
  int disable_uvw = 0;
  int use_tgt_nodes = 0;

  // PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_HEXA8;
  PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_TETRA4;
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


  // printf("separation_x : %12.5e\n", separation_x);

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


  // const double xmin = -length/2;
  // const double ymin = -length/2;
  // const double zmin = -length/2;

  /*
   *  Source cube
   */
  // int case_type = 0; // Random cube mesh
  int case_type = 1; // Helice configuration
  int n_mesh = -1;
  // int n_mesh = 2;
  // int n_mesh = 5;
  if(case_type == 1) {
    n_mesh = 1 + 3*4;
  } else {
    n_mesh = 5;
  }

  int          **n_cell         = malloc(n_mesh * sizeof(int          *));
  int          **n_face         = malloc(n_mesh * sizeof(int          *));
  int          **n_vtx          = malloc(n_mesh * sizeof(int          *));

  int          **n_cell_without_ext = malloc(n_mesh * sizeof(int          *));
  int          **n_face_without_ext = malloc(n_mesh * sizeof(int          *));
  int          **n_vtx_without_ext  = malloc(n_mesh * sizeof(int          *));

  int          **n_cell_ext     = malloc(n_mesh * sizeof(int          *));
  int          **n_face_ext     = malloc(n_mesh * sizeof(int          *));
  int          **n_vtx_ext      = malloc(n_mesh * sizeof(int          *));
  int         ***cell_face_idx  = malloc(n_mesh * sizeof(int         **));
  int         ***cell_face      = malloc(n_mesh * sizeof(int         **));
  int         ***face_vtx_idx   = malloc(n_mesh * sizeof(int         **));
  int         ***face_vtx       = malloc(n_mesh * sizeof(int         **));
  double      ***vtx_coord      = malloc(n_mesh * sizeof(double      **));
  PDM_g_num_t ***cell_ln_to_gn  = malloc(n_mesh * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***face_ln_to_gn  = malloc(n_mesh * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***vtx_ln_to_gn   = malloc(n_mesh * sizeof(PDM_g_num_t **));

  int         ***face_group_idx = malloc(n_mesh * sizeof(int         **));
  int         ***face_group     = malloc(n_mesh * sizeof(int         **));
  PDM_g_num_t ***group_ln_to_gn = malloc(n_mesh * sizeof(PDM_g_num_t **));

  int          **n_extract_cell            = malloc(n_mesh * sizeof(int          *));
  int          **n_extract_face            = malloc(n_mesh * sizeof(int          *));
  PDM_g_num_t ***extract_center_ln_to_gn   = malloc(n_mesh * sizeof(PDM_g_num_t **));
  double      ***extract_center_coord      = malloc(n_mesh * sizeof(double      **));
  PDM_g_num_t ***extract_face_bnd_ln_to_gn = malloc(n_mesh * sizeof(PDM_g_num_t **));
  double      ***extract_face_bnd_coord    = malloc(n_mesh * sizeof(double      **));

  int          **n_external_face         = malloc(n_mesh * sizeof(int          *));
  int          **n_external_vtx          = malloc(n_mesh * sizeof(int          *));
  PDM_g_num_t ***external_face_ln_to_gn  = malloc(n_mesh * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***external_vtx_ln_to_gn   = malloc(n_mesh * sizeof(PDM_g_num_t **));
  double      ***external_vtx_coord      = malloc(n_mesh * sizeof(double      **));
  int         ***external_face_vtx_idx   = malloc(n_mesh * sizeof(int         **));
  int         ***external_face_vtx       = malloc(n_mesh * sizeof(int         **));

  double ***cell_center = malloc(n_mesh * sizeof(double **));
  double ***cell_volume = malloc(n_mesh * sizeof(double **));
  double ***face_center = malloc(n_mesh * sizeof(double **));


  double ***field           = malloc(n_mesh * sizeof(double **));
  double ***grad_field      = malloc(n_mesh * sizeof(double **));
  double ***blk_interp_from = malloc(n_mesh * sizeof(double **));
  double ***blk_interp_vol  = malloc(n_mesh * sizeof(double **));
  double ***cell_nat        = malloc(n_mesh * sizeof(double **));
  double ***mask            = malloc(n_mesh * sizeof(double **));
  int    ***elmt_cross_surf = malloc(n_mesh * sizeof(int    **));

  const int n_timer = 8;
  double cpu_time_max[n_timer];
  double cpu_time_min[n_timer];
  double cpu_time_sum[n_timer];
  double delta_t     [n_timer];
  double t1, t2;

  PDM_MPI_Barrier(comm);

  t1 = PDM_MPI_Wtime();
  int n_group = 6;
  int ldeform = deform;
  int lextract_center_depth = 0;
  int lextract_bnd_faces    = 0;
  PDM_g_num_t n_cell_tot = 0;
  PDM_g_num_t n_interp_loc = 0;
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {

    if(i_mesh == 0){
      ldeform = 0;
    } else {
      ldeform = 1;
    }

    double i_rand_max = 1. / ((double) RAND_MAX);
    double lseparation_x = 0.;
    double lseparation_y = 0.;
    double lseparation_z = 0.;
    double llength       = length;
    double lxmin         = 0.;
    double lymin         = 0.;
    double lzmin         = 0.;

    if(case_type == 0) {
      srand(i_mesh);
      if(i_mesh > 0 ) {
        lseparation_x = rand() * i_rand_max * 0.33;
        lseparation_y = rand() * i_rand_max * 0.33;
        lseparation_z = rand() * i_rand_max * 0.33;

        lextract_center_depth = 1;
        llength = rand() * i_rand_max /3.;
        lxmin   = -llength/2;
        lymin   = -llength/2;
        lzmin   = -llength/2;
      } else {
        lxmin         = -llength/2;
        lymin         = -llength/2;
        lzmin         = -llength/2;
      }
    } else {
      srand(i_mesh);
      lxmin         = -llength/2;
      lymin         = -llength/2;
      lzmin         = -llength/2;

      if(i_mesh > 0 ) {
        lextract_center_depth = 1;
        llength = length/8;

        lxmin   = -llength/2;
        lymin   = +llength*2;
        lzmin   = -llength/2;

        ldeform = (i_mesh+1);
      }
    }

    // Rand a number around 20% for n_vtx_seg
    PDM_g_num_t percent = (PDM_g_num_t) ceil( (double) n_vtx_seg * (25. / 100.));
    PDM_g_num_t n_vtx_add_rand = rand() % ( percent ) - percent / 2;
    PDM_g_num_t n_vtx_seg_rand = n_vtx_seg + n_vtx_add_rand;

    if(i_rank == 0) {
      printf("lseparation = %12.5e / %12.5e / %12.5e \n", lseparation_x, lseparation_y, lseparation_z);
      printf("lxmin/lymin/lzmin = %12.5e / %12.5e / %12.5e \n", lxmin, lymin, lzmin);
      printf("llength     = %12.5e  \n", llength);
      printf("n_vtx_add_rand = "PDM_FMT_G_NUM" \n", n_vtx_add_rand);
      printf("n_vtx_seg_rand = "PDM_FMT_G_NUM" \n", n_vtx_seg_rand);
      printf("percent        = "PDM_FMT_G_NUM" \n", percent);
    }

    n_cell_tot += (n_vtx_seg_rand-1) * (n_vtx_seg_rand-1) * (n_vtx_seg_rand-1);


    /* Generate */
    _cube_mesh (comm,
                n_part,
                part_method,
                n_vtx_seg_rand,
                lxmin+lseparation_x*llength,
                lymin+lseparation_y*llength,
                lzmin+lseparation_z*llength,
                llength,
                ldeform,//deform
                extension_depth_src,
                elt_type,
                &n_cell        [i_mesh],
                &n_face        [i_mesh],
                &n_vtx         [i_mesh],
                &n_cell_ext    [i_mesh],
                &n_face_ext    [i_mesh],
                &n_vtx_ext     [i_mesh],
                &cell_face_idx [i_mesh],
                &cell_face     [i_mesh],
                &face_vtx_idx  [i_mesh],
                &face_vtx      [i_mesh],
                &vtx_coord     [i_mesh],
                &cell_ln_to_gn [i_mesh],
                &face_ln_to_gn [i_mesh],
                &vtx_ln_to_gn  [i_mesh],
                &face_group_idx[i_mesh],
                &face_group    [i_mesh],
                &group_ln_to_gn[i_mesh]);

    /* Compute geometry + extraction */
    _compute_centers(n_part,
                     n_cell[i_mesh],
                     n_face[i_mesh],
                     n_vtx[i_mesh],
                     cell_face_idx[i_mesh],
                     cell_face[i_mesh],
                     face_vtx_idx[i_mesh],
                     face_vtx[i_mesh],
                     vtx_coord[i_mesh],
                     &cell_center[i_mesh],
                     &cell_volume[i_mesh],
                     &face_center[i_mesh]);

    /* Init field */
    _init_fields(n_part,
                 n_cell          [i_mesh],
                 cell_center     [i_mesh],
                 &field          [i_mesh],
                 &grad_field     [i_mesh],
                 &blk_interp_from[i_mesh],
                 &blk_interp_vol [i_mesh],
                 &cell_nat       [i_mesh]);

    /* Prepare extract */
    _prepare_target_cloud(n_part,
                          n_group,
                          n_cell                   [i_mesh],
                          n_face                   [i_mesh],
                          n_vtx                    [i_mesh],
                          cell_ln_to_gn            [i_mesh],
                          face_ln_to_gn            [i_mesh],
                          vtx_ln_to_gn             [i_mesh],
                          group_ln_to_gn           [i_mesh],
                          cell_face_idx            [i_mesh],
                          cell_face                [i_mesh],
                          face_vtx_idx             [i_mesh],
                          face_vtx                 [i_mesh],
                          face_group_idx           [i_mesh],
                          face_group               [i_mesh],
                          vtx_coord                [i_mesh],
                          cell_center              [i_mesh],
                          face_center              [i_mesh],
                          cell_nat                 [i_mesh],
                          lextract_center_depth,
                          lextract_bnd_faces,
                          &n_extract_cell           [i_mesh],
                          &n_extract_face           [i_mesh],
                          &extract_center_ln_to_gn  [i_mesh],
                          &extract_center_coord     [i_mesh],
                          &extract_face_bnd_ln_to_gn[i_mesh],
                          &extract_face_bnd_coord   [i_mesh]);

    n_cell_without_ext[i_mesh] = malloc(n_part * sizeof(int));
    n_face_without_ext[i_mesh] = malloc(n_part * sizeof(int));
    n_vtx_without_ext [i_mesh] = malloc(n_part * sizeof(int));
    for(int i_part = 0; i_part < n_part; ++i_part) {
      n_cell_without_ext[i_mesh][i_part] = n_cell[i_mesh][i_part] - n_cell_ext[i_mesh][i_part];
      n_face_without_ext[i_mesh][i_part] = n_face[i_mesh][i_part] - n_face_ext[i_mesh][i_part];
      n_vtx_without_ext [i_mesh][i_part] = n_vtx [i_mesh][i_part] - n_vtx_ext [i_mesh][i_part];
      n_interp_loc                      += n_extract_cell[i_mesh][i_part];
    }

    /* Prepare extract */
    _prepare_external_faces(comm,
                            n_part,
                            n_group,
                            n_face                   [i_mesh],
                            n_vtx                    [i_mesh],
                            face_ln_to_gn            [i_mesh],
                            vtx_ln_to_gn             [i_mesh],
                            group_ln_to_gn           [i_mesh],
                            face_vtx_idx             [i_mesh],
                            face_vtx                 [i_mesh],
                            face_group_idx           [i_mesh],
                            face_group               [i_mesh],
                            vtx_coord                [i_mesh],
                            &n_external_face         [i_mesh],
                            &n_external_vtx          [i_mesh],
                            &external_face_ln_to_gn  [i_mesh],
                            &external_vtx_ln_to_gn   [i_mesh],
                            &external_vtx_coord      [i_mesh],
                            &external_face_vtx_idx   [i_mesh],
                            &external_face_vtx       [i_mesh]);

    /*
     * Export vtk extract surface
     */
    if(0 == 1) {
      for(int i_part = 0; i_part < n_part; ++i_part) {
        char filename[999];
        sprintf(filename, "external_faces_%i_%i_%i.vtk", i_mesh, i_part, i_rank);
        PDM_vtk_write_polydata(filename,
                               n_external_vtx          [i_mesh][i_part],
                               external_vtx_coord      [i_mesh][i_part],
                               external_vtx_ln_to_gn   [i_mesh][i_part],
                               n_external_face         [i_mesh][i_part],
                               external_face_vtx_idx   [i_mesh][i_part],
                               external_face_vtx       [i_mesh][i_part],
                               external_face_ln_to_gn  [i_mesh][i_part],
                               NULL);
      }
    }

    // /* Visu */
    // char filename[999];
    // sprintf(filename, "mesh_%i", i_mesh);
    // _visu (filename,
    //        n_part,
    //        n_cell_without_ext[i_mesh],
    //        n_face_without_ext[i_mesh],
    //        n_vtx_without_ext [i_mesh],
    //        cell_face_idx     [i_mesh],
    //        cell_face         [i_mesh],
    //        cell_ln_to_gn     [i_mesh],
    //        face_vtx_idx      [i_mesh],
    //        face_vtx          [i_mesh],
    //        face_ln_to_gn     [i_mesh],
    //        vtx_coord         [i_mesh],
    //        vtx_ln_to_gn      [i_mesh],
    //        field             [i_mesh],
    //        blk_interp_from   [i_mesh],
    //        cell_nat          [i_mesh]);
  }
  t2 = PDM_MPI_Wtime();
  delta_t[0] = t2 - t1;


  /*
   * First step : mask
   */
  for(int i_cloud = 0; i_cloud < n_mesh; ++i_cloud) {
    mask           [i_cloud] = malloc(n_part * sizeof(double *));
    elmt_cross_surf[i_cloud] = malloc(n_part * sizeof(int    *));
    for(int i_part = 0; i_part < n_part; ++i_part) {
      mask           [i_cloud][i_part] = malloc(n_cell[i_cloud][i_part] * sizeof(double));
      elmt_cross_surf[i_cloud][i_part] = malloc(n_cell[i_cloud][i_part] * sizeof(int));
      for(int i = 0; i < n_cell[i_cloud][i_part]; ++i) {
        mask           [i_cloud][i_part][i] = 0.; // 0 -> pas masqué
        elmt_cross_surf[i_cloud][i_part][i] = 0.; // Pas d'interaction avec une surface
      }
    }
  }



  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();

  PDM_inside_cloud_surf_t **ics = (PDM_inside_cloud_surf_t **) malloc( n_mesh * sizeof(PDM_inside_cloud_surf_t *));
  for(int i_mesh = 1; i_mesh < n_mesh; ++i_mesh) { // Because first mesh mask all we skip it !!!
    ics[i_mesh] = PDM_inside_cloud_surf_create(n_mesh, PDM_OWNERSHIP_USER, comm);

    for(int i_cloud = 0; i_cloud < n_mesh; ++i_cloud) {
      PDM_inside_cloud_surf_n_part_cloud_set(ics[i_mesh], i_cloud, n_part);
      if(i_mesh ==  i_cloud) {
        for(int i_part = 0; i_part < n_part; ++i_part) {
          PDM_inside_cloud_surf_cloud_set(ics[i_mesh], i_cloud, i_part, 0, NULL, NULL);
        }
      } else {
        for(int i_part = 0; i_part < n_part; ++i_part) {
          PDM_inside_cloud_surf_cloud_set(ics[i_mesh],
                                          i_cloud,
                                          i_part,
                                          n_vtx       [i_cloud][i_part],
                                          vtx_coord   [i_cloud][i_part],
                                          vtx_ln_to_gn[i_cloud][i_part]);
        }
      }
    }

    PDM_inside_cloud_surf_surf_mesh_global_data_set(ics[i_mesh], n_part);

    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_inside_cloud_surf_surf_mesh_part_set(ics[i_mesh],
                                               i_part,
                                               n_external_face       [i_mesh][i_part],
                                               external_face_vtx_idx [i_mesh][i_part],
                                               external_face_vtx     [i_mesh][i_part],
                                               external_face_ln_to_gn[i_mesh][i_part],
                                               n_external_vtx        [i_mesh][i_part],
                                               external_vtx_coord    [i_mesh][i_part],
                                               external_vtx_ln_to_gn [i_mesh][i_part]);
    }

    // PDM_inside_cloud_surf_compute(ics[i_mesh]);
    PDM_inside_cloud_surf_compute_optim(ics[i_mesh]);

    int mask_type = 1;
    for(int i_cloud = 0; i_cloud < n_mesh; ++i_cloud) {
      for(int i_part = 0; i_part < n_part; ++i_part) {
        int *is_inside = NULL;
        PDM_inside_cloud_surf_get(ics[i_mesh], i_cloud, i_part, &is_inside);
        if(i_mesh ==  i_cloud) {
          free(is_inside);
          continue;
        }

        /* Translate information into a single field at cell */
        for(int i_cell = 0; i_cell < n_cell[i_cloud][i_part]; ++i_cell) {

          int one_is_interior = 0;
          int one_is_exterior = 0;
          double cell_is_completely_inside = 1.;
          double face_is_inside = 0.;
          for(int idx_face = cell_face_idx[i_cloud][i_part][i_cell]; idx_face < cell_face_idx[i_cloud][i_part][i_cell+1]; ++idx_face) {
            int i_face = PDM_ABS(cell_face[i_cloud][i_part][idx_face]) - 1;
            for(int idx_vtx = face_vtx_idx[i_cloud][i_part][i_face]; idx_vtx < face_vtx_idx[i_cloud][i_part][i_face+1]; ++idx_vtx) {
              int i_vtx = face_vtx[i_cloud][i_part][idx_vtx]-1;
              // face_is_inside = PDM_MAX(face_is_inside, is_inside[i_vtx]);
              face_is_inside = PDM_MAX(face_is_inside, is_inside[i_vtx]);
              if(is_inside[i_vtx] == 0) {
                cell_is_completely_inside = 0;
                one_is_exterior = 1;
              } else {
                one_is_interior = 1;
              }
            }
          }
          // A voir : On prends le plus grand recouvrement ou le plus petit
          if(mask_type == 0) {
            mask[i_cloud][i_part][i_cell] = PDM_MAX(mask[i_cloud][i_part][i_cell],  face_is_inside);
          } else {
            mask[i_cloud][i_part][i_cell] = PDM_MAX(mask[i_cloud][i_part][i_cell], cell_is_completely_inside);
          }

          if(one_is_interior == 1 && one_is_exterior == 1) {
            elmt_cross_surf[i_cloud][i_part][i_cell] = 1;
          }

        }

        free(is_inside);
      }
    }


    PDM_inside_cloud_surf_free(ics[i_mesh]);
  }
  free(ics);
  t2 = PDM_MPI_Wtime();
  delta_t[1] = t2 - t1;

  /*
   * Prepare localisation
   */
  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();
  int n_point_cloud = 1;  // 1 : adjacent cell center | 2 : center_interface
  int n_tot_cloud = n_point_cloud * n_mesh;
  PDM_mesh_location_t **mesh_loc = (PDM_mesh_location_t **) malloc( n_mesh * sizeof(PDM_mesh_location_t *));

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {

    mesh_loc[i_mesh] = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED,
                                                 n_tot_cloud, // Number of cloud
                                                 comm,
                                                 PDM_OWNERSHIP_KEEP);

    PDM_mesh_location_reverse_results_enable(mesh_loc[i_mesh]);

    /* Set target point cloud */
    for(int i_cloud = 0; i_cloud < n_tot_cloud; ++i_cloud) {
      PDM_mesh_location_n_part_cloud_set (mesh_loc[i_mesh], i_cloud, n_part);
    }

    /* Set target point cloud */
    for(int i_cloud = 0; i_cloud < n_mesh; ++i_cloud) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        // printf("n_extract_cell         [i_cloud][i_part] = %i \n", n_extract_cell         [i_cloud][i_part]);
        if(i_mesh != i_cloud) {
          PDM_mesh_location_cloud_set(mesh_loc[i_mesh],
                                      i_cloud,
                                      i_part,
                                      n_extract_cell         [i_cloud][i_part],
                                      extract_center_coord   [i_cloud][i_part],
                                      extract_center_ln_to_gn[i_cloud][i_part]);
        } else {
          // printf("Not provided cloud \n");
          PDM_mesh_location_cloud_set(mesh_loc[i_mesh],
                                      i_cloud,
                                      i_part,
                                      0,
                                      NULL,
                                      NULL);
        }
      }
    }

    /* Set source mesh */
    PDM_mesh_location_mesh_global_data_set(mesh_loc[i_mesh], n_part);

    int **is_elmt_select_by_user = malloc(n_part * sizeof(int *));

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_mesh_location_part_set (mesh_loc[i_mesh],
                                  i_part,
                                  n_cell_without_ext[i_mesh][i_part],
                                  cell_face_idx     [i_mesh][i_part],
                                  cell_face         [i_mesh][i_part],
                                  cell_ln_to_gn     [i_mesh][i_part],
                                  n_face_without_ext[i_mesh][i_part],
                                  face_vtx_idx      [i_mesh][i_part],
                                  face_vtx          [i_mesh][i_part],
                                  face_ln_to_gn     [i_mesh][i_part],
                                  n_vtx_without_ext [i_mesh][i_part],
                                  vtx_coord         [i_mesh][i_part],
                                  vtx_ln_to_gn      [i_mesh][i_part]);

      is_elmt_select_by_user[i_part] = malloc(n_cell_without_ext[i_mesh][i_part] * sizeof(int));
      for(int i_cell = 0; i_cell < n_cell_without_ext[i_mesh][i_part]; ++i_cell) {
        if(mask[i_mesh][i_part][i_cell] > 0.5) {
          is_elmt_select_by_user[i_part][i_cell] = 0;
        } else {
          is_elmt_select_by_user[i_part][i_cell] = 1;
        }
      }

      // PDM_mesh_location_user_extract_set(mesh_loc[i_mesh],
      //                                    i_part,
      //                                    is_elmt_select_by_user[i_part]);
      PDM_mesh_location_user_extract_set(mesh_loc[i_mesh],
                                         i_part,
                                         elmt_cross_surf[i_mesh][i_part]);
    }

    /* Set location parameters */
    PDM_mesh_location_tolerance_set(mesh_loc[i_mesh], tolerance);
    PDM_mesh_location_method_set(mesh_loc[i_mesh], loc_method);

    char *env_var = NULL;
    env_var = getenv ("PDM_MESH_LOCATION_NEW");
    int algo = 1;
    if (env_var != NULL) {
      algo = atoi(env_var);
    }
    if(algo == 0) {
      PDM_mesh_location_compute (mesh_loc[i_mesh]);
    } else if(algo == 1) {
      // PDM_mesh_location_compute_optim2(mesh_loc[i_mesh]);
      PDM_mesh_location_compute_optim(mesh_loc[i_mesh]);
    }

    PDM_mesh_location_dump_times (mesh_loc[i_mesh]);

    for (int i_part = 0; i_part < n_part; i_part++) {
      free(is_elmt_select_by_user[i_part]);
    }
    free(is_elmt_select_by_user);

    PDM_MPI_Barrier(comm);
    // if (i_rank == 0) {
    //   printf("OK! :D");
    // }
    // PDM_MPI_Finalize();
    // return 0;


  }
  t2 = PDM_MPI_Wtime();
  delta_t[2] = t2 - t1;


  /*
   *  Prepare all part_to_part
   */
  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();
  int         ****mesh_pts_idx            = malloc(n_mesh * sizeof(int         ***));
  PDM_g_num_t ****mesh_pts_gnum           = malloc(n_mesh * sizeof(PDM_g_num_t ***));
  double      ****points_coords           = malloc(n_mesh * sizeof(double      ***));
  double      ****points_uvw              = malloc(n_mesh * sizeof(double      ***));
  int         ****points_weights_idx      = malloc(n_mesh * sizeof(int         ***));
  double      ****points_weights          = malloc(n_mesh * sizeof(double      ***));
  double      ****points_dist2            = malloc(n_mesh * sizeof(double      ***));
  double      ****points_projected_coords = malloc(n_mesh * sizeof(double      ***));

  PDM_part_to_part_t ***ptp = malloc(n_mesh * sizeof(PDM_part_to_part_t **));
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {

    mesh_pts_idx           [i_mesh] = malloc(n_tot_cloud * sizeof(int         ***));
    mesh_pts_gnum          [i_mesh] = malloc(n_tot_cloud * sizeof(PDM_g_num_t ***));
    points_coords          [i_mesh] = malloc(n_tot_cloud * sizeof(double      ***));
    points_uvw             [i_mesh] = malloc(n_tot_cloud * sizeof(double      ***));
    points_weights_idx     [i_mesh] = malloc(n_tot_cloud * sizeof(int         ***));
    points_weights         [i_mesh] = malloc(n_tot_cloud * sizeof(double      ***));
    points_dist2           [i_mesh] = malloc(n_tot_cloud * sizeof(double      ***));
    points_projected_coords[i_mesh] = malloc(n_tot_cloud * sizeof(double      ***));

    ptp[i_mesh] = malloc(n_tot_cloud * sizeof(PDM_part_to_part_t *));
    for(int i_cloud = 0; i_cloud < n_tot_cloud; ++i_cloud) {

      mesh_pts_idx           [i_mesh][i_cloud] = malloc(n_part * sizeof(int         *));
      mesh_pts_gnum          [i_mesh][i_cloud] = malloc(n_part * sizeof(PDM_g_num_t *));
      points_coords          [i_mesh][i_cloud] = malloc(n_part * sizeof(double      *));
      points_uvw             [i_mesh][i_cloud] = malloc(n_part * sizeof(double      *));
      points_weights_idx     [i_mesh][i_cloud] = malloc(n_part * sizeof(int         *));
      points_weights         [i_mesh][i_cloud] = malloc(n_part * sizeof(double      *));
      points_dist2           [i_mesh][i_cloud] = malloc(n_part * sizeof(double      *));
      points_projected_coords[i_mesh][i_cloud] = malloc(n_part * sizeof(double      *));

      // Partition de maillage !!!
      for (int i_part = 0; i_part < n_part; i_part++) {

        PDM_mesh_location_points_in_elt_get (mesh_loc[i_mesh],
                                             i_part,
                                             i_cloud,
                                             &mesh_pts_idx           [i_mesh][i_cloud][i_part],
                                             &mesh_pts_gnum          [i_mesh][i_cloud][i_part],
                                             &points_coords          [i_mesh][i_cloud][i_part],
                                             &points_uvw             [i_mesh][i_cloud][i_part],
                                             &points_weights_idx     [i_mesh][i_cloud][i_part],
                                             &points_weights         [i_mesh][i_cloud][i_part],
                                             &points_dist2           [i_mesh][i_cloud][i_part],
                                             &points_projected_coords[i_mesh][i_cloud][i_part]);

      }

      /*
       * Creation du part_to_part pour chaque mesh / cloud
       */
      assert(n_tot_cloud == n_mesh);
      // ptp[i_mesh][i_cloud] = PDM_part_to_part_create ((const PDM_g_num_t **) cell_ln_to_gn[i_mesh],
      //                                                                        n_cell_without_ext[i_mesh],
      //                                                                        n_part,
      //                                                 (const PDM_g_num_t **) cell_ln_to_gn[i_cloud],
      //                                                                        n_cell       [i_cloud],
      //                                                                        n_part,
      //                                                 (const int         **) mesh_pts_idx [i_mesh][i_cloud],
      //                                                 (const PDM_g_num_t **) mesh_pts_gnum[i_mesh][i_cloud],
      //                                                                        comm);
      PDM_mesh_location_part_to_part_get(mesh_loc[i_mesh],
                                         i_cloud,
                                         &ptp[i_mesh][i_cloud],
                                         PDM_OWNERSHIP_USER);

    }
  }
  t2 = PDM_MPI_Wtime();
  delta_t[3] = t2 - t1;


  /*
   *  Exchange data between src and target
   */
  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();
  int n_data_send = 5; // Value / blk_interp / volume / dist / cell_nat (Reminder if dist > 0 -> Extrapolate )
  double      ****send_buffer  = malloc(n_mesh * sizeof(double ***));
  double      ****recv_buffer  = malloc(n_mesh * sizeof(double ***));
  int           **request_exch = malloc(n_mesh * sizeof(int      *));
  // int           **recv_request = malloc(n_mesh * sizeof(int      *));
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {

    send_buffer [i_mesh] = malloc(n_tot_cloud * sizeof(double **));
    recv_buffer [i_mesh] = malloc(n_tot_cloud * sizeof(double **));
    request_exch[i_mesh] = malloc(n_tot_cloud * sizeof(int      ));
    // recv_request[i_mesh] = malloc(n_tot_cloud * sizeof(int      ));

    for(int i_cloud = 0; i_cloud < n_tot_cloud; ++i_cloud) {

      send_buffer [i_mesh][i_cloud] = malloc(n_part * sizeof(double *));

      for (int i_part = 0; i_part < n_part; i_part++) {

        int _n_cell_without_ext = n_cell_without_ext[i_mesh][i_part];

        // src_tgt_field[i_part] = malloc(n_tgt_to_interp * sizeof(double));

        int *_mesh_pts_idx = mesh_pts_idx [i_mesh][i_cloud][i_part];
        int n_tgt_to_interp = _mesh_pts_idx[_n_cell_without_ext];

        send_buffer[i_mesh][i_cloud][i_part] = malloc( n_data_send * n_tgt_to_interp * sizeof(double));

        double *_send_buffer = send_buffer [i_mesh][i_cloud][i_part];
        double *_pts_dist    = points_dist2[i_mesh][i_cloud][i_part];

        for (int i = 0; i < _n_cell_without_ext; i++) {
          for (int j = _mesh_pts_idx[i]; j < _mesh_pts_idx[i+1]; j++) {

            // Do interpolation
            double tgt_x = points_coords[i_mesh][i_cloud][i_part][3*j  ];
            double tgt_y = points_coords[i_mesh][i_cloud][i_part][3*j+1];
            double tgt_z = points_coords[i_mesh][i_cloud][i_part][3*j+2];

            double xc = cell_center[i_mesh][i_part][3*i  ];
            double yc = cell_center[i_mesh][i_part][3*i+1];
            double zc = cell_center[i_mesh][i_part][3*i+2];

            double val   = field     [i_mesh][i_part][i    ];
            double gvalx = grad_field[i_mesh][i_part][3*i  ];
            double gvaly = grad_field[i_mesh][i_part][3*i+1];
            double gvalz = grad_field[i_mesh][i_part][3*i+2];

            // Do interpolation
            // src_tgt_field[i_part][j] = (double) j; // val + (xc - tgt_x) * gvalx + (yc - tgt_y) * gvaly + (zc - tgt_z) * gvalz;
            // src_tgt_field[i_part][j] = val + (tgt_x - xc) * gvalx + (tgt_y - yc) * gvaly + (tgt_z - zc) * gvalz;

            _send_buffer[n_data_send*j  ] = val + (tgt_x - xc) * gvalx + (tgt_y - yc) * gvaly + (tgt_z - zc) * gvalz;;
            _send_buffer[n_data_send*j+1] = i_mesh;
            _send_buffer[n_data_send*j+2] = cell_volume[i_mesh][i_part][i];
            _send_buffer[n_data_send*j+3] = _pts_dist[j];
            _send_buffer[n_data_send*j+4] = cell_nat[i_mesh][i_part][i];

          }
        }
      }

      /*
       * Exchange - Can be done with explicit P2P
       */
       PDM_part_to_part_iexch(ptp[i_mesh][i_cloud],
                              PDM_MPI_COMM_KIND_P2P,
                              PDM_STRIDE_CST_INTERLACED,
                              PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                              n_data_send,
                              sizeof(double),
              (const int  **) NULL,
              (const void **) send_buffer[i_mesh][i_cloud],
                              NULL,
                   (void ***) &recv_buffer[i_mesh][i_cloud],
                              &request_exch[i_mesh][i_cloud]);

    }
  }
  t2 = PDM_MPI_Wtime();
  delta_t[4] = t2 - t1;


  /*
   * Receive data
   */
  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    for(int i_cloud = 0; i_cloud < n_tot_cloud; ++i_cloud) {
      // We can also do iexch_test and use the first request
      PDM_part_to_part_iexch_wait(ptp[i_mesh][i_cloud], request_exch[i_mesh][i_cloud]);

      /*
       * Post-treatment
       */
      int  *n_ref_tgt;
      int **ref_tgt;
      PDM_part_to_part_ref_lnum2_get (ptp[i_mesh][i_cloud], &n_ref_tgt, &ref_tgt);

      int         **gnum1_come_from_idx = NULL;
      PDM_g_num_t **gnum1_come_from     = NULL;
      PDM_part_to_part_gnum1_come_from_get(ptp[i_mesh][i_cloud], &gnum1_come_from_idx, &gnum1_come_from);

      for(int i_part = 0; i_part < n_part; ++i_part) {

        // int _n_cell_without_ext = n_cell_without_ext[i_mesh][i_part];

        double *_recv_buffer     = recv_buffer     [i_mesh][i_cloud][i_part];

        // CAUTION HERE : We fill the cloud !!!!
        double *_field           = field           [i_cloud][i_part];
        double *_blk_interp_from = blk_interp_from [i_cloud][i_part];
        double *_blk_interp_vol  = blk_interp_vol  [i_cloud][i_part];
        // double *_cell_nat        = cell_nat        [i_cloud][i_part];

        for(int i = 0; i < n_ref_tgt[i_part]; ++i) {
          int i_tgt = ref_tgt[i_part][i]-1;
          for(int j = gnum1_come_from_idx[i_part][i]; j < gnum1_come_from_idx[i_part][i+1]; ++j) {

            double lvol      = _recv_buffer[n_data_send*j+2];
            double ldist     = _recv_buffer[n_data_send*j+3];
            double lcell_nat = _recv_buffer[n_data_send*j+4];
            // printf("dist = %12.5e | vol = %12.5e \n", dist, vol);

            // On doit pouvoir gerer l'extrapolation également comme pour la condition du volume
            if(lcell_nat < 0 && ldist < 0 && lvol < _blk_interp_vol[i_tgt]) {
              _field          [i_tgt] = _recv_buffer[n_data_send*j  ];
              // printf("i_tgt = %i | fld = %12.5e | interp_from = %12.5e \n", i_tgt, _recv_buffer[n_data_send*j  ], _recv_buffer[n_data_send*j+1]);
              _blk_interp_from[i_tgt] = _recv_buffer[n_data_send*j+1];
              _blk_interp_vol [i_tgt] = PDM_MIN(_blk_interp_vol [i_tgt], lvol);
            }
          }
        }
      }
    }
  }
  t2 = PDM_MPI_Wtime();
  delta_t[5] = t2 - t1;

  /*
   * Cleaning
   */
  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    for(int i_cloud = 0; i_cloud < n_tot_cloud; ++i_cloud) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        free(send_buffer[i_mesh][i_cloud][i_part]);
        free(recv_buffer[i_mesh][i_cloud][i_part]);
      }
      free(send_buffer [i_mesh][i_cloud]);
      free(recv_buffer [i_mesh][i_cloud]);
    }
    free(send_buffer [i_mesh]);
    free(recv_buffer [i_mesh]);
    free(request_exch[i_mesh]);
    // free(recv_request[i_mesh]);
  }
  free(send_buffer);
  free(recv_buffer);
  free(request_exch);
  // free(recv_request);

  /*
   * Free ptp
   */
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    for(int i_cloud = 0; i_cloud < n_tot_cloud; ++i_cloud) {
      PDM_part_to_part_free(ptp[i_mesh][i_cloud]);
    }
    free(ptp[i_mesh]);
  }
  free(ptp);


  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    PDM_mesh_location_free(mesh_loc[i_mesh]);
  }
  t2 = PDM_MPI_Wtime();
  delta_t[6] = t2 - t1;

  /* Visu */
  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();
  if(post) {
    for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
      char filename[999];
      sprintf(filename, "mesh_%i", i_mesh);
      _visu (filename,
             n_part,
             n_cell_without_ext[i_mesh],
             n_face_without_ext[i_mesh],
             n_vtx_without_ext [i_mesh],
             cell_face_idx     [i_mesh],
             cell_face         [i_mesh],
             cell_ln_to_gn     [i_mesh],
             face_vtx_idx      [i_mesh],
             face_vtx          [i_mesh],
             face_ln_to_gn     [i_mesh],
             vtx_coord         [i_mesh],
             vtx_ln_to_gn      [i_mesh],
             field             [i_mesh],
             blk_interp_from   [i_mesh],
             cell_nat          [i_mesh],
             mask              [i_mesh]);
    }
  }
  t2 = PDM_MPI_Wtime();
  delta_t[7] = t2 - t1;


  PDM_MPI_Allreduce (delta_t, cpu_time_max, n_timer, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
  PDM_MPI_Allreduce (delta_t, cpu_time_min, n_timer, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce (delta_t, cpu_time_sum, n_timer, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  PDM_g_num_t n_interp_tot = -1;
  PDM_MPI_Allreduce (&n_interp_loc, &n_interp_tot, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);


  if(i_rank == 0) {
    printf("Total number of cells         : "PDM_FMT_G_NUM" | Mean by proc = %i\n", n_cell_tot  , (int) n_cell_tot/n_rank);
    printf("Total number of interpolation : "PDM_FMT_G_NUM" | Mean by proc = %i\n", n_interp_tot, (int) n_interp_tot/n_rank);
    printf("duration min/max : Generate cube   = %12.5e %12.5e | Adim (cell) = %12.5e | Adim (interp)  = %12.5e \n", cpu_time_min[0], cpu_time_max[0], cpu_time_sum[0]/((double) n_cell_tot), cpu_time_sum[0]/((double) n_interp_tot));
    printf("duration min/max : Ray traicing    = %12.5e %12.5e | Adim (cell) = %12.5e | Adim (interp)  = %12.5e \n", cpu_time_min[1], cpu_time_max[1], cpu_time_sum[1]/((double) n_cell_tot), cpu_time_sum[1]/((double) n_interp_tot));
    printf("duration min/max : Localisation    = %12.5e %12.5e | Adim (cell) = %12.5e | Adim (interp)  = %12.5e \n", cpu_time_min[2], cpu_time_max[2], cpu_time_sum[2]/((double) n_cell_tot), cpu_time_sum[2]/((double) n_interp_tot));
    printf("duration min/max : ptp create      = %12.5e %12.5e | Adim (cell) = %12.5e | Adim (interp)  = %12.5e \n", cpu_time_min[3], cpu_time_max[3], cpu_time_sum[3]/((double) n_cell_tot), cpu_time_sum[3]/((double) n_interp_tot));
    printf("duration min/max : ptp iexch       = %12.5e %12.5e | Adim (cell) = %12.5e | Adim (interp)  = %12.5e \n", cpu_time_min[4], cpu_time_max[4], cpu_time_sum[4]/((double) n_cell_tot), cpu_time_sum[4]/((double) n_interp_tot));
    printf("duration min/max : ptp wait + post = %12.5e %12.5e | Adim (cell) = %12.5e | Adim (interp)  = %12.5e \n", cpu_time_min[5], cpu_time_max[5], cpu_time_sum[5]/((double) n_cell_tot), cpu_time_sum[5]/((double) n_interp_tot));
    printf("duration min/max : Cleaning        = %12.5e %12.5e | Adim (cell) = %12.5e | Adim (interp)  = %12.5e \n", cpu_time_min[6], cpu_time_max[6], cpu_time_sum[6]/((double) n_cell_tot), cpu_time_sum[6]/((double) n_interp_tot));
    printf("duration min/max : Visualisation   = %12.5e %12.5e | Adim (cell) = %12.5e | Adim (interp)  = %12.5e \n", cpu_time_min[7], cpu_time_max[7], cpu_time_sum[7]/((double) n_cell_tot), cpu_time_sum[7]/((double) n_interp_tot));
  }

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    for(int i_cloud = 0; i_cloud < n_tot_cloud; ++i_cloud) {
      free(mesh_pts_idx           [i_mesh][i_cloud]);
      free(mesh_pts_gnum          [i_mesh][i_cloud]);
      free(points_coords          [i_mesh][i_cloud]);
      free(points_uvw             [i_mesh][i_cloud]);
      free(points_weights_idx     [i_mesh][i_cloud]);
      free(points_weights         [i_mesh][i_cloud]);
      free(points_dist2           [i_mesh][i_cloud]);
      free(points_projected_coords[i_mesh][i_cloud]);
    }
    free(mesh_pts_idx           [i_mesh]);
    free(mesh_pts_gnum          [i_mesh]);
    free(points_coords          [i_mesh]);
    free(points_uvw             [i_mesh]);
    free(points_weights_idx     [i_mesh]);
    free(points_weights         [i_mesh]);
    free(points_dist2           [i_mesh]);
    free(points_projected_coords[i_mesh]);
  }
  free(mesh_pts_idx           );
  free(mesh_pts_gnum          );
  free(points_coords          );
  free(points_uvw             );
  free(points_weights_idx     );
  free(points_weights         );
  free(points_dist2           );
  free(points_projected_coords);

  free(mesh_loc);

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      free(cell_face_idx  [i_mesh][i_part]);
      free(cell_face      [i_mesh][i_part]);
      free(face_vtx_idx   [i_mesh][i_part]);
      free(face_vtx       [i_mesh][i_part]);
      free(vtx_coord      [i_mesh][i_part]);
      free(cell_ln_to_gn  [i_mesh][i_part]);
      free(face_ln_to_gn  [i_mesh][i_part]);
      free(vtx_ln_to_gn   [i_mesh][i_part]);
      free(face_group_idx [i_mesh][i_part]);
      free(face_group     [i_mesh][i_part]);
      free(group_ln_to_gn [i_mesh][i_part]);
      free(cell_center    [i_mesh][i_part]);
      free(cell_volume    [i_mesh][i_part]);
      free(face_center    [i_mesh][i_part]);
      free(field          [i_mesh][i_part]);
      free(grad_field     [i_mesh][i_part]);
      free(blk_interp_from[i_mesh][i_part]);
      free(blk_interp_vol [i_mesh][i_part]);
      free(cell_nat       [i_mesh][i_part]);

      free(external_face_ln_to_gn[i_mesh][i_part]);
      free(external_vtx_ln_to_gn [i_mesh][i_part]);
      free(external_vtx_coord    [i_mesh][i_part]);
      free(external_face_vtx_idx [i_mesh][i_part]);
      free(external_face_vtx     [i_mesh][i_part]);

      if(extract_center_ln_to_gn[i_mesh][i_part] != NULL) {
        free(extract_center_ln_to_gn[i_mesh][i_part]);
        free(extract_center_coord[i_mesh][i_part]);
      }
      if(extract_face_bnd_ln_to_gn[i_mesh][i_part] != NULL) {
        free(extract_face_bnd_ln_to_gn[i_mesh][i_part]);
        free(extract_face_bnd_coord[i_mesh][i_part]);
      }

    }
    free(n_cell         [i_mesh]);
    free(n_face         [i_mesh]);
    free(n_vtx          [i_mesh]);

    free(n_cell_without_ext[i_mesh]);
    free(n_face_without_ext[i_mesh]);
    free(n_vtx_without_ext [i_mesh]);


    free(n_cell_ext     [i_mesh]);
    free(n_face_ext     [i_mesh]);
    free(n_vtx_ext      [i_mesh]);
    free(cell_face_idx  [i_mesh]);
    free(cell_face      [i_mesh]);
    free(face_vtx_idx   [i_mesh]);
    free(face_vtx       [i_mesh]);
    free(vtx_coord      [i_mesh]);
    free(cell_ln_to_gn  [i_mesh]);
    free(face_ln_to_gn  [i_mesh]);
    free(vtx_ln_to_gn   [i_mesh]);
    free(face_group_idx [i_mesh]);
    free(face_group     [i_mesh]);
    free(group_ln_to_gn [i_mesh]);
    free(cell_center    [i_mesh]);
    free(cell_volume    [i_mesh]);
    free(face_center    [i_mesh]);
    free(field          [i_mesh]);
    free(grad_field     [i_mesh]);
    free(blk_interp_from[i_mesh]);
    free(blk_interp_vol [i_mesh]);
    free(cell_nat       [i_mesh]);

    free(external_face_ln_to_gn[i_mesh]);
    free(external_vtx_ln_to_gn [i_mesh]);
    free(external_vtx_coord    [i_mesh]);
    free(external_face_vtx_idx [i_mesh]);
    free(external_face_vtx     [i_mesh]);

    free(n_extract_cell[i_mesh]);
    free(n_extract_face[i_mesh]);
    free(extract_center_ln_to_gn  [i_mesh]);
    free(extract_center_coord     [i_mesh]);
    free(extract_face_bnd_ln_to_gn[i_mesh]);
    free(extract_face_bnd_coord   [i_mesh]);

    free(n_external_face[i_mesh]);
    free(n_external_vtx [i_mesh]);
  }

  for(int i_cloud = 0; i_cloud < n_mesh; ++i_cloud) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(mask           [i_cloud][i_part]);
      free(elmt_cross_surf[i_cloud][i_part]);
    }
    free(mask           [i_cloud]);
    free(elmt_cross_surf[i_cloud]);
  }
  free(mask);
  free(elmt_cross_surf);

  free(n_cell        );
  free(n_face        );
  free(n_vtx         );

  free(n_cell_without_ext);
  free(n_face_without_ext);
  free(n_vtx_without_ext );

  free(n_cell_ext    );
  free(n_face_ext    );
  free(n_vtx_ext     );
  free(cell_face_idx );
  free(cell_face     );
  free(face_vtx_idx  );
  free(face_vtx      );
  free(vtx_coord     );
  free(cell_ln_to_gn );
  free(face_ln_to_gn );
  free(vtx_ln_to_gn  );
  free(face_group_idx);
  free(face_group    );
  free(group_ln_to_gn);
  free(cell_center);
  free(cell_volume);
  free(face_center);
  free(field          );
  free(grad_field     );
  free(blk_interp_from);
  free(blk_interp_vol);
  free(cell_nat);

  free(n_extract_cell);
  free(n_extract_face);
  free(n_external_face);
  free(n_external_vtx );
  free(extract_center_ln_to_gn  );
  free(extract_center_coord     );
  free(extract_face_bnd_ln_to_gn);
  free(extract_face_bnd_coord   );

  free(external_face_ln_to_gn);
  free(external_vtx_ln_to_gn );
  free(external_vtx_coord    );
  free(external_face_vtx_idx );
  free(external_face_vtx     );

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}

