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
#include "pdm_dcube_gen.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include "pdm_array.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"

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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_vtx_seg,
           double        *length,
           double        *separation_x,
           double        *separation_y,
           double        *separation_z,
           int           *deform,
           double        *tolerance,
           double        *marge,
           int           *n_part,
           PDM_g_num_t   *n_pts,
           int           *post,
           int           *part_method,
           PDM_mesh_location_method_t *loc_method,
           int           *disable_uvw,
           int           *use_tgt_nodes)
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

static PDM_part_t *
_cube_mesh
(
 const int              n_part,
 const PDM_part_split_t part_method,
 const PDM_g_num_t      n_vtx_seg,
 const double           xmin,
 const double           ymin,
 const double           zmin,
 const double           length,
 const int              deform
 )
{
  PDM_dcube_t *dcube = PDM_dcube_gen_init (PDM_MPI_COMM_WORLD,
                                           n_vtx_seg,
                                           length,
                                           xmin,
                                           ymin,
                                           zmin,
                                           PDM_OWNERSHIP_KEEP);

  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  int          n_face_group;
  PDM_g_num_t *dface_cell = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx = NULL;
  double      *dvtx_coord = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group = NULL;
  int          dface_vtx_l;
  int          dface_group_l;


  PDM_dcube_gen_dim_get (dcube,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtx_l,
                         &dface_group_l);

  PDM_dcube_gen_data_get (dcube,
                          &dface_cell,
                          &dface_vtx_idx,
                          &dface_vtx,
                          &dvtx_coord,
                          &dface_group_idx,
                          &dface_group);

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
  // int ppart_id = 0;
  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_part_t *ppart = PDM_part_create (PDM_MPI_COMM_WORLD,
                                       part_method,
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

  PDM_dcube_gen_free(dcube);

  return ppart;
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

  PDM_g_num_t n_vtx_seg    = 10;
  double      length       = 1.;
  double      separation_x = 2.;
  double      separation_y = 0.;
  double      separation_z = 0.;
  int         deform       = 0;
  double      tolerance    = 1e-6;
  double      marge        = 0.;
  int         n_part       = 1;
  int         post         = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t part_method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t part_method  = PDM_PART_SPLIT_PTSCOTCH;
#else
  PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;
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
              &use_tgt_nodes);


  /*
   *  Init
   */
  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);


  const double xmin = 0;
  const double ymin = 0;
  const double zmin = 0;

  /*
   *  Source cube
   */
  PDM_part_t *ppart_src = _cube_mesh (n_part,
                                      part_method,
                                      n_vtx_seg,
                                      xmin,
                                      ymin,
                                      zmin,
                                      length,
                              0);//deform);

  /*
   *  Target cube
   */
  PDM_part_t *ppart_tgt = _cube_mesh (n_part,
                                      part_method,
                                      n_vtx_seg,
                                      xmin + separation_x*length,
                                      ymin + separation_y*length,
                                      zmin + separation_z*length,
                                      length,
                                      deform);

  /*
   *  Mesh location structure initialization
   */
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED,
                                                            1,
                                                            PDM_MPI_COMM_WORLD);


  /* Set target point cloud */
  PDM_mesh_location_n_part_cloud_set (mesh_loc,
                                      0,
                                      n_part);

  /*double **cell_volume = malloc (sizeof(double *) * n_part);
    double **cell_center = malloc (sizeof(double *) * n_part);*/

  int *n_tgt = malloc (sizeof(double *) * n_part);
  PDM_g_num_t **tgt_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  double      **tgt_coord = malloc (sizeof(double *)      * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
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

    PDM_part_part_dim_get (ppart_tgt,
                           ipart,
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

    int         *cell_tag;
    int         *cell_face_idx;
    int         *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int         *face_tag;
    int         *face_cell;
    int         *face_vtx_idx;
    int         *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int         *face_part_boundProcIdx;
    int         *face_part_boundPartIdx;
    int         *face_part_bound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *face_group_idx;
    int         *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart_tgt,
                           ipart,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_boundProcIdx,
                           &face_part_boundPartIdx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    if (use_tgt_nodes) {
      n_tgt[ipart] = n_vtx;
      tgt_g_num[ipart] = malloc (sizeof(PDM_g_num_t) * n_vtx);
      memcpy (tgt_g_num[ipart], vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_vtx);

      tgt_coord[ipart] = malloc (sizeof(double) * n_vtx * 3);
      memcpy (tgt_coord[ipart], vtx, sizeof(double) * n_vtx * 3);
    }
    else {
      n_tgt[ipart] = n_cell;
      tgt_g_num[ipart] = malloc (sizeof(PDM_g_num_t) * n_cell);
      memcpy (tgt_g_num[ipart], cell_ln_to_gn, sizeof(PDM_g_num_t) * n_cell);

      const int is_oriented = 0;
      double *cell_volume = malloc(sizeof(double) * n_cell);
      tgt_coord[ipart] = malloc(sizeof(double) * 3 * n_cell);

      PDM_geom_elem_polyhedra_properties (is_oriented,
                                          n_cell,
                                          n_face,
                                          face_vtx_idx,
                                          face_vtx,
                                          cell_face_idx,
                                          cell_face,
                                          n_vtx,
                                          vtx,
                                          cell_volume,
                                          tgt_coord[ipart],
                                          NULL,
                                          NULL);
      free (cell_volume);
    }

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

  int *n_src = malloc (sizeof(double *) * n_part);
  PDM_g_num_t **src_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {

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

    PDM_part_part_dim_get (ppart_src,
                           ipart,
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

    int         *cell_tag;
    int         *cell_face_idx;
    int         *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int         *face_tag;
    int         *face_cell;
    int         *face_vtx_idx;
    int         *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int         *face_part_boundProcIdx;
    int         *face_part_boundPartIdx;
    int         *face_part_bound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *face_group_idx;
    int         *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart_src,
                           ipart,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_boundProcIdx,
                           &face_part_boundPartIdx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    n_src[ipart] = n_cell;
    src_g_num[ipart] = malloc (sizeof(PDM_g_num_t) * n_cell);
    memcpy (src_g_num[ipart], cell_ln_to_gn, sizeof(PDM_g_num_t) * n_cell);

    PDM_mesh_location_part_set (mesh_loc,
                                ipart,
                                n_cell,
                                cell_face_idx,
                                cell_face,
                                cell_ln_to_gn,
                                n_face,
                                face_vtx_idx,
                                face_vtx,
                                face_ln_to_gn,
                                n_vtx,
                                vtx,
                                vtx_ln_to_gn);
  }


  /* Set location parameters */
  PDM_mesh_location_tolerance_set (mesh_loc,
                                   tolerance);

  PDM_mesh_location_method_set (mesh_loc,
                                loc_method);

  if (disable_uvw) {
    PDM_mesh_location_disable_uvw_computation (mesh_loc);
  }

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
  PDM_g_num_t **tgt_location = malloc (sizeof(PDM_g_num_t *) * n_part);
  double **tgt_proj_coord = malloc (sizeof(double *) * n_part);

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
  PDM_MPI_Allreduce (&n_wrong, &g_n_wrong, 1, PDM_MPI_INT, PDM_MPI_SUM, PDM_MPI_COMM_WORLD);
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

      for (int i = 0; i < n_src[ipart]; i++) {

        PDM_g_num_t ck = (src_g_num[ipart][i] - 1) / (n_cell_seg * n_cell_seg);
        PDM_g_num_t ci = (src_g_num[ipart][i] - 1) % n_cell_seg;
        PDM_g_num_t cj = (src_g_num[ipart][i] - 1 - ck*n_cell_seg*n_cell_seg) / n_cell_seg;

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


    PDM_MPI_Allreduce (&n_wrong, &g_n_wrong, 1, PDM_MPI_INT, PDM_MPI_SUM, PDM_MPI_COMM_WORLD);
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



  for (int ipart = 0; ipart < n_part; ipart++) {
    /*free (cell_center[ipart]);
      free (cell_volume[ipart]);*/
    free (tgt_g_num[ipart]);
    free (tgt_coord[ipart]);
    free (tgt_location[ipart]);
    free (tgt_proj_coord[ipart]);
    free (src_g_num[ipart]);
  }
  free (n_tgt);
  free (n_src);
  /*free (cell_center);
    free (cell_volume);*/
  free (tgt_g_num);
  free (tgt_coord);
  free (tgt_location);
  free (tgt_proj_coord);
  free (src_g_num);

  PDM_mesh_location_free (mesh_loc,
                          0);

  PDM_part_free (ppart_src);
  PDM_part_free (ppart_tgt);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return g_n_wrong;
}

