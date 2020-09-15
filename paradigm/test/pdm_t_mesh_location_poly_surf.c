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
#include "pdm_poly_surf_gen.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

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
     "  -pt-scocth       Call PT-Scotch.\n\n"
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
           double        *tolerance,
           double        *marge,
           int           *n_part,
           PDM_g_num_t   *n_pts,
           int           *post,
           int           *have_random,
           int           *init_random,
           int           *part_method,
           PDM_mesh_location_method_t *loc_method)
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

    else if (strcmp (argv[i], "-no_random") == 0) {
      *have_random = 0;
    }

    else if (strcmp (argv[i], "-random_init") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *init_random = atoi(argv[i]);
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
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_random_cloud
(
 const int      n_pts,
 const double   xyz_min[3],
 const double   xyz_max[3],
 const int      n_procs,
 const int      my_rank,
 double       **coord,
 int           *n_pts_l
 )
{
  double length[3] = {xyz_max[0] - xyz_min[0],
                      xyz_max[1] - xyz_min[1],
                      xyz_max[2] - xyz_min[2]};


  *n_pts_l = (int) (n_pts/n_procs);
  if (my_rank < n_pts%n_procs) {
    *n_pts_l += 1;
  }

  *coord = malloc (sizeof(double) * 3 * (*n_pts_l));
  double *_coord = *coord;
  double x;
  int idx = 0;
  for (PDM_g_num_t i = 0; i < n_procs*(*n_pts_l); i++) {
    for (int idim = 0; idim < 3; idim++) {
      x = xyz_min[idim] + length[idim] * (double) rand() / ((double) RAND_MAX);
      if (i%n_procs == my_rank) {
        _coord[idx++] = x;
      }
    }
  }
}




static void _z_elevation (const int  n_pts,
                          double    *coord)
{
  for (int i = 0; i < n_pts; i++) {
    double x = coord[3*i] - 0.5;
    double y = coord[3*i+1] - 0.5;
    double r = sqrt (x*x + y*y);
    PDM_UNUSED (r);
    coord[3*i+2] = 0.;//0.2 * cos(12*r) * exp(-5*r*r) + 0.3*(y + x);
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

  PDM_g_num_t n_vtx_seg   = 10;
  double      length      = 1.;
  double      tolerance   = 1e-6;
  double      marge       = 0.;
  int         n_part      = 1;
  int         post        = 0;
  int         have_random = 1;
  int         init_random = 0;
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

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_seg,
              &length,
              &tolerance,
              &marge,
              &n_part,
              &n_pts,
              &post,
              &have_random,
              &init_random,
              (int *) &part_method,
              &loc_method);


  /*
   *  Init
   */

  int my_rank;
  int n_ranks;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &my_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_ranks);


  /*
   *  Create distributed surface mesh
   */

  if (my_rank == 0) {
    printf("-- Build surface mesh\n");
    fflush(stdout);
  }

  const double xmin = 0.;
  const double ymin = 0.;
  const double xmax = xmin + length;
  const double ymax = ymin + length;

  const PDM_g_num_t nx = n_vtx_seg;
  const PDM_g_num_t ny = n_vtx_seg;


  PDM_g_num_t  nGFace;
  PDM_g_num_t  nGVtx;
  PDM_g_num_t  nGEdge;
  int          dn_vtx;
  double      *dvtx_coord    = NULL;
  int          dn_face;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *dFaceEdge     = NULL;
  int          dNEdge;
  PDM_g_num_t *dEdgeVtx      = NULL;
  PDM_g_num_t *dEdgeFace     = NULL;
  int          nEdgeGroup;
  int         *dEdgeGroupIdx = NULL;
  PDM_g_num_t *dEdgeGroup    = NULL;

  PDM_poly_surf_gen (PDM_MPI_COMM_WORLD,
                     xmin,
                     xmax,
                     ymin,
                     ymax,
                     have_random,
                     init_random,
                     nx,
                     ny,
                     &nGFace,
                     &nGVtx,
                     &nGEdge,
                     &dn_vtx,
                     &dvtx_coord,
                     &dn_face,
                     &dface_vtx_idx,
                     &dface_vtx,
                     &dFaceEdge,
                     &dNEdge,
                     &dEdgeVtx,
                     &dEdgeFace,
                     &nEdgeGroup,
                     &dEdgeGroupIdx,
                     &dEdgeGroup);

  _z_elevation (dn_vtx,
                dvtx_coord);

  /*
   *  Create mesh partitions
   */

  if (my_rank == 0) {
    printf("-- Part\n");
    fflush(stdout);
  }

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc (dn_face*sizeof(int));
  int *dEdgeVtxIdx = (int *) malloc ((dNEdge+1)*sizeof(int));

  dEdgeVtxIdx[0] = 0;
  for (int i = 0; i < dNEdge; i++) {
    dEdgeVtxIdx[i+1] = 2 + dEdgeVtxIdx[i];
  }

  int ppart_id;
  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_part_create (&ppart_id,
                   PDM_MPI_COMM_WORLD,
                   part_method,
                   "PDM_PART_RENUM_CELL_NONE",
                   "PDM_PART_RENUM_FACE_NONE",
                   n_property_cell,
                   renum_properties_cell,
                   n_property_face,
                   renum_properties_face,
                   n_part,
                   dn_face,
                   dNEdge,
                   dn_vtx,
                   nEdgeGroup,
                   NULL,
                   NULL,
                   NULL,
                   NULL,
                   have_dcell_part,
                   dcell_part,
                   dEdgeFace,
                   dEdgeVtxIdx,
                   dEdgeVtx,
                   NULL,
                   dvtx_coord,
                   NULL,
                   dEdgeGroupIdx,
                   dEdgeGroup);

  free (dcell_part);



  /************************
   *
   * Point cloud definition
   *
   ************************/
  if (my_rank == 0) {
    printf("-- Point cloud\n");
    fflush(stdout);
  }

  int n_pts_l;
  double *pts_coords = NULL;

  marge *= length;
  double xyz_min[3] = {-marge, -marge, -marge};
  double xyz_max[3] = {length + marge, length + marge, marge};
  _random_cloud (n_pts,
                 xyz_min,
                 xyz_max,
                 n_ranks,
                 my_rank,
                 &pts_coords,
                 &n_pts_l);

  _z_elevation (n_pts_l,
                pts_coords);


  /* Point cloud global numbering */
  int id_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD);

  double *char_length = malloc(sizeof(double) * n_pts_l);

  for (int i = 0; i < n_pts_l; i++) {
    char_length[i] = length * 1.e-6;
  }

  PDM_gnum_set_from_coords (id_gnum, 0, n_pts_l, pts_coords, char_length);

  PDM_gnum_compute (id_gnum);

  PDM_g_num_t *pts_gnum = PDM_gnum_get(id_gnum, 0);

  PDM_gnum_free (id_gnum, 1);
  free (char_length);








  /************************
   *
   * Mesh location struct initializaiton
   *
   ************************/

  int id_loc = PDM_mesh_location_create (PDM_MESH_NATURE_SURFACE_MESH,//???
                                         1,//const int n_point_cloud,
                                         PDM_MPI_COMM_WORLD);

  /* Set point cloud(s) */
  PDM_mesh_location_n_part_cloud_set (id_loc,
                                      0,//i_point_cloud,
                                      1);//n_part

  PDM_mesh_location_cloud_set (id_loc,
                               0,//i_point_cloud,
                               0,//i_part,
                               n_pts_l,
                               pts_coords,
                               pts_gnum);

  PDM_mesh_location_mesh_global_data_set (id_loc,
                                          n_part);

  /* Set mesh */
  for (int ipart = 0; ipart < n_part; ipart++) {

    int n_face;
    int n_edge;
    int n_edge_part_bound;
    int n_vtx;
    int n_proc;
    int n_t_part;
    int s_face_edge;
    int s_edge_vtx;
    int s_edge_group;
    int n_edge_group2;

    PDM_part_part_dim_get (ppart_id,
                          ipart,
                          &n_face,
                          &n_edge,
                          &n_edge_part_bound,
                          &n_vtx,
                          &n_proc,
                          &n_t_part,
                          &s_face_edge,
                          &s_edge_vtx,
                          &s_edge_group,
                          &n_edge_group2);

    int         *face_tag;
    int         *face_edge_idx;
    int         *face_edge;
    PDM_g_num_t *face_ln_to_gn;
    int         *edge_tag;
    int         *edge_face;
    int         *edge_vtx_idx;
    int         *edge_vtx;
    PDM_g_num_t *edge_ln_to_gn;
    int         *edge_part_boundProcIdx;
    int         *edge_part_boundPartIdx;
    int         *edge_part_bound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *edge_group_idx;
    int         *edge_group;
    PDM_g_num_t *edge_group_ln_to_gn;

    PDM_part_part_val_get (ppart_id,
                           ipart,
                           &face_tag,
                           &face_edge_idx,
                           &face_edge,
                           &face_ln_to_gn,
                           &edge_tag,
                           &edge_face,
                           &edge_vtx_idx,
                           &edge_vtx,
                           &edge_ln_to_gn,
                           &edge_part_boundProcIdx,
                           &edge_part_boundPartIdx,
                           &edge_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &edge_group_idx,
                           &edge_group,
                           &edge_group_ln_to_gn);

    PDM_mesh_location_part_set_2d (id_loc,
                                   ipart,
                                   n_face,
                                   face_edge_idx,
                                   face_edge,
                                   face_ln_to_gn,
                                   n_edge,
                                   edge_vtx_idx,
                                   edge_vtx,
                                   edge_ln_to_gn,
                                   n_vtx,
                                   vtx,
                                   vtx_ln_to_gn);
  }




  /* Set location parameters */
  PDM_mesh_location_tolerance_set (id_loc,
                                   tolerance);

  PDM_mesh_location_method_set (id_loc,
                                loc_method);


  /*
   * Compute location
   */
  if (my_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  PDM_mesh_location_compute (id_loc);

  PDM_mesh_location_dump_times (id_loc);





  /*
   * Check results
   */
  if (my_rank == 0) {
    printf("-- Check\n");
    fflush(stdout);
  }

  int p_n_points;
  double      *p_coords      = NULL;
  PDM_g_num_t *p_gnum        = NULL;
  PDM_g_num_t *p_location    = NULL;
  int         *p_weights_idx = NULL;
  double      *p_weights     = NULL;
  PDM_mesh_location_get (id_loc,
                         0,//i_point_cloud,
                         0,//i_part,
                         &p_n_points,
                         &p_coords,
                         &p_gnum,
                         &p_location,
                         &p_weights_idx,
                         &p_weights);
  if (0) {
    for (int ipt = 0; ipt < n_pts_l; ipt++) {
      printf("Point ("PDM_FMT_G_NUM") (%f %f %f), location = ("PDM_FMT_G_NUM"), weights =",
             p_gnum[ipt], p_coords[3*ipt], p_coords[3*ipt+1], p_coords[3*ipt+2], p_location[ipt]);
      for (int i = p_weights_idx[ipt]; i < p_weights_idx[ipt+1]; i++) {
        printf(" %f", p_weights[i]);
      }
      printf("\n");
    }
  }





  /*
   * Finalize
   */
  PDM_mesh_location_free (id_loc,
                          0);

  PDM_part_free (ppart_id);


  free (dvtx_coord);
  free (dface_vtx_idx);
  free (dface_vtx);
  free (dFaceEdge);
  free (dEdgeVtx);
  free (dEdgeFace);
  free (dEdgeGroupIdx);
  free (dEdgeGroup);
  free (dEdgeVtxIdx);

  PDM_MPI_Finalize();

  if (my_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
