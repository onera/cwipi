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
           int           *rotation,
           double        *tolerance,
           double        *marge,
           int           *n_part,
           PDM_g_num_t   *n_pts,
           int           *post,
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
    else if (strcmp(argv[i], "-rot") == 0) {
        *rotation = 1;
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

  PDM_g_num_t n_vtx_seg = 10;
  double      length    = 1.;
  int         rotation  = 0;
  double      tolerance = 1e-6;
  double      marge     = 0.;
  int         n_part    = 1;
  int         post      = 0;
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

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &rotation,
             &tolerance,
             &marge,
             &n_part,
             &n_pts,
             &post,
             (int *) &part_method,
             &loc_method);


  /*
   *  Init
   */

  struct timeval t_elaps_debut;

  int my_rank;
  int n_procs;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &my_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_procs);

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

  /*
   *  Create distributed cube
   */

  int          id;

  const double xmin = 0;
  const double ymin = 0;
  const double zmin = 0;

  /*const double xmax = xmin + length;
    const double ymax = ymin + length;
    const double zmax = zmin + length;*/

  if (my_rank == 0) {
    printf("-- Build cube\n");
    fflush(stdout);
  }

  PDM_dcube_gen_init(&id,
                     PDM_MPI_COMM_WORLD,
                     n_vtx_seg,
                     length,
                     xmin,
                     ymin,
                     zmin);

  PDM_dcube_gen_dim_get(id,
                        &n_face_group,
                        &dn_cell,
                        &dn_face,
                        &dn_vtx,
                        &dface_vtx_l,
                        &dface_group_l);

  PDM_dcube_gen_data_get(id,
                         &dface_cell,
                         &dface_vtx_idx,
                         &dface_vtx,
                         &dvtx_coord,
                         &dface_group_idx,
                         &dface_group);

  if (rotation) {
    _rotate (dn_vtx,
             dvtx_coord);
  } 
  int ppart_id = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  if (my_rank == 0) {
    printf("-- Part\n");
    fflush(stdout);
  }

  PDM_part_create(&ppart_id,
                  PDM_MPI_COMM_WORLD,
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
  double xyz_max[3] = {length + marge, length + marge, length + marge};
  _random_cloud (n_pts,
                 xyz_min,
                 xyz_max,
                 n_procs,
                 my_rank,
                 &pts_coords,
                 &n_pts_l);

  if (rotation) {
    _rotate (n_pts_l,
             pts_coords);
  } 

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

  int id_loc = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED,//???
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

    PDM_part_part_dim_get(ppart_id,
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

    PDM_part_part_val_get (ppart_id,
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

    PDM_mesh_location_part_set (id_loc,
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
  PDM_mesh_location_tolerance_set (id_loc,
                                   tolerance);

  PDM_mesh_location_method_set (id_loc,
                                loc_method);


  /* Compute location */
  if (my_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  PDM_mesh_location_compute (id_loc);

  PDM_mesh_location_dump_times (id_loc);

  PDM_g_num_t *p_location    = NULL;
  int         *p_weights_idx = NULL;
  double      *p_weights     = NULL;
  double      *p_proj_coord  = NULL;
  PDM_mesh_location_get (id_loc,
                         0,//i_point_cloud,
                         0,//i_part,
                         &p_location,
                         &p_weights_idx,
                         &p_weights,
                         &p_proj_coord);


  /* Check results */
if (!rotation) {
  if (my_rank == 0) {
    printf("-- Check\n");
    fflush(stdout);
  }

  const double location_tolerance = 1.e-6;

  const PDM_g_num_t n_cell_seg = n_vtx_seg - 1;
  const double cell_side = length / ((double) n_cell_seg);

  for (int ipt = 0; ipt < n_pts_l; ipt++) {
    double *p = pts_coords + 3*ipt;

    int i = (int) floor (p[0] / cell_side);
    int j = (int) floor (p[1] / cell_side);
    int k = (int) floor (p[2] / cell_side);

    PDM_g_num_t box_gnum = 1 + i + n_cell_seg*(j + n_cell_seg*k);

    if (p[0] < -tolerance || p[0] > length + tolerance ||
        p[1] < -tolerance || p[1] > length + tolerance ||
        p[2] < -tolerance || p[2] > length + tolerance) {
      box_gnum = -1;
    }

    //printf("%d: ("PDM_FMT_G_NUM") | ("PDM_FMT_G_NUM")\n", ipt, p_location[ipt], box_gnum);
    if (p_location[ipt] != box_gnum) {
      double *cp = p_proj_coord + 3*ipt;
      printf("%d ("PDM_FMT_G_NUM") (%.15lf %.15lf %.15lf): ("PDM_FMT_G_NUM") | ("PDM_FMT_G_NUM") proj : (%.15lf %.15lf %.15lf)\n",
             ipt, pts_gnum[ipt],
             p[0], p[1], p[2],
             p_location[ipt], box_gnum,
             cp[0], cp[1], cp[2]);
      printf("weights =");
      for (int l = p_weights_idx[ipt]; l < p_weights_idx[ipt+1]; l++) {
        printf(" %f", p_weights[l]);
      }
      printf("\n");

      //-->>
      double cell_min[3] = {cell_side * i,     cell_side * j,     cell_side * k};
      double cell_max[3] = {cell_side * (i+1), cell_side * (j+1), cell_side * (k+1)};
      /*printf("cell min = (%.15lf %.15lf %.15lf)\ncell max = (%.15lf %.15lf %.15lf)\n",
             cell_min[0], cell_min[1], cell_min[2],
             cell_max[0], cell_max[1], cell_max[2]);*/

      double dist = HUGE_VAL;
      for (int idim = 0; idim < 3; idim++) {
        double _dist1 = PDM_ABS (p[idim] - cell_min[idim]);
        double _dist2 = PDM_ABS (p[idim] - cell_max[idim]);
        double _dist = PDM_MIN (_dist1, _dist2);
        dist = PDM_MIN (dist, _dist);
      }
      printf("distance = %e\n\n", dist);
      assert (dist < location_tolerance);
      //<<--
    }
  }
}

  PDM_mesh_location_free (id_loc,
                          0);

  PDM_part_free (ppart_id);

  free (pts_coords);
  free (pts_gnum);

  PDM_dcube_gen_free (id);

  PDM_MPI_Finalize();

  if (my_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
