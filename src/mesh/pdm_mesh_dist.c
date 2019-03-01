/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_mesh_dist.h"
#include "pdm_mesh_nodal.h"
#include "pdm_surf_mesh.h"
#include "pdm_handles.h"
#include "pdm_octree.h"
#include "pdm_dbbtree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_triangle.h"
#include "pdm_polygon.h"
#include "pdm_timer.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define NTIMER 7
  
/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \enum _ol_timer_step_t
 *
 */
 
typedef enum {

    BEGIN                         = 0,
    UPPER_BOUND_DIST              = 1,        
    CANDIDATE_SELECTION           = 2,
    LOAD_BALANCING_ELEM_DIST      = 3,
    COMPUTE_ELEM_DIST             = 4,
    RESULT_TRANSMISSION           = 5,        
    END                           = 6,        

} _ol_timer_step_t;

  
/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 * 
 */

typedef struct {

  int           n_part;
  int          *n_points;
  double      **coords;
  PDM_g_num_t **gnum;
  double      **dist;
  double      **proj;
  PDM_g_num_t **closest_elt_gnum;
  
} _points_cloud_t;

/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 * 
 */

typedef struct {

  int  n_point_cloud; /*!< Number of point clouds */
  PDM_MPI_Comm comm;  /*!< MPI communicator */ 

  PDM_mesh_nature_t mesh_nature;  /*!< Nature of the mesh */

  PDM_surf_mesh_t *surf_mesh;  /*!< Surface mesh pointer */
  
  int  mesh_nodal_id;  /*!< Surface mesh identifier */

  _points_cloud_t *points_cloud; /*!< Point clouds */

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */
  
  double times_cpu[NTIMER];     /*!< CPU time */
  
  double times_cpu_u[NTIMER];  /*!< User CPU time */
  
  double times_cpu_s[NTIMER];  /*!< System CPU time */

  
} _PDM_dist_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_dists   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/
  
/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _PDM_dist_t *
_get_from_id
(
 int  id
)
{
  _PDM_dist_t *dist = (_PDM_dist_t *) PDM_Handles_get (_dists, id);
    
  if (dist == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_dist error : Bad identifier\n");
  }

  return dist;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure to compute distance to a mesh nodal
 *
 * \param [in]   mesh_nature    Nature of the mesh
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 *
 * \return     Identifier
 */

int
PDM_mesh_dist_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Comm comm
)
{
  if (_dists == NULL) {
    _dists = PDM_Handles_create (4);
  }

  _PDM_dist_t *dist = (_PDM_dist_t *) malloc(sizeof(_PDM_dist_t));

  int id = PDM_Handles_store (_dists, dist);

  dist->mesh_nature = mesh_nature;
  dist->mesh_nodal_id = -1;
  dist->surf_mesh = NULL;
  dist->n_point_cloud = n_point_cloud;
  dist->comm = comm;
  dist->points_cloud =
    (_points_cloud_t*) malloc (sizeof(_points_cloud_t) * n_point_cloud);

  for (int i = 0; i <  n_point_cloud; i++) {
    dist->points_cloud[i].n_part = -1;
    dist->points_cloud[i].n_points = NULL;
    dist->points_cloud[i].coords = NULL;
    dist->points_cloud[i].gnum = NULL;
    dist->points_cloud[i].dist = NULL;
    dist->points_cloud[i].proj = NULL;
    dist->points_cloud[i].closest_elt_gnum = NULL;
  }

  dist->timer = PDM_timer_create ();
  
  for (int i = 0; i < NTIMER; i++) {
    dist->times_elapsed[i] = 0.;
    dist->times_cpu[i] = 0.;
    dist->times_cpu_u[i] = 0.;
    dist->times_cpu_s[i] = 0.;
  }
  
  return id;
}


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_mesh_dist_n_part_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_part
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  dist->points_cloud[i_point_cloud].n_part = n_part;
  dist->points_cloud[i_point_cloud].n_points =
    realloc(dist->points_cloud[i_point_cloud].n_points, n_part * sizeof(int));
  dist->points_cloud[i_point_cloud].coords =
    realloc(dist->points_cloud[i_point_cloud].coords, n_part * sizeof(double *));
  dist->points_cloud[i_point_cloud].gnum =
    realloc(dist->points_cloud[i_point_cloud].gnum, n_part * sizeof(PDM_g_num_t *));
  dist->points_cloud[i_point_cloud].dist =
    realloc(dist->points_cloud[i_point_cloud].dist, n_part * sizeof(double *));
  dist->points_cloud[i_point_cloud].proj =
    realloc(dist->points_cloud[i_point_cloud].proj, n_part * sizeof(double *));
  dist->points_cloud[i_point_cloud].closest_elt_gnum =
    realloc(dist->points_cloud[i_point_cloud].closest_elt_gnum,
            n_part * sizeof(PDM_g_num_t * ));
  
  for (int i = 0; i < n_part; i++) {
    dist->points_cloud[i_point_cloud].n_points[i] = -1;
    dist->points_cloud[i_point_cloud].coords[i] = NULL;
    dist->points_cloud[i_point_cloud].gnum[i] = NULL;
    dist->points_cloud[i_point_cloud].dist[i] = NULL;
    dist->points_cloud[i_point_cloud].proj[i] = NULL;
    dist->points_cloud[i_point_cloud].closest_elt_gnum[i] = NULL;
  }
}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number 
 *
 */

void
PDM_mesh_dist_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
 const int          n_points,
       double      *coords,
       PDM_g_num_t *gnum
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  dist->points_cloud[i_point_cloud].n_points[i_part] = n_points;
  dist->points_cloud[i_point_cloud].coords[i_part] = coords;
  dist->points_cloud[i_point_cloud].gnum[i_part] = gnum;
}


/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   id             Identifier
 * \param [in]   mesh_nodal_id  Mesh nodal identifier
 *
 */

void
PDM_mesh_dist_nodal_mesh_set
(
 const int  id,
 const int  mesh_nodal_id
)
{
  _PDM_dist_t *dist = _get_from_id (id);
  dist->mesh_nodal_id = mesh_nodal_id;
}


/**
 *
 * \brief Set global data of a surface mesh
 *
 * \param [in]   id             Identifier
 * \param [in]   n_g_face       Global number of faces
 * \param [in]   n_g_vtx        Global number of vertices
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_mesh_dist_surf_mesh_global_data_set
(
 const int         id,
 const PDM_g_num_t n_g_face,
 const PDM_g_num_t n_g_vtx,
 const int         n_part
)
{
  
  _PDM_dist_t *dist = _get_from_id (id);

  if (dist->surf_mesh != NULL) {
    dist->surf_mesh = PDM_surf_mesh_free (dist->surf_mesh);
  }

  dist->surf_mesh =
    PDM_surf_mesh_create (n_g_face, n_g_vtx, n_part, dist->comm);
                                       
}


/**
 *
 * \brief Set a part of a surface mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define  
 * \param [in]   n_face        Number of faces                     
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering 
 * \param [in]   n_vtx         Number of vertices              
 * \param [in]   coords        Coordinates       
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering 
 *
 */

void
PDM_mesh_dist_surf_mesh_part_set
(
 const int          id,
 const int          i_part,
 const int          n_face,
 const int         *face_vtx_idx,
 const int         *face_vtx,
 const PDM_g_num_t *face_ln_to_gn,
 const int          n_vtx, 
 const double      *coords,
 const PDM_g_num_t *vtx_ln_to_gn
)
{
  
  _PDM_dist_t *dist = _get_from_id (id);

  PDM_surf_mesh_part_input (dist->surf_mesh,
                            i_part,
                            n_face,
                            face_vtx_idx,
                            face_vtx,
                            face_ln_to_gn,
                            n_vtx,
                            coords,
                            vtx_ln_to_gn);
}

  
/**
 *
 * \brief Set a point cloud with initial distance
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   initial_dist    Initial distance  
 * \param [in]   coords          Point coordinates
 *
 */

/* void */
/* PDM_mesh_dist_cloud_with_initial_set */
/* ( */
/*  const int          id, */
/*  const int          i_point_cloud, */
/*  const int          i_part, */
/*  const int          n_points, */
/*  const double      *initial_dist, */
/*  const double      *coords */
/* ) */
/* { */
/* } */


/**
 *
 * \brief Process merge points
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_mesh_dist_process
(
 const int id
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  const int n_point_cloud = dist->n_point_cloud;
  const int mesh_id = dist->mesh_nodal_id;
  PDM_MPI_Comm comm = dist->comm;

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;
  
  PDM_timer_hang_on(dist->timer);  
  dist->times_elapsed[BEGIN] = PDM_timer_elapsed(dist->timer);
  dist->times_cpu[BEGIN]     = PDM_timer_cpu(dist->timer);
  dist->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(dist->timer);
  dist->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);
  
  /* 
   * For each cloud
   */

  for (int i_point_cloud = 0; i_point_cloud < n_point_cloud; i_point_cloud++) {

    _points_cloud_t *pt_cloud = &(dist->points_cloud[i_point_cloud]);
    const int n_part = pt_cloud->n_part; 
    
    /********************************************************************************* 
     * 
     * Compute the upper bound distance. It is the distance from the closest vertex
     *      - Store mesh vetices in a octree
     *      - Compute the closest distance of points to vertices stored in the octree
     *
     *********************************************************************************/

    PDM_timer_hang_on(dist->timer);  
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);
  
    const double tolerance = 1e-4;
    const int depth_max = 1000;
    const int points_in_leaf_max = 4; 

    int n_part_mesh = 0;
    if (dist->mesh_nodal_id != -1) {
      n_part_mesh = PDM_Mesh_nodal_n_part_get (mesh_id);
    }
    else if (dist->surf_mesh != NULL) {
      n_part_mesh = PDM_surf_mesh_n_part_get (dist->surf_mesh);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0,
                "PDM_mesh_dist error : The surface mesh is not defined. To do that : \n"
                "        Call PDM_mesh_dist_nodal_mesh_set or\n" 
                "        Call PDM_mesh_dist_surf_mesh_global_data_set +"
                " PDM_mesh_dist_surf_mesh_part_set\n");
    }
    
    int octree_id = PDM_octree_create (n_part_mesh,
                                       depth_max, 
                                       points_in_leaf_max,
                                       tolerance,
                                       comm);
    
    for (int i_part = 0; i_part < n_part_mesh; i_part++) {

      int n_vertices = 0;
      const double *vertices_coords = NULL;
      const PDM_g_num_t *vertices_gnum = NULL;

      if (dist->mesh_nodal_id != -1) {
        n_vertices      = PDM_Mesh_nodal_n_vertices_get (mesh_id, i_part);
        vertices_coords = PDM_Mesh_nodal_vertices_get (mesh_id, i_part);
        vertices_gnum   = PDM_Mesh_nodal_vertices_g_num_get (mesh_id, i_part);
      }
      else if (dist->surf_mesh != NULL) {
        n_vertices      = PDM_surf_mesh_part_n_vtx_get (dist->surf_mesh, i_part);
        vertices_coords = PDM_surf_mesh_part_vtx_get (dist->surf_mesh, i_part);
        vertices_gnum   = PDM_surf_mesh_part_vtx_g_num_get (dist->surf_mesh, i_part);
      }
      else {
        PDM_error(__FILE__, __LINE__, 0,
                  "PDM_mesh_dist error : The surface mesh is not defined. To do that : \n"
                  "        Call PDM_mesh_dist_nodal_mesh_set or\n" 
                  "        Call PDM_mesh_dist_surf_mesh_global_data_set +"
                  " PDM_mesh_dist_surf_mesh_part_set\n");
      }
      
      PDM_octree_point_cloud_set (octree_id, i_part, n_vertices,
                                  vertices_coords, vertices_gnum);

    }
    
    /*
     * Build octree
     */

    PDM_octree_build (octree_id);

    /*
     * Concatenation of the partitions
     */

    int n_pts_rank = 0;

    for (int i_part = 0; i_part < n_part; i_part++) {
      n_pts_rank += pt_cloud->n_points[i_part];
    }
    
    double *pts_rank = malloc (sizeof(double) * n_pts_rank * 3);
    PDM_g_num_t *pts_g_num_rank = malloc (sizeof(PDM_g_num_t) * n_pts_rank);

    n_pts_rank = 0;
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < pt_cloud->n_points[i_part]; i++) {
        for (int k = 0; k < 3; k++) {
          pts_rank[3*(n_pts_rank + i) + k] = pt_cloud->coords[i_part][3*i+k];
          pts_g_num_rank[n_pts_rank + i] = pt_cloud->gnum[i_part][i];
        }
      }
      n_pts_rank += pt_cloud->n_points[i_part];
    }

    /*
     * Look for closest surface mesh vertices
     */

    PDM_g_num_t *closest_vertices_gnum = malloc (sizeof(PDM_g_num_t) * n_pts_rank);

    double *closest_vertices_dist2 =  malloc (sizeof(double) * n_pts_rank);
    
    PDM_octree_closest_point (octree_id, 
                              n_pts_rank,
                              pts_rank,
                              pts_g_num_rank,
                              closest_vertices_gnum,
                              closest_vertices_dist2);

    free (closest_vertices_gnum);
    
    PDM_octree_free (octree_id);
    
    PDM_timer_hang_on(dist->timer);  
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[UPPER_BOUND_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[UPPER_BOUND_DIST]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[UPPER_BOUND_DIST]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[UPPER_BOUND_DIST]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

    
    /********************************************************************************* 
     * 
     * Compute bounding box structure to find candidates closest 
     *     to the upper bound distance
     *
     *********************************************************************************/

    PDM_timer_hang_on(dist->timer);  
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    PDM_dbbtree_t *dbbt = PDM_dbbtree_create (dist->comm, 3);

          int          *nElts   = malloc (sizeof(int) * n_part_mesh);
    const double      **extents = malloc (sizeof(double *) * n_part_mesh);
    const PDM_g_num_t **gNum    = malloc (sizeof(PDM_g_num_t *) * n_part_mesh);


    if (dist->mesh_nodal_id != -1) {
    
/*     int PDM_Mesh_nodal_n_blocks_get */
/* ( */
/* const int   idx */
/* ); */

/* int * */
/* PDM_Mesh_nodal_blocks_id_get */
/* ( */
/* const int   idx */
/* ); */

 
/* PDM_Mesh_nodal_elt_t */
/* PDM_Mesh_nodal_block_type_get */
/* ( */
/* const int   idx, */
/* const int   id_block      */
/* ); */

/* void */
/* PDM_Mesh_nodal_block_std_get  */
/* (    */
/* const int            idx, */
/* const int            id_block,      */
/* const int            id_part,  */
/*       PDM_l_num_t  **connec    */
/* );  */

/* int */
/* PDM_Mesh_nodal_block_n_elt_get  */
/* (    */
/* const int            idx, */
/* const int            id_block,      */
/* const int            id_part  */
/* ); */
 
/* PDM_g_num_t * */
/* PDM_Mesh_nodal_block_g_num_get  */
/* (    */
/* const int            idx, */
/* const int            id_block,      */
/* const int            id_part  */
/* );  */


/* void */
/* PDM_Mesh_nodal_block_poly2d_get  */
/* ( */
/*  const int          idx, */
/*  const int          id_block,  */
/*  const int          id_part,  */
/*        PDM_l_num_t  **connec_idx,    */
/*        PDM_l_num_t  **connec */
/* );  */

    }
    else if (dist->surf_mesh != NULL) {
      PDM_surf_mesh_compute_faceExtentsMesh (dist->surf_mesh, 1e-4);
      for (int i_part = 0; i_part < n_part_mesh; i_part++) {
        nElts[i_part] = PDM_surf_mesh_part_n_face_get (dist->surf_mesh,
                                                       i_part);
        
        gNum[i_part] = PDM_surf_mesh_part_edge_g_num_get (dist->surf_mesh,
                                                          i_part);

        extents[i_part] = PDM_surf_mesh_part_extents_get (dist->surf_mesh,
                                                          i_part);

      }
    }
    
    else {
      PDM_error(__FILE__, __LINE__, 0,
                "PDM_mesh_dist error : The surface mesh is not defined. To do that : \n"
                "        Call PDM_mesh_dist_nodal_mesh_set or\n" 
                "        Call PDM_mesh_dist_surf_mesh_global_data_set +"
                " PDM_mesh_dist_surf_mesh_part_set\n");
    }
 
    /* PDM_box_set_t  *surf_mesh_boxes = PDM_dbbtree_boxes_set (dbbt, */
    /*                                                          n_part_mesh, */
    /*                                                          nElts, */
    /*                                                          extents, */
    /*                                                          gNum); */
    /* 
     * Find elements closer than closest_vertices_dist2 distance
     */

    int         *box_index; 
    PDM_g_num_t *box_g_num;

    PDM_dbbtree_closest_upper_bound_dist_boxes_get (dbbt,
                                                    n_pts_rank,        
                                                    pts_rank,
                                                    pts_g_num_rank,
                                                    closest_vertices_dist2,
                                                    &box_index,  
                                                    &box_g_num);

    free (pts_rank);
    free (closest_vertices_dist2);
    free (pts_g_num_rank);

    PDM_dbbtree_free (dbbt);

    PDM_timer_hang_on(dist->timer);  
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[CANDIDATE_SELECTION] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[CANDIDATE_SELECTION]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[CANDIDATE_SELECTION]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[CANDIDATE_SELECTION]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

    /*********************************************************************************
     * 
     * Load balancing of elementary computations (distance from a point to an element)
     *
     *********************************************************************************/

    PDM_timer_hang_on(dist->timer);  
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    double *box_n = malloc (sizeof(double) * n_pts_rank);
    int *i_box_n = malloc (sizeof(int) * n_pts_rank);

    for (int i = 0; i < n_pts_rank; i++) {
      box_n[i] = (double) (box_index[i+1] - box_index[i]);
    }

    for (int i = 0; i < n_pts_rank; i++) {
      i_box_n[i] = box_index[i+1] - box_index[i];
    }

    PDM_part_to_block_t *ptb_vtx =
      PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                PDM_PART_TO_BLOCK_POST_MERGE,
                                1.,
                                &pts_g_num_rank,
                                &box_n,  
                                &n_pts_rank,
                                1,  
                                comm);

    int *stride = malloc (sizeof(int) * n_pts_rank);
    for (int i = 0; i < n_pts_rank; i++) {
      stride[i] = 3;
    }

    int *block_stride;
    double *block_pts;
    PDM_part_to_block_exch (ptb_vtx,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &stride,
                            (void **) &pts_rank,
                            &block_stride,
                            (void **) &block_pts);

    free (block_stride);

    for (int i = 0; i < n_pts_rank; i++) {
      stride[i] = 1;
    }

    int *block_g_num_stride;
    PDM_g_num_t *block_g_num;
    PDM_part_to_block_exch (ptb_vtx,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            1,
                            &i_box_n,
                            (void **) &box_g_num,
                            &block_g_num_stride,
                            (void **) &block_g_num);

    free (i_box_n);


    /*
     * Receive of needed elements
     *    - PDM_part_to_part function have to be write
     *    - This step is realised by a couple of  PDM_part_to_block and  PDM_block_to_part
     *
     */

    /* part to block */

    const PDM_g_num_t **gnum_face_mesh = malloc (sizeof(PDM_g_num_t *) * n_part_mesh);
    int *n_face_mesh = malloc (sizeof(int *) * n_part_mesh);
    
    for (int i = 0; i < n_part_mesh; i++) {
      n_face_mesh[i] = PDM_surf_mesh_part_n_face_get (dist->surf_mesh, i);
      gnum_face_mesh[i] = PDM_surf_mesh_part_face_g_num_get (dist->surf_mesh, i);      
    }

    PDM_part_to_block_t *ptb_elt =
      PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                PDM_PART_TO_BLOCK_POST_CLEANUP,
                                1.,
                                (PDM_g_num_t **) gnum_face_mesh,
                                NULL,  
                                n_face_mesh,
                                n_part_mesh,  
                                comm);

    double **coords_face_mesh = malloc (sizeof(double *) * n_part_mesh);
    int **coords_face_mesh_n = malloc (sizeof(int *) * n_part_mesh);
    
    for (int i = 0; i < n_part_mesh; i++) {
      coords_face_mesh[i] = NULL;
      const int *part_face_vtx     =
        PDM_surf_mesh_part_face_vtx_get (dist->surf_mesh, i);
      const int *part_face_vtx_idx =
        PDM_surf_mesh_part_face_vtx_idx_get (dist->surf_mesh, i);
      const double *part_vtx = PDM_surf_mesh_part_vtx_get (dist->surf_mesh, i);

      coords_face_mesh_n[i] = malloc (sizeof(int) * n_face_mesh[i]);
      coords_face_mesh[i] = malloc(sizeof(double) * 3 * part_face_vtx_idx[n_face_mesh[i]]);
      
      for (int j = 0; j < n_face_mesh[i]; j++) {
        coords_face_mesh_n[i][j] = (part_face_vtx_idx[j+1] - part_face_vtx_idx[j]) * 3;
      }

      int idx = 0;
      for (int j = 0; j < part_face_vtx_idx[n_face_mesh[i]]; j++) {
        int _vtx = part_face_vtx[j] - 1;
        for (int k = 0; k < 3; k++) {
          coords_face_mesh[i][idx++] = part_vtx[3*_vtx+k];
        }
      }

    }

    int *block_coords_face_mesh_n = NULL;
    double *block_coords_face_mesh = NULL;
    
    PDM_part_to_block_exch (ptb_elt,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            -1,
                            coords_face_mesh_n,
                            (void **) coords_face_mesh,
                            &block_coords_face_mesh_n,
                            (void **) &block_coords_face_mesh);


    /* block to part */
    
    PDM_g_num_t *block_face_distrib_idx = PDM_part_to_block_distrib_index_get (ptb_elt);

    int block_g_num_n = 0;
    for (int i = 0; i < n_pts_rank; i++) {
      block_g_num_n += block_g_num_stride[i]; 
    }
    
    PDM_block_to_part_t *btp = PDM_block_to_part_create (block_face_distrib_idx,
                                                         (const PDM_g_num_t **) &block_g_num,
                                                         &block_g_num_n,
                                                         1,
                                                         comm);

    int un = 1;
    int *part_coords_vtx_face_n = malloc (sizeof(int) * block_g_num_n);
    PDM_block_to_part_exch (btp,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            &un,
                            block_coords_face_mesh_n,
                            NULL,
                            (void **) &part_coords_vtx_face_n);

    int *part_coords_vtx_face_idx = malloc (sizeof(int) * (block_g_num_n+1));
    part_coords_vtx_face_idx[0] = 0;

    for (int i = 0; i < block_g_num_n; i++) {
      part_coords_vtx_face_idx[i+1] = part_coords_vtx_face_idx[i] +
                                      part_coords_vtx_face_n[i]/3;
    }
    
    double *part_coords_vtx_face = malloc (sizeof(double) * 3 *
                                           part_coords_vtx_face_idx[block_g_num_n]);
    PDM_block_to_part_exch (btp,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            block_coords_face_mesh_n,
                            block_coords_face_mesh,
                            &part_coords_vtx_face_n,
                            (void **) &part_coords_vtx_face);

    /* free */ 

    for (int i = 0; i < n_part_mesh; i++) {
      free (coords_face_mesh[i]);
      free (coords_face_mesh_n[i]);
    }
    free (coords_face_mesh);
    free (gnum_face_mesh);
    free (n_face_mesh);

    PDM_part_to_block_free (ptb_elt);
    
    PDM_timer_hang_on(dist->timer);  
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[LOAD_BALANCING_ELEM_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[LOAD_BALANCING_ELEM_DIST]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[LOAD_BALANCING_ELEM_DIST]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[LOAD_BALANCING_ELEM_DIST]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

    /****************************************************************************
     * 
     * compute distance min per points 
     *
     ***************************************************************************/

    PDM_timer_hang_on(dist->timer);  
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    double *block_closest_dist = malloc (sizeof(double) * n_pts_rank);
    double *block_closest_proj = malloc (sizeof(double) * 3 * n_pts_rank);
    PDM_g_num_t *block_closest_gnum = malloc (sizeof(PDM_g_num_t) *  n_pts_rank);

    int idx = 0;
    for (int i = 0; i < n_pts_rank; i++) {
      double *_pt_coords = block_pts + 3*i;
      double *_block_closest_proj = block_closest_proj + 3*i;
      double *_block_closest_dist = block_closest_dist + i;
      _block_closest_dist[0] = HUGE_VAL;
      
      PDM_g_num_t *_block_closest_gnum = block_closest_gnum + i; 
      
      for (int j = 0; j < block_g_num_stride[i]; j++) {
        int n_vtx_elt = (part_coords_vtx_face_idx[idx+1] - part_coords_vtx_face_idx[idx])/3;
        double *_coords_face_elt = part_coords_vtx_face + part_coords_vtx_face_idx[idx];
        PDM_g_num_t face_g_num = block_g_num[idx];

        double closestPoint[3];
        double minDist2;

        if (n_vtx_elt == 3) {

          PDM_triangle_status_t status =
            PDM_triangle_evaluate_position (_pt_coords,
                                            _coords_face_elt,
                                            closestPoint,
                                            &minDist2,
                                            NULL);
          if (status == PDM_TRIANGLE_DEGENERATED) {
            idx += 1;
            continue;
          }
        }

        else {
          PDM_polygon_status_t status =
            PDM_polygon_evaluate_position (_pt_coords,
                                           n_vtx_elt,
                                           _coords_face_elt,
                                           closestPoint,
                                           &minDist2);
                                         
          if (status == PDM_POLYGON_DEGENERATED) {
            idx += 1;
            continue;
          }

        }

        if (minDist2 < _block_closest_dist[0]) {
          _block_closest_dist[0] = minDist2;
          for (int k = 0; k < 3; k++) {
            _block_closest_proj[k+1] = closestPoint[k];
          }
          _block_closest_gnum[0] = face_g_num;
        }
        
        idx += 1;
      }
    }

    /* Free */
    
    free (part_coords_vtx_face_n);
    free (part_coords_vtx_face_idx);
    free (part_coords_vtx_face);
    
    PDM_block_to_part_free (btp);

    PDM_timer_hang_on(dist->timer);  
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[COMPUTE_ELEM_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[COMPUTE_ELEM_DIST]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[COMPUTE_ELEM_DIST]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[COMPUTE_ELEM_DIST]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

    /**********************************************************************************
     * 
     * Transfer results 
     *
     **********************************************************************************/

    PDM_timer_hang_on(dist->timer);  
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    PDM_g_num_t *block_vtx_distrib_idx = PDM_part_to_block_distrib_index_get (ptb_vtx);

    PDM_block_to_part_t *btp_vtx = PDM_block_to_part_create (block_vtx_distrib_idx,
                                                             (const long **) pt_cloud->gnum,
                                                             pt_cloud->n_points,
                                                             n_part,
                                                             comm);

    for (int i = 0; i < n_part; i++) {
      int npts = pt_cloud->n_points[i];
      pt_cloud->dist[i] = malloc (sizeof(double) * npts);
      pt_cloud->proj[i] = malloc (sizeof(double) * npts * 3);
      pt_cloud->closest_elt_gnum[i] = malloc (sizeof(PDM_g_num_t) * npts);
    }

    PDM_block_to_part_exch (btp,
                            sizeof(double),
                            PDM_STRIDE_CST,
                            &un,
                            block_closest_dist,
                            NULL,
                            (void **) pt_cloud->dist);

    int three = 3;
    PDM_block_to_part_exch (btp,
                            sizeof(double),
                            PDM_STRIDE_CST,
                            &three,
                            block_closest_proj,
                            NULL,
                            (void **) pt_cloud->proj);


    PDM_block_to_part_exch (btp,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST,
                            &un,
                            block_closest_gnum,
                            NULL,
                            (void **) pt_cloud->closest_elt_gnum);

    PDM_block_to_part_free (btp_vtx);
    PDM_part_to_block_free (ptb_vtx);

    free (block_closest_proj);
    free (block_closest_dist);
    free (block_closest_gnum);

    PDM_timer_hang_on(dist->timer);  
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[RESULT_TRANSMISSION] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[RESULT_TRANSMISSION]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[RESULT_TRANSMISSION]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[RESULT_TRANSMISSION]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);
  }

  PDM_timer_hang_on(dist->timer);  
  dist->times_elapsed[END] = PDM_timer_elapsed(dist->timer);
  dist->times_cpu[END]     = PDM_timer_cpu(dist->timer);
  dist->times_cpu_u[END]   = PDM_timer_cpu_user(dist->timer);
  dist->times_cpu_s[END]   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);
  
}


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   id                Identifier
 * \param [in]   i_point_cloud     Current cloud
 * \param [in]   i_part            Index of partition of the cloud
 * \param [out]  distance          Distance
 * \param [out]  projected         Projected point coordinates
 * \param [out]  closest_elt_g_num Global number of the closest element
 *
 */

void
PDM_mesh_dist_get
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
       double      **distance,
       double      **projected,
       PDM_g_num_t **closest_elt_gnum
)
{
 _PDM_dist_t *dist = _get_from_id (id);

 *distance = dist->points_cloud[i_point_cloud].dist[i_part];
 *projected = dist->points_cloud[i_point_cloud].proj[i_part];
 *closest_elt_gnum = dist->points_cloud[i_point_cloud].closest_elt_gnum[i_part];
}


/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  id       Identifier
 * \param [in]  partial  if partial is equal to 0, all data are removed. 
 *                       Otherwise, results are kept. 
 *
 */

void
PDM_mesh_dist_free
(
 const int id,
 const int partial
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  if (!partial) {
    for (int i_point_cloud = 0; i_point_cloud < dist->n_point_cloud; i_point_cloud++) {
      for (int i = 0; i < (dist->points_cloud[i_point_cloud]).n_part; i++) {
        free (dist->points_cloud[i_point_cloud].dist[i]);
        free (dist->points_cloud[i_point_cloud].proj[i]);
        free (dist->points_cloud[i_point_cloud].closest_elt_gnum[i]); 
      }
    }
  }
  
  for (int i_point_cloud = 0; i_point_cloud < dist->n_point_cloud; i_point_cloud++) {
    free (dist->points_cloud[i_point_cloud].n_points);
    free (dist->points_cloud[i_point_cloud].coords);
    free (dist->points_cloud[i_point_cloud].gnum);
    free (dist->points_cloud[i_point_cloud].dist);
    free (dist->points_cloud[i_point_cloud].proj);
    free (dist->points_cloud[i_point_cloud].closest_elt_gnum);
  }
  
  free (dist->points_cloud);

  PDM_timer_free(dist->timer);

  if (dist->surf_mesh != NULL) {
    PDM_surf_mesh_free (dist->surf_mesh);
  }
  
  free (dist);
  
  PDM_Handles_handle_free (_dists, id, PDM_FALSE);
  
  const int n_dists = PDM_Handles_n_get (_dists);
  
  if (n_dists == 0) {
    _dists = PDM_Handles_free (_dists);
  }

}


/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_mesh_dump_times
(
 const int id
)
{
  _PDM_dist_t *dist = _get_from_id (id);
  double t1 = dist->times_elapsed[END] - dist->times_elapsed[BEGIN];
  double t2 = dist->times_cpu[END] - dist->times_cpu[BEGIN];
  
  PDM_printf( "distance times ALL (elapsed and cpu) : %12.5es %12.5es\n", t1, t2);
  PDM_printf( "distance times UPPER_BOUND_DIST (elapsed and cpu) : %12.5es %12.5es\n",
              dist->times_elapsed[UPPER_BOUND_DIST],
              dist->times_cpu[UPPER_BOUND_DIST]);
  PDM_printf( "distance times CANDIDATE_SELECTION (elapsed and cpu) : %12.5es %12.5es\n",
              dist->times_elapsed[CANDIDATE_SELECTION],
              dist->times_cpu[CANDIDATE_SELECTION]);
  PDM_printf( "distance times LOAD_BALANCING_ELEM_DIST (elapsed and cpu) : %12.5es %12.5es\n",
              dist->times_elapsed[LOAD_BALANCING_ELEM_DIST],
              dist->times_cpu[LOAD_BALANCING_ELEM_DIST]);
  PDM_printf( "distance times RESULT_TRANSMISSION (elapsed and cpu) : %12.5es %12.5es\n",
              dist->times_elapsed[RESULT_TRANSMISSION],
              dist->times_cpu[RESULT_TRANSMISSION]);

}

  
#ifdef	__cplusplus
}
#endif
