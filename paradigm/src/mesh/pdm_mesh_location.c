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
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"
#include "pdm_surf_mesh.h"
#include "pdm_handles.h"
#include "pdm_dbbtree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_triangle.h"
#include "pdm_polygon.h"
#include "pdm_timer.h"
#include "pdm_hash_tab.h"
#include "pdm_mesh_location.h"
#include "pdm_point_location.h"

#include "pdm_binary_search.h"
#include "pdm_para_octree.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
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

  BEGIN                        = 0,
  BUILD_BOUNDING_BOXES         = 1,
  SEARCH_CANDIDATES            = 2,
  LOAD_BALANCING               = 3,
  COMPUTE_ELEMENTARY_LOCATIONS = 4,
  MERGE_LOCATION_DATA          = 5,
  END                          = 6,

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
  PDM_g_num_t **location;
  double      **uvw;
  int         **weights_idx;
  double      **weights; /*!< Barycentric coordinates */
  double      **projected_coords;

} _point_cloud_t;


/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  int  n_point_cloud; /*!< Number of point clouds */
  PDM_MPI_Comm comm;  /*!< MPI communicator */

  PDM_mesh_nature_t mesh_nature;  /*!< Nature of the mesh */

  int  shared_nodal;   /*!< 1 if mesh nodal is shared, 0 otherwise */
  int  mesh_nodal_id;  /*!< Mesh identifier */
  int _mesh_nodal_id;
  PDM_l_num_t **face_vtx_n; /* Mandatory to build mesh nodal */
  PDM_l_num_t **cell_face_n; /* Mandatory to build mesh nodal */

  _point_cloud_t *point_clouds; /*!< Point clouds */

  double tolerance;

  PDM_mesh_location_method_t method;

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */


} _PDM_location_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_locations   = NULL;

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

static _PDM_location_t *
_get_from_id
(
 int  id
 )
{
  _PDM_location_t *location = (_PDM_location_t *) PDM_Handles_get (_locations, id);

  if (location == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_location error : Bad identifier\n");
  }

  return location;
}



/*
 *
 * Redistribute evenly across all ranks the elementary location operation to perform
 *
 */

static void
_redistribute_elementary_location
(
 _PDM_location_t       *location,
 int                    n_elt,
 PDM_g_num_t            elt_g_num[],
 const int              pts_idx[],
 const PDM_g_num_t      pts_g_num[],
 const double           pts_coord[],
 int                   *r_n_elt,
 int                    r_type_idx[],
 PDM_g_num_t          **r_elt_g_num,
 int                  **r_vtx_idx,
 double               **r_vtx_coord,
 int                  **r_pts_idx,
 PDM_g_num_t          **r_pts_g_num,
 double               **r_pts_coord,
 PDM_l_num_t          **r_poly3d_face_idx,
 PDM_l_num_t          **r_face_vtx_idx,
 PDM_l_num_t          **r_face_vtx,
 int                  **r_face_orientation
 )
{
  const int order = 1;

  int n_blocks   = PDM_Mesh_nodal_n_blocks_get (location->mesh_nodal_id);
  int n_parts    = PDM_Mesh_nodal_n_part_get (location->mesh_nodal_id);
  int *blocks_id = PDM_Mesh_nodal_blocks_id_get (location->mesh_nodal_id);


  int my_rank, n_ranks;
  PDM_MPI_Comm_rank (location->comm, &my_rank);
  PDM_MPI_Comm_size (location->comm, &n_ranks);

  /*
   * Get element type and number of points to locate per element
   */
  PDM_Mesh_nodal_elt_t *elt_type = malloc (sizeof(PDM_Mesh_nodal_elt_t) * n_elt);
  int *n_pts_per_elt = malloc (sizeof(int) * n_elt);
  int n_poly3d = 0;

  int ielt = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {
    int id_block = blocks_id[iblock];
    PDM_Mesh_nodal_elt_t block_type = PDM_Mesh_nodal_block_type_get (location->mesh_nodal_id,
                                                                     id_block);

    for (int ipart = 0; ipart < n_parts; ipart++) {
      int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                       id_block,
                                                       ipart);
      if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        n_poly3d += n_elt_part;
      }

      for (int i = 0; i < n_elt_part; i++) {
        elt_type[ielt] = block_type;
        n_pts_per_elt[ielt] = pts_idx[ielt+1] - pts_idx[ielt];
        ielt++;
      }
    } // End of loop on parts
  } // End of loop on nodal blocks


  /*
   * Get number of vertices per element and face connectivity for polyhedra
   */
  ielt = 0;
  int *n_vtx_per_elt = malloc (sizeof(int) * n_elt);
  int *n_face_per_elt = malloc (sizeof(int) * n_poly3d);
  PDM_g_num_t *poly3d_g_num = malloc (sizeof(PDM_g_num_t) * n_poly3d);

  PDM_l_num_t *connec_idx = NULL;
  PDM_l_num_t *connec     = NULL;

  PDM_l_num_t  n_face;
  PDM_l_num_t *face_vtx_idx  = NULL;
  PDM_l_num_t *face_vtx      = NULL;
  PDM_l_num_t *cell_face_idx = NULL;
  PDM_l_num_t *cell_face     = NULL;

  int ipoly = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {
    int id_block = blocks_id[iblock];

    /* Polyhedra */
    if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                         id_block,
                                                         ipart);

        PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (location->mesh_nodal_id,
                                                          id_block,
                                                          ipart,
                                                          &connec_idx,
                                                          &connec);

        PDM_Mesh_nodal_block_poly3d_get (location->mesh_nodal_id,
                                         id_block,
                                         ipart,
                                         &n_face,
                                         &face_vtx_idx,
                                         &face_vtx,
                                         &cell_face_idx,
                                         &cell_face);

        for (int i = 0; i < n_elt_part; i++) {
          poly3d_g_num[ipoly] = elt_g_num[ielt];
          n_vtx_per_elt[ielt++] = connec_idx[i+1] - connec_idx[i];
          n_face_per_elt[ipoly++] = cell_face_idx[i+1] - cell_face_idx[i];
        }
      }
    }

    /* Polygons */
    else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                         id_block,
                                                         ipart);
        PDM_Mesh_nodal_block_poly2d_get (location->mesh_nodal_id,
                                         id_block,
                                         ipart,
                                         &connec_idx,
                                         &connec);

        for (int i = 0; i < n_elt_part; i++) {
          n_vtx_per_elt[ielt++] = connec_idx[i+1] - connec_idx[i];
        }
      }
    }

    /* Standard elements */
    else {
      PDM_Mesh_nodal_elt_t std_type = PDM_Mesh_nodal_block_type_get (location->mesh_nodal_id,
                                                                     id_block);
      int n_vtx = PDM_Mesh_nodal_n_vertices_element (std_type,
                                                     order);

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                         id_block,
                                                         ipart);

        for (int i = 0; i < n_elt_part; i++) {
          n_vtx_per_elt[ielt++] = n_vtx;
        }
      }
    }

  } // End of loop on nodal blocks



  /*
   * Get local face-vtx connectivity of polyhedra
   */
  PDM_l_num_t *local_face_vtx   = NULL;
  int         *n_vtx_per_face   = NULL;
  int         *face_orientation = NULL;
  if (n_poly3d > 0) {
    /* Total number of faces */
    int n_face_tot = 0;
    int n_face_max = 0;
    for (ipoly = 0; ipoly < n_poly3d; ipoly++) {
      n_face_tot += n_face_per_elt[ipoly];
      n_face_max = PDM_MAX (n_face_max, n_face_per_elt[ipoly]);
    }

    size_t s_local_face_vtx = 5 * n_face_tot;
    local_face_vtx = malloc (sizeof(PDM_l_num_t) * s_local_face_vtx);
    face_orientation = malloc (sizeof(int) * n_face_tot);
    n_vtx_per_face = malloc (sizeof(int) * n_face_tot);


    int idx_face = 0;
    int idx_vtx = 0;
    ipoly = 0;
    for (int iblock = 0; iblock < n_blocks; iblock++) {
      int id_block = blocks_id[iblock];

      /* Polyhedra */
      if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        for (int ipart = 0; ipart < n_parts; ipart++) {
          int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                           id_block,
                                                           ipart);

          PDM_Mesh_nodal_block_poly3d_get (location->mesh_nodal_id,
                                           id_block,
                                           ipart,
                                           &n_face,
                                           &face_vtx_idx,
                                           &face_vtx,
                                           &cell_face_idx,
                                           &cell_face);

          PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (location->mesh_nodal_id,
                                                            id_block,
                                                            ipart,
                                                            &connec_idx,
                                                            &connec);

          for (int i = 0; i < n_elt_part; i++) {
            int _n_vtx = connec_idx[i+1] - connec_idx[i];
            int _n_face = cell_face_idx[i+1] - cell_face_idx[i];

            for (int iface = 0; iface < _n_face; iface++) {
              int _iface = cell_face[cell_face_idx[i] + iface];

              if (_iface < 0) {
                _iface = -_iface - 1;
                face_orientation[idx_face] = -1;
              } else {
                _iface = _iface - 1;
                face_orientation[idx_face] = 1;
              }

              int n_vtx_face = face_vtx_idx[_iface+1] - face_vtx_idx[_iface];

              for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
                int _ivtx = PDM_binary_search_int (face_vtx[face_vtx_idx[_iface] + ivtx],
                                                   connec + connec_idx[i],
                                                   _n_vtx);
                assert (_ivtx >= 0);

                if ((int) s_local_face_vtx <= idx_vtx) {
                  s_local_face_vtx = PDM_MAX ((int) (2*s_local_face_vtx), idx_vtx);
                  local_face_vtx = realloc (local_face_vtx,
                                            sizeof(PDM_l_num_t) * s_local_face_vtx);
                }
                local_face_vtx[idx_vtx++] = _ivtx + 1;
              }

              n_vtx_per_face[idx_face++] = n_vtx_face;
            }

            ipoly++;
          } // End of loop on elements of current part

        } // End of loop on parts
      }
    } // End of loop on nodal blocks

    if (idx_vtx < (int) s_local_face_vtx) {
      local_face_vtx = realloc (local_face_vtx, sizeof(PDM_l_num_t) * idx_vtx);
    }
  }



  /*
   * Get vertex coordinates
   */
  /* Total number of vertices */
  int n_vtx_tot = 0;
  for (ielt = 0; ielt < n_elt; ielt++) {
    n_vtx_tot += n_vtx_per_elt[ielt];
  }

  double *vtx_coord = malloc (sizeof(double) * 3 * n_vtx_tot);
  double *_vtx_coord = vtx_coord;

  ipoly = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {
    int id_block = blocks_id[iblock];

    /* Polyhedra */
    if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                         id_block,
                                                         ipart);

        const double *vtx_coord_part = PDM_Mesh_nodal_vertices_get (location->mesh_nodal_id,
                                                                    ipart);

        PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (location->mesh_nodal_id,
                                                          id_block,
                                                          ipart,
                                                          &connec_idx,
                                                          &connec);

        for (int i = 0; i < n_elt_part; i++) {
          double xyz_min[3] = {HUGE_VAL,
                               HUGE_VAL,
                               HUGE_VAL};

          double xyz_max[3] = {-HUGE_VAL,
                               -HUGE_VAL,
                               -HUGE_VAL};

          for (int j = connec_idx[i]; j < connec_idx[i+1]; j++) {
            int ivtx = connec[j] - 1;
            for (int k = 0; k < 3; k++) {
              _vtx_coord[k] = vtx_coord_part[3*ivtx + k];
              xyz_min[k] = PDM_MIN (xyz_min[k], _vtx_coord[k]);
              xyz_max[k] = PDM_MAX (xyz_max[k], _vtx_coord[k]);
            }
            _vtx_coord += 3;
          }
        }
      } // End of loop on parts

    }

    /* Polygons */
    else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                         id_block,
                                                         ipart);

        const double *vtx_coord_part = PDM_Mesh_nodal_vertices_get (location->mesh_nodal_id,
                                                                    ipart);

        PDM_Mesh_nodal_block_poly2d_get (location->mesh_nodal_id,
                                         id_block,
                                         ipart,
                                         &connec_idx,
                                         &connec);

        for (int i = 0; i < n_elt_part; i++) {
          for (int j = connec_idx[i]; j < connec_idx[i+1]; j++) {
            int ivtx = connec[j] - 1;
            for (int k = 0; k < 3; k++) {
              _vtx_coord[k] = vtx_coord_part[3*ivtx + k];
            }
            _vtx_coord += 3;
          }
        }
      } // End of loop on parts

    }

    /* Standard elements */
    else {

      PDM_Mesh_nodal_elt_t std_type = PDM_Mesh_nodal_block_type_get (location->mesh_nodal_id,
                                                                     id_block);
      int n_vtx = PDM_Mesh_nodal_n_vertices_element (std_type,
                                                     order);


      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                         id_block,
                                                         ipart);

        const double *vtx_coord_part = PDM_Mesh_nodal_vertices_get (location->mesh_nodal_id,
                                                                    ipart);

        PDM_Mesh_nodal_block_std_get (location->mesh_nodal_id,
                                      id_block,
                                      ipart,
                                      &connec);

        for (int i = 0; i < n_elt_part; i++) {
          for (int j = 0; j < n_vtx; j++) {
            int ivtx = connec[n_vtx*i + j] - 1;
            for (int k = 0; k < 3; k++) {
              _vtx_coord[k] = vtx_coord_part[3*ivtx + k];
            }
            _vtx_coord += 3;
          }
        }
      } // End of loop on parts

    }

  } // End of loop on nodal blocks




  /* Compute elements weights for an even redistribution */
  double *elt_weight = malloc (sizeof(double) * n_elt);
  for (ielt = 0; ielt < n_elt; ielt++) {
    //elt_weight[ielt] = (double) n_pts_per_elt[ielt];
    elt_weight[ielt] = (double) n_pts_per_elt[ielt] * n_vtx_per_elt[ielt];
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &elt_g_num,
                                                       &elt_weight,
                                                       &n_elt,
                                                       1,
                                                       location->comm);
  free (elt_weight);

  *r_n_elt = PDM_part_to_block_n_elt_block_get (ptb);
  PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);


  /*
   * Exchange points to locate
   */

  /* Global number */
  int *block_n_pts_per_elt = NULL;
  PDM_g_num_t *block_pts_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &n_pts_per_elt,
                          (void **) &pts_g_num,
                          &block_n_pts_per_elt,
                          (void **) &block_pts_g_num);


  /* Coordinates */
  for (ielt = 0; ielt < n_elt; ielt++) {
    n_pts_per_elt[ielt] *= 3;
  }

  int *block_stride = NULL;
  double *block_pts_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          1,
                          &n_pts_per_elt,
                          (void **) &pts_coord,
                          &block_stride,
                          (void **) &block_pts_coord);
  free (block_stride);
  free (n_pts_per_elt);



  /*
   * Exchange elements
   */
  int *part_stride = malloc (sizeof(int) * n_elt);
  for (ielt = 0; ielt < n_elt; ielt++) {
    part_stride[ielt] = 1;
  }

  /* Type */
  PDM_Mesh_nodal_elt_t *block_elt_type = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_Mesh_nodal_elt_t),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &elt_type,
                          &block_stride,
                          (void **) &block_elt_type);
  free (block_stride);
  free (elt_type);

  /* Global number */
  PDM_g_num_t *block_elt_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &elt_g_num,
                          &block_stride,
                          (void **) &block_elt_g_num);
  free (block_stride);
  free (part_stride);

  /* Coordinates of vertices */
  for (ielt = 0; ielt < n_elt; ielt++) {
    n_vtx_per_elt[ielt] *= 3;
  }

  int *block_n_vtx_per_elt = NULL;
  double *block_vtx_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          1,
                          &n_vtx_per_elt,
                          (void **) &vtx_coord,
                          &block_n_vtx_per_elt,
                          (void **) &block_vtx_coord);
  free (vtx_coord);

  for (ielt = 0; ielt < *r_n_elt; ielt++) {
    block_n_vtx_per_elt[ielt] /= 3;
  }

  /*
   * Polyhedra
   */

  PDM_part_to_block_t *ptb_poly3d = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                               PDM_PART_TO_BLOCK_POST_MERGE,
                                                               1.,
                                                               &poly3d_g_num,
                                                               block_distrib_idx,
                                                               &n_poly3d,
                                                               1,
                                                               location->comm);

  int r_n_poly3d = PDM_part_to_block_n_elt_block_get (ptb_poly3d);

  /* Number of faces per polyhedron */
  part_stride = malloc (sizeof(int) * n_poly3d);
  for (ipoly = 0; ipoly < n_poly3d; ipoly++) {
    part_stride[ipoly] = 1;
  }

  int *r_n_face_per_elt = NULL;
  PDM_part_to_block_exch (ptb_poly3d,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &n_face_per_elt,
                          &block_stride,
                          (void **) &r_n_face_per_elt);
  free (part_stride);
  free (block_stride);

  *r_poly3d_face_idx = malloc (sizeof(PDM_l_num_t) * (r_n_poly3d + 1));
  (*r_poly3d_face_idx)[0] = 0;
  for (ipoly = 0; ipoly < r_n_poly3d; ipoly++) {
    (*r_poly3d_face_idx)[ipoly+1] = (*r_poly3d_face_idx)[ipoly] + r_n_face_per_elt[ipoly];
  }
  free (r_n_face_per_elt);

  int r_n_face = (*r_poly3d_face_idx)[r_n_poly3d];

  /* Face orientation */
  PDM_part_to_block_exch (ptb_poly3d,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          1,
                          &n_face_per_elt,
                          (void **) &face_orientation,
                          &block_stride,
                          (void **) r_face_orientation);
  free (block_stride);
  free (face_orientation);





  /* Number of vertices per face */
  int *r_n_vtx_per_face = NULL;
  PDM_part_to_block_exch (ptb_poly3d,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          1,
                          &n_face_per_elt,
                          (void **) &n_vtx_per_face,
                          &block_stride,
                          (void **) &r_n_vtx_per_face);
  free (block_stride);


  /* Face-vtx connectivity */
  *r_face_vtx_idx = malloc (sizeof(PDM_l_num_t) * (r_n_face + 1));
  (*r_face_vtx_idx)[0] = 0;
  for (int iface = 0; iface < r_n_face; iface++) {
    (*r_face_vtx_idx)[iface+1] = (*r_face_vtx_idx)[iface] + r_n_vtx_per_face[iface];
  }
  free (r_n_vtx_per_face);

  int idx_face = 0;
  for (ipoly = 0; ipoly < n_poly3d; ipoly++) {
    n_vtx_per_elt[ipoly] = 0;

    for (int iface = 0; iface < n_face_per_elt[ipoly]; iface++) {
      n_vtx_per_elt[ipoly] += n_vtx_per_face[idx_face++];
    }
  }
  free (n_vtx_per_face);
  free (n_face_per_elt);

  PDM_part_to_block_exch (ptb_poly3d,
                          sizeof(PDM_l_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &n_vtx_per_elt,
                          (void **) &local_face_vtx,
                          &block_stride,
                          (void **) r_face_vtx);
  free (block_stride);
  free (n_vtx_per_elt);
  free (local_face_vtx);

  free (poly3d_g_num);
  PDM_part_to_block_free (ptb);
  PDM_part_to_block_free (ptb_poly3d);




  /* Sort elements by type */
  int type_count[PDM_MESH_NODAL_N_ELEMENT_TYPES] = {0};
  for (ielt = 0; ielt < *r_n_elt; ielt++) {
    type_count[block_elt_type[ielt]]++;
  }


  r_type_idx[PDM_MESH_NODAL_POINT] = 0;
  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_POINT;
       type < PDM_MESH_NODAL_N_ELEMENT_TYPES;
       type++) {
    r_type_idx[type+1] = r_type_idx[type] + type_count[type];
    type_count[type] = 0;
  }


  *r_elt_g_num = malloc (sizeof(PDM_g_num_t) * (*r_n_elt));
  *r_vtx_idx = malloc (sizeof(int) * (*r_n_elt + 1));
  (*r_vtx_idx)[0] = 0;
  *r_pts_idx = malloc (sizeof(int) * (*r_n_elt + 1));
  (*r_pts_idx)[0] = 0;

  for (ielt = 0; ielt < *r_n_elt; ielt++) {
    PDM_Mesh_nodal_elt_t type = block_elt_type[ielt];
    int _ielt = r_type_idx[type] + type_count[type]++;

    (*r_elt_g_num)[_ielt] = block_elt_g_num[ielt];
    (*r_vtx_idx)[_ielt+1] = block_n_vtx_per_elt[ielt];
    (*r_pts_idx)[_ielt+1] = block_n_pts_per_elt[ielt];
  }
  free (block_elt_g_num);

  for (ielt = 0; ielt < *r_n_elt; ielt++) {
    (*r_vtx_idx)[ielt+1] += (*r_vtx_idx)[ielt];
    (*r_pts_idx)[ielt+1] += (*r_pts_idx)[ielt];
  }


  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_POINT;
       type < PDM_MESH_NODAL_N_ELEMENT_TYPES;
       type++) {
    type_count[type] = 0;
  }


  *r_vtx_coord = malloc (sizeof(double) * (*r_vtx_idx)[*r_n_elt] * 3);
  *r_pts_coord = malloc (sizeof(double) * (*r_pts_idx)[*r_n_elt] * 3);
  *r_pts_g_num = malloc (sizeof(PDM_g_num_t) * (*r_pts_idx)[*r_n_elt]);

  int idx_vtx = 0;
  int idx_pts = 0;
  for (ielt = 0; ielt < *r_n_elt; ielt++) {
    PDM_Mesh_nodal_elt_t type = block_elt_type[ielt];
    int _ielt = r_type_idx[type] + type_count[type]++;

    for (int i = 0; i < block_n_vtx_per_elt[ielt]; i++) {
      int ipt = (*r_vtx_idx)[_ielt] + i;
      for (int j = 0; j < 3; j++) {
        (*r_vtx_coord)[3*ipt + j] = block_vtx_coord[3*idx_vtx + j];
      }
      idx_vtx++;
    }

    for (int i = 0; i < block_n_pts_per_elt[ielt]; i++) {
      int ipt = (*r_pts_idx)[_ielt] + i;
      (*r_pts_g_num)[ipt] = block_pts_g_num[idx_pts];
      for (int j = 0; j < 3; j++) {
        (*r_pts_coord)[3*ipt + j] = block_pts_coord[3*idx_pts + j];
      }
      idx_pts++;
    }
  }

  free (block_n_vtx_per_elt);
  free (block_n_pts_per_elt);
  free (block_vtx_coord);
  free (block_pts_g_num);
  free (block_pts_coord);
  free (block_elt_type);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure to compute the location of point clouds inta a mesh
 *
 * \param [in]   mesh_nature    Nature of the mesh
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 *
 * \return     Identifier
 *
 */

int
PDM_mesh_location_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Comm comm
 )
{
  if (_locations == NULL) {
    _locations = PDM_Handles_create (4);
  }

  _PDM_location_t *location = (_PDM_location_t *) malloc(sizeof(_PDM_location_t));

  int id = PDM_Handles_store (_locations, location);

  location->n_point_cloud = n_point_cloud;
  location->comm = comm;
  location->mesh_nature = mesh_nature;

  location->shared_nodal   = 0;
  location->mesh_nodal_id  = -1;
  location->_mesh_nodal_id = -1;

  location->point_clouds =
    (_point_cloud_t*) malloc (sizeof(_point_cloud_t) * n_point_cloud);

  for (int i = 0; i <  n_point_cloud; i++) {
    location->point_clouds[i].n_part = -1;
    location->point_clouds[i].n_points = NULL;
    location->point_clouds[i].coords = NULL;
    location->point_clouds[i].gnum = NULL;
    location->point_clouds[i].location = NULL;
    location->point_clouds[i].uvw = NULL;
    location->point_clouds[i].weights = NULL;
    location->point_clouds[i].weights_idx = NULL;
  }

  location->tolerance = 0.;

  location->method = PDM_MESH_LOCATION_OCTREE;

  location->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    location->times_elapsed[i] = 0.;
    location->times_cpu[i]     = 0.;
    location->times_cpu_u[i]   = 0.;
    location->times_cpu_s[i]   = 0.;
  }

  return id;

}

void
PDM_mesh_location_create_cf
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Fint comm,
 int *id
 )
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(comm);

  *id = PDM_mesh_location_create (mesh_nature, n_point_cloud, _comm);

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
PDM_mesh_location_n_part_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_part
 )
{
  _PDM_location_t *location = _get_from_id (id);

  location->point_clouds[i_point_cloud].n_part = n_part;
  location->point_clouds[i_point_cloud].n_points =
    realloc(location->point_clouds[i_point_cloud].n_points, n_part * sizeof(int));
  location->point_clouds[i_point_cloud].coords =
    realloc(location->point_clouds[i_point_cloud].coords,
            n_part * sizeof(double *));
  location->point_clouds[i_point_cloud].gnum =
    realloc(location->point_clouds[i_point_cloud].gnum,
            n_part * sizeof(PDM_g_num_t *));

  for (int i = 0; i < n_part; i++) {
    location->point_clouds[i_point_cloud].n_points[i] = -1;
    location->point_clouds[i_point_cloud].coords[i] = NULL;
    location->point_clouds[i_point_cloud].gnum[i] = NULL;
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
PDM_mesh_location_cloud_set
(
 const int    id,
 const int    i_point_cloud,
 const int    i_part,
 const int    n_points,
 double      *coords,
 PDM_g_num_t *gnum
 )
{
  _PDM_location_t *location = _get_from_id (id);

  location->point_clouds[i_point_cloud].n_points[i_part] = n_points;
  location->point_clouds[i_point_cloud].coords[i_part] = coords;
  location->point_clouds[i_point_cloud].gnum[i_part] = gnum;

}


/**
 *
 * \brief Get a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [out]   n_points        Number of points
 * \param [out]   coords          Point coordinates
 * \param [out]   gnum            Point global number
 *
 */

void
PDM_mesh_location_cloud_get
(
 const int           id,
 const int           i_point_cloud,
 const int           i_part,
       int          *n_points,
       double      **coords,
       PDM_g_num_t **gnum
)
{
  _PDM_location_t *mesh_location = _get_from_id (id);

  assert (mesh_location->point_clouds != NULL);
  assert (i_point_cloud < mesh_location->n_point_cloud);

  _point_cloud_t *pcloud = mesh_location->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);

  *n_points        = pcloud->n_points[i_part];
  *coords          = pcloud->coords[i_part];
  *gnum            = pcloud->gnum[i_part];
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
PDM_mesh_location_shared_nodal_mesh_set
(
 const int  id,
 const int  mesh_nodal_id
 )
{
  _PDM_location_t *location = _get_from_id (id);

  location->mesh_nodal_id = mesh_nodal_id;
  location->shared_nodal = 1;
}


/**
 *
 * \brief Set global data of a mesh
 *
 * \param [in]   id             Identifier
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_mesh_location_mesh_global_data_set
(
 const int  id,
 const int  n_part
 )
{
  _PDM_location_t *location = _get_from_id (id);

  if ((location->shared_nodal == 0) && (location->mesh_nodal_id != -1)) {
    PDM_Mesh_nodal_free (location->mesh_nodal_id);
  }

  location->mesh_nodal_id = PDM_Mesh_nodal_create (n_part, location->comm);

  location->face_vtx_n  = malloc(sizeof(PDM_l_num_t *) * n_part);
  location->cell_face_n = malloc(sizeof(PDM_l_num_t *) * n_part);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    location->face_vtx_n [i_part] = NULL;
    location->cell_face_n[i_part] = NULL;
  }


}


/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_face_idx Index in the cell -> face connectivity
 * \param [in]   cell_face     cell -> face connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
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
PDM_mesh_location_part_set
(
 const int          id,
 const int          i_part,
 const int          n_cell,
 const int         *cell_face_idx,
 const int         *cell_face,
 const PDM_g_num_t *cell_ln_to_gn,
 const int          n_face,
 const int         *face_vtx_idx,
 const int         *face_vtx,
 const PDM_g_num_t *face_ln_to_gn,
 const int          n_vtx,
 const double      *coords,
 const PDM_g_num_t *vtx_ln_to_gn
 )
{
  _PDM_location_t *location = _get_from_id (id);

  PDM_UNUSED(face_ln_to_gn);

  /*
   * Creation de mesh nodal
   */

  PDM_Mesh_nodal_coord_set (location->mesh_nodal_id,
                            i_part,
                            n_vtx,
                            coords,
                            vtx_ln_to_gn);



  location->face_vtx_n[i_part]  = malloc (sizeof(PDM_l_num_t) * n_face);
  location->cell_face_n[i_part] = malloc (sizeof(PDM_l_num_t) * n_cell);

  for (int i = 0; i < n_face; i++) {
    location->face_vtx_n[i_part][i] = face_vtx_idx[i+1] - face_vtx_idx[i];
  }

  for (int i = 0; i < n_cell; i++) {
    location->cell_face_n[i_part][i] = cell_face_idx[i+1] - cell_face_idx[i];
  }

  PDM_Mesh_nodal_cell3d_cellface_add (location->mesh_nodal_id,
                                      i_part,
                                      n_cell,
                                      n_face,
                                      face_vtx_idx,
                                      location->face_vtx_n[i_part],
                                      face_vtx,
                                      cell_face_idx,
                                      location->cell_face_n[i_part],
                                      cell_face,
                                      cell_ln_to_gn);
}



/**
 *
 * \brief Set a part of a mesh (2d version)
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_edge_idx Index in the cell -> edge connectivity
 * \param [in]   cell_edge     cell -> edge connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
 * \param [in]   n_edge        Number of edges
 * \param [in]   edge_vtx_idx  Index in the edge -> vertex connectivity
 * \param [in]   edge_vtx      edge -> vertex connectivity
 * \param [in]   edge_ln_to_gn Local edge numbering to global edge numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */

void
PDM_mesh_location_part_set_2d
(
 const int          id,
 const int          i_part,
 const int          n_cell,
 const int         *cell_edge_idx,
 const int         *cell_edge,
 const PDM_g_num_t *cell_ln_to_gn,
 const int          n_edge,
 const int         *edge_vtx_idx,
 const int         *edge_vtx,
 const PDM_g_num_t *edge_ln_to_gn,
 const int          n_vtx,
 const double      *coords,
 const PDM_g_num_t *vtx_ln_to_gn
 )
{
  _PDM_location_t *location = _get_from_id (id);

  PDM_UNUSED (edge_ln_to_gn);

  /*
   * Creation de mesh nodal
   */

  PDM_Mesh_nodal_coord_set (location->mesh_nodal_id,
                            i_part,
                            n_vtx,
                            coords,
                            vtx_ln_to_gn);

  location->face_vtx_n[i_part]  = malloc (sizeof(PDM_l_num_t) * n_edge);
  location->cell_face_n[i_part] = malloc (sizeof(PDM_l_num_t) * n_cell);

  PDM_l_num_t *edge_vtx_nb  = location->face_vtx_n[i_part];
  PDM_l_num_t *cell_edge_nb = location->cell_face_n[i_part];

  for (int i = 0; i < n_edge; i++) {
    edge_vtx_nb[i] = edge_vtx_idx[i+1] - edge_vtx_idx[i];
  }

  for (int i = 0; i < n_cell; i++) {
    cell_edge_nb[i] = cell_edge_idx[i+1] - cell_edge_idx[i];
  }


  PDM_Mesh_nodal_cell2d_celledge_add (location->mesh_nodal_id,
                                      i_part,
                                      n_cell,
                                      n_edge,
                                      edge_vtx_idx,
                                      edge_vtx_nb,
                                      edge_vtx,
                                      cell_edge_idx,
                                      cell_edge_nb,
                                      cell_edge,
                                      cell_ln_to_gn);
}





/**
 *
 * \brief Set the tolerance for bounding boxes
 *
 * \param [in]   id              Identifier
 * \param [in]   tol             Tolerance
 *
 */

void
PDM_mesh_location_tolerance_set
(
 const int    id,
 const double tol
 )
{
  _PDM_location_t *location = _get_from_id (id);

  location->tolerance = tol;
}


/**
 *
 * \brief Set the method for computing location
 *
 * \param [in]   id              Identifier
 * \param [in]   method          Method
 *
 */

void
PDM_mesh_location_method_set
(
 const int                        id,
 const PDM_mesh_location_method_t method
 )
{
  _PDM_location_t *location = _get_from_id (id);

  location->method = method;
}




/**
 *
 * \brief Get mesh location
 *
 * \param [in]   id                    Identifier
 * \param [in]   i_point_cloud         Current cloud
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  n_points              Number of points in point cloud
 * \param [out]  coord                 Coordinates of points in point cloud
 * \param [out]  g_num                 Global numbers of points in point cloud
 * \param [out]  location              The global number of the closest element if the point is located,
 *                                     -1 otherwise
 *
 */

void
PDM_mesh_location_get
(
 const int     id,
 const int     i_point_cloud,
 const int     i_part,
 PDM_g_num_t **location,
 int         **weights_idx,
 double      **weights,
 double      **projected_coord
 )
{
  _PDM_location_t *mesh_location = _get_from_id (id);

  assert (mesh_location->point_clouds != NULL);
  assert (i_point_cloud < mesh_location->n_point_cloud);

  _point_cloud_t *pcloud = mesh_location->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);

  *location        = pcloud->location[i_part];
  *weights_idx     = pcloud->weights_idx[i_part];
  *weights         = pcloud->weights[i_part];
  *projected_coord = pcloud->projected_coords[i_part];
}




/**
 *
 * \brief Free a locationd mesh structure
 *
 * \param [in]  id       Identifier
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 */

void
PDM_mesh_location_free
(
 const int id,
 const int partial
 )
{
  _PDM_location_t *location = _get_from_id (id);

  /* Free point clouds */
  if (location->point_clouds != NULL) {
    for (int icloud = 0; icloud < location->n_point_cloud; icloud++) {
      _point_cloud_t *pcloud = location->point_clouds + icloud;

      if (pcloud->n_points != NULL) {
        free (pcloud->n_points);
      }

      if (pcloud->coords != NULL) {
        free (pcloud->coords);
      }

      if (pcloud->gnum != NULL) {
        free (pcloud->gnum);
      }

      if (pcloud->location != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (!partial) {
            if (pcloud->location[ipart] != NULL) {
              free (pcloud->location[ipart]);
            }
          }
        }
        free (pcloud->location);
      }

      if (pcloud->uvw != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (!partial) {
            if (pcloud->uvw[ipart] != NULL) {
              free (pcloud->uvw[ipart]);
            }
          }
        }
        free (pcloud->uvw);
      }

      if (pcloud->weights_idx != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (!partial) {
            if (pcloud->weights_idx[ipart] != NULL) {
              free (pcloud->weights_idx[ipart]);
            }
          }
        }
        free (pcloud->weights_idx);
      }

      if (pcloud->weights != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (!partial) {
            if (pcloud->weights[ipart] != NULL) {
              free (pcloud->weights[ipart]);
            }
          }
        }
        free (pcloud->weights);
      }

      if (pcloud->projected_coords != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (!partial) {
            if (pcloud->projected_coords[ipart] != NULL) {
              free (pcloud->projected_coords[ipart]);
            }
          }
        }
        free (pcloud->projected_coords);
      }

    }
    free (location->point_clouds);
  }

  /* Free mesh nodal */
  //PDM_Mesh_nodal_partial_free (location->mesh_nodal_id);?

  if (!location->shared_nodal) {

    int _n_part = PDM_Mesh_nodal_n_part_get(location->mesh_nodal_id);

    PDM_Mesh_nodal_free (location->mesh_nodal_id);

    if(location->cell_face_n != NULL){
      for (int i = 0; i< _n_part; i++) {
        if(location->cell_face_n[i] != NULL) {
          free(location->cell_face_n[i]);
        }
      }
      free (location->cell_face_n);
    }

    if(location->face_vtx_n != NULL){
      for (int i = 0; i< _n_part; i++) {
        if(location->face_vtx_n[i] != NULL) {
          free(location->face_vtx_n[i]);
        }
      }
      free (location->face_vtx_n);
    }
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
PDM_mesh_location_dump_times
(
 const int id
 )
{
  _PDM_location_t *location = _get_from_id (id);

  double t1 = location->times_elapsed[END] - location->times_elapsed[BEGIN];
  double t2 = location->times_cpu[END] - location->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, location->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, location->comm);

  double t_elaps_max[NTIMER];
  PDM_MPI_Allreduce (location->times_elapsed,
                     t_elaps_max,
                     NTIMER,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     location->comm);

  double t_cpu_max[NTIMER];
  PDM_MPI_Allreduce (location->times_cpu,
                     t_cpu_max, NTIMER,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     location->comm);

  int rank;
  PDM_MPI_Comm_rank (location->comm, &rank);

  if (rank == 0) {

    PDM_printf( "mesh_location timer : all (elapsed and cpu) :                                   "
                " %12.5es %12.5es\n",
                t1max, t2max);

    PDM_printf( "mesh_location timer : build bounding boxes (elapsed and cpu) :                  "
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_BOUNDING_BOXES],
                t_cpu_max[BUILD_BOUNDING_BOXES]);

    PDM_printf( "mesh_location timer : build aux. struct + search candidates (elapsed and cpu) : "
                " %12.5es %12.5es\n",
                t_elaps_max[SEARCH_CANDIDATES],
                t_cpu_max[SEARCH_CANDIDATES]);

    PDM_printf( "mesh_location timer : load balancing (elapsed and cpu) :                        "
                " %12.5es %12.5es\n",
                t_elaps_max[LOAD_BALANCING],
                t_cpu_max[LOAD_BALANCING]);

    PDM_printf( "mesh_location timer : compute elementary locations (elapsed and cpu) :          "
                " %12.5es %12.5es\n",
                t_elaps_max[COMPUTE_ELEMENTARY_LOCATIONS],
                t_cpu_max[COMPUTE_ELEMENTARY_LOCATIONS]);

    PDM_printf( "mesh_location timer : merge location data (elapsed and cpu) :                   "
                " %12.5es %12.5es\n",
                t_elaps_max[MERGE_LOCATION_DATA],
                t_cpu_max[MERGE_LOCATION_DATA]);
  }
}



/**
 *
 * \brief Compute point location
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_mesh_location_compute
(
 const int id
 )
{
  const double eps_dist = 1.e-10;
  const double tolerance = 1e-6;

  const int DEBUG = 0;
  const int dim = 3;

  const int octree_depth_max = 15;//?
  const int octree_points_in_leaf_max = 1;
  const int octree_build_leaf_neighbours = 0;
  int octree_id;

  _PDM_location_t *location = _get_from_id (id);

  int my_rank;
  PDM_MPI_Comm_rank (location->comm, &my_rank);

  int n_procs;
  PDM_MPI_Comm_size (location->comm, &n_procs);

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  location->times_elapsed[BEGIN] = PDM_timer_elapsed(location->timer);
  location->times_cpu[BEGIN]     = PDM_timer_cpu(location->timer);
  location->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(location->timer);
  location->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(location->timer);

  b_t_elapsed = location->times_elapsed[BEGIN];
  b_t_cpu     = location->times_cpu[BEGIN];
  b_t_cpu_u   = location->times_cpu_u[BEGIN];
  b_t_cpu_s   = location->times_cpu_s[BEGIN];
  PDM_timer_resume(location->timer);

  /*
   * Build the bounding boxes of mesh elements
   */
  int n_blocks = PDM_Mesh_nodal_n_blocks_get (location->mesh_nodal_id);
  int n_parts  = PDM_Mesh_nodal_n_part_get (location->mesh_nodal_id);
  int *blocks_id = PDM_Mesh_nodal_blocks_id_get (location->mesh_nodal_id);

  int n_boxes = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {
    n_boxes += PDM_Mesh_nodal_n_cell_get (location->mesh_nodal_id,
                                          ipart);
  }



  PDM_g_num_t *box_g_num   = malloc (sizeof(PDM_g_num_t) * n_boxes);
  double      *box_extents = malloc (sizeof(double)      * n_boxes * 6);

  int ibox = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {

    int id_block = blocks_id[iblock];

    for (int ipart = 0; ipart < n_parts; ipart++) {
      /* get element extents */
      PDM_Mesh_nodal_compute_cell_extents (location->mesh_nodal_id,
                                           id_block,
                                           ipart,
                                           location->tolerance,
                                           box_extents + 6*ibox);

      /* get elements gnum */
      PDM_g_num_t *_gnum = PDM_Mesh_nodal_g_num_get (location->mesh_nodal_id,
                                                     id_block,
                                                     ipart);

      int n_elt = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                  id_block,
                                                  ipart);

      for (int ielt = 0; ielt < n_elt; ielt++) {
        box_g_num[ibox] = _gnum[ielt];
        ibox++;
      }
    }
  }

  PDM_timer_hang_on(location->timer);
  e_t_elapsed = PDM_timer_elapsed(location->timer);
  e_t_cpu     = PDM_timer_cpu(location->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

  location->times_elapsed[BUILD_BOUNDING_BOXES] += e_t_elapsed - b_t_elapsed;
  location->times_cpu[BUILD_BOUNDING_BOXES]     += e_t_cpu - b_t_cpu;
  location->times_cpu_u[BUILD_BOUNDING_BOXES]   += e_t_cpu_u - b_t_cpu_u;
  location->times_cpu_s[BUILD_BOUNDING_BOXES]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;
  PDM_timer_resume(location->timer);




  PDM_dbbtree_t *dbbt = NULL;
  if (location->method == PDM_MESH_LOCATION_DBBTREE) {

    /* Compute local extents */
    double my_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    for (int i = 0; i < n_boxes; i++) {
      for (int j = 0; j < 3; j++) {
        my_extents[j]   = PDM_MIN (my_extents[j],   box_extents[6*i + j]);
        my_extents[j+3] = PDM_MAX (my_extents[j+3], box_extents[6*i + 3 + j]);
      }
    }

    /* Compute global extents */
    double global_extents[6];
    PDM_MPI_Allreduce (my_extents,   global_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, location->comm);
    PDM_MPI_Allreduce (my_extents+3, global_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, location->comm);

    /* Break symmetry */
    double max_range = 0.;
    for (int i = 0; i < 3; i++) {
      max_range = PDM_MAX (max_range, global_extents[i+3] - global_extents[i]);
    }
    for (int i = 0; i < 3; i++) {
      global_extents[i]   -= max_range * 1.1e-3;
      global_extents[i+3] += max_range * 1.0e-3;
    }

    dbbt = PDM_dbbtree_create (location->comm, dim, global_extents);

    PDM_dbbtree_boxes_set (dbbt,
                           1,//const int n_part,
                           &n_boxes,
                           (const double **) (&box_extents),
                           (const PDM_g_num_t **) (&box_g_num));

    PDM_timer_hang_on(location->timer);
    e_t_elapsed = PDM_timer_elapsed(location->timer);
    e_t_cpu     = PDM_timer_cpu(location->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

    location->times_elapsed[SEARCH_CANDIDATES] += e_t_elapsed - b_t_elapsed;
    location->times_cpu[SEARCH_CANDIDATES]     += e_t_cpu - b_t_cpu;
    location->times_cpu_u[SEARCH_CANDIDATES]   += e_t_cpu_u - b_t_cpu_u;
    location->times_cpu_s[SEARCH_CANDIDATES]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(location->timer);
  }


  /*
   * Locate points
   */
  int         *pts_idx   = NULL;
  PDM_g_num_t *pts_g_num = NULL;
  double      *pts_coord = NULL;

  for (int icloud = 0; icloud < location->n_point_cloud; icloud++) {

    PDM_timer_hang_on(location->timer);
    b_t_elapsed = PDM_timer_elapsed(location->timer);
    b_t_cpu     = PDM_timer_cpu(location->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(location->timer);
    PDM_timer_resume(location->timer);


    _point_cloud_t *pcloud = location->point_clouds + icloud;


    /*
     * Concatenate point cloud partitions
     */
    int n_pts_pcloud = 0;
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      n_pts_pcloud += pcloud->n_points[ipart];
    }
    PDM_g_num_t *pcloud_g_num = malloc (sizeof(PDM_g_num_t) * n_pts_pcloud);
    double      *pcloud_coord = malloc (sizeof(double)      * n_pts_pcloud * dim);
    int idx = 0;
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {

        pcloud_g_num[idx] = pcloud->gnum[ipart][ipt];

        for (int idim = 0; idim < dim; idim++) {
          pcloud_coord[dim*idx + idim] = pcloud->coords[ipart][dim*ipt + idim];
        }

        idx++;
      }
    }


    /*
     * Get points inside bounding boxes of elements
     */

    switch (location->method) {

    case PDM_MESH_LOCATION_OCTREE:
      /* Create octree structure */
      octree_id = PDM_para_octree_create (1,
                                          octree_depth_max,
                                          octree_points_in_leaf_max,
                                          octree_build_leaf_neighbours,
                                          location->comm);

      /* Set octree point cloud */
      PDM_para_octree_point_cloud_set (octree_id,
                                       0,
                                       n_pts_pcloud,
                                       pcloud_coord,
                                       pcloud_g_num);

      /* Build parallel octree */
      PDM_para_octree_build (octree_id, NULL);
      //PDM_para_octree_dump (octree_id);
      if (DEBUG) {
        PDM_para_octree_dump_times (octree_id);
      }


      /* Locate points inside boxes */
      PDM_para_octree_points_inside_boxes (octree_id,
                                           n_boxes,
                                           box_extents,
                                           box_g_num,
                                           &pts_idx,
                                           &pts_g_num,
                                           &pts_coord);

      /* Free octree */
      PDM_para_octree_free (octree_id);
      break;

    case PDM_MESH_LOCATION_DBBTREE:
      PDM_dbbtree_points_inside_boxes (dbbt,
                                       n_pts_pcloud,
                                       pcloud_g_num,
                                       pcloud_coord,
                                       n_boxes,
                                       box_g_num,
                                       &pts_idx,
                                       &pts_g_num,
                                       &pts_coord);
      break;

    default:
      printf("Error: unknown location method %d\n", location->method);
      assert (1 == 0);

    }
    free (pcloud_coord);


    if (0) {//DEBUG) {
      printf("\n[%d] --- Pts in box ---\n", my_rank);
      for (ibox = 0; ibox < n_boxes; ibox++) {

        if (pts_idx[ibox+1] <= pts_idx[ibox]) {
          continue;
        }

        printf("[%d] %d ("PDM_FMT_G_NUM"): ", my_rank, ibox, box_g_num[ibox]);
        for (int i = pts_idx[ibox]; i < pts_idx[ibox+1]; i++) {
          /*printf("((%ld); %f %f %f) ",
            pts_g_num[i], pts_coord[dim*i], pts_coord[dim*i+1], pts_coord[dim*i+2]);*/
          printf("("PDM_FMT_G_NUM") ", pts_g_num[i]);
        }
        printf("\n");
      }
      printf("[%d] ------------------\n\n\n", my_rank);
    }

    PDM_timer_hang_on(location->timer);
    e_t_elapsed = PDM_timer_elapsed(location->timer);
    e_t_cpu     = PDM_timer_cpu(location->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

    location->times_elapsed[SEARCH_CANDIDATES] += e_t_elapsed - b_t_elapsed;
    location->times_cpu[SEARCH_CANDIDATES]     += e_t_cpu - b_t_cpu;
    location->times_cpu_u[SEARCH_CANDIDATES]   += e_t_cpu_u - b_t_cpu_u;
    location->times_cpu_s[SEARCH_CANDIDATES]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(location->timer);




    /*
     * Load balancing: redistribute evenly elementary location operations
     */
    int redistrib_n_elt = 0;

    PDM_Mesh_nodal_elt_t *redistrib_elt_type  = NULL;
    PDM_g_num_t          *redistrib_elt_g_num = NULL;
    int                  *redistrib_vtx_idx   = NULL;
    double               *redistrib_vtx_coord = NULL;
    int                  *redistrib_pts_idx   = NULL;
    PDM_g_num_t          *redistrib_pts_g_num = NULL;
    double               *redistrib_pts_coord = NULL;

    PDM_l_num_t *redistrib_poly3d_face_idx    = NULL;
    PDM_l_num_t *redistrib_face_vtx_idx       = NULL;
    PDM_l_num_t *redistrib_face_vtx           = NULL;
    int         *redistrib_face_orientation   = NULL;

    int redistrib_type_idx[PDM_MESH_NODAL_N_ELEMENT_TYPES + 1];

    _redistribute_elementary_location (location,
                                       n_boxes,
                                       box_g_num,
                                       pts_idx,
                                       pts_g_num,
                                       pts_coord,
                                       &redistrib_n_elt,
                                       redistrib_type_idx,
                                       &redistrib_elt_g_num,
                                       &redistrib_vtx_idx,
                                       &redistrib_vtx_coord,
                                       &redistrib_pts_idx,
                                       &redistrib_pts_g_num,
                                       &redistrib_pts_coord,
                                       &redistrib_poly3d_face_idx,
                                       &redistrib_face_vtx_idx,
                                       &redistrib_face_vtx,
                                       &redistrib_face_orientation);
    free (pts_idx);
    free (pts_g_num);
    free (pts_coord);


    PDM_timer_hang_on(location->timer);
    e_t_elapsed = PDM_timer_elapsed(location->timer);
    e_t_cpu     = PDM_timer_cpu(location->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

    location->times_elapsed[LOAD_BALANCING] += e_t_elapsed - b_t_elapsed;
    location->times_cpu[LOAD_BALANCING]     += e_t_cpu - b_t_cpu;
    location->times_cpu_u[LOAD_BALANCING]   += e_t_cpu_u - b_t_cpu_u;
    location->times_cpu_s[LOAD_BALANCING]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(location->timer);


    /*
     * Location in elements
     */
    int n_pts = redistrib_pts_idx[redistrib_n_elt];
    PDM_g_num_t *redistrib_pts_location = malloc (sizeof(PDM_g_num_t) * n_pts);
    for (ibox = 0; ibox < redistrib_n_elt; ibox++) {
      for (int i = redistrib_pts_idx[ibox]; i < redistrib_pts_idx[ibox+1]; i++) {
        redistrib_pts_location[i] = redistrib_elt_g_num[ibox];
      }
    }
    free (redistrib_elt_g_num);
    free (redistrib_elt_type);

    PDM_Mesh_nodal_elt_t *redistrib_pts_elt_type = malloc (sizeof(PDM_Mesh_nodal_elt_t) * n_pts);
    for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_POINT;
         type < PDM_MESH_NODAL_N_ELEMENT_TYPES;
         type++) {
      for (ibox = redistrib_type_idx[type]; ibox < redistrib_type_idx[type+1]; ibox++) {
        for (int i = redistrib_pts_idx[ibox]; i < redistrib_pts_idx[ibox+1]; i++) {
          redistrib_pts_elt_type[i] = type;
        }
      }
    }

    double *distance        = NULL;
    double *projected_coord = NULL;
    int    *weights_idx     = NULL;
    double *weights         = NULL;

    PDM_point_location_nodal (redistrib_type_idx,
                              redistrib_vtx_idx,
                              redistrib_vtx_coord,
                              redistrib_poly3d_face_idx,
                              redistrib_face_vtx_idx,
                              redistrib_face_vtx,
                              redistrib_face_orientation,
                              redistrib_pts_idx,
                              redistrib_pts_coord,
                              tolerance,
                              &distance,
                              &projected_coord,
                              &weights_idx,
                              &weights);

    if (DEBUG) {
      for (int i = 0; i < n_pts; i++) {
        printf("Point gnum = ("PDM_FMT_G_NUM")\n", redistrib_pts_g_num[i]);
        printf("\t        coords = (%f, %f, %f)\n",
               redistrib_pts_coord[dim*i],
               redistrib_pts_coord[dim*i+1],
               redistrib_pts_coord[dim*i+2]);
        printf("\tlocation = ("PDM_FMT_G_NUM")\n", redistrib_pts_location[i]);
        printf("\t  proj. coords = (%f, %f, %f)\n",
               projected_coord[dim*i],
               projected_coord[dim*i+1],
               projected_coord[dim*i+2]);
        printf("\tdistance = %f\n", distance[i]);
        printf("\t weights =");
        for (int j = weights_idx[i]; j < weights_idx[i+1]; j++) {
          printf(" %f", (float) weights[j]);
        }
        printf("\n");
      }
    }
    free (redistrib_vtx_coord);
    free (redistrib_vtx_idx);
    free (redistrib_pts_coord);
    free (redistrib_pts_idx);
    free (redistrib_poly3d_face_idx);
    free (redistrib_face_vtx_idx);
    free (redistrib_face_vtx);
    free (redistrib_face_orientation);

    PDM_timer_hang_on(location->timer);
    e_t_elapsed = PDM_timer_elapsed(location->timer);
    e_t_cpu     = PDM_timer_cpu(location->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

    location->times_elapsed[COMPUTE_ELEMENTARY_LOCATIONS] += e_t_elapsed - b_t_elapsed;
    location->times_cpu[COMPUTE_ELEMENTARY_LOCATIONS]     += e_t_cpu - b_t_cpu;
    location->times_cpu_u[COMPUTE_ELEMENTARY_LOCATIONS]   += e_t_cpu_u - b_t_cpu_u;
    location->times_cpu_s[COMPUTE_ELEMENTARY_LOCATIONS]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(location->timer);



    /*
     * Merge location data
     */
    PDM_g_num_t *pcloud_location       = NULL;
    int         *pcloud_weights_stride = NULL;
    int         *pcloud_weights_idx    = NULL;
    double      *pcloud_weights        = NULL;
    double      *pcloud_proj_coord     = NULL;

    /*
     *   1) Part-to-block
     */
    PDM_part_to_block_t *ptb2 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                          &pcloud_g_num,
                                                          NULL,
                                                          &n_pts_pcloud,
                                                          1,
                                                          location->comm);

    PDM_g_num_t *block_g_num2 = PDM_part_to_block_block_gnum_get (ptb2);
    PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb2);

    PDM_part_to_block_t *ptb1 = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           &redistrib_pts_g_num,
                                                           block_distrib_idx,
                                                           &n_pts,
                                                           1,
                                                           location->comm);
    free (redistrib_pts_g_num);

    const int n_pts_block1 = PDM_part_to_block_n_elt_block_get (ptb1);

    int *part_stride = malloc (sizeof(int) * n_pts);
    for (int i = 0; i < n_pts; i++) {
      part_stride[i] = 1;
    }
    int *block_stride = NULL;

    /* Exchange location */
    PDM_g_num_t *block_location1 = NULL;
    PDM_part_to_block_exch (ptb1,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            1,
                            &part_stride,
                            (void **) &redistrib_pts_location,
                            &block_stride,
                            (void **) &block_location1);
    free (redistrib_pts_location);
    free (block_stride);

    /* Exchange element type */
    PDM_Mesh_nodal_elt_t *block_elt_type = NULL;
    PDM_part_to_block_exch (ptb1,
                            sizeof(PDM_Mesh_nodal_elt_t),
                            PDM_STRIDE_VAR,
                            1,
                            &part_stride,
                            (void **) &redistrib_pts_elt_type,
                            &block_stride,
                            (void **) &block_elt_type);
    free (redistrib_pts_elt_type);
    free (block_stride);

    /* Exchange distance */
    double *block_distance = NULL;
    PDM_part_to_block_exch (ptb1,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &part_stride,
                            (void **) &distance,
                            &block_stride,
                            (void **) &block_distance);
    free (distance);
    free (block_stride);


    /* Exchange weights */
    int *weights_stride = malloc (sizeof(int) * n_pts);
    for (int i = 0; i < n_pts; i++) {
      weights_stride[i] = weights_idx[i+1] - weights_idx[i];
    }
    free (weights_idx);

    double *block_weights1 = NULL;
    PDM_part_to_block_exch (ptb1,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &weights_stride,
                            (void **) &weights,
                            &block_stride,
                            (void **) &block_weights1);
    free (block_stride);
    free (weights);

    /* Exchange weights stride */
    int *block_n_vtx_elt    = NULL;
    int *block_n_candidates = NULL;
    PDM_part_to_block_exch (ptb1,
                            sizeof(int),
                            PDM_STRIDE_VAR,
                            1,
                            &part_stride,
                            (void **) &weights_stride,
                            &block_n_candidates,
                            (void **) &block_n_vtx_elt);
    free (weights_stride);

    /* Exchange projected coords */
    for (int i = 0; i < n_pts; i++) {
      part_stride[i] = 3;
    }
    double *block_proj_coord1 = NULL;
    PDM_part_to_block_exch (ptb1,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &part_stride,
                            (void **) &projected_coord,
                            &block_stride,
                            (void **) &block_proj_coord1);
    free (projected_coord);
    free (block_stride);
    free (part_stride);


    /*
     *   2) Among candidate elements, keep closest one for each point (set location to -1 if no candidate -> unlocated point)
     */
    PDM_g_num_t *block_g_num1 = PDM_part_to_block_block_gnum_get (ptb1);

    const int n_pts_block2 = PDM_part_to_block_n_elt_block_get (ptb2);

    PDM_g_num_t *block_location2 = malloc (sizeof(PDM_g_num_t) * n_pts_block2);
    double *block_proj_coord2 = malloc (sizeof(double) * n_pts_block2 * 3);
    int *block_weights_stride2 = malloc (sizeof(int) * n_pts_block2);
    for (int i = 0; i < n_pts_block2; i++) {
      block_location2[i] = -1;
      block_weights_stride2[i] = 0;
    }
    for (int i = 0; i < 3*n_pts_block2; i++) {
      block_proj_coord2[i] = 0.;
    }

    int *idx_min = malloc (sizeof(int) * n_pts_block2);
    idx = 0;
    int idw = 0;
    size_t s_weights = 0;
    for (int i = 0; i < n_pts_block1; i++) {
      idx_min[i] = idx;

      if (block_n_candidates[i] > 1) {
        double min_dist = HUGE_VAL;
        PDM_Mesh_nodal_elt_t type_min = PDM_MESH_NODAL_N_ELEMENT_TYPES;

        for (int j = idx; j < idx + block_n_candidates[i]; j++) {

          if (block_distance[j] < min_dist - eps_dist ||
              (block_distance[j] < min_dist + eps_dist &&
               type_min > block_elt_type[j])) {
            min_dist = block_distance[j];
            type_min = block_elt_type[j];
            idx_min[i] = j;
          }

          idw += block_n_vtx_elt[j];
        }
      }

      idx += block_n_candidates[i];

      int ipt = block_g_num1[i] - 1 - block_distrib_idx[my_rank];
      assert (block_g_num2[ipt] == block_g_num1[i]);

      block_location2[ipt] = block_location1[idx_min[i]];
      block_weights_stride2[ipt] = block_n_vtx_elt[idx_min[i]];
      for (int j = 0; j < 3; j++) {
        block_proj_coord2[3*ipt + j] = block_proj_coord1[3*idx_min[i] + j];
      }
      s_weights += block_weights_stride2[ipt];
    }

    int *block_weights_idx1 = malloc (sizeof(int) * (idx + 1));
    block_weights_idx1[0] = 0;
    for (int i = 0; i < idx; i++) {
      block_weights_idx1[i+1] = block_weights_idx1[i] + block_n_vtx_elt[i];
    }


    double *block_weights2 = malloc (sizeof(double) * s_weights);
    idx = 0;
    for (int i = 0; i < n_pts_block1; i++) {
      int ipt = block_g_num1[i] - 1 - block_distrib_idx[my_rank];

      for (int j = 0; j < block_weights_stride2[ipt]; j++) {
        block_weights2[idx++] = block_weights1[block_weights_idx1[idx_min[i]] + j];
      }
    }
    free (idx_min);
    free (block_n_candidates);
    free (block_n_vtx_elt);
    free (block_weights_idx1);
    free (block_distance);
    free (block_elt_type);
    free (block_location1);
    free (block_weights1);
    free (block_proj_coord1);

    /*
     *   3) Block-to-part
     */
    int one = 1;
    int three = 3;
    PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                         (const PDM_g_num_t **) &pcloud_g_num,
                                                         &n_pts_pcloud,
                                                         1,
                                                         location->comm);
    free (pcloud_g_num);

    /* Exchange location */
    pcloud_location = malloc (sizeof(PDM_g_num_t) * n_pts_pcloud);
    PDM_block_to_part_exch (btp,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST,
                            &one,
                            block_location2,
                            NULL,
                            (void **) &pcloud_location);
    free (block_location2);

    /* Exchange projected coords */
    pcloud_proj_coord = malloc (sizeof(double) * n_pts_pcloud * 3);
    PDM_block_to_part_exch (btp,
                            sizeof(double),
                            PDM_STRIDE_CST,
                            &three,
                            block_proj_coord2,
                            NULL,
                            (void **) &pcloud_proj_coord);
    free (block_proj_coord2);

    /* Exchange weights stride */
    pcloud_weights_stride = malloc (sizeof(int) * n_pts_pcloud);
    PDM_block_to_part_exch (btp,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            &one,
                            (void *) block_weights_stride2,
                            NULL,
                            (void **) &pcloud_weights_stride);

    pcloud_weights_idx = malloc (sizeof(int) * (n_pts_pcloud + 1));
    pcloud_weights_idx[0] = 0;
    for (int i = 0; i < n_pts_pcloud; i++) {
      pcloud_weights_idx[i+1] = pcloud_weights_idx[i] + pcloud_weights_stride[i];
    }

    /* Exchange weights */
    pcloud_weights = malloc (sizeof(double) * pcloud_weights_idx[n_pts_pcloud]);
    PDM_block_to_part_exch (btp,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            block_weights_stride2,
                            (void *) block_weights2,
                            &pcloud_weights_stride,
                            (void **) &pcloud_weights);
    free (block_weights2);
    free (block_weights_stride2);

    PDM_part_to_block_free (ptb1);
    PDM_part_to_block_free (ptb2);
    PDM_block_to_part_free (btp);


    /*
     * Conform to original partitioning of current point cloud
     */
    assert (pcloud->location == NULL);

    pcloud->location         = malloc (sizeof(PDM_g_num_t *) * pcloud->n_part);
    pcloud->weights_idx      = malloc (sizeof(int *)         * pcloud->n_part);
    pcloud->weights          = malloc (sizeof(double *)      * pcloud->n_part);
    pcloud->projected_coords = malloc (sizeof(double *)      * pcloud->n_part);

    idx = 0;
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      pcloud->location[ipart]         = malloc (sizeof(PDM_g_num_t) * pcloud->n_points[ipart]);
      pcloud->projected_coords[ipart] = malloc (sizeof(double) * pcloud->n_points[ipart] * 3);
      pcloud->weights_idx[ipart]      = malloc (sizeof(int) * (pcloud->n_points[ipart] + 1));

      pcloud->weights_idx[ipart][0] = 0;
      int idx_tmp = idx;
      for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {
        pcloud->location[ipart][ipt] = pcloud_location[idx];

        for (int j = 0; j < 3; j++) {
          pcloud->projected_coords[ipart][3*ipt + j] = pcloud_proj_coord[3*idx + j];
        }

        pcloud->weights_idx[ipart][ipt+1] = pcloud->weights_idx[ipart][ipt] + pcloud_weights_stride[idx];

        idx++;
      }

      pcloud->weights[ipart] = malloc (sizeof(double) * pcloud->weights_idx[ipart][pcloud->n_points[ipart]]);

      idx = idx_tmp;
      for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {
        for (int j = 0; j < pcloud_weights_stride[idx]; j++) {
          pcloud->weights[ipart][pcloud->weights_idx[ipart][ipt] + j] =
            pcloud_weights[pcloud_weights_idx[idx] + j];
        }
        idx++;
      }
    }
    free (pcloud_location);
    free (pcloud_weights_stride);
    free (pcloud_weights_idx);
    free (pcloud_weights);
    free (pcloud_proj_coord);


    PDM_timer_hang_on(location->timer);
    e_t_elapsed = PDM_timer_elapsed(location->timer);
    e_t_cpu     = PDM_timer_cpu(location->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

    location->times_elapsed[MERGE_LOCATION_DATA] += e_t_elapsed - b_t_elapsed;
    location->times_cpu[MERGE_LOCATION_DATA]     += e_t_cpu - b_t_cpu;
    location->times_cpu_u[MERGE_LOCATION_DATA]   += e_t_cpu_u - b_t_cpu_u;
    location->times_cpu_s[MERGE_LOCATION_DATA]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(location->timer);

  } // Loop over point clouds

  free (box_g_num);
  free (box_extents);


  if (dbbt != NULL) {
    PDM_dbbtree_free (dbbt);
  }


  PDM_timer_hang_on(location->timer);

  location->times_elapsed[END] = PDM_timer_elapsed(location->timer);
  location->times_cpu[END]     = PDM_timer_cpu(location->timer);
  location->times_cpu_u[END]   = PDM_timer_cpu_user(location->timer);
  location->times_cpu_s[END]   = PDM_timer_cpu_sys(location->timer);

  b_t_elapsed = location->times_elapsed[END];
  b_t_cpu     = location->times_cpu[END];
  b_t_cpu_u   = location->times_cpu_u[END];
  b_t_cpu_s   = location->times_cpu_s[END];
  PDM_timer_resume(location->timer);
}



int
PDM_mesh_location_mesh_nodal_id_get
(
 const int id
 )
{
  _PDM_location_t *location = _get_from_id (id);

  return location->mesh_nodal_id;
}

#ifdef	__cplusplus
}
#endif
