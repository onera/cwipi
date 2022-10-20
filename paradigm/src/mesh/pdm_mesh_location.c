/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

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
#include "pdm_dbbtree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_part.h"
#include "pdm_mesh_location.h"
#include "pdm_mesh_location_priv.h"
#include "pdm_point_location.h"
#include "pdm_ho_location.h"
#include "pdm_array.h"
#include "pdm_distrib.h"
#include "pdm_doctree.h"
// #include "pdm_mpi_priv.h"

#include "pdm_binary_search.h"
#include "pdm_para_octree.h"
#include "pdm_gnum.h"
#include "pdm_sort.h"
#include "pdm_logging.h"
#include "pdm_writer.h"
#include "pdm_vtk.h"
#include "pdm_extract_part.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_gnum_location.h"
#include "pdm_order.h"
#include "pdm_unique.h"

#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_mesh_nodal_priv.h"

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

//#define NTIMER_MESH_LOCATION 15

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \enum _ml_timer_step_t
 *
 */

typedef enum {

  BEGIN                        = 0,
  BUILD_BOUNDING_BOXES         = 1,
  STORE_CONNECTIVITY           = 2,
  SEARCH_CANDIDATES            = 3,
  LOAD_BALANCING               = 4,
  COMPUTE_ELEMENTARY_LOCATIONS = 5,
  MERGE_LOCATION_DATA_STEP1    = 6,
  MERGE_LOCATION_DATA_STEP2    = 7,
  MERGE_LOCATION_DATA_STEP3    = 8,
  COMPRESS_LOCATION_DATA       = 9,
  REVERSE_LOCATION_DATA        = 10,
  REVERSE_LOCATION_DATA_PTB    = 11,
  REVERSE_LOCATION_DATA_BTP1   = 12,
  REVERSE_LOCATION_DATA_BTP2   = 13,
  REVERSE_LOCATION_DATA_BTP3   = 14,
  REVERSE_LOCATION_DATA_UVW    = 15,
  END                          = 16

} _ml_timer_step_t;


/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

static
void
end_timer_and_print(const char* msg, PDM_MPI_Comm comm, double t1){

  double t2 = PDM_MPI_Wtime();

  double delta_t = t2 - t1;
  double delta_max;
  double delta_min;

  PDM_MPI_Allreduce (&delta_t,
                     &delta_max,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     comm);

  PDM_MPI_Allreduce (&delta_t,
                     &delta_min,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MIN,
                     comm);

  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);
  if(i_rank == 0) {
    printf("[%i] %s : duration min/max -> %12.5e %12.5e \n", n_rank, msg, delta_min, delta_max);
  }
}

struct _pdm_mpi_double_int_t {
  double val;
  int    rank;
};
typedef struct _pdm_mpi_double_int_t PDM_MPI_double_int_t;


// static
// void
// end_timer_and_log(const char* msg, PDM_MPI_Comm comm, double t1){

//   int n_rank;
//   int i_rank;

//   PDM_MPI_Comm_size(comm, &n_rank);
//   PDM_MPI_Comm_rank(comm, &i_rank);

//   double t2 = PDM_MPI_Wtime();

//   double delta_t = t2 - t1;

//   PDM_MPI_double_int_t l_info;

//   l_info.val  = delta_t;
//   l_info.rank = i_rank;

//   PDM_MPI_double_int_t g_max_info;
//   PDM_MPI_double_int_t g_min_info;


//   PDM_MPI_Allreduce (&l_info,
//                      &g_max_info,
//                      1,
//                      PDM_MPI_DOUBLE_INT,
//                      PDM_MPI_MAXLOC,
//                      comm);

//   PDM_MPI_Allreduce (&l_info,
//                      &g_min_info,
//                      1,
//                      PDM_MPI_DOUBLE_INT,
//                      PDM_MPI_MINLOC,
//                      comm);

//   log_trace("[%i] %s : duration min/max -> %12.5e [on rank = %i] %12.5e [on rank = %i] \n",
//            n_rank, msg, g_min_info.val, g_min_info.rank, g_max_info.val, g_max_info.rank);
//   if(i_rank == 0) {
//     printf("[%i] %s : duration min/max -> %12.5e [on rank = %i] %12.5e [on rank = %i] \n",
//            n_rank, msg, g_min_info.val, g_min_info.rank, g_max_info.val, g_max_info.rank);
//   }
// }


static
void
end_timer_and_log_from_dt(const char* msg, PDM_MPI_Comm comm, double delta_t){

  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_MPI_double_int_t l_info;

  l_info.val  = delta_t;
  l_info.rank = i_rank;

  PDM_MPI_double_int_t g_max_info;
  PDM_MPI_double_int_t g_min_info;


  PDM_MPI_Allreduce (&l_info,
                     &g_max_info,
                     1,
                     PDM_MPI_DOUBLE_INT,
                     PDM_MPI_MAXLOC,
                     comm);

  PDM_MPI_Allreduce (&l_info,
                     &g_min_info,
                     1,
                     PDM_MPI_DOUBLE_INT,
                     PDM_MPI_MINLOC,
                     comm);

  log_trace("[%i] %s : duration min/max -> %12.5e [on rank = %i] %12.5e [on rank = %i] \n",
           n_rank, msg, g_min_info.val, g_min_info.rank, g_max_info.val, g_max_info.rank);
  if(i_rank == 0) {
    printf("[%i] %s : duration min/max -> %12.5e [on rank = %i] %12.5e [on rank = %i] \n",
           n_rank, msg, g_min_info.val, g_min_info.rank, g_max_info.val, g_max_info.rank);
  }
}
/*
 *
 * Redistribute evenly across all ranks the elementary location operation to perform
 *
 */

static void
_redistribute_elementary_location
(
 PDM_mesh_location_t   *ml,
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

  int n_blocks   = 0;
  int n_parts    = 0;
  int *blocks_id = NULL;

  if (ml->mesh_nodal != NULL) {
    n_blocks   = PDM_Mesh_nodal_n_blocks_get (ml->mesh_nodal);
    n_parts    = PDM_Mesh_nodal_n_part_get   (ml->mesh_nodal);
    blocks_id = PDM_Mesh_nodal_blocks_id_get(ml->mesh_nodal);
  }

  int my_rank, n_ranks;
  PDM_MPI_Comm_rank (ml->comm, &my_rank);
  PDM_MPI_Comm_size (ml->comm, &n_ranks);

  /*
   * Get element type and number of points to locate per element
   */
  PDM_Mesh_nodal_elt_t *elt_type = malloc (sizeof(PDM_Mesh_nodal_elt_t) * n_elt);
  int *n_pts_per_elt = malloc (sizeof(int) * n_elt);
  int n_poly3d = 0;

  int ielt = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {
    int id_block = blocks_id[iblock];
    PDM_Mesh_nodal_elt_t block_type = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                     id_block);

    for (int ipart = 0; ipart < n_parts; ipart++) {
      int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
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

  /* Compute elements weights for an even redistribution */
  double *elt_weight = malloc (sizeof(double) * n_elt);
  for (ielt = 0; ielt < n_elt; ielt++) {
    elt_weight[ielt] = (double) n_pts_per_elt[ielt];
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &elt_g_num,
                                                       &elt_weight,
                                                       &n_elt,
                                                       1,
                                                       ml->comm);
  free (elt_weight);

  *r_n_elt = PDM_part_to_block_n_elt_block_get (ptb);
  PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);

  double **elt_g_weight = PDM_part_to_block_global_weight_get (ptb);
  double *_elt_g_weight = elt_g_weight[0];

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
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                         id_block,
                                                         ipart);

        PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (ml->mesh_nodal,
                                                          id_block,
                                                          ipart,
                                                          &connec_idx,
                                                          &connec);

        PDM_Mesh_nodal_block_poly3d_get (ml->mesh_nodal,
                                         id_block,
                                         ipart,
                                         &n_face,
                                         &face_vtx_idx,
                                         &face_vtx,
                                         &cell_face_idx,
                                         &cell_face);

        for (int i = 0; i < n_elt_part; i++) {
          poly3d_g_num[ipoly] = elt_g_num[ielt];
          if (_elt_g_weight[ielt] > 0) {
            n_vtx_per_elt[ielt++] = connec_idx[i+1] - connec_idx[i];
            n_face_per_elt[ipoly++] = cell_face_idx[i+1] - cell_face_idx[i];
          } else {
            n_vtx_per_elt[ielt++] = 0;
            n_face_per_elt[ipoly++] = 0;
          }
        }
      }
    }

    /* Polygons */
    else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                         id_block,
                                                         ipart);
        PDM_Mesh_nodal_block_poly2d_get (ml->mesh_nodal,
                                         id_block,
                                         ipart,
                                         &connec_idx,
                                         &connec);

        for (int i = 0; i < n_elt_part; i++) {
          if (_elt_g_weight[ielt] > 0) {
            n_vtx_per_elt[ielt++] = connec_idx[i+1] - connec_idx[i];
          } else {
            n_vtx_per_elt[ielt++] = 0;
          }
        }
      }
    }

    /* Standard elements */
    else {
      PDM_Mesh_nodal_elt_t std_type = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                     id_block);
      int n_vtx = PDM_Mesh_nodal_n_vertices_element (std_type,
                                                     order);

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                         id_block,
                                                         ipart);

        for (int i = 0; i < n_elt_part; i++) {
          if (_elt_g_weight[ielt] > 0) {
            n_vtx_per_elt[ielt++] = n_vtx;
          } else {
            n_vtx_per_elt[ielt++] = 0;
          }
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
          int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                           id_block,
                                                           ipart);

          PDM_Mesh_nodal_block_poly3d_get (ml->mesh_nodal,
                                           id_block,
                                           ipart,
                                           &n_face,
                                           &face_vtx_idx,
                                           &face_vtx,
                                           &cell_face_idx,
                                           &cell_face);

          PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (ml->mesh_nodal,
                                                            id_block,
                                                            ipart,
                                                            &connec_idx,
                                                            &connec);

          for (int i = 0; i < n_elt_part; i++) {
            if (n_face_per_elt[ipoly] == 0) {
              ipoly++;
              continue;
            }
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
  ielt = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {
    int id_block = blocks_id[iblock];

    /* Polyhedra */
    if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                         id_block,
                                                         ipart);

        const double *vtx_coord_part = PDM_Mesh_nodal_vertices_get (ml->mesh_nodal,
                                                                    ipart);

        PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (ml->mesh_nodal,
                                                          id_block,
                                                          ipart,
                                                          &connec_idx,
                                                          &connec);

        for (int i = 0; i < n_elt_part; i++) {
          if (n_vtx_per_elt[ielt++] > 0) {
            for (int j = connec_idx[i]; j < connec_idx[i+1]; j++) {
              int ivtx = connec[j] - 1;
              for (int k = 0; k < 3; k++) {
                _vtx_coord[k] = vtx_coord_part[3*ivtx + k];
              }
              _vtx_coord += 3;
            }
          }
        }
      } // End of loop on parts

    }

    /* Polygons */
    else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                         id_block,
                                                         ipart);

        const double *vtx_coord_part = PDM_Mesh_nodal_vertices_get (ml->mesh_nodal,
                                                                    ipart);

        PDM_Mesh_nodal_block_poly2d_get (ml->mesh_nodal,
                                         id_block,
                                         ipart,
                                         &connec_idx,
                                         &connec);

        for (int i = 0; i < n_elt_part; i++) {
          if (n_vtx_per_elt[ielt++] > 0) {
            for (int j = connec_idx[i]; j < connec_idx[i+1]; j++) {
              int ivtx = connec[j] - 1;
              for (int k = 0; k < 3; k++) {
                _vtx_coord[k] = vtx_coord_part[3*ivtx + k];
              }
              _vtx_coord += 3;
            }
          }
        }
      } // End of loop on parts

    }

    /* Standard elements */
    else {

      PDM_Mesh_nodal_elt_t std_type = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                     id_block);
      int n_vtx = PDM_Mesh_nodal_n_vertices_element (std_type,
                                                     order);


      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt_part = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                         id_block,
                                                         ipart);

        const double *vtx_coord_part = PDM_Mesh_nodal_vertices_get (ml->mesh_nodal,
                                                                    ipart);

        PDM_Mesh_nodal_block_std_get (ml->mesh_nodal,
                                      id_block,
                                      ipart,
                                      &connec);

        for (int i = 0; i < n_elt_part; i++) {
          if (n_vtx_per_elt[ielt++] > 0) {
            for (int j = 0; j < n_vtx; j++) {
              int ivtx = connec[n_vtx*i + j] - 1;
              for (int k = 0; k < 3; k++) {
                _vtx_coord[k] = vtx_coord_part[3*ivtx + k];
              }
              _vtx_coord += 3;
            }
          }
        }
      } // End of loop on parts

    }

  } // End of loop on nodal blocks



  /*
   * Exchange points to locate
   */

  /* Global number */
  int *block_n_pts_per_elt = NULL;
  PDM_g_num_t *block_pts_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &n_pts_per_elt,
                          (void **) &pts_g_num,
                          &block_n_pts_per_elt,
                          (void **) &block_pts_g_num);


  /* Coordinates */
  /*for (ielt = 0; ielt < n_elt; ielt++) {
    n_pts_per_elt[ielt] *= 3;
    }*/

  int *block_stride = NULL;
  double *block_pts_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          3*sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
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
  int *part_stride = PDM_array_const_int(n_elt, 1);

  /* Type */
  PDM_Mesh_nodal_elt_t *block_elt_type = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_Mesh_nodal_elt_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &part_stride,
                          (void **) &elt_type,
                          &block_stride,
                          (void **) &block_elt_type);
  free (block_stride);
  free (part_stride);
  free (elt_type);

  /* Global number */
  PDM_g_num_t *block_elt_g_num = PDM_part_to_block_block_gnum_get (ptb);

  /* Coordinates of vertices */
  /*for (ielt = 0; ielt < n_elt; ielt++) {
    n_vtx_per_elt[ielt] *= 3;
    }*/

  int *block_n_vtx_per_elt = NULL;
  double *block_vtx_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          3*sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &n_vtx_per_elt,
                          (void **) &vtx_coord,
                          &block_n_vtx_per_elt,
                          (void **) &block_vtx_coord);
  free (vtx_coord);

  /*for (ielt = 0; ielt < *r_n_elt; ielt++) {
    block_n_vtx_per_elt[ielt] /= 3;
    }*/

  /*
   * Polyhedra
   */

  PDM_part_to_block_t *ptb_poly3d = PDM_part_to_block_create_from_distrib (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                               PDM_PART_TO_BLOCK_POST_MERGE,
                                                               1.,
                                                               &poly3d_g_num,
                                                               block_distrib_idx,
                                                               &n_poly3d,
                                                               1,
                                                               ml->comm);

  int r_n_poly3d = PDM_part_to_block_n_elt_block_get (ptb_poly3d);

  /* Number of faces per polyhedron */
  part_stride = PDM_array_const_int(n_poly3d, 1);

  int *r_n_face_per_elt = NULL;
  PDM_part_to_block_exch (ptb_poly3d,
                          sizeof(int),
                          PDM_STRIDE_VAR_INTERLACED,
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
                          PDM_STRIDE_VAR_INTERLACED,
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
                          PDM_STRIDE_VAR_INTERLACED,
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
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &n_vtx_per_elt,
                          (void **) &local_face_vtx,
                          &block_stride,
                          (void **) r_face_vtx);
  free (block_stride);
  free (n_vtx_per_elt);
  free (local_face_vtx);

  free (poly3d_g_num);
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
  PDM_part_to_block_free (ptb);

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

  if (0) {
    printf("r_pts_idx[%d] = %d\n", *r_n_elt, (*r_pts_idx)[*r_n_elt]);
    printf("r_vtx_idx[%d] = %d\n", *r_n_elt, (*r_vtx_idx)[*r_n_elt]);
  }

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







static void
_extract_selected_mesh_elements
(
 PDM_mesh_location_t   *ml,
 int                  **n_select_elt,
 int                 ***select_elt_l_num,
 PDM_g_num_t           *select_elt_parent_g_num,
 PDM_g_num_t           *select_elt_g_num,
 const int              pts_idx[],
 const PDM_g_num_t      pts_g_num[],
 const double           pts_coord[],
 int                   *r_n_elt,
 int                    r_type_idx[],
 PDM_g_num_t          **r_elt_parent_g_num,
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

  int n_block   = 0;
  int n_part    = 0;
  int *block_id = NULL;

  if (ml->mesh_nodal != NULL) {
    n_block   = PDM_Mesh_nodal_n_blocks_get (ml->mesh_nodal);
    n_part    = PDM_Mesh_nodal_n_part_get   (ml->mesh_nodal);
    block_id = PDM_Mesh_nodal_blocks_id_get(ml->mesh_nodal);
  }

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (ml->comm, &i_rank);
  PDM_MPI_Comm_size (ml->comm, &n_rank);


  int n_elt = 0;
  for (int iblock = 0; iblock < n_block; iblock++) {
    for (int ipart = 0; ipart < n_part; ipart++) {
      n_elt += n_select_elt[iblock][ipart];
    }
  }

  /*
   * Get element type and number of points to locate per element
   */
  PDM_Mesh_nodal_elt_t *elt_type = malloc (sizeof(PDM_Mesh_nodal_elt_t) * n_elt);
  int *n_pts_per_elt = malloc (sizeof(int) * n_elt);
  int n_poly3d = 0;

  int ielt = 0;
  for (int iblock = 0; iblock < n_block; iblock++) {
    int id_block = block_id[iblock];
    PDM_Mesh_nodal_elt_t block_type = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                     id_block);

    for (int ipart = 0; ipart < n_part; ipart++) {

      if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        n_poly3d += n_select_elt[iblock][ipart];
      }

      for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
        elt_type[ielt] = block_type;
        n_pts_per_elt[ielt] = pts_idx[ielt+1] - pts_idx[ielt];
        ielt++;
      }
    } // End of loop on parts
  } // End of loop on nodal blocks


  /* Compute elements weights for an even redistribution */
  double *elt_weight = malloc (sizeof(double) * n_elt);
  for (ielt = 0; ielt < n_elt; ielt++) {
    elt_weight[ielt] = (double) n_pts_per_elt[ielt];
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       // PDM_PART_TO_BLOCK_POST_MERGE,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                       1.,
                                                       &select_elt_g_num,
                                                       &elt_weight,
                                                       &n_elt,
                                                       1,
                                                       ml->comm);
  free (elt_weight);

  *r_n_elt = PDM_part_to_block_n_elt_block_get (ptb);
  PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);

  double **elt_g_weight = PDM_part_to_block_global_weight_get (ptb);
  double *_elt_g_weight = elt_g_weight[0];

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
  for (int iblock = 0; iblock < n_block; iblock++) {
    int id_block = block_id[iblock];

    /* Polyhedra */
    if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

      for (int ipart = 0; ipart < n_part; ipart++) {
        PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (ml->mesh_nodal,
                                                          id_block,
                                                          ipart,
                                                          &connec_idx,
                                                          &connec);

        PDM_Mesh_nodal_block_poly3d_get (ml->mesh_nodal,
                                         id_block,
                                         ipart,
                                         &n_face,
                                         &face_vtx_idx,
                                         &face_vtx,
                                         &cell_face_idx,
                                         &cell_face);

        for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
          int l_num = select_elt_l_num[iblock][ipart][i];
          poly3d_g_num[ipoly] = select_elt_g_num[ielt];
          if (_elt_g_weight[ielt] > 0) {
            n_vtx_per_elt[ielt++] = connec_idx[l_num + 1] - connec_idx[l_num];
            n_face_per_elt[ipoly++] = cell_face_idx[l_num + 1] - cell_face_idx[l_num];
          } else {
            n_vtx_per_elt[ielt++] = 0;
            n_face_per_elt[ipoly++] = 0;
          }
        }
      }
    }

    /* Polygons */
    else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

      for (int ipart = 0; ipart < n_part; ipart++) {
        PDM_Mesh_nodal_block_poly2d_get (ml->mesh_nodal,
                                         id_block,
                                         ipart,
                                         &connec_idx,
                                         &connec);

        for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
          if (_elt_g_weight[ielt] > 0) {
            int l_num = select_elt_l_num[iblock][ipart][i];
            n_vtx_per_elt[ielt++] = connec_idx[l_num + 1] - connec_idx[l_num];
          } else {
            n_vtx_per_elt[ielt++] = 0;
          }
        }
      }
    }

    /* Standard elements */
    else {
      PDM_Mesh_nodal_elt_t std_type = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                     id_block);
      int n_vtx = PDM_Mesh_nodal_n_vertices_element (std_type,
                                                     order);

      for (int ipart = 0; ipart < n_part; ipart++) {
        for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
          if (_elt_g_weight[ielt] > 0) {
            n_vtx_per_elt[ielt++] = n_vtx;
          } else {
            n_vtx_per_elt[ielt++] = 0;
          }
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
    for (int iblock = 0; iblock < n_block; iblock++) {
      int id_block = block_id[iblock];

      /* Polyhedra */
      if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        for (int ipart = 0; ipart < n_part; ipart++) {
          PDM_Mesh_nodal_block_poly3d_get (ml->mesh_nodal,
                                           id_block,
                                           ipart,
                                           &n_face,
                                           &face_vtx_idx,
                                           &face_vtx,
                                           &cell_face_idx,
                                           &cell_face);

          PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (ml->mesh_nodal,
                                                            id_block,
                                                            ipart,
                                                            &connec_idx,
                                                            &connec);

          for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
            if (n_face_per_elt[ipoly] == 0) {
              ipoly++;
              continue;
            }
            int l_num = select_elt_l_num[iblock][ipart][i];
            int _n_vtx = connec_idx[l_num + 1] - connec_idx[l_num];
            int _n_face = cell_face_idx[l_num + 1] - cell_face_idx[l_num];

            for (int iface = 0; iface < _n_face; iface++) {
              int _iface = cell_face[cell_face_idx[l_num] + iface];

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
                                                   connec + connec_idx[l_num],
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
  ielt = 0;
  for (int iblock = 0; iblock < n_block; iblock++) {
    int id_block = block_id[iblock];

    /* Polyhedra */
    if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

      for (int ipart = 0; ipart < n_part; ipart++) {
        const double *vtx_coord_part = PDM_Mesh_nodal_vertices_get (ml->mesh_nodal,
                                                                    ipart);

        PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (ml->mesh_nodal,
                                                          id_block,
                                                          ipart,
                                                          &connec_idx,
                                                          &connec);

        for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
          if (n_vtx_per_elt[ielt++] > 0) {
            int l_num = select_elt_l_num[iblock][ipart][i];

            for (int j = connec_idx[l_num]; j < connec_idx[l_num+1]; j++) {
              int ivtx = connec[j] - 1;
              for (int k = 0; k < 3; k++) {
                _vtx_coord[k] = vtx_coord_part[3*ivtx + k];
              }
              _vtx_coord += 3;
            }
          }
        }
      } // End of loop on parts

    }

    /* Polygons */
    else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

      for (int ipart = 0; ipart < n_part; ipart++) {
        const double *vtx_coord_part = PDM_Mesh_nodal_vertices_get (ml->mesh_nodal,
                                                                    ipart);

        PDM_Mesh_nodal_block_poly2d_get (ml->mesh_nodal,
                                         id_block,
                                         ipart,
                                         &connec_idx,
                                         &connec);

        for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
          if (n_vtx_per_elt[ielt++] > 0) {
            int l_num = select_elt_l_num[iblock][ipart][i];

            for (int j = connec_idx[l_num]; j < connec_idx[l_num+1]; j++) {
              int ivtx = connec[j] - 1;
              for (int k = 0; k < 3; k++) {
                _vtx_coord[k] = vtx_coord_part[3*ivtx + k];
              }
              _vtx_coord += 3;
            }
          }
        }
      } // End of loop on parts

    }

    /* Standard elements */
    else {

      PDM_Mesh_nodal_elt_t std_type = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                     id_block);
      int n_vtx = PDM_Mesh_nodal_n_vertices_element (std_type,
                                                     order);


      for (int ipart = 0; ipart < n_part; ipart++) {
        const double *vtx_coord_part = PDM_Mesh_nodal_vertices_get (ml->mesh_nodal,
                                                                    ipart);

        PDM_Mesh_nodal_block_std_get (ml->mesh_nodal,
                                      id_block,
                                      ipart,
                                      &connec);

        for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
          if (n_vtx_per_elt[ielt++] > 0) {
            int l_num = select_elt_l_num[iblock][ipart][i];

            for (int j = 0; j < n_vtx; j++) {
              int ivtx = connec[n_vtx*l_num + j] - 1;
              for (int k = 0; k < 3; k++) {
                _vtx_coord[k] = vtx_coord_part[3*ivtx + k];
              }
              _vtx_coord += 3;
            }
          }
        }
      } // End of loop on parts

    }

  } // End of loop on nodal blocks



  /*
   * Exchange points to locate
   */

  /* Global number */
  int *block_n_pts_per_elt = NULL;
  PDM_g_num_t *block_pts_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &n_pts_per_elt,
                          (void **) &pts_g_num,
                          &block_n_pts_per_elt,
                          (void **) &block_pts_g_num);


  /* Coordinates */
  int *block_stride = NULL;
  double *block_pts_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          3*sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
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
  int *part_stride = PDM_array_const_int (n_elt, 1);

  /* Type */
  PDM_Mesh_nodal_elt_t *block_elt_type = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_Mesh_nodal_elt_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                          (void **) &elt_type,
                          NULL,
                          (void **) &block_elt_type);
  free (elt_type);

  /* Parent global number */
  PDM_g_num_t *block_elt_parent_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                          (void **) &select_elt_parent_g_num,
                          NULL,
                          (void **) &block_elt_parent_g_num);
  free (part_stride);

  /* Global number */
  PDM_g_num_t *block_elt_g_num = PDM_part_to_block_block_gnum_get (ptb);

  /* Coordinates of vertices */
  int *block_n_vtx_per_elt = NULL;
  double *block_vtx_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          3*sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &n_vtx_per_elt,
                          (void **) &vtx_coord,
                          &block_n_vtx_per_elt,
                          (void **) &block_vtx_coord);
  free (vtx_coord);

  /*
   * Polyhedra
   */

  PDM_part_to_block_t *ptb_poly3d = PDM_part_to_block_create_from_distrib (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                               PDM_PART_TO_BLOCK_POST_MERGE,
                                                               1.,
                                                               &poly3d_g_num,
                                                               block_distrib_idx,
                                                               &n_poly3d,
                                                               1,
                                                               ml->comm);

  int r_n_poly3d = PDM_part_to_block_n_elt_block_get (ptb_poly3d);

  /* Number of faces per polyhedron */
  part_stride = PDM_array_const_int(n_poly3d, 1);

  int *r_n_face_per_elt = NULL;
  PDM_part_to_block_exch (ptb_poly3d,
                          sizeof(int),
                          PDM_STRIDE_VAR_INTERLACED,
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
                          PDM_STRIDE_VAR_INTERLACED,
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
                          PDM_STRIDE_VAR_INTERLACED,
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
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &n_vtx_per_elt,
                          (void **) &local_face_vtx,
                          &block_stride,
                          (void **) r_face_vtx);
  free (block_stride);
  free (n_vtx_per_elt);
  free (local_face_vtx);

  free (poly3d_g_num);
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


  *r_elt_parent_g_num = malloc (sizeof(PDM_g_num_t) * (*r_n_elt));
  *r_elt_g_num = malloc (sizeof(PDM_g_num_t) * (*r_n_elt));
  *r_vtx_idx = malloc (sizeof(int) * (*r_n_elt + 1));
  (*r_vtx_idx)[0] = 0;
  *r_pts_idx = malloc (sizeof(int) * (*r_n_elt + 1));
  (*r_pts_idx)[0] = 0;

  for (ielt = 0; ielt < *r_n_elt; ielt++) {
    PDM_Mesh_nodal_elt_t type = block_elt_type[ielt];
    int _ielt = r_type_idx[type] + type_count[type]++;

    (*r_elt_parent_g_num)[_ielt] = block_elt_parent_g_num[ielt];
    (*r_elt_g_num)[_ielt] = block_elt_g_num[ielt];
    (*r_vtx_idx)[_ielt+1] = block_n_vtx_per_elt[ielt];
    (*r_pts_idx)[_ielt+1] = block_n_pts_per_elt[ielt];
  }
  free (block_elt_parent_g_num);
  PDM_part_to_block_free (ptb);

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

  if (0) {
    printf("r_pts_idx[%d] = %d\n", *r_n_elt, (*r_pts_idx)[*r_n_elt]);
    printf("r_vtx_idx[%d] = %d\n", *r_n_elt, (*r_vtx_idx)[*r_n_elt]);
  }

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

/**
 *
 * \brief Compute point location
 *
 * \param [in]   ml  Pointer to \ref PDM_mesh_location object
 *
 */

static void
_store_cell_vtx
(
  PDM_mesh_location_t  *ml
)
{

  int n_blocks = 0;
  int n_parts  = 0;
  int *blocks_id = NULL;

  if (ml->mesh_nodal != NULL) {
    n_blocks = PDM_Mesh_nodal_n_blocks_get (ml->mesh_nodal);
    n_parts  = PDM_Mesh_nodal_n_part_get (ml->mesh_nodal);
    blocks_id = PDM_Mesh_nodal_blocks_id_get (ml->mesh_nodal);
  }

  int **n_vtx_per_elt = malloc (sizeof(int *) * n_parts);

  for (int ipart = 0; ipart < n_parts; ipart++) {
    int n_elt = PDM_Mesh_nodal_n_cell_get (ml->mesh_nodal,
                                           ipart);

    n_vtx_per_elt[ipart] = malloc (sizeof(int) * n_elt);
    // ml->cell_vtx_idx[ipart] = malloc (sizeof(int) * (n_elt+1));
    // ml->cell_vtx_idx[ipart][0] = 0;
  }

  PDM_g_num_t n_g_cell = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {
    int n_elt = PDM_Mesh_nodal_n_cell_get (ml->mesh_nodal,
                                           ipart);
    int ielt = 0;
    for (int iblock = 0; iblock < n_blocks; iblock++) {
      int id_block = blocks_id[iblock];

      PDM_Mesh_nodal_elt_t t_elt = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                  id_block);

      int *parent_num = PDM_Mesh_nodal_block_parent_num_get (ml->mesh_nodal,
                                                             id_block,
                                                             ipart);

      PDM_g_num_t *g_num = PDM_Mesh_nodal_g_num_get (ml->mesh_nodal,
                                                      id_block,
                                                      ipart);

      int n_elt_block = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                        id_block,
                                                        ipart);

      for (int i = 0; i < n_elt_block; i++) {
        n_g_cell = PDM_MAX (n_g_cell, g_num[i]);
      }
      int n_vtx_elt = 0;
      switch (t_elt) {
      case PDM_MESH_NODAL_POINT:
        n_vtx_elt = 1;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_BAR2:
        n_vtx_elt = 2;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_TRIA3:
        n_vtx_elt = 3;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_QUAD4:
      case PDM_MESH_NODAL_TETRA4:
        n_vtx_elt = 4;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_PYRAMID5:
        n_vtx_elt = 5;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_PRISM6:
        n_vtx_elt = 6;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_HEXA8:
        n_vtx_elt = 8;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_POLY_2D:
        {
          int *connec_idx;
          int *connec;
          PDM_Mesh_nodal_block_poly2d_get (ml->mesh_nodal,
                                           id_block,
                                           ipart,
                                          &connec_idx,
                                          &connec);
          if (parent_num != NULL) {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][parent_num[i]] = connec_idx[i+1] - connec_idx[i];
            }
          }
          else {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][ielt++] = connec_idx[i+1] - connec_idx[i];
            }
          }
          break;
        }
      case PDM_MESH_NODAL_POLY_3D:
        {
          int *connec_idx;
          int *connec;
          PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (ml->mesh_nodal,
                                                            id_block,
                                                            ipart,
                                                            &connec_idx,
                                                            &connec);
          if (parent_num != NULL) {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][parent_num[i]] = connec_idx[i+1] - connec_idx[i];
            }
          }
          else {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][ielt++] = connec_idx[i+1] - connec_idx[i];
            }
          }
          break;
        }
      default :
        PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_location error : Bad Element type\n");
      }
    }
    
    ml->cell_vtx_idx[ipart] = malloc (sizeof(int) * (n_elt+1));
    ml->cell_vtx_idx[ipart][0] = 0;
    for (int i = 0; i < n_elt; i++) {
      ml->cell_vtx_idx[ipart][i+1] = ml->cell_vtx_idx[ipart][i] + n_vtx_per_elt[ipart][i];
    }
    ml->cell_vtx[ipart] = malloc (sizeof(int) * ml->cell_vtx_idx[ipart][n_elt]);
  }

  PDM_g_num_t _n_g_cell = 0;
  PDM_MPI_Allreduce (&n_g_cell, &_n_g_cell, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, ml->comm);
  n_g_cell = _n_g_cell;

  for (int ipart = 0; ipart < n_parts; ipart++) {
    int ielt = 0;
    for (int iblock = 0; iblock < n_blocks; iblock++) {
      int id_block = blocks_id[iblock];

      PDM_Mesh_nodal_elt_t t_elt = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                  id_block);

      int *parent_num = PDM_Mesh_nodal_block_parent_num_get (ml->mesh_nodal,
                                                             id_block,
                                                             ipart);

      int n_elt_block = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                        id_block,
                                                        ipart);

      int n_vtx_elt = 0;
      switch (t_elt) {
      case PDM_MESH_NODAL_POINT:
      case PDM_MESH_NODAL_BAR2:
      case PDM_MESH_NODAL_TRIA3:
      case PDM_MESH_NODAL_QUAD4:
      case PDM_MESH_NODAL_TETRA4:
      case PDM_MESH_NODAL_PYRAMID5:
      case PDM_MESH_NODAL_PRISM6:
      case PDM_MESH_NODAL_HEXA8: {
        int *connec;
        PDM_Mesh_nodal_block_std_get (ml->mesh_nodal,
                                      id_block,
                                      ipart,
                                      &connec);
        if (parent_num != NULL) {
          int idx2 = 0;
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][parent_num[i]];
            int idx = ml->cell_vtx_idx[ipart][parent_num[i]];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx2++];
            }
          }
        }
        else {
          int idx2 = 0;
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][ielt];
            int idx = ml->cell_vtx_idx[ipart][ielt++];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx2++];
            }
          }
        }
        break;
      }
      case PDM_MESH_NODAL_POLY_2D: {
        int *connec_idx;
        int *connec;
        PDM_Mesh_nodal_block_poly2d_get (ml->mesh_nodal,
                                         id_block,
                                         ipart,
                                        &connec_idx,
                                        &connec);
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][parent_num[i]];
            int idx = ml->cell_vtx_idx[ipart][parent_num[i]];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][ielt];
            int idx = ml->cell_vtx_idx[ipart][ielt++];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        break;
      }
      case PDM_MESH_NODAL_POLY_3D:{
        int *connec_idx;
        int *connec;
        PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (ml->mesh_nodal,
                                                          id_block,
                                                          ipart,
                                                          &connec_idx,
                                                          &connec);
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][parent_num[i]];
            int idx = ml->cell_vtx_idx[ipart][parent_num[i]];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][ielt];
            int idx = ml->cell_vtx_idx[ipart][ielt++];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        break;
      }
      default :
        PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_location error : Bad Element type\n");
      }
    }

    // int n_elt = PDM_Mesh_nodal_n_cell_get (ml->mesh_nodal,
    //                                        ipart);
    // PDM_log_trace_array_int(n_vtx_per_elt[ipart],
    //                         n_elt,
    //                         "n_vtx_per_elt : ");
    // PDM_log_trace_connectivity_int(ml->cell_vtx_idx[ipart],
    //                                ml->cell_vtx[ipart],
    //                                n_elt,
    //                                "ml->cell_vtx : ");
  }

  for (int ipart = 0; ipart < n_parts; ipart++) {
    free (n_vtx_per_elt[ipart]);
  }
  free (n_vtx_per_elt);
}



static void _pmesh_nodal_elmts_free
(
 PDM_part_mesh_nodal_elmts_t *pmne,
 int                          hack
 )
{
  free(pmne->sections_id);
  if (hack) {
    for (int i = 0; i < pmne->n_section_std; i++) {
      free(pmne->sections_std[i]);
    }
  }
  free(pmne->sections_std   );
  free(pmne->sections_poly2d);
  free(pmne->sections_poly3d);
  free(pmne->n_elmts);
  free(pmne);
}

static PDM_part_mesh_nodal_elmts_t *
_mesh_nodal_to_pmesh_nodal_elmts
(
 PDM_Mesh_nodal_t *mesh_nodal,
 PDM_MPI_Comm     location_comm
)
{
  PDM_part_mesh_nodal_elmts_t *pmne = NULL;

  int  n_block   = 0;
  int  n_part    = 0;
  int *blocks_id = NULL;

  if (mesh_nodal != NULL) {
    n_block   = PDM_Mesh_nodal_n_blocks_get (mesh_nodal);
    n_part    = PDM_Mesh_nodal_n_part_get   (mesh_nodal);
    blocks_id = PDM_Mesh_nodal_blocks_id_get(mesh_nodal);
  }

  /* Infer mesh dimension from mesh_nodal */
  int has_dim[4] = {0, 0, 0, 0};
  for (int iblock = 0; iblock < n_block; iblock++) {
    int id_block = blocks_id[iblock];

    PDM_Mesh_nodal_elt_t t_elt = PDM_Mesh_nodal_block_type_get(mesh_nodal,
                                                               id_block);

    int dim = PDM_Mesh_nodal_elt_dim_get(t_elt);
    has_dim[dim] = 1;
  }
  // PDM_log_trace_array_int(has_dim, 4, "has_dim : ");

  int mesh_dimension = 0;
  for (int i = 3; i >= 0; i--) {
    if (has_dim[i] != 0) {
      for (int j = 0; j < i; j++) {
        assert(has_dim[j] == 0);
      }

      mesh_dimension = i;
    }
  }


  int max_mesh_dimension = 0;
  PDM_MPI_Allreduce (&mesh_dimension, &max_mesh_dimension, 1, PDM_MPI_INT, PDM_MPI_MAX, location_comm);

  /* Gather block types (check coherence and add info for procs without mesh_nodal) */
  int max_n_block = 0;
  PDM_MPI_Allreduce (&n_block, &max_n_block, 1, PDM_MPI_INT, PDM_MPI_MAX, location_comm);

  PDM_g_num_t *section_gnum = malloc(sizeof(PDM_g_num_t) * max_n_block);
  int         *section_id   = malloc(sizeof(int        ) * max_n_block);
  int         *section_type = malloc(sizeof(int        ) * max_n_block);
  if (mesh_nodal == NULL) {
    for (int iblock = 0; iblock < max_n_block; iblock++) {
      section_gnum[iblock] = iblock + 1;
      section_id  [iblock] = -1;
      section_type[iblock] = -1;
    }
  }
  else {
    assert(n_block == max_n_block);

    for (int iblock = 0; iblock < n_block; iblock++) {
      int id_block = blocks_id[iblock];
      section_gnum[iblock] = iblock + 1;
      section_id  [iblock] = id_block;
      section_type[iblock] = (int) PDM_Mesh_nodal_block_type_get(mesh_nodal,
                                                                 id_block);
    }
  }

  // log_trace("before ptb\n");
  // PDM_log_trace_array_int(section_id,   max_n_block, "section_id   : ");
  // PDM_log_trace_array_int(section_type, max_n_block, "section_type : ");

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &section_gnum,
                                                      NULL,
                                                      &max_n_block,
                                                      1,
                                                      location_comm);
  free(section_gnum);

  int *part_stride  = PDM_array_const_int(max_n_block, 1);
  int *block_stride = NULL;
  int *block_section_type = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &part_stride,
               (void **) &section_type,
                         &block_stride,
               (void  *) &block_section_type);
  free(section_type);
  free(block_stride);

  int *block_section_id = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &part_stride,
               (void **) &section_id,
                         &block_stride,
               (void  *) &block_section_id);
  free(section_id);
  free(part_stride);

  /* Check coherence (merge) */
  int block_n_section = PDM_part_to_block_n_elt_block_get(ptb);
  int idx_read = 0;
  for (int isection = 0; isection < block_n_section; isection++) {
    assert(block_stride[isection] > 0);
    // log_trace("isection = %d:\n", isection);
    // PDM_log_trace_array_int(block_section_type + idx_read, block_stride[isection], "  type : ");
    // PDM_log_trace_array_int(block_section_id   + idx_read, block_stride[isection], "  id   : ");
    block_section_type[isection] = block_section_type[idx_read];
    block_section_id  [isection] = block_section_id  [idx_read++];
    for (int i = 1; i < block_stride[isection]; i++) {
      int type = block_section_type[idx_read];
      int id   = block_section_id  [idx_read++];
      if (type != -1) {
        if (block_section_type[isection] == -1) {
          block_section_type[isection] = type;
        }
        else {
          assert(type == block_section_type[isection]);
        }
      }

      if (id != -1) {
        if (block_section_id[isection] == -1) {
          block_section_id[isection] = id;
        }
        else {
          assert(id == block_section_id[isection]);
        }
      }

      assert(block_section_type[isection] != -1);
      assert(block_section_id  [isection] != -1);
    }
  }
  free(block_stride);
  block_section_type = realloc(block_section_type, sizeof(int) * block_n_section);
  block_section_id   = realloc(block_section_id,   sizeof(int) * block_n_section);

  int **psection_type = NULL;
  PDM_part_to_block_reverse_exch(ptb,
                                 sizeof(int),
                                 PDM_STRIDE_CST_INTERLACED,
                                 1,
                                 NULL,
                      (void   *) block_section_type,
                                 NULL,
                      (void ***) &psection_type);
  free(block_section_type);


  int **psection_id = NULL;
  PDM_part_to_block_reverse_exch(ptb,
                                 sizeof(int),
                                 PDM_STRIDE_CST_INTERLACED,
                                 1,
                                 NULL,
                      (void   *) block_section_id,
                                 NULL,
                      (void ***) &psection_id);
  free(block_section_id);

  PDM_part_to_block_free(ptb);

  section_type = psection_type[0];
  free(psection_type);

  section_id = psection_id[0];
  free(psection_id);

  // log_trace("after ptb\n");
  // PDM_log_trace_array_int(section_id,   max_n_block, "section_id   : ");
  // PDM_log_trace_array_int(section_type, max_n_block, "section_type : ");


  // log_trace("mesh_dimension = %d", mesh_dimension);

  pmne = PDM_part_mesh_nodal_elmts_create(max_mesh_dimension,
                                          n_part,
                                          location_comm);

  pmne->n_section        = max_n_block;//n_block;
  pmne->n_section_std    = 0;
  pmne->n_section_poly2d = 0;
  pmne->n_section_poly3d = 0;

  pmne->sections_id = section_id;//blocks_id; // Legit???

  for (int iblock = 0; iblock < max_n_block; iblock++) {

    // int id_block = blocks_id[iblock];
    int id_block = section_id[iblock];

    if (id_block < PDM_BLOCK_ID_BLOCK_POLY2D) {
      pmne->n_section_std++;
    }
    else if (id_block < PDM_BLOCK_ID_BLOCK_POLY3D) {
      id_block -= PDM_BLOCK_ID_BLOCK_POLY2D;
      pmne->n_section_poly2d++;
    }
    else {
      id_block -= PDM_BLOCK_ID_BLOCK_POLY3D;
      pmne->n_section_poly3d++;
    }

  }


  pmne->sections_std    = malloc(sizeof(PDM_Mesh_nodal_block_std_t    *) * pmne->n_section_std   );
  pmne->sections_poly2d = malloc(sizeof(PDM_Mesh_nodal_block_poly2d_t *) * pmne->n_section_poly2d);
  pmne->sections_poly3d = malloc(sizeof(PDM_Mesh_nodal_block_poly3d_t *) * pmne->n_section_poly3d);

  pmne->n_section_std    = 0;
  pmne->n_section_poly2d = 0;
  pmne->n_section_poly3d = 0;

  if (mesh_nodal == NULL) {

    // Add empty nodal blocks with appropriate types
    int *order = malloc(sizeof(int) * max_n_block);
    for (int i = 0; i < max_n_block; i++) {
      order[i] = i;
    }
    PDM_sort_int(section_id, order, max_n_block);

    for (int isection = 0; isection < max_n_block; isection++) {
      PDM_Mesh_nodal_elt_t type = (PDM_Mesh_nodal_elt_t) section_type[order[isection]];
      // int id = section_id[isection];

      if (type == PDM_MESH_NODAL_POLY_2D) {
        pmne->sections_poly2d[pmne->n_section_poly2d++] = NULL;
      }
      else if (type == PDM_MESH_NODAL_POLY_3D) {
        pmne->sections_poly3d[pmne->n_section_poly3d++] = NULL;
      }
      else {
        pmne->sections_std[pmne->n_section_std++] = malloc(sizeof(PDM_Mesh_nodal_block_std_t));
        PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[pmne->n_section_std-1];
        block->t_elt = type;
      }

    }


    for (int i = 0; i < n_part; i++) {
      pmne->n_elmts[i] = 0;
    }
    free(order);
  }

  else {

    for (int iblock = 0; iblock < n_block; iblock++) {

      int id_block = blocks_id[iblock];

      if (id_block < PDM_BLOCK_ID_BLOCK_POLY2D) {
        pmne->sections_std[pmne->n_section_std++] = mesh_nodal->blocks_std[id_block];
      }
      else if (id_block < PDM_BLOCK_ID_BLOCK_POLY3D) {
        id_block -= PDM_BLOCK_ID_BLOCK_POLY2D;
        pmne->sections_poly2d[pmne->n_section_poly2d++] = mesh_nodal->blocks_poly2d[id_block];
      }
      else {
        id_block -= PDM_BLOCK_ID_BLOCK_POLY3D;
        pmne->sections_poly3d[pmne->n_section_poly3d++] = mesh_nodal->blocks_poly3d[id_block];
      }

    }

    for (int i = 0; i < n_part; i++) {
      pmne->n_elmts[i] = PDM_Mesh_nodal_n_cell_get(mesh_nodal,
                                                   i);
    }
  }

  free(section_type);


  return pmne;
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
 * \param [in]   owner          OwnerShip
 *
 * \return     Pointer to \ref PDM_mesh_location object
 *
 */

PDM_mesh_location_t*
PDM_mesh_location_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int               n_point_cloud,
 const PDM_MPI_Comm      comm,
 const PDM_ownership_t   owner
)
{

  PDM_mesh_location_t *ml = (PDM_mesh_location_t *) malloc(sizeof(PDM_mesh_location_t));

  ml->n_point_cloud = n_point_cloud;
  ml->comm = comm;
  ml->mesh_nature = mesh_nature;

  ml->shared_nodal = 0;
  ml->mesh_nodal   = NULL;
  ml->_mesh_nodal  = NULL;

  ml->use_user_extract = 0;
  ml->is_elmt_select_by_user = NULL;

  ml->point_clouds =
    (_point_cloud_t*) malloc (sizeof(_point_cloud_t) * n_point_cloud);

  ml->ptp = malloc(sizeof(PDM_part_to_part_t *) * n_point_cloud);

  ml->ptp_ownership = malloc(sizeof(PDM_ownership_t) * n_point_cloud);

  for (int i = 0; i <  n_point_cloud; i++) {
    ml->point_clouds[i].n_part = -1;
    ml->point_clouds[i].n_points = NULL;
    ml->point_clouds[i].coords = NULL;
    ml->point_clouds[i].gnum = NULL;
    ml->point_clouds[i].location = NULL;
    ml->point_clouds[i].uvw = NULL;
    ml->point_clouds[i].weights = NULL;
    ml->point_clouds[i].weights_idx = NULL;
    ml->point_clouds[i].projected_coords = NULL;
    ml->point_clouds[i].n_located = NULL;
    ml->point_clouds[i].n_un_located = NULL;
    ml->point_clouds[i].located = NULL;
    ml->point_clouds[i].un_located = NULL;

    ml->ptp[i] = NULL;
    ml->ptp_ownership[i] = PDM_OWNERSHIP_KEEP;
  }

  ml->points_in_elements = NULL;

  ml->reverse_result = 0;
  ml->tolerance = 0.;

  ml->method = PDM_MESH_LOCATION_OCTREE;
  // ml->method = PDM_MESH_LOCATION_DBBTREE;

  ml->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER_MESH_LOCATION; i++) {
    ml->times_elapsed[i] = 0.;
    ml->times_cpu[i]     = 0.;
    ml->times_cpu_u[i]   = 0.;
    ml->times_cpu_s[i]   = 0.;
  }

  ml->reverse_result = 0;

  ml->owner = owner;

  ml->tag_unlocated_get = 0;
  ml->tag_located_get = 0;
  ml->tag_point_location_get = 0;
  ml->tag_points_in_elt_get = 0;
  ml->tag_cell_vtx_get = 0;

  return ml;

}



/**
 *
 * \brief Get the number of located points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The number of located points
 *
 */

int
PDM_mesh_location_n_located_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);
  return pcloud->n_located[i_part];
}


/**
 *
 * \brief Get the number of unlocated points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The number of unlocated points
 *
 */
int
PDM_mesh_location_n_unlocated_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);
  return pcloud->n_un_located[i_part];
}


/**
 *
 * \brief Get the list of unlocated points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The list of unlocated points
 *
 */
int *
PDM_mesh_location_unlocated_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);
  ml->tag_unlocated_get = 1;
  return pcloud->un_located[i_part];
}


/**
 *
 * \brief Get the list of located points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The list of located points
 *
 */
int *
PDM_mesh_location_located_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);
  ml->tag_located_get = 1;
  return pcloud->located[i_part];
}

/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_mesh_location_n_part_cloud_set
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  n_part
)
{

  ml->point_clouds[i_point_cloud].n_part = n_part;
  ml->point_clouds[i_point_cloud].n_points =
    realloc(ml->point_clouds[i_point_cloud].n_points, n_part * sizeof(int));
  ml->point_clouds[i_point_cloud].coords =
    realloc(ml->point_clouds[i_point_cloud].coords,
            n_part * sizeof(double *));
  ml->point_clouds[i_point_cloud].gnum =
    realloc(ml->point_clouds[i_point_cloud].gnum,
            n_part * sizeof(PDM_g_num_t *));

  for (int i = 0; i < n_part; i++) {
    ml->point_clouds[i_point_cloud].n_points[i] = -1;
    ml->point_clouds[i_point_cloud].coords[i] = NULL;
    ml->point_clouds[i_point_cloud].gnum[i] = NULL;
  }

}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
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
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part,
 const int                  n_points,
       double              *coords,
       PDM_g_num_t         *gnum
)
{

  ml->point_clouds[i_point_cloud].n_points[i_part] = n_points;
  ml->point_clouds[i_point_cloud].coords[i_part] = coords;
  ml->point_clouds[i_point_cloud].gnum[i_part] = gnum;

}


/**
 *
 * \brief Get a point cloud
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [out]  n_points        Number of points
 * \param [out]  coords          Point coordinates
 * \param [out]  gnum            Point global number
 *
 */
void
PDM_mesh_location_cloud_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_point_cloud,
 const int                   i_part,
       int                  *n_points,
       double              **coords,
       PDM_g_num_t         **gnum
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);

  *n_points        = pcloud->n_points[i_part];
  *coords          = pcloud->coords[i_part];
  *gnum            = pcloud->gnum[i_part];
}



/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   id             Pointer to \ref PDM_mesh_location object
 * \param [in]   mesh_nodal  Mesh nodal Pointer to \ref PDM_mesh_location object
 *
 */

void
PDM_mesh_location_shared_nodal_mesh_set
(
 PDM_mesh_location_t *ml,
 PDM_Mesh_nodal_t    *mesh_nodal
)
{
  ml->mesh_nodal = mesh_nodal;
  ml->shared_nodal = 1;

  int n_part = 0;
  if (mesh_nodal != NULL) {
    n_part = PDM_Mesh_nodal_n_part_get(mesh_nodal);
  }

  ml->cell_vtx_idx = malloc(sizeof(PDM_l_num_t *) * n_part);
  ml->cell_vtx     = malloc(sizeof(PDM_l_num_t *) * n_part);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    ml->cell_vtx_idx[i_part] = NULL;
    ml->cell_vtx    [i_part] = NULL;
  }

}


/**
 *
 * \brief Set global data of a mesh
 *
 * \param [in]   id             Pointer to \ref PDM_mesh_location object
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_mesh_location_mesh_global_data_set
(
       PDM_mesh_location_t *ml,
 const int                  n_part
)
{

  assert (ml->shared_nodal == 0);

  if (ml->shared_nodal == 0) {
    if(ml->mesh_nodal != NULL) {
      PDM_Mesh_nodal_free (ml->mesh_nodal);
    }
    ml->mesh_nodal = PDM_Mesh_nodal_create (n_part, ml->comm);
  }

  if(ml->shared_nodal == 0) {
    ml->face_vtx_n   = malloc(sizeof(PDM_l_num_t *) * n_part);
    ml->cell_face_n  = malloc(sizeof(PDM_l_num_t *) * n_part);
  }

  ml->cell_vtx_idx = malloc(sizeof(PDM_l_num_t *) * n_part);
  ml->cell_vtx     = malloc(sizeof(PDM_l_num_t *) * n_part);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    if(ml->shared_nodal == 0) {
      ml->face_vtx_n  [i_part] = NULL;
      ml->cell_face_n [i_part] = NULL;
    }
    ml->cell_vtx_idx[i_part] = NULL;
    ml->cell_vtx    [i_part] = NULL;
  }

  ml->is_elmt_select_by_user = malloc(sizeof(int *) * n_part);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    ml->is_elmt_select_by_user[i_part] = NULL;
  }

}


/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id            Pointer to \ref PDM_mesh_location object
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
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                  n_cell,
 const int                 *cell_face_idx,
 const int                 *cell_face,
 const PDM_g_num_t         *cell_ln_to_gn,
 const int                  n_face,
 const int                 *face_vtx_idx,
 const int                 *face_vtx,
 const PDM_g_num_t         *face_ln_to_gn,
 const int                  n_vtx,
 const double              *coords,
 const PDM_g_num_t         *vtx_ln_to_gn
)
{

  PDM_UNUSED(face_ln_to_gn);

  /*
   * Creation de mesh nodal
   */

  PDM_Mesh_nodal_coord_set (ml->mesh_nodal,
                            i_part,
                            n_vtx,
                            coords,
                            vtx_ln_to_gn,
                            PDM_OWNERSHIP_USER);



  ml->face_vtx_n[i_part]  = malloc (sizeof(PDM_l_num_t) * n_face);
  ml->cell_face_n[i_part] = malloc (sizeof(PDM_l_num_t) * n_cell);

  for (int i = 0; i < n_face; i++) {
    ml->face_vtx_n[i_part][i] = face_vtx_idx[i+1] - face_vtx_idx[i];
  }

  for (int i = 0; i < n_cell; i++) {
    ml->cell_face_n[i_part][i] = cell_face_idx[i+1] - cell_face_idx[i];
  }

  PDM_Mesh_nodal_cell3d_cellface_add(ml->mesh_nodal,
                                     i_part,
                                     n_cell,
                                     n_face,
                                     face_vtx_idx,
                                     ml->face_vtx_n[i_part],
                                     face_vtx,
                                     face_ln_to_gn,
                                     cell_face_idx,
                                     ml->cell_face_n[i_part],
                                     cell_face,
                                     cell_ln_to_gn,
                                     PDM_OWNERSHIP_KEEP);
}

/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id                     Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part                 Partition to define
 * \param [in]   n_cell                 Number of cells
 * \param [in]   is_elmt_select_by_user Flag to determine if user want or no to extract current cell
 *
 */
void
PDM_mesh_location_user_extract_set
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                 *is_elmt_select_by_user
)
{
  ml->is_elmt_select_by_user[i_part] = (int *) is_elmt_select_by_user;
}

/**
 *
 * \brief Set a part of a mesh (2d version)
 *
 * \param [in]   id            Pointer to \ref PDM_mesh_location object
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
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                  n_cell,
 const int                 *cell_edge_idx,
 const int                 *cell_edge,
 const PDM_g_num_t         *cell_ln_to_gn,
 const int                  n_edge,
 const int                 *edge_vtx_idx,
 const int                 *edge_vtx,
 const PDM_g_num_t         *edge_ln_to_gn,
 const int                  n_vtx,
 const double              *coords,
 const PDM_g_num_t         *vtx_ln_to_gn
)
{

  PDM_UNUSED (edge_ln_to_gn);

  /*
   * Creation de mesh nodal
   */

  PDM_Mesh_nodal_coord_set (ml->mesh_nodal,
                            i_part,
                            n_vtx,
                            coords,
                            vtx_ln_to_gn,
                            PDM_OWNERSHIP_USER);

  ml->face_vtx_n[i_part]  = malloc (sizeof(PDM_l_num_t) * n_edge);
  ml->cell_face_n[i_part] = malloc (sizeof(PDM_l_num_t) * n_cell);

  PDM_l_num_t *edge_vtx_nb  = ml->face_vtx_n[i_part];
  PDM_l_num_t *cell_edge_nb = ml->cell_face_n[i_part];

  for (int i = 0; i < n_edge; i++) {
    edge_vtx_nb[i] = edge_vtx_idx[i+1] - edge_vtx_idx[i];
  }

  for (int i = 0; i < n_cell; i++) {
    cell_edge_nb[i] = cell_edge_idx[i+1] - cell_edge_idx[i];
  }


  PDM_Mesh_nodal_cell2d_celledge_add (ml->mesh_nodal,
                                      i_part,
                                      n_cell,
                                      n_edge,
                                      edge_vtx_idx,
                                      edge_vtx_nb,
                                      edge_vtx,
                                      cell_edge_idx,
                                      cell_edge_nb,
                                      cell_edge,
                                      cell_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);
}





/**
 *
 * \brief Set the tolerance for bounding boxes
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   tol             Tolerance
 *
 */
void
PDM_mesh_location_tolerance_set
(
       PDM_mesh_location_t *ml,
 const double               tol
)
{

  ml->tolerance = tol;
}


/**
 *
 * \brief Set the method for computing location
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   method          Method
 *
 */
void
PDM_mesh_location_method_set
(
       PDM_mesh_location_t        *ml,
 const PDM_mesh_location_method_t  method
)
{

  ml->method = method;
}


/**
 *
 * \brief Get point location
 *
 * \param [in]   id                    Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud         Current cloud
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  n_points              Number of points in point cloud
 * \param [out]  coord                 Coordinates of points in point cloud
 * \param [out]  location              The global number of the closest element if the point is located,
 *                                     -1 otherwise
 *
 */

void
PDM_mesh_location_point_location_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_point_cloud,
 const int                   i_part,
       PDM_g_num_t         **location,
       double              **dist2,
       double              **projected_coord
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);

  *location        = pcloud->location[i_part];
  // TODO :Leak in python
  // *weights_idx     = pcloud->weights_idx[i_part];
  // *weights         = pcloud->weights[i_part];
  *projected_coord = pcloud->projected_coords[i_part];
  *dist2           = pcloud->dist2[i_part];

  ml->tag_point_location_get = 1;
}


/**
 *
 * \brief get cell vertex connectivity
 *
 * \param [in]   id                    Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  cell_vtx_idx          Index in (size = n_elt + 1)
 * \param [out]  cell_vtx              Cell vertex connectivity
 *
 */
void
PDM_mesh_location_cell_vertex_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_part,
       int                 **cell_vtx_idx,
       int                 **cell_vtx
)
{

  if (!ml->reverse_result) {
    PDM_error(__FILE__, __LINE__, 0, 
      "PDM_mesh_location_cell_vertex_get : Reverse results computation is disable, "
      "call PDM_mesh_location_reverse_results_enable to enable its\n");
  }

  *cell_vtx_idx = ml->cell_vtx_idx[i_part];
  *cell_vtx     = ml->cell_vtx[i_part];

  ml->tag_cell_vtx_get = 1;

}


/**
 *
 * \brief Get point list located in elements
 *
 * \param [in]   id                      Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part                  Index of partition of the mesh
 * \param [in]   i_point_cloud           Index of cloud
 * \param [out]  elt_pts_inside_idx      Points index (size = n_elt + 1)
 * \param [out]  points_gnum             Points global number
 * \param [out]  points_coords           Points coordinates
 * \param [out]  points_uvw              Points parametric coordinates in elements
 * \param [out]  points_weights_idx      Interpolation weights index (size = elt_pts_inside_idx[n_elt] + 1)
 * \param [out]  points_weights          Interpolation weights
 * \param [out]  points_dist2            Distance element-points (dist < 0 if the point is inside)
 * \param [out]  points_projected_coords Point projection on element if the point is outside
 *
 */

void
PDM_mesh_location_points_in_elt_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_part,
 const int                   i_point_cloud,
       int                 **elt_pts_inside_idx,
       PDM_g_num_t         **points_gnum,
       double              **points_coords,
       double              **points_uvw,
       int                 **points_weights_idx,
       double              **points_weights,
       double              **points_dist2,
       double              **points_projected_coords
)
{

  if (!ml->reverse_result) {
    PDM_error(__FILE__, __LINE__, 0, 
      "PDM_mesh_location_points_in_elt_get : Reverse results computation is disable, "
      "call PDM_mesh_location_reverse_results_enable to enable its\n");
  }

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  assert (ml->points_in_elements != NULL);
  assert (i_part < ml->points_in_elements->n_part);

  _points_in_element_t *_points_in_elements = ml->points_in_elements + i_point_cloud;

  *elt_pts_inside_idx      = _points_in_elements->pts_inside_idx[i_part];
  *points_gnum             = _points_in_elements->gnum[i_part];
  *points_coords           = _points_in_elements->coords[i_part];
  *points_uvw              = _points_in_elements->uvw[i_part];
  *points_weights_idx      = _points_in_elements->weights_idx[i_part];
  *points_weights          = _points_in_elements->weights[i_part];
  *points_dist2            = _points_in_elements->dist2[i_part];
  *points_projected_coords = _points_in_elements->projected_coords[i_part];

  ml->tag_points_in_elt_get = 1;
}


/**
 *
 * \brief Free a mesh location structure
 *
 * \param [in]  id       Pointer to \ref PDM_mesh_location object
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 */
void
PDM_mesh_location_free
(
 PDM_mesh_location_t  *ml
)
{

  /* Free point clouds */

  if (ml->point_clouds != NULL) {
    for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {
      _point_cloud_t *pcloud = ml->point_clouds + icloud;

      if (ml->reverse_result && ml->points_in_elements != NULL) {
        _points_in_element_t *_points_in_elements = ml->points_in_elements + icloud;
        if(( ml->owner == PDM_OWNERSHIP_KEEP ) ||
           ( ml->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !ml->tag_points_in_elt_get)) {
          for (int i_part = 0; i_part < _points_in_elements->n_part; ++i_part) {
            free (_points_in_elements->pts_inside_idx[i_part]);
            free (_points_in_elements->gnum[i_part]);
            free (_points_in_elements->uvw[i_part]);
            free (_points_in_elements->coords[i_part]);
            free (_points_in_elements->projected_coords[i_part]);
            free (_points_in_elements->weights_idx[i_part]);
            free (_points_in_elements->weights[i_part]);
            free (_points_in_elements->dist2[i_part]);
          }
        }
        free (_points_in_elements->pts_inside_idx);
        free (_points_in_elements->n_elts);
        free (_points_in_elements->gnum);
        free (_points_in_elements->uvw);
        free (_points_in_elements->coords);
        free (_points_in_elements->projected_coords);
        free (_points_in_elements->weights_idx);
        free (_points_in_elements->weights);
        free (_points_in_elements->dist2);
        // free (_points_in_elements);
      }

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
        if(( ml->owner == PDM_OWNERSHIP_KEEP ) ||
           ( ml->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !ml->tag_point_location_get)) {
          for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
            if (pcloud->location[ipart] != NULL) {
              free (pcloud->location[ipart]);
            }
            if (pcloud->dist2[ipart] != NULL) {
              free (pcloud->dist2[ipart]);
            }
            if (pcloud->projected_coords[ipart] != NULL) {
              free (pcloud->projected_coords[ipart]);
            }
          }
        }
        free (pcloud->location);
        free (pcloud->dist2);
        free (pcloud->projected_coords);
      }

      if (pcloud->uvw != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->uvw[ipart] != NULL) {
            free (pcloud->uvw[ipart]);
          }
        }
        free (pcloud->uvw);
      }

      if (pcloud->weights_idx != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->weights_idx[ipart] != NULL) {
            free (pcloud->weights_idx[ipart]);
          }
        }
        free (pcloud->weights_idx);
      }

      if (pcloud->weights != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->weights[ipart] != NULL) {
            free (pcloud->weights[ipart]);
          }
        }
        free (pcloud->weights);
      }


      if (pcloud->n_located != NULL) {
        free (pcloud->n_located);
      }

      if (pcloud->n_un_located != NULL) {
        free (pcloud->n_un_located);
      }

      if (pcloud->located != NULL) {

        if(( ml->owner == PDM_OWNERSHIP_KEEP ) ||
           ( ml->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !ml->tag_located_get)) {
          for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
            if (pcloud->located[ipart] != NULL) {
              free (pcloud->located[ipart]);
            }
          }
        }
        free (pcloud->located);
        pcloud->located = NULL;
      }

      if (pcloud->un_located != NULL) {
        if(( ml->owner == PDM_OWNERSHIP_KEEP ) ||
           ( ml->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !ml->tag_unlocated_get)) {
          for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
            if (pcloud->un_located[ipart] != NULL) {
              free (pcloud->un_located[ipart]);
            }
          }
        }
        free (pcloud->un_located);
        pcloud->un_located = NULL;
      }

    }
    if (ml->points_in_elements != NULL) {
      free (ml->points_in_elements);
    }
    ml->points_in_elements = NULL;

    free (ml->point_clouds);
    ml->point_clouds = NULL;
  }

  /* Free mesh nodal */
  //PDM_Mesh_nodal_partial_free (ml->mesh_nodal);?

  if (ml->mesh_nodal != NULL) {
    int _n_part = PDM_Mesh_nodal_n_part_get(ml->mesh_nodal);

    if (ml->cell_vtx_idx != NULL) {
      for (int i = 0; i< _n_part; i++) {
        if(( ml->owner == PDM_OWNERSHIP_KEEP ) ||
           ( ml->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !ml->tag_cell_vtx_get)) {
          if(ml->cell_vtx_idx[i] != NULL) {
            free(ml->cell_vtx[i]);
            free(ml->cell_vtx_idx[i]);
          }
        }
      }
      free(ml->cell_vtx);
      free(ml->cell_vtx_idx);
    }

    ml->cell_vtx_idx = NULL;
    ml->cell_vtx = NULL;

    if (!ml->shared_nodal) {

      PDM_Mesh_nodal_free (ml->mesh_nodal);

      if(ml->cell_face_n != NULL){
        for (int i = 0; i< _n_part; i++) {
          if(ml->cell_face_n[i] != NULL) {
            free(ml->cell_face_n[i]);
          }
        }
        free (ml->cell_face_n);
      }

      if(ml->face_vtx_n != NULL){
        for (int i = 0; i< _n_part; i++) {
          if(ml->face_vtx_n[i] != NULL) {
            free(ml->face_vtx_n[i]);
          }
        }
        free (ml->face_vtx_n);
      }
    }
  }

  if(ml->is_elmt_select_by_user != NULL) {
    free(ml->is_elmt_select_by_user);
  }


  for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {
    if (ml->ptp_ownership[icloud] == PDM_OWNERSHIP_KEEP) {
      PDM_part_to_part_free(ml->ptp[icloud]);
      ml->ptp[icloud] = NULL;
    }
  }
  free(ml->ptp);
  free(ml->ptp_ownership);

  PDM_timer_free(ml->timer);

  free(ml);
}

/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Pointer to \ref PDM_mesh_location object
 *
 */
void
PDM_mesh_location_dump_times
(
PDM_mesh_location_t *ml
)
{

  double t1 = ml->times_elapsed[END] - ml->times_elapsed[BEGIN];
  double t2 = ml->times_cpu[END] - ml->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

  double t_elaps_max[NTIMER_MESH_LOCATION];
  PDM_MPI_Allreduce (ml->times_elapsed,
                     t_elaps_max,
                     NTIMER_MESH_LOCATION,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     ml->comm);

  double t_cpu_max[NTIMER_MESH_LOCATION];
  PDM_MPI_Allreduce (ml->times_cpu,
                     t_cpu_max, NTIMER_MESH_LOCATION,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     ml->comm);

  int rank;
  PDM_MPI_Comm_rank (ml->comm, &rank);

  if (rank == 0) {

    PDM_printf( "mesh_location timer : all (elapsed and cpu) :                                   "
                " %12.5es %12.5es\n",
                t1max, t2max);

    PDM_printf( "mesh_location timer : build bounding boxes + extract mesh (elapsed and cpu) :   "
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_BOUNDING_BOXES],
                t_cpu_max[BUILD_BOUNDING_BOXES]);

    PDM_printf( "mesh_location timer : Store connectivity (elapsed and cpu) :                    "
                " %12.5es %12.5es\n",
                t_elaps_max[STORE_CONNECTIVITY],
                t_cpu_max[STORE_CONNECTIVITY]);

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

    PDM_printf( "mesh_location timer : merge location data - step 1 (elapsed and cpu) :          "
                " %12.5es %12.5es\n",
                t_elaps_max[MERGE_LOCATION_DATA_STEP1],
                t_cpu_max[MERGE_LOCATION_DATA_STEP1]);

    PDM_printf( "mesh_location timer : merge location data - step 2 (elapsed and cpu) :          "
                " %12.5es %12.5es\n",
                t_elaps_max[MERGE_LOCATION_DATA_STEP2],
                t_cpu_max[MERGE_LOCATION_DATA_STEP2]);

    PDM_printf( "mesh_location timer : merge location data - step 3 (elapsed and cpu) :          "
                " %12.5es %12.5es\n",
                t_elaps_max[MERGE_LOCATION_DATA_STEP3],
                t_cpu_max[MERGE_LOCATION_DATA_STEP3]);

    PDM_printf( "mesh_location timer : compress location data (elapsed and cpu) :                "
                " %12.5es %12.5es\n",
                t_elaps_max[COMPRESS_LOCATION_DATA],
                t_cpu_max[COMPRESS_LOCATION_DATA]);

    // PDM_printf( "mesh_location timer : reverse location data (elapsed and cpu) :                 "
    //             " %12.5es %12.5es\n",
    //             t_elaps_max[REVERSE_LOCATION_DATA],
    //             t_cpu_max[REVERSE_LOCATION_DATA]);

    PDM_printf( "mesh_location timer : reverse location data - ptb (elapsed and cpu) :           "
                " %12.5es %12.5es\n",
                t_elaps_max[REVERSE_LOCATION_DATA_PTB],
                t_cpu_max[REVERSE_LOCATION_DATA_PTB]);

    PDM_printf( "mesh_location timer : reverse location data - btp step 1 (elapsed and cpu) :    "
                " %12.5es %12.5es\n",
                t_elaps_max[REVERSE_LOCATION_DATA_BTP1],
                t_cpu_max[REVERSE_LOCATION_DATA_BTP1]);

    PDM_printf( "mesh_location timer : reverse location data - btp step 2 (elapsed and cpu) :    "
                " %12.5es %12.5es\n",
                t_elaps_max[REVERSE_LOCATION_DATA_BTP2],
                t_cpu_max[REVERSE_LOCATION_DATA_BTP2]);

    PDM_printf( "mesh_location timer : reverse location data - btp step 3 (elapsed and cpu) :    "
                " %12.5es %12.5es\n",
                t_elaps_max[REVERSE_LOCATION_DATA_BTP3],
                t_cpu_max[REVERSE_LOCATION_DATA_BTP3]);

    PDM_printf( "mesh_location timer : reverse location data - uvw (elapsed and cpu) :           "
                " %12.5es %12.5es\n",
                t_elaps_max[REVERSE_LOCATION_DATA_UVW],
                t_cpu_max[REVERSE_LOCATION_DATA_UVW]);
  }
}



static void _export_boxes
(
 const char        *filename,
 const int          n_box,
 const double      *box_extents,
 const PDM_g_num_t *box_g_num
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nboxes\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_box);
  for (int i = 0; i < n_box; i++) {
    const double *e = box_extents + 6*i;
    fprintf(f, "%f %f %f\n", e[0], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[5]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[5]);
  }

  fprintf(f, "CELLS %d %d\n", n_box, 9*n_box);
  for (int i = 0; i < n_box; i++) {
    int j = 8*i;
    fprintf(f, "8 %d %d %d %d %d %d %d %d\n", j, j+1, j+2, j+3, j+4, j+5, j+6, j+7);
  }

  fprintf(f, "CELL_TYPES %d\n", n_box);
  for (int i = 0; i < n_box; i++) {
    fprintf(f, "12\n");
  }

  fprintf(f, "CELL_DATA %d\n", n_box);
  fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
  for (int i = 0; i < n_box; i++) {
    fprintf(f, ""PDM_FMT_G_NUM"\n", box_g_num[i]);
  }

  fclose(f);
}


static void _export_point_cloud
(
 char         *filename,
 int           n_part,
 int          *n_pts,
 double      **coord,
 PDM_g_num_t **g_num,
 PDM_g_num_t **parent_g_num
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\noctree points\nASCII\nDATASET UNSTRUCTURED_GRID\n");

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

  if (parent_g_num != NULL) {
    fprintf(f, "FIELD FieldData 1\n");
    fprintf(f, "parent_gnum 1 %d int\n", n_pts_t);
    for (int ipart = 0; ipart < n_part; ipart++) {
      for (int i = 0; i < n_pts[ipart]; i++) {
        fprintf(f, ""PDM_FMT_G_NUM"\n", parent_g_num[ipart][i]);
      }
    }
  }

  fclose(f);
}


static void
_point_cloud_extract_selection
(
 PDM_MPI_Comm            comm,
 const _point_cloud_t   *parent_cloud,
 int                    *n_select_pts,
 int                   **select_pts_l_num,
 PDM_g_num_t          ***select_pts_parent_g_num,
 PDM_g_num_t          ***select_pts_g_num,
 double               ***select_pts_coord
 )
{
  const int n_part = parent_cloud->n_part;

  *select_pts_parent_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  *select_pts_g_num        = malloc (sizeof(PDM_g_num_t *) * n_part);
  *select_pts_coord        = malloc (sizeof(double *)      * n_part);

  PDM_gen_gnum_t *gen_gnum = PDM_gnum_create (3,
                                              n_part,
                                              PDM_FALSE,
                                              0.,
                                              comm,
                                              PDM_OWNERSHIP_USER);


  for (int ipart = 0; ipart < n_part; ipart++) {

    (*select_pts_parent_g_num)[ipart] = malloc (sizeof(PDM_g_num_t *) * n_select_pts[ipart]);
    (*select_pts_g_num)[ipart] = NULL;
    (*select_pts_coord)[ipart] = malloc (sizeof(double) * n_select_pts[ipart] * 3);

    for (int i = 0; i < n_select_pts[ipart]; i++) {
      int j = select_pts_l_num[ipart][i];

      (*select_pts_parent_g_num)[ipart][i] = parent_cloud->gnum[ipart][j];
      for (int k = 0; k < 3; k++) {
        (*select_pts_coord)[ipart][3*i + k] = parent_cloud->coords[ipart][3*j + k];
      }
    }

    PDM_gnum_set_from_parents (gen_gnum,
                               ipart,
                               n_select_pts[ipart],
                               (*select_pts_parent_g_num)[ipart]);
  }

  /* Generate a new global numbering for selected points */
  PDM_gnum_compute (gen_gnum);

  for (int ipart = 0; ipart < n_part; ipart++) {
    (*select_pts_g_num)[ipart] = PDM_gnum_get (gen_gnum, ipart);
  }

  /* Free memory */
  PDM_gnum_free (gen_gnum);
}



/**
 *
 * \brief Compute point location
 *
 * \param [in]   id  Pointer to \ref PDM_mesh_location object
 *
 */

void
PDM_mesh_location_compute
(
PDM_mesh_location_t        *ml
)
{
  const double eps_dist = 1.e-10;
  const double tolerance = 1e-6;

  const int dbg_enabled = 0;
  const int log_timer   = 0;
  const int dim = 3;

  const int octree_depth_max = 31;
  const int octree_points_in_leaf_max = 1;
  const int octree_build_leaf_neighbours = 0;
  PDM_para_octree_t *octree = NULL;

  const int VISU = 0;
  int allow_extraction = 1;
  float extraction_threshold = 0.5; // max size ratio between extracted and original meshes

  //
  // TODO: - ok - Conserver que "allow_extraction = 1 et donc suprrimer la variable et les else"
  //       - ok - Mettre l'extraction mesh dans une fonction dédiée
  //       - Faire une extraction par paire nuage/maillage ?
  //       - Voir l'utilite d'un appel extract_part
  //       - part_to_block_geom sur select_pts (pcloud_*) (normalement ca casse pas tout, ajout d'un part_to_part à ma fin ?)
  //       - part_to_block_geom sur select_boxes (prendre le centre ? )
  //       - "Store cell vertex connectivity" nécessaire si calcul des poids : voir si on peut pas optimiser
  //         (concatenation des connectivites cell->vtx des blocs en un seul tableau)
  //       - Faire un dbbtree par paire nuage/maillage ?
  //       - PDM_para_octree_points_inside_boxes : 
  //            * Faire un resultat suivant les boites intersectant localement l'octree
  //            * Ajout d'un part_to_block geometrique pour redistribuer les boites de maniere equilibrée avec le nombre de points inclus en poids (box_set)
  //       - PDM_dbbtree_points_inside_boxes :
  //            * Mettre le résultat sous forme de points
  //            * Repasser sur une vision blocs par boites pour la suite de l'algo
  //       - Dans le case ajouter l'octree clone de dbbtree
  //       - Dans le case ajouter la collision d'arbres
  //       - 2eme extraction : ok (simplification supprimer les else car fait tout le temps)
  //       - Redistribution (part_to_block + block_to_part) des parent gnum des pts dans la distribution des éléments (A simplifier, a revoir ?)
  //       - PDM_point_location_nodal - pas de modification pour l'instant
  //       - Merge location data : il y a du boulot !
  //          - 1 : premier part_to_block (echanger seulement la distance pour le prendre min)                 
  //          - 2 : Ajouter un retour qui indique la vraie cellule dan l'etat précédent
  //          - Depuis l'etat des cellules pour le pdm_point_location_nodale
  //          -   part_to_part vers le ragement initial des points
  //          -   part_to_part vers le ragement initial des cellules
  //       - Plusieurs nuages de points : recouvrir les communications finales par le prochain nuage

  char *env_var = getenv ("ALLOW_EXTRACTION");
  if (env_var != NULL) {
    allow_extraction = atoi(env_var);
  }
  // int use_extracted_pts  = allow_extraction;
  int *use_extracted_pts  = malloc(ml->n_point_cloud * sizeof(int));
  for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {
    use_extracted_pts[icloud] = allow_extraction;
  }
  int use_extracted_mesh = allow_extraction;

  int          **n_select_pts = NULL;
  PDM_g_num_t ***select_pts_parent_g_num = NULL;
  PDM_g_num_t ***select_pts_g_num = NULL;
  double      ***select_pts_coord = NULL;

  int my_rank;
  PDM_MPI_Comm_rank (ml->comm, &my_rank);

  int n_procs;
  PDM_MPI_Comm_size (ml->comm, &n_procs);

  // int USE_OCTREE_BTSHARED = 0;
  // env_var = getenv ("USE_OCTREE_BTSHARED");
  // if (env_var != NULL) {
  //   USE_OCTREE_BTSHARED = atoi(env_var);
  // }
  // if (dbg_enabled && my_rank == 0 && ml->method == PDM_MESH_LOCATION_OCTREE)
  //   printf("USE_OCTREE_BTSHARED = %d\n", USE_OCTREE_BTSHARED);

  // int USE_OCTREE_COPIES = 1;
  // env_var = getenv ("USE_OCTREE_COPIES");
  // if (env_var != NULL) {
  //   USE_OCTREE_COPIES = atoi(env_var);
  // }
  // if (dbg_enabled && my_rank == 0) // && ml->method == PDM_MESH_LOCATION_OCTREE)
  //   printf("USE_OCTREE_COPIES = %d\n", USE_OCTREE_COPIES);

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  PDM_MPI_Barrier (ml->comm);
  ml->times_elapsed[BEGIN] = PDM_timer_elapsed(ml->timer);
  ml->times_cpu[BEGIN]     = PDM_timer_cpu(ml->timer);
  ml->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(ml->timer);
  ml->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(ml->timer);

  b_t_elapsed = ml->times_elapsed[BEGIN];
  b_t_cpu     = ml->times_cpu[BEGIN];
  b_t_cpu_u   = ml->times_cpu_u[BEGIN];
  b_t_cpu_s   = ml->times_cpu_s[BEGIN];
  PDM_timer_resume(ml->timer);

  /*
   * Build the bounding boxes of mesh elements
   */

  int n_blocks = 0;
  int n_parts  = 0;
  int *blocks_id = NULL;

  if (ml->mesh_nodal != NULL) {
    n_blocks = PDM_Mesh_nodal_n_blocks_get (ml->mesh_nodal);
    n_parts  = PDM_Mesh_nodal_n_part_get (ml->mesh_nodal);
    blocks_id = PDM_Mesh_nodal_blocks_id_get (ml->mesh_nodal);
  }

  int n_boxes = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {
    n_boxes += PDM_Mesh_nodal_n_cell_get (ml->mesh_nodal,
                                          ipart);
  }

  PDM_g_num_t *box_g_num   = malloc (sizeof(PDM_g_num_t) * n_boxes);
  double      *box_extents = malloc (sizeof(double)      * n_boxes * 6);

  int n_select_boxes = n_boxes;
  PDM_g_num_t *select_box_parent_g_num = NULL;
  PDM_g_num_t *select_box_g_num   = box_g_num;
  double      *select_box_extents = box_extents;

  int  **n_select_elt = NULL;
  int ***select_elt_l_num = NULL;

  int ibox = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {

    int id_block = blocks_id[iblock];

    for (int ipart = 0; ipart < n_parts; ipart++) {
      /* get element extents */
      PDM_Mesh_nodal_compute_cell_extents (ml->mesh_nodal,
                                           id_block,
                                           ipart,
                                           ml->tolerance,
                                           box_extents + 6*ibox);

      /* get elements gnum */
      PDM_g_num_t *_gnum = PDM_Mesh_nodal_g_num_get (ml->mesh_nodal,
                                                     id_block,
                                                     ipart);

      int n_elt = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                  id_block,
                                                  ipart);

      for (int ielt = 0; ielt < n_elt; ielt++) {
        box_g_num[ibox] = _gnum[ielt];
        ibox++;
      }
    }
  }

  if (VISU) {
    char filename[999];
    sprintf(filename, "mesh_boxes_%3.3d.vtk", my_rank);
    _export_boxes (filename, n_boxes, box_extents, box_g_num);


    PDM_writer_t *cs = PDM_writer_create("Ensight",
                                         PDM_WRITER_FMT_ASCII,
                                         PDM_WRITER_TOPO_CST,
                                         PDM_WRITER_OFF,
                                         "mesh_location",
                                         "source_mesh",
                                         PDM_MPI_COMM_WORLD,
                                         PDM_IO_KIND_MPI_SIMPLE,
                                         1.,
                                         NULL);

    int id_geom = PDM_writer_geom_create_from_mesh_nodal (cs,
                                                          "source_mesh_geom",
                                                          ml->mesh_nodal);

    PDM_writer_step_beg(cs, 0.);

    PDM_writer_geom_write(cs,
                          id_geom);

    PDM_writer_step_end(cs);
    PDM_writer_free(cs);
  }


  if (allow_extraction) {

    /*
     *  Compute global extents of source mesh
     */

    double mesh_extents[6] = { HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                              -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    for (int i = 0; i < n_boxes; i++) {
      for (int j = 0; j < 3; j++) {
        mesh_extents[j]   = PDM_MIN (mesh_extents[j],   box_extents[6*i + j]);
        mesh_extents[j+3] = PDM_MAX (mesh_extents[j+3], box_extents[6*i + j + 3]);
      }
    }

    double g_mesh_extents[6];
    PDM_MPI_Allreduce (mesh_extents,   g_mesh_extents,   3,
                       PDM_MPI_DOUBLE, PDM_MPI_MIN, ml->comm);
    PDM_MPI_Allreduce (mesh_extents+3, g_mesh_extents+3, 3,
                       PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

    if (dbg_enabled && my_rank == 0) {
      printf("g_mesh_extents = %f %f %f / %f %f %f\n",
             g_mesh_extents[0], g_mesh_extents[1], g_mesh_extents[2],
             g_mesh_extents[3], g_mesh_extents[4], g_mesh_extents[5]);
    }

    /*
     *  Extract points that intersect the source mesh global extents
     *  Brute force (could be accelerated using octree)
     */

    n_select_pts = malloc (sizeof(int *)  * ml->n_point_cloud);
    select_pts_parent_g_num = malloc (sizeof(PDM_g_num_t **) * ml->n_point_cloud);
    select_pts_g_num = malloc (sizeof(PDM_g_num_t **) * ml->n_point_cloud);
    select_pts_coord = malloc (sizeof(double **) * ml->n_point_cloud);

    for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {

      _point_cloud_t *pcloud = ml->point_clouds + icloud;

      n_select_pts[icloud] = malloc (sizeof(int) * pcloud->n_part);
      int **select_pts_l_num = malloc (sizeof(int *) * pcloud->n_part);

      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        n_select_pts[icloud][ipart] = 0;
        select_pts_l_num[ipart] = malloc (sizeof(int) * pcloud->n_points[ipart]);

        for (int i = 0; i < pcloud->n_points[ipart]; i++) {
          int inside = 1;
          for (int j = 0; j < 3; j++) {
            if (pcloud->coords[ipart][3*i+j] < g_mesh_extents[j] ||
                pcloud->coords[ipart][3*i+j] > g_mesh_extents[j+3]) {
              inside = 0;
              break;
            }
          }

          if (inside) {
            select_pts_l_num[ipart][n_select_pts[icloud][ipart]++] = i;// +1?
          }
        }

        select_pts_l_num[ipart] = realloc (select_pts_l_num[ipart],
                                           sizeof(int) * n_select_pts[icloud][ipart]);
      } // End of loop on parts


      PDM_g_num_t l_n_pts[2] = {0, 0};
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        l_n_pts[0] += pcloud->n_points[ipart];
        l_n_pts[1] += n_select_pts[icloud][ipart];
      }

      PDM_g_num_t g_n_pts[2];
      PDM_MPI_Allreduce (l_n_pts, g_n_pts, 2, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ml->comm);

      use_extracted_pts[icloud] = (g_n_pts[1] < extraction_threshold * g_n_pts[0]);

      if (use_extracted_pts[icloud]) {
        _point_cloud_extract_selection (ml->comm,
                                        pcloud,
                                        n_select_pts[icloud],
                                        select_pts_l_num,
                                        &(select_pts_parent_g_num[icloud]),
                                        &(select_pts_g_num[icloud]),
                                        &(select_pts_coord[icloud]));
      } else {
        free (n_select_pts[icloud]);
        n_select_pts[icloud] = pcloud->n_points;
        select_pts_parent_g_num[icloud] = pcloud->gnum;
        select_pts_g_num[icloud] = pcloud->gnum;
        select_pts_coord[icloud] = pcloud->coords;
      }

      if (allow_extraction) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          free (select_pts_l_num[ipart]);
        }
        free (select_pts_l_num);
      }

      if (VISU) {
        char filename[999];

        sprintf(filename, "parent_cloud_%d_%3.3d.vtk", icloud, my_rank);
        _export_point_cloud (filename,
                             pcloud->n_part,
                             pcloud->n_points,
                             pcloud->coords,
                             pcloud->gnum,
                             NULL);

        sprintf(filename, "extracted_cloud_%d_%3.3d.vtk", icloud, my_rank);
        _export_point_cloud (filename,
                             pcloud->n_part,
                             n_select_pts[icloud],
                             select_pts_coord[icloud],
                             select_pts_g_num[icloud],
                             select_pts_parent_g_num[icloud]);
      }



      if (dbg_enabled && my_rank == 0 && use_extracted_pts[icloud]) {
        printf("extract "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM" pts from cloud #%d ("PDM_FMT_G_NUM"%%)\n",
               g_n_pts[1], g_n_pts[0], icloud, (100 * g_n_pts[1]) / g_n_pts[0]);
      }

    } // End loop on point clouds




    /*
     *  Compute global extents of extracted point clouds
     */
    double pts_extents[6] = { HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                             -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {

      _point_cloud_t *pcloud = ml->point_clouds + icloud;

      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        for (int i = 0; i < n_select_pts[icloud][ipart]; i++) {
          for (int j = 0; j < 3; j++) {
            pts_extents[j]   = PDM_MIN (pts_extents[j],
                                        select_pts_coord[icloud][ipart][3*i + j]);
            pts_extents[j+3] = PDM_MAX (pts_extents[j+3],
                                        select_pts_coord[icloud][ipart][3*i + j]);
          }
        }
      }
    }

    double g_pts_extents[6];
    PDM_MPI_Allreduce (pts_extents,   g_pts_extents,   3,
                       PDM_MPI_DOUBLE, PDM_MPI_MIN, ml->comm);
    PDM_MPI_Allreduce (pts_extents+3, g_pts_extents+3, 3,
                       PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

    if (dbg_enabled && my_rank == 0) {
      printf("g_pts_extents = %f %f %f / %f %f %f\n",
             g_pts_extents[0], g_pts_extents[1], g_pts_extents[2],
             g_pts_extents[3], g_pts_extents[4], g_pts_extents[5]);
    }


    /*
     *  Select elements whose bounding box intersect
     *  the global extents of the extracted point clouds
     *  Brute force (could be accelerated using bbtree)
     */

    n_select_elt     = malloc (n_blocks * sizeof(int  *));
    select_elt_l_num = malloc (n_blocks * sizeof(int **));

    ibox = 0;
    for (int iblock = 0; iblock < n_blocks; iblock++) {

      n_select_elt    [iblock] = malloc (n_parts * sizeof(int  ));
      select_elt_l_num[iblock] = malloc (n_parts * sizeof(int *));

      int id_block = blocks_id[iblock];

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int part_n_elt = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                         id_block,
                                                         ipart);

        n_select_elt[iblock][ipart] = 0;
        select_elt_l_num[iblock][ipart] = malloc (sizeof(int) * part_n_elt);

        for (int ielt = 0; ielt < part_n_elt; ielt++) {

          double *box_min = box_extents + 6*ibox;
          double *box_max = box_min + 3;

          int intersect = 1;
          for (int j = 0; j < 3; j++) {
            if (box_min[j] > g_pts_extents[j+3] ||
                box_max[j] < g_pts_extents[j]) {
              intersect = 0;
              break;
            }
          }

          if (intersect) {
            select_elt_l_num[iblock][ipart][n_select_elt[iblock][ipart]++] = ielt;// +1?
          }

          ibox++;
        }

        select_elt_l_num[iblock][ipart] = realloc (select_elt_l_num[iblock][ipart],
                                                   sizeof(int) * n_select_elt[iblock][ipart]);
      } // End of loop on parts
    } // End of loop on nodal blocks

    PDM_g_num_t l_n_elt[2];

    l_n_elt[0] = (PDM_g_num_t) n_boxes;
    l_n_elt[1] = 0;
    for (int iblock = 0; iblock < n_blocks; iblock++) {
      for (int ipart = 0; ipart < n_parts; ipart++) {
        l_n_elt[1] += n_select_elt[iblock][ipart];
      }
    }

    PDM_g_num_t g_n_elt[2];
    PDM_MPI_Allreduce (l_n_elt, g_n_elt, 2, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ml->comm);

    use_extracted_mesh = (g_n_elt[1] < extraction_threshold * g_n_elt[0]);


    if (dbg_enabled && my_rank == 0 && use_extracted_mesh) {
      printf("extract "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM" mesh elts ("PDM_FMT_G_NUM"%%)\n",
             g_n_elt[1], g_n_elt[0], (100 * g_n_elt[1]) / g_n_elt[0]);
    }


    if (use_extracted_mesh) {
      n_select_boxes = 0;
      for (int iblock = 0; iblock < n_blocks; iblock++) {
        for (int ipart = 0; ipart < n_parts; ipart++) {
          n_select_boxes += n_select_elt[iblock][ipart];
        }
      }

      select_box_parent_g_num = malloc (sizeof(PDM_g_num_t) * n_select_boxes);
      select_box_extents = malloc (sizeof(double)      * n_select_boxes * 6);
      ibox = 0;
      int idx_elt = 0;
      for (int iblock = 0; iblock < n_blocks; iblock++) {
        int id_block = blocks_id[iblock];

        for (int ipart = 0; ipart < n_parts; ipart++) {
          int part_n_elt = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                           id_block,
                                                           ipart);

          for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
            int ielt = idx_elt + select_elt_l_num[iblock][ipart][i];
            select_box_parent_g_num[ibox] = box_g_num[ielt];
            for (int j = 0; j < 6; j++) {
              select_box_extents[6*ibox + j] = box_extents[6*ielt + j];
            }
            ibox++;
          }

          idx_elt += part_n_elt;
        } // End of loop on parts
      } // End of loop on blocks

      /*
       *  Generate a new global numbering for selected boxes
       */
      PDM_gen_gnum_t *gen_gnum_boxes = PDM_gnum_create (3,
                                                        1,
                                                        PDM_FALSE,
                                                        0.,
                                                        ml->comm,
                                                        PDM_OWNERSHIP_USER);
      PDM_gnum_set_from_parents (gen_gnum_boxes,
                                 0,
                                 n_select_boxes,
                                 select_box_parent_g_num);

      PDM_gnum_compute (gen_gnum_boxes);

      select_box_g_num = PDM_gnum_get (gen_gnum_boxes, 0);

      PDM_gnum_free (gen_gnum_boxes);


      if (VISU) {
        char filename[999];
        sprintf(filename, "extracted_mesh_boxes_%3.3d.vtk", my_rank);
        _export_boxes (filename, n_select_boxes, select_box_extents, select_box_g_num);
      }

    }
  }


  if (0) {
    PDM_g_num_t ln_elt = n_select_boxes;
    PDM_g_num_t gn_elt;
    PDM_MPI_Allreduce (&ln_elt, &gn_elt, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ml->comm);
    if (my_rank == 0) {
      printf("gn_elt = "PDM_FMT_G_NUM"\n", gn_elt);
    }
  }


  PDM_MPI_Barrier (ml->comm);
  PDM_timer_hang_on(ml->timer);
  e_t_elapsed = PDM_timer_elapsed(ml->timer);
  e_t_cpu     = PDM_timer_cpu(ml->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

  ml->times_elapsed[STORE_CONNECTIVITY] += e_t_elapsed - b_t_elapsed;
  ml->times_cpu[STORE_CONNECTIVITY]     += e_t_cpu     - b_t_cpu;
  ml->times_cpu_u[STORE_CONNECTIVITY]   += e_t_cpu_u   - b_t_cpu_u;
  ml->times_cpu_s[STORE_CONNECTIVITY]   += e_t_cpu_s   - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;
  PDM_timer_resume(ml->timer);

  /*
   *  Store cell vertex connectivity
   *  ------------------------------
   */

  int **n_vtx_per_elt = malloc (sizeof(int *) * n_parts);

  for (int ipart = 0; ipart < n_parts; ipart++) {
    int n_elt = PDM_Mesh_nodal_n_cell_get (ml->mesh_nodal,
                                           ipart);

    n_vtx_per_elt[ipart] = malloc (sizeof(int) * n_elt);
    // ml->cell_vtx_idx[ipart] = malloc (sizeof(int) * (n_elt+1));
    // ml->cell_vtx_idx[ipart][0] = 0;
  }

  PDM_g_num_t n_g_cell = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {
    int n_elt = PDM_Mesh_nodal_n_cell_get (ml->mesh_nodal,
                                           ipart);
    int ielt = 0;
    for (int iblock = 0; iblock < n_blocks; iblock++) {
      int id_block = blocks_id[iblock];

      PDM_Mesh_nodal_elt_t t_elt = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                  id_block);

      int *parent_num = PDM_Mesh_nodal_block_parent_num_get (ml->mesh_nodal,
                                                             id_block,
                                                             ipart);

      PDM_g_num_t *g_num = PDM_Mesh_nodal_g_num_get (ml->mesh_nodal,
                                                      id_block,
                                                      ipart);

      int n_elt_block = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                        id_block,
                                                        ipart);

      for (int i = 0; i < n_elt_block; i++) {
        n_g_cell = PDM_MAX (n_g_cell, g_num[i]);
      }
      int n_vtx_elt = 0;
      switch (t_elt) {
      case PDM_MESH_NODAL_POINT:
        n_vtx_elt = 1;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_BAR2:
        n_vtx_elt = 2;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_TRIA3:
        n_vtx_elt = 3;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_QUAD4:
      case PDM_MESH_NODAL_TETRA4:
        n_vtx_elt = 4;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_PYRAMID5:
        n_vtx_elt = 5;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_PRISM6:
        n_vtx_elt = 6;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_HEXA8:
        n_vtx_elt = 8;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_POLY_2D:
        {
          int *connec_idx;
          int *connec;
          PDM_Mesh_nodal_block_poly2d_get (ml->mesh_nodal,
                                           id_block,
                                           ipart,
                                          &connec_idx,
                                          &connec);
          if (parent_num != NULL) {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][parent_num[i]] = connec_idx[i+1] - connec_idx[i];
            }
          }
          else {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][ielt++] = connec_idx[i+1] - connec_idx[i];
            }
          }
          break;
        }
      case PDM_MESH_NODAL_POLY_3D:
        {
          int *connec_idx;
          int *connec;
          PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (ml->mesh_nodal,
                                                            id_block,
                                                            ipart,
                                                            &connec_idx,
                                                            &connec);
          if (parent_num != NULL) {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][parent_num[i]] = connec_idx[i+1] - connec_idx[i];
            }
          }
          else {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][ielt++] = connec_idx[i+1] - connec_idx[i];
            }
          }
          break;
        }
      default :
        PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_location error : Bad Element type\n");
      }
    }
    
    ml->cell_vtx_idx[ipart] = malloc (sizeof(int) * (n_elt+1));
    ml->cell_vtx_idx[ipart][0] = 0;
    for (int i = 0; i < n_elt; i++) {
      ml->cell_vtx_idx[ipart][i+1] = ml->cell_vtx_idx[ipart][i] + n_vtx_per_elt[ipart][i];
    }
    ml->cell_vtx[ipart] = malloc (sizeof(int) * ml->cell_vtx_idx[ipart][n_elt]);
  }

  PDM_g_num_t _n_g_cell = 0;
  PDM_MPI_Allreduce (&n_g_cell, &_n_g_cell, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, ml->comm);
  n_g_cell = _n_g_cell;

  for (int ipart = 0; ipart < n_parts; ipart++) {
    int ielt = 0;
    for (int iblock = 0; iblock < n_blocks; iblock++) {
      int id_block = blocks_id[iblock];

      PDM_Mesh_nodal_elt_t t_elt = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                  id_block);

      int *parent_num = PDM_Mesh_nodal_block_parent_num_get (ml->mesh_nodal,
                                                             id_block,
                                                             ipart);

      int n_elt_block = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                        id_block,
                                                        ipart);

      int n_vtx_elt = 0;
      switch (t_elt) {
      case PDM_MESH_NODAL_POINT:
      case PDM_MESH_NODAL_BAR2:
      case PDM_MESH_NODAL_TRIA3:
      case PDM_MESH_NODAL_QUAD4:
      case PDM_MESH_NODAL_TETRA4:
      case PDM_MESH_NODAL_PYRAMID5:
      case PDM_MESH_NODAL_PRISM6:
      case PDM_MESH_NODAL_HEXA8: {
        int *connec;
        PDM_Mesh_nodal_block_std_get (ml->mesh_nodal,
                                      id_block,
                                      ipart,
                                      &connec);
        if (parent_num != NULL) {
          int idx2 = 0;
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][parent_num[i]];
            int idx = ml->cell_vtx_idx[ipart][parent_num[i]];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx2++];
            }
          }
        }
        else {
          int idx2 = 0;
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][ielt];
            int idx = ml->cell_vtx_idx[ipart][ielt++];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx2++];
            }
          }
        }
        break;
      }
      case PDM_MESH_NODAL_POLY_2D: {
        int *connec_idx;
        int *connec;
        PDM_Mesh_nodal_block_poly2d_get (ml->mesh_nodal,
                                         id_block,
                                         ipart,
                                        &connec_idx,
                                        &connec);
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][parent_num[i]];
            int idx = ml->cell_vtx_idx[ipart][parent_num[i]];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][ielt];
            int idx = ml->cell_vtx_idx[ipart][ielt++];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        break;
      }
      case PDM_MESH_NODAL_POLY_3D:{
        int *connec_idx;
        int *connec;
        PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (ml->mesh_nodal,
                                                          id_block,
                                                          ipart,
                                                          &connec_idx,
                                                          &connec);
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][parent_num[i]];
            int idx = ml->cell_vtx_idx[ipart][parent_num[i]];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][ielt];
            int idx = ml->cell_vtx_idx[ipart][ielt++];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        break;
      }
      default :
        PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_location error : Bad Element type\n");
      }
    }
  }

  for (int ipart = 0; ipart < n_parts; ipart++) {
    free (n_vtx_per_elt[ipart]);
  }
  free (n_vtx_per_elt);


  PDM_MPI_Barrier (ml->comm);
  PDM_timer_hang_on(ml->timer);
  e_t_elapsed = PDM_timer_elapsed(ml->timer);
  e_t_cpu     = PDM_timer_cpu(ml->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

  ml->times_elapsed[BUILD_BOUNDING_BOXES] += e_t_elapsed - b_t_elapsed;
  ml->times_cpu[BUILD_BOUNDING_BOXES]     += e_t_cpu - b_t_cpu;
  ml->times_cpu_u[BUILD_BOUNDING_BOXES]   += e_t_cpu_u - b_t_cpu_u;
  ml->times_cpu_s[BUILD_BOUNDING_BOXES]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;
  PDM_timer_resume(ml->timer);

  /*
   *  Location
   *  --------
   */

  PDM_dbbtree_t *dbbt = NULL;
  PDM_box_set_t *box_set = NULL;
  if (ml->method == PDM_MESH_LOCATION_DBBTREE) {

    /* Compute local extents */
    double l_extents[6] = { HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                           -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    for (int i = 0; i < n_select_boxes; i++) {
      for (int j = 0; j < 3; j++) {
        l_extents[j]   = PDM_MIN (l_extents[j],   select_box_extents[6*i + j]);
        l_extents[j+3] = PDM_MAX (l_extents[j+3], select_box_extents[6*i + 3 + j]);
      }
    }

    /* Compute global extents */
    double g_extents[6];
    PDM_MPI_Allreduce (l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, ml->comm);
    PDM_MPI_Allreduce (l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

    /* Break symmetry */
    double max_range = 0.;
    for (int i = 0; i < 3; i++) {
      max_range = PDM_MAX (max_range, g_extents[i+3] - g_extents[i]);
    }
    for (int i = 0; i < 3; i++) {
      g_extents[i]   -= max_range * 1.1e-3;
      g_extents[i+3] += max_range * 1.0e-3;
    }

    dbbt = PDM_dbbtree_create (ml->comm, dim, g_extents);

    box_set = PDM_dbbtree_boxes_set (dbbt,
                                     1,
                                     &n_select_boxes,
                                     (const double **) &select_box_extents,
                                     (const PDM_g_num_t **) &select_box_g_num);

  }


  PDM_MPI_Barrier (ml->comm);
  PDM_timer_hang_on(ml->timer);
  e_t_elapsed = PDM_timer_elapsed(ml->timer);
  e_t_cpu     = PDM_timer_cpu(ml->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

  ml->times_elapsed[SEARCH_CANDIDATES] += e_t_elapsed - b_t_elapsed;
  ml->times_cpu[SEARCH_CANDIDATES]     += e_t_cpu - b_t_cpu;
  ml->times_cpu_u[SEARCH_CANDIDATES]   += e_t_cpu_u - b_t_cpu_u;
  ml->times_cpu_s[SEARCH_CANDIDATES]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;
  PDM_timer_resume(ml->timer);


  /*
   * Locate points
   */

  int         *pts_idx   = NULL;
  PDM_g_num_t *pts_g_num = NULL;
  double      *pts_coord = NULL;

  if (ml->reverse_result) {
    ml->points_in_elements = malloc (sizeof(_points_in_element_t) * ml->n_point_cloud);
  }

  for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    b_t_elapsed = PDM_timer_elapsed(ml->timer);
    b_t_cpu     = PDM_timer_cpu(ml->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);
    PDM_timer_resume(ml->timer);

    int loc_use_extracted_mesh = use_extracted_mesh;

    _point_cloud_t *pcloud = ml->point_clouds + icloud;

    /*
     * Concatenate point cloud partitions
     */

    int n_pts_pcloud = 0;
    PDM_g_num_t *pcloud_parent_g_num = NULL;
    PDM_g_num_t *pcloud_g_num = NULL;
    double      *pcloud_coord = NULL;

    if (allow_extraction) {
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        n_pts_pcloud += n_select_pts[icloud][ipart];
      }
      pcloud_g_num = malloc (sizeof(PDM_g_num_t) * n_pts_pcloud);
      pcloud_coord = malloc (sizeof(double)      * n_pts_pcloud * dim);
      int idx = 0;
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        for (int ipt = 0; ipt < n_select_pts[icloud][ipart]; ipt++) {

          pcloud_g_num[idx] = select_pts_g_num[icloud][ipart][ipt];

          for (int idim = 0; idim < dim; idim++) {
            pcloud_coord[dim*idx + idim] = select_pts_coord[icloud][ipart][dim*ipt + idim];
          }

          idx++;
        }
      }

      if (use_extracted_pts[icloud]) {
        pcloud_parent_g_num = malloc (sizeof(PDM_g_num_t) * n_pts_pcloud);
        idx = 0;
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          for (int ipt = 0; ipt < n_select_pts[icloud][ipart]; ipt++) {
            pcloud_parent_g_num[idx++] = select_pts_parent_g_num[icloud][ipart][ipt];
          }

          free (select_pts_parent_g_num[icloud][ipart]);
          free (select_pts_g_num[icloud][ipart]);
          free (select_pts_coord[icloud][ipart]);
        }
        free (n_select_pts[icloud]);
        free (select_pts_parent_g_num[icloud]);
        free (select_pts_g_num[icloud]);
        free (select_pts_coord[icloud]);
      }
      else {
        pcloud_parent_g_num = pcloud_g_num;
      }
    }

    else {
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        n_pts_pcloud += pcloud->n_points[ipart];
      }
      pcloud_g_num = malloc (sizeof(PDM_g_num_t) * n_pts_pcloud);
      pcloud_coord = malloc (sizeof(double)      * n_pts_pcloud * dim);
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
      pcloud_parent_g_num = pcloud_g_num;
    }

    if (0) {
      PDM_g_num_t ln_pts = (PDM_g_num_t) n_pts_pcloud;
      if (0) printf("[%6d] n_pts_pcloud = %d, ln_pts = "PDM_FMT_G_NUM"\n", my_rank, n_pts_pcloud, ln_pts);
      PDM_g_num_t gn_pts;
      PDM_MPI_Allreduce (&ln_pts, &gn_pts, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ml->comm);
      if (my_rank == 0) {
        printf("cloud %d : gn_pts = "PDM_FMT_G_NUM"\n", icloud, gn_pts);
      }
    }

    /*
     * Get points inside bounding boxes of elements
     */
    double t1all = PDM_MPI_Wtime();
    char *env_var_oct = getenv ("OCTREE_SHARED");
    int use_shared_octree = 0;
    if (env_var_oct != NULL) {
      use_shared_octree = atoi(env_var_oct);
    }

    switch (ml->method) {

    case PDM_MESH_LOCATION_OCTREE: {
      /* Create octree structure */
      octree = PDM_para_octree_create (1,
                                       octree_depth_max,
                                       octree_points_in_leaf_max,
                                       octree_build_leaf_neighbours,
                                       ml->comm);

      /* Set octree point cloud */
      PDM_para_octree_point_cloud_set (octree,
                                       0,
                                       n_pts_pcloud,
                                       pcloud_coord,
                                       pcloud_g_num);

      /* Build parallel octree */
      PDM_MPI_Barrier(ml->comm);
      double t1 = PDM_MPI_Wtime();

      if(use_shared_octree == 0) {
        PDM_para_octree_build (octree, NULL);
      } else {
        PDM_para_octree_build_shared (octree, NULL);
      }
      if (1) {
        end_timer_and_print("PDM_para_octree_build ", ml->comm, t1);
      }
      // PDM_para_octree_dump (octree);
      // if (dbg_enabled) {
        // PDM_para_octree_dump_times (octree);
      // }
      // PDM_para_octree_export_vtk(octree,"octree_debug_ml_");

      /* Locate points inside boxes */
      // PDM_MPI_Barrier (ml->comm);
      t1 = PDM_MPI_Wtime();
      if(use_shared_octree == 0) {
        PDM_para_octree_points_inside_boxes (octree,
                                             n_select_boxes,
                                             select_box_extents,
                                             select_box_g_num,
                                             &pts_idx,
                                             &pts_g_num,
                                             &pts_coord);
      } else {
       PDM_para_octree_points_inside_boxes_shared (octree,
                                                   n_select_boxes,
                                                   select_box_extents,
                                                   select_box_g_num,
                                                   &pts_idx,
                                                   &pts_g_num,
                                                   &pts_coord);
      }
      if (1) {
        end_timer_and_print("PDM_para_octree_points_inside_boxes ", ml->comm, t1);
      }

      t1 = PDM_MPI_Wtime();
      /* Free octree */
      PDM_para_octree_free (octree);
      if (1) {
        end_timer_and_print("PDM_para_octree_free ", ml->comm, t1);
      }
      break;
     }
    case PDM_MESH_LOCATION_DBBTREE: {
      if (dbg_enabled) printf("[%d] n_pts_pcloud = %d, n_select_boxes = %d\n", my_rank, n_pts_pcloud, n_select_boxes);//
      // if (USE_OCTREE_COPIES) {
      double t1dbtree = PDM_MPI_Wtime();
      if(use_shared_octree == 0) {
        PDM_dbbtree_points_inside_boxes (dbbt,
                                         n_pts_pcloud,
                                         pcloud_g_num,
                                         pcloud_coord,
                                         n_select_boxes,
                                         select_box_g_num,
                                         &pts_idx,
                                         &pts_g_num,
                                         &pts_coord,
                                         0);
      } else {
        // PDM_MPI_Barrier (ml->comm);
        PDM_dbbtree_points_inside_boxes_shared(dbbt,
                                               n_pts_pcloud,
                                               pcloud_g_num,
                                               pcloud_coord,
                                               n_select_boxes,
                                               select_box_g_num,
                                               &pts_idx,
                                               &pts_g_num,
                                               &pts_coord,
                                               0);
      }
      if (0) {
        end_timer_and_print("PDM_dbbtree_points_inside_boxes ", ml->comm, t1dbtree);
      }
      break;
     }
    case PDM_MESH_LOCATION_DOCTREE: {

      PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_KDTREE;
      // PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_OCTREE;
      PDM_doctree_t *doct = PDM_doctree_create(ml->comm,
                                               3,
                                               1,
                                               NULL, // global_extents
                                               local_tree_kind);

      int *init_location_pts = NULL;
      if(0) {
        init_location_pts = malloc(3 * n_pts_pcloud * sizeof(int));
        for(int i = 0; i < n_pts_pcloud; ++i) {
          init_location_pts[3*i  ] = my_rank;
          init_location_pts[3*i+1] = 0;
          init_location_pts[3*i+2] = i;
        }
      }

      PDM_doctree_point_set(doct,
                            0,
                            n_pts_pcloud,
                            init_location_pts,
                            pcloud_g_num,
                            pcloud_coord);


      int *init_location_box = malloc(3 * n_select_boxes * sizeof(int));
      for(int i = 0; i < n_select_boxes; ++i) {
        init_location_box[3*i  ] = my_rank;
        init_location_box[3*i+1] = 0;
        init_location_box[3*i+2] = i;
      }

      PDM_doctree_solicitation_set(doct,
                                   PDM_TREE_SOLICITATION_BOXES_POINTS,
                                   1,
                                   &n_select_boxes,
                                   &init_location_box,
                                   &select_box_g_num,
                                   &select_box_extents);

      double t1 = PDM_MPI_Wtime();
      PDM_doctree_build(doct);
      if (log_timer) {
        log_trace("PDM_doctree_build : %12.5e \n", PDM_MPI_Wtime()-t1);
      }

      t1 = PDM_MPI_Wtime();
      PDM_doctree_results_in_orig_frame_get(doct,
                                            n_select_boxes,
                                            select_box_g_num,
                                            &pts_idx,
                                            &pts_g_num,
                                            &pts_coord);
      if (log_timer) {
        log_trace("PDM_doctree_results_in_orig_frame_get : %12.5e \n", PDM_MPI_Wtime()-t1);
      }

      // PDM_log_trace_connectivity_long(pts_idx, pts_g_num, n_select_boxes, "pts_g_num : ");
      // PDM_log_trace_array_long(select_box_g_num, n_select_boxes, "select_box_g_num : ");
      PDM_doctree_dump_times(doct);
      PDM_doctree_free(doct);
      if(init_location_pts != NULL) {
        free(init_location_pts);
      }
      free(init_location_box);

      break;
    }
    default:
      printf("Error: unknown location method %d\n", ml->method);
      assert (1 == 0);

    }
    free (pcloud_coord);

    if (dbg_enabled) {
      printf("n_select_boxes = %i \n", n_select_boxes);
      printf("\n[%d] --- Pts in box ---\n", my_rank);
      for (ibox = 0; ibox < n_select_boxes; ibox++) {

        if (pts_idx[ibox+1] <= pts_idx[ibox]) {
          continue;
        }

        printf("[%d] %d ("PDM_FMT_G_NUM"): ", my_rank, ibox, select_box_g_num[ibox]);
        //printf("[%d] %d ("PDM_FMT_G_NUM"): ", my_rank, ibox, select_box_parent_g_num[ibox]);
        for (int i = pts_idx[ibox]; i < pts_idx[ibox+1]; i++) {
          /*printf("((%ld); %f %f %f) ",
            pts_g_num[i], pts_coord[dim*i], pts_coord[dim*i+1], pts_coord[dim*i+2]);*/
          printf("("PDM_FMT_G_NUM") ", pts_g_num[i]);
        }
        printf("\n");
      }
      printf("[%d] ------------------\n\n\n", my_rank);
    }

    if (0) {
      for (ibox = 0; ibox < n_select_boxes; ibox++) {
        if (select_box_g_num[ibox] == 2793384) {
          printf("[%d] box "PDM_FMT_G_NUM", %d point(s) inside:", my_rank, select_box_g_num[ibox], pts_idx[ibox+1] - pts_idx[ibox]);
          for (int i = pts_idx[ibox]; i < pts_idx[ibox+1]; i++) {
            printf(" "PDM_FMT_G_NUM, pts_g_num[i]);
          }
          printf("\n");
        }
      }
    }
    end_timer_and_print("Search all  ", ml->comm, t1all);

    // PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[SEARCH_CANDIDATES] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu[SEARCH_CANDIDATES]     += e_t_cpu - b_t_cpu;
    ml->times_cpu_u[SEARCH_CANDIDATES]   += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s[SEARCH_CANDIDATES]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    /*
     * 2nd extraction: Remove elements that cannot contain any target point
     */
    int          **loc_n_select_elt            = NULL;
    int         ***loc_select_elt_l_num        = NULL;
    PDM_g_num_t   *loc_select_box_parent_g_num = NULL;
    PDM_g_num_t   *loc_select_box_g_num        = NULL;
    // printf("------------------------ icloud = %i loc_use_extracted_mesh = %i \n", icloud, use_extracted_mesh);
    if (loc_use_extracted_mesh) {

      loc_n_select_elt            = malloc(n_blocks       * sizeof(int  *     ));
      loc_select_elt_l_num        = malloc(n_blocks       * sizeof(int **     ));
      loc_select_box_parent_g_num = malloc(n_select_boxes * sizeof(PDM_g_num_t));
      loc_select_box_g_num        = malloc(n_select_boxes * sizeof(PDM_g_num_t));
      // printf("n_select_boxes = %i \n", n_select_boxes);

      int idx1 = 0;
      ibox = 0;
      for (int iblock = 0; iblock < n_blocks; iblock++) {

        loc_n_select_elt    [iblock] = malloc (n_parts * sizeof(int  ));
        loc_select_elt_l_num[iblock] = malloc (n_parts * sizeof(int *));

        // PDM_log_trace_array_long(select_box_parent_g_num, n_select_boxes, "select_box_parent_g_num :: ");

        for (int ipart = 0; ipart < n_parts; ipart++) {

          loc_n_select_elt    [iblock][ipart] = 0;
          loc_select_elt_l_num[iblock][ipart] = malloc (sizeof(int) * n_select_elt[iblock][ipart]);

          int idx2 = 0;
          // printf("n_select_elt[iblock][ipart] = %i \n", n_select_elt[iblock][ipart]);
          for (int ielt = 0; ielt < n_select_elt[iblock][ipart]; ielt++) {
            // printf("ielt = %i | ibox = %i - %i %i \n", ielt, ibox, pts_idx[ibox], pts_idx[ibox+1]);
            if (pts_idx[ibox] < pts_idx[ibox+1]) {
              loc_select_elt_l_num[iblock][ipart][idx2] = select_elt_l_num[iblock][ipart][ielt];

              loc_select_box_parent_g_num[idx1] = select_box_parent_g_num[ibox];
              loc_select_box_g_num       [idx1] = select_box_g_num       [ibox];
              pts_idx[idx1+1] = pts_idx[ibox+1];
              idx2++;
              idx1++;
            }
            ibox++;
          }
          // realloc?
          loc_n_select_elt[iblock][ipart] = idx2;

          loc_select_elt_l_num[iblock][ipart] = realloc(loc_select_elt_l_num[iblock][ipart], loc_n_select_elt[iblock][ipart] * sizeof(int));

        } // End of loop on parts
      } // End of loop on nodal blocks

      if (dbg_enabled) printf("[%4d] before : %8d, after : %8d\n", my_rank, n_select_boxes, idx1);

      /* Realloc */
      loc_select_box_parent_g_num = realloc(loc_select_box_parent_g_num, idx1 * sizeof(PDM_g_num_t));
      loc_select_box_g_num        = realloc(loc_select_box_g_num       , idx1 * sizeof(PDM_g_num_t));
    }

    else if (allow_extraction) {
      loc_use_extracted_mesh = 1;

      loc_n_select_elt            = malloc(n_blocks       * sizeof(int  *     ));
      loc_select_elt_l_num        = malloc(n_blocks       * sizeof(int **     ));
      loc_select_box_parent_g_num = malloc(n_select_boxes * sizeof(PDM_g_num_t));

      int idx1 = 0;
      ibox = 0;
      for (int iblock = 0; iblock < n_blocks; iblock++) {

        loc_n_select_elt    [iblock] = malloc (n_parts * sizeof(int  ));
        loc_select_elt_l_num[iblock] = malloc (n_parts * sizeof(int *));

        int id_block = blocks_id[iblock];
        for (int ipart = 0; ipart < n_parts; ipart++) {
          int part_n_elt = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                           id_block,
                                                           ipart);

          loc_select_elt_l_num[iblock][ipart] = malloc (sizeof(int) * part_n_elt);

          int idx2 = 0;
          for (int ielt = 0; ielt < part_n_elt; ielt++) {
            if (pts_idx[ibox] < pts_idx[ibox+1]) {
              loc_select_elt_l_num[iblock][ipart][idx2] = ielt;
              loc_select_box_parent_g_num[idx1] = box_g_num[ibox];
              pts_idx[idx1+1] = pts_idx[ibox+1];
              idx2++;
              idx1++;
            }
            ibox++;
          }
          // realloc?
          loc_n_select_elt[iblock][ipart] = idx2;

          loc_select_elt_l_num[iblock][ipart] = realloc(loc_select_elt_l_num[iblock][ipart], idx2 * sizeof(int));

        } // End of loop on parts
      } // End of loop on nodal blocks

      if (dbg_enabled) printf("[%4d] before : %8d, after : %8d\n", my_rank, n_select_boxes, idx1);
      int loc_n_select_boxes = idx1;

      loc_select_box_parent_g_num = realloc(loc_select_box_parent_g_num, idx1 * sizeof(PDM_g_num_t));

      /*
       *  Generate a new global numbering for selected boxes
       */
      PDM_gen_gnum_t *gen_gnum_boxes = PDM_gnum_create (3,
                                                        1,
                                                        PDM_FALSE,
                                                        0.,
                                                        ml->comm,
                                                        PDM_OWNERSHIP_USER);
      PDM_gnum_set_from_parents (gen_gnum_boxes,
                                 0,
                                 loc_n_select_boxes,
                                 loc_select_box_parent_g_num);

      PDM_gnum_compute (gen_gnum_boxes);

      loc_select_box_g_num = PDM_gnum_get (gen_gnum_boxes, 0);

      PDM_gnum_free (gen_gnum_boxes);
    } /* End if extracted_mesh */


    /*
     * Load balancing: redistribute evenly elementary location operations
     */
    int redistrib_n_elt = 0;

    PDM_g_num_t *redistrib_elt_parent_g_num = NULL;
    PDM_g_num_t *redistrib_elt_g_num        = NULL;
    int         *redistrib_vtx_idx          = NULL;
    double      *redistrib_vtx_coord        = NULL;
    int         *redistrib_pts_idx          = NULL;
    PDM_g_num_t *redistrib_pts_parent_g_num = NULL;
    PDM_g_num_t *redistrib_pts_g_num        = NULL;
    double      *redistrib_pts_coord        = NULL;
    PDM_l_num_t *redistrib_poly3d_face_idx  = NULL;
    PDM_l_num_t *redistrib_face_vtx_idx     = NULL;
    PDM_l_num_t *redistrib_face_vtx         = NULL;
    int         *redistrib_face_orientation = NULL;

    int redistrib_type_idx[PDM_MESH_NODAL_N_ELEMENT_TYPES + 1];

    if (loc_use_extracted_mesh) {
      _extract_selected_mesh_elements (ml,
                                       loc_n_select_elt,
                                       loc_select_elt_l_num,
                                       loc_select_box_parent_g_num,
                                       loc_select_box_g_num,
                                       pts_idx,
                                       pts_g_num,
                                       pts_coord,
                                       &redistrib_n_elt,
                                       redistrib_type_idx,
                                       &redistrib_elt_parent_g_num,
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
      // free (select_box_parent_g_num);
    }
    else {
      _redistribute_elementary_location (ml,
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
    }
    free (pts_idx);
    free (pts_g_num);
    free (pts_coord);

    if (allow_extraction) {
      for (int iblock = 0; iblock < n_blocks; iblock++) {
        for (int ipart = 0; ipart < n_parts; ipart++) {
          free (loc_select_elt_l_num[iblock][ipart]);
        }
        free (loc_n_select_elt[iblock]);
        free (loc_select_elt_l_num[iblock]);
      }
      free (loc_n_select_elt);
      free (loc_select_elt_l_num);

      free(loc_select_box_parent_g_num);
      free(loc_select_box_g_num);
    }


    int n_pts = redistrib_pts_idx[redistrib_n_elt];

    if (use_extracted_pts[icloud]) {
      if (0) {
        printf("[%d] redistrib_pts_g_num = ", my_rank);
        for (int i = 0; i < n_pts; i++) {
          printf(PDM_FMT_G_NUM" ", redistrib_pts_g_num[i]);
        }
        printf("\n");
      }

      /*for (int i = 0; i < redistrib_n_elt; i++) {
        if (redistrib_elt_parent_g_num[i] == 657) {
          printf("elt ("PDM_FMT_G_NUM") : %d points\n", redistrib_elt_parent_g_num[i], redistrib_pts_idx[i+1] - redistrib_pts_idx[i]);
        }
        }*/

      // substitute redistrib_pts_g_num with redistrib_pts_parent_g_num
      PDM_part_to_block_t *ptb_parent =
        PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                  PDM_PART_TO_BLOCK_POST_CLEANUP,
                                  1.,
                                  &pcloud_g_num,
                                  NULL,
                                  &n_pts_pcloud,
                                  1,
                                  ml->comm);
      PDM_g_num_t *block_parent_distrib_idx = PDM_part_to_block_distrib_index_get (ptb_parent);

      int *_block_stride = NULL;
      PDM_g_num_t *block_pcloud_parent_gnum = NULL;
      PDM_part_to_block_exch (ptb_parent,
                              sizeof(PDM_g_num_t),
                              PDM_STRIDE_CST_INTERLACED,
                              1,
                              NULL,
                              (void **) &pcloud_parent_g_num,
                              &_block_stride,
                              (void **) &block_pcloud_parent_gnum);
      free (pcloud_parent_g_num);

      if (0) {
        int n_pts_a = PDM_part_to_block_n_elt_block_get (ptb_parent);
        PDM_g_num_t *block_g_num_a = PDM_part_to_block_block_gnum_get (ptb_parent);
        printf("[%d] block_pcloud_parent_gnum :\n", my_rank);
        for (int i = 0; i < n_pts_a; i++) {
          printf("  "PDM_FMT_G_NUM" --> "PDM_FMT_G_NUM"\n", block_g_num_a[i], block_pcloud_parent_gnum[i]);
        }
      }

      PDM_block_to_part_t *btp_parent =
        PDM_block_to_part_create (block_parent_distrib_idx,
                                  (const PDM_g_num_t **) &redistrib_pts_g_num,
                                  &n_pts,
                                  1,
                                  ml->comm);

      int _one = 1;
      redistrib_pts_parent_g_num = malloc (sizeof(PDM_g_num_t) * n_pts);
      PDM_block_to_part_exch_in_place (btp_parent,
                              sizeof(PDM_g_num_t),
                              PDM_STRIDE_CST_INTERLACED,
                              &_one,
                              (void *) block_pcloud_parent_gnum,
                              NULL,
                              (void **) &redistrib_pts_parent_g_num);
      free (block_pcloud_parent_gnum);
      if (0) {
        printf("[%d] redistrib_pts_parent_g_num = ", my_rank);
        for (int i = 0; i < n_pts; i++) {
          printf(PDM_FMT_G_NUM" ", redistrib_pts_parent_g_num[i]);
        }
        printf("\n");
      }

      ptb_parent = PDM_part_to_block_free (ptb_parent);
      btp_parent = PDM_block_to_part_free (btp_parent);
    }
    else {
      redistrib_pts_parent_g_num = redistrib_pts_g_num;
    }

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[LOAD_BALANCING] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu[LOAD_BALANCING]     += e_t_cpu - b_t_cpu;
    ml->times_cpu_u[LOAD_BALANCING]   += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s[LOAD_BALANCING]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);


    /*
     * Location in elements
     */
    PDM_g_num_t *redistrib_pts_location = malloc (sizeof(PDM_g_num_t) * n_pts);
    if (loc_use_extracted_mesh) {
      for (ibox = 0; ibox < redistrib_n_elt; ibox++) {
        for (int i = redistrib_pts_idx[ibox]; i < redistrib_pts_idx[ibox+1]; i++) {
          redistrib_pts_location[i] = redistrib_elt_parent_g_num[ibox];
        }
      }
      // free (redistrib_elt_parent_g_num);
    }
    else {
      for (ibox = 0; ibox < redistrib_n_elt; ibox++) {
        for (int i = redistrib_pts_idx[ibox]; i < redistrib_pts_idx[ibox+1]; i++) {
          redistrib_pts_location[i] = redistrib_elt_g_num[ibox];
        }
      }
    }
    // free (redistrib_elt_g_num);

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

    // if (dbg_enabled) {
    //   for (int k = 0; k < redistrib_type_idx[PDM_MESH_NODAL_N_ELEMENT_TYPES]; k++) {
    //     if (use_extracted_mesh) {
    //       log_trace("Elt gnum = "PDM_FMT_G_NUM"\n", redistrib_elt_parent_g_num[k]);
    //     } else {
    //       log_trace("Elt gnum = "PDM_FMT_G_NUM"\n", redistrib_elt_g_num[k]);
    //     }
    //     for (int i = redistrib_pts_idx[k]; i < redistrib_pts_idx[k+1]; i++) {
    //       log_trace("Point gnum = ("PDM_FMT_G_NUM")\n", redistrib_pts_g_num[i]);
    //       log_trace("\t        coords = (%f, %f, %f)\n",
    //                 redistrib_pts_coord[dim*i],
    //                 redistrib_pts_coord[dim*i+1],
    //                 redistrib_pts_coord[dim*i+2]);
    //       log_trace("\tlocation = ("PDM_FMT_G_NUM")\n", redistrib_pts_location[i]);
    //       log_trace("\t  proj. coords = (%f, %f, %f)\n",
    //                 projected_coord[dim*i],
    //                 projected_coord[dim*i+1],
    //                 projected_coord[dim*i+2]);
    //       log_trace("\tdistance = %f\n", distance[i]);
    //       log_trace("\t weights =");
    //       for (int j = weights_idx[i]; j < weights_idx[i+1]; j++) {
    //         log_trace(" %f", (float) weights[j]);
    //       }
    //       log_trace("\n");
    //     }
    //   }
    // }
    free (redistrib_elt_g_num);
    free (redistrib_vtx_coord);
    free (redistrib_vtx_idx);
    free (redistrib_pts_coord);
    free (redistrib_pts_idx);
    free (redistrib_poly3d_face_idx);
    free (redistrib_face_vtx_idx);
    free (redistrib_face_vtx);
    free (redistrib_face_orientation);
    if (loc_use_extracted_mesh) {
      free (redistrib_elt_parent_g_num);
    }

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[COMPUTE_ELEMENTARY_LOCATIONS] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu[COMPUTE_ELEMENTARY_LOCATIONS]     += e_t_cpu - b_t_cpu;
    ml->times_cpu_u[COMPUTE_ELEMENTARY_LOCATIONS]   += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s[COMPUTE_ELEMENTARY_LOCATIONS]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    /*
     * Merge location data
     */
    /*
     *   1) Part-to-block
     */
    PDM_g_num_t *block_parent_distrib_idx =
      PDM_compute_uniform_entity_distribution_from_partition (ml->comm,
                                                              pcloud->n_part,
                                                              pcloud->n_points,
                                                              (const PDM_g_num_t **) pcloud->gnum);

    PDM_part_to_block_t *ptb1 = PDM_part_to_block_create_from_distrib (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           &redistrib_pts_parent_g_num,
                                                           block_parent_distrib_idx,
                                                           &n_pts,
                                                           1,
                                                           ml->comm);
    PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb1);
    if (redistrib_pts_parent_g_num != redistrib_pts_g_num) {
      free (redistrib_pts_parent_g_num);
    }
    free (redistrib_pts_g_num);

    int *part_stride = PDM_array_const_int(n_pts, 1);

    int *block_stride = NULL;

    /* Exchange location (synchronous) */
    PDM_g_num_t *block_location1 = NULL;
    PDM_part_to_block_exch (ptb1,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
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
                            PDM_STRIDE_VAR_INTERLACED,
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
                            PDM_STRIDE_VAR_INTERLACED,
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
                            PDM_STRIDE_VAR_INTERLACED,
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
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            &part_stride,
                            (void **) &weights_stride,
                            &block_n_candidates,
                            (void **) &block_n_vtx_elt);
    free (weights_stride);

    /* Exchange projected coords */
    /*for (int i = 0; i < n_pts; i++) {
      part_stride[i] = 3;
      }*/
    double *block_proj_coord1 = NULL;
    PDM_part_to_block_exch (ptb1,
                            3*sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            &part_stride,
                            (void **) &projected_coord,
                            &block_stride,
                            (void **) &block_proj_coord1);
    free (projected_coord);
    free (block_stride);
    free (part_stride);

    /* Exchange location ( Async )*/
    // PDM_g_num_t *block_location1 = NULL;
    // int id1 = PDM_part_to_block_async_exch (ptb1,
    //                         sizeof(PDM_g_num_t),
    //                         PDM_STRIDE_VAR,
    //                         1,
    //                         &part_stride,
    //                         (void **) &redistrib_pts_location);
    // free (redistrib_pts_location);

    // /* Exchange element type */
    // PDM_Mesh_nodal_elt_t *block_elt_type = NULL;
    // int id2 = PDM_part_to_block_async_exch (ptb1,
    //                         sizeof(PDM_Mesh_nodal_elt_t),
    //                         PDM_STRIDE_VAR,
    //                         1,
    //                         &part_stride,
    //                         (void **) &redistrib_pts_elt_type);

    // free (redistrib_pts_elt_type);

    // /* Exchange distance */
    // double *block_distance = NULL;
    // int id3 = PDM_part_to_block_async_exch (ptb1,
    //                         sizeof(double),
    //                         PDM_STRIDE_VAR,
    //                         1,
    //                         &part_stride,
    //                         (void **) &distance);
    // free (distance);

    // /* Exchange weights */
    // int *weights_stride = malloc (sizeof(int) * n_pts);
    // for (int i = 0; i < n_pts; i++) {
    //   weights_stride[i] = weights_idx[i+1] - weights_idx[i];
    // }
    // free (weights_idx);

    // double *block_weights1 = NULL;
    // int id4 = PDM_part_to_block_async_exch (ptb1,
    //                         sizeof(double),
    //                         PDM_STRIDE_VAR,
    //                         1,
    //                         &weights_stride,
    //                         (void **) &weights);
    // free (weights);

    // /* Exchange weights stride */
    // int *block_n_vtx_elt    = NULL;
    // int *block_n_candidates = NULL;
    // int id5 = PDM_part_to_block_async_exch (ptb1,
    //                         sizeof(int),
    //                         PDM_STRIDE_VAR,
    //                         1,
    //                         &part_stride,
    //                         (void **) &weights_stride);
    // free (weights_stride);

    // /* Exchange projected coords */
    // for (int i = 0; i < n_pts; i++) {
    //   part_stride[i] = 3;
    // }
    // double *block_proj_coord1 = NULL;
    // int id6 = PDM_part_to_block_async_exch (ptb1,
    //                         sizeof(double),
    //                         PDM_STRIDE_VAR,
    //                         1,
    //                         &part_stride,
    //                         (void **) &projected_coord);
    // free (projected_coord);
    // free (part_stride);


    // PDM_part_to_block_async_wait(ptb1, id1);
    // PDM_part_to_block_asyn_post_treatment(ptb1, id1, &block_stride, (void **) &block_location1);
    // free (block_stride);

    // PDM_part_to_block_async_wait(ptb1, id2);
    // PDM_part_to_block_asyn_post_treatment(ptb1, id2, &block_stride, (void **) &block_elt_type);
    // free (block_stride);

    // PDM_part_to_block_async_wait(ptb1, id3);
    // PDM_part_to_block_asyn_post_treatment(ptb1, id3, &block_stride, (void **) &block_distance);
    // free (block_stride);

    // PDM_part_to_block_async_wait(ptb1, id4);
    // PDM_part_to_block_asyn_post_treatment(ptb1, id4, &block_stride, (void **) &block_weights1);
    // free (block_stride);

    // PDM_part_to_block_async_wait(ptb1, id5);
    // PDM_part_to_block_asyn_post_treatment(ptb1, id5, &block_n_candidates, (void **) &block_n_vtx_elt);

    // PDM_part_to_block_async_wait(ptb1, id6);
    // PDM_part_to_block_asyn_post_treatment(ptb1, id6, &block_stride, (void **) &block_proj_coord1);
    // free (block_stride);

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[MERGE_LOCATION_DATA_STEP1] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu[MERGE_LOCATION_DATA_STEP1]     += e_t_cpu - b_t_cpu;
    ml->times_cpu_u[MERGE_LOCATION_DATA_STEP1]   += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s[MERGE_LOCATION_DATA_STEP1]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);


    /*
     *   2) Among candidate elements, keep closest one for each point (set location to -1 if no candidate -> unlocated point)
     */
    PDM_g_num_t *block_g_num1 = PDM_part_to_block_block_gnum_get (ptb1);

    const int n_pts_block1 = PDM_part_to_block_n_elt_block_get (ptb1);

    int n_pts_block2 = (int) (block_parent_distrib_idx[my_rank+1] -
                              block_parent_distrib_idx[my_rank]);

    PDM_g_num_t *block_location2 = malloc (sizeof(PDM_g_num_t) * n_pts_block2);
    double *block_proj_coord2 = malloc (sizeof(double) * n_pts_block2 * 3);
    double *block_dist2 = malloc (sizeof(double) * n_pts_block2);

    int *block_weights_stride2 = malloc (sizeof(int) * n_pts_block2);
    for (int i = 0; i < n_pts_block2; i++) {
      block_location2[i] = -1;
      block_dist2[i] = -1.;
      block_weights_stride2[i] = 0;
    }
    for (int i = 0; i < 3*n_pts_block2; i++) {
      block_proj_coord2[i] = 0.;
    }

    int *idx_min = malloc (sizeof(int) * n_pts_block2);
    int idx = 0;
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

        }
      }

      idx += block_n_candidates[i];

      int ipt = block_g_num1[i] - 1 - block_distrib_idx[my_rank];

      block_location2[ipt] = block_location1[idx_min[i]];
      block_weights_stride2[ipt] = block_n_vtx_elt[idx_min[i]];
      block_dist2[ipt] = block_distance[idx_min[i]];

      for (int j = 0; j < 3; j++) {
        block_proj_coord2[3*ipt + j] = block_proj_coord1[3*idx_min[i] + j];
      }
      s_weights += block_weights_stride2[ipt];
    }

    int *block_weights_idx1 = PDM_array_new_idx_from_sizes_int(block_n_vtx_elt, idx);


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

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[MERGE_LOCATION_DATA_STEP2] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu[MERGE_LOCATION_DATA_STEP2]     += e_t_cpu - b_t_cpu;
    ml->times_cpu_u[MERGE_LOCATION_DATA_STEP2]   += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s[MERGE_LOCATION_DATA_STEP2]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    /*
     *   3) Block-to-part
     */
    int one = 1;
    int three = 3;
    PDM_block_to_part_t *btp = PDM_block_to_part_create (block_parent_distrib_idx,
                                                         (const PDM_g_num_t **) pcloud->gnum,
                                                         pcloud->n_points,
                                                         pcloud->n_part,
                                                         ml->comm);
    free (pcloud_g_num);

    pcloud->location         = malloc (sizeof(PDM_g_num_t *) * pcloud->n_part);
    pcloud->dist2            = malloc (sizeof(double *)      * pcloud->n_part);
    pcloud->weights_idx      = malloc (sizeof(int *)         * pcloud->n_part);
    pcloud->weights          = malloc (sizeof(double *)      * pcloud->n_part);
    pcloud->projected_coords = malloc (sizeof(double *)      * pcloud->n_part);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      pcloud->location[ipart]         = malloc (sizeof(PDM_g_num_t) * pcloud->n_points[ipart]);
      pcloud->dist2[ipart]            = malloc (sizeof(double) * pcloud->n_points[ipart]);
      pcloud->projected_coords[ipart] = malloc (sizeof(double) * pcloud->n_points[ipart] * 3);
      pcloud->weights_idx[ipart]      = malloc (sizeof(int) * (pcloud->n_points[ipart] + 1));
    }

    /* Exchange location */
    PDM_block_to_part_exch_in_place (btp,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST_INTERLACED,
                            &one,
                            block_location2,
                            NULL,
                            (void **) pcloud->location);
    free (block_location2);

    /* Exchange distance */
    PDM_block_to_part_exch_in_place (btp,
                            sizeof(double),
                            PDM_STRIDE_CST_INTERLACED,
                            &one,
                            block_dist2,
                            NULL,
                            (void **) pcloud->dist2);
    free (block_dist2);

    /* Exchange projected coords */
    PDM_block_to_part_exch_in_place (btp,
                            sizeof(double),
                            PDM_STRIDE_CST_INTERLACED,
                            &three,
                            block_proj_coord2,
                            NULL,
                            (void **) pcloud->projected_coords);
    free (block_proj_coord2);

    /* Exchange weights stride */
    int **_weights_stride = malloc (sizeof(int *) * pcloud->n_part);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      _weights_stride[ipart] = malloc (sizeof(int) * pcloud->n_points[ipart]);
    }
    PDM_block_to_part_exch_in_place (btp,
                            sizeof(int),
                            PDM_STRIDE_CST_INTERLACED,
                            &one,
                            (void *) block_weights_stride2,
                            NULL,
                            (void **) _weights_stride);

    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      PDM_array_idx_from_sizes_int (_weights_stride[ipart],
                                    pcloud->n_points[ipart],
                                    pcloud->weights_idx[ipart]);
      pcloud->weights[ipart] =
        malloc (sizeof(double) * pcloud->weights_idx[ipart][pcloud->n_points[ipart]]);
    }

    /* Exchange weights */
    PDM_block_to_part_exch_in_place (btp,
                            sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            block_weights_stride2,
                            (void *) block_weights2,
                            _weights_stride,
                            (void **) pcloud->weights);
    free (block_weights2);
    free (block_weights_stride2);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      free (_weights_stride[ipart]);
    }
    free (_weights_stride);

    free (block_parent_distrib_idx);
    PDM_part_to_block_free (ptb1);
    PDM_block_to_part_free (btp);

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[MERGE_LOCATION_DATA_STEP3] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu[MERGE_LOCATION_DATA_STEP3]     += e_t_cpu - b_t_cpu;
    ml->times_cpu_u[MERGE_LOCATION_DATA_STEP3]   += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s[MERGE_LOCATION_DATA_STEP3]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    /*
     *  Compress results : Get points in elements
     *  ----------------------------------------
     */

    pcloud->n_located    = malloc (sizeof(int) * pcloud->n_part);
    pcloud->n_un_located = malloc (sizeof(int) * pcloud->n_part);
    pcloud->located      = malloc (sizeof(int *) * pcloud->n_part);
    pcloud->un_located   = malloc (sizeof(int *) * pcloud->n_part);

    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      pcloud->n_located[ipart]    = 0;
      pcloud->n_un_located[ipart] = 0;

      pcloud->located[ipart]    = malloc (sizeof(int) * pcloud->n_points[ipart]);
      pcloud->un_located[ipart] = malloc (sizeof(int) * pcloud->n_points[ipart]);

      for (int j = 0; j < pcloud->n_points[ipart]; j++) {
        if (pcloud->location[ipart][j] > 0) {
          int idx2 = pcloud->n_located[ipart];

          pcloud->location[ipart][idx2] = pcloud->location[ipart][j];
          pcloud->dist2[ipart][idx2]    = pcloud->dist2[ipart][j];
          for (int k = 0; k < 3; ++k) {
            pcloud->projected_coords[ipart][3*idx2+k] = pcloud->projected_coords[ipart][3*j+k];
          }

          int n_weights =  pcloud->weights_idx[ipart][j+1] - pcloud->weights_idx[ipart][j];

          for (int k = 0; k < n_weights; k++) {
            pcloud->weights[ipart][pcloud->weights_idx[ipart][idx2] +k] = pcloud->weights[ipart][pcloud->weights_idx[ipart][j] +k];   
          }

          pcloud->weights_idx[ipart][idx2+1] = pcloud->weights_idx[ipart][idx2] + n_weights; 

          pcloud->located[ipart][idx2] = j+1;
          pcloud->n_located[ipart]++;
        }
        else {
          pcloud->un_located[ipart][pcloud->n_un_located[ipart]++] = j+1;
        }
      }

      pcloud->located[ipart]    = realloc (pcloud->located[ipart],
                                          sizeof(int) * pcloud->n_located[ipart]);
      pcloud->un_located[ipart] = realloc (pcloud->un_located[ipart],
                                          sizeof(int) * pcloud->n_un_located[ipart]);

      pcloud->location[ipart]   = realloc (pcloud->location[ipart],
                                          sizeof(PDM_g_num_t) * pcloud->n_located[ipart]);
      pcloud->dist2[ipart]      = realloc (pcloud->dist2[ipart],
                                          sizeof(double) * pcloud->n_located[ipart]);

      pcloud->weights_idx[ipart] = realloc (pcloud->weights_idx[ipart],
                                           sizeof(int) * (pcloud->n_located[ipart]+1));

      pcloud->weights[ipart]     = realloc (pcloud->weights[ipart],
                                            sizeof(double) * pcloud->weights_idx[ipart][pcloud->n_located[ipart]]);

      pcloud->projected_coords[ipart] = realloc (pcloud->projected_coords[ipart],
                                          sizeof(double) * 3 * pcloud->n_located[ipart]);

    }

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[COMPRESS_LOCATION_DATA] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu[COMPRESS_LOCATION_DATA]     += e_t_cpu - b_t_cpu;
    ml->times_cpu_u[COMPRESS_LOCATION_DATA]   += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s[COMPRESS_LOCATION_DATA]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    if (ml->reverse_result) {

      /*
       *  Reverse results : Get points in elements
       *  ----------------------------------------
       */

      double      **elt_weight     = malloc (sizeof(double      *) * pcloud->n_part);
      int         **ptb_stride     = malloc (sizeof(int         *) * pcloud->n_part);
      int         **weight_stride  = malloc (sizeof(int         *) * pcloud->n_part);
      PDM_g_num_t **_gnum_points   = malloc (sizeof(PDM_g_num_t *) * pcloud->n_part);
      double      **_coords_points = malloc (sizeof(double      *) * pcloud->n_part);

      /* Part to block on elements from location array (be carreful : partial part to block) */
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        elt_weight[ipart]     = malloc (sizeof(double     ) *     pcloud->n_located[ipart]);
        ptb_stride[ipart]     = malloc (sizeof(int        ) *     pcloud->n_located[ipart]);
        weight_stride[ipart]  = malloc (sizeof(int        ) *     pcloud->n_located[ipart]);
        _gnum_points[ipart]   = malloc (sizeof(PDM_g_num_t) *     pcloud->n_located[ipart]);
        _coords_points[ipart] = malloc (sizeof(double     ) * 3 * pcloud->n_located[ipart]);

        for (int j = 0; j < pcloud->n_located[ipart]; j++) {
          elt_weight[ipart][j]    = 1.;
          ptb_stride[ipart][j]   = 1;
          weight_stride[ipart][j] = pcloud->weights_idx[ipart][j+1] -
                                    pcloud->weights_idx[ipart][j];
          _gnum_points[ipart][j] = pcloud->gnum[ipart][pcloud->located[ipart][j]-1];
          for (int k = 0; k < 3; k++) {
            _coords_points[ipart][3*j+k] = pcloud->coords[ipart][3*(pcloud->located[ipart][j]-1)+k];
          }
        }
      }

      PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                          pcloud->location,
                                                          elt_weight,
                                                          pcloud->n_located,
                                                          pcloud->n_part,
                                                          ml->comm);

      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        free (elt_weight[ipart]);
      }
      free (elt_weight);

      int block_n_elt = PDM_part_to_block_n_elt_block_get (ptb);
      block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);

      PDM_g_num_t *block_g_num_location = PDM_part_to_block_block_gnum_get (ptb);

      /*
       * Exchange data
       */

      /* _gnum_points */

      int *block_pts_per_elt_n = NULL;
      PDM_g_num_t *block_pts_g_num = NULL;
      PDM_part_to_block_exch (ptb,
                              sizeof(PDM_g_num_t),
                              PDM_STRIDE_VAR_INTERLACED,
                              1,
                              ptb_stride,
                              (void **) _gnum_points,
                              &block_pts_per_elt_n,
                              (void **) &block_pts_g_num);

      // int *block_pts_per_elt_idx = malloc(sizeof(int) * (block_n_elt + 1));
      // block_pts_per_elt_idx[0] = 0;
      // for (int i = 0; i < block_n_elt; i++) {
      //   block_pts_per_elt_idx[i + 1] = block_pts_per_elt_idx[i] + block_pts_per_elt_n[i];
      // }

      // printf("block_pts_g_num : \n");
      // for (int i = 0; i< block_n_elt; i++) {
      //   printf("i %d %d : ", my_rank, block_pts_per_elt_n[i]);
      //   for (int j = block_pts_per_elt_idx[i]; j< block_pts_per_elt_idx[i+1]; j++) {
      //     printf(PDM_FMT_G_NUM, block_pts_g_num[j]);
      //   }
      //   printf("\n");
      // }
      // free(block_pts_per_elt_idx);

      // printf("_block_distrib_idx :");
      // for (int i = 0; i< n_procs + 1; i++) {
      //   printf(PDM_FMT_G_NUM,  block_distrib_idx[i]);
      // }
      // printf("\n");

      free (block_pts_per_elt_n);
      block_pts_per_elt_n = NULL;

      /* _weight_stride */

      int *block_pts_weights_stride = NULL;
      PDM_part_to_block_exch (ptb,
                              sizeof(int),
                              PDM_STRIDE_VAR_INTERLACED,
                              1,
                              ptb_stride,
                              (void **) weight_stride,
                              &block_pts_per_elt_n,
                              (void **) &block_pts_weights_stride);

      free (block_pts_per_elt_n);
      block_pts_per_elt_n = NULL;

      /* dist2 */

      double *block_pts_dist2 = NULL;
      PDM_part_to_block_exch (ptb,
                              sizeof(double),
                              PDM_STRIDE_VAR_INTERLACED,
                              1,
                              ptb_stride,
                              (void **) pcloud->dist2,
                              &block_pts_per_elt_n,
                              (void **) &block_pts_dist2);

      //free (block_pts_per_elt_n);
      //block_pts_per_elt_n = NULL;

      /* _coords_points */

      int *block_pts_per_elt_n2;
      double *block_pts_coords = NULL;
      PDM_part_to_block_exch (ptb,
                              3*sizeof(double),
                              PDM_STRIDE_VAR_INTERLACED,
                              1,
                              ptb_stride,
                              (void **) _coords_points,
                              &block_pts_per_elt_n2,
                              (void **) &block_pts_coords);

      free (block_pts_per_elt_n2);
      block_pts_per_elt_n2 = NULL;

      /* _proj */

      double *block_pts_proj = NULL;
      PDM_part_to_block_exch (ptb,
                              3*sizeof(double),
                              PDM_STRIDE_VAR_INTERLACED,
                              1,
                              ptb_stride,
                              (void **) pcloud->projected_coords,
                              &block_pts_per_elt_n2,
                              (void **) &block_pts_proj);

      free (block_pts_per_elt_n2);
      block_pts_per_elt_n2 = NULL;

      /* _weight */

      double *block_pts_weights = NULL;
      PDM_part_to_block_exch (ptb,
                              sizeof(double),
                              PDM_STRIDE_VAR_INTERLACED,
                              1,
                              weight_stride,
                              (void **) pcloud->weights,
                              &block_pts_per_elt_n2,
                              (void **) &block_pts_weights);

      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        free (ptb_stride[ipart]);
        free (weight_stride[ipart]);
        free (_gnum_points[ipart]);
        free (_coords_points[ipart]);
      }

      free (ptb_stride);
      free (weight_stride);
      free (_gnum_points);
      free (_coords_points);

      PDM_MPI_Barrier (ml->comm);
      PDM_timer_hang_on(ml->timer);
      e_t_elapsed = PDM_timer_elapsed(ml->timer);
      e_t_cpu     = PDM_timer_cpu(ml->timer);
      e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
      e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

      ml->times_elapsed[REVERSE_LOCATION_DATA_PTB] += e_t_elapsed - b_t_elapsed;
      ml->times_cpu[REVERSE_LOCATION_DATA_PTB]     += e_t_cpu - b_t_cpu;
      ml->times_cpu_u[REVERSE_LOCATION_DATA_PTB]   += e_t_cpu_u - b_t_cpu_u;
      ml->times_cpu_s[REVERSE_LOCATION_DATA_PTB]   += e_t_cpu_s - b_t_cpu_s;

      b_t_elapsed = e_t_elapsed;
      b_t_cpu     = e_t_cpu;
      b_t_cpu_u   = e_t_cpu_u;
      b_t_cpu_s   = e_t_cpu_s;
      PDM_timer_resume(ml->timer);

      /* Block to part */
      int n_part_nodal = 0;
      if (ml->mesh_nodal != NULL) {
        n_part_nodal =PDM_Mesh_nodal_n_part_get (ml->mesh_nodal);
      }
      PDM_g_num_t **numabs_nodal = malloc (sizeof(PDM_g_num_t *) * n_part_nodal);
      int *n_elt_nodal = malloc (sizeof(int) * n_part_nodal);

      for (int i_part = 0; i_part < n_part_nodal; i_part++) {
        n_elt_nodal[i_part]= PDM_Mesh_nodal_n_cell_get (ml->mesh_nodal,
                                                        i_part);
        numabs_nodal[i_part] =  PDM_Mesh_nodal_g_num_get_from_part (ml->mesh_nodal,
                                                                    i_part);
      }

      _points_in_element_t *_points_in_elements = ml->points_in_elements + icloud;

      btp = PDM_block_to_part_create_from_sparse_block (block_g_num_location,
                                                        block_n_elt,
                                 (const PDM_g_num_t **) numabs_nodal,
                                                        n_elt_nodal,
                                                        n_part_nodal,
                                                        ml->comm);

      PDM_MPI_Barrier (ml->comm);
      PDM_timer_hang_on(ml->timer);
      e_t_elapsed = PDM_timer_elapsed(ml->timer);
      e_t_cpu     = PDM_timer_cpu(ml->timer);
      e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
      e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

      ml->times_elapsed[REVERSE_LOCATION_DATA_BTP1] += e_t_elapsed - b_t_elapsed;
      ml->times_cpu[REVERSE_LOCATION_DATA_BTP1]     += e_t_cpu - b_t_cpu;
      ml->times_cpu_u[REVERSE_LOCATION_DATA_BTP1]   += e_t_cpu_u - b_t_cpu_u;
      ml->times_cpu_s[REVERSE_LOCATION_DATA_BTP1]   += e_t_cpu_s - b_t_cpu_s;

      b_t_elapsed = e_t_elapsed;
      b_t_cpu     = e_t_cpu;
      b_t_cpu_u   = e_t_cpu_u;
      b_t_cpu_s   = e_t_cpu_s;
      PDM_timer_resume(ml->timer);

      PDM_part_to_block_free (ptb);

      // PDM_log_trace_array_long(numabs_nodal[0], n_elt_nodal[0], "numabs_nodal :: ");
      /* _gnum_points  _points_in_elements->gnum */

      int **pts_in_elt_n;
      PDM_block_to_part_exch (btp,
                               sizeof(PDM_g_num_t),
                               PDM_STRIDE_VAR_INTERLACED,
                               block_pts_per_elt_n,
                               (void *) block_pts_g_num,
                               &pts_in_elt_n,
                               (void ***) &(_points_in_elements->gnum));

      _points_in_elements->n_part = n_part_nodal; // To see with Eric
      _points_in_elements->n_elts = malloc (sizeof(int) * n_part_nodal);
      _points_in_elements->pts_inside_idx = malloc(sizeof(int *) * n_part_nodal);
      for (int i_part = 0; i_part < n_part_nodal; i_part++) {
        _points_in_elements->pts_inside_idx[i_part] = malloc(sizeof(int) * (n_elt_nodal[i_part]+1));
        _points_in_elements->pts_inside_idx[i_part][0] = 0;
        for (int i = 0; i < n_elt_nodal[i_part]; i++) {
          _points_in_elements->pts_inside_idx[i_part][i+1] = _points_in_elements->pts_inside_idx[i_part][i] +
                                                             pts_in_elt_n[i_part][i];
        }

        _points_in_elements->n_elts[i_part] = n_elt_nodal[i_part]; // Not sure See with Eric
        free (pts_in_elt_n[i_part]);
      }
      free (pts_in_elt_n);
      free (block_pts_g_num);

      /* _weight_stride block_pts_weights_stride*/

      int **pts_weights_stride;
      PDM_block_to_part_exch (btp,
                               sizeof(int),
                               PDM_STRIDE_VAR_INTERLACED,
                               block_pts_per_elt_n,
                               (void *) block_pts_weights_stride,
                               &pts_in_elt_n,
                               (void ***) &(pts_weights_stride));

      _points_in_elements->weights_idx = malloc(sizeof(int *) * n_part_nodal);
      for (int i_part = 0; i_part < n_part_nodal; i_part++) {
        _points_in_elements->weights_idx[i_part] =
             malloc(sizeof(int) * (_points_in_elements->pts_inside_idx[i_part][n_elt_nodal[i_part]]+1));
        _points_in_elements->weights_idx[i_part][0] = 0;
        for (int i = 0; i < _points_in_elements->pts_inside_idx[i_part][n_elt_nodal[i_part]]; i++) {
          _points_in_elements->weights_idx[i_part][i+1] = _points_in_elements->weights_idx[i_part][i] +
                                                          pts_weights_stride[i_part][i];
        }
        free (pts_in_elt_n[i_part]);
        free (pts_weights_stride[i_part]);
      }
      free (pts_in_elt_n);
      free (pts_weights_stride);

      /* dist2 block_pts_dist2 */

      PDM_block_to_part_exch (btp,
                               sizeof(double),
                               PDM_STRIDE_VAR_INTERLACED,
                               block_pts_per_elt_n,
                               (void *) block_pts_dist2,
                               &pts_in_elt_n,
                               (void ***) &(_points_in_elements->dist2));

      for (int i_part = 0; i_part < n_part_nodal; i_part++) {
        free (pts_in_elt_n[i_part]);
      }
      free (pts_in_elt_n);
      free (block_pts_dist2);

      /* _coords_points */

      /*for (int i = 0; i < block_n_elt; i++) {
        block_pts_per_elt_n[i] *= 3;
        }*/

      // PDM_log_trace_array_int(block_pts_per_elt_n , block_n_elt, "block_pts_per_elt_n  AFTER :: ");
      // PDM_log_trace_array_int(block_pts_per_elt_n2, block_n_elt, "block_pts_per_elt_n2 AFTER :: ");
      PDM_block_to_part_exch (btp,
                               3*sizeof(double),
                               PDM_STRIDE_VAR_INTERLACED,
                               block_pts_per_elt_n,
                               (void *) block_pts_coords,
                               &pts_in_elt_n,
                               (void ***) &(_points_in_elements->coords));

      for (int i_part = 0; i_part < n_part_nodal; i_part++) {
        free (pts_in_elt_n[i_part]);
      }
      free (pts_in_elt_n);
      free (block_pts_coords);

      /* _proj */

      PDM_block_to_part_exch (btp,
                               3*sizeof(double),
                               PDM_STRIDE_VAR_INTERLACED,
                               block_pts_per_elt_n,
                               (void *) block_pts_proj,
                               &pts_in_elt_n,
                               (void ***) &(_points_in_elements->projected_coords));

      for (int i_part = 0; i_part < n_part_nodal; i_part++) {
        free (pts_in_elt_n[i_part]);
      }
      free (pts_in_elt_n);
      free (block_pts_proj);

      /* _weight */

      int k5 = 0;
      int *block_pts_weights_stride_sum = malloc (sizeof(int) * block_n_elt);
      for (int k1 = 0; k1 < block_n_elt; k1++) {
        block_pts_weights_stride_sum[k1] = 0;
        for (int k4 = 0; k4 < block_pts_per_elt_n[k1]; k4++) {
          block_pts_weights_stride_sum[k1] += block_pts_weights_stride[k5++];
        }
      }

      PDM_block_to_part_exch (btp,
                               sizeof(double),
                               PDM_STRIDE_VAR_INTERLACED,
                               block_pts_weights_stride_sum,
                               (void *) block_pts_weights,
                               &pts_in_elt_n,
                               (void ***) &(_points_in_elements->weights));
      free (block_pts_weights_stride);
      free (block_pts_weights_stride_sum);

      int elt_n_pts_max    = 0;
      int elt_n_weight_max = 0;

      for (int i_part = 0; i_part < n_part_nodal; i_part++) {
        for (int i = 0; i < n_elt_nodal[i_part]; i++) {
          int elt_i_pt = _points_in_elements->pts_inside_idx[i_part][i];

          int elt_n_pts = _points_in_elements->pts_inside_idx[i_part][i+1] - elt_i_pt;

          elt_n_pts_max = PDM_MAX(elt_n_pts_max, elt_n_pts);

          for (int j = _points_in_elements->pts_inside_idx[i_part][i]; 
                   j < _points_in_elements->pts_inside_idx[i_part][i+1]; j++) {

            elt_n_weight_max = PDM_MAX (elt_n_weight_max, 
                                        _points_in_elements->weights_idx[i_part][j+1] -
                                        _points_in_elements->weights_idx[i_part][j]); 
          }    
        }
      }

      double *tmp_elt_coords           = malloc (sizeof(double     ) * 3 * elt_n_pts_max);
      double *tmp_elt_projected_coords = malloc (sizeof(double     ) * 3 * elt_n_pts_max);
      double *tmp_elt_dist2            = malloc (sizeof(double     ) *     elt_n_pts_max); 
      int    *tmp_elt_weights_idx      = malloc (sizeof(int        ) *     elt_n_pts_max);
      double *tmp_elt_weights          = malloc (sizeof(double     ) *     elt_n_weight_max * elt_n_pts_max);

      int    *order                    = malloc (sizeof(int        ) *     elt_n_pts_max);

      if (1 == 0) {

        for (int i_part = 0; i_part < n_part_nodal; i_part++) {

          for (int i = 0; i < n_elt_nodal[i_part]; i++) {
            printf("elt 0: %d\n", i); 

            int elt_i_pt = _points_in_elements->pts_inside_idx[i_part][i];

            int elt_n_pts = _points_in_elements->pts_inside_idx[i_part][i+1] - elt_i_pt;

            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              printf (" "PDM_FMT_G_NUM"", _points_in_elements->gnum[i_part][j]);
            }
            printf("\n");
      
            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              printf (" %12.5e  %12.5e  %12.5e /", _points_in_elements->coords[i_part][3*j], _points_in_elements->coords[i_part][3*j+1], _points_in_elements->coords[i_part][3*j+2]);
            }
            printf("\n");

            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              printf (" %12.5e  %12.5e  %12.5e /", _points_in_elements->projected_coords[i_part][3*j], _points_in_elements->projected_coords[i_part][3*j+1], _points_in_elements->projected_coords[i_part][3*j+2]);
            }
            printf("\n");

            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              printf (" %12.5e", _points_in_elements->dist2[i_part][j]);
            }
            printf("\n");

            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              printf (" %d", _points_in_elements->weights_idx[i_part][j]);
            }
            printf("\n");
      
            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              for (int k = _points_in_elements->weights_idx[i_part][j]; k < _points_in_elements->weights_idx[i_part][j+1]; k++) {
                printf (" %12.5e", _points_in_elements->weights[i_part][k]);
              }
              printf("/");
            }
            printf("\n");

          }
        }
      }

      PDM_MPI_Barrier (ml->comm);
      PDM_timer_hang_on(ml->timer);
      e_t_elapsed = PDM_timer_elapsed(ml->timer);
      e_t_cpu     = PDM_timer_cpu(ml->timer);
      e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
      e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

      ml->times_elapsed[REVERSE_LOCATION_DATA_BTP2] += e_t_elapsed - b_t_elapsed;
      ml->times_cpu[REVERSE_LOCATION_DATA_BTP2]     += e_t_cpu - b_t_cpu;
      ml->times_cpu_u[REVERSE_LOCATION_DATA_BTP2]   += e_t_cpu_u - b_t_cpu_u;
      ml->times_cpu_s[REVERSE_LOCATION_DATA_BTP2]   += e_t_cpu_s - b_t_cpu_s;

      b_t_elapsed = e_t_elapsed;
      b_t_cpu     = e_t_cpu;
      b_t_cpu_u   = e_t_cpu_u;
      b_t_cpu_s   = e_t_cpu_s;
      PDM_timer_resume(ml->timer);

      for (int i_part = 0; i_part < n_part_nodal; i_part++) {

        int idx1       = 0;
        int idx_weight = 0;

        for (int i = 0; i < n_elt_nodal[i_part]; i++) {

          int elt_i_pt = _points_in_elements->pts_inside_idx[i_part][i];

          int elt_n_pts = _points_in_elements->pts_inside_idx[i_part][i+1] - elt_i_pt;

          _points_in_elements->pts_inside_idx[i_part][i] = idx1;

          for (int j = 0; j < elt_n_pts; j++) {  
            order[j] = j;
          }

          int n_vtx_elt = -1;
          if (elt_n_pts > 0) {
            n_vtx_elt = _points_in_elements->weights_idx[i_part][elt_i_pt + 1] - _points_in_elements->weights_idx[i_part][elt_i_pt ];
          }

          PDM_g_num_t *_elt_gnum             = _points_in_elements->gnum[i_part]             +     elt_i_pt;
          double      *_elt_coords           = _points_in_elements->coords[i_part]           + 3 * elt_i_pt;
          double      *_elt_projected_coords = _points_in_elements->projected_coords[i_part] + 3 * elt_i_pt;
          double      *_elt_dist2            = _points_in_elements->dist2[i_part]            +     elt_i_pt; 
          int         *_elt_weights_idx      = _points_in_elements->weights_idx[i_part]      +     elt_i_pt;
          double      *_elt_weights          = _points_in_elements->weights[i_part]          +     _elt_weights_idx[0];

          memcpy (tmp_elt_coords,           _elt_coords,           sizeof(double) * 3 *         elt_n_pts);
          memcpy (tmp_elt_projected_coords, _elt_projected_coords, sizeof(double) * 3 *         elt_n_pts);
          memcpy (tmp_elt_dist2,            _elt_dist2,            sizeof(double) *             elt_n_pts);
          memcpy (tmp_elt_weights_idx,      _elt_weights_idx,      sizeof(int)    *             elt_n_pts);
          memcpy (tmp_elt_weights,          _elt_weights,          sizeof(double) * n_vtx_elt * elt_n_pts);

          PDM_sort_long (_elt_gnum, order, elt_n_pts);

          PDM_g_num_t pre = -1;

          for (int j = 0; j < elt_n_pts; j++) {
            if (pre != _elt_gnum[j]) {
              _points_in_elements->gnum[i_part][idx1]          = _elt_gnum[j];
              _points_in_elements->dist2[i_part][idx1]         = tmp_elt_dist2[order[j]];
              _points_in_elements->weights_idx[i_part][idx1+1] = _points_in_elements->weights_idx[i_part][idx1] + n_vtx_elt; 
              for (int k = 0; k < 3; k++) {
                _points_in_elements->coords[i_part][3*idx1+k]           = tmp_elt_coords[3*order[j]+k];
                _points_in_elements->projected_coords[i_part][3*idx1+k] = tmp_elt_projected_coords[3*order[j]+k];
              }

              for (int k = 0; k < n_vtx_elt; k++) {
                _points_in_elements->weights[i_part][idx_weight++] = tmp_elt_weights[n_vtx_elt*order[j]+k]; 
              }

              pre = _elt_gnum[j];
              idx1++;
            }
          }      
        }

        _points_in_elements->pts_inside_idx[i_part][n_elt_nodal[i_part]] = idx1;

        _points_in_elements->gnum[i_part]             = realloc (_points_in_elements->gnum[i_part], 
                                                                  sizeof(PDM_g_num_t) * _points_in_elements->pts_inside_idx[i_part][n_elt_nodal[i_part]]);

        _points_in_elements->dist2[i_part]            = realloc (_points_in_elements->dist2[i_part], 
                                                                  sizeof(double) * _points_in_elements->pts_inside_idx[i_part][n_elt_nodal[i_part]]);

        _points_in_elements->weights_idx[i_part]      = realloc (_points_in_elements->weights_idx[i_part], 
                                                                  sizeof(int) * (_points_in_elements->pts_inside_idx[i_part][n_elt_nodal[i_part]] + 1)); 

        _points_in_elements->coords[i_part]           = realloc (_points_in_elements->coords[i_part], 
                                                                  sizeof(double) * 3 * _points_in_elements->pts_inside_idx[i_part][n_elt_nodal[i_part]]);

        _points_in_elements->projected_coords[i_part] = realloc (_points_in_elements->projected_coords[i_part], 
                                                                  sizeof(double) * 3 * _points_in_elements->pts_inside_idx[i_part][n_elt_nodal[i_part]]);

        _points_in_elements->weights[i_part]          = realloc (_points_in_elements->weights[i_part], 
                                                                  sizeof(double) * (_points_in_elements->weights_idx[i_part][_points_in_elements->pts_inside_idx[i_part][n_elt_nodal[i_part]]]+1)); 

      }

      for (int i_part = 0; i_part < n_part_nodal; i_part++) {
        free (pts_in_elt_n[i_part]);
      }

      if (1 == 0) {
        for (int i_part = 0; i_part < n_part_nodal; i_part++) {

          for (int i = 0; i < n_elt_nodal[i_part]; i++) {
            printf("elt 2 : %d\n", i); 

            int elt_i_pt = _points_in_elements->pts_inside_idx[i_part][i];

            int elt_n_pts = _points_in_elements->pts_inside_idx[i_part][i+1] - elt_i_pt;

            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              printf (" "PDM_FMT_G_NUM"", _points_in_elements->gnum[i_part][j]);
            }
            printf("\n");
      
            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              printf (" %12.5e  %12.5e  %12.5e /", _points_in_elements->coords[i_part][3*j], _points_in_elements->coords[i_part][3*j+1], _points_in_elements->coords[i_part][3*j+2]);
            }
            printf("\n");

            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              printf (" %12.5e  %12.5e  %12.5e /", _points_in_elements->projected_coords[i_part][3*j], _points_in_elements->projected_coords[i_part][3*j+1], _points_in_elements->projected_coords[i_part][3*j+2]);
            }
            printf("\n");

            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              printf (" %12.5e", _points_in_elements->dist2[i_part][j]);
            }
            printf("\n");

            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              printf (" %d", _points_in_elements->weights_idx[i_part][j]);
            }
            printf("\n");
      
            for (int j = elt_i_pt; j < elt_i_pt + elt_n_pts; j++) {
              for (int k = _points_in_elements->weights_idx[i_part][j]; k < _points_in_elements->weights_idx[i_part][j+1]; k++) {
                printf (" %12.5e", _points_in_elements->weights[i_part][k]);
              }
              printf("/");
            }
            printf("\n");

          }
        }
      }

      PDM_MPI_Barrier (ml->comm);
      PDM_timer_hang_on(ml->timer);
      e_t_elapsed = PDM_timer_elapsed(ml->timer);
      e_t_cpu     = PDM_timer_cpu(ml->timer);
      e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
      e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

      ml->times_elapsed[REVERSE_LOCATION_DATA_BTP3] += e_t_elapsed - b_t_elapsed;
      ml->times_cpu[REVERSE_LOCATION_DATA_BTP3]     += e_t_cpu - b_t_cpu;
      ml->times_cpu_u[REVERSE_LOCATION_DATA_BTP3]   += e_t_cpu_u - b_t_cpu_u;
      ml->times_cpu_s[REVERSE_LOCATION_DATA_BTP3]   += e_t_cpu_s - b_t_cpu_s;

      b_t_elapsed = e_t_elapsed;
      b_t_cpu     = e_t_cpu;
      b_t_cpu_u   = e_t_cpu_u;
      b_t_cpu_s   = e_t_cpu_s;
      PDM_timer_resume(ml->timer);

      /* Compute uvw */
    
      _points_in_elements->uvw = malloc (sizeof(double *) * n_part_nodal);

      for (int i_part = 0; i_part < n_part_nodal; i_part++) {
        _points_in_elements->uvw[i_part] = malloc(sizeof(double) * 3 * _points_in_elements->pts_inside_idx[i_part][n_elt_nodal[i_part]]);
      }

      const double newton_tol = 1.e-6;

      // Allocation uvw
      double *_cell_coords = malloc (sizeof(double) * 8 * 3);
      double *_cell_coords_ijk = malloc (sizeof(double) * 8 * 3);
      for (int ipart = 0; ipart < n_part_nodal; ipart++) {
        int ielt = -1;
        const double *coords_vtx = PDM_Mesh_nodal_vertices_get (ml->mesh_nodal,
                                                                ipart);
        for (int iblock = 0; iblock < n_blocks; iblock++) {
          int id_block = blocks_id[iblock];

          PDM_Mesh_nodal_elt_t t_elt = PDM_Mesh_nodal_block_type_get (ml->mesh_nodal,
                                                                      id_block);

          if (t_elt != PDM_MESH_NODAL_PYRAMID5 &&
              t_elt != PDM_MESH_NODAL_PRISM6   &&
              t_elt != PDM_MESH_NODAL_HEXA8) {
            continue;
          }
          int *parent_num = PDM_Mesh_nodal_block_parent_num_get (ml->mesh_nodal,
                                                                 id_block,
                                                                 ipart);

          int n_elt_block = PDM_Mesh_nodal_block_n_elt_get (ml->mesh_nodal,
                                                            id_block,
                                                            ipart);

          int *connec;
          PDM_Mesh_nodal_block_std_get (ml->mesh_nodal,
                                        id_block,
                                        ipart,
                                        &connec);

          int n_vtx = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, 1);

          // PDM_log_trace_array_double(_points_in_elements->projected_coords[ipart],  3 * _points_in_elements->pts_inside_idx[ipart][n_elt_block], "projected_coords");
          // PDM_log_trace_array_double(_points_in_elements->coords[ipart], 3 * _points_in_elements->pts_inside_idx[ipart][n_elt_block], "coords");

          for (int i = 0; i < n_elt_block; i++) {
            int order_ijk = 0;
            for (int k = 0; k < n_vtx; k++) {
              int ivtx = connec[i*n_vtx + k] - 1;
              for (int j = 0; j < 3; j++) {
                _cell_coords[3*k+j] = coords_vtx[3*ivtx+j];
              }
            }

            double _projected_coords[3]; // Ignored in this case

            int ielt_parent = ++ielt;
            if (parent_num != NULL) {
              ielt_parent = parent_num[i];
            }

            //int idx_pts_elts = _points_in_elements->pts_inside_idx[ipart][ielt_parent];
            //int n_pts_elts = _points_in_elements->pts_inside_idx[ipart][ielt_parent+1] - idx_pts_elts;
            //for (int k1 = 0; k1 < n_pts_elts; k1++) {
            for (int ipt = _points_in_elements->pts_inside_idx[ipart][ielt_parent];
                 ipt < _points_in_elements->pts_inside_idx[ipart][ielt_parent+1];
                 ipt++) {

              const double *_point_coords;
              double *_point_uvw = _points_in_elements->uvw[ipart] + 3*ipt;

              //if (_points_in_elements->dist2[ipart][idx_pts_elts + k1] >= 0.) {
              if (_points_in_elements->dist2[ipart][ipt] >= 0.) {
                _point_coords = _points_in_elements->projected_coords[ipart] + 3*ipt;
                //&(_points_in_elements->projected_coords[ipart][3*(idx_pts_elts + k1)]);
              }
              else {
                _point_coords = _points_in_elements->coords[ipart] + 3*ipt;
                //&(_points_in_elements->coords[ipart][3*(idx_pts_elts + k1)]);
              }

              PDM_bool_t stat = PDM_point_location_compute_uvw (t_elt,
                                                                _point_coords,
                                                                _cell_coords,
                                                                newton_tol,
                                                                _point_uvw);
              //&(_points_in_elements->uvw[ipart][3 * (idx_pts_elts + k1)]));

              int uvw_ok = 0;
              if (stat) {
                uvw_ok = (_point_uvw[0] > -newton_tol && _point_uvw[0] < 1+newton_tol &&
                          _point_uvw[1] > -newton_tol && _point_uvw[1] < 1+newton_tol &&
                          _point_uvw[2] > -newton_tol && _point_uvw[2] < 1+newton_tol);
              }

              if (!uvw_ok) {
                /* Newton failed, try subdivision method */
                if (!order_ijk) {
                  /* Get cell coord in ijk order */
                  order_ijk = 1;
                  if (t_elt == PDM_MESH_NODAL_PRISM6 ||
                      t_elt == PDM_MESH_NODAL_TETRA4 ||
                      t_elt == PDM_MESH_NODAL_TRIA3  ||
                      t_elt == PDM_MESH_NODAL_BAR2) {
                    memcpy(_cell_coords_ijk, _cell_coords, sizeof(double) * n_vtx * 3);
                  } else {
                    _cell_coords_ijk[ 0] = _cell_coords[ 0];
                    _cell_coords_ijk[ 1] = _cell_coords[ 1];
                    _cell_coords_ijk[ 2] = _cell_coords[ 2];

                    _cell_coords_ijk[ 3] = _cell_coords[ 3];
                    _cell_coords_ijk[ 4] = _cell_coords[ 4];
                    _cell_coords_ijk[ 5] = _cell_coords[ 5];

                    _cell_coords_ijk[ 6] = _cell_coords[ 9];
                    _cell_coords_ijk[ 7] = _cell_coords[10];
                    _cell_coords_ijk[ 8] = _cell_coords[11];

                    _cell_coords_ijk[ 9] = _cell_coords[ 6];
                    _cell_coords_ijk[10] = _cell_coords[ 7];
                    _cell_coords_ijk[11] = _cell_coords[ 8];

                    if (t_elt == PDM_MESH_NODAL_PYRAMID5 ||
                        t_elt == PDM_MESH_NODAL_HEXA8 ) {
                      _cell_coords_ijk[12] = _cell_coords[12];
                      _cell_coords_ijk[13] = _cell_coords[13];
                      _cell_coords_ijk[14] = _cell_coords[14];

                      if (t_elt == PDM_MESH_NODAL_HEXA8) {
                        _cell_coords_ijk[15] = _cell_coords[15];
                        _cell_coords_ijk[16] = _cell_coords[16];
                        _cell_coords_ijk[17] = _cell_coords[17];

                        _cell_coords_ijk[18] = _cell_coords[21];
                        _cell_coords_ijk[19] = _cell_coords[22];
                        _cell_coords_ijk[20] = _cell_coords[23];

                        _cell_coords_ijk[21] = _cell_coords[18];
                        _cell_coords_ijk[22] = _cell_coords[19];
                        _cell_coords_ijk[23] = _cell_coords[20];
                      }
                    }
                  }
                }

                PDM_Mesh_nodal_elt_t t_elt_ho;
                switch (t_elt) {
                  case PDM_MESH_NODAL_BAR2:
                    t_elt_ho = PDM_MESH_NODAL_BARHO;
                    break;

                  case PDM_MESH_NODAL_TRIA3:
                    t_elt_ho = PDM_MESH_NODAL_TRIAHO;
                    break;

                  case PDM_MESH_NODAL_QUAD4:
                    t_elt_ho = PDM_MESH_NODAL_QUADHO;
                    break;

                  case PDM_MESH_NODAL_TETRA4:
                    t_elt_ho = PDM_MESH_NODAL_TETRAHO;
                    break;

                  case PDM_MESH_NODAL_PYRAMID5:
                    t_elt_ho = PDM_MESH_NODAL_PYRAMIDHO;
                    break;

                  case PDM_MESH_NODAL_PRISM6:
                    t_elt_ho = PDM_MESH_NODAL_PRISMHO;
                    break;

                  case PDM_MESH_NODAL_HEXA8:
                    t_elt_ho = PDM_MESH_NODAL_HEXAHO;
                    break;

                  default:
                    t_elt_ho = t_elt;
                }

                PDM_ho_location (t_elt_ho,
                                 1,
                                 n_vtx,
                                 _cell_coords_ijk,
                                 _point_coords,
                                 _projected_coords,
                                 _point_uvw);
                //&(_points_in_elements->uvw[ipart][3 * (idx_pts_elts + k1)]));
              }

              // Check uvw
              if (!(_point_uvw[0] > -newton_tol && _point_uvw[0] < 1+newton_tol &&
                    _point_uvw[1] > -newton_tol && _point_uvw[1] < 1+newton_tol &&
                    _point_uvw[2] > -newton_tol && _point_uvw[2] < 1+newton_tol)) {//_points_in_elements->dist2[ipart][ipt] <= 0.) {
                log_trace("\n\nt_elt = %d\n", (int) t_elt);
                log_trace("_point_coords = %20.16f %20.16f %20.16f\n",
                          _point_coords[0], _point_coords[1], _point_coords[2]);
                log_trace("_cell_coords =\n");
                for (int k = 0; k < n_vtx; k++) {
                  log_trace("%20.16f %20.16f %20.16f\n",
                            _cell_coords[3*k+0], _cell_coords[3*k+1], _cell_coords[3*k+2]);
                }
                log_trace("_point_uvw = %20.16f %20.16f %20.16f\n",
                          _point_uvw[0], _point_uvw[1], _point_uvw[2]);
              }
              assert(_point_uvw[0] > -newton_tol && _point_uvw[0] < 1+newton_tol &&
                     _point_uvw[1] > -newton_tol && _point_uvw[1] < 1+newton_tol &&
                     _point_uvw[2] > -newton_tol && _point_uvw[2] < 1+newton_tol);
            }
          }
        }
      }
      free (_cell_coords_ijk);
      free (_cell_coords);

      free (pts_in_elt_n);
      free (block_pts_weights);

      free (block_pts_per_elt_n);
      free (block_pts_per_elt_n2);
      free (numabs_nodal);
      free (n_elt_nodal);
      PDM_block_to_part_free (btp);
      free (order);
      free (tmp_elt_coords);
      free (tmp_elt_projected_coords);
      free (tmp_elt_dist2);
      free (tmp_elt_weights_idx);
      free (tmp_elt_weights);

      PDM_MPI_Barrier (ml->comm);
      PDM_timer_hang_on(ml->timer);
      e_t_elapsed = PDM_timer_elapsed(ml->timer);
      e_t_cpu     = PDM_timer_cpu(ml->timer);
      e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
      e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

      /*ml->times_elapsed[REVERSE_LOCATION_DATA] += e_t_elapsed - b_t_elapsed;
      ml->times_cpu[REVERSE_LOCATION_DATA]     += e_t_cpu - b_t_cpu;
      ml->times_cpu_u[REVERSE_LOCATION_DATA]   += e_t_cpu_u - b_t_cpu_u;
      ml->times_cpu_s[REVERSE_LOCATION_DATA]   += e_t_cpu_s - b_t_cpu_s;*/
      ml->times_elapsed[REVERSE_LOCATION_DATA_UVW] += e_t_elapsed - b_t_elapsed;
      ml->times_cpu[REVERSE_LOCATION_DATA_UVW]     += e_t_cpu - b_t_cpu;
      ml->times_cpu_u[REVERSE_LOCATION_DATA_UVW]   += e_t_cpu_u - b_t_cpu_u;
      ml->times_cpu_s[REVERSE_LOCATION_DATA_UVW]   += e_t_cpu_s - b_t_cpu_s;
      ml->times_elapsed[REVERSE_LOCATION_DATA] = ml->times_elapsed[REVERSE_LOCATION_DATA_PTB] + 
                                                 ml->times_elapsed[REVERSE_LOCATION_DATA_BTP1] + 
                                                 ml->times_elapsed[REVERSE_LOCATION_DATA_BTP2] + 
                                                 ml->times_elapsed[REVERSE_LOCATION_DATA_BTP3] + 
                                                 ml->times_elapsed[REVERSE_LOCATION_DATA_UVW];
      ml->times_cpu[REVERSE_LOCATION_DATA] = ml->times_cpu[REVERSE_LOCATION_DATA_PTB] +
                                             ml->times_cpu[REVERSE_LOCATION_DATA_BTP1] + 
                                             ml->times_cpu[REVERSE_LOCATION_DATA_BTP2] + 
                                             ml->times_cpu[REVERSE_LOCATION_DATA_BTP3] + 
                                             ml->times_cpu[REVERSE_LOCATION_DATA_UVW];
      ml->times_cpu_u[REVERSE_LOCATION_DATA] = ml->times_cpu_u[REVERSE_LOCATION_DATA_PTB] + 
                                               ml->times_cpu_u[REVERSE_LOCATION_DATA_BTP1] + 
                                               ml->times_cpu_u[REVERSE_LOCATION_DATA_BTP2] + 
                                               ml->times_cpu_u[REVERSE_LOCATION_DATA_BTP3] + 
                                               ml->times_cpu_u[REVERSE_LOCATION_DATA_UVW];
      ml->times_cpu_s[REVERSE_LOCATION_DATA] = ml->times_cpu_s[REVERSE_LOCATION_DATA_PTB] + 
                                               ml->times_cpu_s[REVERSE_LOCATION_DATA_BTP1] + 
                                               ml->times_cpu_s[REVERSE_LOCATION_DATA_BTP2] + 
                                               ml->times_cpu_s[REVERSE_LOCATION_DATA_BTP3] + 
                                               ml->times_cpu_s[REVERSE_LOCATION_DATA_UVW];

      b_t_elapsed = e_t_elapsed;
      b_t_cpu     = e_t_cpu;
      b_t_cpu_u   = e_t_cpu_u;
      b_t_cpu_s   = e_t_cpu_s;
      PDM_timer_resume(ml->timer);
    }

  } // Loop over point clouds


  if (use_extracted_mesh) {
    free(select_box_parent_g_num);
  }

  if(allow_extraction) {
    for (int iblock = 0; iblock < n_blocks; iblock++) {
      for (int ipart = 0; ipart < n_parts; ipart++) {
        free (select_elt_l_num[iblock][ipart]);
      }
      free (n_select_elt[iblock]);
      free (select_elt_l_num[iblock]);
    }
  }
  free (n_select_elt);
  free (select_elt_l_num);
  free(use_extracted_pts);

  if (allow_extraction) {
    free (n_select_pts);
    free (select_pts_parent_g_num);
    free (select_pts_g_num);
    free (select_pts_coord);
  }

  if (select_box_g_num != box_g_num)     free (select_box_g_num);
  if (select_box_extents != box_extents) free (select_box_extents);
  free (box_g_num);
  free (box_extents);


  if (ml->method == PDM_MESH_LOCATION_DBBTREE) {
    PDM_dbbtree_free (dbbt);
    PDM_box_set_destroy (&box_set);
  }

  PDM_timer_hang_on(ml->timer);

  ml->times_elapsed[END] = PDM_timer_elapsed(ml->timer);
  ml->times_cpu[END]     = PDM_timer_cpu(ml->timer);
  ml->times_cpu_u[END]   = PDM_timer_cpu_user(ml->timer);
  ml->times_cpu_s[END]   = PDM_timer_cpu_sys(ml->timer);

  b_t_elapsed = ml->times_elapsed[END];
  b_t_cpu     = ml->times_cpu[END];
  b_t_cpu_u   = ml->times_cpu_u[END];
  b_t_cpu_s   = ml->times_cpu_s[END];
  PDM_timer_resume(ml->timer);
}

/**
 * NOTES:
 * - _points_in_element_t : separate point clouds?
 * - origin frame for elements (one per cloud?)
 *   -> _store_cell_vtx_selected (one per cloud for PDM_mesh_location_cell_vertex_get?)
 */

 static void
_dump_point_cloud
(
 const char     *name,
 PDM_MPI_Comm    comm,
 int             n_part,
 int            *pn_pts,
 double        **ppts_coord,
 PDM_g_num_t   **ppts_ln_to_gn,
 const int       n_field,
 const char    **field_name,
 double       ***field_value
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t **g_num = ppts_ln_to_gn;

  if (1) {
    PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(3,
                                               n_part,
                                               PDM_FALSE,
                                               1.,
                                               comm,
                                               PDM_OWNERSHIP_USER);

    for (int ipart = 0; ipart < n_part; ipart++) {
      PDM_gnum_set_from_parents(gen_gnum,
                                ipart,
                                pn_pts[ipart],
                                ppts_ln_to_gn[ipart]);
    }

    PDM_gnum_compute(gen_gnum);

    g_num = malloc(sizeof(PDM_g_num_t *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      g_num[ipart] = PDM_gnum_get(gen_gnum,
                                  ipart);
    }

    PDM_gnum_free(gen_gnum);
  }

  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        name,
                                        name,
                                        comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       name,
                                       n_part);

  int id_var_num_part = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");

  int id_var[n_field];
  for (int i = 0; i < n_field; i++) {
    id_var[i] = PDM_writer_var_create(wrt,
                                      PDM_WRITER_OFF,
                                      PDM_WRITER_VAR_SCALAR,
                                      PDM_WRITER_VAR_ELEMENTS,
                                      field_name[i]);
  }

  PDM_writer_step_beg(wrt, 0.);

  int id_block = PDM_writer_geom_bloc_add(wrt,
                                          id_geom,
                                          PDM_WRITER_POINT,
                                          PDM_OWNERSHIP_USER);

  PDM_real_t **val_num_part = malloc(sizeof(PDM_real_t  *) * n_part);
  int **connec = malloc(sizeof(int *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    // PDM_log_trace_array_long(ppts_ln_to_gn[ipart],
    //                          pn_pts[ipart],
    //                          "ppts_ln_to_gn[ipart] : ");

    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              ipart,
                              pn_pts[ipart],
                              ppts_coord[ipart],
                              g_num[ipart],
                              PDM_OWNERSHIP_USER);

    connec[ipart] = malloc(sizeof(int) * pn_pts[ipart]);
    for (int i = 0; i < pn_pts[ipart]; i++) {
      connec[ipart][i] = i+1;
    }
    PDM_writer_geom_bloc_std_set(wrt,
                                 id_geom,
                                 id_block,
                                 ipart,
                                 pn_pts[ipart],
                                 connec[ipart],
                                 g_num[ipart]);

    val_num_part[ipart] = malloc(sizeof(PDM_real_t) * pn_pts[ipart]);
    for (int i = 0; i < pn_pts[ipart]; i++) {
      val_num_part[ipart][i] = n_part*i_rank + ipart;
    }

    PDM_writer_var_set(wrt,
                       id_var_num_part,
                       id_geom,
                       ipart,
                       val_num_part[ipart]);

    for (int i = 0; i < n_field; i++) {
      PDM_writer_var_set(wrt,
                         id_var[i],
                         id_geom,
                         ipart,
                         field_value[i][ipart]);
    }
  }

  PDM_writer_geom_write(wrt,
                        id_geom);

  PDM_writer_var_write(wrt,
                       id_var_num_part);
  PDM_writer_var_free(wrt,
                      id_var_num_part);

  for (int i = 0; i < n_field; i++) {
    PDM_writer_var_write(wrt,
                         id_var[i]);
    PDM_writer_var_free(wrt,
                        id_var[i]);
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(val_num_part[ipart]);
    free(connec[ipart]);
    if (g_num != ppts_ln_to_gn) {
      free(g_num[ipart]);
    }
  }
  free(val_num_part);
  free(connec);
  if (g_num != ppts_ln_to_gn) {
    free(g_num);
  }

  PDM_writer_step_end(wrt);

  PDM_writer_free(wrt);
}



static void
_dump_pmne
(
      PDM_MPI_Comm                  comm,
const char                         *prefix,
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           n_part,
      PDM_g_num_t                 **pelt_ln_to_gn,
      int                          *pn_vtx,
      double                      **pvtx_coord
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  char filename[999];

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int ipart = 0; ipart < n_part; ipart++) {

    for (int isection = 0; isection < n_section; isection++) {

      sprintf(filename, "%s_part_%d_section_%d_%3.3d.vtk",
              prefix, ipart, isection, i_rank);

      int id_section = sections_id[isection];
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne,
                                                                             id_section);

      int n_elt = PDM_part_mesh_nodal_elmts_block_n_elt_get(pmne,
                                                             id_section,
                                                             ipart);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 ipart);

      // if (parent_num != NULL) {
      //   log_trace("section %d (type %d), ", isection, t_elt);
      //   PDM_log_trace_array_int(parent_num, n_elt, "parent_num : ");
      // }


      PDM_g_num_t *gnum = malloc(sizeof(PDM_g_num_t) * n_elt);
      for (int i = 0; i < n_elt; i++) {
        int icell = i;
        if (parent_num != NULL) {
          icell = parent_num[i];
        }
        gnum[i] = pelt_ln_to_gn[ipart][icell];
      }


      if (t_elt == PDM_MESH_NODAL_POLY_2D) {

        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_elmts_block_poly2d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec_idx,
                                                   &connec);

        PDM_vtk_write_polydata(filename,
                               pn_vtx[ipart],
                               pvtx_coord[ipart],
                               NULL,
                               n_elt,
                               connec_idx,
                               connec,
                               gnum,
                               NULL);

      }

      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {

        int  n_face;
        int *face_vtx_idx;
        int *face_vtx;
        int *cell_face_idx;
        int *cell_face;
        PDM_g_num_t *face_ln_to_gn    = NULL;
        int         *_parent_num      = NULL;
        PDM_g_num_t *elt_ln_to_gn     = NULL;
        PDM_g_num_t *parent_elt_g_num = NULL;
        PDM_part_mesh_nodal_elmts_block_poly3d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &n_face,
                                                   &face_ln_to_gn,
                                                   &face_vtx_idx,
                                                   &face_vtx,
                                                   &elt_ln_to_gn,
                                                   &cell_face_idx,
                                                   &cell_face,
                                                   &_parent_num,
                                                   &parent_elt_g_num);

        // PDM_log_trace_connectivity_int(cell_face_idx,
        //                                cell_face,
        //                                n_elt,
        //                                "cell_face : ");

        PDM_vtk_write_polydata(filename,
                               pn_vtx[ipart],
                               pvtx_coord[ipart],
                               NULL,
                               n_face,
                               face_vtx_idx,
                               face_vtx,
                               face_ln_to_gn,
                               NULL);

      }

      else {
        int         *elmt_vtx                 = NULL;
        int         *_parent_num              = NULL;
        PDM_g_num_t *numabs                   = NULL;
        PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
        PDM_part_mesh_nodal_elmts_block_std_get(pmne,
                                                id_section,
                                                ipart,
                                                &elmt_vtx,
                                                &numabs,
                                                &_parent_num,
                                                &parent_entitity_ln_to_gn);

        const char *field_name[] = {"parent_num"};
        PDM_vtk_write_std_elements(filename,
                                   pn_vtx[ipart],
                                   pvtx_coord[ipart],
                                   NULL,
                                   t_elt,
                                   n_elt,
                                   elmt_vtx,
                                   gnum,
                                   0,
                                   field_name,
                    (const int **) &parent_num);
      }

      free(gnum);
    }
  }
}



/**
 *
 * \brief Compute point location optim
 *
 * \param [in]   id  Pointer to \ref PDM_mesh_location object
 *
 */

void
PDM_mesh_location_compute_optim
(
 PDM_mesh_location_t        *ml
)
{
  /*
   *  Parameters
   */
  const double tolerance = 1e-6;
  float extraction_threshold = 0.5; // max size ratio between extracted and original meshes

  const int full_async  = 0;
  const int dbg_enabled = 0;
  const int log_timer   = 0;
  const int dim = 3;

  /* Octree parameters */
  const int octree_depth_max             = 31;
  const int octree_points_in_leaf_max    = 1;
  const int octree_build_leaf_neighbours = 0;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (ml->comm, &i_rank);
  PDM_MPI_Comm_size (ml->comm, &n_rank);


  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  // PDM_MPI_Barrier (ml->comm);
  ml->times_elapsed[BEGIN] = PDM_timer_elapsed (ml->timer);
  ml->times_cpu    [BEGIN] = PDM_timer_cpu     (ml->timer);
  ml->times_cpu_u  [BEGIN] = PDM_timer_cpu_user(ml->timer);
  ml->times_cpu_s  [BEGIN] = PDM_timer_cpu_sys (ml->timer);

  b_t_elapsed = ml->times_elapsed[BEGIN];
  b_t_cpu     = ml->times_cpu    [BEGIN];
  b_t_cpu_u   = ml->times_cpu_u  [BEGIN];
  b_t_cpu_s   = ml->times_cpu_s  [BEGIN];
  PDM_timer_resume(ml->timer);


  ml->points_in_elements = malloc(sizeof(_points_in_element_t) * ml->n_point_cloud);

  /*
   *  Compute global extents of source mesh
   *  -------------------------------------
   */

  /* Create dummy part_mesh_nodal */

  PDM_part_mesh_nodal_elmts_t *pmne = _mesh_nodal_to_pmesh_nodal_elmts(ml->mesh_nodal, ml->comm);

  /* Build the bounding boxes of all mesh elements
    (concatenate nodal blocks for each part) */
  int n_block = 0;
  int n_part  = 0;
  int *blocks_id = NULL;

  if (ml->mesh_nodal != NULL) {
    n_block   = PDM_Mesh_nodal_n_blocks_get (ml->mesh_nodal);
    n_part    = PDM_Mesh_nodal_n_part_get   (ml->mesh_nodal);
    blocks_id = PDM_Mesh_nodal_blocks_id_get(ml->mesh_nodal);
  }

  int *req_pts_proj_coord = PDM_array_const_int(ml->n_point_cloud, -1);
  int *req_pts_dist2      = PDM_array_const_int(ml->n_point_cloud, -1);


  int                   *pn_elt      = malloc(sizeof(int                   ) * n_part);
  PDM_g_num_t          **elt_g_num   = malloc(sizeof(PDM_g_num_t          *) * n_part);
  double               **elt_extents = malloc(sizeof(double               *) * n_part);
  PDM_Mesh_nodal_elt_t **elt_type    = malloc(sizeof(PDM_Mesh_nodal_elt_t *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_elt = PDM_Mesh_nodal_n_cell_get(ml->mesh_nodal,
                                          ipart);

    pn_elt     [ipart] = n_elt;
    elt_g_num  [ipart] = malloc(sizeof(PDM_g_num_t         ) * n_elt);
    elt_extents[ipart] = malloc(sizeof(double              ) * n_elt * 6);
    elt_type   [ipart] = malloc(sizeof(PDM_Mesh_nodal_elt_t) * n_elt);
    int idx = 0;
    for (int iblock = 0; iblock < n_block; iblock++) {
      int id_block = blocks_id[iblock];
      int n_elt_in_block = PDM_Mesh_nodal_block_n_elt_get(ml->mesh_nodal,
                                                          id_block,
                                                          ipart);
      double *_extents = malloc(sizeof(double) * n_elt_in_block * 6);
      PDM_Mesh_nodal_compute_cell_extents(ml->mesh_nodal,
                                          id_block,
                                          ipart,
                                          ml->tolerance,
                                          _extents);


      PDM_g_num_t *_elt_g_num = PDM_Mesh_nodal_g_num_get(ml->mesh_nodal,
                                                         id_block,
                                                         ipart);

      PDM_Mesh_nodal_elt_t t_elt = PDM_Mesh_nodal_block_type_get(ml->mesh_nodal,
                                                                 id_block);

      int *parent_num = PDM_Mesh_nodal_block_parent_num_get(ml->mesh_nodal,
                                                            id_block,
                                                            ipart);

      for (int ielt = 0; ielt < n_elt_in_block; ielt++) {
        idx = ielt;
        if (parent_num) {
          idx = parent_num[ielt];
        }
        elt_g_num[ipart][idx] = _elt_g_num[ielt];
        elt_type [ipart][idx] = t_elt;
        memcpy(elt_extents[ipart] + 6*idx, _extents + 6*ielt, sizeof(double)*6);
      }
      free(_extents);
    }
  }

  if (dbg_enabled) {
    int     *pn_vtx     = malloc(sizeof(int     ) * n_part);
    double **pvtx_coord = malloc(sizeof(double *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      pvtx_coord[ipart] = (double *) PDM_Mesh_nodal_vertices_get(ml->mesh_nodal,
                                                                 ipart);

      PDM_g_num_t *pvtx_ln_to_gn = NULL;
      pn_vtx[ipart] = PDM_Mesh_nodal_vertices_ln_to_gn_get(ml->mesh_nodal,
                                                           ipart,
                                                           &pvtx_ln_to_gn);
    }
    _dump_pmne(ml->comm,
               "init_pmne",
               pmne,
               n_part,
               elt_g_num,
               pn_vtx,
               pvtx_coord);
    free(pn_vtx    );
    free(pvtx_coord);
  }

  if (1 && dbg_enabled) {
    PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                          PDM_WRITER_FMT_BIN,
                                          PDM_WRITER_TOPO_CST,
                                          PDM_WRITER_OFF,
                                          "mesh_location__source",
                                          "source_mesh",
                                          ml->comm,
                                          PDM_IO_KIND_MPI_SIMPLE,
                                          1.,
                                          NULL);

    int id_geom = PDM_writer_geom_create_from_mesh_nodal(wrt,
                                                         "source_mesh_geom",
                                                         ml->mesh_nodal);

    int id_var_num_part = PDM_writer_var_create(wrt,
                                                PDM_WRITER_OFF,
                                                PDM_WRITER_VAR_SCALAR,
                                                PDM_WRITER_VAR_ELEMENTS,
                                                "num_part");

    int id_var_gnum = PDM_writer_var_create(wrt,
                                            PDM_WRITER_OFF,
                                            PDM_WRITER_VAR_SCALAR,
                                            PDM_WRITER_VAR_ELEMENTS,
                                            "elt_gnum");

    PDM_writer_step_beg(wrt, 0.);

    PDM_writer_geom_write(wrt,
                          id_geom);

    PDM_real_t **val_num_part = malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_gnum     = malloc (sizeof(PDM_real_t *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      int n_elt = PDM_Mesh_nodal_n_cell_get(ml->mesh_nodal,
                                            ipart);
      val_num_part[ipart] = malloc(sizeof(PDM_real_t) * n_elt);
      val_gnum    [ipart] = malloc(sizeof(PDM_real_t) * n_elt);
      for (int i = 0; i < n_elt; i++) {
        val_num_part[ipart][i] = (PDM_real_t) (i_rank*n_part + ipart);
        val_gnum    [ipart][i] = (PDM_real_t) elt_g_num[ipart][i];
      }

      PDM_writer_var_set(wrt,
                         id_var_num_part,
                         id_geom,
                         ipart,
                         val_num_part[ipart]);

      PDM_writer_var_set(wrt,
                         id_var_gnum,
                         id_geom,
                         ipart,
                         val_gnum[ipart]);
    }
    PDM_writer_var_write(wrt,
                         id_var_num_part);
    PDM_writer_var_free(wrt,
                        id_var_num_part);

    PDM_writer_var_write(wrt,
                         id_var_gnum);
    PDM_writer_var_free(wrt,
                        id_var_gnum);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free(val_num_part[ipart]);
      free(val_gnum    [ipart]);
    }
    free(val_num_part);
    free(val_gnum);

    PDM_writer_step_end(wrt);
    PDM_writer_free(wrt);
  }

  double l_mesh_extents[6] = {  HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                               -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};

  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_elt = PDM_Mesh_nodal_n_cell_get(ml->mesh_nodal, ipart);
    for (int i = 0; i < n_elt; i++) {
      for (int j = 0; j < 3; j++) {
        l_mesh_extents[j  ] = PDM_MIN(l_mesh_extents[j  ], elt_extents[ipart][6*i+j  ]);
        l_mesh_extents[j+3] = PDM_MAX(l_mesh_extents[j+3], elt_extents[ipart][6*i+j+3]);
      }
    }
  }

  double g_mesh_extents[6];
  PDM_MPI_Allreduce(l_mesh_extents,   g_mesh_extents,   3,
                    PDM_MPI_DOUBLE, PDM_MPI_MIN, ml->comm);
  PDM_MPI_Allreduce(l_mesh_extents+3, g_mesh_extents+3, 3,
                    PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

  if (dbg_enabled && i_rank == 0) {
    printf("g_mesh_extents = %f %f %f / %f %f %f\n",
           g_mesh_extents[0], g_mesh_extents[1], g_mesh_extents[2],
           g_mesh_extents[3], g_mesh_extents[4], g_mesh_extents[5]);
  }

  // PDM_MPI_Barrier (ml->comm);
  PDM_timer_hang_on(ml->timer);
  e_t_elapsed = PDM_timer_elapsed (ml->timer);
  e_t_cpu     = PDM_timer_cpu     (ml->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys (ml->timer);

  ml->times_elapsed[BUILD_BOUNDING_BOXES] += e_t_elapsed - b_t_elapsed;
  ml->times_cpu    [BUILD_BOUNDING_BOXES] += e_t_cpu     - b_t_cpu;
  ml->times_cpu_u  [BUILD_BOUNDING_BOXES] += e_t_cpu_u   - b_t_cpu_u;
  ml->times_cpu_s  [BUILD_BOUNDING_BOXES] += e_t_cpu_s   - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;
  PDM_timer_resume(ml->timer);

  /*
  *  Store cell vertex connectivity
  *  ------------------------------
  */

  _store_cell_vtx(ml);

  // PDM_MPI_Barrier (ml->comm);
  PDM_timer_hang_on(ml->timer);
  e_t_elapsed = PDM_timer_elapsed (ml->timer);
  e_t_cpu     = PDM_timer_cpu     (ml->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys (ml->timer);

  ml->times_elapsed[STORE_CONNECTIVITY] += e_t_elapsed - b_t_elapsed;
  ml->times_cpu    [STORE_CONNECTIVITY] += e_t_cpu     - b_t_cpu;
  ml->times_cpu_u  [STORE_CONNECTIVITY] += e_t_cpu_u   - b_t_cpu_u;
  ml->times_cpu_s  [STORE_CONNECTIVITY] += e_t_cpu_s   - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;
  PDM_timer_resume(ml->timer);

  double dt_search_all = 0.;
  double dt_extract_all = 0.;

  double dt_step1_all = 0.;
  double dt_step2_all = 0.;
  double dt_step3_all = 0.;
  double dt_step4_all = 0.;
  double dt_step5_all = 0.;
  double dt_step6_all = 0.;

  // PDM_mpi_comm_kind_t comm_kind = PDM_MPI_COMM_KIND_P2P;
  PDM_mpi_comm_kind_t comm_kind = PDM_MPI_COMM_KIND_COLLECTIVE;

  /* Big loop on point clouds */
  for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {

    // PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    b_t_elapsed = PDM_timer_elapsed(ml->timer);
    b_t_cpu     = PDM_timer_cpu(ml->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);
    PDM_timer_resume(ml->timer);

    if (dbg_enabled) {
      log_trace("Point cloud %d\n", icloud);
    }

    double t1 = PDM_MPI_Wtime();

    _point_cloud_t *pcloud = ml->point_clouds + icloud;

    if (dbg_enabled) {
      char name[999];
      sprintf(name, "mesh_location_point_cloud_%d", icloud);
      _dump_point_cloud(name,
                        ml->comm,
                        pcloud->n_part,
                        pcloud->n_points,
                        pcloud->coords,
                        pcloud->gnum,
                        0,
                        NULL,
                        NULL);
    }

    /*
     *  Extract points that intersect the source mesh global extents
     *  Brute force (could be accelerated using octree)
     */
    int  *n_select_pts     = malloc(sizeof(int  ) * pcloud->n_part);
    int **select_pts_l_num = malloc(sizeof(int *) * pcloud->n_part);

    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {

      n_select_pts    [ipart] = 0;
      select_pts_l_num[ipart] = malloc (sizeof(int) * pcloud->n_points[ipart]);

      for (int i = 0; i < pcloud->n_points[ipart]; i++) {
        int inside = 1;
        for (int j = 0; j < 3; j++) {
          if (pcloud->coords[ipart][3*i+j] < g_mesh_extents[j  ] ||
              pcloud->coords[ipart][3*i+j] > g_mesh_extents[j+3]) {
            inside = 0;
            break;
          }
        }

        if (inside) {
          select_pts_l_num[ipart][n_select_pts[ipart]++] = i;
        }
      } // End of loop on current parition's points

      select_pts_l_num[ipart] = realloc(select_pts_l_num[ipart],
                                        sizeof(int) * n_select_pts[ipart]);

    } // End of loop on current point cloud's partitions

    /* Compute global number of selected points */
    PDM_g_num_t l_n_pts[2] = {0, 0};
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      l_n_pts[0] += pcloud->n_points[ipart];
      l_n_pts[1] += n_select_pts    [ipart];
    }

    PDM_g_num_t g_n_pts[2];
    PDM_MPI_Allreduce(l_n_pts, g_n_pts, 2, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ml->comm);


    if (g_n_pts[1] == 0) {
      /* We extracted zero points, take a shortcut */
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        free(select_pts_l_num[ipart]);
      }
      free(n_select_pts);
      free(select_pts_l_num);


      /* Result from mesh PoV */
      if (ml->reverse_result) {
        _points_in_element_t *pts_in_elt = ml->points_in_elements + icloud;

        pts_in_elt->n_part = n_part;
        pts_in_elt->n_elts = malloc(sizeof(int) * n_part);
        memcpy(pts_in_elt->n_elts, pn_elt, sizeof(int) * n_part);

        pts_in_elt->pts_inside_idx   = malloc(sizeof(int         *) * n_part);
        pts_in_elt->gnum             = malloc(sizeof(PDM_g_num_t *) * n_part);
        pts_in_elt->coords           = malloc(sizeof(double      *) * n_part);
        pts_in_elt->uvw              = malloc(sizeof(double      *) * n_part);
        pts_in_elt->projected_coords = malloc(sizeof(double      *) * n_part);
        pts_in_elt->weights_idx      = malloc(sizeof(int         *) * n_part);
        pts_in_elt->weights          = malloc(sizeof(double      *) * n_part);
        pts_in_elt->dist2            = malloc(sizeof(double      *) * n_part);

        int **pts_in_elt_triplet_idx = malloc(sizeof(int *) * n_part);
        int **pts_in_elt_triplet     = malloc(sizeof(int *) * n_part);
        for (int ipart = 0; ipart < n_part; ipart++) {
          pts_in_elt->pts_inside_idx  [ipart] = PDM_array_zeros_int(pn_elt[ipart] + 1);
          pts_in_elt->gnum            [ipart] = malloc(sizeof(PDM_g_num_t) * 0);
          pts_in_elt->coords          [ipart] = malloc(sizeof(double     ) * 0);
          pts_in_elt->uvw             [ipart] = malloc(sizeof(double     ) * 0);
          pts_in_elt->projected_coords[ipart] = malloc(sizeof(double     ) * 0);
          pts_in_elt->weights_idx     [ipart] = malloc(sizeof(int        ) * 0);
          pts_in_elt->weights         [ipart] = malloc(sizeof(double     ) * 0);
          pts_in_elt->dist2           [ipart] = malloc(sizeof(double     ) * 0);

          pts_in_elt_triplet_idx[ipart] = PDM_array_zeros_int(1);
          pts_in_elt_triplet    [ipart] = malloc(sizeof(int) * 0);
        }


        ml->ptp[icloud] = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) elt_g_num,
                                                                    (const int          *) pn_elt,
                                                                    n_part,
                                                                    (const int          *) pcloud->n_points,
                                                                    pcloud->n_part,
                                                                    (const int         **) pts_in_elt->pts_inside_idx, // size = n_elt
                                                                    (const int         **) pts_in_elt_triplet_idx,     // size = pts_inside_idx[n_elt]
                                                                    (const int         **) pts_in_elt_triplet,
                                                                    ml->comm);
        for (int ipart = 0; ipart < n_part; ipart++) {
          free(pts_in_elt_triplet_idx[ipart]);
          free(pts_in_elt_triplet    [ipart]);
        }
        free(pts_in_elt_triplet_idx);
        free(pts_in_elt_triplet    );
      }

      /* Result from point cloud PoV */
      pcloud->n_located        = malloc(sizeof(int          ) * pcloud->n_part);
      pcloud->n_un_located     = malloc(sizeof(int          ) * pcloud->n_part);
      pcloud->located          = malloc(sizeof(int         *) * pcloud->n_part);
      pcloud->un_located       = malloc(sizeof(int         *) * pcloud->n_part);
      pcloud->location         = malloc(sizeof(PDM_g_num_t *) * pcloud->n_part);
      pcloud->projected_coords = malloc(sizeof(double      *) * pcloud->n_part);
      pcloud->dist2            = malloc(sizeof(double      *) * pcloud->n_part);
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        pcloud->n_located   [ipart] = 0;
        pcloud->n_un_located[ipart] = pcloud->n_points[ipart];

        pcloud->located   [ipart] = malloc(sizeof(int) * pcloud->n_located   [ipart]);
        pcloud->un_located[ipart] = malloc(sizeof(int) * pcloud->n_un_located[ipart]);

        for (int i = 0; i < pcloud->n_points[ipart]; i++) {
          pcloud->un_located[ipart][i] = i+1;
        }

        pcloud->location        [ipart] = malloc(sizeof(PDM_g_num_t) * 0);
        pcloud->projected_coords[ipart] = malloc(sizeof(double     ) * 0);
        pcloud->dist2           [ipart] = malloc(sizeof(double     ) * 0);
      }


      continue; // move on the next point cloud
    }



    int use_extracted_pts = (g_n_pts[1] < extraction_threshold * g_n_pts[0]);

    PDM_g_num_t **select_pts_g_num_user    = NULL;
    double      **select_pts_coord         = NULL;
    int         **select_pts_init_location = NULL;
    if (use_extracted_pts) {
      if (dbg_enabled) {
        log_trace("point cloud extraction %d / %d\n", l_n_pts[1], l_n_pts[0]);
      }

      select_pts_g_num_user    = malloc(pcloud->n_part * sizeof(PDM_g_num_t * ));
      select_pts_coord         = malloc(pcloud->n_part * sizeof(double      * ));
      select_pts_init_location = malloc(pcloud->n_part * sizeof(int         * ));

      // Just extract gnum
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        select_pts_g_num_user   [ipart] = malloc(    n_select_pts[ipart] * sizeof(PDM_g_num_t));
        select_pts_coord        [ipart] = malloc(3 * n_select_pts[ipart] * sizeof(double     ));
        select_pts_init_location[ipart] = malloc(3 * n_select_pts[ipart] * sizeof(int        ));

        for (int i = 0; i < n_select_pts[ipart]; i++) {
          int j = select_pts_l_num[ipart][i];

          select_pts_init_location[ipart][3*i  ] = i_rank;
          select_pts_init_location[ipart][3*i+1] = ipart;
          select_pts_init_location[ipart][3*i+2] = j;

          select_pts_g_num_user[ipart][i] = pcloud->gnum[ipart][j];
          for (int k = 0; k < 3; k++) {
            select_pts_coord[ipart][3*i + k] = pcloud->coords[ipart][3*j + k];
          }
        }
      }
    }
    else {
      /* We keep the whole point cloud */
      if (dbg_enabled) {
        log_trace("no point cloud extraction\n");
      }
      free(n_select_pts);
      n_select_pts            = pcloud->n_points;
      select_pts_g_num_user   = pcloud->gnum;
      select_pts_coord        = pcloud->coords;

      select_pts_init_location = malloc(pcloud->n_part * sizeof(int         * ));
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        select_pts_init_location[ipart] = malloc(3 * n_select_pts[ipart] * sizeof(int        ));
        for (int i = 0; i < n_select_pts[ipart]; i++) {
          select_pts_init_location[ipart][3*i  ] = i_rank;
          select_pts_init_location[ipart][3*i+1] = ipart;
          select_pts_init_location[ipart][3*i+2] = i;
        }
      }
    }

    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      free(select_pts_l_num[ipart]);
    }
    free(select_pts_l_num);


    /*
     *  Compute global extents of extracted point cloud
     */
    double l_pts_extents[6] = {
      HUGE_VAL, HUGE_VAL, HUGE_VAL,
      -HUGE_VAL, -HUGE_VAL, -HUGE_VAL
    };

    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      for (int i = 0; i < n_select_pts[ipart]; i++) {
        for (int j = 0; j < 3; j++) {
          l_pts_extents[j  ] = PDM_MIN(l_pts_extents[j  ], select_pts_coord[ipart][3*i + j]);
          l_pts_extents[j+3] = PDM_MAX(l_pts_extents[j+3], select_pts_coord[ipart][3*i + j]);
        }
      }
    }

    double g_pts_extents[6];
    PDM_MPI_Allreduce(l_pts_extents,   g_pts_extents,   3,
                      PDM_MPI_DOUBLE, PDM_MPI_MIN, ml->comm);
    PDM_MPI_Allreduce(l_pts_extents+3, g_pts_extents+3, 3,
                      PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

    if (dbg_enabled && i_rank == 0) {
      printf("  g_pts_extents = %f %f %f / %f %f %f\n",
             g_pts_extents[0], g_pts_extents[1], g_pts_extents[2],
             g_pts_extents[3], g_pts_extents[4], g_pts_extents[5]);
    }



    /*
     *  Select elements whose bounding box intersect
     *  the global extents of the extracted point clouds
     *  Brute force (could be accelerated using bbtree)
     */
    int  *n_select_elt     = malloc(sizeof(int  ) * n_part);
    int **select_elt_l_num = malloc(sizeof(int *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      int n_elt = PDM_Mesh_nodal_n_cell_get(ml->mesh_nodal,
                                            ipart);

      n_select_elt[ipart] = 0;
      select_elt_l_num[ipart] = malloc(sizeof(int) * n_elt);

      for (int ielt = 0; ielt < n_elt; ielt++) {

        double *box_min = elt_extents[ipart] + 6*ielt;
        double *box_max = box_min + 3;

        int intersect = 1;
        for (int j = 0; j < 3; j++) {
          if (box_min[j] > g_pts_extents[j+3] ||
              box_max[j] < g_pts_extents[j]) {
            intersect = 0;
            break;
          }
        }

        // if (intersect &&
        //     ml->is_elmt_select_by_user[ipart]       !=  NULL &&
        //     ml->is_elmt_select_by_user[ipart][ielt] == 1) {
        //   select_elt_l_num[ipart][n_select_elt[ipart]++] = ielt;
        // }
        if (intersect) {
          select_elt_l_num[ipart][n_select_elt[ipart]++] = ielt;
        }

      } // End of loop on current part's boxes

      select_elt_l_num[ipart] = realloc(select_elt_l_num[ipart],
                                        sizeof(int) * n_select_elt[ipart]);

    } // End of loop on mesh parts

    /* Compute global number of selected elements */
    PDM_g_num_t l_n_elt[2] = {0, 0};
    for (int ipart = 0; ipart < n_part; ipart++) {
      int n_elt = PDM_Mesh_nodal_n_cell_get(ml->mesh_nodal,
                                            ipart);
      l_n_elt[0] += n_elt;
      l_n_elt[1] += n_select_elt[ipart];
    }

    PDM_g_num_t g_n_elt[2];
    PDM_MPI_Allreduce(l_n_elt, g_n_elt, 2, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ml->comm);

    int use_extracted_mesh = (g_n_elt[1] < extraction_threshold * g_n_elt[0]);

    int         **select_elt_init_location_user = NULL;
    double      **select_elt_extents            = NULL;
    PDM_g_num_t **select_elt_g_num_user         = NULL;

    if (use_extracted_mesh) {
      if (dbg_enabled) {
        // log_trace("mesh extraction %d / %d\n", l_n_elt[1], l_n_elt[0]);
        log_trace("mesh extraction %d / %d\n", g_n_elt[1], g_n_elt[0]);
      }

      select_elt_init_location_user = malloc(sizeof(int         *) * n_part);
      select_elt_extents            = malloc(sizeof(double      *) * n_part);
      select_elt_g_num_user         = malloc(sizeof(PDM_g_num_t *) * n_part);
      for (int ipart = 0; ipart < n_part; ipart++) {
        select_elt_init_location_user[ipart] = malloc(sizeof(int        ) * n_select_elt[ipart] * 3);
        select_elt_extents           [ipart] = malloc(sizeof(double     ) * n_select_elt[ipart] * 6);
        select_elt_g_num_user        [ipart] = malloc(sizeof(PDM_g_num_t) * n_select_elt[ipart]);
        for (int i = 0; i < n_select_elt[ipart]; i++) {
          int elt_id = select_elt_l_num[ipart][i];

          memcpy(select_elt_extents[ipart] + 6*i,
                 elt_extents       [ipart] + 6*elt_id,
                 sizeof(double) * 6);

          select_elt_g_num_user[ipart][i] = elt_g_num[ipart][elt_id];

          select_elt_init_location_user[ipart][3*i  ] = i_rank;
          select_elt_init_location_user[ipart][3*i+1] = ipart;
          select_elt_init_location_user[ipart][3*i+2] = elt_id;

        }
      }
    }
    else {
      if (dbg_enabled) {
        log_trace("no mesh extraction\n");
      }
      free(n_select_elt);
      n_select_elt          = pn_elt;
      select_elt_extents    = elt_extents;
      select_elt_g_num_user = elt_g_num;

      select_elt_init_location_user = malloc(sizeof(int*) * n_part);
      for (int ipart = 0; ipart < n_part; ipart++) {
        select_elt_init_location_user[ipart] = malloc(sizeof(int) * n_select_elt[ipart] * 3);
        for (int i = 0; i < n_select_elt[ipart]; i++) {
          select_elt_init_location_user[ipart][3*i  ] = i_rank;
          select_elt_init_location_user[ipart][3*i+1] = ipart;
          select_elt_init_location_user[ipart][3*i+2] = i;
        }
      }
    }

    /*
     *  Redistribute evenly the selected points (ptb_geom)
     */
    int **weight      = malloc(sizeof(int *) * pcloud->n_part);
    int **pstride_one = malloc(sizeof(int *) * pcloud->n_part);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      weight     [ipart] = PDM_array_const_int(n_select_pts[ipart], 1);
      pstride_one[ipart] = PDM_array_const_int(n_select_pts[ipart], 1);
    }

    /* Use same extents for Hilbert encoding of pts and boxes?? */
    PDM_part_to_block_t *ptb_pts = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                 PDM_PART_TO_BLOCK_POST_MERGE, // TODO: Merge
                                                                 1.,
                                                                 PDM_PART_GEOM_HILBERT,
                                                                 select_pts_coord,
                                                                 select_pts_g_num_user,
                                                                 weight,
                                                                 n_select_pts,
                                                                 pcloud->n_part,
                                                                 ml->comm);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      free(weight[ipart]);
    }
    free(weight);

    int dn_pts = PDM_part_to_block_n_elt_block_get(ptb_pts);
    PDM_g_num_t *dpts_g_num_user = PDM_part_to_block_block_gnum_get(ptb_pts);


    PDM_g_num_t *distrib_pts = PDM_part_to_block_distrib_index_get(ptb_pts);
    PDM_g_num_t *dpts_g_num_geom = malloc(dn_pts * sizeof(PDM_g_num_t));

    for(int i = 0; i < dn_pts; ++i) {
      dpts_g_num_geom[i] = distrib_pts[i_rank] + i + 1;
    }

    if (dbg_enabled) {
      PDM_log_trace_array_long(dpts_g_num_user, dn_pts, "dpts_g_num_user        : ");
    }

    /* Exchange coordinates (do this with abstract distrib?) */
    int    *blk_coord_n       = NULL;
    double *tmp_blk_pts_coord = NULL;
    int request_pts_coord = -1;
    PDM_part_to_block_iexch(ptb_pts,
                            PDM_MPI_COMM_KIND_COLLECTIVE,
                            3 * sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            pstride_one,
                  (void **) select_pts_coord,
                            &blk_coord_n,
                  (void **) &tmp_blk_pts_coord,
                            &request_pts_coord);

    /*
     * Exchange init_location of pts --> Is almost the selected with proc / part info
     */
    int *dpts_init_location_pts_n = NULL;
    int *dpts_init_location_pts   = NULL;
    PDM_part_to_block_exch(ptb_pts,
                           3 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           pstride_one,
                 (void **) select_pts_init_location,
                           &dpts_init_location_pts_n,
                 (void **) &dpts_init_location_pts);

    PDM_part_to_block_iexch_wait(ptb_pts, request_pts_coord);

    double *dpts_coord = malloc(3 * dn_pts * sizeof(double));
    int idx_read  = 0;
    for(int i = 0; i < dn_pts; ++i) {
      dpts_coord[3*i  ] = tmp_blk_pts_coord[3*idx_read  ];
      dpts_coord[3*i+1] = tmp_blk_pts_coord[3*idx_read+1];
      dpts_coord[3*i+2] = tmp_blk_pts_coord[3*idx_read+2];

      idx_read += blk_coord_n[i];
    }
    free(blk_coord_n);
    free(tmp_blk_pts_coord);

    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      free(pstride_one[ipart]);
    }
    free(pstride_one);

    if (dbg_enabled) {
      char filename[999];
      sprintf(filename, "mesh_location_dpts_%d_%3.3d.vtk", icloud, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                dn_pts,
                                dpts_coord,
                                dpts_g_num_user,
                                NULL);
    }

    /*
     *  Redistribute evenly the selected boxes (ptb_geom)
     */
    double **select_box_center = malloc(sizeof(double *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      select_box_center[ipart] = malloc(sizeof(double) * n_select_elt[ipart] * 3);
      for (int i = 0; i < n_select_elt[ipart]; i++) {
        for (int j = 0; j < 3; j++) {
          select_box_center[ipart][3*i+j] = 0.5*(select_elt_extents[ipart][6*i+j  ] +
                                                 select_elt_extents[ipart][6*i+j+3]);
        }
      }
    }

    weight = malloc(sizeof(int *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      weight[ipart] = PDM_array_const_int(n_select_elt[ipart], 1);
    }

    if(0){
      char filename[999];
      sprintf(filename, "mesh_location_extract_boxes_%d_%3.3d.vtk", icloud, i_rank);
      PDM_vtk_write_boxes(filename,
                          n_select_elt[0],
                          select_elt_extents[0],
                          select_elt_g_num_user[0]);
    }

    PDM_part_to_block_t *ptb_elt = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                 PDM_PART_TO_BLOCK_POST_CLEANUP, // A faire en merge
                                                                 1.,
                                                                 PDM_PART_GEOM_HILBERT,
                                                                 select_box_center,
                                                                 select_elt_g_num_user,
                                                                 weight,
                                                                 n_select_elt,
                                                                 n_part,
                                                                 ml->comm);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free(select_box_center[ipart]);
      free(weight[ipart]);
    }
    free(select_box_center);
    free(weight);


    int dn_elt1 = PDM_part_to_block_n_elt_block_get(ptb_elt);
    PDM_g_num_t *delt_g_num_user = PDM_part_to_block_block_gnum_get(ptb_elt);
    if (dbg_enabled) {
      PDM_log_trace_array_long(delt_g_num_user, dn_elt1, "delt_g_num_user : ");
    }

    PDM_g_num_t *distrib_elt1 = PDM_part_to_block_distrib_index_get(ptb_elt);
    PDM_g_num_t *delmt_g_num_geom = malloc(dn_elt1 * sizeof(PDM_g_num_t));

    for(int i = 0; i < dn_elt1; ++i) {
      delmt_g_num_geom[i] = distrib_elt1[i_rank] + i + 1;
    }


    /* Exchange extents (do this with abstract distrib?) */
    double *delt_extents1 = NULL;
    int request_elt_extents = -1;
    PDM_part_to_block_iexch(ptb_elt,
                            PDM_MPI_COMM_KIND_COLLECTIVE,
                            6 * sizeof(double),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                  (void **) select_elt_extents,
                            NULL,
                  (void **) &delt_extents1,
                            &request_elt_extents);

    /*
     * TODO - Exchange box_init_location en merge  !!! Attention on doit merger les blocks mais garder le non mergé pour le renvoie
     */
    double *delt_init_location_user = NULL;
    int request_elt_init_location = -1;
    PDM_part_to_block_iexch(ptb_elt,
                            PDM_MPI_COMM_KIND_COLLECTIVE,
                            3 * sizeof(int),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                  (void **) select_elt_init_location_user,
                            NULL,
                  (void **) &delt_init_location_user,
                            &request_elt_init_location);
    PDM_part_to_block_iexch_wait(ptb_elt, request_elt_extents);
    PDM_part_to_block_iexch_wait(ptb_elt, request_elt_init_location);


    if (dbg_enabled) {
      char filename[999];
      sprintf(filename, "mesh_location_dboxes_%d_%3.3d.vtk", icloud, i_rank);
      PDM_vtk_write_boxes(filename,
                          dn_elt1,
                          delt_extents1,
                          delt_g_num_user);
    }

    /*
     *  Location : search candidates
     *  ----------------------------
     */
    double t2 = PDM_MPI_Wtime();
    if (log_timer) {
      log_trace("Extract + hilbert  = %12.5e \n", t2 -t1);
    }

    dt_extract_all += t2 -t1;

    /*
     * Get points inside bounding boxes of elements
     */

    int         *box_pts_idx   = NULL;
    PDM_g_num_t *box_pts_g_num = NULL;
    double      *box_pts_coord = NULL;

    int dn_elt2 = 0;
    PDM_g_num_t *delt_g_num_geom2 = NULL;
    int         *delt_pts_n2      = NULL;
    PDM_g_num_t *delt_pts_g_num_geom  = NULL;
    double      *delt_pts_coord2  = NULL;

    char *env_var_oct = getenv("OCTREE_SHARED");
    int use_shared_tree = 0;
    if (env_var_oct != NULL) {
      use_shared_tree = atoi(env_var_oct);
    }

    switch (ml->method) {
      case PDM_MESH_LOCATION_OCTREE: {
        /* Create octree structure */
        PDM_para_octree_t *octree = PDM_para_octree_create(1,
                                                           octree_depth_max,
                                                           octree_points_in_leaf_max,
                                                           octree_build_leaf_neighbours,
                                                           ml->comm);

        /* Set octree point cloud */
        PDM_para_octree_point_cloud_set(octree,
                                        0,
                                        dn_pts,
                                        dpts_coord,
                                        dpts_g_num_geom);

        /* Build parallel octree */
        // PDM_MPI_Barrier(ml->comm);

        if (use_shared_tree == 0) {
          PDM_para_octree_build(octree, NULL);
        }
        else {
          PDM_para_octree_build_shared(octree, NULL);
        }

        if (use_shared_tree == 0) {
          PDM_part_to_block_t *ptb_pib = NULL;
          PDM_para_octree_points_inside_boxes_block_frame(octree,
                                                          dn_elt1,
                                                          delt_extents1,
                                                          delmt_g_num_geom,
                                                          &ptb_pib,
                                                          &delt_pts_n2,
                                                          &delt_pts_g_num_geom,
                                                          &delt_pts_coord2);

          dn_elt2 = PDM_part_to_block_n_elt_block_get(ptb_pib);
          PDM_g_num_t *_g_num = PDM_part_to_block_block_gnum_get(ptb_pib);

          delt_g_num_geom2 = malloc(sizeof(PDM_g_num_t) * dn_elt2);
          memcpy(delt_g_num_geom2, _g_num, sizeof(PDM_g_num_t) * dn_elt2);

          PDM_part_to_block_free(ptb_pib);

        }
        else {
          PDM_part_to_block_t *ptb_pib = NULL;
          PDM_para_octree_points_inside_boxes_shared_block_frame(octree,
                                                                 dn_elt1,
                                                                 delt_extents1,
                                                                 delmt_g_num_geom,
                                                                 &ptb_pib,
                                                                 &delt_pts_n2,
                                                                 &delt_pts_g_num_geom,
                                                                 &delt_pts_coord2);
          dn_elt2 = PDM_part_to_block_n_elt_block_get(ptb_pib);
          PDM_g_num_t *_g_num = PDM_part_to_block_block_gnum_get(ptb_pib);

          delt_g_num_geom2 = malloc(sizeof(PDM_g_num_t) * dn_elt2);
          memcpy(delt_g_num_geom2, _g_num, sizeof(PDM_g_num_t) * dn_elt2);

          PDM_part_to_block_free(ptb_pib);
        }

        PDM_para_octree_free(octree);
        break;
      }
      case PDM_MESH_LOCATION_DOCTREE: {

        // PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_KDTREE;
        PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_OCTREE;
        PDM_doctree_t *doct = PDM_doctree_create(ml->comm,
                                                 3,
                                                 1,
                                                 NULL, // global_extents
                                                 local_tree_kind);

        /* pass abstract distrib of pts to doctree? */
        int *init_location_pts = NULL;
        if (0) {
          init_location_pts = malloc(3 * dn_pts * sizeof(int));
          for(int i = 0; i < dn_pts; ++i) {
            init_location_pts[3*i  ] = i_rank;
            init_location_pts[3*i+1] = 0;
            init_location_pts[3*i+2] = i;
          }
        }
        PDM_doctree_point_set(doct,
                              0,
                              dn_pts,
                              init_location_pts,
                              dpts_g_num_geom,
                              dpts_coord);

        /* pass abstract distrib of boxes to doctree? */
        /* or get init location from ad_elt? */
        int *init_location_box = malloc(3 * dn_elt1 * sizeof(int));
        for(int i = 0; i < dn_elt1; ++i) {
          init_location_box[3*i  ] = i_rank;
          init_location_box[3*i+1] = 0;
          init_location_box[3*i+2] = i;
        }

        PDM_doctree_solicitation_set(doct,
                                     PDM_TREE_SOLICITATION_BOXES_POINTS,
                                     1,
                                     &dn_elt1,
                                     &init_location_box,
                                     &delmt_g_num_geom,
                                     &delt_extents1);

        PDM_doctree_build(doct);

        PDM_doctree_results_in_block_frame_get(doct,
                                               &dn_elt2,
                                               &delt_g_num_geom2,
                                               &delt_pts_n2,
                                               &delt_pts_g_num_geom,
                                               &delt_pts_coord2,
                                               PDM_OWNERSHIP_USER);


        PDM_doctree_dump_times(doct);
        PDM_doctree_free(doct);
        if (init_location_pts != NULL) {
          free(init_location_pts);
        }
        free(init_location_box);
        break;
      }
      case PDM_MESH_LOCATION_DBBTREE: {
        /* Compute local extents */
        double l_extents[6] = {
          HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
          -HUGE_VAL, -HUGE_VAL, -HUGE_VAL
        };

        for (int ipart = 0; ipart < n_part; ipart++) {
          for (int i = 0; i < n_select_elt[ipart]; i++) {
            for (int j = 0; j < 3; j++) {
              l_extents[j  ] = PDM_MIN(l_extents[j  ], select_elt_extents[ipart][6*i+j]);
              l_extents[j+3] = PDM_MAX(l_extents[j+3], select_elt_extents[ipart][6*i+3+j]);
            }
          }
        }

        /* Compute global extents */
        double g_extents[6];
        PDM_MPI_Allreduce(l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, ml->comm);
        PDM_MPI_Allreduce(l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

        /* Break symmetry */
        double max_range = 1e-12;
        for (int i = 0; i < 3; i++) {
          max_range = PDM_MAX(max_range, g_extents[i+3] - g_extents[i]);
        }
        for (int i = 0; i < 3; i++) {
          g_extents[i]   -= max_range * 1.1e-3;
          g_extents[i+3] += max_range * 1.0e-3;
        }

        PDM_dbbtree_t *dbbt = PDM_dbbtree_create(ml->comm, dim, g_extents);

        PDM_box_set_t *box_set = PDM_dbbtree_boxes_set(dbbt,
                                                       1,
                                                       &dn_elt1,
                                (const double      **) &delt_extents1,
                                (const PDM_g_num_t **) &delmt_g_num_geom);


        if (use_shared_tree == 0) {

          PDM_part_to_block_t *ptb_pib = NULL;
          PDM_dbbtree_points_inside_boxes_block_frame(dbbt,
                                                      dn_pts,
                                                      dpts_g_num_geom,
                                                      dpts_coord,
                                                      &ptb_pib,
                                                      &delt_pts_n2,
                                                      &delt_pts_g_num_geom,
                                                      &delt_pts_coord2,
                                                      0);

          dn_elt2 = PDM_part_to_block_n_elt_block_get(ptb_pib);
          PDM_g_num_t *_g_num = PDM_part_to_block_block_gnum_get(ptb_pib);

          delt_g_num_geom2 = malloc(sizeof(PDM_g_num_t) * dn_elt2);
          memcpy(delt_g_num_geom2, _g_num, sizeof(PDM_g_num_t) * dn_elt2);

          PDM_part_to_block_free(ptb_pib);

        }
        else {
          abort();
          // PDM_MPI_Barrier (ml->comm);
          PDM_dbbtree_points_inside_boxes_shared(dbbt,
                                                 dn_pts,
                                                 dpts_g_num_geom,
                                                 dpts_coord,
                                                 dn_elt1,
                                                 delmt_g_num_geom, // Attention faire une distribution part_to_bloc_geom dans le cas octree
                                                 &box_pts_idx,
                                                 &box_pts_g_num,
                                                 &box_pts_coord,
                                                 0);
        }

        PDM_dbbtree_free(dbbt);
        PDM_box_set_destroy(&box_set);
        break;
      }
      default: {
        PDM_error(__FILE__, __LINE__, 0,
                  "PDM_mesh_location : unknown location method %d\n", (int) ml->method);
      }
    }
    free(delt_extents1);
    free(dpts_coord);
    free(delmt_g_num_geom);
    free(dpts_g_num_geom);


    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed (ml->timer);
    e_t_cpu     = PDM_timer_cpu     (ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys (ml->timer);

    ml->times_elapsed[SEARCH_CANDIDATES] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [SEARCH_CANDIDATES] += e_t_cpu     - b_t_cpu;
    ml->times_cpu_u  [SEARCH_CANDIDATES] += e_t_cpu_u   - b_t_cpu_u;
    ml->times_cpu_s  [SEARCH_CANDIDATES] += e_t_cpu_s   - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);
    PDM_MPI_Barrier (ml->comm);
    if (log_timer) {
      log_trace("Search  = %12.5e \n", PDM_MPI_Wtime() - t2);
    }

    dt_search_all += PDM_MPI_Wtime() - t2;


    if (dbg_enabled) {
      log_trace("before compression\n");
      PDM_log_trace_array_long(delt_g_num_geom2, dn_elt2, "delt_g_num_geom2 : ");
    }

    /* TODO: maybe take into account element type to yield better weights ?? */
    if (dbg_enabled) {
      int *_idx = PDM_array_new_idx_from_sizes_int(delt_pts_n2, dn_elt2);
      PDM_log_trace_connectivity_long(_idx,
                                      delt_pts_g_num_geom,
                                      dn_elt2,
                                      "delt_pts_g_num_geom : ");
      free(_idx);
    }

    /* Compress block (remove elements with zero candidate point) -> frame (2) */
    int tmp_dn_elt2 = 0;
    for (int i = 0; i < dn_elt2; i++) {
      if (delt_pts_n2[i] > 0) {
        delt_pts_n2[tmp_dn_elt2] = delt_pts_n2[i];
        delt_g_num_geom2[tmp_dn_elt2] = delt_g_num_geom2[i];
        tmp_dn_elt2++;
      }
    }
    dn_elt2 = tmp_dn_elt2;

    delt_g_num_geom2 = realloc(delt_g_num_geom2, sizeof(PDM_g_num_t) * dn_elt2);
    delt_pts_n2      = realloc(delt_pts_n2     , sizeof(int        ) * dn_elt2);

    if (dbg_enabled) {
      log_trace("after compression\n");
      PDM_log_trace_array_long(delt_g_num_geom2, dn_elt2, "delt_g_num_geom2 : ");
      int total_weight = 0;
      for (int i = 0; i < dn_elt2; i++) {
        total_weight += delt_pts_n2[i];
      }
      log_trace("total weight = %d\n", total_weight);
    }

    int *delt_pts_idx2 = PDM_array_new_idx_from_sizes_int(delt_pts_n2, dn_elt2);
    free(delt_pts_n2);

    /*
     * Comment faire le lien avec le extract part ?
     *    delt_g_num_geom2 = Numero geometrique mais requilibré fonction de la soliciation
     *    Je pense que c'est presque un block_to_block ...
     */

    PDM_block_to_part_t* btp_elmt_geom_to_elmt_user = PDM_block_to_part_create(distrib_elt1,
                                                        (const PDM_g_num_t **) &delt_g_num_geom2,
                                                                               &dn_elt2,
                                                                               1,
                                                                               ml->comm);

    /*
     * Exchange gnum
     */
    PDM_g_num_t **tmp_delt_parent_g_num2 = NULL;
    int stride_one = 1;
    PDM_block_to_part_exch(btp_elmt_geom_to_elmt_user,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           &stride_one,
                           delt_g_num_user,
                           NULL,
            (void ***)    &tmp_delt_parent_g_num2);
    PDM_g_num_t *delt_parent_g_num2 = tmp_delt_parent_g_num2[0];
    free(tmp_delt_parent_g_num2);

    // TODO : Adpat when we merge elmt and take : delt_init_location_user_unified
    int **tmp_delt_init_location2 = NULL;
    PDM_block_to_part_exch(btp_elmt_geom_to_elmt_user,
                           3 * sizeof(int),
                           PDM_STRIDE_CST_INTERLACED,
                           &stride_one,
                           delt_init_location_user,
                           NULL,
         (void ***)        &tmp_delt_init_location2);
    int *delt_init_location2 = tmp_delt_init_location2[0];
    free(tmp_delt_init_location2);
    free(delt_init_location_user);

    PDM_block_to_part_free(btp_elmt_geom_to_elmt_user);

    if(dbg_enabled) {
      PDM_log_trace_array_long(delt_parent_g_num2, dn_elt2, "delt_parent_g_num2 ::");
      PDM_log_trace_array_int (delt_init_location2, 3 * dn_elt2, "delt_init_location2 ::");
    }

    /* Extract partition associated to current frame (2) */
    if (dbg_enabled) {
      log_trace(">> PDM_extract_part_create\n");
    }
    // TODO: proper get
    int mesh_dimension = pmne->mesh_dimension;

    PDM_extract_part_t *extrp = PDM_extract_part_create(mesh_dimension,
                                                        n_part,
                                                        1,                                 // n_part_out
                                                        PDM_EXTRACT_PART_KIND_FROM_TARGET,
                                                        PDM_SPLIT_DUAL_WITH_PTSCOTCH,      // unused in this case
                                                        PDM_FALSE,                         // compute_child_gnum
                                                        PDM_OWNERSHIP_KEEP,
                                                        ml->comm);

    PDM_extract_part_part_nodal_set(extrp, pmne);

    /* Set vtx_coord */
    for (int ipart = 0; ipart < n_part; ipart++) {
      const double *pvtx_coord = PDM_Mesh_nodal_vertices_get(ml->mesh_nodal,
                                                             ipart);

      PDM_g_num_t *pvtx_ln_to_gn = NULL;
      const int pn_vtx = PDM_Mesh_nodal_vertices_ln_to_gn_get(ml->mesh_nodal,
                                                              ipart,
                                                              &pvtx_ln_to_gn);

      int n_cell = 0;
      int n_face = 0;
      int n_edge = 0;
      PDM_g_num_t *cell_g_num = NULL;
      PDM_g_num_t *face_g_num = NULL;
      PDM_g_num_t *edge_g_num = NULL;
      switch (mesh_dimension) {
        case 3: {
          n_cell = pn_elt[ipart];
          cell_g_num = elt_g_num[ipart];
          break;
        }
        case 2: {
          n_face = pn_elt[ipart];
          face_g_num = elt_g_num[ipart];
          break;
        }
        case 1: {
          n_edge = pn_elt[ipart];
          edge_g_num = elt_g_num[ipart];
          break;
        }
        default:
        PDM_error(__FILE__, __LINE__, 0, "incorrect mesh_dimension %d\n", mesh_dimension);
      }

      PDM_extract_part_part_set(extrp,
                                ipart,
                                n_cell,
                                n_face,
                                n_edge,
                                pn_vtx,
                                NULL, // pcell_face_idx[ipart],
                                NULL, // pcell_face[ipart],
                                NULL, // pface_edge_idx[ipart],
                                NULL, // pface_edge[ipart],
                                NULL, // pedge_vtx[ipart],
                                NULL, // pface_vtx_idx[ipart],
                                NULL, // pface_vtx[ipart],
                                cell_g_num,
                                face_g_num,
                                edge_g_num,
                                pvtx_ln_to_gn,
                     (double *) pvtx_coord);
    }


    PDM_mesh_entities_t entity_type = PDM_MESH_ENTITY_CELL;
    switch (mesh_dimension) {
      case 3: {
        entity_type = PDM_MESH_ENTITY_CELL;
        break;
      }
      case 2: {
        entity_type = PDM_MESH_ENTITY_FACE;
        break;
      }
      case 1: {
        entity_type = PDM_MESH_ENTITY_EDGE;
        break;
      }
      default:
      PDM_error(__FILE__, __LINE__, 0, "incorrect mesh_dimension %d\n", mesh_dimension);
    }
    PDM_extract_part_target_set(extrp,
                                0,
                                dn_elt2,
                                delt_parent_g_num2,
                                delt_init_location2);

    if (dbg_enabled) {
      log_trace(">> PDM_extract_part_compute\n");
    }
    PDM_extract_part_compute(extrp);
    if (dbg_enabled) {
      log_trace("Yeah :D\n");
    }

    PDM_part_mesh_nodal_elmts_t *extract_pmne = NULL;
    PDM_extract_part_part_mesh_nodal_get(extrp,
                                         &extract_pmne,
                                         PDM_OWNERSHIP_USER);

    PDM_part_to_part_t *ptp_elt = NULL;
    PDM_extract_part_part_to_part_get(extrp,
                                      entity_type,
                                      &ptp_elt,
                                      PDM_OWNERSHIP_USER);

    int          pextract_n_elt        = 0;
    int          pextract_n_vtx        = 0;
    double      *pextract_vtx_coord    = NULL;

    pextract_n_elt = PDM_extract_part_n_entity_get(extrp,
                                                   0,
                                                   entity_type);
    assert(pextract_n_elt == dn_elt2);

    pextract_n_vtx = PDM_extract_part_n_entity_get(extrp,
                                                   0,
                                                   PDM_MESH_ENTITY_VERTEX);

    PDM_extract_part_vtx_coord_get(extrp,
                                   0,
                                   &pextract_vtx_coord,
                                   PDM_OWNERSHIP_KEEP);

    if (dbg_enabled) {
      _dump_pmne(ml->comm,
                 "extract_pmne",
                 extract_pmne,
                 1,
                 &delt_g_num_geom2,
                 &pextract_n_vtx,
                 &pextract_vtx_coord);
    }


    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed (ml->timer);
    e_t_cpu     = PDM_timer_cpu     (ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys (ml->timer);

    ml->times_elapsed[LOAD_BALANCING] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [LOAD_BALANCING] += e_t_cpu     - b_t_cpu;
    ml->times_cpu_u  [LOAD_BALANCING] += e_t_cpu_u   - b_t_cpu_u;
    ml->times_cpu_s  [LOAD_BALANCING] += e_t_cpu_s   - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);
    // PDM_MPI_Barrier (ml->comm);

    /* Perform elementary point locations */
    double **pelt_pts_distance2   = NULL;
    double **pelt_pts_proj_coord2 = NULL;
    int    **pelt_pts_weight_idx2 = NULL;
    double **pelt_pts_weight2     = NULL;
    double **pelt_pts_uvw2        = NULL;

    PDM_point_location_nodal2(extract_pmne,
                              1,
            (const double **) &pextract_vtx_coord,
            (const int    **) &delt_pts_idx2,
            (const double **) &delt_pts_coord2,
                              tolerance,
                              &pelt_pts_distance2,
                              &pelt_pts_proj_coord2,
                              &pelt_pts_weight_idx2,
                              &pelt_pts_weight2,
                              &pelt_pts_uvw2);

    double *delt_pts_distance2   = pelt_pts_distance2  [0];
    double *delt_pts_proj_coord2 = pelt_pts_proj_coord2[0];
    int    *delt_pts_weight_idx2 = pelt_pts_weight_idx2[0];
    double *delt_pts_weight2     = pelt_pts_weight2    [0];
    double *delt_pts_uvw2        = pelt_pts_uvw2       [0];
    free(pelt_pts_distance2  );
    free(pelt_pts_proj_coord2);
    free(pelt_pts_weight_idx2);
    free(pelt_pts_weight2    );
    free(pelt_pts_uvw2       );

    PDM_part_mesh_nodal_elmts_free(extract_pmne);
    PDM_extract_part_free(extrp);
    free(delt_init_location2);

    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[COMPUTE_ELEMENTARY_LOCATIONS] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [COMPUTE_ELEMENTARY_LOCATIONS] += e_t_cpu - b_t_cpu;
    ml->times_cpu_u  [COMPUTE_ELEMENTARY_LOCATIONS] += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s  [COMPUTE_ELEMENTARY_LOCATIONS] += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);
    // PDM_MPI_Barrier (ml->comm);

    double beg_step1 = PDM_MPI_Wtime();

    /*
     *  Pass in block frame for points to select closest element
     */
    int n_pts2 = delt_pts_idx2[dn_elt2];
    int         *part_stride = PDM_array_const_int(n_pts2, 1);
    PDM_g_num_t *part_elt_id = malloc(sizeof(PDM_g_num_t) * n_pts2);
    // double      *part_weight = malloc(sizeof(double     ) * n_pts2);
    for (int ielt = 0; ielt < dn_elt2; ielt++) {
      for (int i = delt_pts_idx2[ielt]; i < delt_pts_idx2[ielt+1]; i++) {
        part_elt_id[i] = delt_g_num_geom2[ielt];
        // part_weight[i] = 1.;
      }
    }

    PDM_g_num_t *pts_ln_to_gn = malloc(sizeof(PDM_g_num_t) * n_pts2);
    int *pts_unique_order = malloc(sizeof(int) * n_pts2);
    memcpy(pts_ln_to_gn, delt_pts_g_num_geom, sizeof(PDM_g_num_t) * n_pts2);
    int n_pts_unique = PDM_inplace_unique_long2(pts_ln_to_gn,
                                                pts_unique_order,
                                                0,
                                                n_pts2-1);
    pts_ln_to_gn = realloc(pts_ln_to_gn, sizeof(PDM_g_num_t) * n_pts_unique);
    if (dbg_enabled) {
      log_trace("%d unique pts / %d\n", n_pts_unique, n_pts2);
    }

    PDM_g_num_t *local_pts_elt_g_num = malloc(sizeof(PDM_g_num_t) * n_pts_unique);
    double      *local_pts_elt_dist2 = malloc(sizeof(double     ) * n_pts_unique);
    double      *part_weight         = malloc(sizeof(double     ) * n_pts_unique);
    for (int i = 0; i < n_pts_unique; i++) {
      local_pts_elt_dist2[i] = HUGE_VAL;
      part_weight[i] = 1.;
    }

    for (int ielt = 0; ielt < dn_elt2; ielt++) {
      for (int i = delt_pts_idx2[ielt]; i < delt_pts_idx2[ielt+1]; i++) {
        int pt_id = pts_unique_order[i];
        if (delt_pts_distance2[i] < local_pts_elt_dist2[pt_id]) {
          local_pts_elt_g_num[pt_id] = delt_g_num_geom2[ielt];
          local_pts_elt_dist2[pt_id] = delt_pts_distance2[i];
        }
      }
    }


    if (dbg_enabled) {
      for (int i = 0; i < n_pts2; i++) {
        log_trace("pt "PDM_FMT_G_NUM" (%f %f %f) : elt "PDM_FMT_G_NUM", dist = %e\n",
                  delt_pts_g_num_geom[i],
                  delt_pts_coord2[3*i], delt_pts_coord2[3*i+1], delt_pts_coord2[3*i+2],
                  part_elt_id[i], delt_pts_distance2[i]);
      }
    }

    if (dbg_enabled) {
      for (int i = 0; i < n_pts_unique; i++) {
        log_trace("  pt "PDM_FMT_G_NUM" : elt "PDM_FMT_G_NUM", dist = %e\n",
                  pts_ln_to_gn[i], local_pts_elt_g_num[i], local_pts_elt_dist2[i]);
      }
    }

    // PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
    //                                                     PDM_PART_TO_BLOCK_POST_MERGE,
    //                                                     1.,
    //                                                     &delt_pts_g_num_geom,
    //                                                     &part_weight,
    //                                                     &n_pts2,
    //                                                     1,
    //                                                     ml->comm);
    PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                        1.,
                                                        &pts_ln_to_gn,
                                                        &part_weight,
                                                        &n_pts_unique,
                                                        1,
                                                        ml->comm);
    free(part_weight);

    int    *block_pts_elt_n     = NULL;
    double *block_pts_elt_dist2 = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(double),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           &part_stride,
                 // (void **) &delt_pts_distance2,
                 (void **) &local_pts_elt_dist2,
                           &block_pts_elt_n,
                 (void **) &block_pts_elt_dist2);
    free(block_pts_elt_n);
    free(local_pts_elt_dist2);

    PDM_g_num_t *block_pts_elt_id = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           &part_stride,
                 // (void **) &part_elt_id,
                 (void **) &local_pts_elt_g_num,
                           &block_pts_elt_n,
                 (void **) &block_pts_elt_id);
    free(part_elt_id);
    free(part_stride);
    free(local_pts_elt_g_num);


    /* Pick closest elt for each point in current block */
    int block_n_pts = PDM_part_to_block_n_elt_block_get(ptb);
    PDM_g_num_t *block_pts_g_num = PDM_part_to_block_block_gnum_get(ptb);

    if (dbg_enabled) {
      int k = 0;
      log_trace("block_pts_elt_id (pre) :\n");
      for (int i = 0; i < block_n_pts; i++) {
        log_trace(" "PDM_FMT_G_NUM" [%d] -> ", block_pts_g_num[i], i);
        for (int j = 0; j < block_pts_elt_n[i]; j++) {
          log_trace("(%d, %f) ", block_pts_elt_id[k], block_pts_elt_dist2[k]);
          k++;
        }
        log_trace("\n");
      }
    }

    int idx = 0;
    for (int i = 0; i < block_n_pts; i++) {
      double min_dist2 = HUGE_VAL;
      for (int j = 0; j < block_pts_elt_n[i]; j++) {
        if (block_pts_elt_dist2[idx] < min_dist2) {
          min_dist2           = block_pts_elt_dist2[idx];
          block_pts_elt_id[i] = block_pts_elt_id   [idx];
        }
        idx++;
      }
    }

    if (dbg_enabled) {
      int *block_pts_elt_idx = PDM_array_new_idx_from_sizes_int(block_pts_elt_n,
                                                                block_n_pts);
      PDM_log_trace_connectivity_long(block_pts_elt_idx,
                                      block_pts_elt_id,
                                      block_n_pts,
                                      "block_pts_elt_id (post) : ");
      free(block_pts_elt_idx);
    }
    free(block_pts_elt_dist2);


    /* Send back to elements */
    PDM_g_num_t **tmp_part_elt_id = NULL;
    PDM_part_to_block_reverse_exch(ptb,
                                   sizeof(PDM_g_num_t),
                                   PDM_STRIDE_CST_INTERLACED,
                                   1,
                                   block_pts_elt_n,
                        (void   *) block_pts_elt_id,
                                   NULL,
                        (void ***) &tmp_part_elt_id);
    part_elt_id = tmp_part_elt_id[0];
    free(tmp_part_elt_id);
    free(block_pts_elt_n);
    free(block_pts_elt_id);
    PDM_part_to_block_free(ptb);

    // if (dbg_enabled) {
    //   PDM_log_trace_connectivity_long(delt_pts_idx2,
    //                                   part_elt_id,
    //                                   dn_elt2,
    //                                   "part_elt_id : ");
    // }

    /* Compress (get rid of false positives) */
    int         *final_elt_pts_n          = PDM_array_zeros_int(dn_elt2);
    PDM_g_num_t *final_elt_pts_g_num_geom = malloc(sizeof(PDM_g_num_t) * n_pts2                      );
    double      *final_elt_pts_coord      = malloc(sizeof(double     ) * n_pts2 * 3                  );
    double      *final_elt_pts_distance   = malloc(sizeof(double     ) * n_pts2                      );
    double      *final_elt_pts_proj_coord = malloc(sizeof(double     ) * n_pts2 * 3                  );
    int         *final_elt_pts_weight_idx = malloc(sizeof(int        ) * (n_pts2+1)                  );
    double      *final_elt_pts_weight     = malloc(sizeof(double     ) * delt_pts_weight_idx2[n_pts2]);
    double      *final_elt_pts_uvw        = malloc(sizeof(double     ) * n_pts2 * 3                  );

    final_elt_pts_weight_idx[0] = 0;
    idx = 0;
    for (int ielt = 0; ielt < dn_elt2; ielt++) {

      int idx_pt = delt_pts_idx2[ielt];
      int n_pts  = delt_pts_idx2[ielt+1] - idx_pt;
      // PDM_g_num_t *_elt_id       = part_elt_id            + idx_pt;
      PDM_g_num_t *_parent_g_num = delt_pts_g_num_geom    + idx_pt;
      double      *_dist2        = delt_pts_distance2     + idx_pt;
      double      *_coord        = delt_pts_coord2        + idx_pt*3;
      double      *_proj         = delt_pts_proj_coord2   + idx_pt*3;
      int         *_weight_idx   = delt_pts_weight_idx2   + idx_pt;
      double      *_uvw          = delt_pts_uvw2          + idx_pt*3;

      for (int i = 0; i < n_pts; i++) {

        // if (_elt_id[i] == delt_g_num_geom2[ielt]) {
        int pt_id = pts_unique_order[idx_pt+i];
        if (part_elt_id[pt_id] == delt_g_num_geom2[ielt]) {

          final_elt_pts_n[ielt]++;
          final_elt_pts_g_num_geom[idx] = _parent_g_num[i];
          final_elt_pts_distance  [idx] = _dist2[i];
          memcpy(final_elt_pts_coord      + 3*idx, _coord + 3*i, sizeof(double)*3);
          memcpy(final_elt_pts_proj_coord + 3*idx, _proj  + 3*i, sizeof(double)*3);

          int _weight_n = _weight_idx[i+1] - _weight_idx[i];
          double *_weight = delt_pts_weight2 + _weight_idx[i];
          final_elt_pts_weight_idx[idx+1] = final_elt_pts_weight_idx[idx] + _weight_n;
          memcpy(final_elt_pts_weight + final_elt_pts_weight_idx[idx],
                 _weight,
                 sizeof(double) * _weight_n);

          memcpy(final_elt_pts_uvw + 3*idx, _uvw + 3*i, sizeof(double)*3);

          idx++;

        }
      } // End of loop on current elt's pts
    } // End of loop on elts in frame 2
    free(delt_pts_idx2         );
    free(delt_pts_coord2       );
    free(delt_pts_distance2    );
    free(delt_pts_proj_coord2  );
    free(delt_pts_weight_idx2  );
    free(delt_pts_weight2      );
    free(delt_pts_uvw2         );
    free(part_elt_id           );
    free(delt_g_num_geom2);
    free(delt_pts_g_num_geom);

    free(pts_unique_order);
    free(pts_ln_to_gn);

    int final_n_pts = idx;
    final_elt_pts_g_num_geom = realloc(final_elt_pts_g_num_geom, sizeof(PDM_g_num_t) * final_n_pts);
    final_elt_pts_coord      = realloc(final_elt_pts_coord     , sizeof(double     ) * final_n_pts*3);
    final_elt_pts_distance   = realloc(final_elt_pts_distance  , sizeof(double     ) * final_n_pts);
    final_elt_pts_proj_coord = realloc(final_elt_pts_proj_coord, sizeof(double     ) * final_n_pts*3);
    final_elt_pts_weight_idx = realloc(final_elt_pts_weight_idx, sizeof(int        ) * (final_n_pts+1));
    final_elt_pts_weight     = realloc(final_elt_pts_weight    , sizeof(double     ) * final_elt_pts_weight_idx[final_n_pts]);
    final_elt_pts_uvw        = realloc(final_elt_pts_uvw       , sizeof(double     ) * final_n_pts*3);

    // PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[MERGE_LOCATION_DATA_STEP1] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [MERGE_LOCATION_DATA_STEP1] += e_t_cpu - b_t_cpu;
    ml->times_cpu_u  [MERGE_LOCATION_DATA_STEP1] += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s  [MERGE_LOCATION_DATA_STEP1] += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    dt_step1_all += PDM_MPI_Wtime() - beg_step1;
    double beg_step2 = PDM_MPI_Wtime();
    // PDM_MPI_Barrier (ml->comm);

    /*
     * Update results in user frame
     */
    PDM_block_to_part_t* btp_pts_gnum_geom_to_user = PDM_block_to_part_create(distrib_pts,
                                                      (const PDM_g_num_t **)  &final_elt_pts_g_num_geom,
                                                                              &final_n_pts,
                                                                              1,
                                                                              ml->comm);
    /*
     * Exchange gnum
     */
    PDM_g_num_t **tmp_final_elt_pts_g_num = NULL;
    PDM_block_to_part_exch(btp_pts_gnum_geom_to_user,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           &stride_one,
                           dpts_g_num_user,
                           NULL,
            (void ***)    &tmp_final_elt_pts_g_num);
    PDM_g_num_t *final_elt_pts_g_num = tmp_final_elt_pts_g_num[0];
    free(tmp_final_elt_pts_g_num);
    free(final_elt_pts_g_num_geom); // No longer used

    /*
     * Exchange init_location
     */
    int **tmp_final_elt_pts_triplet_n = NULL;
    int **tmp_final_elt_pts_triplet   = NULL;
    PDM_block_to_part_exch(btp_pts_gnum_geom_to_user,
                           3 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           dpts_init_location_pts_n,
                           dpts_init_location_pts,
                          &tmp_final_elt_pts_triplet_n,
            (void ***)    &tmp_final_elt_pts_triplet);
    int *final_elt_pts_triplet_n = tmp_final_elt_pts_triplet_n[0];
    int *final_elt_pts_triplet   = tmp_final_elt_pts_triplet  [0];
    free(tmp_final_elt_pts_triplet_n);
    free(tmp_final_elt_pts_triplet);

    free(dpts_init_location_pts_n);
    free(dpts_init_location_pts);

    idx_read = 0;
    int *final_elt_pts_triplet_idx = malloc((dn_elt2+1) * sizeof(int));
    final_elt_pts_triplet_idx[0] = 0;
    for(int i = 0; i < dn_elt2; ++i) {
      final_elt_pts_triplet_idx[i+1] = final_elt_pts_triplet_idx[i];
      int n_pts_in_elmt = final_elt_pts_n[i];
      for(int j = 0; j < n_pts_in_elmt; ++j) {
        final_elt_pts_triplet_idx[i+1] += final_elt_pts_triplet_n[idx_read++];
      }
    }
    assert(idx_read == final_n_pts);

    PDM_block_to_part_free(btp_pts_gnum_geom_to_user);

    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[MERGE_LOCATION_DATA_STEP2] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [MERGE_LOCATION_DATA_STEP2] += e_t_cpu - b_t_cpu;
    ml->times_cpu_u  [MERGE_LOCATION_DATA_STEP2] += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s  [MERGE_LOCATION_DATA_STEP2] += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);
    // PDM_MPI_Barrier (ml->comm);

    dt_step2_all += PDM_MPI_Wtime() - beg_step2;
    double beg_step3 = PDM_MPI_Wtime();

    /*
     *  Transfer location data from elt (current frame) to elt (user frame)
     */
    int   request_pts_in_elt_n   = -1;
    int   request_pts_g_num      = -1;
    int   request_pts_triplet    = -1;
    int   request_pts_triplet_n  = -1;
    int   request_pts_dist2      = -1;
          request_pts_coord      = -1;
    int   request_pts_proj_coord = -1;
    int   request_pts_uvw        = -1;
    int   request_pts_weight     = -1;
    int **pts_in_elt_n           = NULL;
    int **stride_pts             = NULL;
    int **stride_pts_triplet     = NULL;
    int **stride_pts_weight      = NULL;

    int *final_elt_pts_triplet_stride = NULL;
    int **pts_in_elt_triplet_n   = NULL;
    int **pts_in_elt_triplet     = NULL;

    _points_in_element_t *pts_in_elt = NULL;

    if (ml->reverse_result) {
      pts_in_elt = ml->points_in_elements + icloud;

      pts_in_elt->n_part = n_part;  // foireux?
      pts_in_elt->n_elts = malloc(sizeof(int) * n_part);
      memcpy(pts_in_elt->n_elts, pn_elt, sizeof(int) * n_part);

      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_CST_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             sizeof(int),
                             NULL,
            (const void  **) &final_elt_pts_n,
                             NULL,
            (      void ***) &pts_in_elt_n,
                             &request_pts_in_elt_n);
      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             sizeof(PDM_g_num_t),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_g_num,
                             &stride_pts,
            (      void ***) &pts_in_elt->gnum,
                             &request_pts_g_num);

      /* Echanger triplets des pts pour construire le ptp elt_user <-> pts_user */
      final_elt_pts_triplet_stride = malloc(sizeof(int) * dn_elt2);
      for (int i = 0; i < dn_elt2; i++) {
        final_elt_pts_triplet_stride[i] = final_elt_pts_triplet_idx[i+1] - final_elt_pts_triplet_idx[i];
      }

      // Exchange all triplets
      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             3*sizeof(int),
            (const int   **) &final_elt_pts_triplet_stride,
            (const void  **) &final_elt_pts_triplet,
                             &stride_pts_triplet,
            (      void ***) &pts_in_elt_triplet,
                             &request_pts_triplet);

      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_g_num);


      // Exchange number of triplet per point
      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             sizeof(int),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_triplet_n,
                             &stride_pts,
            (      void ***) &pts_in_elt_triplet_n,
                             &request_pts_triplet_n);

      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             sizeof(double),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_distance,
                             &stride_pts,
            (      void ***) &pts_in_elt->dist2,
                             &request_pts_dist2);

      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             3*sizeof(double),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_coord,
                             &stride_pts,
            (      void ***) &pts_in_elt->coords,
                             &request_pts_coord);

      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             3*sizeof(double),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_proj_coord,
                             &stride_pts,
            (      void ***) &pts_in_elt->projected_coords,
                             &request_pts_proj_coord);

      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             3*sizeof(double),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_uvw,
                             &stride_pts,
            (      void ***) &pts_in_elt->uvw,
                             &request_pts_uvw);
    }
    // free(final_elt_pts_triplet_n);

    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[REVERSE_LOCATION_DATA_BTP1] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [REVERSE_LOCATION_DATA_BTP1] += e_t_cpu - b_t_cpu;
    ml->times_cpu_u  [REVERSE_LOCATION_DATA_BTP1] += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s  [REVERSE_LOCATION_DATA_BTP1] += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    // end_timer_and_log("Step 2 : ", ml->comm, t1_step2);
    // PDM_MPI_Barrier (ml->comm);

    dt_step3_all += PDM_MPI_Wtime() - beg_step3;
    double beg_step4 = PDM_MPI_Wtime();


    int *final_elt_pts_idx = PDM_array_new_idx_from_sizes_int(final_elt_pts_n, dn_elt2);
    // free(final_elt_pts_n);

    int *elt_pts_weight_stride = NULL;
    if (ml->reverse_result) {
      elt_pts_weight_stride = malloc(sizeof(int) * dn_elt2);
      for (int ielt = 0; ielt < dn_elt2; ielt++) {
        int head = final_elt_pts_idx[ielt  ];
        int tail = final_elt_pts_idx[ielt+1];

        elt_pts_weight_stride[ielt] = final_elt_pts_weight_idx[tail] - final_elt_pts_weight_idx[head];
      }
      if (dbg_enabled) {
        PDM_log_trace_array_int(elt_pts_weight_stride, dn_elt2, "elt_pts_weight_stride : ");
      }
    }

    if (ml->reverse_result) {
      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             sizeof(double),
            (const int   **) &elt_pts_weight_stride,
            (const void  **) &final_elt_pts_weight,
                             &stride_pts_weight,
            (      void ***) &pts_in_elt->weights,
                             &request_pts_weight);


      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_in_elt_n);


      int  *n_ref_elt = NULL;
      int **ref_elt   = NULL;
      PDM_part_to_part_ref_lnum2_get(ptp_elt,
                                     &n_ref_elt,
                                     &ref_elt);

      int  *n_unref_elt = NULL;
      int **unref_elt   = NULL;
      PDM_part_to_part_unref_lnum2_get(ptp_elt,
                                       &n_unref_elt,
                                       &unref_elt);

      pts_in_elt->pts_inside_idx = malloc(sizeof(int *) * n_part);
      for (int ipart = 0; ipart < n_part; ipart++) {
        pts_in_elt->pts_inside_idx[ipart] = malloc(sizeof(int) * (pn_elt[ipart]+1));
        pts_in_elt->pts_inside_idx[ipart][0] = 0;

        for (int i = 0; i < n_ref_elt[ipart]; i++) {
          pts_in_elt->pts_inside_idx[ipart][ref_elt[ipart][i]] = pts_in_elt_n[ipart][i];
        }
        free(pts_in_elt_n[ipart]);

        for (int i = 0; i < n_unref_elt[ipart]; i++) {
          pts_in_elt->pts_inside_idx[ipart][unref_elt[ipart][i]] = 0;
        }

        for (int i = 0; i < pn_elt[ipart]; i++) {
          pts_in_elt->pts_inside_idx[ipart][i+1] += pts_in_elt->pts_inside_idx[ipart][i];
        }
      }
      free(pts_in_elt_n);

      pts_in_elt->weights_idx = malloc(sizeof(int *) * n_part);
      for (int ipart = 0; ipart < n_part; ipart++) {
        pts_in_elt->weights_idx[ipart] = malloc(sizeof(int) * (pts_in_elt->pts_inside_idx[ipart][pn_elt[ipart]] + 1));
        pts_in_elt->weights_idx[ipart][0] = 0;

        for (int ielt = 0; ielt < pn_elt[ipart]; ielt++) {
          int n_vtx = ml->cell_vtx_idx[ipart][ielt+1] - ml->cell_vtx_idx[ipart][ielt];
          for (int idx_pt = pts_in_elt->pts_inside_idx[ipart][ielt]; idx_pt < pts_in_elt->pts_inside_idx[ipart][ielt+1]; idx_pt++) {
            pts_in_elt->weights_idx[ipart][idx_pt+1] = pts_in_elt->weights_idx[ipart][idx_pt] + n_vtx;
          }
        }
      }

      if (!full_async) {
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_dist2);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_coord);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_proj_coord);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_weight);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_uvw);
        // for (int ipart = 0; ipart < n_part; ipart++) {
        //   free(stride_pts_weight[ipart]);
        // }
        // free(stride_pts_weight);
      }

      free(elt_pts_weight_stride);
    }

    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[REVERSE_LOCATION_DATA_BTP2] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [REVERSE_LOCATION_DATA_BTP2] += e_t_cpu - b_t_cpu;
    ml->times_cpu_u  [REVERSE_LOCATION_DATA_BTP2] += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s  [REVERSE_LOCATION_DATA_BTP2] += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    // end_timer_and_log("Step 3 : ", ml->comm, t1_step3);
    // PDM_MPI_Barrier (ml->comm);
    // double t1_step4 = PDM_MPI_Wtime();

    dt_step4_all += PDM_MPI_Wtime() - beg_step4;
    double beg_step5 = PDM_MPI_Wtime();


    /*
     *  Create ptp to exchange data between elt and points (both in user frame)
     */
    if (ml->reverse_result) {

      double t1_wait = PDM_MPI_Wtime();
      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_triplet);
      free(final_elt_pts_triplet_stride);
      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_triplet_n);
      free(final_elt_pts_triplet_n);
      if (log_timer) {
        log_trace("wait = %12.5e \n", PDM_MPI_Wtime() - t1_wait);
      }

      int  *n_ref_elt = NULL;
      int **ref_elt   = NULL;
      PDM_part_to_part_ref_lnum2_get(ptp_elt,
                                     &n_ref_elt,
                                     &ref_elt);

      int  *n_unref_elt = NULL;
      int **unref_elt   = NULL;
      PDM_part_to_part_unref_lnum2_get(ptp_elt,
                                       &n_unref_elt,
                                       &unref_elt);

      int **pts_in_elt_triplet_idx = malloc(sizeof(int *) * n_part);
      for (int ipart = 0; ipart < n_part; ipart++) {
        free(stride_pts_triplet[ipart]);

        int _n_pts = pts_in_elt->pts_inside_idx[ipart][pn_elt[ipart]];

        // PDM_log_trace_array_int(pts_in_elt_triplet_n[ipart],
        //                         _n_pts,
        //                         "pts_in_elt_triplet_n : ");

        pts_in_elt_triplet_idx[ipart] = malloc(sizeof(int) * (_n_pts + 1));
        pts_in_elt_triplet_idx[ipart][0] = 0;
        int max_triplet_n = 0;
        for (int i = 0; i < _n_pts; i++) {
          max_triplet_n = PDM_MAX(max_triplet_n, pts_in_elt_triplet_n[ipart][i]);

          pts_in_elt_triplet_idx[ipart][i+1] = pts_in_elt_triplet_idx[ipart][i] +
          pts_in_elt_triplet_n[ipart][i]*3;
        }

        if (0) {
          int *_idx = malloc(sizeof(int) * (_n_pts + 1));
          for (int i = 0; i <= _n_pts; i++) {
            _idx[i] = pts_in_elt_triplet_idx[ipart][i]/3;
          }

          PDM_log_trace_array_long(pts_in_elt->gnum[ipart],
                                   _n_pts,
                                   "pts_in_elt_gnum : ");

          PDM_log_trace_graph_nuplet_int(_idx,
                                         pts_in_elt_triplet[ipart],
                                         3,
                                         _n_pts,
                                         "pts_in_elt_triplet      : ");
          free(_idx);
        }



        /* Lexicographic sort each point's triplets */
        double t1_sort = PDM_MPI_Wtime();
        int *order = malloc(sizeof(int) * max_triplet_n);
        for (int i = 0; i < _n_pts; i++) {
          PDM_order_lnum_s(pts_in_elt_triplet[ipart] + pts_in_elt_triplet_idx[ipart][i],
                           3,
                           order,
                           pts_in_elt_triplet_n[ipart][i]);
        }
        free(order);
        if (log_timer) {
          log_trace("sort = %12.5e \n", PDM_MPI_Wtime() - t1_sort);
        }


        free(pts_in_elt_triplet_n[ipart]);

      }
      free(stride_pts_triplet);
      free(pts_in_elt_triplet_n);
      free(final_elt_pts_triplet);
      free(final_elt_pts_triplet_idx);
      free(final_elt_pts_idx);
      free(final_elt_pts_g_num);

      // PDM_part_to_part_free(ptp_elt);

      dt_step5_all += PDM_MPI_Wtime() - beg_step5;
      double beg_step6 = PDM_MPI_Wtime();

      // log_trace("PDM_part_to_part_create_from_num2_triplet start \n");
      ml->ptp[icloud] = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) elt_g_num,
                                                                  (const int          *) pn_elt,
                                                                  n_part,
                                                                  (const int          *) pcloud->n_points,
                                                                  pcloud->n_part,
                                                                  (const int         **) pts_in_elt->pts_inside_idx, // size = n_elt
                                                                  (const int         **) pts_in_elt_triplet_idx,     // size = pts_inside_idx[n_elt]
                                                                  (const int         **) pts_in_elt_triplet,
                                                                  ml->comm);
      // log_trace("PDM_part_to_part_create_from_num2_triplet end \n");

      for (int ipart = 0; ipart < n_part; ipart++) {
        free(pts_in_elt_triplet    [ipart]);
        free(pts_in_elt_triplet_idx[ipart]);
      }
      free(pts_in_elt_triplet);
      free(pts_in_elt_triplet_idx);

      dt_step6_all += PDM_MPI_Wtime() - beg_step6;
    }

    // for (int ipart = 0; ipart < n_part; ipart++) {
    //   free(stride_pts[ipart]);
    // }
    // free(stride_pts);

    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[REVERSE_LOCATION_DATA_BTP3] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [REVERSE_LOCATION_DATA_BTP3] += e_t_cpu - b_t_cpu;
    ml->times_cpu_u  [REVERSE_LOCATION_DATA_BTP3] += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s  [REVERSE_LOCATION_DATA_BTP3] += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    // end_timer_and_log("Step 4 : ", ml->comm, t1_step4);
    // PDM_MPI_Barrier (ml->comm);

    // TODO: ownership on located/unlocated??
    int  *n_located = NULL;
    int **located   = NULL;
    PDM_part_to_part_ref_lnum2_get(ml->ptp[icloud],
                                   &n_located,
                                   &located);

    int  *n_unlocated = NULL;
    int **unlocated   = NULL;
    PDM_part_to_part_unref_lnum2_get(ml->ptp[icloud],
                                     &n_unlocated,
                                     &unlocated);

    pcloud->n_located    = malloc(sizeof(int  ) * pcloud->n_part);
    pcloud->n_un_located = malloc(sizeof(int  ) * pcloud->n_part);
    pcloud->located      = malloc(sizeof(int *) * pcloud->n_part);
    pcloud->un_located   = malloc(sizeof(int *) * pcloud->n_part);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      pcloud->n_located   [ipart] = n_located  [ipart];
      pcloud->n_un_located[ipart] = n_unlocated[ipart];

      pcloud->located   [ipart] = malloc(sizeof(int) * n_located  [ipart]);
      pcloud->un_located[ipart] = malloc(sizeof(int) * n_unlocated[ipart]);

      memcpy(pcloud->located   [ipart], located  [ipart], sizeof(int) * n_located  [ipart]);
      memcpy(pcloud->un_located[ipart], unlocated[ipart], sizeof(int) * n_unlocated[ipart]);
    }

    // TODO: ownership on gnum1_come_from(_idx)??
    int         **gnum1_come_from_idx = NULL;
    PDM_g_num_t **gnum1_come_from     = NULL;
    PDM_part_to_part_gnum1_come_from_get(ml->ptp[icloud],
                                         &gnum1_come_from_idx,
                                         &gnum1_come_from);

    pcloud->location = malloc(sizeof(PDM_g_num_t *) * pcloud->n_part);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      pcloud->location[ipart] = malloc(sizeof(PDM_g_num_t) * n_located[ipart]);
      for (int i = 0; i < n_located[ipart]; i++) {

        if (gnum1_come_from_idx[ipart][i+1] != gnum1_come_from_idx[ipart][i] + 1) {
          int ipt = pcloud->located[ipart][i] - 1;
          log_trace("point "PDM_FMT_G_NUM" has locations ",
                    pcloud->gnum[ipart][ipt]);
          PDM_log_trace_array_long(gnum1_come_from[ipart] + gnum1_come_from_idx[ipart][i],
                                   gnum1_come_from_idx[ipart][i+1] - gnum1_come_from_idx[ipart][i],
                                   "");
        }
        assert(gnum1_come_from_idx[ipart][i+1] == gnum1_come_from_idx[ipart][i] + 1);

        pcloud->location[ipart][i] = gnum1_come_from[ipart][i];
      }

      if (dbg_enabled) {
        PDM_log_trace_array_long(pcloud->location[ipart], n_located[ipart], "pcloud->location[ipart] : ");
      }
    }


    if (dbg_enabled) {
      double **is_located = malloc(sizeof(double) * pcloud->n_part);
      double **location   = malloc(sizeof(double) * pcloud->n_part);
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        is_located[ipart] = malloc(sizeof(double) * pcloud->n_points[ipart]);
        location  [ipart] = malloc(sizeof(double) * pcloud->n_points[ipart]);
        for (int i = 0; i < pcloud->n_points[ipart]; i++) {
          is_located[ipart][i] = -1;
        }
        for (int i = 0; i < n_located[ipart]; i++) {
          is_located[ipart][located[ipart][i]-1] = 1;
          location  [ipart][located[ipart][i]-1] = (double) pcloud->location[ipart][i];
        }
        for (int i = 0; i < n_unlocated[ipart]; i++) {
          is_located[ipart][unlocated[ipart][i]-1] = 0;
          location  [ipart][unlocated[ipart][i]-1] = -1;
        }
      }

      const char  *field_name[]   = {"is_located", "location"};
      double     **field_value[2] = {is_located, location};

      char name[999];
      sprintf(name, "mesh_location_point_cloud_%d_loc", icloud);
      _dump_point_cloud(name,
                        ml->comm,
                        pcloud->n_part,
                        pcloud->n_points,
                        pcloud->coords,
                        pcloud->gnum,
                        2,
                        field_name,
                        field_value);

      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        free(is_located[ipart]);
        free(location  [ipart]);
      }
      free(is_located);
      free(location);
    }

    if (full_async) {
      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_proj_coord);
    }
    /* Exchange other fields (dist2, uvw, weights(_idx), proj_coord) */
    /* Really necessary from the points' PoV ?? */
    int **tmp_projected_coords_n = NULL;
    PDM_part_to_part_iexch(ml->ptp[icloud],
                           comm_kind,
                           PDM_STRIDE_CST_INTERLACED,
                           PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                           1,
                           3*sizeof(double),
                           NULL,
          (const void  **) pts_in_elt->projected_coords,
                           &tmp_projected_coords_n,
          (      void ***) &pcloud->projected_coords,
                           &req_pts_proj_coord[icloud]);

    // int request_pts_dist2;
    if (full_async) {
      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_dist2);
    }

    int **tmp_pts_dist2_n = NULL;
    PDM_part_to_part_iexch(ml->ptp[icloud],
                           comm_kind,
                           PDM_STRIDE_CST_INTERLACED,
                           PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                           1,
                           sizeof(double),
                           NULL,
          (const void  **) pts_in_elt->dist2,
                           &tmp_pts_dist2_n,
          (      void ***) &pcloud->dist2,
                           &req_pts_dist2[icloud]);

    // PDM_part_to_part_iexch_wait(ml->ptp[icloud], request_pts_proj_coord);
    // PDM_part_to_part_iexch_wait(ml->ptp[icloud], request_pts_dist2);

    if (ptp_elt != NULL) {
      if (full_async) {
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_uvw);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_weight);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_coord);
      }
      for (int ipart = 0; ipart < n_part; ipart++) {
        free(stride_pts_weight[ipart]);
        free(stride_pts       [ipart]);
      }
      free(stride_pts_weight);
      free(stride_pts       );
      PDM_part_to_part_free(ptp_elt);
    }


    // PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[REVERSE_LOCATION_DATA_UVW] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [REVERSE_LOCATION_DATA_UVW] += e_t_cpu - b_t_cpu;
    ml->times_cpu_u  [REVERSE_LOCATION_DATA_UVW] += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s  [REVERSE_LOCATION_DATA_UVW] += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    /* Free memory */

    /* See if stuff can be freed earlier.... */
    free(delt_parent_g_num2      );
    free(final_elt_pts_coord     );
    free(final_elt_pts_distance  );
    free(final_elt_pts_proj_coord);
    free(final_elt_pts_weight_idx);
    free(final_elt_pts_weight    );
    free(final_elt_pts_uvw       );
    free(final_elt_pts_n         );


    if (use_extracted_pts) {
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        free(select_pts_g_num_user[ipart]);
        free(select_pts_coord     [ipart]);
      }
      free(select_pts_g_num_user);
      free(select_pts_coord);
      free(n_select_pts);
    }
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      free(select_pts_init_location[ipart]);
    }
    free(select_pts_init_location);


    if (use_extracted_mesh) {
      for (int ipart = 0; ipart < n_part; ipart++) {
        free(select_elt_extents   [ipart]);
        free(select_elt_g_num_user[ipart]);
      }
      free(select_elt_extents);
      free(select_elt_g_num_user);
      free(n_select_elt);
    }
    for (int ipart = 0; ipart < n_part; ipart++) {
      free(select_elt_l_num[ipart]);
      free(select_elt_init_location_user[ipart]);
    }
    free(select_elt_l_num);
    free(select_elt_init_location_user);

    PDM_part_to_block_free(ptb_pts);
    PDM_part_to_block_free(ptb_elt);

  } /* End icloud */

  _pmesh_nodal_elmts_free(pmne,
                          (ml->mesh_nodal == NULL));

  for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {
    if (ml->ptp[icloud] != NULL) {
      if (req_pts_proj_coord[icloud] != -1) {
        // log_trace("out wait req_pts_proj_coord = %i \n", req_pts_proj_coord[icloud]);
        PDM_part_to_part_iexch_wait(ml->ptp[icloud], req_pts_proj_coord[icloud]);
        // log_trace("out wait req_pts_proj_coord END = %i \n", req_pts_proj_coord[icloud]);
      }

      if (req_pts_dist2[icloud] != -1) {
        // log_trace("out wait req_pts_dist2 = %i \n", req_pts_dist2[icloud]);
        PDM_part_to_part_iexch_wait(ml->ptp[icloud], req_pts_dist2[icloud]);
        // log_trace("out wait req_pts_dist2 END = %i \n", req_pts_dist2[icloud]);
      }
    }
  }
  free(req_pts_proj_coord);
  free(req_pts_dist2);

  if (log_timer) {
    log_trace("dt_extract_all = %12.5e / dt_search_all = %12.5e \n", dt_extract_all, dt_search_all);
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(elt_extents[ipart]);
    free(elt_g_num  [ipart]);
    free(elt_type   [ipart]);
  }
  free(elt_extents);
  free(elt_g_num);
  free(elt_type);
  free(pn_elt);

  PDM_timer_hang_on(ml->timer);
  ml->times_elapsed[END] = PDM_timer_elapsed (ml->timer);
  ml->times_cpu    [END] = PDM_timer_cpu     (ml->timer);
  ml->times_cpu_u  [END] = PDM_timer_cpu_user(ml->timer);
  ml->times_cpu_s  [END] = PDM_timer_cpu_sys (ml->timer);
  PDM_timer_resume(ml->timer);

  if (log_timer) {
    end_timer_and_log_from_dt("Step 1 : ", ml->comm, dt_step1_all);
    end_timer_and_log_from_dt("Step 2 : ", ml->comm, dt_step2_all);
    end_timer_and_log_from_dt("Step 3 : ", ml->comm, dt_step3_all);
    end_timer_and_log_from_dt("Step 4 : ", ml->comm, dt_step4_all);
    end_timer_and_log_from_dt("Step 5 : ", ml->comm, dt_step5_all);
    end_timer_and_log_from_dt("Step 6 : ", ml->comm, dt_step6_all);
  }

  // exit(1);
}

/**
 *
 * \brief Get the number of cells
 *
 * \param [in]  id       Pointer to \ref PDM_mesh_location object
 * \param [in]  i_part   Index of partition of the mesh
 *
 * \return Number of cells
 */

int
PDM_mesh_location_n_cell_get
(
       PDM_mesh_location_t *ml,
 const int                  i_part
)
{
  return ml->points_in_elements[0].n_elts[i_part];
}


PDM_Mesh_nodal_t*
PDM_mesh_location_mesh_nodal_get
(
 PDM_mesh_location_t *ml
)
{
  return ml->mesh_nodal;
}


/**
 * Enable reverse results computation (To call PDM_mesh_location_points_in_elt_get)
 */

void
PDM_mesh_location_reverse_results_enable
(
PDM_mesh_location_t *ml
)
{
  ml->reverse_result = 1;
}



/**
 * \brief Get part_to_part object to exchange data between
 * the source mesh and a target point cloud (both in user frame)
 *
 * \param [in ] ml         Pointer to \ref PDM_mesh_location_t object
 * \param [in ] icloud     Point cloud ID
 * \param [out] ptp        Pointer to \ref PDM_part_to_part_t object
 * \param [in ] ownership  Ownership for ptp
 *
 */

void
PDM_mesh_location_part_to_part_get
(
       PDM_mesh_location_t  *ml,
 const int                   icloud,
       PDM_part_to_part_t  **ptp,
       PDM_ownership_t       ownership
 )
{
  assert(icloud < ml->n_point_cloud);

  *ptp = ml->ptp[icloud];
  ml->ptp_ownership[icloud] = ownership;
}


#ifdef	__cplusplus
}
#endif
