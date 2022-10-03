/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi.h"
#include "pdm_octree.h"
#include "pdm_octree_seq.h"
#include "pdm_array.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_points_merge.h"
#include "pdm_points_merge_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/


/*============================================================================
 * Global variable
 *============================================================================*/

static const double _default_eps = 1e-9;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Search a point
 *
 * \param [in]   ppart_id        ppart identifier
 *
 */

static int
_intersect_extents
(
const double *first_extents,
const double *second_extents
)
{
  int intersect = 0;

  for (int i = 0; i < 3; i++) {
    if ((first_extents[i] >= second_extents[i]) &&
        (first_extents[i] <= second_extents[i+3])) {
      intersect += 1;
    }
    else if ((first_extents[i+3] >= second_extents[i]) &&
            (first_extents[i+3] <= second_extents[i+3])) {
      intersect += 1;
    }
    else if ((first_extents[i] <= second_extents[i]) &&
             (first_extents[i+3] >= second_extents[i+3])) {
      intersect += 1;
    }
  }

  return (intersect == 3);
}


/**
 *
 * \brief Search a point in local partitions
 *
 * \param [in]   ppart_id        ppart identifier
 *
 */

static void
_search_local_couple
(
int **local_couple,
int  *n_couple,
int  *s_couple,
const int  point_cloud,
const int  point_idx,
const double *point_coords,
const double *point_box,
const int search_cloud,
PDM_octree_seq_t *octree,
const int associated_octree_node_id,
const double *associated_coords,
const double *associated_char_length,
const double tolerance
)
{
  int node_id = associated_octree_node_id;
  const double *coords = associated_coords;
  const double *char_length = associated_char_length;

  if (PDM_octree_seq_leaf_is (octree, node_id)) {

    int *points_clouds_id;
    int *point_indexes;
    int n_candidates = PDM_octree_seq_n_points_get (octree, node_id);
    PDM_octree_seq_points_get (octree, node_id,
                               &points_clouds_id, &point_indexes);

    for (int i = 0; i < n_candidates; i++) {
      double dist2;
      double tol;

      if (point_box != NULL) {
        const double *coords_candidate = associated_coords + 3 * point_indexes[i];
        double char_length_candidate = associated_char_length[point_indexes[i]];
        dist2 = (coords_candidate[0] - point_coords[0]) *
                (coords_candidate[0] - point_coords[0]) +
                (coords_candidate[1] - point_coords[1]) *
                (coords_candidate[1] - point_coords[1]) +
                (coords_candidate[2] - point_coords[2]) *
                (coords_candidate[2] - point_coords[2]);
        double tol1 = char_length_candidate * tolerance *
                      char_length_candidate * tolerance;

        double tol2 = (point_box[0] - point_coords[0]) * (point_box[0] - point_coords[0]);

        tol = PDM_MIN (tol1, tol2);

      }
      else {
        const double *coords_candidate = associated_coords + 3 * point_indexes[i];
        dist2 = (coords_candidate[0] - point_coords[0]) *
                (coords_candidate[0] - point_coords[0]) +
                (coords_candidate[1] - point_coords[1]) *
                (coords_candidate[1] - point_coords[1]) +
                (coords_candidate[2] - point_coords[2]) *
                (coords_candidate[2] - point_coords[2]);
        tol = _default_eps * tolerance *
              _default_eps * tolerance;
      }

      if (dist2 <= tol) {
        if (*n_couple >= *s_couple) {
          if (*s_couple == 0) {
            *s_couple = 4;
          }
          else {
            *s_couple *= 2;
          }
          *local_couple = realloc(*local_couple, sizeof(int) * (*s_couple) * 4);
        }

        int _n_couple = *n_couple;
        int *_local_couple = *local_couple;

        if (!((search_cloud == point_cloud) && (point_indexes[i] == point_idx))) {
          _local_couple[4*_n_couple]   = search_cloud;
          _local_couple[4*_n_couple+1] = point_indexes[i];
          _local_couple[4*_n_couple+2] = point_cloud;
          _local_couple[4*_n_couple+3] = point_idx;
          *n_couple += 1;
        }
      }
    }

  }

  else {
    for (int i = 0; i < 8; i++) {
      const int node_child =
            PDM_octree_seq_children_get (octree, node_id,
                                         (PDM_octree_seq_child_t) i);
      if (node_child != -1) {
        if (point_box != NULL) {
          if (_intersect_extents (PDM_octree_seq_node_extents_get (octree, node_child),
                                  point_box)) {

            _search_local_couple (local_couple, n_couple, s_couple, point_cloud,
                                  point_idx, point_coords, point_box, search_cloud,
                                  octree, node_child, coords, char_length, tolerance);
          }
        }
        else {
          double _extents[6] = {point_coords[0] - _default_eps,
                                point_coords[1] - _default_eps,
                                point_coords[2] - _default_eps,
                                point_coords[0] + _default_eps,
                                point_coords[1] + _default_eps,
                                point_coords[2] + _default_eps};

          if (_intersect_extents (PDM_octree_seq_node_extents_get (octree, node_child),
                                  _extents)) {

            _search_local_couple (local_couple, n_couple, s_couple, point_cloud,
                                  point_idx, point_coords, point_box, search_cloud,
                                  octree, node_child, coords, char_length, tolerance);
          }
        }
      }
    }
  }
}


/**
 *
 * \brief Search a point in distant partitions
 *
 * \param [in]   ppart_id        ppart identifier
 *
 */

static void
_search_distant_couple
(
int           *n_fusion_from_proc,
int          **distant_couple,
int           *n_couple,
int           *s_couple,
const int      point_proc,
const int      point_cloud,
const int      point_idx,
const double  *point_coords,
const double  *point_box,
PDM_octree_t *octree,
const int      associated_octree_node_id,
const double **associated_coords,
const double **associated_char_length,
const double   tolerance
)
{

  int node_id = associated_octree_node_id;

  assert(node_id != -1);

  if (PDM_octree_leaf_is (octree, node_id)) {

    int *points_clouds_id;
    int *point_indexes;
    int n_candidates = PDM_octree_n_points_get (octree, node_id);
    PDM_octree_points_get (octree, node_id,
                               &points_clouds_id, &point_indexes);

    for (int i = 0; i < n_candidates; i++) {
      double dist2;
      double tol;

      if (point_box != NULL) {
        const double *coords_candidate =
                      associated_coords[points_clouds_id[i]] + 3 * point_indexes[i];
        const double char_length_candidate =
                      associated_char_length[points_clouds_id[i]][point_indexes[i]];
        dist2 = (coords_candidate[0] - point_coords[0]) *
                (coords_candidate[0] - point_coords[0]) +
                (coords_candidate[1] - point_coords[1]) *
                (coords_candidate[1] - point_coords[1]) +
                (coords_candidate[2] - point_coords[2]) *
                (coords_candidate[2] - point_coords[2]);
        double tol1 = char_length_candidate * tolerance *
                      char_length_candidate * tolerance;

        double tol2 = (point_box[0] - point_coords[0]) * (point_box[0] - point_coords[0]);

        tol = PDM_MIN (tol1, tol2);
        if (0 == 1) {
          if (dist2 <= tol) {
            printf ("point local : %d\n", point_indexes[i]);
            printf ("point distant : %d %d\n", point_proc, point_idx) ;
            printf("dist2 tol : %12.5e %12.5e\n", dist2, tol);
            if (dist2 <= tol) printf("**** Fusion *****\n");
            printf("%12.5e %12.5e\n", coords_candidate[0], point_coords[0]);
            printf("%12.5e %12.5e\n", coords_candidate[1], point_coords[1]);
            printf("%12.5e %12.5e\n", coords_candidate[2], point_coords[2]);
          }
        }
      }
      else {
        const double *coords_candidate =
                      associated_coords[points_clouds_id[i]] + 3 * point_indexes[i];
        dist2 = (coords_candidate[0] - point_coords[0]) *
                (coords_candidate[0] - point_coords[0]) +
                (coords_candidate[1] - point_coords[1]) *
                (coords_candidate[1] - point_coords[1]) +
                (coords_candidate[2] - point_coords[2]) *
                (coords_candidate[2] - point_coords[2]);
        tol = _default_eps * tolerance *
              _default_eps * tolerance;
        if (1 == 0) {
          if (dist2 <= tol) {
            printf("dist2 tol : %12.5e %12.5e\n", dist2, tol);
            printf("%12.5e %12.5e\n", coords_candidate[0], point_coords[0]);
            printf("%12.5e %12.5e\n", coords_candidate[1], point_coords[1]);
            printf("%12.5e %12.5e\n", coords_candidate[2], point_coords[2]);
          }
        }
      }

      if (dist2 <= tol) {
        if (*n_couple >= *s_couple) {
          if (*s_couple == 0) {
            *s_couple = 4;
          }
          else {
            *s_couple *= 2;
          }
          *distant_couple = realloc(*distant_couple, sizeof(int) * (*s_couple) * 5);
        }

        int _n_couple = *n_couple;
        int *_distant_couple = *distant_couple;

        _distant_couple[5*_n_couple  ] = points_clouds_id[i];
        _distant_couple[5*_n_couple+1] = point_indexes[i];
        _distant_couple[5*_n_couple+2] = point_proc;
        _distant_couple[5*_n_couple+3] = point_cloud;
        _distant_couple[5*_n_couple+4] = point_idx;
        *n_couple += 1;
        n_fusion_from_proc[point_proc] += 1;

      }
    }
  }

  else {
    for (int i = 0; i < 8; i++) {
      const int node_child =
            PDM_octree_children_get (octree, node_id,
                                     (PDM_octree_child_t) i);
      if (node_child != -1) {
        if (point_box != NULL) {
          if (_intersect_extents (PDM_octree_node_extents_get (octree, node_child),
                                  point_box)) {

            _search_distant_couple (n_fusion_from_proc,
                                    distant_couple, n_couple, s_couple,
                                    point_proc, point_cloud, point_idx,
                                    point_coords, point_box,
                                    octree,
                                    node_child,
                                    associated_coords,
                                    associated_char_length,
                                    tolerance);

          }
        }
        else {
          double _extents[6] = {point_coords[0] - _default_eps,
                                point_coords[1] - _default_eps,
                                point_coords[2] - _default_eps,
                                point_coords[0] + _default_eps,
                                point_coords[1] + _default_eps,
                                point_coords[2] + _default_eps};

          if (_intersect_extents (PDM_octree_node_extents_get (octree, node_child),
                                  _extents)) {

            _search_distant_couple (n_fusion_from_proc,
                                    distant_couple, n_couple, s_couple,
                                    point_proc, point_cloud, point_idx,
                                    point_coords, point_box,
                                    octree,
                                    node_child,
                                    associated_coords,
                                    associated_char_length,
                                    tolerance);
          }
        }
      }
    }
  }
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a points merge structure
 *
 * \param [in]   n_point_cloud      Number of point cloud
 * \param [in]   tolerance          Relative geometric tolerance
 * \param [in]   comm               MPI communicator
 * \param [in]   owner              Ownership
 *
 * \return     Pointer to \ref PDM_points_merge object
 */

PDM_points_merge_t *
PDM_points_merge_create
(
 const int             n_point_cloud,
 const double          tolerance,
 const PDM_MPI_Comm    comm,
 const PDM_ownership_t owner
)
{
  PDM_points_merge_t *pm = (PDM_points_merge_t *) malloc(sizeof(PDM_points_merge_t));

  pm->comm              = comm;
  pm->owner             = owner;
  pm->results_is_getted = PDM_FALSE;
  pm->tolerance         = tolerance;
  pm->n_point_clouds    = n_point_cloud;
  pm->n_points          = malloc (sizeof(int     ) * n_point_cloud);
  pm->point_clouds      = malloc (sizeof(double *) * n_point_cloud);
  pm->char_length       = malloc (sizeof(double *) * n_point_cloud);
  pm->octree            = NULL;
  pm->candidates_idx    = malloc (sizeof(int *) * n_point_cloud);
  pm->candidates_desc   = malloc (sizeof(int *) * n_point_cloud);

  for (int i = 0; i < n_point_cloud; i++) {
    pm->candidates_idx [i] = NULL;
    pm->point_clouds   [i] = NULL;
    pm->char_length    [i] = NULL;
    pm->candidates_desc[i] = NULL;
  }

  pm->depth_max = 31;
  pm->points_in_leaf_max = 4;

  pm->octree = PDM_octree_create (n_point_cloud, pm->depth_max,
                                  pm->points_in_leaf_max, tolerance, comm);

  return pm;

}


/**
 *
 * \brief Free a \ref PDM_points_merge object
 *
 * \param [in]   pm             Pointer to \ref PDM_points_merge object
 *
 */

void
PDM_points_merge_free
(
 PDM_points_merge_t *pm
)
{

  if(( pm->owner == PDM_OWNERSHIP_KEEP ) ||
     ( pm->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !pm->results_is_getted)){
    for (int i = 0; i < pm->n_point_clouds; i++) {
      if (pm->candidates_idx[i] != NULL) {
        free (pm->candidates_idx[i]);
      }
      if (pm->candidates_desc[i] != NULL) {
        free (pm->candidates_desc[i]);
      }
    }
  }

  free (pm->candidates_idx);
  free (pm->candidates_desc);
  free (pm->point_clouds);
  free (pm->char_length);
  free (pm->n_points);

  free (pm);
}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   pm             Pointer to \ref PDM_points_merge object
 * \param [in]   i_point_cloud  Index of point cloud
 * \param [in]   n_points       Number of points
 * \param [in]   coords         Point coordinates
 * \param [in]   char_length    Characteristic length (or NULL)
 *
 */

void
PDM_points_merge_cloud_set
(
       PDM_points_merge_t *pm,
 const int                 i_point_cloud,
 const int                 n_points,
 const double             *coords,
 const double             *char_length
)
{

  pm->char_length[i_point_cloud] = char_length;
  pm->point_clouds[i_point_cloud] = coords;
  pm->n_points[i_point_cloud] = n_points;

  PDM_octree_point_cloud_set (pm->octree, i_point_cloud, n_points, coords, NULL);

}


/**
 *
 * \brief Process merge points
 *
 * \param [in]   pm             Pointer to \ref PDM_points_merge object
 *
 */

void
PDM_points_merge_process
(
 PDM_points_merge_t *pm
)
{
  PDM_octree_build (pm->octree);

  int n_rank;
  PDM_MPI_Comm_size(pm->comm , &n_rank);

  int i_rank;
  PDM_MPI_Comm_rank(pm->comm , &i_rank);

  int *local_couple = NULL;
  int n_local_couple = 0;
  int s_local_couple = 0;

  /*
   * Local fusion
   */

  double *point_box = NULL;
  double _point_box[6];

  if (pm->char_length != NULL) {
    point_box = _point_box;
  }

  pm->max_n_points = 0;
  for (int i = 0; i < pm->n_point_clouds; i++) {
    pm->max_n_points = PDM_MAX (pm->max_n_points, pm->n_points[i]);
  }
  for (int i = 0; i < pm->n_point_clouds; i++) {
    PDM_octree_seq_t *octree_seq = PDM_octree_seq_create (1,
                                                          pm->depth_max,
                                                          pm->points_in_leaf_max,
                                                          pm->tolerance);

    PDM_octree_seq_point_cloud_set (octree_seq, 0,
                                    pm->n_points[i], pm->point_clouds[i]);
    PDM_octree_seq_build (octree_seq);

    const int root_id = PDM_octree_seq_root_node_id_get (octree_seq);

    const double *_char_length = NULL;
    if (pm->char_length != NULL) {
      _char_length = pm->char_length[i];
    }

    for (int j = i; j < pm->n_point_clouds; j++) {


      for (int k = 0; k < pm->n_points[j]; k++) {
        const double *_coord = pm->point_clouds[j] + 3 * k;
        if (point_box != NULL) {
          double char_length_point = pm->char_length[j][k];
          double tolerance = pm->tolerance;
          point_box[0] = _coord[0] - tolerance * char_length_point;
          point_box[1] = _coord[1] - tolerance * char_length_point;
          point_box[2] = _coord[2] - tolerance * char_length_point;
          point_box[3] = _coord[0] + tolerance * char_length_point;
          point_box[4] = _coord[1] + tolerance * char_length_point;
          point_box[5] = _coord[2] + tolerance * char_length_point;
        }

        _search_local_couple (&local_couple, &n_local_couple, &s_local_couple,
                              j, k, _coord, point_box, i, octree_seq,
                              root_id, pm->point_clouds[i], _char_length,
                              pm->tolerance);

      }
    }

    PDM_octree_seq_free (octree_seq);

  }

  /*
   * Distant fusion
   *   - Send/recv points between processes candidates
   *   - Search points candidates in the octree
   *   - Update distant table couple
   */

  double *extents_proc;
  int    *used_ranks;

  const int n_used_ranks = PDM_octree_processes_extents_get (pm->octree,
                                                             &used_ranks, &extents_proc);

  int s_tmp_store = sizeof(int) * pm->max_n_points;
  int n_tmp_store = 0;

  int *tmp_store = malloc (sizeof(int) * s_tmp_store * 3);

//  printf ("Extents proc \n");
//  for (int k = 0; k < n_rank ; k++) {
//    const double *_extents_proc = extents_proc + k * 6;
//    printf("%d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",k, _extents_proc[0],
//            _extents_proc[1],
//            _extents_proc[2],
//            _extents_proc[3],
//            _extents_proc[4],
//            _extents_proc[5]
//           );
//  }

  for (int i_cloud = 0; i_cloud < pm->n_point_clouds; i_cloud++) {
    int n_points = pm->n_points[i_cloud];
    const double *_coord = pm->point_clouds[i_cloud];
    for (int i = 0; i < n_points; i++) {
      const double *__coord = _coord + 3 * i;
      double box[6];

      if (pm->char_length != NULL) {
        box[0] = __coord[0] - pm->char_length[i_cloud][i] * pm->tolerance;
        box[1] = __coord[1] - pm->char_length[i_cloud][i] * pm->tolerance;
        box[2] = __coord[2] - pm->char_length[i_cloud][i] * pm->tolerance;
        box[3] = __coord[0] + pm->char_length[i_cloud][i] * pm->tolerance;
        box[4] = __coord[1] + pm->char_length[i_cloud][i] * pm->tolerance;
        box[5] = __coord[2] + pm->char_length[i_cloud][i] * pm->tolerance;
      }
      else {
        box[0] = __coord[0] - _default_eps;
        box[1] = __coord[1] - _default_eps;
        box[2] = __coord[2] - _default_eps;
        box[3] = __coord[0] + _default_eps;
        box[4] = __coord[1] + _default_eps;
        box[5] = __coord[2] + _default_eps;
      }

      for (int k = 0; k < n_used_ranks ; k++) {
        const int curr_rank  = used_ranks[k];
        const double *_extents_proc = extents_proc + k * 6;
        if (curr_rank != i_rank) {

          if (_intersect_extents(box, _extents_proc)) {
            if (n_tmp_store >= s_tmp_store) {
              s_tmp_store *= 2;
              tmp_store = realloc (tmp_store, sizeof(int) * s_tmp_store * 3);
            }
            tmp_store[3*n_tmp_store]   = curr_rank;
            tmp_store[3*n_tmp_store+1] = i_cloud;
            tmp_store[3*n_tmp_store+2] = i;
            n_tmp_store += 1;
          }
        }
      }
    }
  }

  int *val_send_n = PDM_array_zeros_int(n_rank);

  for (int i = 0; i < n_tmp_store; i++) {
    val_send_n[tmp_store[3*i]]++;
  }

  int *val_recv_n = malloc (sizeof(int)*n_rank);
  PDM_MPI_Alltoall (val_send_n, 1, PDM_MPI_INT, val_recv_n, 1, PDM_MPI_INT, pm->comm);

  // Envoi des points + char length en option sur les autres procs (test bounding box)

  int *val_send_idx = malloc (sizeof(int)*(n_rank+1));
  int *val_recv_idx = malloc (sizeof(int)*(n_rank+1));

  int _stride = 3 * 8 + 4 + 4; /* Coords + icloud + ipoint */
  if (pm->char_length != NULL) {
    _stride += 8; /* char_length */
  }

  for (int i = 0; i < n_rank; i++) {
    val_send_n[i] *= _stride;
    val_recv_n[i] *= _stride;
  }

  val_send_idx[0] = 0;
  val_recv_idx[0] = 0;
  for (int i = 0; i < n_rank; i++) {
    val_send_idx[i+1] = val_send_idx[i] + val_send_n[i];
    val_recv_idx[i+1] = val_recv_idx[i] + val_recv_n[i];
    val_send_n[i] = 0;
  }

  unsigned char *val_send =
        malloc (sizeof(unsigned char) * val_send_idx[n_rank]);
  unsigned char *val_recv =
        malloc (sizeof(unsigned char) * val_recv_idx[n_rank]);

  for (int i = 0; i < n_tmp_store; i++) {
    int iproc   = tmp_store[3*i];
    int i_cloud = tmp_store[3*i+1];
    int i_point = tmp_store[3*i+2];

    double *_coord = (double *) pm->point_clouds[i_cloud] + 3 * i_point;
    double *_tmp_val_double = (double *) (val_send + val_send_idx[iproc] + val_send_n[iproc]);

    _tmp_val_double[0] = _coord[0];
    _tmp_val_double[1] = _coord[1];
    _tmp_val_double[2] = _coord[2];

    val_send_n[iproc] += 24;

    if (pm->char_length != NULL) {
      double _char_length = pm->char_length[i_cloud][i_point];
      _tmp_val_double[3] = _char_length;
      val_send_n[iproc] += 8;
    }

    int *_tmp_val_int = (int *) (val_send + val_send_idx[iproc] + val_send_n[iproc]);

    _tmp_val_int[0] = i_cloud;
    _tmp_val_int[1] = i_point;

    val_send_n[iproc] += 8;

  }

  free (tmp_store);

  PDM_MPI_Alltoallv(val_send, val_send_n, val_send_idx, PDM_MPI_UNSIGNED_CHAR,
                    val_recv, val_recv_n, val_recv_idx, PDM_MPI_UNSIGNED_CHAR,
                    pm->comm);

  free (val_send);
  free (val_send_idx);

  int *n_fusion_from_proc = val_send_n;

  PDM_array_reset_int(n_fusion_from_proc, n_rank, 0);

  int *distant_couple  = NULL;
  int n_distant_couple = 0;
  int s_distant_couple = 0;

  double distant_coord[3];

  unsigned char *_tmp_recv = val_recv;

  for (int i = 0; i < n_rank; i++) {
    int _deb = val_recv_idx[i] / _stride;
    int _end = _deb + val_recv_n[i] / _stride;

    for (int j = _deb; j < _end; j++) {
      distant_coord[0] = *((double *) _tmp_recv);
      _tmp_recv += 8;
      distant_coord[1] = *((double *) _tmp_recv);
      _tmp_recv += 8;
      distant_coord[2] = *((double *) _tmp_recv);
      _tmp_recv += 8;
      double _char_length = -1;
      if (pm->char_length != NULL) {
        _char_length  = *((double *) _tmp_recv);
        _tmp_recv += 8;
        point_box[0] = distant_coord[0] - _char_length * pm->tolerance;
        point_box[1] = distant_coord[1] - _char_length * pm->tolerance;
        point_box[2] = distant_coord[2] - _char_length * pm->tolerance;
        point_box[3] = distant_coord[0] + _char_length * pm->tolerance;
        point_box[4] = distant_coord[1] + _char_length * pm->tolerance;
        point_box[5] = distant_coord[2] + _char_length * pm->tolerance;
      }

      int distant_cloud = *((int *) _tmp_recv);
      _tmp_recv += 4;
      int distant_point = *((int *) _tmp_recv);
      _tmp_recv += 4;

      int root_id = PDM_octree_root_node_id_get (pm->octree);

      _search_distant_couple (n_fusion_from_proc,
                              &distant_couple,
                              &n_distant_couple,
                              &s_distant_couple,
                              i,
                              distant_cloud,
                              distant_point,
                              distant_coord,
                              point_box,
                              pm->octree,
                              root_id,
                              pm->point_clouds,
                              pm->char_length,
                              pm->tolerance);

    }
  }

  free (val_recv);
  free (val_recv_idx);

  /*
   *  Check if the number of couple is coherent between other processes
   *
   */

  int *n_fusion_with_proc = val_recv_n;

  PDM_MPI_Alltoall (n_fusion_from_proc, 1, PDM_MPI_INT,
                    n_fusion_with_proc, 1, PDM_MPI_INT,
                    pm->comm);

  for (int i = 0; i < n_rank; i++) {
    assert (n_fusion_with_proc[i] == n_fusion_from_proc[i]);
  }

  free (n_fusion_from_proc);
  free (n_fusion_with_proc);

  /*
   * Build candidates_idx and candidates_desc arrays
   * from distant_couples and local_couples arrays
   *
   */

  for (int i = 0; i < pm->n_point_clouds; i++) {
    pm->candidates_idx[i] = malloc(sizeof(int) * (pm->n_points[i] + 1));
    pm->candidates_desc[i] = NULL;
  }

  for (int i = 0; i < pm->n_point_clouds; i++) {
    int *_candidates_idx = pm->candidates_idx[i];
    int _n_points = pm->n_points[i];

    PDM_array_reset_int(_candidates_idx, _n_points + 1, 0);
  }

  for (int i = 0; i < n_local_couple; i++) {
    int first_cloud  = local_couple[4*i    ];
    int first_index  = local_couple[4*i + 1];
    int second_cloud = local_couple[4*i + 2];
    int second_index = local_couple[4*i + 3];

    pm->candidates_idx[first_cloud][first_index+1]++;
    pm->candidates_idx[second_cloud][second_index+1]++;

  }

  for (int i = 0; i < n_distant_couple; i++) {
    int local_cloud = distant_couple[5*i    ];
    int local_index = distant_couple[5*i + 1];
//    int point_proc  = distant_couple[5*i + 2];
//    int point_cloud = distant_couple[5*i + 3];
//    int point_idx   = distant_couple[5*i + 4];

    pm->candidates_idx[local_cloud][local_index+1]++;

  }

  int **candidates_n = malloc (sizeof(int*) * pm->n_point_clouds);

  for (int i = 0; i < pm->n_point_clouds; i++) {
    int *_candidates_idx = pm->candidates_idx[i];
    int _n_points = pm->n_points[i];

    candidates_n[i] = malloc(sizeof(int) * _n_points);
    int *_candidates_n = candidates_n[i];

    for (int j = 0; j < _n_points; j++) {
      _candidates_idx[j+1] += _candidates_idx[j];
      _candidates_n[j] = 0;
    }

    pm->candidates_desc[i] = malloc(sizeof(int)*_candidates_idx[_n_points]*3);

  }

  int iproc;
  PDM_MPI_Comm_rank (pm->comm, &iproc);

  for (int i = 0; i < n_local_couple; i++) {
    int first_cloud  = local_couple[4*i    ];
    int first_index  = local_couple[4*i + 1];
    int second_cloud = local_couple[4*i + 2];
    int second_index = local_couple[4*i + 3];

    int idx = pm->candidates_idx[first_cloud][first_index] +
              candidates_n[first_cloud][first_index];

    pm->candidates_desc[first_cloud][3*idx]     = iproc;
    pm->candidates_desc[first_cloud][3*idx + 1] = second_cloud;
    pm->candidates_desc[first_cloud][3*idx + 2] = second_index;

    candidates_n[first_cloud][first_index]++;

    idx = pm->candidates_idx[second_cloud][second_index] +
               candidates_n[second_cloud][second_index];

    pm->candidates_desc[second_cloud][3*idx]     = iproc;
    pm->candidates_desc[second_cloud][3*idx + 1] = first_cloud;
    pm->candidates_desc[second_cloud][3*idx + 2] = first_index;

    candidates_n[second_cloud][second_index]++;

  }

  for (int i = 0; i < n_distant_couple; i++) {
    int local_cloud = distant_couple[5*i    ];
    int local_index = distant_couple[5*i + 1];
    int point_proc  = distant_couple[5*i + 2];
    int point_cloud = distant_couple[5*i + 3];
    int point_idx   = distant_couple[5*i + 4];

    int idx = pm->candidates_idx[local_cloud][local_index] +
              candidates_n[local_cloud][local_index];

    pm->candidates_desc[local_cloud][3*idx]     = point_proc;
    pm->candidates_desc[local_cloud][3*idx + 1] = point_cloud;
    pm->candidates_desc[local_cloud][3*idx + 2] = point_idx;

    candidates_n[local_cloud][local_index]++;

  }

  /* Free local data */

  for (int i = 0; i < pm->n_point_clouds; i++) {
    free (candidates_n[i]);
  }

  free (candidates_n);
  free (distant_couple);
  free (local_couple);

  PDM_octree_free (pm->octree);

}


/**
 *
 * \brief Get candidates to merge for each point
 *
 * \param [in]   pm              Pointer to \ref PDM_points_merge object
 * \param [in]   i_point_cloud   Current cloud
 * \param [out]  candidates_idx  Indexes of candidate for each current cloud point
 *                               (size = number of points in the current cloud + 1)
 * \param [out]  candidates_desc Candidates description (process,
 *                                                       cloud in the process,
 *                                                       point in the cloud)
 *
 */

void
PDM_points_merge_candidates_get
(
       PDM_points_merge_t  *pm,
 const int                  i_point_cloud,
       int                **candidates_idx,
       int                **candidates_desc
)
{

  assert(pm->candidates_idx  != NULL);
  assert(pm->candidates_desc != NULL);

  *candidates_idx  = pm->candidates_idx [i_point_cloud];
  *candidates_desc = pm->candidates_desc[i_point_cloud];

  pm->results_is_getted = PDM_TRUE;

  if (0 == 1) {
    printf("candidates : \n");
    for (int i = 0; i < pm->n_points[i_point_cloud]; i++) {
      if (pm->candidates_idx[i_point_cloud][i+1] > pm->candidates_idx[i_point_cloud][i]) {
      printf("-- %d %d ", i_point_cloud, i);
      for (int j = pm->candidates_idx[i_point_cloud][i];
               j <  pm->candidates_idx[i_point_cloud][i+1]; j++) {

        printf(" : %d", pm->candidates_desc[i_point_cloud][3*j]);
        printf(" %d", pm->candidates_desc[i_point_cloud][3*j+1]);
        printf(" %d", pm->candidates_desc[i_point_cloud][3*j+2]);
      }
      printf("\n");
      }
    }
  }
}

/**
 *
 * \brief Get size of the resulting array
 *
 * \param [in]   pm                Pointer to \ref PDM_points_merge object
 * \param [in]   i_point_cloud     Current cloud
 * \param [out]  n_point_cloud     Number of points in the current cloud
 * \param [out]  n_candidates_desc Size of candidates_desc = candidates_idx[n_point_cloud+1]
 *
 */
void
PDM_points_merge_candidates_size_get
(
       PDM_points_merge_t *pm,
 const int                i_point_cloud,
       int               *n_point_cloud,
       int               *n_candidates_desc
)
{

  assert(pm->candidates_idx  != NULL);
  assert(pm->candidates_desc != NULL);

  *n_point_cloud     = pm->n_points[i_point_cloud];
  *n_candidates_desc = pm->candidates_idx[i_point_cloud][pm->n_points[i_point_cloud]]; // ou x3 ?

}
