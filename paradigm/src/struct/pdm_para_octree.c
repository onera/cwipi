
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <sys/resource.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_morton.h"
#include "pdm_handles.h"
#include "pdm_para_octree.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define NTIMER 11
//#define NGB_ON_THE_FLY 1

const int NGB_ON_THE_FLY = 1;

/*============================================================================
 * Type definitions
 *============================================================================*/


/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                         = 0,
  BUILD_ORDER_POINTS            = 1,
  BUILD_BLOCK_PARTITION         = 2,
  BUILD_LOCAL_NODES             = 3,
  BUILD_LOCAL_NEIGHBOURS_STEP1  = 4,
  BUILD_LOCAL_NEIGHBOURS_STEP2  = 5,
  BUILD_LOCAL_NEIGHBOURS_STEP3  = 6,
  BUILD_LOCAL_NEIGHBOURS        = 7,
  BUILD_DISTANT_NEIGHBOURS      = 8,
  BUILD_TOTAL                   = 9,
  END                           = 10,

} _ol_timer_step_t;


/**
 * \struct _heap_t
 * \brief  Heap used to recursively subdivide nodes
 *
 */

typedef struct  {

  int   top;                  /*!< Top of head  */
  int   size;                 /*!< Size of heap */
  PDM_morton_code_t *codes;   /*!< Morton codes */
  int *range;                 /*!< Points range */
  int *n_points;              /*!< Points number */
  int   max_top;

} _heap_t;


/**
 * \struct _l_octant_t
 * \brief  Define a list of octants
 *
 */

typedef struct  {

  int   n_nodes;                 /*!< Current number of nodes in octree */
  int   n_nodes_max;             /*!< Maximum number of nodes in octree */

  PDM_morton_code_t *codes;        /*!< Morton codes */

  int  *n_points;          /*!< Number of points in octant*/
  int  *range;             /*!< Start index of point list for each octant */

  int   *neighbour_idx;
  int   *neighbours;               /*!< rank + id_node size = 2 * n_nodes */
  int   dim;

} _l_octant_t;


/**
 * \struct _octree_t
 * \brief  Define an octree
 *
 */

typedef struct  {

  double  global_extents[6];     /*!< Extents of current process */
  int     depth_max;             /*!< Maximum depth of the three */
  int     points_in_leaf_max;    /*!< Maximum number of points in a leaf */
  double      s[3];           /*!< Translation for the normalization */
  double      d[3];           /*!< Dilatation for the normalization */

  int     n_point_clouds;        /*!< Number of point cloud */

  PDM_g_num_t        t_n_points;     /*!< total number of points */
  int                n_points;       /*!< Number of points in each cloud */
  double            *points;         /*!< Point coordinates */
  int               *points_icloud;  /*!< Point cloud */
  PDM_g_num_t       *points_gnum;    /*!< Point global number */
  PDM_morton_code_t *points_code;    /*!< Morton codes */

  PDM_morton_code_t *rank_octants_index;
  _l_octant_t *octants;       /*!< list of octants */

  PDM_MPI_Comm comm;           /*!< MPI communicator */
  int   dim;                     /*!< Dimension */

  int  n_part_boundary_elt;    /*!< Number of partitioning boundary element */
  int *part_boundary_elt_idx; /*!< Index for part_boundary_elt (size=\ref n_part_boundary_elt + 1 */
  int *part_boundary_elt;     /*!< Partitioning boundary elements description (proc number + element number) */

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */

  int neighboursToBuild;

  int  n_connected;
  int *connected_idx;

} _octree_t;



/**
 * \struct _neighbours_tmp_t
 * \brief  Define a temporary neighbour structure
 *
 */


typedef struct  {

  int n_neighbour[6];     /*!< Number of neighbours in the arrays  */
  int s_neighbour[6];     /*!< Size of arrays */
  int *neighbours[6];     /*!< Arrays */

} _neighbours_tmp_t;



/**
 * \struct _min_heap_t
 * \brief  Binary heap used as (min-)priority queue
 *
 */


typedef struct {

  int          size;
  int          count;
  int         *lnum;
  PDM_g_num_t *gnum;
  double      *priority;

} _min_heap_t;



/*static const int _2d_sibling_neighbours[4][4] = {{-1, 2,-1, 1},
						 {-1, 3, 0,-1},
						 { 0,-1,-1, 3},
						 { 1,-1, 2,-1}};*/

static const int _3d_sibling_neighbours[8][6] = {{-1, 4,-1, 2,-1, 1},
						 {-1, 5,-1, 3, 0,-1},
						 {-1, 6, 0,-1,-1, 3},
						 {-1, 7, 1,-1, 2,-1},
						 { 0,-1,-1, 6,-1, 5},
						 { 1,-1,-1, 7, 4,-1},
						 { 2,-1, 4,-1,-1, 7},
						 { 3,-1, 5,-1, 6,-1}};

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_octrees    = NULL;

//static const double _eps_default  = 1.e-12;

static const PDM_morton_int_t max_morton_level = 15;
//static const int max_morton_level = 2;

/*============================================================================
 * Private function definitions
 *============================================================================*/
static int
_binary_search
(
 const int  elem,
 const int *array,
 const int  n,
 int       *in_array
 )
{
  int l = 0;
  int r = n;

  *in_array = 0;

  if (n < 1)
    return 0;

  while (l + 1 < r) {
    int m = l + (r - l)/2;

    if (elem < array[m])
      r = m;
    else
      l = m;
  }

  if (array[l] == elem) {
    *in_array = 1;
    return l;

  } else if (array[l] < elem)
    return l + 1;

  else
    return l;
}


static int
_binary_search_long
(
 const PDM_g_num_t  elem,
 const PDM_g_num_t *array,
 const int          n,
 int               *in_array
 )
{
  int l = 0;
  int r = n;

  *in_array = 0;

  if (n < 1)
    return 0;

  while (l + 1 < r) {
    int m = l + (r - l)/2;

    if (elem < array[m])
      r = m;
    else
      l = m;
  }

  if (array[l] == elem) {
    *in_array = 1;
    return l;

  } else if (array[l] < elem)
    return l + 1;

  else
    return l;
}



static _min_heap_t *
_min_heap_create
(
 const int size
 )
{
  _min_heap_t *h = malloc (sizeof(_min_heap_t));
  h->size        = size;
  h->count       = 0;
  h->lnum        = malloc (sizeof(int)         * size);
  h->gnum        = malloc (sizeof(PDM_g_num_t) * size);
  h->priority    = malloc (sizeof(double)      * size);

  for (int i = 0; i < h->size; i++)
    h->priority[i] = HUGE_VAL;

  return h;
}


static void
_min_heap_reset
(
 _min_heap_t *h
 )
{
  h->count = 0;
}


static void
_min_heap_free
(
 _min_heap_t *h
 )
{
  free (h->lnum);
  free (h->gnum);
  free (h->priority);
  free (h);
}



static void
_min_heap_swap
(
 _min_heap_t *h,
 const int    a,
 const int    b
 )
{
  int l          = h->lnum[a];
  PDM_g_num_t g  = h->gnum[a];
  double p       = h->priority[a];

  h->lnum[a]     = h->lnum[b];
  h->gnum[a]     = h->gnum[b];
  h->priority[a] = h->priority[b];

  h->lnum[b]     = l;
  h->gnum[b]     = g;
  h->priority[b] = p;
}



static void
_min_heap_push
(
 _min_heap_t       *h,
 const int          lnum,
 const PDM_g_num_t  gnum,
 const double       priority
 )
{
  int i = 0;
  /* make sure the heap is large enough to contain the new element */
  if (h->count >= h->size) {
    h->size = (h->size) ? 2*h->size : 10;
    h->lnum =     (int *)         realloc (h->lnum,     sizeof(int)         * h->size);
    h->gnum =     (PDM_g_num_t *) realloc (h->gnum,     sizeof(PDM_g_num_t) * h->size);
    h->priority = (double *)      realloc (h->priority, sizeof(double)      * h->size);
    for (i = h->count+1; i < h->size; i++)
      h->priority[i] = HUGE_VAL;
  }

  i = h->count;
  h->count++;
  h->lnum[i]     = lnum;
  h->gnum[i]     = gnum;
  h->priority[i] = priority;

  int parent = (i - 1)/2;
  while (i > 0 && h->priority[parent] > priority) {
    _min_heap_swap (h, parent, i);
    i = parent;
    parent = (i - 1)/2;
  }
}



static void
_min_heap_heapify
(
 _min_heap_t *h,
 const int    i
 )
{
  int l = 2*i+1; // left child node
  int r = 2*i+2; // right child node
  int s = i; // node with smallest priority

  if (l < h->count && h->priority[l] < h->priority[i])
    s = l;

  if (r < h->count && h->priority[r] < h->priority[s])
    s = r;

  if (s != i) {
    _min_heap_swap (h, i, s);
    _min_heap_heapify (h, s);
  }
}


static int
_min_heap_pop
(
 _min_heap_t *h,
 int         *lnum,
 PDM_g_num_t *gnum,
 double      *priority
 )
{
  if (h->count < 1)
    return 0;

  *lnum     = h->lnum[0];
  *gnum     = h->gnum[0];
  *priority = h->priority[0];

  h->count--;
  h->lnum[0]     = h->lnum[h->count];
  h->gnum[0]     = h->gnum[h->count];
  h->priority[0] = h->priority[h->count];

  _min_heap_heapify (h, 0);

  return 1;
}




inline static double
_octant_min_dist2
(
 const int          dim,
 PDM_morton_code_t  code,
 const double      *d,
 const double      *s,
 const double      *coords
 )
{
  double min_dist2 = 0., delta = 0.;
  //  double side = 1./pow(2, code.L);
  double side = 1./(double)(1 << code.L);

  for (int i = 0; i < dim; i++) {
    double x = coords[i];
    double xmin = s[i] + d[i] * side * code.X[i];
    double xmax = xmin + d[i] * side;

    if (x > xmax) {
      delta = x - xmax;
      min_dist2 += delta * delta;
    } else if (x < xmin) {
      delta = x - xmin;
      min_dist2 += delta * delta;
    }
  }

  return min_dist2;
}


inline static double
_octant_min_dist2_normalized
(
 const int          dim,
 PDM_morton_code_t  code,
 const double      *d,
 const double      *coords
 )
{
  double min_dist2 = 0., delta = 0.;
  //  double side = 1./pow(2, code.L);
  double side = 1./(double)(1 << code.L);

  for (int i = 0; i < dim; i++) {
    double x = coords[i];
    double xmin = side * code.X[i];
    double xmax = xmin + side;

    if (x > xmax) {
      delta = d[i] * (x - xmax);
      min_dist2 += delta * delta;
    } else if (x < xmin) {
      delta = d[i] * (x - xmin);
      min_dist2 += delta * delta;
    }
  }

  return min_dist2;
}




inline static int
_get_octant_part_id
(
 const _octree_t *octree,
 const int        id_octant
 )
{
  for (int i_part = 0; i_part < octree->n_connected; i_part++) {
    if (octree->connected_idx[i_part] <= id_octant &&
        id_octant < octree->connected_idx[i_part+1]) {
      return i_part;
    }
  }

  return -1;
}


/**
 *
 * \brief Create a heap
 *
 * \param [in]  Size    Size of heap
 *
 * \return   a new heap
 *
 */

static _heap_t *
_heap_create
(
 const int size
 )
{
  _heap_t *heap = malloc(sizeof(_heap_t));
  heap->top = 0;
  heap->max_top = 0;
  heap->size = size;
  heap->codes = malloc(sizeof(PDM_morton_code_t) * size);
  heap->range =  malloc(sizeof(int) * size);
  heap->n_points =  malloc(sizeof(int) * size);
  return heap;
}


/**
 *
 * \brief Free a heap
 *
 * \param [in]  heap   Heap to free
 *
 * \return NULL
 *
 */

static _heap_t *
_heap_free
(
 _heap_t *heap
 )
{
  free (heap->codes);
  free (heap->range);
  free (heap->n_points);
  free (heap);
  return NULL;
}


/**
 *
 * \brief Push a new element in the heap
 *
 * \param [inout]  heap      Heap
 * \param [in]     code      Morton code
 * \param [in]     range     Range
 * \param [in]     n_points  Number of points
 *
 * \return  1 if pushed 0 otherwise
 *
 */

static int
_heap_push
(
 _heap_t *heap,
 const PDM_morton_code_t code,
 const int range,
 const int n_points
 )
{
  if (heap->top >= heap->size) {
    return 0;
  }
  int idx = heap->top;
  PDM_morton_copy (code, &(heap->codes[idx]));
  heap->range[idx] = range;
  heap->n_points[idx] = n_points;
  heap->top++;
  heap->max_top = PDM_MAX(heap->top, heap->max_top);
  return 1;
}


/**
 *
 * \brief Pull top element of the heap
 *
 * \param [inout]  heap      Heap
 * \param [out]    code      Morton code
 * \param [out]    range     Range
 * \param [out]    n_points  Number of points
 *
 * \return  1 if pulled 0 otherwise
 *
 */

static int
_heap_pull
(
 _heap_t *heap,
 PDM_morton_code_t *code,
 int *range,
 int *n_points
 )
{
  heap->top--;
  if (heap->top < 0) {
    return 0;
  }
  int idx = heap->top;
  PDM_morton_copy (heap->codes[idx], code);
  *range = heap->range[idx];
  *n_points = heap->n_points[idx];
  return 1;
}

/**
 *
 * \brief Neighbour
 *
 * \param [inout]   octants     Octants
 *
 * \return NULL
 *
 */

static PDM_para_octree_direction_t
_inv_direction
(
 PDM_para_octree_direction_t direc
 )
{
  if (direc == PDM_BOTTOM) {
    return PDM_UP;
  }
  else if (direc == PDM_UP) {
    return PDM_BOTTOM;
  }
  else if (direc == PDM_SOUTH) {
    return PDM_NORTH;
  }
  else if (direc == PDM_NORTH) {
    return PDM_SOUTH;
  }
  else if (direc == PDM_WEST) {
    return PDM_EAST;
  }
  else if (direc == PDM_EAST) {
    return PDM_WEST;
  }
  else {
    abort();
  }

}

/**
 *
 * \brief Neighbour
 *
 * \param [inout]   octants     Octants
 *
 * \return neighbour or NULL ()
 *
 */

static PDM_morton_code_t *
_neighbour
(
 PDM_morton_code_t code,
 PDM_para_octree_direction_t direction
 )
{
  const int dim = direction / 2;
  const int _direction = 2 * (direction % 2) - 1;

  PDM_morton_code_t *neighbour = NULL;

  if (((_direction > 0) && (code.X[dim] < (unsigned int) ((1 << code.L) - 1))) ||
      ((_direction < 0) && (code.X[dim] > 0))) {

    neighbour = malloc(sizeof(PDM_morton_code_t));

    neighbour->L = code.L;
    neighbour->X[0] = code.X[0];
    neighbour->X[1] = code.X[1];
    neighbour->X[2] = code.X[2];

    neighbour->X[dim] = code.X[dim] + _direction;
  }

  return neighbour;
}


/**
 *
 * \brief Free octants
 *
 * \param [inout]   octants     Octants
 *
 * \return NULL
 *
 */

static void
_octants_purge
(
 _l_octant_t *octants
 )
{
  octants->n_nodes_max = 0;
  octants->n_nodes     = 0;

  if (octants->codes != NULL) {
    free (octants->codes);
  }

  if (octants->n_points != NULL) {
    free (octants->n_points);
  }

  if (octants->range != NULL) {
    free (octants->range);
  }

  if (octants->neighbour_idx != NULL) {
    free (octants->neighbour_idx);
  }

  if (octants->neighbours != NULL) {
    free (octants->neighbours);
  }
}

/**
 *
 * \brief Free octants
 *
 * \param [inout]   octants     Octants
 *
 * \return NULL
 *
 */

static _l_octant_t *
_octants_free
(
 _l_octant_t *octants
 )
{

  _octants_purge (octants);

  free(octants);
  return NULL;
}


/**
 *
 * \brief Initialize list of octants
 *
 * \param [inout]   octants     Octants
 * \param [in]      octant_dim  Dimension of an octant
 * \param [in]      init_size   Initial size of octants
 *
 */

static void
_octants_init
(
 _l_octant_t *octants,
 const int   octant_dim,
 const int   init_size
 )
{
  octants->n_nodes_max = init_size;
  octants->n_nodes     = 0;

  octants->codes    = malloc (sizeof(PDM_morton_code_t) * octants->n_nodes_max);
  octants->n_points = malloc (sizeof(int) * octants->n_nodes_max);
  octants->range    = malloc (sizeof(int) * (octants->n_nodes_max+1));

  octants->neighbour_idx = NULL;
  octants->neighbours    = NULL;
  octants->dim = octant_dim;
}


/**
 *
 * \brief Check size of the size of a list of octants
 *
 * \param [in]   octants       Octants
 * \param [in]   n_free_node   Number of required fre nodes
 *
 */

static int
_octants_check_alloc
(
 _l_octant_t *octants,
 const int n_free_node
 )
{
  int is_realloc = 0;
  if (octants->n_nodes + n_free_node > octants->n_nodes_max) {

    octants->n_nodes_max *= 2;

    octants->codes    = realloc (octants->codes,
                                 sizeof(PDM_morton_code_t) * octants->n_nodes_max);
    octants->n_points = realloc (octants->n_points,
                                 sizeof(int) * octants->n_nodes_max);
    octants->range = realloc (octants->range,
                              sizeof(int) * (octants->n_nodes_max+1));
    octants->neighbour_idx = NULL;
    octants->neighbours    = NULL;

    is_realloc = 1;
  }
  return is_realloc;
}


/**
 *
 * \brief Push back a octant to a list of octants
 *
 * \param [in]   octants     Octants
 *
 */

static void
_octants_push_back
(
 _l_octant_t *octants,
 const PDM_morton_code_t code,
 const int n_points,
 const int range
 )
{

  _octants_check_alloc (octants, 1);

  const int idx = octants->n_nodes;

  PDM_morton_copy (code, octants->codes + idx);

  octants->n_points[idx] = n_points;

  octants->range[idx] = range;

  octants->n_nodes += 1;

}


/**
 *
 * \brief Push front a octant to a list of octants
 *
 * \param [in]   octants     Octants
 *
 */

static void
_octants_push_front
(
 _l_octant_t *octants,
 const PDM_morton_code_t code,
 const int n_points,
 const int range
 )
{

  _octants_check_alloc (octants, 1);

  for (int i = octants->n_nodes; i > 0; i--) {

    PDM_morton_copy (octants->codes[i - 1], octants->codes + i);

    octants->n_points[i] =  octants->n_points[i-1];

    octants->range[i] = octants->range[i-1];
  }

  const int idx = 0;

  PDM_morton_copy (code, octants->codes + idx);

  octants->n_points[idx] = n_points;

  octants->range[idx] = range;

  octants->n_nodes += 1;

}


/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _octree_t *
_get_from_id
(
 int  id
 )
{
  _octree_t *octree = (_octree_t *) PDM_Handles_get (_octrees, id);

  if (octree == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_octree error : Bad identifier\n");
  }

  return octree;
}


/**
 *
 * \brief Build minimal octree between two octants
 *
 * \param [in]     octree    Current octree
 * \param [in]     code      Morton code
 * \param [inout]  extents   Extents associated to the Morton code
 *
 */

/* static void */
/* g_extents */
/* ( */
/*  _octree_t *octree, */
/*  PDM_morton_code_t code, */
/*  double    extents[] */
/* ) */
/* { */
/*   for (int i = 0; i < octree->dim; i++) { */
/*     extents[i] = */
/*       ((double) code.X[i]/((double) pow(2,code.L)))* octree->d[i] + octree->s[i]; */
/*     extents[octree->dim + i] = */
/*       (((double) code.X[i] + 1)/((double) pow(2,code.L))) * octree->d[i] + octree->s[i]; */
/*   } */
/* } */

/**
 *
 * \brief Remove duplicates
 *
 * \param [in]  octants List of octants
 *
 * \return octants without duplicates
 *
 */

static _l_octant_t *
_remove_duplicates
(
 _l_octant_t *octants
 )
{
  PDM_morton_code_t *_codes = octants->codes;
  _l_octant_t *r_octants = malloc(sizeof(_l_octant_t));

  const int dim = 3;

  _octants_init (r_octants, dim, octants->n_nodes);

  PDM_morton_code_t prev_code;

  prev_code.L = -1;
  prev_code.X[0] = 0;
  prev_code.X[1] = 0;
  prev_code.X[2] = 0;

  for (int i = 0; i < octants->n_nodes; i++) {

    if (_codes[i].L == prev_code.L) {
      if ((prev_code.X[0] == _codes[i].X[0]) &&
          (prev_code.X[1] == _codes[i].X[1]) &&
          (prev_code.X[2] == _codes[i].X[2])) {

        break;
      }
    }

    prev_code.L    = _codes[i].L;
    prev_code.X[0] = _codes[i].X[0];
    prev_code.X[1] = _codes[i].X[1];
    prev_code.X[2] = _codes[i].X[2];

    _octants_push_back (r_octants,
                        _codes[i],
                        0,
                        0);

  }

  return r_octants;
}

/**
 *
 * \brief Removing overlaps from a sorted lis of octants
 *
 * \param [inout]  octants A lis of octants
 *
 */

static _l_octant_t *
_linearize
(
 _l_octant_t *octants
 )
{
  PDM_morton_code_t *_codes = octants->codes;
  _l_octant_t *r_octants = malloc(sizeof(_l_octant_t));

  const int dim = 3;

  _octants_init (r_octants, dim, octants->n_nodes);

  if  (octants->n_nodes > 0) {
    for (int i = 0; i < octants->n_nodes - 1; i++) {

      if (!PDM_morton_ancestor_is (_codes[i], _codes[i+1])) {
        _octants_push_back (r_octants,
                            _codes[i],
                            0,
                            0);
      }

    }

    _octants_push_back (r_octants,
                        _codes[octants->n_nodes-1],
                        0,
                        0);
  }

  return r_octants;
}


/**
 *
 * \brief Constructing a minimal linear octree between two octants
 *
 * \param [in]  a     Morton code a
 * \param [in]  b     Morton code b
 *
 * \return octants The minimal linear octree between a and b or NULL if a >= b
 *
 */

static _l_octant_t *
_complete_region
(
 PDM_morton_code_t a,
 PDM_morton_code_t b
 )
{
  const int dim = 3;

  _l_octant_t *r_octants = NULL;

  if (PDM_morton_a_gt_b (b, a)) {

    _l_octant_t *w_octants = malloc(sizeof(_l_octant_t));
    r_octants = malloc(sizeof(_l_octant_t));

    _octants_init (w_octants, dim, 4);
    _octants_init (r_octants, dim, 4);

    /* printf("_complete_region\n"); */

    /* printf("a_d\n"); */
    /* PDM_morton_dump(3, a); */
    /* printf("a_f\n"); */

    /* printf("b_d\n"); */
    /* PDM_morton_dump(3, b); */
    /* printf("b_f\n"); */

    PDM_morton_code_t nca;
    PDM_morton_nearest_common_ancestor (a, b, &nca);

    const int n_child = 8;

    int  size = PDM_morton_max_level * 8;
    _heap_t *heap = _heap_create (size);

    PDM_morton_code_t children[8];
    PDM_morton_get_children(dim,
                            nca,
                            children);

    for (int i = n_child - 1; i >= 0; i--) {
      int is_pushed = _heap_push (heap,
                                  children[i],
                                  0,
                                  0);
      if (!is_pushed) {
        printf ("Internal error PDM_para_octree 1 : heap is full\n");
        exit(1);
      }
    }

    PDM_morton_code_t code;
    int range;
    int n_points;

    while (_heap_pull (heap, &code, &range, &n_points)) {
      if (PDM_morton_a_gt_b (code, a) &&
          PDM_morton_a_gt_b (b, code) &&
          !PDM_morton_ancestor_is (code, b)) {

        _octants_push_back (r_octants,
                            code,
                            0,
                            0);
      }

      else if ((PDM_morton_ancestor_is (code, b) ||
                PDM_morton_ancestor_is (code, a)) &&
               !((code.X[0] == a.X[0]) &&
                 (code.X[1] == a.X[1]) &&
                 (code.X[2] == a.X[2]) &&
                 (code.L == a.L)) &&
               !((code.X[0] == b.X[0]) &&
                 (code.X[1] == b.X[1]) &&
                 (code.X[2] == b.X[2]) &&
                 (code.L == b.L))) {

        PDM_morton_get_children(dim,
                                code,
                                children);

        for (int i = n_child - 1; i >= 0; i--) {
          int is_pushed = _heap_push (heap,
                                      children[i],
                                      0,
                                      0);

          if (!is_pushed) {
            printf ("Internal error PDM_para_octree 2 : heap is full\n");
            exit(1);
          }
        }
      }
    }

    _octants_free (w_octants);

    _heap_free (heap);
  }

  return r_octants;
}



/**
 *
 * \brief Redistribute octants
 *
 * \param [inout]  L             Distributed list of octants
 * \param [in]     morton_index  Morton index
 * \param [in]     comm          MPI communicator
 *
 */

static void
_distribute_octants
(
 _l_octant_t       *L,
 PDM_morton_code_t *morton_index,
 PDM_MPI_Comm       comm
 )
{
  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  int *send_count = malloc(sizeof(int) * n_ranks);
  size_t *send_shift = malloc(sizeof(size_t) * (n_ranks+1));

  int *recv_count = malloc(sizeof(int) * n_ranks);
  size_t *recv_shift = malloc(sizeof(size_t) * (n_ranks+1));

  for (int i = 0; i < n_ranks; i++) {
    send_count[i] = 0;
  }

  int irank = 0;
  for (int i = 0; i < L->n_nodes; i++) {
    if (PDM_morton_a_ge_b (L->codes[i], morton_index[irank+1])) {

      irank += 1 + PDM_morton_binary_search (n_ranks - (irank + 1),
                                             L->codes[i],
                                             morton_index + irank + 1);
    }
    send_count[irank] += L->dim + 1;
  }

  /* Exchange number of coords to send to each process */

  PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT,
                   recv_count, 1, PDM_MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
  }

  /* Build send and receive buffers */

  PDM_morton_int_t *send_codes =
    malloc (send_shift[n_ranks] * sizeof(PDM_morton_int_t));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  irank = 0;
  for (int i = 0; i < L->n_nodes; i++) {

    if (PDM_morton_a_ge_b (L->codes[i], morton_index[irank+1])) {

      irank += 1 + PDM_morton_binary_search(n_ranks - (irank + 1),
                                            L->codes[i],
                                            morton_index + irank + 1);
    }

    int shift = send_shift[irank] + send_count[irank];
    send_codes[shift++] = L->codes[i].L;

    for (int j = 0; j < L->dim; j++) {
      send_codes[shift++] = L->codes[i].X[j];
    }

    send_count[irank] += L->dim + 1;
  }

  PDM_morton_int_t * recv_codes = malloc (recv_shift[n_ranks] * sizeof(PDM_morton_int_t));

  /* - exchange codes between processes */

  PDM_MPI_Alltoallv_l(send_codes, send_count, send_shift, PDM_MPI_UNSIGNED,
                      recv_codes, recv_count, recv_shift, PDM_MPI_UNSIGNED,
                      comm);

  free (send_codes);
  free (send_count);
  free (send_shift);
  free (recv_count);

  /* - tri des codes recus */

  const int _dim = L->dim;
  _octants_purge (L);
  _octants_init (L, _dim, recv_shift[n_ranks]/(_dim + 1));

  size_t idx = 0;
  for (size_t i = 0; i < recv_shift[n_ranks]/4; i++) {
    PDM_morton_code_t _code;
    _code.L = recv_codes[idx++];
    for (int j = 0; j < L->dim; j++) {
      _code.X[j] = recv_codes[idx++];
    }
    _octants_push_back (L,
                        _code,
                        0,
                        0);
  }

  free (recv_shift);
  free (recv_codes);

  PDM_morton_local_sort (L->n_nodes, L->codes);

}


/**
 *
 * \brief Constructing a complete linear octree from partial set of octants
 *
 * \param [in]  L     Distributed list of octants
 * \param [in]  comm  MPI Communicator
 *
 * \return octants The complete linear octree
 *
 */

static _l_octant_t *
_complete_octree
(
 _l_octant_t *L,
 PDM_MPI_Comm comm
 )
{
  const int dim = 3;

  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank (comm, &rank);

  /* Remove duplicates */

  _l_octant_t *L1 = _remove_duplicates (L);

  /* Linearize */

  _l_octant_t *L2 = _linearize (L1);

  _octants_free (L1);

  PDM_g_num_t _n_nodes_global = 0;
  PDM_g_num_t _n_nodes_local =  L2->n_nodes;

  PDM_MPI_Allreduce (&_n_nodes_local, &_n_nodes_global, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

  _l_octant_t *R = malloc(sizeof(_l_octant_t));
  _octants_init (R, dim, PDM_MAX(L2->n_nodes, 1));

  if (_n_nodes_global > 0) {

    PDM_morton_code_t *L2_morton_index = malloc(sizeof(PDM_morton_code_t) * (n_ranks + 1));

    int *order = malloc (sizeof(int) * L2->n_nodes);
    int *weight = malloc (sizeof(int) * L2->n_nodes);

    for (int i = 0; i < L2->n_nodes; i++) {
      weight[i] = 1;
      order[i] = i;
    }

    PDM_morton_int_t max_level = 0;
    for (int i = 0; i < L2->n_nodes; i++) {
      max_level = PDM_MAX (L2->codes[i].L, max_level);
    }

    PDM_morton_int_t max_max_level;
    PDM_MPI_Allreduce(&max_level, &max_max_level, 1,
                      PDM_MPI_UNSIGNED, PDM_MPI_MAX, comm);

    /* printf("n L2 : %d %d\n" , rank, L2->n_nodes); */
    /* printf("\n-------------\nL2 avant d\n"); */
    /* for (int i = 0; i < L2->n_nodes; i++) { */
    /*   PDM_morton_dump (3, L2->codes[i]); */
    /* } */
    PDM_morton_ordered_build_rank_index (dim,
                                         max_max_level,
                                         L2->n_nodes,
                                         L2->codes,
                                         weight,
                                         L2_morton_index,
                                         comm);

    /* printf("\nL2 morton index d\n"); */
    /* for (int i = 0; i < n_ranks + 1; i++) { */
    /*   PDM_morton_dump (3,  L2_morton_index[i]); */
    /* } */
    /* printf("L2 morton, index f\n"); */

    free(weight);
    free(order);

    _distribute_octants (L2, L2_morton_index, comm);

    free (L2_morton_index);

    /* printf("n L2 : %d %d\n" , rank, L2->n_nodes); */
    /* printf("\nL2 d\n"); */
    /* for (int i = 0; i < L2->n_nodes; i++) { */
    /*   PDM_morton_dump (3, L2->codes[i]); */
    /* } */
    /* printf("L2 f\n--------------\n"); */

    /* PDM_MPI_Barrier(comm); */
    /* exit(1); */

    int *rank_n_nodes = malloc (sizeof(int) * n_ranks);

    PDM_MPI_Allgather(&L2->n_nodes, 1, PDM_MPI_INT,
                      rank_n_nodes, 1, PDM_MPI_INT,
                      comm);

    int first_rank = 0;
    while (first_rank < n_ranks-1
           && rank_n_nodes[first_rank] == 0) {
      first_rank++;
    }

    int last_rank = n_ranks-1;
    while (last_rank > 0
           && rank_n_nodes[last_rank] == 0) {
      last_rank--;
    }

    int next_rank = rank + 1;
    while (next_rank < n_ranks-1
           && rank_n_nodes[next_rank] == 0) {
      next_rank++;
    }

    int prev_rank = rank - 1;
    while (prev_rank > 0
           && rank_n_nodes[prev_rank] == 0) {
      prev_rank--;
    }

    /* printf ("[%d] first last next prev : %d %d %d %d\n", rank, */
    /*         first_rank, last_rank, next_rank, prev_rank); */

    if (rank == first_rank && rank_n_nodes[rank] > 0) {
      PDM_morton_code_t root_DFD;

      root_DFD.L = max_morton_level;
      root_DFD.X[0] = 0;
      root_DFD.X[1] = 0;
      root_DFD.X[2] = 0;

      PDM_morton_code_t FINA;

      PDM_morton_nearest_common_ancestor (root_DFD,
                                          L2->codes[0],
                                          &FINA);

      PDM_morton_code_t child[8];
      PDM_morton_get_children(dim,
                              FINA,
                              child);

      _octants_push_front (L2,
                           child[0],
                           0,
                           0);
    }


    if (rank == last_rank && rank_n_nodes[rank] > 0) {
      PDM_morton_code_t root_DLD;

      root_DLD.L = max_morton_level;
      root_DLD.X[0] = (1u << max_morton_level) - 1u;
      root_DLD.X[1] = (1u << max_morton_level) - 1u;
      root_DLD.X[2] = (1u << max_morton_level) - 1u;

      PDM_morton_code_t FINA;
      PDM_morton_nearest_common_ancestor (root_DLD,
                                          L2->codes[L2->n_nodes -1],
                                          &FINA);

      PDM_morton_code_t child[8];
      PDM_morton_get_children(dim,
                              FINA,
                              child);

      _octants_push_back (L2,
                          child[7],
                          0,
                          0);

    }

    unsigned int sbuff[4];
    unsigned int rbuff[4];
    PDM_MPI_Request srequest;
    PDM_MPI_Request rrequest;

    if (rank < last_rank && rank_n_nodes[rank] > 0) {

      PDM_MPI_Irecv ((void *) rbuff,
                     4,
                     PDM_MPI_UNSIGNED,
                     next_rank,
                     0,
                     comm,
                     &rrequest);

    }

    if (rank > first_rank && rank_n_nodes[rank] > 0) {

      assert (L2->n_nodes > 0);
      sbuff[0] = L2->codes[0].L;
      sbuff[1] = L2->codes[0].X[0];
      sbuff[2] = L2->codes[0].X[1];
      sbuff[3] = L2->codes[0].X[2];

      PDM_MPI_Issend ((void *) sbuff,
                      4,
                      PDM_MPI_UNSIGNED,
                      prev_rank,
                      0,
                      comm,
                      &srequest);


    }

    if (rank < last_rank && rank_n_nodes[rank] > 0) {

      PDM_MPI_Wait (&rrequest);
      PDM_morton_code_t code;

      code.L = rbuff[0];
      code.X[0] = rbuff[1];
      code.X[1] = rbuff[2];
      code.X[2] = rbuff[3];

      _octants_push_back (L2,
                          code,
                          0,
                          0);

    }

    if (rank > first_rank && rank_n_nodes[rank] > 0) {

      PDM_MPI_Wait (&srequest);

    }

    for (int i = 0; i < L2->n_nodes - 1; i++) {
      _l_octant_t *A = _complete_region (L2->codes[i], L2->codes[i+1]);

      _octants_push_back (R,
                          L2->codes[i],
                          0,
                          0);

      if (A != NULL) {
        for (int j = 0; j < A->n_nodes; j++) {
          _octants_push_back (R,
                              A->codes[j],
                              0,
                              0);
        }
        _octants_free (A);
      }
    }

    if (rank == last_rank  && rank_n_nodes[rank] > 0) {
      _octants_push_back (R,
                          L2->codes[L2->n_nodes-1],
                          0,
                          0);
    }

    _octants_free (L2);

    free (rank_n_nodes);
  }

  else {
    if (rank == n_ranks - 1) {

      PDM_morton_code_t _code;
      _code.L = 0;
      _code.X[0] = 0;
      _code.X[1] = 0;
      _code.X[2] = 0;

      _octants_push_back (R,
                          _code,
                          0,
                          0);

    }

  }

  /* printf("fin complete_octree\n"); */
  /* fflush(stdout); */

  return R;
}



/**
 *
 * \brief Distribute points
 *
 * \param [in]   id                 Identifier
 *
 */

static void
_distribute_points
(
 int *n_points,
 double **points,
 int **points_icloud,
 PDM_g_num_t **points_gnum,
 PDM_morton_code_t **points_code,
 PDM_morton_code_t *morton_index,
 const PDM_MPI_Comm comm,
 const int dim,
 const PDM_morton_int_t max_level,
 const double *global_extents
 )
{
  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  int _n_points = *n_points;

  double *__points = *points;
  int *__points_icloud = *points_icloud;
  PDM_g_num_t *__points_gnum = *points_gnum;
  PDM_morton_code_t *__points_code = *points_code;

  int *c_rank = malloc (_n_points * sizeof(int));

  for (int i = 0; i < _n_points; i++) {
    size_t _c_rank = PDM_morton_quantile_search((size_t) n_ranks,
                                                __points_code[i],
                                                morton_index);
    c_rank[i] = (int) _c_rank;
  }

  int *send_count = malloc (n_ranks * sizeof (int));
  int *recv_count = malloc (n_ranks * sizeof (int));
  int *send_shift = malloc ((n_ranks + 1) * sizeof (int));
  int *recv_shift = malloc ((n_ranks + 1) * sizeof (int));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  for (int i = 0; i < _n_points; i++) {
    send_count[c_rank[i]] += dim;
  }

  /* Exchange number of coords to send to each process */

  PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT,
                   recv_count, 1, PDM_MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
  }

  /* Build send and receive buffers */

  double *send_coords = malloc (send_shift[n_ranks] * sizeof(double));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++)
    send_count[rank_id] = 0;

  for (int i = 0; i < _n_points; i++) {
    int rank_id = c_rank[i];
    int shift = send_shift[rank_id] + send_count[rank_id];
    for (int j = 0; j < dim; j++)
      send_coords[shift + j] = __points[i*dim + j];
    send_count[rank_id] += dim;
  }

  double *recv_coords = malloc (recv_shift[n_ranks] * sizeof(double));

  /* Exchange coords between processes */

  PDM_MPI_Alltoallv(send_coords, send_count, send_shift, PDM_MPI_DOUBLE,
                    recv_coords, recv_count, recv_shift, PDM_MPI_DOUBLE,
                    comm);

  free(send_coords);

  /* Build send and receive buffers */

  for (int rank_id = 0; rank_id < n_ranks + 1; rank_id++) {
    send_shift[rank_id] = send_shift[rank_id]/dim;
    recv_shift[rank_id] = recv_shift[rank_id]/dim;
  }

  int *send_points_icloud = malloc (send_shift[n_ranks] * sizeof(int));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    recv_count[rank_id] = recv_count[rank_id]/dim;
    send_count[rank_id] = 0;
  }

  for (int i = 0; i < _n_points; i++) {
    int rank_id = c_rank[i];
    int shift = send_shift[rank_id] + send_count[rank_id];
    send_points_icloud[shift] = __points_icloud[i];
    send_count[rank_id] += 1;
  }

  int *recv_points_icloud = malloc (recv_shift[n_ranks] * sizeof(int));

  /* Exchange points_icloud between processes */

  PDM_MPI_Alltoallv(send_points_icloud, send_count, send_shift, PDM_MPI_INT,
                    recv_points_icloud, recv_count, recv_shift, PDM_MPI_INT,
                    comm);

  free(send_points_icloud);


  /* Build send and receive buffers : points_gnum*/

  PDM_g_num_t *send_points_gnum =
    malloc (send_shift[n_ranks] * sizeof(PDM_g_num_t));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  for (int i = 0; i < _n_points; i++) {
    int rank_id = c_rank[i];
    int shift = send_shift[rank_id] + send_count[rank_id];
    send_points_gnum[shift] = __points_gnum[i];
    send_count[rank_id] += 1;
  }

  free (c_rank);

  PDM_g_num_t *recv_points_gnum =
    malloc (recv_shift[n_ranks] * sizeof(PDM_g_num_t));

  /* Exchange points_gnum between processes */

  PDM_MPI_Alltoallv(send_points_gnum, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                    recv_points_gnum, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                    comm);

  free(send_points_gnum);

  _n_points = recv_shift[n_ranks];

  free (send_count);
  free (recv_count);
  free (send_shift);
  free (recv_shift);

  __points = realloc (__points, sizeof(double)  * 3 * _n_points);

  __points_icloud =
    realloc (__points_icloud, sizeof(int) * _n_points);

  __points_gnum =
    realloc (__points_gnum, sizeof(PDM_g_num_t) * _n_points);

  /* Re-encode points */

  __points_code = realloc (__points_code,
                           sizeof(PDM_morton_code_t) * _n_points);

  double d[3];
  double s[3];

  PDM_morton_encode_coords(dim,
                           max_level,
                           global_extents,
                           _n_points,
                           recv_coords,
                           __points_code,
                           d,
                           s);

  int *order = malloc (sizeof(int) * _n_points);

  for (int i = 0; i < _n_points; i++) {
    order[i] = i;
  }

  PDM_morton_local_order(_n_points, __points_code, order);

  for (int i = 0; i < _n_points; i++) {
    __points_icloud[i] = recv_points_icloud[order[i]];
    __points_gnum[i] = recv_points_gnum[order[i]];
    for (int j = 0; j < dim; j++) {
      __points[dim*i+j] = recv_coords[dim*order[i]+j];
    }
  }

  free (recv_points_icloud);
  free (recv_points_gnum);
  free (recv_coords);

  PDM_morton_code_t *_points_code =
    malloc (sizeof(PDM_morton_code_t) * _n_points);

  for (int i = 0; i < _n_points; i++) {
    _points_code[i].L = __points_code[order[i]].L;
    _points_code[i].X[0] = __points_code[order[i]].X[0];
    _points_code[i].X[1] = __points_code[order[i]].X[1];
    _points_code[i].X[2] = __points_code[order[i]].X[2];
  }

  free (__points_code);
  free (order);

  *points_code = _points_code;

  *points = __points;
  *points_icloud = __points_icloud;
  *points_gnum = __points_gnum;

  *n_points = _n_points;
}


/**
 *
 * \brief Partitioning octants into large contiguous blocks. The list of octants
 *        is redistributed
 *
 * \param [in]  octant_list  a list of distributed octants,
 *                           octant_list is not redistributed at the end
 *
 * \return block_octants  A list of distributed blocks
 *
 */

static _l_octant_t *
_block_partition
(
 _l_octant_t *octant_list,
 const PDM_MPI_Comm comm,
 PDM_morton_code_t **G_morton_index
 )
{

  /* Complete region */

  _l_octant_t *T = NULL;

  int max_level = -1;
  int min_level = 32;

  int comm_rank;
  PDM_MPI_Comm_rank(comm, &comm_rank);

  /* printf("\n_block_partition : octant_list %d : d\n", comm_rank); */
  /* if (octant_list!=NULL) { */
  /*   printf("\n_nodes : %d\n", octant_list->n_nodes); */
  /*   for (int i = 0; i < octant_list->n_nodes; i++) { */
  /*     PDM_morton_dump (3, octant_list->codes[i]); */
  /*   } */
  /* } */
  /* else { */
  /*   printf ("octant_list NULL\n"); */
  /* } */

  /* printf("_block_partition :octant_list f\n\n"); */

  if (octant_list->n_nodes > 1 ) {

    T = _complete_region (octant_list->codes[0],
                          octant_list->codes[octant_list->n_nodes-1]);

    if (T !=NULL) {
      for (int i = 0; i < T->n_nodes; i++) {
        max_level = PDM_MAX ((int) T->codes[i].L, max_level);
        min_level = PDM_MIN ((int) T->codes[i].L, min_level);
      }
    }
  }

  /* printf("\n_block_partition : complete_region %d :d\n", comm_rank); */
  /* if (T!=NULL) { */
  /*   printf("\n_nodes : %d\n", T->n_nodes); */
  /*   for (int i = 0; i < T->n_nodes; i++) { */
  /*     PDM_morton_dump (3, T->codes[i]); */
  /*   } */
  /* } */
  /* else { */
  /*   printf ("T NULL\n"); */
  /* } */
  /* printf("_block_partition : complete_region f\n\n"); */

  int max_max_level;
  PDM_MPI_Allreduce(&max_level, &max_max_level, 1,
                    PDM_MPI_INT, PDM_MPI_MAX, comm);

  /* Intialize C with coarse octants */

  _l_octant_t *C = malloc(sizeof(_l_octant_t));

  if (T != NULL) {

    _octants_init (C, T->dim, T->n_nodes);

    for (int i = 0; i < T->n_nodes; i++) {

      if ( (int) T->codes[i].L <= min_level) {
        _octants_push_back (C,
                            T->codes[i],
                            T->n_points[i],
                            T->range[i]);
      }
    }
  }

  else {
    _octants_init (C, octant_list->dim, 1);
  }

  /* Complete octree */

  /* printf("\n_block_partition : before complete_octree %d %d : d\n", comm_rank, C->n_nodes); */
  /* for (int i = 0; i < C->n_nodes; i++) { */
  /*   PDM_morton_dump (3, C->codes[i]); */
  /* } */
  /* printf("_block_partition : before complete_octree f\n\n"); */

  _l_octant_t *G = _complete_octree (C, comm);


  double vol = 0;
  for (int i = 0; i < G->n_nodes; i++) {
    double _side = 1. / (double) (1 << G->codes[i].L);
    vol += _side*_side*_side;
    G->range[i+1] =
      G->range[i] +
      G->n_points[i];
  }
  double total_vol;
  PDM_MPI_Allreduce(&vol, &total_vol, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  if ( (PDM_ABS(total_vol - 1.)>= 1e-15)) {
    printf("Erreur volume different de 1 apres complete_octree : %12.5e\n", total_vol);
    for (int i = 0; i < G->n_nodes; i++) {
      PDM_morton_dump (3, G->codes[i]);
    }
  }

  assert (PDM_ABS(total_vol - 1.) < 1e-15);


  /* printf("\n_block_partition : after complete_octree %d %d : d\n", comm_rank, C->n_nodes); */
  /* for (int i = 0; i < G->n_nodes; i++) { */
  /*   PDM_morton_dump (3, G->codes[i]); */
  /* } */
  /* printf("_block_partition : after complete_octree f\n\n"); */

  /* PDM_MPI_Barrier (comm); */
  /* exit(1); */
  _octants_free (C);

  if (T != NULL) {
    _octants_free (T);
  }

  /*
   * Compute weight
   */

  /* - exchange codes to ranks (weight per rank)*/

  int n_ranks;
  PDM_MPI_Comm_size(comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank(comm, &rank);

  PDM_morton_int_t *code_buff = malloc (sizeof(PDM_morton_int_t) * (G->dim + 1));
  PDM_morton_int_t *rank_buff = malloc (sizeof(PDM_morton_int_t) * n_ranks * (G->dim + 1));
  int *n_nodes_rank = malloc (sizeof(int) * n_ranks);

  for (int i = 0; i < G->dim + 1; i++) {
    code_buff[i] = 0;
  }

  if ( G->n_nodes > 0) {
    code_buff[0] = G->codes[0].L;
    for (int i = 0; i < G->dim; i++) {
      code_buff[i+1] =  G->codes[0].X[i];
    }
  }

  PDM_MPI_Allgather (&(G->n_nodes), 1, PDM_MPI_INT,
                     n_nodes_rank,  1, PDM_MPI_INT,
                     comm);

  PDM_MPI_Allgather (code_buff, G->dim + 1, PDM_MPI_UNSIGNED,
                     rank_buff, G->dim + 1, PDM_MPI_UNSIGNED,
                     comm);

  int n_active_ranks = 0;
  for (int i = 0; i < n_ranks; i++) {
    if (n_nodes_rank[i] > 0) {
      n_active_ranks++;
    }
  }

  //assert (n_active_ranks > 0);

  PDM_morton_code_t *rank_codes = malloc (sizeof(PDM_morton_code_t) * n_active_ranks);
  int *active_ranks = malloc (sizeof(int) * n_active_ranks);

  n_active_ranks = 0;
  for (int i = 0; i < n_ranks; i++) {
    if (n_nodes_rank[i] > 0) {
      active_ranks[n_active_ranks] = i;
      rank_codes[n_active_ranks].L = rank_buff[(G->dim + 1) * i];
      for (int j = 0; j < G->dim; j++) {
        rank_codes[n_active_ranks].X[j] = rank_buff[(G->dim + 1) * i + j + 1];
      }
      n_active_ranks++;
    }
  }

  free (code_buff);
  free (rank_buff);

  int *send_count = malloc(sizeof(int) * n_ranks);
  int *send_shift = malloc(sizeof(int) * (n_ranks+1));

  int *recv_count = malloc(sizeof(int) * n_ranks);
  int *recv_shift = malloc(sizeof(int) * (n_ranks+1));

  int irank = 0;
  for (int i = 0; i < n_ranks; i++) {
    send_count[i] = 0;
  }

  /* printf("rank codes deb\n"); */
  /* for (int i = 0; i < n_active_ranks; i++) { */
  /*   PDM_morton_dump(3, rank_codes[i]); */
  /* } */
  /* printf("rank codes fin\n"); */

  for (int i = 0; i < octant_list->n_nodes; i++) {

    if (irank < (n_active_ranks - 1)) {
      if (PDM_morton_a_ge_b (octant_list->codes[i], rank_codes[irank+1])) {

        irank += 1 + PDM_morton_binary_search(n_active_ranks - (irank + 1),
                                              octant_list->codes[i],
                                              rank_codes + irank + 1);
      }
    }
    send_count[active_ranks[irank]] += octant_list->dim + 2;
  }

  /* Exchange number of coords to send to each process */

  PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT,
                   recv_count, 1, PDM_MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
  }

  /* Build send and receive buffers */

  PDM_morton_int_t *send_codes =
    malloc (send_shift[n_ranks] * sizeof(PDM_morton_int_t));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  irank = 0;
  for (int i = 0; i < octant_list->n_nodes; i++) {

    if (irank < (n_active_ranks - 1)) {
      if (PDM_morton_a_ge_b (octant_list->codes[i], rank_codes[irank+1])) {

        irank += 1 + PDM_morton_binary_search(n_active_ranks - (irank + 1),
                                              octant_list->codes[i],
                                              rank_codes + irank + 1);
      }
    }

    int shift = send_shift[active_ranks[irank]] + send_count[active_ranks[irank]];

    assert(octant_list->n_points[i] >= 0);

    send_codes[shift++] = octant_list->codes[i].L;

    for (int j = 0; j < octant_list->dim; j++) {
      send_codes[shift++] = octant_list->codes[i].X[j];
    }

    send_codes[shift++] = (PDM_morton_int_t) octant_list->n_points[i];

    send_count[active_ranks[irank]] += octant_list->dim + 2;
  }

  free (rank_codes);

  PDM_morton_int_t *recv_codes =
    malloc (recv_shift[n_ranks] * sizeof(PDM_morton_int_t));

  /* - exchange codes between processes */

  PDM_MPI_Alltoallv(send_codes, send_count, send_shift, PDM_MPI_UNSIGNED,
                    recv_codes, recv_count, recv_shift, PDM_MPI_UNSIGNED,
                    comm);


  free (send_codes);
  free (send_count);
  free (send_shift);
  free (recv_count);

  const int _stride = octant_list->dim + 2;
  const int n_recv_codes = recv_shift[n_ranks] / _stride;

  free (recv_shift);

  int *weight = malloc (sizeof(int) * G->n_nodes);

  for (int i = 0; i < G->n_nodes; i++) {
    weight[i] = 0;
  }

  /* - compute weight of each cell */

  for (int i = 0; i < n_recv_codes; i++) {

    PDM_morton_code_t code;

    code.L = recv_codes[i*_stride];

    for (int j = 0; j < octant_list->dim; j++) {
      code.X[j] = recv_codes[i*_stride+j+1];
    }

    int G_node =  PDM_morton_binary_search(G->n_nodes,
                                           code,
                                           G->codes);

    weight[G_node] +=  recv_codes[i*_stride + 1 + octant_list->dim];
  }

  free (recv_codes);

  /*
   * Load balancing G from weight
   */

  int *order = malloc (sizeof(int) * G->n_nodes);

  for (int i = 0; i <  G->n_nodes; i++) {
    order[i] = i;
  }

  *G_morton_index = malloc(sizeof(PDM_morton_code_t) * (n_ranks + 1));
  PDM_morton_code_t *_G_morton_index = *G_morton_index;

  PDM_morton_ordered_build_rank_index (octant_list->dim,
                                       max_max_level,
                                       G->n_nodes,
                                       G->codes,
                                       weight,
                                       _G_morton_index,
                                       comm);

  free (order);
  free (weight);

  /* printf("\nblock_partition avant d\n"); */
  /* for (int i = 0; i <  G->n_nodes; i++) { */
  /*   PDM_morton_dump (3, G->codes[i]); */
  /* } */
  /* printf("block_partition avant f\n\n"); */

  _distribute_octants (G, _G_morton_index, comm);

  /*
   * Redistribute octant list from coarse load balancing
   */

  _distribute_octants (octant_list, _G_morton_index, comm);

  free (active_ranks);
  free (n_nodes_rank);
  return G;

}





static void
_compute_connected_parts
(
 _octree_t         *octree,
 _neighbours_tmp_t *neighbours
 )
{
  const int n_direction = (int) PDM_N_DIRECTION;

  int s_connected = 3;
  octree->connected_idx = malloc (sizeof(int) * s_connected);
  octree->connected_idx[0] = 0;

  int *visited = malloc (sizeof(int) * octree->octants->n_nodes);
  for (int i = 0; i < octree->octants->n_nodes; i++)
    visited[i] = 0;

  int *stack = malloc (sizeof(int) * octree->octants->n_nodes);
  int pos_stack = 0;

  int max = 0;
  while (max < octree->octants->n_nodes) {

    stack[pos_stack++] = max;
    visited[max] = 1;

    while (pos_stack > 0) {
      int i_node = stack[--pos_stack];

      for (int j = 0; j < n_direction; j++) {
        for (int k = 0; k < neighbours[i_node].n_neighbour[j]; k++) {
          int i_ngb = neighbours[i_node].neighbours[j][k];
          if (i_ngb >= 0) {
            if (!visited[i_ngb]) {
              stack[pos_stack++] = i_ngb;
              visited[i_ngb] = 1;
              max = PDM_MAX (max, i_ngb);
            }
          }
        }
      }

    }

    max++;
    if (s_connected <= octree->n_connected) {
      s_connected *= 2;
      octree->connected_idx = realloc (octree->connected_idx, sizeof(int) * s_connected);
    }
    octree->connected_idx[++octree->n_connected] = max;
  }
  free (visited);
  free (stack);

  octree->connected_idx = realloc (octree->connected_idx,
				   sizeof(int) * (octree->n_connected+1));
}









/**
 *
 * \brief Compute neighbours
 *
 * \param [inout]   octree
 * \param [in]      b_t_elapsed
 * \param [in]      b_t_cpu
 * \param [in]      b_t_cpu_u
 * \param [in]      b_t_cpu_s
 *
 */

static void
_compute_neighbours
(
 _octree_t *octree,
 double   b_t_elapsed,
 double   b_t_cpu,
 double   b_t_cpu_u,
 double   b_t_cpu_s
 )
{
  double   e_t_elapsed;
  double   e_t_cpu;
  double   e_t_cpu_u;
  double   e_t_cpu_s;

  const int n_direction = (int) PDM_N_DIRECTION;

  int n_ranks;
  PDM_MPI_Comm_size (octree->comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  _neighbours_tmp_t *neighbours_tmp = malloc (sizeof(_neighbours_tmp_t) * octree->octants->n_nodes);
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j =  0; j < n_direction; j++) {
      neighbours_tmp[i].n_neighbour[j] = 0;
      neighbours_tmp[i].s_neighbour[j] = 1;
      neighbours_tmp[i].neighbours[j] = malloc (sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
      neighbours_tmp[i].neighbours[j][0] = 0;
    }
  }

  /* Boucle sur les noeuds : */
  size_t start_intersect, end_intersect;
  size_t start_intersect_quantile, end_intersect_quantile;

  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 1; j < n_direction; j+=2) {
      PDM_morton_code_t *neighbour_code =
        _neighbour (octree->octants->codes[i], (PDM_para_octree_direction_t)j);
      PDM_para_octree_direction_t inv_j =
        _inv_direction((PDM_para_octree_direction_t) j);

      if (neighbour_code != NULL) {

        if (octree->rank_octants_index != NULL) {

          PDM_morton_quantile_intersect (n_ranks,
                                         *neighbour_code,
                                         octree->rank_octants_index,
                                         &start_intersect_quantile,
                                         &end_intersect_quantile);
        }
        else {
          start_intersect_quantile = 0;
          end_intersect_quantile = 1;
        }

        for (int neighbour_rank = start_intersect_quantile;
             neighbour_rank < (int) end_intersect_quantile; neighbour_rank++) {

          if (neighbour_rank == rank) {

            PDM_morton_list_intersect (octree->octants->n_nodes - (i+1),
                                       *neighbour_code,
                                       octree->octants->codes + i + 1,
                                       &start_intersect,
                                       &end_intersect);

            for (int k = start_intersect; k < (int) end_intersect; k++) {
              int idx = k + i + 1;
              PDM_morton_code_t *neighbour_neighbour_code =
                _neighbour (octree->octants->codes[idx], inv_j);

              assert (neighbour_neighbour_code != NULL);

              if (PDM_morton_ancestor_is (octree->octants->codes[i], *neighbour_neighbour_code) ||
                  PDM_morton_ancestor_is (*neighbour_neighbour_code, octree->octants->codes[i])) {

                if (neighbours_tmp[i].n_neighbour[j] >= neighbours_tmp[i].s_neighbour[j]) {
                  neighbours_tmp[i].s_neighbour[j] *= 2;
                  neighbours_tmp[i].neighbours[j] =
                    realloc (neighbours_tmp[i].neighbours[j],
                             sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
                }
                neighbours_tmp[i].neighbours[j][neighbours_tmp[i].n_neighbour[j]++] = idx;

                if (neighbours_tmp[idx].n_neighbour[inv_j] >= neighbours_tmp[idx].s_neighbour[inv_j]) {
                  neighbours_tmp[idx].s_neighbour[inv_j] *= 2;
                  neighbours_tmp[idx].neighbours[inv_j] =
                    realloc (neighbours_tmp[idx].neighbours[inv_j],
                             sizeof(int) * neighbours_tmp[idx].s_neighbour[inv_j]);
                }
                neighbours_tmp[idx].neighbours[inv_j][neighbours_tmp[idx].n_neighbour[inv_j]++] = i;
              }

              free (neighbour_neighbour_code);
            }

          }

          else {
            if (neighbours_tmp[i].n_neighbour[j] >= neighbours_tmp[i].s_neighbour[j]) {
              neighbours_tmp[i].s_neighbour[j] *= 2;
              neighbours_tmp[i].neighbours[j] = realloc (neighbours_tmp[i].neighbours[j],
                                                         sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
            }
            neighbours_tmp[i].neighbours[j][neighbours_tmp[i].n_neighbour[j]++] = - (neighbour_rank + 1);
          }
        }
        free (neighbour_code);
      }
    }
  }

  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP1] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP1]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP1]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP1]   += e_t_cpu_s - b_t_cpu_s;

  PDM_timer_resume(octree->timer);

  PDM_timer_hang_on(octree->timer);
  b_t_elapsed = PDM_timer_elapsed(octree->timer);
  b_t_cpu     = PDM_timer_cpu(octree->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);
  PDM_timer_resume(octree->timer);

  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j+=2) {
      PDM_morton_code_t *neighbour_code =
        _neighbour (octree->octants->codes[i], (PDM_para_octree_direction_t) j);

      if (neighbour_code != NULL) {

        if (octree->rank_octants_index != NULL) {

          PDM_morton_quantile_intersect (n_ranks,
                                         *neighbour_code,
                                         octree->rank_octants_index,
                                         &start_intersect_quantile,
                                         &end_intersect_quantile);
        }
        else {
          start_intersect_quantile = 0;
          end_intersect_quantile = 1;
        }

        for (int neighbour_rank = start_intersect_quantile;
             neighbour_rank < (int) end_intersect_quantile; neighbour_rank++) {

          if (neighbour_rank != rank) {
            if (neighbours_tmp[i].n_neighbour[j] >= neighbours_tmp[i].s_neighbour[j]) {
              neighbours_tmp[i].s_neighbour[j] *= 2;
              neighbours_tmp[i].neighbours[j] = realloc (neighbours_tmp[i].neighbours[j],
                                                         sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
            }
            neighbours_tmp[i].neighbours[j][neighbours_tmp[i].n_neighbour[j]++] = - (neighbour_rank + 1);
          }
        }
        free (neighbour_code);
      }
    }
  }


  /*************************************************************************
   *
   * Compute connected parts
   *
   *************************************************************************/
  _compute_connected_parts (octree, neighbours_tmp);


  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP2] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP2]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP2]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP2]   += e_t_cpu_s - b_t_cpu_s;

  PDM_timer_resume(octree->timer);

  PDM_timer_hang_on(octree->timer);
  b_t_elapsed = PDM_timer_elapsed(octree->timer);
  b_t_cpu     = PDM_timer_cpu(octree->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);
  PDM_timer_resume(octree->timer);

  /*************************************************************************
   *
   * Build parallel partition boundary
   *
   *************************************************************************/

  int FALSE_NEIGHBOUR = -1;
  if (n_ranks > 1) {
    const int n_quantile =  n_ranks * n_direction;

    int *neighbour_rank_n = malloc (sizeof(int) * n_quantile);
    int *neighbour_rank_idx = malloc (sizeof(int) * (n_quantile + 1));

    for (int i = 0; i < n_quantile; i++) {
      neighbour_rank_n[i] = 0;
    }

    /* Premiere boucle pour compter */

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (int j = 0; j < n_direction; j++) {
        for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
          if (neighbours_tmp[i].neighbours[j][k] < 0) {
            neighbour_rank_n[-(neighbours_tmp[i].neighbours[j][k] + 1)*n_direction +j]++;
          }
        }
      }
    }

    int max_node_dir = -1;
    neighbour_rank_idx[0] = 0;
    for (int i = 0; i < n_quantile; i++) {
      neighbour_rank_idx[i+1] = neighbour_rank_idx[i] + neighbour_rank_n[i];
      max_node_dir = PDM_MAX (max_node_dir, neighbour_rank_n[i]);
      neighbour_rank_n[i] = 0;
    }

    /* Allocation */

    int *neighbour_rank_node_id = malloc (sizeof(int) * neighbour_rank_idx[n_quantile]);
    PDM_morton_code_t *neighbour_rank_code = malloc (sizeof(PDM_morton_code_t) *
                                                     neighbour_rank_idx[n_quantile]);
    int *neighbour_rank_node_k = malloc (sizeof(int) * neighbour_rank_idx[n_quantile]);
    int *neighbour_rank_node_part = malloc (sizeof(int) * neighbour_rank_idx[n_quantile]);//

    /* Deuxieme boucle pour stocker avec tri suivant la direction */

    for (int i_part = 0; i_part < octree->n_connected; i_part++) {
      for (int i = octree->connected_idx[i_part]; i < octree->connected_idx[i_part+1]; i++) {
        for (int j = 0; j < n_direction; j++) {
          for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
            if (neighbours_tmp[i].neighbours[j][k] < 0) {
              int index = -(neighbours_tmp[i].neighbours[j][k] + 1)*n_direction +j;
              int index2 = neighbour_rank_idx[index] + neighbour_rank_n[index];

              neighbour_rank_node_id[index2] = i;
              PDM_morton_copy (octree->octants->codes[i],
                               neighbour_rank_code + index2);
              neighbour_rank_node_k[index2] = k;
              neighbour_rank_node_part[index2] = i_part;

              neighbour_rank_n[index]++;
            }
          }
        }
      }
    }


    /* Tri des codes pour chaque direction de chaque rang */
    int *order = malloc (sizeof(int) * max_node_dir);
    int *tmp_node_id = malloc (sizeof(int) * max_node_dir);
    PDM_morton_code_t *tmp_code = malloc (sizeof(PDM_morton_code_t) * max_node_dir);
    int *tmp_node_k = malloc (sizeof(int) * max_node_dir);
    int *tmp_node_part = malloc (sizeof(int) * max_node_dir);

    for (int i = 0; i < n_ranks; i++) {
      for (int j = 0; j < n_direction; j++) {
        PDM_morton_local_order (neighbour_rank_n[n_direction * i + j],
                                neighbour_rank_code + neighbour_rank_idx[n_direction * i + j],
                                order);
        int idx1 = 0;
        for (int k = neighbour_rank_idx[n_direction * i + j];
             k < neighbour_rank_idx[n_direction * i + j + 1];
             k++) {
          PDM_morton_copy (neighbour_rank_code[k], tmp_code + idx1);
          tmp_node_id[idx1]   = neighbour_rank_node_id[k];
          tmp_node_k[idx1]    = neighbour_rank_node_k[k];
          tmp_node_part[idx1] = neighbour_rank_node_part[k];
          idx1 += 1;
        }

        idx1 = 0;
        for (int k = neighbour_rank_idx[n_direction * i + j];
             k < neighbour_rank_idx[n_direction * i + j + 1];
             k++) {
          PDM_morton_copy (tmp_code[order[idx1]], neighbour_rank_code + k );
          neighbour_rank_node_id[k]   = tmp_node_id[order[idx1]];
          neighbour_rank_node_k[k]    = tmp_node_k[order[idx1]];
          neighbour_rank_node_part[k] = tmp_node_part[order[idx1]];
          idx1 += 1;
        }
      }
    }

    free (tmp_code);
    free (order);
    free (tmp_node_id);
    free (tmp_node_k);
    free (tmp_node_part);


    /* Envoi / reception (Les donnees recues sont triees) */

    int *recv_neighbour_rank_n = malloc (sizeof(int) * n_quantile);

    for (int i = 0; i < n_quantile; i++) {
      recv_neighbour_rank_n[i] = 0;
    }

    PDM_MPI_Request *recv_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);
    PDM_MPI_Request *send_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);

    int *used_ranks = malloc (sizeof(int) * n_ranks);

    PDM_MPI_Alltoall (neighbour_rank_n, n_direction, PDM_MPI_INT,
                      recv_neighbour_rank_n, n_direction, PDM_MPI_INT,
                      octree->comm);

    int *recv_neighbour_rank_idx = malloc (sizeof(int) * (n_direction * n_ranks + 1));
    recv_neighbour_rank_idx[0] = 0;


    for (int i = 0; i <  n_direction * n_ranks; i++)
      recv_neighbour_rank_idx[i+1] = recv_neighbour_rank_idx[i] + recv_neighbour_rank_n[i];



    int *recv_neighbour_rank_node_id   = malloc (sizeof(int) * recv_neighbour_rank_idx[n_quantile]);
    int *recv_neighbour_rank_node_part = malloc (sizeof(int) * recv_neighbour_rank_idx[n_quantile]);
    PDM_morton_code_t *recv_neighbour_rank_code =
      malloc (sizeof(PDM_morton_code_t) * recv_neighbour_rank_idx[n_quantile]);


    unsigned int *_neighbour_rank_code =
      malloc (sizeof(unsigned int) * 4 * neighbour_rank_idx[n_quantile]);
    unsigned int *_recv_neighbour_rank_code =
      malloc (sizeof(unsigned int) * 4 * recv_neighbour_rank_idx[n_quantile]);

    int idx = 0;
    for (int i = 0; i < neighbour_rank_idx[n_quantile]; i++) {
      _neighbour_rank_code[idx++] = neighbour_rank_code[i].L;
      for (int j = 0; j < 3; j++) {
        _neighbour_rank_code[idx++] = neighbour_rank_code[i].X[j];
      }
    }

    int *rank_neighbour_rank_n = malloc (sizeof(int) * n_ranks);
    int *rank_neighbour_rank_idx = malloc (sizeof(int) * (n_ranks + 1));
    int *rank_recv_neighbour_rank_n = malloc (sizeof(int) * n_ranks);
    int *rank_recv_neighbour_rank_idx = malloc (sizeof(int) * (n_ranks + 1));

    rank_neighbour_rank_idx[0] = 0;
    rank_recv_neighbour_rank_idx[0] = 0;

    for (int i = 0; i < n_ranks; i++) {
      rank_neighbour_rank_n[i] = 0;
      rank_recv_neighbour_rank_n[i] = 0;
    }

    for (int i = 0; i < n_ranks; i++) {
      for (int j = 0; j < n_direction; j++) {
        rank_neighbour_rank_n[i] += neighbour_rank_n[i*n_direction+j];
        rank_recv_neighbour_rank_n[i] += recv_neighbour_rank_n[i*n_direction+j];
      }
      rank_neighbour_rank_idx[i+1] = rank_neighbour_rank_n[i] + rank_neighbour_rank_idx[i];
      rank_recv_neighbour_rank_idx[i+1] = rank_recv_neighbour_rank_n[i] + rank_recv_neighbour_rank_idx[i];
    }

    PDM_MPI_Alltoallv (neighbour_rank_node_id,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_INT,
                       recv_neighbour_rank_node_id,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_INT,
                       octree->comm);

    PDM_MPI_Alltoallv (neighbour_rank_node_part,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_INT,
                       recv_neighbour_rank_node_part,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_INT,
                       octree->comm);

    for (int i = 0; i < n_ranks; i++) {
      rank_neighbour_rank_n[i] *= 4;
      rank_recv_neighbour_rank_n[i] *= 4;
      rank_neighbour_rank_idx[i+1] *= 4;
      rank_recv_neighbour_rank_idx[i+1] *= 4;
    }


    PDM_MPI_Alltoallv (_neighbour_rank_code,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_UNSIGNED,
                       _recv_neighbour_rank_code,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_UNSIGNED,
                       octree->comm);


    free (_neighbour_rank_code);

    free (rank_neighbour_rank_n);
    free (rank_neighbour_rank_idx);
    free (rank_recv_neighbour_rank_n);
    free (rank_recv_neighbour_rank_idx);

    idx = 0;
    for (int i = 0; i < recv_neighbour_rank_idx[n_quantile]; i++) {
      recv_neighbour_rank_code[i].L = _recv_neighbour_rank_code[idx++];
      for (int j = 0; j < 3; j++) {
        recv_neighbour_rank_code[i].X[j] = _recv_neighbour_rank_code[idx++];
      }
    }
    free (_recv_neighbour_rank_code);

    free (recv_request);
    free (send_request);

    free (used_ranks);


    octree->n_part_boundary_elt = neighbour_rank_idx[n_quantile];
    octree->part_boundary_elt_idx = malloc (sizeof(int) * (octree->n_part_boundary_elt + 1));

    int s_part_boundary_elt = 2 * 3 * neighbour_rank_idx[n_quantile];
    octree->part_boundary_elt = malloc (sizeof(int) * s_part_boundary_elt);

    int n_part_boundary_elt = 0;

    idx = 0;
    int idx_part_boundary_elt = 0;
    for (int i = 0; i <= octree->n_part_boundary_elt; i++)
      octree->part_boundary_elt_idx[i] = 0;


    FALSE_NEIGHBOUR = -(octree->n_part_boundary_elt + 1);

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (int j = 0; j < n_direction; j++) {
        for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
          if (neighbours_tmp[i].neighbours[j][k] < 0)
            neighbours_tmp[i].neighbours[j][k] = FALSE_NEIGHBOUR;
        }
      }
    }



    for (int i = 0; i < n_ranks; i++ ) {
      for (PDM_para_octree_direction_t j = PDM_BOTTOM; j < PDM_N_DIRECTION; j++) {
        PDM_para_octree_direction_t inv_j = _inv_direction(j);

        int idx_recv = n_direction * i + inv_j;

        int idx_candidate = recv_neighbour_rank_idx[idx_recv];
        int n_candidate = recv_neighbour_rank_idx[idx_recv+1] - idx_candidate;

        if (n_candidate > 0) {

          for (int k = neighbour_rank_idx[i * n_direction + j];
               k < neighbour_rank_idx[i * n_direction + j + 1]; k++) {
            PDM_morton_code_t *neighbour_code = _neighbour (neighbour_rank_code[k], j);

            PDM_morton_list_intersect (n_candidate,
                                       *neighbour_code,
                                       recv_neighbour_rank_code + idx_candidate,
                                       &start_intersect,
                                       &end_intersect);

            int n_intersect_neighbours = 0;

            if (end_intersect > start_intersect) {
              for (int k1 = start_intersect; k1 < (int) end_intersect; k1++) {
                int k2 = idx_candidate + k1;

                PDM_morton_code_t *neighbour_neighbour_code =
                  _neighbour (recv_neighbour_rank_code[k2], inv_j);

                assert (neighbour_neighbour_code != NULL);

                if (PDM_morton_ancestor_is (neighbour_rank_code[k], *neighbour_neighbour_code) ||
                    PDM_morton_ancestor_is (*neighbour_neighbour_code, neighbour_rank_code[k])) {
                  n_intersect_neighbours++;

                  if ((s_part_boundary_elt - idx_part_boundary_elt) <= 3) {
                    s_part_boundary_elt *= 2;
                    octree->part_boundary_elt = realloc (octree->part_boundary_elt,
                                                         sizeof(int) * s_part_boundary_elt);
                  }
                  octree->part_boundary_elt_idx[n_part_boundary_elt+1]++;
                  octree->part_boundary_elt[idx_part_boundary_elt++] = i; // rank
                  octree->part_boundary_elt[idx_part_boundary_elt++] = recv_neighbour_rank_node_id[k2]; // neighbour's local number in rank i
                  octree->part_boundary_elt[idx_part_boundary_elt++] = recv_neighbour_rank_node_part[k2]; // neighbour's part number in rank i
                }

                free (neighbour_neighbour_code);
              }
            }

            int k3 = neighbour_rank_node_k[k];
            int i2 = neighbour_rank_node_id[k];

            if (n_intersect_neighbours > 0) {
              neighbours_tmp[i2].neighbours[j][k3] = -(n_part_boundary_elt+1);

              assert (neighbours_tmp[i2].neighbours[j][k3] != FALSE_NEIGHBOUR);

              n_part_boundary_elt++;
            }
            else {
              neighbours_tmp[i2].neighbours[j][k3] = FALSE_NEIGHBOUR;
            }

            free (neighbour_code);
          }
        }
      }
    }


    free (neighbour_rank_n);
    free (neighbour_rank_idx);
    free (neighbour_rank_node_id);
    free (neighbour_rank_node_k);
    free (neighbour_rank_node_part);
    free (neighbour_rank_code);

    free (recv_neighbour_rank_n);
    free (recv_neighbour_rank_idx);

    free (recv_neighbour_rank_node_id);
    free (recv_neighbour_rank_node_part);
    free (recv_neighbour_rank_code);

    octree->n_part_boundary_elt = n_part_boundary_elt;
  }

  for (int i = 0; i < octree->n_part_boundary_elt; i++) {
    octree->part_boundary_elt_idx[i+1] += octree->part_boundary_elt_idx[i];
  }


  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_DISTANT_NEIGHBOURS] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_DISTANT_NEIGHBOURS]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_DISTANT_NEIGHBOURS]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_DISTANT_NEIGHBOURS]   += e_t_cpu_s - b_t_cpu_s;

  PDM_timer_resume(octree->timer);

  PDM_timer_hang_on(octree->timer);
  b_t_elapsed = PDM_timer_elapsed(octree->timer);
  b_t_cpu     = PDM_timer_cpu(octree->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);
  PDM_timer_resume(octree->timer);



  /*************************************************************************
   *
   * Copy temporary neighbours in the neighbour structure
   *
   *************************************************************************/
  octree->octants->neighbour_idx =
    malloc(sizeof(int) * (n_direction * octree->octants->n_nodes + 1));

  int idx = 0;
  octree->octants->neighbour_idx[0] = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      octree->octants->neighbour_idx[idx+1] =
        octree->octants->neighbour_idx[idx] + neighbours_tmp[i].n_neighbour[j];

      /* account for false distant neighbours */
      for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
        if (neighbours_tmp[i].neighbours[j][k] == FALSE_NEIGHBOUR)
          octree->octants->neighbour_idx[idx+1]--;
      }

      idx += 1;
    }
  }

  octree->octants->neighbours =
    malloc(sizeof(int) *
           octree->octants->neighbour_idx[n_direction * octree->octants->n_nodes]);

  idx = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {

        if (neighbours_tmp[i].neighbours[j][k] != FALSE_NEIGHBOUR)
          octree->octants->neighbours[idx++] = neighbours_tmp[i].neighbours[j][k];

      }
    }
  }

  /* Free temporary arrays */
  /* printf("sortie 2 neighbours_tmp debut\n"); */
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      if (neighbours_tmp[i].neighbours[j] != NULL) {
        free (neighbours_tmp[i].neighbours[j]);
      }
    }
  }
  /* printf("sortie 2 neighbours_tmp fin\n"); */

  free (neighbours_tmp);




  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP3] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP3]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP3]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP3]   += e_t_cpu_s - b_t_cpu_s;


  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS] += octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS] += octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS] += octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS] += octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP3];

  PDM_timer_resume(octree->timer);
}




static void
_finalize_neighbours
(
 _octree_t          *octree,
 _neighbours_tmp_t **ngb_octree,
 double              b_t_elapsed,
 double              b_t_cpu,
 double              b_t_cpu_u,
 double              b_t_cpu_s
 )
{
  double   e_t_elapsed;
  double   e_t_cpu;
  double   e_t_cpu_u;
  double   e_t_cpu_s;

  const int n_direction = (int) PDM_N_DIRECTION;

  int n_ranks;
  PDM_MPI_Comm_size (octree->comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  _neighbours_tmp_t *neighbours_tmp = *ngb_octree;

  /*************************************************************************
   *
   * Compute connected parts
   *
   *************************************************************************/
  _compute_connected_parts (octree,
			    neighbours_tmp);

  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP1] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP1]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP1]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP1]   += e_t_cpu_s - b_t_cpu_s;

  PDM_timer_resume(octree->timer);


  PDM_timer_hang_on(octree->timer);
  b_t_elapsed = PDM_timer_elapsed(octree->timer);
  b_t_cpu     = PDM_timer_cpu(octree->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);
  PDM_timer_resume(octree->timer);


  /*************************************************************************
   *
   * Build parallel partition boundary
   *
   *************************************************************************/
  int FALSE_NEIGHBOUR = -1;
  if (n_ranks > 1) {
    const int n_quantile = n_ranks * n_direction;

    int *neighbour_rank_n   = malloc (sizeof(int) * n_quantile);
    int *neighbour_rank_idx = malloc (sizeof(int) * (n_quantile + 1));

    for (int i = 0; i < n_quantile; i++) {
      neighbour_rank_n[i] = 0;
    }

    /* Premiere boucle pour compter */
    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < PDM_N_DIRECTION; dir++) {

        PDM_morton_code_t *neighbour_code = _neighbour (octree->octants->codes[i], dir);
        if (neighbour_code == NULL) {
          continue;
        }

        // single neighbour with inferior or equal level?...

        size_t start, end;
        PDM_morton_quantile_intersect (n_ranks,
                                       *neighbour_code,
                                       octree->rank_octants_index,
                                       &start,
                                       &end);

        for (int neighbour_rank = start; neighbour_rank < (int) end; neighbour_rank++) {

          if (neighbour_rank == rank) {
            continue;
          }

          if (neighbours_tmp[i].n_neighbour[dir] >= neighbours_tmp[i].s_neighbour[dir]) {
            neighbours_tmp[i].s_neighbour[dir] *= 2;
            neighbours_tmp[i].neighbours[dir] = realloc (neighbours_tmp[i].neighbours[dir],
                                                         sizeof(int) * neighbours_tmp[i].s_neighbour[dir]);
          }
          neighbours_tmp[i].neighbours[dir][neighbours_tmp[i].n_neighbour[dir]++] = - (neighbour_rank + 1);
          neighbour_rank_n[neighbour_rank*n_direction + dir]++;
        }
        free (neighbour_code);
      }
    }

    int max_node_dir = -1;
    neighbour_rank_idx[0] = 0;
    for (int i = 0; i < n_quantile; i++) {
      neighbour_rank_idx[i+1] = neighbour_rank_idx[i] + neighbour_rank_n[i];
      max_node_dir = PDM_MAX (max_node_dir, neighbour_rank_n[i]);
      neighbour_rank_n[i] = 0;
    }



    /* Allocation */
    int *neighbour_rank_node_id = malloc (sizeof(int) * neighbour_rank_idx[n_quantile]);
    PDM_morton_code_t *neighbour_rank_code = malloc (sizeof(PDM_morton_code_t) *
                                                     neighbour_rank_idx[n_quantile]);
    int *neighbour_rank_node_k = malloc (sizeof(int) * neighbour_rank_idx[n_quantile]);
    int *neighbour_rank_node_part = malloc (sizeof(int) * neighbour_rank_idx[n_quantile]);


    /* Deuxieme boucle pour stocker avec tri suivant la direction */
    for (int i_part = 0; i_part < octree->n_connected; i_part++) {
      for (int i = octree->connected_idx[i_part]; i < octree->connected_idx[i_part+1]; i++) {
        for (int j = 0; j < n_direction; j++) {
          for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
            if (neighbours_tmp[i].neighbours[j][k] < 0) {
              int index = -(neighbours_tmp[i].neighbours[j][k] + 1)*n_direction +j;
              int index2 = neighbour_rank_idx[index] + neighbour_rank_n[index];

              neighbour_rank_node_id[index2] = i;
              PDM_morton_copy (octree->octants->codes[i],
                               neighbour_rank_code + index2);
              neighbour_rank_node_k[index2] = k;
              neighbour_rank_node_part[index2] = i_part;

              neighbour_rank_n[index]++;
            }
          }
        }
      }
    }


    /* Tri des codes pour chaque direction de chaque rang */
    int *order = malloc (sizeof(int) * max_node_dir);
    int *tmp_node_id = malloc (sizeof(int) * max_node_dir);
    PDM_morton_code_t *tmp_code = malloc (sizeof(PDM_morton_code_t) * max_node_dir);
    int *tmp_node_k = malloc (sizeof(int) * max_node_dir);
    int *tmp_node_part = malloc (sizeof(int) * max_node_dir);

    for (int i = 0; i < n_ranks; i++) {
      for (int j = 0; j < n_direction; j++) {
        PDM_morton_local_order (neighbour_rank_n[n_direction * i + j],
                                neighbour_rank_code + neighbour_rank_idx[n_direction * i + j],
                                order);
        int idx1 = 0;
        for (int k = neighbour_rank_idx[n_direction * i + j];
             k < neighbour_rank_idx[n_direction * i + j + 1];
             k++) {
          PDM_morton_copy (neighbour_rank_code[k], tmp_code + idx1);
          tmp_node_id[idx1]   = neighbour_rank_node_id[k];
          tmp_node_k[idx1]    = neighbour_rank_node_k[k];
          tmp_node_part[idx1] = neighbour_rank_node_part[k];
          idx1 += 1;
        }

        idx1 = 0;
        for (int k = neighbour_rank_idx[n_direction * i + j];
             k < neighbour_rank_idx[n_direction * i + j + 1];
             k++) {
          PDM_morton_copy (tmp_code[order[idx1]], neighbour_rank_code + k );
          neighbour_rank_node_id[k]   = tmp_node_id[order[idx1]];
          neighbour_rank_node_k[k]    = tmp_node_k[order[idx1]];
          neighbour_rank_node_part[k] = tmp_node_part[order[idx1]];
          idx1 += 1;
        }
      }
    }

    free (tmp_code);
    free (order);
    free (tmp_node_id);
    free (tmp_node_k);
    free (tmp_node_part);

    /* Envoi / reception (Les donnees recues sont triees) */

    int *recv_neighbour_rank_n = malloc (sizeof(int) * n_quantile);

    for (int i = 0; i < n_quantile; i++) {
      recv_neighbour_rank_n[i] = 0;
    }

    PDM_MPI_Request *recv_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);
    PDM_MPI_Request *send_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);

    int *used_ranks = malloc (sizeof(int) * n_ranks);

    PDM_MPI_Alltoall (neighbour_rank_n, n_direction, PDM_MPI_INT,
                      recv_neighbour_rank_n, n_direction, PDM_MPI_INT,
                      octree->comm);

    int *recv_neighbour_rank_idx = malloc (sizeof(int) * (n_direction * n_ranks + 1));
    recv_neighbour_rank_idx[0] = 0;


    for (int i = 0; i <  n_direction * n_ranks; i++)
      recv_neighbour_rank_idx[i+1] = recv_neighbour_rank_idx[i] + recv_neighbour_rank_n[i];



    int *recv_neighbour_rank_node_id   = malloc (sizeof(int) * recv_neighbour_rank_idx[n_quantile]);
    int *recv_neighbour_rank_node_part = malloc (sizeof(int) * recv_neighbour_rank_idx[n_quantile]);
    PDM_morton_code_t *recv_neighbour_rank_code =
      malloc (sizeof(PDM_morton_code_t) * recv_neighbour_rank_idx[n_quantile]);


    unsigned int *_neighbour_rank_code =
      malloc (sizeof(unsigned int) * 4 * neighbour_rank_idx[n_quantile]);
    unsigned int *_recv_neighbour_rank_code =
      malloc (sizeof(unsigned int) * 4 * recv_neighbour_rank_idx[n_quantile]);

    int idx = 0;
    for (int i = 0; i < neighbour_rank_idx[n_quantile]; i++) {
      _neighbour_rank_code[idx++] = neighbour_rank_code[i].L;
      for (int j = 0; j < 3; j++) {
        _neighbour_rank_code[idx++] = neighbour_rank_code[i].X[j];
      }
    }

    int *rank_neighbour_rank_n = malloc (sizeof(int) * n_ranks);
    int *rank_neighbour_rank_idx = malloc (sizeof(int) * (n_ranks + 1));
    int *rank_recv_neighbour_rank_n = malloc (sizeof(int) * n_ranks);
    int *rank_recv_neighbour_rank_idx = malloc (sizeof(int) * (n_ranks + 1));

    rank_neighbour_rank_idx[0] = 0;
    rank_recv_neighbour_rank_idx[0] = 0;

    for (int i = 0; i < n_ranks; i++) {
      rank_neighbour_rank_n[i] = 0;
      rank_recv_neighbour_rank_n[i] = 0;
    }

    for (int i = 0; i < n_ranks; i++) {
      for (int j = 0; j < n_direction; j++) {
        rank_neighbour_rank_n[i] += neighbour_rank_n[i*n_direction+j];
        rank_recv_neighbour_rank_n[i] += recv_neighbour_rank_n[i*n_direction+j];
      }
      rank_neighbour_rank_idx[i+1] = rank_neighbour_rank_n[i] + rank_neighbour_rank_idx[i];
      rank_recv_neighbour_rank_idx[i+1] = rank_recv_neighbour_rank_n[i] + rank_recv_neighbour_rank_idx[i];
    }

    PDM_MPI_Alltoallv (neighbour_rank_node_id,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_INT,
                       recv_neighbour_rank_node_id,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_INT,
                       octree->comm);

    PDM_MPI_Alltoallv (neighbour_rank_node_part,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_INT,
                       recv_neighbour_rank_node_part,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_INT,
                       octree->comm);

    for (int i = 0; i < n_ranks; i++) {
      rank_neighbour_rank_n[i] *= 4;
      rank_recv_neighbour_rank_n[i] *= 4;
      rank_neighbour_rank_idx[i+1] *= 4;
      rank_recv_neighbour_rank_idx[i+1] *= 4;
    }


    PDM_MPI_Alltoallv (_neighbour_rank_code,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_UNSIGNED,
                       _recv_neighbour_rank_code,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_UNSIGNED,
                       octree->comm);


    free (_neighbour_rank_code);

    free (rank_neighbour_rank_n);
    free (rank_neighbour_rank_idx);
    free (rank_recv_neighbour_rank_n);
    free (rank_recv_neighbour_rank_idx);

    idx = 0;
    for (int i = 0; i < recv_neighbour_rank_idx[n_quantile]; i++) {
      recv_neighbour_rank_code[i].L = _recv_neighbour_rank_code[idx++];
      for (int j = 0; j < 3; j++) {
        recv_neighbour_rank_code[i].X[j] = _recv_neighbour_rank_code[idx++];
      }
    }
    free (_recv_neighbour_rank_code);

    free (recv_request);
    free (send_request);

    free (used_ranks);


    octree->n_part_boundary_elt = neighbour_rank_idx[n_quantile];
    octree->part_boundary_elt_idx = malloc (sizeof(int) * (octree->n_part_boundary_elt + 1));

    int s_part_boundary_elt = 2 * 3 * neighbour_rank_idx[n_quantile];
    octree->part_boundary_elt = malloc (sizeof(int) * s_part_boundary_elt);

    int n_part_boundary_elt = 0;

    idx = 0;
    int idx_part_boundary_elt = 0;
    for (int i = 0; i <= octree->n_part_boundary_elt; i++)
      octree->part_boundary_elt_idx[i] = 0;


    FALSE_NEIGHBOUR = -(octree->n_part_boundary_elt + 1);

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (int j = 0; j < n_direction; j++) {
        for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
          if (neighbours_tmp[i].neighbours[j][k] < 0)
            neighbours_tmp[i].neighbours[j][k] = FALSE_NEIGHBOUR;
        }
      }
    }



    for (int i = 0; i < n_ranks; i++ ) {
      for (PDM_para_octree_direction_t j = PDM_BOTTOM; j < PDM_N_DIRECTION; j++) {
        PDM_para_octree_direction_t inv_j = _inv_direction(j);

        int idx_recv = n_direction * i + inv_j;

        int idx_candidate = recv_neighbour_rank_idx[idx_recv];
        int n_candidate = recv_neighbour_rank_idx[idx_recv+1] - idx_candidate;

        if (n_candidate > 0) {

          for (int k = neighbour_rank_idx[i * n_direction + j];
               k < neighbour_rank_idx[i * n_direction + j + 1]; k++) {
            PDM_morton_code_t *neighbour_code = _neighbour (neighbour_rank_code[k], j);

            size_t start_intersect, end_intersect;
            PDM_morton_list_intersect (n_candidate,
                                       *neighbour_code,
                                       recv_neighbour_rank_code + idx_candidate,
                                       &start_intersect,
                                       &end_intersect);

            int n_intersect_neighbours = 0;

            if (end_intersect > start_intersect) {
              for (int k1 = start_intersect; k1 < (int )end_intersect; k1++) {
                int k2 = idx_candidate + k1;

                PDM_morton_code_t *neighbour_neighbour_code =
                  _neighbour (recv_neighbour_rank_code[k2], inv_j);

                assert (neighbour_neighbour_code != NULL);

                if (PDM_morton_ancestor_is (neighbour_rank_code[k], *neighbour_neighbour_code) ||
                    PDM_morton_ancestor_is (*neighbour_neighbour_code, neighbour_rank_code[k])) {
                  n_intersect_neighbours++;

                  if ((s_part_boundary_elt - idx_part_boundary_elt) <= 3) {
                    s_part_boundary_elt *= 2;
                    octree->part_boundary_elt = realloc (octree->part_boundary_elt,
                                                         sizeof(int) * s_part_boundary_elt);
                  }
                  octree->part_boundary_elt_idx[n_part_boundary_elt+1]++;
                  octree->part_boundary_elt[idx_part_boundary_elt++] = i; // rank
                  octree->part_boundary_elt[idx_part_boundary_elt++] = recv_neighbour_rank_node_id[k2]; // neighbour's local number in rank i
                  octree->part_boundary_elt[idx_part_boundary_elt++] = recv_neighbour_rank_node_part[k2]; // neighbour's part number in rank i
                }

                free (neighbour_neighbour_code);
              }
            }

            int k3 = neighbour_rank_node_k[k];
            int i2 = neighbour_rank_node_id[k];

            if (n_intersect_neighbours > 0) {
              neighbours_tmp[i2].neighbours[j][k3] = -(n_part_boundary_elt+1);

              assert (neighbours_tmp[i2].neighbours[j][k3] != FALSE_NEIGHBOUR);

              n_part_boundary_elt++;
            }
            else {
              neighbours_tmp[i2].neighbours[j][k3] = FALSE_NEIGHBOUR;
            }

            free (neighbour_code);
          }
        }
      }
    }

    free (neighbour_rank_n);
    free (neighbour_rank_idx);
    free (neighbour_rank_node_id);
    free (neighbour_rank_node_k);
    free (neighbour_rank_node_part);
    free (neighbour_rank_code);

    free (recv_neighbour_rank_n);
    free (recv_neighbour_rank_idx);

    free (recv_neighbour_rank_node_id);
    free (recv_neighbour_rank_node_part);
    free (recv_neighbour_rank_code);

    octree->n_part_boundary_elt = n_part_boundary_elt;
  }

  for (int i = 0; i < octree->n_part_boundary_elt; i++) {
    octree->part_boundary_elt_idx[i+1] += octree->part_boundary_elt_idx[i];
  }


  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_DISTANT_NEIGHBOURS] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_DISTANT_NEIGHBOURS]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_DISTANT_NEIGHBOURS]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_DISTANT_NEIGHBOURS]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);



  /*************************************************************************
   *
   * Copy temporary neighbours in the neighbour structure
   *
   *************************************************************************/
  octree->octants->neighbour_idx = malloc(sizeof(int) * (n_direction * octree->octants->n_nodes + 1));

  int idx = 0;
  octree->octants->neighbour_idx[0] = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      octree->octants->neighbour_idx[idx+1] =
        octree->octants->neighbour_idx[idx] + neighbours_tmp[i].n_neighbour[j];

      /* account for false distant neighbours */
      for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
        if (neighbours_tmp[i].neighbours[j][k] == FALSE_NEIGHBOUR)
          octree->octants->neighbour_idx[idx+1]--;
      }

      idx += 1;
    }
  }

  octree->octants->neighbours =
    malloc(sizeof(int) *
           octree->octants->neighbour_idx[n_direction * octree->octants->n_nodes]);

  idx = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {

        if (neighbours_tmp[i].neighbours[j][k] != FALSE_NEIGHBOUR)
          octree->octants->neighbours[idx++] = neighbours_tmp[i].neighbours[j][k];

      }
    }
  }

  /* Free temporary arrays */
  /* printf("sortie 2 neighbours_tmp debut\n"); */
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      if (neighbours_tmp[i].neighbours[j] != NULL) {
        free (neighbours_tmp[i].neighbours[j]);
      }
    }
  }
  /* printf("sortie 2 neighbours_tmp fin\n"); */

  free (neighbours_tmp);




  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP3] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP3]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP3]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP3]   += e_t_cpu_s - b_t_cpu_s;


  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS] += octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS] += octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS] += octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS] += octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP3];

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);
}


static void
_check_neighbours_area
(
 const _octree_t *octree
 )
{
  int my_rank, n_ranks;
  PDM_MPI_Comm_rank (octree->comm, &my_rank);
  PDM_MPI_Comm_size (octree->comm, &n_ranks);

  if (my_rank == 0) {
    printf("-- Check neighbours\n");
  }

  _l_octant_t *octants = octree->octants;
  double *area = malloc (sizeof(double) * octants->n_nodes);

  int *rank_ngb_n = malloc (sizeof(double) * n_ranks);
  for (int i = 0; i < n_ranks; i++)
    rank_ngb_n[i] = 0;

  for (int i = 0; i < octants->n_nodes; i++) {
    PDM_morton_code_t code = octants->codes[i];

    /* check neighbours of current octant */
    area[i] = 0;
    for (int j = 0; j < 6; j++) {
      for (int k = octants->neighbour_idx[6*i+j];
           k < octants->neighbour_idx[6*i+j+1]; k++) {
        int ingb = octants->neighbours[k];

        if (ingb < 0) {
          // distant neighbour
          ingb = -(ingb + 1);
          for (int l = octree->part_boundary_elt_idx[ingb];
               l < octree->part_boundary_elt_idx[ingb+1]; l++) {
            int ngb_rank = octree->part_boundary_elt[3*l];
            rank_ngb_n[ngb_rank]++;
          }

        } else {
          // local neighbour
          PDM_morton_code_t ngb_code = octants->codes[ingb];
          //          double side = 1./pow(2, PDM_MAX(code.L, ngb_code.L));
          double side = 1./(double)(1 << PDM_MAX(code.L, ngb_code.L));
          area[i] += side * side;
        }
      }
    }
  }

  // MPI communications to check distant neighbours
  int *rank_ngb_idx = NULL;
  int *rank_ngb_id_level = NULL;
  int *recv_rank_ngb_n = NULL;
  int *recv_rank_ngb_idx = NULL;
  int *recv_rank_ngb_id_level = NULL;

  if (n_ranks > 1) {
    rank_ngb_idx = malloc (sizeof(int) * (n_ranks + 1));
    rank_ngb_idx[0] = 0;
    for (int i = 0; i < n_ranks; i++) {
      rank_ngb_idx[i+1] = rank_ngb_idx[i] + 2*rank_ngb_n[i];
      rank_ngb_n[i] = 0;
    }
    rank_ngb_id_level = malloc (sizeof(int) * rank_ngb_idx[n_ranks]);

    for (int i = 0; i < octants->n_nodes; i++) {
      for (int j = 0; j < 6; j++) {
        for (int k = octants->neighbour_idx[6*i+j];
             k < octants->neighbour_idx[6*i+j+1]; k++) {
          int ingb = octants->neighbours[k];

          if (ingb < 0) {
            ingb = -(ingb + 1);

            for (int l = octree->part_boundary_elt_idx[ingb];
                 l < octree->part_boundary_elt_idx[ingb+1]; l++) {
              int ngb_rank = octree->part_boundary_elt[3*l];
              int ngb_id = octree->part_boundary_elt[3*l+1];
              int idx = rank_ngb_idx[ngb_rank] + rank_ngb_n[ngb_rank];
              rank_ngb_id_level[idx++] = ngb_id;
              rank_ngb_id_level[idx++] = (int) octants->codes[i].L;
              rank_ngb_n[ngb_rank] += 2;
            }
          }
        }
      }
    }

    recv_rank_ngb_n = malloc (sizeof(int) * n_ranks);
    PDM_MPI_Alltoall (rank_ngb_n, 1, PDM_MPI_INT,
                      recv_rank_ngb_n, 1, PDM_MPI_INT,
                      octree->comm);

    recv_rank_ngb_idx = malloc (sizeof(int) * (n_ranks + 1));
    recv_rank_ngb_idx[0] = 0;
    for (int i = 0; i < n_ranks; i++) {
      recv_rank_ngb_idx[i+1] = recv_rank_ngb_idx[i] + recv_rank_ngb_n[i];
    }

    recv_rank_ngb_id_level = malloc (sizeof(int) * recv_rank_ngb_idx[n_ranks]);
    PDM_MPI_Alltoallv (rank_ngb_id_level, rank_ngb_n, rank_ngb_idx, PDM_MPI_INT,
                       recv_rank_ngb_id_level, recv_rank_ngb_n, recv_rank_ngb_idx, PDM_MPI_INT,
                       octree->comm);

    free (rank_ngb_id_level);
    free (rank_ngb_n);
    free (rank_ngb_idx);


    for (int i = 0; i < n_ranks; i++) {
      recv_rank_ngb_n[i] /= 2;
      recv_rank_ngb_idx[i+1] /= 2;
    }
    for (int i = 0; i < n_ranks; i++) {
      for (int j = recv_rank_ngb_idx[i]; j < recv_rank_ngb_idx[i+1]; j++) {
        int id = recv_rank_ngb_id_level[2*j];
        PDM_morton_int_t level = (PDM_morton_int_t) recv_rank_ngb_id_level[2*j+1];
        double side = 1./(double) (1 << PDM_MAX (level, octants->codes[id].L));

        area[id] += side * side;
      }
    }

    free (recv_rank_ngb_id_level);
    free (recv_rank_ngb_n);
    free (recv_rank_ngb_idx);
  }



  for (int i = 0; i < octants->n_nodes; i++) {
    /* compute exact interior surface area of current octant */
    PDM_morton_code_t code = octants->codes[i];
    double side = 1./ (double) (1 << code.L);
    int ndir = 0;
    for (int j = 0; j < 6; j++) {
      PDM_morton_code_t *ngb_code = _neighbour (code, (PDM_para_octree_direction_t) j);
      if (ngb_code != NULL) {
        ndir++;
        free (ngb_code);
      }
    }
    double exact_area = ndir * side * side;

    /* compare with actual area */
    /*printf("node #%d: exact area = %f, actual area = %f, err = %f\n",
      i,
      exact_area,
      area[i],
      PDM_ABS(area[i]/exact_area - 1));*/
    //-->>
    /*if (PDM_ABS(area[i]/exact_area - 1) > 1e-15) {
      printf("[%d] node %d, level %u, area = %f, exact = %f, relative area error = %f\n",
      my_rank, i, octants->codes[i].L, area[i], exact_area, area[i]/exact_area - 1);
      printf("\tneighbours (rank, node_id):\n");
      for (int j = 0; j < 6; j++) {
      printf("\t\tdirection %d:", j);
      for (int k = octants->neighbour_idx[6*i+j];
      k < octants->neighbour_idx[6*i+j+1]; k++) {
      int ingb = octants->neighbours[k];
      if (ingb < 0) {
      ingb = -(ingb + 1);

      for (int l = octree->part_boundary_elt_idx[ingb];
      l < octree->part_boundary_elt_idx[ingb+1]; l++) {
      int ngb_rank = octree->part_boundary_elt[2*l];
      int ngb_id = octree->part_boundary_elt[2*l+1];
      printf(" (%d, %d)", ngb_rank, ngb_id);
      }
      } else {
      printf(" (%d, %d)", my_rank, ingb);
      }
      }
      printf("\n");
      }
      printf("\n");
      }*/
    //<<--
    //assert (PDM_ABS(area[i]/exact_area - 1) < 1e-15);
    assert (PDM_ABS(area[i] - exact_area) <= 1e-15 * exact_area);
  }

  free (area);
}




/*
 * ---> pdm_morton.c
 */
inline static void
_maximal_intersecting_range
(
 const PDM_morton_code_t  lo,
 const PDM_morton_code_t  hi,
 const PDM_morton_code_t *codes,
 const int                n_codes,
 size_t                  *l,
 size_t                  *r
 )
{
  if (PDM_morton_a_gt_b (codes[0], hi)) {
    *l = 0;
    *r = *l;
    return;
  } else if (PDM_morton_a_gt_b (lo, codes[n_codes-1])) {
    *l = n_codes-1;
    *r = *l;
  }


  // get left bound (included)
  size_t _l = 0;
  size_t _r = n_codes;
  while (_l < _r - 1) {
    size_t m = _l + (_r - _l) / 2;
    if (PDM_morton_a_gt_b (codes[m], lo)) {
      _r = m;
    } else {
      _l = m;
    }
  }

  *l = _l;

  // get right bound (excluded)
  _l = _l - 1;
  _r = n_codes - 1;
  while (_l < _r - 1) {
    size_t m = _r - (_r - _l) / 2;
    if (PDM_morton_a_gt_b (hi, codes[m])) {
      _l = m;
    } else {
      _r = m;
    }
  }

  *r = _r + 1;
}






static void
_closest_points_local
(
 const _octree_t  *octree,
 const int         n_closest_points,
 const int         n_pts,
 double           *pts_coord,
 PDM_g_num_t      *pts_g_num, // ONLY FOR DEBUG
 int              *start_leaves,
 int              *start_leaves_idx,
 double           *upper_bound_dist,
 PDM_g_num_t      *local_closest_src_gnum,
 double           *local_closest_src_dist,
 int              *send_count,
 int              *send_shift,
 int             **send_tgt_lnum,
 int             **send_start_leaves,
 int             **send_start_leaves_count,
 int             **send_start_leaves_rank_shift
 )
{
  const int DEBUG = 0;
  const int CHECK_FACE_DIST = 1;
  //const int CHECK_EMPTY_START = 0;
  const int CHECK_CLOSEST = 1;
  const double THRESHOLD_CLOSEST = 0.75; // check closest if min_dist > threshold * upper bound
  const int CHECK_INTERSECT = 1;
  const int NORMALIZE = 1;

  const _l_octant_t *octants = octree->octants;
  const int dim = octree->dim;

  int myRank, lComm;
  PDM_MPI_Comm_rank (octree->comm, &myRank);
  PDM_MPI_Comm_size (octree->comm, &lComm);

  for (int i = 0; i < lComm; i++) {
    send_count[i] = 0;
  }

  //--->>>
  int **tmp_send_tgt_lnum       = malloc (sizeof(int *) * lComm);
  int **tmp_send_tgt_n_leaves   = malloc (sizeof(int *) * lComm);
  int **tmp_send_start_leaves   = malloc (sizeof(int *) * lComm);
  int  *s_tmp_send_start_leaves = malloc (sizeof(int) * lComm);
  int  *n_tmp_send_start_leaves = malloc (sizeof(int) * lComm);

  int **send_to_rank_leaves   = malloc (sizeof(int *) * lComm);
  int  *n_send_to_rank_leaves = malloc (sizeof(int) * lComm);
  int  *s_send_to_rank_leaves = malloc (sizeof(int) * lComm);
  for (int i = 0; i < lComm; i++) {
    s_send_to_rank_leaves[i] = 16;// ?
    if (i != myRank) {
      send_to_rank_leaves[i] = malloc (sizeof(int) * 2 * s_send_to_rank_leaves[i]);

      tmp_send_tgt_lnum[i] = malloc (sizeof(int) * n_pts);// could be smaller and dynamically re-alloc'ed when necessary
      tmp_send_tgt_n_leaves[i] = malloc (sizeof(int) * n_pts);// idem
      s_tmp_send_start_leaves[i] = 16;//?
      tmp_send_start_leaves[i] = malloc (sizeof(int) * s_tmp_send_start_leaves[i]);
    } else {
      s_tmp_send_start_leaves[i] = 0;
    }
    n_tmp_send_start_leaves[i] = 0;
  }
  //<<<---

  int *is_visited_part = malloc (sizeof(int) * octree->n_connected);

  int *is_visited     = malloc (sizeof(int) * octants->n_nodes);
  int *visited_leaves = malloc (sizeof(int) * octants->n_nodes);// could be smaller...
  int n_visited = 0;
  for (int i = 0; i < octants->n_nodes; i++) {
    is_visited[i] = 0;
  }

  /* Min heap used to sort start leaves */
  _min_heap_t *start_heap = _min_heap_create (start_leaves_idx[n_pts]); // smaller size?

  /* Min heap used to visit leaves from neighbour to neighbour */
  _min_heap_t *leaf_heap = _min_heap_create (octants->n_nodes);


  /* Loop over target points */
  for (int i_tgt = 0; i_tgt < n_pts; i_tgt++) {

    /* Init */
    for (int i = 0; i < octree->n_connected; i++) {
      is_visited_part[i] = 0;
    }

    for (int i = 0; i < lComm; i++) {
      n_send_to_rank_leaves[i] = 0;
    }

    for (int i = 0; i < n_visited; i++) {
      is_visited[visited_leaves[i]] = 0;
    }
    n_visited = 0;

    PDM_g_num_t *closest_src_gnum = local_closest_src_gnum + n_closest_points * i_tgt;
    double      *closest_src_dist = local_closest_src_dist + n_closest_points * i_tgt;
    for (int j = 0; j < n_closest_points; j++) {
      closest_src_dist[j] = upper_bound_dist[i_tgt];
      closest_src_gnum[j] = -1; // USEFUL ONLY FOR DEBUG
    }

    double *max_src_dist = closest_src_dist + n_closest_points - 1;


    /* Get current target point data */
    const double *_pt = pts_coord + i_tgt * dim;
    double _ptn[3];
    if (NORMALIZE) {
      for (int i_dim = 0; i_dim < dim; i_dim++) {
        _ptn[i_dim] = (_pt[i_dim] - octree->s[i_dim]) / octree->d[i_dim];
      }
    }

    int n_start_leaves = start_leaves_idx[i_tgt+1] - start_leaves_idx[i_tgt];
    if (DEBUG) {
      printf("\n=== pt ("PDM_FMT_G_NUM") (upper_bound_dist = %f) ===\nstart leaves (%d):\n",
             pts_g_num[i_tgt], upper_bound_dist[i_tgt], n_start_leaves);
    }

    if (n_start_leaves < 1) {
      continue; /* move on to next target point */
    }

    /* Sort start leaves in ascending order of min distance from tgt point */
    _min_heap_reset (start_heap);
    double min_start_dist = *max_src_dist;
    for (int i_start = start_leaves_idx[i_tgt]; i_start < start_leaves_idx[i_tgt+1]; i_start++) {

      int leaf_id = start_leaves[i_start];
      int part_id = _get_octant_part_id (octree,
                                         leaf_id);
      is_visited_part[part_id] = 1;

      double start_dist;
      if (NORMALIZE) {
        start_dist = _octant_min_dist2_normalized (dim,
                                                   octants->codes[leaf_id],
                                                   octree->d,
                                                   _ptn);
      } else {
        start_dist = _octant_min_dist2 (dim,
                                        octants->codes[leaf_id],
                                        octree->d,
                                        octree->s,
                                        _pt);
      }

      min_start_dist = PDM_MIN (min_start_dist, start_dist);

      if (start_dist < *max_src_dist) {
        _min_heap_push (start_heap,
                        leaf_id,
                        0,
                        start_dist);
      }

      if (DEBUG) {
        int i_part = _get_octant_part_id (octree, leaf_id);
        printf("\t%d (part %d): dist = %f\n", leaf_id, i_part, start_dist);
      }
    }
    if (DEBUG) {
      printf("============================\n");
    }


    /* Check whether start_heap is empty */
    //if (CHECK_EMPTY_START) {
    //if (start_heap->count == 0 && n_start_leaves > 0) {
    if (CHECK_CLOSEST && n_start_leaves > 0) {
      if (min_start_dist >= THRESHOLD_CLOSEST * (*max_src_dist)) {
        /* Look for the closest octant (within the search radius), which would have been missed if we returned too early) */
        double            box[2*dim];
        PDM_morton_code_t box_corners[2];
        double            width = sqrt (*max_src_dist);
        double s[3], d[3];
        /* Encode corners of search box */
        for (int i_dim = 0; i_dim < dim; i_dim++) {
          box[i_dim]       = PDM_MAX (octree->global_extents[i_dim],
                                      _pt[i_dim] - width);
          box[dim + i_dim] = PDM_MAX (octree->global_extents[dim + i_dim],
                                      _pt[i_dim] + width);
        }

        PDM_morton_encode_coords (dim,
                                  PDM_morton_max_level,
                                  octree->global_extents,
                                  (size_t) 2,
                                  box,
                                  box_corners,
                                  d,
                                  s);

        /* Loop over visited parts */
        for (int i_part = 0; i_part < octree->n_connected; i_part++) {

          if (!is_visited_part[i_part]) {
            continue;
          }

          /* Get range of octants that intersect the search box in the sense of Morton codes */
          size_t l, r;
          _maximal_intersecting_range (box_corners[0],
                                       box_corners[1],
                                       octants->codes + octree->connected_idx[i_part],
                                       octree->connected_idx[i_part+1] - octree->connected_idx[i_part],
                                       &l,
                                       &r);

          l += octree->connected_idx[i_part];
          r += octree->connected_idx[i_part];

          /* Narrow down this range as much as possible */
          // BIGMIN/LITMAX? ...

          /* Find closest leaf within this range */
          int    closest_id = -1;
          double min_dist = *max_src_dist;
          double dist;
          for (int i = l; i < (int) r; i++) {

            if (CHECK_INTERSECT) {
              /* make sure leaf intersects search box so it's worth computing min distance */
              int intersect = 1;
              PDM_morton_code_t *code = octants->codes + i;

              double min_octant[dim];
              double max_octant[dim];
              double side = 1. / (double) (1 << code->L);

              for (int j = 0; j < dim; j++) {
                min_octant[j] = octree->s[j] + octree-> d[j] * side * code->X[j];
                max_octant[j] = min_octant[j] + octree->d[j] * side;

                if (max_octant[j] < box[j] || min_octant[j] > box[dim + j]) {
                  intersect = 0;
                  break;
                }
              }

              if (intersect) {
                dist = 0.;
                double delta;
                for (int j = 0; j < dim; j++) {
                  double x = _pt[j];

                  if (x > max_octant[j]) {
                    delta = x - max_octant[j];
                    dist += delta * delta;
                  } else if (x < min_octant[j]) {
                    delta = x - min_octant[j];
                    dist += delta * delta;
                  }
                }
              } else {
                dist = HUGE_VAL;
              }

            } else {
              if (NORMALIZE) {
                dist = _octant_min_dist2_normalized (dim,
                                                     octants->codes[i],
                                                     octree->d,
                                                     _ptn);
              } else {
                dist = _octant_min_dist2 (dim,
                                          octants->codes[i],
                                          octree->d,
                                          octree->s,
                                          _pt);
              }
            } // end if (CHECK_INTERSECT)

            if (dist < min_dist) {
              closest_id = i;
              min_dist = dist;
            }
          }

          if (DEBUG) {
            printf("[%d] closest octant: %d, dist = %f / %f\n",
                   myRank, closest_id, min_dist, *max_src_dist);
          }

          /* If the closest octant is within the search radius, push it in start_heap */
          if (closest_id >= 0) {
            _min_heap_push (start_heap,
                            closest_id,
                            0,
                            min_dist);
          }
        } // end loop over visited parts
      } // end if (min_start_dist >= THRESHOLD_CLOSEST * (*max_src_dist))
    } // end if (CHECK_CLOSEST && n_start_leaves > 0)

    /* Loop over (sorted) start leaves */
    int start_id;
    PDM_g_num_t unused;
    double start_dist;
    while (_min_heap_pop (start_heap, &start_id, &unused, &start_dist)) {

      if (DEBUG) {
        printf("tgt point ("PDM_FMT_G_NUM") start leaf %d: dist = %f / %f  (is visited? %d)\n",
               pts_g_num[i_tgt], start_id, start_dist, *max_src_dist, is_visited[start_id]);
      }

      if (start_dist >= *max_src_dist) {
        break;
      }

      if (is_visited[start_id]) {
        continue;
      }

      /* Push start leaf in priority queue */
      _min_heap_reset (leaf_heap);
      _min_heap_push (leaf_heap,
                      start_id,
                      0,
                      start_dist);
      is_visited[start_id] = 1;
      visited_leaves[n_visited++] = start_id;

      /* Visit octree leaves from neighbour to neighbour using priority queue */
      int leaf_id;
      double leaf_dist;
      while (_min_heap_pop (leaf_heap, &leaf_id, &unused, &leaf_dist)) {

        if (DEBUG) {
          printf("tgt point ("PDM_FMT_G_NUM") inspecting leaf %d: dist = %f / %f\n",
                 pts_g_num[i_tgt], leaf_id, leaf_dist, *max_src_dist);
        }

        if (leaf_dist >= *max_src_dist) {
          break;
        }

        /* inspect source points inside popped leaf */
        for (int i = 0; i < octants->n_points[leaf_id]; i++) {

          // get source point coords and gnum
          int i_src = octants->range[leaf_id] + i;
          double *src_pt = octree->points + i_src * dim;
          PDM_g_num_t src_gnum = octree->points_gnum[i_src];

          // compute (squared) distance from target point
          double src_dist = 0;
          for (int j = 0; j < dim; j++) {
            double delta = _pt[j] - src_pt[j];
            src_dist += delta * delta;
          }

          if (DEBUG && pts_g_num[i_tgt] == 38) {
            printf("  src pt ("PDM_FMT_G_NUM") [%f %f %f] at dist %f / %f\n",
                   src_gnum, src_pt[0], src_pt[1], src_pt[2], src_dist, *max_src_dist);
          }

          // insertion sort
          if (src_dist < *max_src_dist) {
            int j = n_closest_points - 1;
            while (j > 0 && src_dist < closest_src_dist[j-1]) {
              closest_src_gnum[j] = closest_src_gnum[j-1];
              closest_src_dist[j] = closest_src_dist[j-1];
              j--;
            }

            closest_src_gnum[j] = src_gnum;
            closest_src_dist[j] = src_dist;

            if (DEBUG) {
              printf("  src pt ("PDM_FMT_G_NUM") [%f %f %f] at dist %f / %f --> insert at pos %d\n",
                     src_gnum, src_pt[0], src_pt[1], src_pt[2], src_dist, *max_src_dist, j);
            }
          }
        } // end loop over source points inside popped leaf


        /* inspect neighbours of popped leaf */
        double side = 1./(double) (1 << octants->codes[leaf_id].L);
        int check_dist;
        for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
          if (CHECK_FACE_DIST) {
            /* coarse test to discard neighbours that are clearly beyond the current search radius */
            int i_dim = dir / 2;
            int sign = dir % 2;

            double x_face = octree->s[i_dim] + octree->d[i_dim] * side * (octants->codes[leaf_id].X[i_dim] + sign);

            if (sign == 0) {
              check_dist = _pt[i_dim] > x_face;
            } else {
              check_dist = _pt[i_dim] < x_face;
            }

            if (check_dist) {
              double dist_face = (_pt[i_dim] -  x_face) * (_pt[i_dim] -  x_face);

              if (dist_face >= *max_src_dist) {
                /* all neighbours in direction dir are at least
                   at a distance from the tgt point equal to dist_face
                   so they need not be inspected */
                continue;
              }
            }
          }

          /* inspect neighbours in direction dir */
          for (int i = octants->neighbour_idx[6*leaf_id + dir];
               i < octants->neighbour_idx[6*leaf_id + dir + 1]; i++) {

            int ngb = octants->neighbours[i];

            if (ngb < 0) {
              // distant neighbour(s)
              ngb = -(ngb + 1);

              for (int j = octree->part_boundary_elt_idx[ngb];
                   j < octree->part_boundary_elt_idx[ngb+1]; j++) {

                int ngb_rank = octree->part_boundary_elt[3*j];
                int ngb_leaf = octree->part_boundary_elt[3*j+1];
                int ngb_part = octree->part_boundary_elt[3*j+2];

                if (n_send_to_rank_leaves[ngb_rank] == 0) {
                  tmp_send_tgt_lnum[ngb_rank][send_count[ngb_rank]] = i_tgt;
                  tmp_send_tgt_n_leaves[ngb_rank][send_count[ngb_rank]] = 0;
                  send_count[ngb_rank]++;
                }

                if (s_send_to_rank_leaves[ngb_rank] <= n_send_to_rank_leaves[ngb_rank]) {
                  s_send_to_rank_leaves[ngb_rank] *= 2;
                  send_to_rank_leaves[ngb_rank] = realloc (send_to_rank_leaves[ngb_rank],
                                                           sizeof(int) * 2 * s_send_to_rank_leaves[ngb_rank]);
                }
                if (DEBUG) {
                  printf("  Send pt ("PDM_FMT_G_NUM") to rank %d, leaf %d (part %d)\n",
                         pts_g_num[i_tgt], ngb_rank, ngb_leaf, ngb_part);
                }
                send_to_rank_leaves[ngb_rank][2*n_send_to_rank_leaves[ngb_rank]]   = ngb_leaf;
                send_to_rank_leaves[ngb_rank][2*n_send_to_rank_leaves[ngb_rank]+1] = ngb_part;
                n_send_to_rank_leaves[ngb_rank]++;

                tmp_send_tgt_n_leaves[ngb_rank][send_count[ngb_rank]-1]++;

              } // end loop over distant neighbours (j)

            } else {
              // local neighbour
              if (is_visited[ngb] == 1) continue;

              is_visited[ngb] = 1;
              visited_leaves[n_visited++] = ngb;

              // compute min dist from target point to current neighbor leaf
              double ngb_min_dist;
              if (NORMALIZE) {
                ngb_min_dist = _octant_min_dist2_normalized (dim,
                                                             octants->codes[ngb],
                                                             octree->d,
                                                             _ptn);
              } else {
                ngb_min_dist  = _octant_min_dist2 (dim,
                                                   octants->codes[ngb],
                                                   octree->d,
                                                   octree->s,
                                                   _pt);
              }

              if (ngb_min_dist < *max_src_dist) {
                // push current neighbour in priority queue
                _min_heap_push (leaf_heap,
                                ngb,
                                0,
                                ngb_min_dist);
              }

            } // end if local/distant neighbour
          } // end loop over neighbours in direction dir (i)
        } // end loop over directions (dir)
      } // end while (_min_heap_pop (leaf_heap))

    } // end loop over (sorted) start leaves

    for (int rank = 0; rank < lComm; rank++) {
      if (rank == myRank) continue;

      if (s_tmp_send_start_leaves[rank] <= n_tmp_send_start_leaves[rank] + n_send_to_rank_leaves[rank]) {
        s_tmp_send_start_leaves[rank] = PDM_MAX (2 * s_tmp_send_start_leaves[rank],
                                                 n_tmp_send_start_leaves[rank] + n_send_to_rank_leaves[rank]);
        tmp_send_start_leaves[rank] = realloc (tmp_send_start_leaves[rank],
                                               sizeof(int) * s_tmp_send_start_leaves[rank]);
      }

      for (int i = 0; i < n_send_to_rank_leaves[rank]; i++) {
        tmp_send_start_leaves[rank][n_tmp_send_start_leaves[rank]++] = send_to_rank_leaves[rank][2*i];
      }
    }

  } // end loop over target points (i_tgt)
  free (is_visited);
  free (is_visited_part);
  free (visited_leaves);
  _min_heap_free (start_heap);
  _min_heap_free (leaf_heap);


  free (n_send_to_rank_leaves);
  free (s_send_to_rank_leaves);
  for (int i = 0; i < lComm; i++) {
    if (i != myRank) {
      free (send_to_rank_leaves[i]);
    }
  }
  free (send_to_rank_leaves);

  /* Manage send_tgt_lnum buffer */
  send_shift[0] = 0;
  for (int i = 0; i < lComm; i++) {
    send_shift[i+1] = send_shift[i] + send_count[i];
  }

  *send_tgt_lnum = malloc (sizeof(int) * send_shift[lComm]);
  int idx = 0;
  for (int i = 0; i < lComm; i++) {
    for (int j = 0; j < send_count[i]; j++) {
      (*send_tgt_lnum)[idx] = tmp_send_tgt_lnum[i][j];
      idx++;
    }

    if (i != myRank) {
      free (tmp_send_tgt_lnum[i]);
    }
  }
  free (tmp_send_tgt_lnum);


  /* Manage send_start_leaves buffer */
  *send_start_leaves_rank_shift = malloc (sizeof(int) * (lComm+1));
  *send_start_leaves_count = malloc (sizeof(int) * send_shift[lComm]);
  (*send_start_leaves_rank_shift)[0] = 0;
  idx = 0;
  for (int i = 0; i < lComm; i++) {
    (*send_start_leaves_rank_shift)[i+1] = (*send_start_leaves_rank_shift)[i] + n_tmp_send_start_leaves[i];
    for (int j = 0; j < send_count[i]; j++) {
      (*send_start_leaves_count)[idx++] = tmp_send_tgt_n_leaves[i][j];
    }
  }

  *send_start_leaves = malloc (sizeof(int) * (*send_start_leaves_rank_shift)[lComm]);
  idx = 0;
  for (int i = 0; i < lComm; i++) {
    for (int j = 0; j < n_tmp_send_start_leaves[i]; j++) {
      (*send_start_leaves)[idx++] = tmp_send_start_leaves[i][j];
    }
  }



  for (int i = 0; i < lComm; i++) {
    if (i != myRank) {
      free (tmp_send_start_leaves[i]);
      free (tmp_send_tgt_n_leaves[i]);
    }
  }
  free (tmp_send_tgt_n_leaves);
  free (s_tmp_send_start_leaves);
  free (n_tmp_send_start_leaves);
  free (tmp_send_start_leaves);
}








/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create an octree structure
 *
 * \param [in]   n_point_cloud          Number of point cloud
 * \param [in]   depth_max              Maximum depth
 * \param [in]   points_in_leaf_max     Maximum points in a leaf
 * \param [in]   build_leaf_neighbours  Build leaf nieghbours (1 = true)
 * \param [in]   comm                   MPI communicator
 *
 * \return     Identifier
 */

int
PDM_para_octree_create
(
 const int n_point_cloud,
 const int depth_max,
 const int points_in_leaf_max,
 const int build_leaf_neighbours,
 const PDM_MPI_Comm comm
 )
{

  if (_octrees == NULL) {
    _octrees = PDM_Handles_create (4);
  }

  _octree_t *octree = (_octree_t *) malloc(sizeof(_octree_t));

  int id = PDM_Handles_store (_octrees, octree);

  octree->dim = 3;

  for (int i = 0; i < octree->dim; i++) {
    octree->global_extents[i]   = -HUGE_VAL;
    octree->global_extents[octree->dim+i] =  HUGE_VAL;
    octree->s[i]         = 0.;
    octree->d[i]         = 0.;
  }

  octree->depth_max = depth_max;
  octree->points_in_leaf_max = points_in_leaf_max;

  octree->n_point_clouds = n_point_cloud;
  octree->t_n_points = 0;
  octree->n_points = 0;
  octree->points = NULL;
  octree->points_icloud = NULL;
  octree->points_gnum = NULL;
  octree->points_code = NULL;

  octree->rank_octants_index = NULL;
  octree->octants = NULL;

  octree->n_part_boundary_elt = 0;
  octree->part_boundary_elt_idx = NULL;
  octree->part_boundary_elt = NULL;

  octree->neighboursToBuild = build_leaf_neighbours;

  octree->comm = comm;

  octree->n_connected = 0;
  octree->connected_idx = NULL;

  octree->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    octree->times_elapsed[i] = 0.;
    octree->times_cpu[i] = 0.;
    octree->times_cpu_u[i] = 0.;
    octree->times_cpu_s[i] = 0.;
  }

  return id;

}


/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_free
(
 const int          id
 )
{
  _octree_t *octree = _get_from_id (id);

  if (octree->points != NULL) {
    free (octree->points);
  }

  if (octree->points_icloud != NULL) {
    free (octree->points_icloud);
  }

  if (octree->points_gnum != NULL) {
    free (octree->points_gnum);
  }

  if (octree->points_code != NULL) {
    free (octree->points_code);
  }

  if (octree->part_boundary_elt_idx != NULL) {
    free (octree->part_boundary_elt_idx);
  }

  if (octree->part_boundary_elt != NULL) {
    free (octree->part_boundary_elt);
  }

  if (octree->rank_octants_index != NULL) {
    free (octree->rank_octants_index);
  }

  if (octree->octants != NULL) {

    if (octree->octants->codes != NULL) {
      free (octree->octants->codes);
    }

    if (octree->octants->range != NULL) {
      free (octree->octants->range);
    }

    if (octree->octants->n_points != NULL) {
      free (octree->octants->n_points);
    }

    if (octree->octants->neighbour_idx != NULL) {
      free (octree->octants->neighbour_idx);
    }

    if (octree->octants->neighbours != NULL) {
      free (octree->octants->neighbours);
    }

    free (octree->octants);
  }

  if (octree->connected_idx != NULL) {
    free (octree->connected_idx);
  }

  PDM_timer_free (octree->timer);

  free (octree);

  PDM_Handles_handle_free (_octrees, id, PDM_FALSE);

  const int n_octrees = PDM_Handles_n_get (_octrees);

  if (n_octrees == 0) {
    _octrees = PDM_Handles_free (_octrees);
  }
}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id                 Identifier
 * \param [in]   i_point_cloud      Number of point cloud
 * \param [in]   n_points           Maximum depth
 * \param [in]   coords             Point coordinates
 * \param [in]   g_num              Point global number or NULL
 *
 */


void
PDM_para_octree_point_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_points,
 const double      *coords,
 const PDM_g_num_t *g_num
 )
{
  _octree_t *octree = _get_from_id (id);

  const int idx = octree->n_points;

  octree->n_points += n_points;
  octree->points =
    realloc (octree->points, octree->n_points * sizeof(double) * octree->dim);
  octree->points_icloud =
    realloc (octree->points_icloud, octree->n_points * sizeof(int));
  octree->points_gnum =
    realloc (octree->points_gnum, octree->n_points * sizeof(PDM_g_num_t));
  octree->points_code =
    realloc (octree->points_code, octree->n_points * sizeof(PDM_morton_code_t));

  for (int i = 0; i < octree->dim * n_points; i++) {
    octree->points[octree->dim*idx + i] = coords[i];
  }

  for (int i = 0; i < n_points; i++) {
    octree->points_gnum[idx + i] = g_num[i];
  }

  for (int i = 0; i < n_points; i++) {
    octree->points_icloud[idx + i] = i_point_cloud;
  }

}


/**
 *
 * \brief Build octree
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_build
(
 const int  id,
 double    *global_extents
 )
{
  _octree_t *octree = _get_from_id (id);

  const int dim = octree->dim;
  const PDM_morton_int_t max_level = PDM_morton_max_level;

  int n_ranks;
  PDM_MPI_Comm_size (octree->comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  octree->times_elapsed[BEGIN] = PDM_timer_elapsed(octree->timer);
  octree->times_cpu[BEGIN]     = PDM_timer_cpu(octree->timer);
  octree->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(octree->timer);
  octree->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(octree->timer);

  b_t_elapsed = octree->times_elapsed[BEGIN];
  b_t_cpu     = octree->times_cpu[BEGIN];
  b_t_cpu_u   = octree->times_cpu_u[BEGIN];
  b_t_cpu_s   = octree->times_cpu_s[BEGIN];
  PDM_timer_resume(octree->timer);

  /*
   * Get coord extents
   */
  if (global_extents != NULL) {
    memcpy (octree->global_extents, global_extents, sizeof(double) * dim * 2);
  }

  else {
    PDM_morton_get_coord_extents(dim,
                                 octree->n_points,
                                 octree->points,
                                 octree->global_extents,
                                 octree->comm);

    /*
     * Dilate extents
     */
    double max_range = 0.;
    for (int i = 0; i < dim; i++) {
      max_range = PDM_MAX (max_range,
                           octree->global_extents[i+dim] - octree->global_extents[i]);
    }

    const double epsilon = 1.e-3 * max_range;

    for (int i = 0; i < dim; i++) {
      octree->global_extents[i]     -= 1.1 * epsilon; // On casse la symetrie !
      octree->global_extents[i+dim] +=       epsilon;
    }
  }

  /*
   * Encode coords
   */

  PDM_morton_encode_coords(dim,
                           max_level,
                           octree->global_extents,
                           octree->n_points,
                           octree->points,
                           octree->points_code,
                           octree->d,
                           octree->s);

  int *order = malloc (sizeof(int) * octree->n_points);

  for (int i = 0; i < octree->n_points; i++) {
    order[i] = i;
  }

  /**************************************
   *
   * Global order of codes and balancing
   *
   **************************************/

  PDM_morton_local_order (octree->n_points,
                          octree->points_code,
                          order);

  if (n_ranks > 1) {

    int *weight = malloc (sizeof(int) * octree->n_points);
    for (int i = 0; i < octree->n_points; i++) {
      weight[i] = 1;
    }

    PDM_morton_code_t *morton_index =
      malloc (sizeof(PDM_morton_code_t) * (n_ranks + 1));

    PDM_morton_build_rank_index(dim,
                                max_level,
                                octree->n_points,
                                octree->points_code,
                                weight,
                                order,
                                morton_index,
                                octree->comm);

    free (weight);
    free (order);

    /* distribute point from morton_index */

    _distribute_points (&octree->n_points,
                        &octree->points,
                        &octree->points_icloud,
                        &octree->points_gnum,
                        &octree->points_code,
                        morton_index,
                        octree->comm,
                        octree->dim,
                        max_level,
                        octree->global_extents);

    free(morton_index);

  }

  else {

    int *_points_icloud = malloc (sizeof(int) * octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {
      _points_icloud[i] =  octree->points_icloud[order[i]];
    }

    free (octree->points_icloud);
    octree->points_icloud = _points_icloud;

    PDM_g_num_t *_points_gnum = malloc (sizeof(PDM_g_num_t) * octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {
      _points_gnum[i] =  octree->points_gnum[order[i]];
    }

    free (octree->points_gnum);
    octree->points_gnum = _points_gnum;

    PDM_morton_code_t *_points_code =
      malloc (sizeof(PDM_morton_code_t) * octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {
      _points_code[i].L = octree->points_code[order[i]].L;
      _points_code[i].X[0] = octree->points_code[order[i]].X[0];
      _points_code[i].X[1] = octree->points_code[order[i]].X[1];
      _points_code[i].X[2] = octree->points_code[order[i]].X[2];
    }

    free (octree->points_code);
    octree->points_code = _points_code;

    double *_points = malloc (sizeof(double) * dim * octree->n_points);
    for (int i = 0; i < octree->n_points; i++) {
      for (int j = 0; j < dim; j++) {
        _points[dim*i+j] = octree->points[dim*order[i]+j];
      }
    }

    free (octree->points);
    octree->points = _points;

    free (order);
  }


  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_ORDER_POINTS] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_ORDER_POINTS]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_ORDER_POINTS]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_ORDER_POINTS]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  if (n_ranks > 1) {

    /*************************************************************************
     *
     * Store points in the octants (leaves) at the maximum depth of the octree
     * to build
     *
     *************************************************************************/

    int chg_code = 1;
    _l_octant_t *point_octants = malloc(sizeof(_l_octant_t));

    int curr_node = -1;

    _octants_init (point_octants, octree->dim, octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {

      PDM_morton_code_t _point_code;
      PDM_morton_copy (octree->points_code[i], &_point_code);

      PDM_morton_assign_level (&_point_code, octree->depth_max);

      if (curr_node != -1) {
        chg_code = !(PDM_morton_a_eq_b (point_octants->codes[curr_node],
                                        _point_code));
      }

      if (chg_code) {

        _octants_check_alloc (point_octants, 1);

        int idx = point_octants->n_nodes;

        curr_node = idx;

        PDM_morton_copy (octree->points_code[i], &(point_octants->codes[idx]));

        point_octants->n_points[idx] = 1;
        point_octants->range[idx] = i;

        point_octants->n_nodes += 1;
      }

      else {
        point_octants->n_points[curr_node] += 1;
      }

    }

    /*************************************************************************
     *
     * Block partition (algo 2 sundar)
     *
     *************************************************************************/

    octree->octants = _block_partition (point_octants,
                                        octree->comm,
                                        &octree->rank_octants_index);

    _octants_free (point_octants);

    /*************************************************************************
     *
     * Redistribute points
     *
     *************************************************************************/

    _distribute_points (&octree->n_points,
                        &octree->points,
                        &octree->points_icloud,
                        &octree->points_gnum,
                        &octree->points_code,
                        octree->rank_octants_index,
                        octree->comm,
                        octree->dim,
                        max_level,
                        octree->global_extents);

    int iblock = 0;
    for (int i = 0; i < octree->n_points; i++) {
      while (!PDM_morton_ancestor_is (octree->octants->codes[iblock],
                                      octree->points_code[i])) {
        iblock++;
      }
      assert (iblock < octree->octants->n_nodes);
      octree->octants->n_points[iblock] += 1;
    }

    octree->octants->range[0] = 0;

    double vol = 0;
    for (int i = 0; i < octree->octants->n_nodes; i++) {
      double side = 1./(double) (1 << octree->octants->codes[i].L);
      vol += (side * side * side);
      octree->octants->range[i+1] =
        octree->octants->range[i] +
        octree->octants->n_points[i];
    }
    double total_vol;
    PDM_MPI_Allreduce(&vol, &total_vol, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, octree->comm);

    if ( (PDM_ABS(total_vol - 1.)>= 1e-15)) {
      printf("Erreur volume different de 1 : %12.5e\n", total_vol);
    }

    assert (PDM_ABS(total_vol - 1.) < 1e-15);

  }

  else {

    octree->octants = malloc(sizeof(_l_octant_t));

    _octants_init (octree->octants, octree->dim, octree->n_points);

    PDM_morton_code_t code;

    code.L = 0;
    code.X[0] = 0;
    code.X[1] = 0;
    code.X[2] = 0;

    _octants_push_back (octree->octants,
                        code,
                        octree->n_points,
                        0);

    octree->octants->range[0] = 0;
    octree->octants->range[1] = octree->n_points;
    octree->octants->n_points[0] = octree->n_points;

  }

  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_BLOCK_PARTITION] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_BLOCK_PARTITION]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_BLOCK_PARTITION]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_BLOCK_PARTITION]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  /*************************************************************************
   *
   * Build local octree
   *
   *************************************************************************/
  const int n_child = 8;
  //const int n_direction = (int) PDM_N_DIRECTION;

  int  size = octree->depth_max * 8;

  //long mem = 0;
  _neighbours_tmp_t *ngb_octree = NULL;
  _neighbours_tmp_t *ngb_heap   = NULL;
  const int n_coarse = octree->octants->n_nodes;
  const int init_s = 1;
  _neighbours_tmp_t parent_ngb;

  size_t s_ngb_octree;
  if (NGB_ON_THE_FLY) {

    if (octree->neighboursToBuild) {
      for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
        parent_ngb.n_neighbour[dir] = 0;
        parent_ngb.s_neighbour[dir] = init_s;
        parent_ngb.neighbours[dir] = malloc (sizeof(int) * parent_ngb.s_neighbour[dir]);
      }

      ngb_heap = malloc (sizeof(_neighbours_tmp_t) * size);
      for (int i = 0; i < size; i++) {
        for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
          ngb_heap[i].neighbours[dir] = NULL;
        }
      }

      s_ngb_octree = octree->octants->n_nodes_max;
      ngb_octree = malloc (sizeof(_neighbours_tmp_t) * s_ngb_octree);
    }
  }

  _heap_t *heap = _heap_create (size);
  for (int i = octree->octants->n_nodes - 1; i >= 0; i--) {
    int is_pushed = _heap_push (heap,
                                octree->octants->codes[i],
                                octree->octants->range[i],
                                octree->octants->n_points[i]);
    if (!is_pushed) {
      printf ("Internal error PDM_para_octree 3 : heap is full\n");
      exit(1);
    }

    if (NGB_ON_THE_FLY) {
      if (octree->neighboursToBuild) {
        int h = heap->top - 1;

        /* init ngb_heap[h] */
        for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
          ngb_heap[h].n_neighbour[dir] = 0;
          ngb_heap[h].s_neighbour[dir] = init_s;
          ngb_heap[h].neighbours[dir] = malloc (sizeof(int) * ngb_heap[h].s_neighbour[dir]);
        }

        if (h == 0) {
          continue;
        }

        for (PDM_para_octree_direction_t dir = PDM_UP; dir < 6; dir+= (PDM_para_octree_direction_t)2) {
          PDM_para_octree_direction_t inv_dir = _inv_direction (dir);

          PDM_morton_code_t *ngb_code = _neighbour (heap->codes[h],
                                                    dir);
          if (ngb_code == NULL) {
            continue;
          }

          size_t start, end;
          PDM_morton_list_intersect (n_coarse - (i+1),
                                     *ngb_code,
                                     octree->octants->codes + (i+1),
                                     &start,
                                     &end);
          free (ngb_code);

          if (start >= end) {
            continue;
          }

          start += i+1;
          end   += i+1;

          size_t tmp = start;
          start = n_coarse - end;
          end = n_coarse - tmp;

          for (int j = start; j < (int) end; j++) {
            PDM_morton_code_t *ngb_ngb_code = _neighbour (heap->codes[j], inv_dir);
            assert (ngb_ngb_code != NULL);

            if (PDM_morton_ancestor_is (heap->codes[h], *ngb_ngb_code) ||
                PDM_morton_ancestor_is (*ngb_ngb_code, heap->codes[h])) {
              /* add -(j+1) to ngb_heap[h].neighbours[dir] */
              if (ngb_heap[h].n_neighbour[dir] >= ngb_heap[h].s_neighbour[dir]) {
                ngb_heap[h].s_neighbour[dir] = PDM_MAX (2*ngb_heap[h].s_neighbour[dir],
                                                        ngb_heap[h].n_neighbour[dir] + 1);
                ngb_heap[h].neighbours[dir] = realloc (ngb_heap[h].neighbours[dir],
                                                       sizeof(int) * ngb_heap[h].s_neighbour[dir]);
              }
              ngb_heap[h].neighbours[dir][ngb_heap[h].n_neighbour[dir]++] = -(j+1);

              /* add -(h+1) to ngb_heap[j].neighbours[inv_dir] */
              if (ngb_heap[j].n_neighbour[inv_dir] >= ngb_heap[j].s_neighbour[inv_dir]) {
                ngb_heap[j].s_neighbour[inv_dir] = PDM_MAX (2*ngb_heap[j].s_neighbour[inv_dir],
                                                            ngb_heap[j].n_neighbour[inv_dir] + 1);
                ngb_heap[j].neighbours[inv_dir] = realloc (ngb_heap[j].neighbours[inv_dir],
                                                           sizeof(int) * ngb_heap[j].s_neighbour[inv_dir]);
              }
              ngb_heap[j].neighbours[inv_dir][ngb_heap[j].n_neighbour[inv_dir]++] = -(h+1);
            }
            free (ngb_ngb_code);
          }
        }
      }
    }
  }


  PDM_morton_code_t code;
  int range;
  int n_points;

  octree->octants->n_nodes = 0;

  while (_heap_pull (heap, &code, &range, &n_points)) {

    /* Add children into the heap*/

    if ((code.L < max_morton_level) && (code.L < max_level) &&
        (n_points > octree->points_in_leaf_max)) {

      if (NGB_ON_THE_FLY) {
        if (octree->neighboursToBuild) {
          int h_parent = heap->top;

          /* copy ngb_heap[h_parent] into parent_ngb
             and remove references to -(h_parent+1) from all neighbours of parent */
          for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
            /* copy ngb_heap[h_parent].neighbours[dir] into parent_ngb.neighbours[dir] */
            parent_ngb.n_neighbour[dir] = ngb_heap[h_parent].n_neighbour[dir];
            if (parent_ngb.n_neighbour[dir] >= parent_ngb.s_neighbour[dir]) {
              /*parent_ngb.s_neighbour[dir] = PDM_MAX (2*parent_ngb.s_neighbour[dir],
                parent_ngb.n_neighbour[dir]);*/
              parent_ngb.s_neighbour[dir] = parent_ngb.n_neighbour[dir];

              parent_ngb.neighbours[dir] = realloc (parent_ngb.neighbours[dir],
                                                    sizeof(int) * parent_ngb.s_neighbour[dir]);
            }
            for (int j = 0; j < parent_ngb.n_neighbour[dir]; j++) {
              parent_ngb.neighbours[dir][j] = ngb_heap[h_parent].neighbours[dir][j];
            }

            /* remove all references to -(h_parent+1) from all neighbours of parent */
            PDM_para_octree_direction_t inv_dir = _inv_direction (dir);
            _neighbours_tmp_t *ngb = NULL;
            for (int j = 0; j < parent_ngb.n_neighbour[dir]; j++) {
              int ingb = parent_ngb.neighbours[dir][j];

              if (ingb < 0) {
                // neighbour in heap
                ngb = ngb_heap - (ingb+1);
              } else {
                // neighbour in octree
                ngb = ngb_octree + ingb;
              }

              int found = 0;
              int pos;
              for (pos = 0; pos < ngb->n_neighbour[inv_dir]; pos++) {
                if (ngb->neighbours[inv_dir][pos] == -(h_parent+1)) {
                  found = 1;
                  break;
                }
              }

              assert (found);

              ngb->n_neighbour[inv_dir]--;
              if (pos != ngb->n_neighbour[inv_dir]) {
                ngb->neighbours[inv_dir][pos] = ngb->neighbours[inv_dir][ngb->n_neighbour[inv_dir]];
              }

            }
          }
        }
      }

      PDM_morton_code_t children[n_child];
      PDM_morton_get_children(dim,
                              code,
                              children);

      int range_children[n_child];
      int n_points_children[n_child];

      for (int i = 0; i < n_child; i++) {
        n_points_children[i] = 0;
      }

      int ichild = 0;
      for (int i = 0; i < n_points; i++) {
        assert ((range + i) < octree->n_points);
        if (!PDM_morton_ancestor_is(code, octree->points_code[range + i])) {
          printf("Erreur : n'est pas un ancetre !!!!!\n");

        }
        assert (PDM_morton_ancestor_is(code, octree->points_code[range + i]));
        while (!PDM_morton_ancestor_is (children[ichild], octree->points_code[range + i])) {
          ichild += 1;
        }
        assert (ichild < n_child);
        n_points_children[ichild] += 1;
      }

      range_children[0] = 0;
      for (int i = 0; i < n_child - 1; i++) {
        range_children[i+1] = range_children[i] + n_points_children[i];
      }

      for (int i = n_child - 1; i >= 0; i--) {
        int is_pushed = _heap_push (heap,
                                    children[i],
                                    range + range_children[i],
                                    n_points_children[i]);
        if (!is_pushed) {
          printf ("Internal error PDM_para_octree 4 : heap is full\n");
          exit(1);
        }

        if (NGB_ON_THE_FLY) {
          if (octree->neighboursToBuild) {
            int h_child = heap->top - 1;

            for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
              /* reset ngb_heap[h_child].neighbours[dir] */
              ngb_heap[h_child].n_neighbour[dir] = 0;
              if (ngb_heap[h_child].neighbours[dir] == NULL) {
                ngb_heap[h_child].s_neighbour[dir] = PDM_MAX (init_s,
                                                              parent_ngb.n_neighbour[dir]);
                ngb_heap[h_child].neighbours[dir] = malloc (sizeof(int) * ngb_heap[h_child].s_neighbour[dir]);
              } else {
                if (parent_ngb.n_neighbour[dir] >= ngb_heap[h_child].s_neighbour[dir]) {
                  /*ngb_heap[h_child].s_neighbour[dir] = PDM_MAX (2*ngb_heap[h_child].s_neighbour[dir],
                    parent_ngb.n_neighbour[dir]);*/
                  ngb_heap[h_child].s_neighbour[dir] = parent_ngb.n_neighbour[dir];

                  ngb_heap[h_child].neighbours[dir] = realloc (ngb_heap[h_child].neighbours[dir],
                                                               sizeof(int) * ngb_heap[h_child].s_neighbour[dir]);
                }
              }


              /* set neighbours */
              int i_sibling = _3d_sibling_neighbours[i][dir];

              if (i_sibling < 0) {
                /* inherit neighbours from parent */
                PDM_para_octree_direction_t inv_dir = _inv_direction (dir);
                _neighbours_tmp_t *ngb = NULL;
                PDM_morton_code_t *ngb_code = NULL;
                for (int j = 0; j < parent_ngb.n_neighbour[dir]; j++) {
                  int ingb = parent_ngb.neighbours[dir][j];

                  if (ingb < 0) {
                    // neighbour in heap
                    ngb = ngb_heap - (ingb+1);
                    ngb_code = heap->codes - (ingb+1);
                  } else {
                    // neighbour in octree
                    ngb = ngb_octree + ingb;
                    ngb_code = octree->octants->codes + ingb;
                  }

                  /* check if current neighbour of parent is also a neighbour of current child */
                  PDM_morton_code_t *ngb_ngb_code = _neighbour (*ngb_code, inv_dir);
                  assert (ngb_ngb_code != NULL);

                  if (PDM_morton_ancestor_is (children[i], *ngb_ngb_code) ||
                      PDM_morton_ancestor_is (*ngb_ngb_code, children[i])) {
                    /* add ingb to ngb_heap[h_child].neighbours[dir] */
                    ngb_heap[h_child].neighbours[dir][ngb_heap[h_child].n_neighbour[dir]++] = ingb;
                    //printf("[%d] append %d to ngb_heap[%d].neighbours[%d] (neighbour of parent)\n", rank, ingb, h_child, dir);

                    /* add -(h_child+1) to ngb->neighbours[inv_dir] */
                    if (ngb->n_neighbour[inv_dir] >= ngb->s_neighbour[inv_dir]) {
                      /*ngb->s_neighbour[inv_dir] = PDM_MAX (2*ngb->s_neighbour[inv_dir],
                        ngb->n_neighbour[inv_dir]);*/
                      ngb->s_neighbour[inv_dir] = ngb->n_neighbour[inv_dir] + 1;

                      ngb->neighbours[inv_dir] = realloc (ngb->neighbours[inv_dir],
                                                          sizeof(int) * ngb->s_neighbour[inv_dir]);
                    }
                    ngb->neighbours[inv_dir][ngb->n_neighbour[inv_dir]++] = -(h_child+1);
                  }
                  free (ngb_ngb_code);
                }
              }
              else {
                /* add sibling neighbour */
                int h_sibling = h_child + i - i_sibling;
                ngb_heap[h_child].neighbours[dir][ngb_heap[h_child].n_neighbour[dir]++] = -(h_sibling+1);
              }

            }
          }
        }
      }
    }

    /* Store the leaf */

    else {
      _octants_push_back (octree->octants, code, n_points, range);

      if (NGB_ON_THE_FLY) {
        if (octree->neighboursToBuild) {
          int i = octree->octants->n_nodes - 1;
          int h = heap->top;

          if (i >= (int) s_ngb_octree) {
            s_ngb_octree = PDM_MAX (2*s_ngb_octree,
                                    (size_t) octree->octants->n_nodes_max);
            ngb_octree = realloc (ngb_octree, sizeof(_neighbours_tmp_t) * s_ngb_octree);
          }

          /* copy ngb_heap[h] into ngb_octree[i],
             and change references to -(h+1) into i for all neighbours */
          for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
            /* copy ngb_heap[h].neighbours[dir] into ngb_octree[i].neighbours[dir] */
            ngb_octree[i].n_neighbour[dir] = ngb_heap[h].n_neighbour[dir];
            //ngb_octree[i].s_neighbour[dir] = ngb_heap[h].s_neighbour[dir];
            ngb_octree[i].s_neighbour[dir] = PDM_MAX (init_s,
                                                      ngb_heap[h].n_neighbour[dir]);

            ngb_octree[i].neighbours[dir] = malloc (sizeof(int) * ngb_octree[i].s_neighbour[dir]);
            for (int j = 0; j < ngb_heap[h].n_neighbour[dir]; j++) {
              ngb_octree[i].neighbours[dir][j] = ngb_heap[h].neighbours[dir][j];
            }

            /* change references to -(h+1) into i for all neighbours */
            PDM_para_octree_direction_t inv_dir = _inv_direction (dir);
            _neighbours_tmp_t *ngb = NULL;
            for (int j = 0; j < ngb_octree[i].n_neighbour[dir]; j++) {
              int ingb = ngb_octree[i].neighbours[dir][j];

              if (ingb < 0) {
                // neighbour in heap
                ngb = ngb_heap - (ingb+1);
              } else {
                // neighbour in octree
                ngb = ngb_octree + ingb;
              }

              int found = 0;
              int pos;
              for (pos = 0; pos < ngb->n_neighbour[inv_dir]; pos++) {
                if (ngb->neighbours[inv_dir][pos] == -(h+1)) {
                  found = 1;
                  break;
                }
              }

              assert (found);

              ngb->neighbours[inv_dir][pos] = i;
            }
          }
        }
      }
    }
  }

  if (NGB_ON_THE_FLY) {
    if (octree->neighboursToBuild) {
      ngb_octree = realloc (ngb_octree, sizeof(_neighbours_tmp_t) * octree->octants->n_nodes);

      for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
        free (parent_ngb.neighbours[dir]);
      }

      for (int i = 0; i < heap->max_top; i++) {
        for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
          if (ngb_heap[i].neighbours[dir] != NULL) {
            free (ngb_heap[i].neighbours[dir]);
          }
        }
      }
      free (ngb_heap);
    }
  }

  double vol = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    double _side = 1./ (double) (1 << octree->octants->codes[i].L);
    vol += (_side * _side * _side);
  }
  double total_vol;
  PDM_MPI_Allreduce(&vol, &total_vol, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, octree->comm);

  assert (PDM_ABS(total_vol - 1.) < 1e-15);

  heap = _heap_free (heap);

  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NODES] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NODES]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NODES]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NODES]   += e_t_cpu_s - b_t_cpu_s;

  PDM_timer_resume(octree->timer);

  PDM_timer_hang_on(octree->timer);
  b_t_elapsed = PDM_timer_elapsed(octree->timer);
  b_t_cpu     = PDM_timer_cpu(octree->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);
  PDM_timer_resume(octree->timer);

  /*************************************************************************
   *
   * Neighbours
   *
   *************************************************************************/
  if (octree->neighboursToBuild) {
    if (NGB_ON_THE_FLY) {
      _finalize_neighbours (octree,
                            &ngb_octree,
                            b_t_elapsed,
                            b_t_cpu,
                            b_t_cpu_u,
                            b_t_cpu_s);
    }
    else {
      _compute_neighbours (octree,
        b_t_elapsed,
        b_t_cpu,
        b_t_cpu_u,
        b_t_cpu_s);
    }

    if (1 == 0) {
      _check_neighbours_area (octree);
    }
  }

  PDM_timer_hang_on(octree->timer);
  octree->times_elapsed[BUILD_TOTAL] = PDM_timer_elapsed(octree->timer);
  octree->times_cpu[BUILD_TOTAL]     = PDM_timer_cpu(octree->timer);
  octree->times_cpu_u[BUILD_TOTAL]   = PDM_timer_cpu_user(octree->timer);
  octree->times_cpu_s[BUILD_TOTAL]   = PDM_timer_cpu_sys(octree->timer);
  PDM_timer_resume(octree->timer);

  octree->times_elapsed[END] = octree->times_elapsed[BUILD_TOTAL];
  octree->times_cpu[END]     = octree->times_cpu[BUILD_TOTAL];
  octree->times_cpu_u[END]   = octree->times_cpu_u[BUILD_TOTAL];
  octree->times_cpu_s[END]   = octree->times_cpu_s[BUILD_TOTAL];

}


/**
 *
 * \brief Get extents
 *
 * \param [in]   id                 Identifier
 *
 * \return     Extents
 *
 */

double *
PDM_para_octree_extents_get
(
 const int  id
 )
{
  _octree_t *octree = _get_from_id (id);

  return octree->global_extents;
}

/**
 *
 * \brief Dump octree
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_dump
(
 const int          id
 )
{
  _octree_t *octree = _get_from_id (id);
  PDM_printf("PDM_dump_para_octree : %d\n",id);

  PDM_printf("  - n_nodes : %d\n", octree->octants->n_nodes);
  printf("  - global_extents :");
  for (int i = 0; i < 2*octree->dim; i++) {
    printf(" %12.5e", octree->global_extents[i]);
  }
  PDM_printf("\n");
  PDM_printf("  - depth_max : %d\n", octree->depth_max);
  PDM_printf("  - points_in_leaf_max : %d\n", octree->points_in_leaf_max);

  PDM_printf("  - s : %12.5e %12.5e %12.5e\n", octree->s[0], octree->s[1], octree->s[2]);
  PDM_printf("  - d : %12.5e %12.5e %12.5e\n", octree->d[0], octree->d[1], octree->d[2]);

  PDM_printf("  - n_point_clouds : %d\n", octree->n_point_clouds);
  PDM_printf("  - t_n_points : "PDM_FMT_G_NUM"\n", octree->t_n_points);
  PDM_printf("  - n_points : %d\n", octree->n_points);
  for (int i = 0; i < octree->n_points; i++) {
    PDM_printf("  %d "PDM_FMT_G_NUM" : %12.5e %12.5e %12.5e / level %u code [%u, %u, %u]\n",
               i, octree->points_gnum[i],
               octree->points[3*i], octree->points[3*i+1], octree->points[3*i+2],
               octree->points_code[i].L,
               octree->points_code[i].X[0],
               octree->points_code[i].X[1],
               octree->points_code[i].X[2]);
  }
  PDM_printf("  - n_nodes : %d\n", octree->octants->n_nodes);
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    PDM_printf("  %d : level %u code [%u, %u, %u], range %d, n_points %d\n",
               i,
               octree->octants->codes[i].L,
               octree->octants->codes[i].X[0],
               octree->octants->codes[i].X[1],
               octree->octants->codes[i].X[2],
               octree->octants->range[i],
               octree->octants->n_points[i]
               );
  }
  if (octree->neighboursToBuild) {
    for (int i = 0; i < octree->octants->n_nodes; i++) {
      PDM_printf("  %d : neighbors\n", i);
      for (int j = 0; j < 6; j++) {
        PDM_printf("    - direction %d : ", j);
        for (int k = octree->octants->neighbour_idx[6*i+j];
             k < octree->octants->neighbour_idx[6*i+j+1]; k++) {
          PDM_printf(" %d",octree->octants->neighbours[k]);
        }
        PDM_printf("\n");
      }
    }
  }
}


/**
 *
 * Look for closest points stored inside an octree
 *
 * \param [in]   id                     Identifier
 * \param [in]   n_closest_points       Number of closest points to find
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [in]   pts_g_num              Point global numbers
 * \param [out]  closest_octree_pt_id   Closest points in octree global number
 * \param [out]  closest_octree_pt_dist Closest points in octree distance
 *
 */


void
PDM_para_octree_closest_point
(
 const int    id,
 const int    n_closest_points,
 const int    n_pts,
 double      *pts,
 PDM_g_num_t *pts_g_num,
 PDM_g_num_t *closest_octree_pt_g_num,
 double      *closest_octree_pt_dist2
 )
{
  const int DEBUG = 0;
  const int DEBUG_FILTER = 0;
  const int DEBUG_MERGE = 0;

  const int COMPUTE_FIRST_UPPER_BOUND = 1;
  const double eps_upper_bound_dist = 1.e-6;

  const _octree_t *octree = _get_from_id (id);
  const _l_octant_t *octants = octree->octants;

  const int dim = octree->dim;
  int _n_closest_points = n_closest_points;
  int _n_pts = n_pts;

  int myRank, lComm;
  PDM_MPI_Comm_rank (octree->comm, &myRank);
  PDM_MPI_Comm_size (octree->comm, &lComm);


  /* /!\ /!\ /!\ Force target points inside octree extents /!\ /!\ /!\ -->> */
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < dim; j++) {
      pts[dim*i+j] = PDM_MAX (pts[dim*i+j], octree->global_extents[j]);
      pts[dim*i+j] = PDM_MIN (pts[dim*i+j], octree->global_extents[dim+j]);
    }
  }
  /* <<-- */


  /* Part-to-block create (only to get block distribution) */
  PDM_part_to_block_t *ptb1 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                        1.,
                                                        &pts_g_num,
                                                        NULL,
                                                        &_n_pts,
                                                        1,
                                                        octree->comm);

  PDM_g_num_t *block_distrib_idx1 = PDM_part_to_block_distrib_index_get (ptb1);
  const int n_pts_block1 = PDM_part_to_block_n_elt_block_get (ptb1);

  PDM_g_num_t *block_closest_src_gnum1 = malloc (sizeof(PDM_g_num_t) * n_pts_block1 * n_closest_points);
  double      *block_closest_src_dist1 = malloc (sizeof(double)      * n_pts_block1 * n_closest_points);

  for (int i = 0; i < n_pts_block1; i++) {
    for (int j = 0; j < n_closest_points; j++) {
      block_closest_src_dist1[n_closest_points*i + j] = HUGE_VAL;
      block_closest_src_gnum1[n_closest_points*i + j] = -1; // USEFUL ONLY FOR DEBUG
    }
  }



  /*************************************************************************
   *
   * Distribute the target points
   *
   *************************************************************************/
  /*   1) Encode the coordinates of every target point */
  PDM_morton_code_t *pts_code = malloc (sizeof(PDM_morton_code_t) * n_pts);
  double d[3], s[3];
  PDM_morton_encode_coords (dim,
                            PDM_morton_max_level,
                            octree->global_extents,
                            (size_t) n_pts,
                            pts,
                            pts_code,
                            d,
                            s);

  /*   2) Use binary search to associate each target point to the appropriate process */
  int *send_count = malloc (sizeof(int) * lComm);
  int *recv_count = malloc (sizeof(int) * lComm);
  for (int i = 0; i < lComm; i++) {
    send_count[i] = 0;
  }

  int *rank_pt = malloc (sizeof(int) * n_pts);
  for (int i = 0; i < n_pts; i++) {
    rank_pt[i] = PDM_morton_binary_search (lComm,
                                           pts_code[i],
                                           octree->rank_octants_index);
    send_count[rank_pt[i]]++;
  }
  free (pts_code);

  /*   3) Exchange send/recv counts */
  PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                    recv_count, 1, PDM_MPI_INT,
                    octree->comm);

  int *send_shift = malloc (sizeof(int) * (lComm+1));
  int *recv_shift = malloc (sizeof(int) * (lComm+1));
  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int i = 0; i < lComm; i++) {
    send_shift[i+1] = send_shift[i] + send_count[i];
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
    send_count[i] = 0;
  }

  /*   4) Fill send buffers */
  PDM_g_num_t *send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[lComm]);
  PDM_g_num_t *recv_g_num = malloc (sizeof(PDM_g_num_t) * recv_shift[lComm]);
  double      *send_coord = malloc (sizeof(double)      * send_shift[lComm]*dim);
  double      *recv_coord = malloc (sizeof(double)      * recv_shift[lComm]*dim);
  for (int i = 0; i < n_pts; i++) {
    int rank = rank_pt[i];
    int k = send_shift[rank] + send_count[rank];
    send_g_num[k] = pts_g_num[i];
    for (int j = 0; j < dim; j++)
      send_coord[dim*k+j] = pts[dim*i+j];

    send_count[rank]++;
  }
  free (rank_pt);

  /*   5) Send gnum buffer */
  PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                     recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                     octree->comm);

  /*   6) Send coord buffer */
  int n_recv_pts = recv_shift[lComm];
  for (int i = 0; i < lComm; i++) {
    send_count[i] *= dim;
    recv_count[i] *= dim;
    send_shift[i+1] *= dim;
    recv_shift[i+1] *= dim;
  }
  PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                     recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                     octree->comm);




  PDM_g_num_t *local_closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_recv_pts * n_closest_points);
  double      *local_closest_src_dist = malloc (sizeof(double)      * n_recv_pts * n_closest_points);
  double      *upper_bound_dist       = malloc (sizeof(double)      * n_recv_pts);

  /* Encode the coordinates of the received target points */
  pts_code = malloc (sizeof(PDM_morton_code_t) * n_recv_pts);
  PDM_morton_encode_coords (dim,
                            PDM_morton_max_level,
                            octree->global_extents,
                            (size_t) n_recv_pts,
                            recv_coord,
                            pts_code,
                            d,
                            s);

  if (COMPUTE_FIRST_UPPER_BOUND && octree->n_points >= n_closest_points) {
    const int window_width = (int) ceil (0.5 * n_closest_points);
    int window_start, window_end;
    /* Inspect a window of src points around each tgt point on the Z-order curve */
    for (int i = 0; i < n_recv_pts; i++) {
      int pos = PDM_morton_binary_search (octree->n_points,
                                          pts_code[i],
                                          octree->points_code);

      if (pos < window_width) {
        window_start = 0;
        window_end   = n_closest_points;
      } else if (pos >= octree->n_points - window_width) {
        window_end   = octree->n_points;
        window_start = window_end - n_closest_points;
      } else {
        window_start = pos - window_width;
        window_end   = window_start + n_closest_points;
      }

      const double *_pt = recv_coord + i * dim;
      double max_dist = 0.;
      for (int i_src = window_start; i_src < window_end; i_src++) {
        double *src_pt = octree->points + i_src * dim;
        double src_dist = 0;
        for (int j = 0; j < dim; j++) {
          double delta = _pt[j] - src_pt[j];
          src_dist += delta * delta;
        }

        max_dist = PDM_MAX (max_dist, src_dist);
      }

      upper_bound_dist[i] = max_dist + eps_upper_bound_dist;
    }

  } else {

    for (int i = 0; i < n_recv_pts; i++) {
      upper_bound_dist[i] = HUGE_VAL;
    }

  }


  /* Find start leaf for each target point */
  int *start_leaves = malloc (sizeof(int) * n_recv_pts);
  for (int i = 0; i < n_recv_pts; i++) {
    start_leaves[i] = PDM_morton_binary_search (octants->n_nodes,
                                                pts_code[i],
                                                octants->codes);
  }
  free (pts_code);

  int *start_leaves_idx = malloc (sizeof(int) * (n_recv_pts+1));
  start_leaves_idx[0] = 0;
  for (int i = 0; i < n_recv_pts; i++) {
    start_leaves_idx[i+1] = start_leaves_idx[i] + 1;
  }


  /* Stuff used for filering data coming from concurrent ranks */
  /*PDM_hash_tab_t *processed_tgt = malloc (sizeof(PDM_hash_tab_t) * octree->n_connected);
    PDM_g_num_t keyMax = n_recv_pts; //?
    for (int i_part = 0; i_part < octree->n_connected; i_part++) {
    processed_tgt[i_part] = PDM_hash_tab_create (PDM_HASH_TAB_KEY_LONG
    &keyMax);
    }*/
  int *s_processed_tgt = malloc (sizeof(int) * octree->n_connected);
  int *n_processed_tgt = malloc (sizeof(int) * octree->n_connected);
  PDM_g_num_t **processed_tgt = malloc (sizeof(PDM_g_num_t*) * octree->n_connected);
  for (int i = 0; i < octree->n_connected; i++) {
    s_processed_tgt[i] = PDM_MAX (128, 2 * n_recv_pts);//?
    n_processed_tgt[i] = 0;
    processed_tgt[i] = malloc (sizeof(PDM_g_num_t) * s_processed_tgt[i]);
  }

  PDM_g_num_t **new_processed_tgt   = malloc (sizeof(PDM_g_num_t *) * octree->n_connected);
  int          *n_new_processed_tgt = malloc (sizeof(int) * octree->n_connected);
  int          *s_new_processed_tgt = malloc (sizeof(int) * octree->n_connected);
  for (int i = 0; i < octree->n_connected; i++) {
    s_new_processed_tgt[i] = PDM_MAX (128, n_recv_pts);//?
    new_processed_tgt[i] = malloc (sizeof(PDM_g_num_t) * s_new_processed_tgt[i]);
  }


  /* Stuff used to merge results received from 'part-to-block'2 */
  _min_heap_t *merge_heap = _min_heap_create (10 * n_closest_points); // used to merge part-to-block results

  double      *tmp_closest_src_dist = malloc (sizeof(double)      * n_closest_points);
  PDM_g_num_t *tmp_closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_closest_points);

  int *block_stride2 = NULL;
  double      *block_closest_src_dist2 = NULL;
  PDM_g_num_t *block_closest_src_gnum2 = NULL;
  double      *block_upper_bound_dist = malloc (sizeof(double) * n_pts_block1);
  int one = 1;


  /* Stuff used to redistribute target points */
  int *send_tgt_lnum = NULL;
  int *send_start_leaves = NULL;
  int *send_start_leaves_count = NULL;
  int *send_start_leaves_rank_shift = NULL;

  // while loop...
  int iteration = 0;
  while (1) {
    //-->>
    iteration++;
    if (DEBUG) {
      printf("\n\n\n[%d] iteration %d\n", myRank, iteration);
    }
    else {
      if (1== 0 && myRank == 0)
        printf("\n\n\niteration %d\n", iteration);
    }
    if (iteration > 10*lComm) break; // emergency exit
    //<<--

    /* Filter 'recv' data: */
    //--->>>
    if (DEBUG) {
      printf("\n\n++++ BEFORE FILTERING ++++\n");
      printf("[%d] processed_pts:\n", myRank);
      for (int i = 0; i < octree->n_connected; i++) {
        printf("\tpart %d:", i);
        for (int j = 0; j < n_processed_tgt[i]; j++) {
          printf(" "PDM_FMT_G_NUM"", processed_tgt[i][j]);
        }
        printf("\n");
      }


      printf("\n[%d] recv_g_num = [\n", myRank);
      for (int i = 0; i < lComm; i++) {
        printf("\t{rank %d:", i);
        for (int j = recv_shift[i]/dim; j < recv_shift[i+1]/dim; j++) {
          printf(" "PDM_FMT_G_NUM"", recv_g_num[j]);
        }
        printf(" }\n");
      }
      printf(" ]\n");


      printf("\n[%d] start_leaves = [\n", myRank);
      for (int i = 0; i < n_recv_pts; i++) {
        printf("  ("PDM_FMT_G_NUM") : {", recv_g_num[i]);
        for (int j = start_leaves_idx[i]; j < start_leaves_idx[i+1]; j++) {
          printf(" %d", start_leaves[j]);
        }
        printf("}\n");
      }
      printf(" ]\n\n\n");
    }
    //<<<---

    if (iteration == 1) {
      // at this point, there is exactly one start leaf per tgt point
      for (int i = 0; i < n_recv_pts; i++) {
        // get part number
        int i_part = _get_octant_part_id (octree,
                                          start_leaves[i]);
        processed_tgt[i_part][n_processed_tgt[i_part]++] = recv_g_num[i];
      }

      // sort in ascending order
      for (int i_part = 0; i_part < octree->n_connected; i_part++) {
        PDM_sort_long (processed_tgt[i_part],
                       NULL,
                       n_processed_tgt[i_part]);
      }

    } else {
      // at this point, there are possibly
      //   1) duplicate tgt points (but with different start leaves)
      //   2) tgt points already processed in some parts of current rank
      int tmp_n = 0;
      PDM_g_num_t *tmp_g_num       = malloc (sizeof(PDM_g_num_t) * n_recv_pts);
      double      *tmp_coord       = malloc (sizeof(double)      * n_recv_pts * dim);
      double      *tmp_upper_bound = malloc (sizeof(double)      * n_recv_pts);

      int  *n_tmp_start_leaves = malloc (sizeof(int) * n_recv_pts);
      int  *s_tmp_start_leaves = malloc (sizeof(int) * n_recv_pts);
      int **tmp_start_leaves   = malloc (sizeof(int *) * n_recv_pts);


      for (int i_part = 0; i_part < octree->n_connected; i_part++) {
        n_new_processed_tgt[i_part] = 0;
      }

      for (int i_tgt = 0; i_tgt < n_recv_pts; i_tgt++) {
        PDM_g_num_t tgt_gnum = recv_g_num[i_tgt];
        if (DEBUG && DEBUG_FILTER) {
          printf("\nFilter: recv pt ("PDM_FMT_G_NUM")\n", tgt_gnum);
        }

        // check whether that tgt point has already been added to tmp_* arrays
        int found_tmp = 0;
        int pos_tmp = _binary_search_long (tgt_gnum,
                                           tmp_g_num,
                                           tmp_n,
                                           &found_tmp);
        if (DEBUG && DEBUG_FILTER) {
          printf(" found_tmp = %d\tpos_tmp = %d\n", found_tmp, pos_tmp);
        }

        if (!found_tmp) {
          pos_tmp = tmp_n;
          // add tgt point to tmp_* arrays
          tmp_g_num[pos_tmp] = tgt_gnum;
          for (int i = 0; i < dim; i++) {
            tmp_coord[dim*pos_tmp + i] = recv_coord[dim*i_tgt + i];
          }

          tmp_upper_bound[pos_tmp] = upper_bound_dist[i_tgt];

          s_tmp_start_leaves[pos_tmp] = 2 * (start_leaves_idx[i_tgt+1] - start_leaves_idx[i_tgt]);
          tmp_start_leaves[pos_tmp] = malloc (sizeof(int) * s_tmp_start_leaves[pos_tmp]);
          n_tmp_start_leaves[pos_tmp] = 0;

          tmp_n++;
        }

        for (int i_start = start_leaves_idx[i_tgt];
             i_start < start_leaves_idx[i_tgt+1]; i_start++) {
          int leaf_id = start_leaves[i_start];
          int leaf_part = _get_octant_part_id (octree,
                                               leaf_id);
          if (DEBUG && DEBUG_FILTER) {
            printf("\tstart leaf id %d (part %d)\n", leaf_id, leaf_part);
          }

          // check whether that tgt point has already been processed in part #leaf_part
          int found = 0;
          _binary_search_long (tgt_gnum,
                               processed_tgt[leaf_part],
                               n_processed_tgt[leaf_part],
                               &found);
          if (DEBUG && DEBUG_FILTER) {
            printf("\t already processed? %d\n", found);
          }

          if (!found) {
            // check whether that start leaf has already been added to tmp_start_leaves[leaf_part]
            int found_leaf = 0;
            int pos_leaf = _binary_search (leaf_id,
                                           tmp_start_leaves[pos_tmp],
                                           n_tmp_start_leaves[pos_tmp],
                                           &found_leaf);
            if (DEBUG && DEBUG_FILTER) {
              printf("\t\tfound_leaf = %d, pos_leaf = %d\n", found_leaf, pos_leaf);
            }

            if (!found_leaf) {
              /* add start leaf to tmp_start_leaves[pos_tmp] */
              // realloc tmp_start_leaves[pos_tmp] if necessary
              if (s_tmp_start_leaves[pos_tmp] <= n_tmp_start_leaves[pos_tmp]) {
                s_tmp_start_leaves[pos_tmp] *= 2;
                tmp_start_leaves[pos_tmp] = realloc (tmp_start_leaves[pos_tmp],
                                                     sizeof(int) * s_tmp_start_leaves[pos_tmp]);
              }

              // insert-sort leaf_id in tmp_start_leaves[pos_tmp]
              for (int i = n_tmp_start_leaves[pos_tmp]; i > pos_leaf; i--) {
                tmp_start_leaves[pos_tmp][i] = tmp_start_leaves[pos_tmp][i-1];
              }
              tmp_start_leaves[pos_tmp][pos_leaf] = leaf_id;
              n_tmp_start_leaves[pos_tmp]++;



              // check whether that tgt point has already been added to new_processed_tgt[leaf_part]
              int found_new = 0;
              int pos_new = _binary_search_long (tgt_gnum,
                                                 new_processed_tgt[leaf_part],
                                                 n_new_processed_tgt[leaf_part],
                                                 &found_new);
              if (DEBUG && DEBUG_FILTER) {
                printf("\t\tfound_new = %d, pos_new = %d\n", found_new, pos_new);
              }

              if (!found_new) {
                /* add tgt point to new_processed_tgt[leaf_part] */
                // realloc new_processed_tgt[leaf_part] if necessary
                if (s_new_processed_tgt[leaf_part] <= n_new_processed_tgt[leaf_part]) {
                  s_new_processed_tgt[leaf_part] *= 2;
                  new_processed_tgt[leaf_part] = realloc (new_processed_tgt[leaf_part],
                                                          sizeof(PDM_g_num_t) * s_new_processed_tgt[leaf_part]);
                }

                // insert-sort tgt_gnum in new_processed_tgt[leaf_part]
                for (int i = n_new_processed_tgt[leaf_part]; i > pos_new; i--) {
                  new_processed_tgt[leaf_part][i] = new_processed_tgt[leaf_part][i-1];
                }
                new_processed_tgt[leaf_part][pos_new] = tgt_gnum;
                n_new_processed_tgt[leaf_part]++;
              }

            } // end if (leaf_id not found in processed_tgt[leaf_part])
          } // end if (tgt_gnum not found in processed_tgt[leaf_part])
        } // end loop over start leaves (i_start)
      } // end loop over received tgt points (i_tgt)

      int k = 0;
      start_leaves_idx[0] = 0;
      for (int i = 0; i < tmp_n; i++) {
        if (n_tmp_start_leaves[i] > 0) {
          start_leaves_idx[k+1] = start_leaves_idx[k] + n_tmp_start_leaves[i];

          recv_g_num[k] = tmp_g_num[i];

          for (int j = 0; j < dim; j++) {
            recv_coord[dim*k+j] = tmp_coord[dim*i+j];
          }

          upper_bound_dist[k] = tmp_upper_bound[i];
          k++;
        }
      }

      if (k < n_recv_pts) {
        recv_g_num       = realloc (recv_g_num,       sizeof(PDM_g_num_t) * k);
        recv_coord       = realloc (recv_coord,       sizeof(double)      * k * dim);
        upper_bound_dist = realloc (upper_bound_dist, sizeof(double)      * k);
      }
      free (tmp_g_num);
      free (tmp_coord);
      free (tmp_upper_bound);

      /* manage start leaves */
      start_leaves = realloc (start_leaves, sizeof(int) * start_leaves_idx[k]);
      k = 0;
      int idx = 0;
      for (int i = 0; i < tmp_n; i++) {
        if (n_tmp_start_leaves[i] > 0) {
          for (int j = 0; j < n_tmp_start_leaves[i]; j++) {
            start_leaves[idx++] = tmp_start_leaves[i][j];
          }
          k++;
        }
      }

      for (int i = 0; i < tmp_n; i++) {
        free (tmp_start_leaves[i]);
      }
      free (tmp_start_leaves);
      free (s_tmp_start_leaves);
      free (n_tmp_start_leaves);

      n_recv_pts = k;


      /* merge new_processed_tgt into processed_tgt */
      for (int i_part = 0; i_part < octree->n_connected; i_part++) {
        // realloc processed_tgt[i_part] if necessary
        if (s_processed_tgt[i_part] <= n_processed_tgt[i_part] + n_new_processed_tgt[i_part]) {
          s_processed_tgt[i_part] = PDM_MAX (2 * s_processed_tgt[i_part],
                                             n_processed_tgt[i_part] + n_new_processed_tgt[i_part]);
          processed_tgt[i_part] = realloc (processed_tgt[i_part],
                                           sizeof(PDM_g_num_t) * s_processed_tgt[i_part]);
        }

        for (int i = 0; i < n_new_processed_tgt[i_part]; i++) {
          int found = 0;
          int pos = _binary_search_long (new_processed_tgt[i_part][i],
                                         processed_tgt[i_part],
                                         n_processed_tgt[i_part],
                                         &found);

          // insert-sort
          for (int j = n_processed_tgt[i_part]; j > pos; j--) {
            processed_tgt[i_part][j] = processed_tgt[i_part][j-1];
          }
          processed_tgt[i_part][pos] = new_processed_tgt[i_part][i];
          n_processed_tgt[i_part]++;
        }
      }

    } // end if/else (iteration == 1)


    //--->>>
    if (DEBUG) {
      printf("++++ AFTER FILTERING ++++\n");
      printf("[%d] processed_pts:\n", myRank);
      for (int i = 0; i < octree->n_connected; i++) {
        printf("\tpart %d:", i);
        for (int j = 0; j < n_processed_tgt[i]; j++) {
          printf(" "PDM_FMT_G_NUM"", processed_tgt[i][j]);
        }
        printf("\n");
      }


      printf("\n[%d] recv_g_num = [", myRank);
      for (int i = 0; i < n_recv_pts; i++) {
        printf(" "PDM_FMT_G_NUM"", recv_g_num[i]);
      }
      printf(" ]\n");


      printf("\n[%d] start_leaves = [\n", myRank);
      for (int i = 0; i < n_recv_pts; i++) {
        printf("  ("PDM_FMT_G_NUM") : {", recv_g_num[i]);
        for (int j = start_leaves_idx[i]; j < start_leaves_idx[i+1]; j++) {
          printf(" %d", start_leaves[j]);
        }
        printf("}\n");
      }
      printf(" ]\n\n\n");
    }
    //<<<---

    /* Search closest src points in local octree */
    _closest_points_local (octree,
                           n_closest_points,
                           n_recv_pts,
                           recv_coord,
                           recv_g_num, // ONLY FOR DEBUG
                           start_leaves,
                           start_leaves_idx,
                           upper_bound_dist,
                           local_closest_src_gnum,
                           local_closest_src_dist,
                           send_count,
                           send_shift,
                           &send_tgt_lnum,
                           &send_start_leaves,
                           &send_start_leaves_count,
                           &send_start_leaves_rank_shift);

    /* Part-to-block exchanges to merge results in block arrays */
    PDM_part_to_block_t *ptb2 = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           &recv_g_num,
                                                           block_distrib_idx1,
                                                           &n_recv_pts,
                                                           1,
                                                           octree->comm);

    int n_pts_block2 = PDM_part_to_block_n_elt_block_get (ptb2);
    PDM_g_num_t *block_tgt_gnum2 = PDM_part_to_block_block_gnum_get (ptb2);

    int *stride2 = malloc (sizeof(int) * n_recv_pts);
    for (int i = 0; i < n_recv_pts; i++) {
      stride2[i] = n_closest_points;
    }

    PDM_part_to_block_exch (ptb2,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &stride2,
                            (void **) &local_closest_src_dist,
                            &block_stride2,
                            (void **) &block_closest_src_dist2);
    free (block_stride2);

    PDM_part_to_block_exch (ptb2,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            1,
                            &stride2,
                            (void **) &local_closest_src_gnum,
                            &block_stride2,
                            (void **) &block_closest_src_gnum2);
    free (stride2);

    /* Merge block data */
    if (DEBUG && DEBUG_MERGE) {
      printf("\n\n- - - Merge - - -");
    }
    int *block_idx2 = malloc (sizeof(int) * (n_pts_block2 + 1));
    block_idx2[0] = 0;
    for (int i = 0; i < n_pts_block2; i++) {
      block_idx2[i+1] = block_idx2[i] + block_stride2[i];
    }

    for (int i = 0; i < n_pts_block2; i++) {
      int id_block1 = (int) n_closest_points * (block_tgt_gnum2[i]-1 - block_distrib_idx1[myRank]);

      if (DEBUG && DEBUG_MERGE) {
        /*printf("\nblock2 i = %d ("PDM_FMT_G_NUM") <%d>\n", i, block_tgt_gnum2[i],
          (int) (block_tgt_gnum2[i]-1 - block_distrib_idx1[myRank]));
          printf("id_block1 = %d / %d\n", id_block1, n_pts_block1*n_closest_points);*/
        printf("\npoint ("PDM_FMT_G_NUM")\n", block_tgt_gnum2[i]);
      }

      int n_procs_pt = block_stride2[i] / n_closest_points;

      int *idx_proc = malloc (sizeof(int) * (n_procs_pt + 1));

      _min_heap_reset (merge_heap);

      for (int j = 0; j < n_closest_points; j++) {
        tmp_closest_src_gnum[j] = block_closest_src_gnum1[id_block1 + j];
        tmp_closest_src_dist[j] = block_closest_src_dist1[id_block1 + j];
      }

      if (DEBUG && DEBUG_MERGE) {
        printf("\tj = %d: ("PDM_FMT_G_NUM", %f)\n",
               -1, tmp_closest_src_gnum[0], tmp_closest_src_dist[0]);
      }

      _min_heap_push (merge_heap,
                      -1,
                      tmp_closest_src_gnum[0],
                      tmp_closest_src_dist[0]);
      idx_proc[0] = 0;


      for (int j = 0; j < n_procs_pt; j++) {
        int k = block_idx2[i] + j*n_closest_points;

        if (DEBUG && DEBUG_MERGE) {
          printf("\tj =  %d: ("PDM_FMT_G_NUM", %f)\n",
                 j, block_closest_src_gnum2[k], block_closest_src_dist2[k]);
        }

        _min_heap_push (merge_heap,
                        j,
                        block_closest_src_gnum2[k],
                        block_closest_src_dist2[k]);
        idx_proc[j+1] = 0;
      }


      int _proc;
      PDM_g_num_t _gnum;
      double _dist;
      for (int j = 0; j < n_closest_points; j++) {
        int popped = _min_heap_pop (merge_heap,
                                    &_proc,
                                    &_gnum,
                                    &_dist);
        assert (popped);
        if (DEBUG && DEBUG_MERGE) {
          printf("\tpopped: "PDM_FMT_G_NUM" %f \t / %f\n", _gnum, _dist, block_closest_src_dist1[id_block1 + j]);
        }

        if (_dist >= block_closest_src_dist1[id_block1 + n_closest_points - 1]) break;

        block_closest_src_dist1[id_block1 + j] = _dist;
        block_closest_src_gnum1[id_block1 + j] = _gnum;

        if (DEBUG && DEBUG_MERGE) {
          printf("\t\t%d/%d: "PDM_FMT_G_NUM", %f\n",
                 j+1, n_closest_points,
                 block_closest_src_gnum1[id_block1 + j], block_closest_src_dist1[id_block1 + j]);
        }

        if (j >= n_closest_points - 1)
          break;

        idx_proc[_proc+1]++;

        if (_proc < 0) {
          _min_heap_push (merge_heap,
                          -1,
                          tmp_closest_src_gnum[idx_proc[0]],
                          tmp_closest_src_dist[idx_proc[0]]);
        } else {
          int k = block_idx2[i] + _proc*n_closest_points + idx_proc[_proc+1];
          _min_heap_push (merge_heap,
                          _proc,
                          block_closest_src_gnum2[k],
                          block_closest_src_dist2[k]);
        }

      }
      free (idx_proc);
    }
    if (DEBUG && DEBUG_MERGE) {
      printf("- - - - - - - - -\n\n");
    }

    free (block_idx2);
    free (block_stride2);
    free (block_closest_src_gnum2);
    free (block_closest_src_dist2);
    ptb2 = PDM_part_to_block_free (ptb2);
    // end merge


    /* Update upper_bound_dist */
    PDM_block_to_part_t *btp2 = PDM_block_to_part_create (block_distrib_idx1,
                                                          (const PDM_g_num_t **) &recv_g_num,
                                                          &n_recv_pts,
                                                          1,
                                                          octree->comm);

    for (int i = 0; i < n_pts_block1; i++) {
      block_upper_bound_dist[i] = block_closest_src_dist1[n_closest_points*(i+1) - 1];
    }

    PDM_block_to_part_exch (btp2,
                            sizeof(double),
                            PDM_STRIDE_CST,
                            &one,
                            block_upper_bound_dist,
                            NULL,
                            (void **) &upper_bound_dist);
    btp2 = PDM_block_to_part_free (btp2);


    /* Redistribute target points for next iteration */
    // send count
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      octree->comm);
    recv_shift[0] = 0;
    for (int i = 0; i < lComm; i++)
      recv_shift[i+1] = recv_shift[i] + recv_count[i];

    n_recv_pts = recv_shift[lComm];


    /* Termination criterion */
    int max_n_recv_pts = 0;
    PDM_MPI_Allreduce (&n_recv_pts, &max_n_recv_pts, 1,
                       PDM_MPI_INT, PDM_MPI_MAX, octree->comm);

    if (max_n_recv_pts == 0) {
      free (send_tgt_lnum);
      free (send_start_leaves);
      free (send_start_leaves_count);
      free (send_start_leaves_rank_shift);
      break;
    }


    // send g_num
    send_g_num = realloc (send_g_num, sizeof(PDM_g_num_t) * send_shift[lComm]);
    for (int i = 0; i < send_shift[lComm]; i++) {
      send_g_num[i] = recv_g_num[send_tgt_lnum[i]];
    }
    recv_g_num = realloc (recv_g_num, sizeof(PDM_g_num_t) * recv_shift[lComm]);
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       octree->comm);

    // send upper bound dist
    double *send_upper_bound_dist = malloc (sizeof(double) * send_shift[lComm]);
    for (int i = 0; i < send_shift[lComm]; i++) {
      send_upper_bound_dist[i] = upper_bound_dist[send_tgt_lnum[i]];
    }
    upper_bound_dist = realloc (upper_bound_dist, sizeof(double) * n_recv_pts);
    PDM_MPI_Alltoallv (send_upper_bound_dist, send_count, send_shift, PDM_MPI_DOUBLE,
                       upper_bound_dist,      recv_count, recv_shift, PDM_MPI_DOUBLE,
                       octree->comm);
    free (send_upper_bound_dist);

    // send start leaves
    int *send_start_leaves_rank_count = malloc (sizeof(int) * lComm);
    int *recv_start_leaves_rank_count = malloc (sizeof(int) * lComm);
    for (int i = 0; i < lComm; i++) {
      send_start_leaves_rank_count[i] =
        send_start_leaves_rank_shift[i+1] - send_start_leaves_rank_shift[i];
    }


    PDM_MPI_Alltoall (send_start_leaves_rank_count, 1, PDM_MPI_INT,
                      recv_start_leaves_rank_count, 1, PDM_MPI_INT,
                      octree->comm);

    int *recv_start_leaves_rank_shift = malloc (sizeof(int) * (lComm+1));
    recv_start_leaves_rank_shift[0] = 0;
    for (int i = 0; i < lComm; i++) {
      recv_start_leaves_rank_shift[i+1] =
        recv_start_leaves_rank_shift[i] + recv_start_leaves_rank_count[i];
    }

    start_leaves = realloc (start_leaves, sizeof(int) * recv_start_leaves_rank_shift[lComm]);
    PDM_MPI_Alltoallv (send_start_leaves,
                       send_start_leaves_rank_count,
                       send_start_leaves_rank_shift,
                       PDM_MPI_INT,
                       start_leaves,
                       recv_start_leaves_rank_count,
                       recv_start_leaves_rank_shift,
                       PDM_MPI_INT,
                       octree->comm);
    free (send_start_leaves);
    free (send_start_leaves_rank_shift);


    int *recv_start_leaves_count = malloc (sizeof(int) * recv_shift[lComm]);
    PDM_MPI_Alltoallv (send_start_leaves_count,
                       send_count,
                       send_shift,
                       PDM_MPI_INT,
                       recv_start_leaves_count,
                       recv_count,
                       recv_shift,
                       PDM_MPI_INT,
                       octree->comm);

    start_leaves_idx = realloc (start_leaves_idx, sizeof(int) * (recv_shift[lComm]+1));
    start_leaves_idx[0] = 0;
    for (int i = 0; i < recv_shift[lComm]; i++) {
      start_leaves_idx[i+1] = start_leaves_idx[i] + recv_start_leaves_count[i];
    }


    free (send_start_leaves_rank_count);
    free (recv_start_leaves_rank_count);
    free (recv_start_leaves_rank_shift);

    free (send_start_leaves_count);
    free (recv_start_leaves_count);


    // send coords
    send_coord = realloc (send_coord, sizeof(double) * send_shift[lComm] * dim);

    for (int i = 0; i < send_shift[lComm]; i++) {
      for (int j = 0; j < dim; j++)
        send_coord[dim*i + j] = recv_coord[dim*send_tgt_lnum[i] + j];
    }
    free (send_tgt_lnum);

    n_recv_pts = recv_shift[lComm];
    for (int i = 0; i < lComm; i++) {
      send_count[i]   *= dim;
      recv_count[i]   *= dim;
      send_shift[i+1] *= dim;
      recv_shift[i+1] *= dim;
    }
    recv_coord = realloc (recv_coord, sizeof(double) * recv_shift[lComm]);
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       octree->comm);


    local_closest_src_gnum = realloc (local_closest_src_gnum, sizeof(PDM_g_num_t) * n_recv_pts * n_closest_points);
    local_closest_src_dist = realloc (local_closest_src_dist, sizeof(double)      * n_recv_pts * n_closest_points);

  } // end while loop

  /* Free stuff */
  free (start_leaves);
  free (start_leaves_idx);
  free (local_closest_src_gnum);
  free (local_closest_src_dist);
  free (upper_bound_dist);

  free (send_coord);
  free (send_g_num);
  free (send_count);
  free (send_shift);

  free (recv_coord);
  free (recv_g_num);
  free (recv_count);
  free (recv_shift);

  _min_heap_free (merge_heap);
  free (tmp_closest_src_dist);
  free (tmp_closest_src_gnum);

  free (block_upper_bound_dist);


  for (int i = 0; i < octree->n_connected; i++) {
    free (processed_tgt[i]);
    free (new_processed_tgt[i]);
  }
  free (s_processed_tgt);
  free (n_processed_tgt);
  free (processed_tgt);
  free (s_new_processed_tgt);
  free (n_new_processed_tgt);
  free (new_processed_tgt);


  /* Final Block-to-part exchanges */
  PDM_block_to_part_t *btp1 = PDM_block_to_part_create (block_distrib_idx1,
                                                        (const PDM_g_num_t **) &pts_g_num,
                                                        &n_pts,
                                                        1,
                                                        octree->comm);
  PDM_part_to_block_free (ptb1);

  PDM_block_to_part_exch (btp1,
                          sizeof(double),
                          PDM_STRIDE_CST,
                          &_n_closest_points,
                          block_closest_src_dist1,
                          NULL,
                          (void **) &closest_octree_pt_dist2);

  PDM_block_to_part_exch (btp1,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          &_n_closest_points,
                          block_closest_src_gnum1,
                          NULL,
                          (void **) &closest_octree_pt_g_num);

  free (block_closest_src_dist1);
  free (block_closest_src_gnum1);

  btp1 = PDM_block_to_part_free (btp1);
}



















static void
_local_search1
(
 const _octree_t *octree,
 const int        n_closest_points,
 const int        n_pts,
 double          *pts_coord,
 PDM_g_num_t     *pts_g_num, // ONLY FOR DEBUG
 int             *start_leaf,
 PDM_g_num_t     *local_closest_src_gnum,
 double          *local_closest_src_dist
 )
{
  const int DEBUG = 0;
  const int CHECK_FACE_DIST = 0;

  const _l_octant_t *octants = octree->octants;
  const int dim = octree->dim;

  int myRank;
  PDM_MPI_Comm_rank (octree->comm, &myRank);

  int *is_visited     = malloc (sizeof(int) * octants->n_nodes);
  int *visited_leaves = malloc (sizeof(int) * octants->n_nodes);// could be smaller...
  int n_visited = 0;
  for (int i = 0; i < octants->n_nodes; i++) {
    is_visited[i] = 0;
  }

  /* Min heap used to visit leaves from neighbour to neighbour */
  _min_heap_t *leaf_heap = _min_heap_create (octants->n_nodes);

  /* Loop over target points */
  for (int i_tgt = 0; i_tgt < n_pts; i_tgt++) {
    /* Init */
    for (int i = 0; i < n_visited; i++) {
      is_visited[visited_leaves[i]] = 0;
    }
    n_visited = 0;

    PDM_g_num_t *closest_src_gnum = local_closest_src_gnum + n_closest_points * i_tgt;
    double      *closest_src_dist = local_closest_src_dist + n_closest_points * i_tgt;
    for (int j = 0; j < n_closest_points; j++) {
      closest_src_dist[j] = HUGE_VAL;
      closest_src_gnum[j] = -1; // USEFUL ONLY FOR DEBUG
    }
    double *max_src_dist = closest_src_dist + n_closest_points - 1;

    /* Get current target point data */
    const double *_pt = pts_coord + i_tgt * dim;

    /* Push start leaf in priority queue */
    int start_id = start_leaf[i_tgt];
    _min_heap_reset (leaf_heap);
    _min_heap_push (leaf_heap,
                    start_id,
                    0,
                    0.);
    is_visited[start_id] = 1;
    visited_leaves[n_visited++] = start_id;

    /* Visit octree leaves from neighbour to neighbour using priority queue */
    int leaf_id;
    double leaf_dist;
    PDM_g_num_t unused;
    while (_min_heap_pop (leaf_heap, &leaf_id, &unused, &leaf_dist)) {

      if (DEBUG) {
        printf("tgt point ("PDM_FMT_G_NUM") inspecting leaf %d: dist = %f / %f\n",
               pts_g_num[i_tgt], leaf_id, leaf_dist, *max_src_dist);
      }

      if (leaf_dist >= *max_src_dist) {
        break;
      }

      /* inspect source points inside popped leaf */
      for (int i = 0; i < octants->n_points[leaf_id]; i++) {

        // get source point coords and gnum
        int i_src = octants->range[leaf_id] + i;
        double *src_pt = octree->points + i_src * dim;
        PDM_g_num_t src_gnum = octree->points_gnum[i_src];

        // compute (squared) distance from target point
        double src_dist = 0;
        for (int j = 0; j < dim; j++) {
          double delta = _pt[j] - src_pt[j];
          src_dist += delta * delta;
        }

        if (DEBUG && pts_g_num[i_tgt] == 38) {
          printf("  src pt ("PDM_FMT_G_NUM") [%f %f %f] at dist %f / %f\n",
                 src_gnum, src_pt[0], src_pt[1], src_pt[2], src_dist, *max_src_dist);
        }

        // insertion sort
        if (src_dist < *max_src_dist) {
          int j = n_closest_points - 1;
          while (j > 0 && src_dist < closest_src_dist[j-1]) {
            closest_src_gnum[j] = closest_src_gnum[j-1];
            closest_src_dist[j] = closest_src_dist[j-1];
            j--;
          }

          closest_src_gnum[j] = src_gnum;
          closest_src_dist[j] = src_dist;

          if (DEBUG) {
            printf("  src pt ("PDM_FMT_G_NUM") [%f %f %f] at dist %f / %f --> insert at pos %d\n",
                   src_gnum, src_pt[0], src_pt[1], src_pt[2], src_dist, *max_src_dist, j);
          }
        }
      } // end loop over source points inside popped leaf


      /* inspect neighbours of popped leaf */
      double side = 1./(double) (1 << octants->codes[leaf_id].L);
      int check_dist;
      for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
        if (CHECK_FACE_DIST) {
          /* coarse test to discard neighbours that are clearly beyond the current search radius */
          int i_dim = dir / 2;
          int sign = dir % 2;

          double x_face = octree->s[i_dim] + octree->d[i_dim] * side * (octants->codes[leaf_id].X[i_dim] + sign);

          if (sign == 0) {
            check_dist = _pt[i_dim] > x_face;
          } else {
            check_dist = _pt[i_dim] < x_face;
          }

          if (check_dist) {
            double dist_face = (_pt[i_dim] -  x_face) * (_pt[i_dim] -  x_face);

            if (dist_face >= *max_src_dist) {
              /* all neighbours in direction dir are at least
                 at a distance from the tgt point equal to dist_face
                 so they need not be inspected */
              continue;
            }
          }
        }

        /* inspect neighbours in direction dir */
        for (int i = octants->neighbour_idx[6*leaf_id + dir];
             i < octants->neighbour_idx[6*leaf_id + dir + 1]; i++) {

          int ngb = octants->neighbours[i];

          if (ngb >= 0) { // ignore distant neighbours
            // local neighbour
            if (is_visited[ngb] == 1) continue;

            is_visited[ngb] = 1;
            visited_leaves[n_visited++] = ngb;

            // compute min dist from target point to current neighbor leaf
            double ngb_min_dist = _octant_min_dist2 (dim,
                                                     octants->codes[ngb],
                                                     octree->d,
                                                     octree->s,
                                                     _pt);

            if (ngb_min_dist < *max_src_dist) {
              // push current neighbour in priority queue
              _min_heap_push (leaf_heap,
                              ngb,
                              0,
                              ngb_min_dist);
            }
          } // end if local neighbour
        } // end loop over neighbours in direction dir (i)
      } // end loop over directions (dir)

    } // end while (_min_heap_pop (leaf_heap))
  } // end loop over target points (i_tgt)
}








static void
_local_search2
(
 const _octree_t *octree,
 const int        n_closest_points,
 const int        n_pts,
 double          *pts_coord,
 PDM_g_num_t     *pts_g_num, // ONLY FOR DEBUG
 double          *upper_bound_dist,
 PDM_g_num_t     *local_closest_src_gnum,
 double          *local_closest_src_dist
 )
{
  const int DEBUG = 0;

  const _l_octant_t *octants = octree->octants;
  const int dim = octree->dim;

  int myRank;
  PDM_MPI_Comm_rank (octree->comm, &myRank);

  /* Min heap used to visit octants in ascending order of distance */
  _min_heap_t *leaf_heap = _min_heap_create (octants->n_nodes);

  /* Encode pts_coords */
  PDM_morton_code_t *pts_code = malloc (sizeof(PDM_morton_code_t) * n_pts);
  double d[3], s[3];
  PDM_morton_encode_coords (dim,
                            PDM_morton_max_level,
                            octree->global_extents,
                            (size_t) n_pts,
                            pts_coord,
                            pts_code,
                            d,
                            s);


  double            box[2*dim];
  PDM_morton_code_t box_corners[2];
  double            extents[2*dim];
  /* Loop over target points */
  for (int i_tgt = 0; i_tgt < n_pts; i_tgt++) {

    /* Init */
    PDM_g_num_t *closest_src_gnum = local_closest_src_gnum + n_closest_points * i_tgt;
    double      *closest_src_dist = local_closest_src_dist + n_closest_points * i_tgt;
    for (int j = 0; j < n_closest_points; j++) {
      closest_src_dist[j] = HUGE_VAL;
      closest_src_gnum[j] = -1; // USEFUL ONLY FOR DEBUG
    }
    double *max_src_dist = closest_src_dist + n_closest_points - 1;

    /* Get current target point data */
    double   dist = sqrt(upper_bound_dist[i_tgt]);
    double *coord = pts_coord + dim * i_tgt;

    /* Encode corners of search box */
    for (int i_dim = 0; i_dim < dim; i_dim++) {
      box[i_dim]       = PDM_MAX (octree->global_extents[i_dim],
                                  coord[i_dim] - dist);
      box[dim + i_dim] = PDM_MAX (octree->global_extents[dim + i_dim],
                                  coord[i_dim] + dist);
    }

    PDM_morton_encode_coords (dim,
                              PDM_morton_max_level,
                              octree->global_extents,
                              (size_t) 2,
                              box,
                              box_corners,
                              d,
                              s);

    /* Find octants that intersect the search box */
    /*   1) Find narrow range of octants that might intersect the box */
    size_t l, r;
    _maximal_intersecting_range (box_corners[0],
                                 box_corners[1],
                                 octants->codes,
                                 octants->n_nodes,
                                 &l,
                                 &r);

    /*   2) Perform intersection test */
    _min_heap_reset (leaf_heap);
    for (int i = l; i < (int) r; i++) {
      int intersect = 1;
      PDM_morton_code_t *code = octants->codes + i;

      double side = 1. / (double) (1 << code->L);

      for (int j = 0; j < dim; j++) {
        extents[j]       = octree->s[j] + octree-> d[j] * side * code->X[j];
        extents[dim + j] = extents[j] + octree->d[j] * side;

        if (extents[dim + j] < box[j] || extents[j] > box[dim + j]) {
          intersect = 0;
          break;
        }
      }

      if (intersect) {
        /* Compute min distance between octant and target point */
        double min_dist = 0., delta;
        for (int j = 0; j < dim; j++) {
          double x = coord[j];

          if (x > extents[dim + j]) {
            delta = x - extents[dim + j];
            min_dist += delta * delta;
          } else if (x < extents[j]) {
            delta = x - extents[j];
            min_dist += delta * delta;
          }
        }

        if (min_dist < upper_bound_dist[i_tgt]) {
          _min_heap_push (leaf_heap,
                          i,
                          0,
                          min_dist);
        }
      } // end if intersect
    } // end loop over candidates for intersection (i)


    /* Visit octants in ascending order of distance */
    int leaf_id;
    double leaf_dist;
    PDM_g_num_t unused;
    while (_min_heap_pop (leaf_heap, &leaf_id, &unused, &leaf_dist)) {

      if (leaf_dist >= *max_src_dist) {
        break;
      }

      /* inspect source points inside popped leaf */
      for (int i = 0; i < octants->n_points[leaf_id]; i++) {

        // get source point coords and gnum
        int i_src = octants->range[leaf_id] + i;
        double *src_pt = octree->points + i_src * dim;
        PDM_g_num_t src_gnum = octree->points_gnum[i_src];

        // compute (squared) distance from target point
        double src_dist = 0;
        for (int j = 0; j < dim; j++) {
          double delta = coord[j] - src_pt[j];
          src_dist += delta * delta;
        }

        if (DEBUG && pts_g_num[i_tgt] == 38) {
          printf("  src pt ("PDM_FMT_G_NUM") [%f %f %f] at dist %f / %f\n",
                 src_gnum, src_pt[0], src_pt[1], src_pt[2], src_dist, *max_src_dist);
        }

        // insertion sort
        if (src_dist < *max_src_dist) {
          int j = n_closest_points - 1;
          while (j > 0 && src_dist < closest_src_dist[j-1]) {
            closest_src_gnum[j] = closest_src_gnum[j-1];
            closest_src_dist[j] = closest_src_dist[j-1];
            j--;
          }

          closest_src_gnum[j] = src_gnum;
          closest_src_dist[j] = src_dist;

          if (DEBUG) {
            printf("  src pt ("PDM_FMT_G_NUM") [%f %f %f] at dist %f / %f --> insert at pos %d\n",
                   src_gnum, src_pt[0], src_pt[1], src_pt[2], src_dist, *max_src_dist, j);
          }
        }
      } // end loop over source points inside popped leaf
    } // end while (_min_heap_pop (leaf_heap))


  } // end loop over tgt points (i_tgt)
}



/*
 * More robust/exhaustive (but much slower) implementaion
 */
void
PDM_para_octree_closest_point2
(
 const int    id,
 const int    n_closest_points,
 const int    n_pts,
 double      *pts,
 PDM_g_num_t *pts_g_num,
 PDM_g_num_t *closest_octree_pt_g_num,
 double      *closest_octree_pt_dist2
 )
{
  const int DEBUG = 0;
  const int DEBUG_MERGE = 0;

  const _octree_t *octree = _get_from_id (id);
  const _l_octant_t *octants = octree->octants;

  const int dim = octree->dim;
  int _n_closest_points = n_closest_points;
  int _n_pts = n_pts;

  int myRank, lComm;
  PDM_MPI_Comm_rank (octree->comm, &myRank);
  PDM_MPI_Comm_size (octree->comm, &lComm);


  /* Force target points inside octree extents /!\ /!\ /!\ -->> */
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < dim; j++) {
      pts[dim*i+j] = PDM_MAX (pts[dim*i+j], octree->global_extents[j]);
      pts[dim*i+j] = PDM_MIN (pts[dim*i+j], octree->global_extents[dim+j]);
    }
  }
  /* <<-- */


  /* Part-to-block create (only to get block distribution) */
  PDM_part_to_block_t *ptb1 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_NOTHING,
                                                        1.,
                                                        &pts_g_num,
                                                        NULL,
                                                        &_n_pts,
                                                        1,
                                                        octree->comm);

  PDM_g_num_t *block_distrib_idx1 = PDM_part_to_block_distrib_index_get (ptb1);
  const int n_pts_block1 = PDM_part_to_block_n_elt_block_get (ptb1);

  PDM_g_num_t *block_closest_src_gnum1 = malloc (sizeof(PDM_g_num_t) * n_pts_block1 * n_closest_points);
  double      *block_closest_src_dist1 = malloc (sizeof(double)      * n_pts_block1 * n_closest_points);

  for (int i = 0; i < n_pts_block1; i++) {
    for (int j = 0; j < n_closest_points; j++) {
      block_closest_src_dist1[n_closest_points*i + j] = HUGE_VAL;
      block_closest_src_gnum1[n_closest_points*i + j] = -1; // USEFUL ONLY FOR DEBUG
      /*block_closest_src_dist1[n_closest_points*i + j] = 1234567890;
        block_closest_src_gnum1[n_closest_points*i + j] = 1; // USEFUL ONLY FOR DEBUG*/
    }
  }



  /*************************************************************************
   *
   * Distribute the target points
   *
   *************************************************************************/
  /*   1) Encode the coordinates of every target point */
  PDM_morton_code_t *pts_code = malloc (sizeof(PDM_morton_code_t) * n_pts);
  double d[3], s[3];
  PDM_morton_encode_coords (dim,
                            PDM_morton_max_level,
                            octree->global_extents,
                            (size_t) n_pts,
                            pts,
                            pts_code,
                            d,
                            s);

  /*   2) Use binary search to associate each target point to the appropriate process */
  int *send_count = malloc (sizeof(int) * lComm);
  int *recv_count = malloc (sizeof(int) * lComm);
  for (int i = 0; i < lComm; i++) {
    send_count[i] = 0;
  }

  int *rank_pt = malloc (sizeof(int) * n_pts);
  for (int i = 0; i < n_pts; i++) {
    rank_pt[i] = PDM_morton_binary_search (lComm,
                                           pts_code[i],
                                           octree->rank_octants_index);
    send_count[rank_pt[i]]++;
  }
  free (pts_code);

  /*   3) Exchange send/recv counts */
  PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                    recv_count, 1, PDM_MPI_INT,
                    octree->comm);

  int *send_shift = malloc (sizeof(int) * (lComm+1));
  int *recv_shift = malloc (sizeof(int) * (lComm+1));
  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int i = 0; i < lComm; i++) {
    send_shift[i+1] = send_shift[i] + send_count[i];
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
    send_count[i] = 0;
  }

  /*   4) Fill send buffers */
  PDM_g_num_t *send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[lComm]);
  PDM_g_num_t *recv_g_num = malloc (sizeof(PDM_g_num_t) * recv_shift[lComm]);
  double      *send_coord = malloc (sizeof(double)      * send_shift[lComm]*dim);
  double      *recv_coord = malloc (sizeof(double)      * recv_shift[lComm]*dim);
  for (int i = 0; i < n_pts; i++) {
    int rank = rank_pt[i];
    int k = send_shift[rank] + send_count[rank];
    send_g_num[k] = pts_g_num[i];
    for (int j = 0; j < dim; j++)
      send_coord[dim*k+j] = pts[dim*i+j];

    send_count[rank]++;
  }
  free (rank_pt);

  /*   5) Send gnum buffer */
  PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                     recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                     octree->comm);

  /*   6) Send coord buffer */
  int n_recv_pts = recv_shift[lComm];
  for (int i = 0; i < lComm; i++) {
    send_count[i] *= dim;
    recv_count[i] *= dim;
    send_shift[i+1] *= dim;
    recv_shift[i+1] *= dim;
  }
  PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                     recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                     octree->comm);



  /*****************
   *
   * 1st phase
   *
   *****************/
  if (myRank == 0) {
    printf("  First phase\n");
  }
  /* Encode the coordinates of the received target points */
  pts_code = malloc (sizeof(PDM_morton_code_t) * n_recv_pts);
  PDM_morton_encode_coords (dim,
                            PDM_morton_max_level,
                            octree->global_extents,
                            (size_t) n_recv_pts,
                            recv_coord,
                            pts_code,
                            d,
                            s);

  /* Find start leaf for each target point */
  int *start_leaf = malloc (sizeof(int) * n_recv_pts);
  for (int i = 0; i < n_recv_pts; i++) {
    start_leaf[i] = PDM_morton_binary_search (octants->n_nodes,
                                              pts_code[i],
                                              octants->codes);
  }
  free (pts_code);

  /* First local search to get an upper bound for distance */
  PDM_g_num_t *local_closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_recv_pts * n_closest_points);
  double      *local_closest_src_dist = malloc (sizeof(double)      * n_recv_pts * n_closest_points);

  _local_search1 (octree,
                  n_closest_points,
                  n_recv_pts,
                  recv_coord,
                  recv_g_num, // ONLY FOR DEBUG
                  start_leaf,
                  local_closest_src_gnum,
                  local_closest_src_dist);
  free (start_leaf);

  /* /!\ if octree->n_connected > 1 */
  int send_to_self = (octree->n_connected > 1);
  // TO DO: check if the search box intersects ANY unvisited part of current rank,
  //        if so, send to self, else do not



  /* Part-to-block results of 1st phase */
  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                         PDM_PART_TO_BLOCK_POST_NOTHING,
                                                         1.,
                                                         &recv_g_num,
                                                         block_distrib_idx1,
                                                         &n_recv_pts,
                                                         1,
                                                         octree->comm);

  int *block_stride2 = NULL;
  PDM_part_to_block_exch (ptb2,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          _n_closest_points,
                          NULL,
                          (void **) &local_closest_src_gnum,
                          &block_stride2,
                          (void **) &block_closest_src_gnum1);
  free (block_stride2);

  PDM_part_to_block_exch (ptb2,
                          sizeof(double),
                          PDM_STRIDE_CST,
                          _n_closest_points,
                          NULL,
                          (void **) &local_closest_src_dist,
                          &block_stride2,
                          (void **) &block_closest_src_dist1);
  free (block_stride2);

  ptb2 = PDM_part_to_block_free (ptb2);



  /*****************
   *
   * 2nd phase
   *
   *****************/
  if (myRank == 0) {
    printf("  Second phase\n");
  }
  int **send_lnum = malloc (sizeof(int *) * lComm);
  for (int i = 0; i < lComm; i++) {
    send_count[i] = 0;
    send_lnum[i] = malloc (sizeof(int) * n_recv_pts); // worst case, could be smaller and re-allocated when necessary...
  }

  double            box[2*dim];
  PDM_morton_code_t box_corners[2];
  for (int i_tgt = 0; i_tgt < n_recv_pts; i_tgt++) {
    double dist = sqrt(local_closest_src_dist[n_closest_points * (i_tgt+1) - 1]);
    double *coord = recv_coord + dim * i_tgt;

    /* Encode corners of search box */
    for (int i_dim = 0; i_dim < dim; i_dim++) {
      box[i_dim]       = PDM_MAX (octree->global_extents[i_dim],
                                  coord[i_dim] - dist);
      box[dim + i_dim] = PDM_MAX (octree->global_extents[dim + i_dim],
                                  coord[i_dim] + dist);
    }

    PDM_morton_encode_coords (dim,
                              PDM_morton_max_level,
                              octree->global_extents,
                              (size_t) 2,
                              box,
                              box_corners,
                              d,
                              s);

    /* Find ranks that intersect the search box */
    size_t l, r, unused;
    PDM_morton_quantile_intersect (lComm,
                                   box_corners[0],
                                   octree->rank_octants_index,
                                   &l,
                                   &unused);

    PDM_morton_quantile_intersect (lComm,
                                   box_corners[1],
                                   octree->rank_octants_index,
                                   &unused,
                                   &r);

    for (int i_rank = l; i_rank < (int) r; i_rank++) {
      if (i_rank == myRank && !send_to_self) {
        continue;
      }

      send_lnum[i_rank][send_count[i_rank]++] = i_tgt;
    }

  } // end loop over tgt points (i_tgt)


  /* Exchange send/recv count  */
  PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                    recv_count, 1, PDM_MPI_INT,
                    octree->comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int i = 0; i < lComm; i++) {
    send_shift[i+1] = send_shift[i] + send_count[i];
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
  }

  /* Fill send buffers */
  double *send_upper_bound_dist = malloc (sizeof(double) * send_shift[lComm]);
  send_g_num = realloc (send_g_num, sizeof(PDM_g_num_t) * send_shift[lComm]);
  send_coord = realloc (send_coord, sizeof(double) * send_shift[lComm] * dim);

  for (int i_rank = 0; i_rank < lComm; i_rank++) {
    for (int i = 0; i < send_count[i_rank]; i++) {
      int i_tgt = send_lnum[i_rank][i];

      double   dist = local_closest_src_dist[n_closest_points * (i_tgt+1) - 1];
      double *coord = recv_coord + dim * i_tgt;

      int j = send_shift[i_rank] + i;

      send_g_num[j] = recv_g_num[i_tgt];

      for (int k = 0; k < dim; k++) {
        send_coord[dim*j+k] = coord[k];
      }

      send_upper_bound_dist[j] = dist;
    }
  } // end loop over ranks (i_rank)
  free (send_lnum);

  double *recv_upper_bound_dist = malloc (sizeof(double) * recv_shift[lComm]);
  recv_g_num = realloc (recv_g_num, sizeof(PDM_g_num_t) * recv_shift[lComm]);
  recv_coord = realloc (recv_coord, sizeof(double) * recv_shift[lComm] * dim);

  /* Send gnum buffer */
  PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                     recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                     octree->comm);
  free (send_g_num);

  /* Send upper_bound_dist buffer */
  PDM_MPI_Alltoallv (send_upper_bound_dist, send_count, send_shift, PDM_MPI_DOUBLE,
                     recv_upper_bound_dist, recv_count, recv_shift, PDM_MPI_DOUBLE,
                     octree->comm);
  free (send_upper_bound_dist);

  /* Send coord buffer */
  n_recv_pts = recv_shift[lComm];
  for (int i = 0; i < lComm; i++) {
    send_count[i] *= dim;
    recv_count[i] *= dim;
    send_shift[i+1] *= dim;
    recv_shift[i+1] *= dim;
  }
  PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                     recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                     octree->comm);
  free (send_coord);
  free (send_count);
  free (send_shift);



  /* Second local search */
  local_closest_src_gnum = realloc (local_closest_src_gnum,
                                    sizeof(PDM_g_num_t) * n_recv_pts * n_closest_points);
  local_closest_src_dist = realloc (local_closest_src_dist,
                                    sizeof(double)      * n_recv_pts * n_closest_points);
  _local_search2 (octree,
                  n_closest_points,
                  n_recv_pts,
                  recv_coord,
                  recv_g_num, // ONLY FOR DEBUG
                  recv_upper_bound_dist,
                  local_closest_src_gnum,
                  local_closest_src_dist);

  free (recv_coord);
  free (recv_upper_bound_dist);
  free (recv_count);
  free (recv_shift);






  /* Part-to-block exchanges */
  int *block_stride3 = NULL;
  double      *block_closest_src_dist3 = NULL;
  PDM_g_num_t *block_closest_src_gnum3 = NULL;

  PDM_part_to_block_t *ptb3 = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                         PDM_PART_TO_BLOCK_POST_MERGE,
                                                         1.,
                                                         &recv_g_num,
                                                         block_distrib_idx1,
                                                         &n_recv_pts,
                                                         1,
                                                         octree->comm);

  int n_pts_block3 = PDM_part_to_block_n_elt_block_get (ptb3);
  PDM_g_num_t *block_tgt_gnum3 = PDM_part_to_block_block_gnum_get (ptb3);

  int *stride3 = malloc (sizeof(int) * n_recv_pts);
  for (int i = 0; i < n_recv_pts; i++) {
    stride3[i] = n_closest_points;
  }

  PDM_part_to_block_exch (ptb3,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          1,
                          &stride3,
                          (void **) &local_closest_src_dist,
                          &block_stride3,
                          (void **) &block_closest_src_dist3);
  free (block_stride3);

  PDM_part_to_block_exch (ptb3,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &stride3,
                          (void **) &local_closest_src_gnum,
                          &block_stride3,
                          (void **) &block_closest_src_gnum3);
  free (stride3);


  /* Merge block data */
  _min_heap_t *merge_heap = _min_heap_create (1024); // used to merge part-to-block results

  double      *tmp_closest_src_dist = malloc (sizeof(double)      * n_closest_points);
  PDM_g_num_t *tmp_closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_closest_points);

  int *block_idx3 = malloc (sizeof(int) * (n_pts_block3 + 1));
  block_idx3[0] = 0;
  for (int i = 0; i < n_pts_block3; i++) {
    block_idx3[i+1] = block_idx3[i] + block_stride3[i];
  }

  for (int i = 0; i < n_pts_block3; i++) {
    int id_block1 = (int) n_closest_points * (block_tgt_gnum3[i]-1 - block_distrib_idx1[myRank]);

    if (DEBUG && DEBUG_MERGE) {
      printf("\npoint ("PDM_FMT_G_NUM")\n", block_tgt_gnum3[i]);
    }

    int n_procs_pt = block_stride3[i] / n_closest_points;

    int *idx_proc = malloc (sizeof(int) * (n_procs_pt + 1));

    _min_heap_reset (merge_heap);

    for (int j = 0; j < n_closest_points; j++) {
      tmp_closest_src_gnum[j] = block_closest_src_gnum1[id_block1 + j];
      tmp_closest_src_dist[j] = block_closest_src_dist1[id_block1 + j];
    }

    if (DEBUG && DEBUG_MERGE) {
      printf("\tj = %d: ("PDM_FMT_G_NUM", %f)\n",
             -1, tmp_closest_src_gnum[0], tmp_closest_src_dist[0]);
    }

    _min_heap_push (merge_heap,
                    -1,
                    tmp_closest_src_gnum[0],
                    tmp_closest_src_dist[0]);
    idx_proc[0] = 0;


    for (int j = 0; j < n_procs_pt; j++) {
      int k = block_idx3[i] + j*n_closest_points;

      if (DEBUG && DEBUG_MERGE) {
        printf("\tj =  %d: ("PDM_FMT_G_NUM", %f)\n",
               j, block_closest_src_gnum3[k], block_closest_src_dist3[k]);
      }

      _min_heap_push (merge_heap,
                      j,
                      block_closest_src_gnum3[k],
                      block_closest_src_dist3[k]);
      idx_proc[j+1] = 0;
    }


    int _proc;
    PDM_g_num_t _gnum;
    double _dist;
    for (int j = 0; j < n_closest_points; j++) {
      int popped = _min_heap_pop (merge_heap,
                                  &_proc,
                                  &_gnum,
                                  &_dist);
      assert (popped);
      if (DEBUG && DEBUG_MERGE) {
        printf("\tpopped: "PDM_FMT_G_NUM" %f \t / %f\n", _gnum, _dist, block_closest_src_dist1[id_block1 + j]);
      }

      if (_dist >= block_closest_src_dist1[id_block1 + n_closest_points - 1]) break;

      block_closest_src_dist1[id_block1 + j] = _dist;
      block_closest_src_gnum1[id_block1 + j] = _gnum;

      if (DEBUG && DEBUG_MERGE) {
        printf("\t\t%d/%d: "PDM_FMT_G_NUM", %f\n",
               j+1, n_closest_points,
               block_closest_src_gnum1[id_block1 + j], block_closest_src_dist1[id_block1 + j]);
      }

      if (j >= n_closest_points - 1)
        break;

      idx_proc[_proc+1]++;

      if (_proc < 0) {
        _min_heap_push (merge_heap,
                        -1,
                        tmp_closest_src_gnum[idx_proc[0]],
                        tmp_closest_src_dist[idx_proc[0]]);
      } else {
        int k = block_idx3[i] + _proc*n_closest_points + idx_proc[_proc+1];
        _min_heap_push (merge_heap,
                        _proc,
                        block_closest_src_gnum3[k],
                        block_closest_src_dist3[k]);
      }

    }
    free (idx_proc);
  }
  if (DEBUG && DEBUG_MERGE) {
    printf("- - - - - - - - -\n\n");
  }

  free (block_idx3);
  free (block_stride3);
  free (block_closest_src_gnum3);
  free (block_closest_src_dist3);
  ptb3 = PDM_part_to_block_free (ptb3);
  free (recv_g_num);


  /* Final Block-to-part exchanges */
  PDM_block_to_part_t *btp1 = PDM_block_to_part_create (block_distrib_idx1,
                                                        (const PDM_g_num_t **) &pts_g_num,
                                                        &n_pts,
                                                        1,
                                                        octree->comm);
  PDM_part_to_block_free (ptb1);

  PDM_block_to_part_exch (btp1,
                          sizeof(double),
                          PDM_STRIDE_CST,
                          &_n_closest_points,
                          block_closest_src_dist1,
                          NULL,
                          (void **) &closest_octree_pt_dist2);

  PDM_block_to_part_exch (btp1,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          &_n_closest_points,
                          block_closest_src_gnum1,
                          NULL,
                          (void **) &closest_octree_pt_g_num);

  free (block_closest_src_dist1);
  free (block_closest_src_gnum1);

  btp1 = PDM_block_to_part_free (btp1);


  free (local_closest_src_gnum);
  free (local_closest_src_dist);

}




static void
_my_dump_times
(
 const _octree_t *octree
 )
{
  double tot = octree->times_elapsed[END] - octree->times_elapsed[BEGIN];

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  PDM_printf("[%d] Elapsed = %12.5es\n", rank, tot);

  PDM_printf("[%d]   total              = %12.5es (%d%%)\n",
	     rank,
	     octree->times_elapsed[BUILD_TOTAL],
	     (int) (100 * octree->times_elapsed[BUILD_TOTAL] / tot));

  PDM_printf("[%d]   order points       = %12.5es (%d%%)\n",
	     rank,
	     octree->times_elapsed[BUILD_ORDER_POINTS],
	     (int) (100 * octree->times_elapsed[BUILD_ORDER_POINTS] / tot));

  PDM_printf("[%d]   block partition    = %12.5es (%d%%)\n",
	     rank,
	     octree->times_elapsed[BUILD_BLOCK_PARTITION],
	     (int) (100 * octree->times_elapsed[BUILD_BLOCK_PARTITION] / tot));

  if (1 == 0) {
    PDM_printf("[%d]   local nodes        = %12.5es (%d%%)\n",
               rank,
               octree->times_elapsed[BUILD_LOCAL_NODES],
               (int) (100 * octree->times_elapsed[BUILD_LOCAL_NODES] / tot));

    PDM_printf("[%d]   local neighbours   = %12.5es (%d%%)\n",
               rank,
               octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS],
	     (int) (100 * octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS] / tot));
  }
  double t = octree->times_elapsed[BUILD_LOCAL_NODES] + octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS];
  PDM_printf("[%d]   local nodes + ngbs = %12.5es (%d%%)\n",
	     rank,
	     t,
	     (int) (100 * t / tot));
  //#endif

  PDM_printf("[%d]   distant neighbours = %12.5es (%d%%)\n",
    rank,
    octree->times_elapsed[BUILD_DISTANT_NEIGHBOURS],
    (int) (100 * octree->times_elapsed[BUILD_DISTANT_NEIGHBOURS] / tot));
}

/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_para_octree_dump_times
(
  const int id
)
{
  _octree_t *octree = _get_from_id (id);
  if (1 == 0) {
    _my_dump_times (octree);
    return;
  }

  double t1 = octree->times_elapsed[END] - octree->times_elapsed[BEGIN];
  double t2 = octree->times_cpu[END] - octree->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  double t_elaps_max[NTIMER];
  PDM_MPI_Allreduce (octree->times_elapsed, t_elaps_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  double t_cpu_max[NTIMER];
  PDM_MPI_Allreduce (octree->times_cpu, t_cpu_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  if (rank == 0) {

    PDM_printf( "PDM_para_octree timer : all (elapsed and cpu)                                           :"
		" %12.5es %12.5es\n",
                t1max, t2max);
    PDM_printf( "PDM_para_octree timer : build octree : total (elapsed and cpu)                          :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_TOTAL],
                t_cpu_max[BUILD_TOTAL]);
    PDM_printf( "PDM_para_octree timer : build octree : step order points (elapsed and cpu)              :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_ORDER_POINTS],
                t_cpu_max[BUILD_ORDER_POINTS]);
    PDM_printf( "PDM_para_octree timer : build octree : step block partition (elapsed and cpu)           :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_BLOCK_PARTITION],
                t_cpu_max[BUILD_BLOCK_PARTITION]);
    PDM_printf( "PDM_para_octree timer : build octree : step local nodes (elapsed and cpu)               :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NODES],
                t_cpu_max[BUILD_LOCAL_NODES]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours (elapsed and cpu)          :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours - step 1 (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS_STEP1],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS_STEP1]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours - step 2 (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS_STEP2],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS_STEP2]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours - step 3 (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS_STEP3],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS_STEP3]);
    PDM_printf( "PDM_para_octree timer : build octree : step distant neighbours (elapsed and cpu)        :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_DISTANT_NEIGHBOURS],
                t_cpu_max[BUILD_DISTANT_NEIGHBOURS]);

}

}





//---->>>
static
void write_octree_octants
(
 const int   id,
 const char *filename
 )
{
  const int HEXA = 1;

  _octree_t *octree = _get_from_id (id);
  _l_octant_t *octants = octree->octants;

  int myRank;
  PDM_MPI_Comm_rank (octree->comm, &myRank);

  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "octree_octants_%4.4d\n", myRank);
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*octants->n_nodes);

  double ext[6];
  int ii = 0;
  for (int inode = 0; inode < octants->n_nodes; inode++) {
    double side = 1./ (double) (1 << octants->codes[inode].L);
    for (int idim = 0; idim < 3; idim++) {
      double co = octants->codes[inode].X[idim];
      ext[idim] = octree->s[idim] + octree->d[idim]*co*side;
      ext[3+idim] = ext[idim] + octree->d[idim]*side;
    }
    for (int k = 0; k < 2; k++) {
      for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
          if (HEXA == 1) {
            ii = (1-j)*i + j*(1-i);
          } else {
            ii = i;
	  }
          fprintf(f, "%f %f %f\n", ext[3*ii], ext[3*j+1], ext[3*k+2]);
        }
      }
    }
  }

  fprintf(f, "CELLS %d %d\n", octants->n_nodes, 9*octants->n_nodes);
  for (int i = 0; i < octants->n_nodes; i++) {
    fprintf(f, "8 ");
    for (int j = 0; j < 8; j++) {
      fprintf(f, "%d ", 8*i+j);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELL_TYPES %d\n", octants->n_nodes);
  for (int i = 0; i < octants->n_nodes; i++) {
    fprintf(f, "%d\n", 11+HEXA);
  }

  fprintf(f, "CELL_DATA %d\n", octants->n_nodes);
  if (octree->neighboursToBuild) {
    fprintf(f, "SCALARS connected_comp int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < octree->n_connected; i++) {
      for (int j = octree->connected_idx[i]; j < octree->connected_idx[i+1]; j++) {
        fprintf(f, "%d\n", i);
      }
    }
  } else {
    fprintf(f, "SCALARS order int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < octants->n_nodes; i++) {
      fprintf(f, "%d\n", i);
    }
  }

  fclose(f);
}




static
void write_octree_points
(
 const int   id,
 const char *filename,
 const int   use_gnum
 )
{
  _octree_t *octree = _get_from_id (id);

  int myRank;
  PDM_MPI_Comm_rank (octree->comm, &myRank);

  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "octree_points_%4.4d\n", myRank);
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", octree->n_points);
  for (int i = 0; i < octree->n_points; i++) {
    for (int j = 0; j < octree->dim; j++) {
      fprintf(f, "%f ", octree->points[octree->dim*i+j]);
    }
    for (int j = 0; j < 3 - octree->dim; j++) {
      fprintf(f, "0 ");
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", octree->n_points, 2*octree->n_points);
  for (int i = 0; i < octree->n_points; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", octree->n_points);
  for (int i = 0; i < octree->n_points; i++) {
    fprintf(f, "1\n");
  }

  fprintf(f, "CELL_DATA %d\n", octree->n_points);
  if (use_gnum) {
    fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < octree->n_points; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", octree->points_gnum[i]);
    }
  } else {
    fprintf(f, "SCALARS lnum int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < octree->n_points; i++) {
      fprintf(f, "%d\n", i);
    }
  }

  fclose(f);
}


static
void write_boxes
(
 const PDM_MPI_Comm comm,
 const int          n_boxes,
 const double       boxes_extents[],
 const PDM_g_num_t  boxes_g_num[],
 const char        *filename
 )
{
  const int HEXA = 1;

  int myRank;
  PDM_MPI_Comm_rank (comm, &myRank);

  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "boxes_%4.4d\n", myRank);
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_boxes);

  const double *ext;
  int ii = 0;
  for (int ibox = 0; ibox < n_boxes; ibox++) {
    ext = boxes_extents + 6*ibox;

    for (int k = 0; k < 2; k++) {
      for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
          if (HEXA) {
            ii = (1-j)*i + j*(1-i);
          } else {
            ii = i;
	  }
          fprintf(f, "%f %f %f\n", ext[3*ii], ext[3*j+1], ext[3*k+2]);
        }
      }
    }
    /*fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+0], boxes_extents[6*ibox+1], boxes_extents[6*ibox+2]);//0
    fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+3], boxes_extents[6*ibox+1], boxes_extents[6*ibox+2]);//1
    if (HEXA) {
      fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+3], boxes_extents[6*ibox+4], boxes_extents[6*ibox+2]);//2
      fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+0], boxes_extents[6*ibox+4], boxes_extents[6*ibox+2]);//3
    } else {
      fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+0], boxes_extents[6*ibox+4], boxes_extents[6*ibox+2]);//2
      fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+3], boxes_extents[6*ibox+4], boxes_extents[6*ibox+2]);//3
    }

    fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+0], boxes_extents[6*ibox+1], boxes_extents[6*ibox+5]);//4
    fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+3], boxes_extents[6*ibox+1], boxes_extents[6*ibox+5]);//5
    if (HEXA) {
      fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+3], boxes_extents[6*ibox+4], boxes_extents[6*ibox+5]);//6
      fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+0], boxes_extents[6*ibox+4], boxes_extents[6*ibox+5]);//7
    } else {
      fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+0], boxes_extents[6*ibox+4], boxes_extents[6*ibox+5]);//6
      fprintf(f, "%f %f %f\n", boxes_extents[6*ibox+3], boxes_extents[6*ibox+4], boxes_extents[6*ibox+5]);//7
      }*/
  }

  fprintf(f, "CELLS %d %d\n", n_boxes, 9*n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    fprintf(f, "8 ");
    for (int j = 0; j < 8; j++) {
      fprintf(f, "%d ", 8*i+j);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELL_TYPES %d\n", n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    fprintf(f, "%d\n", 11+HEXA);
  }

  fprintf(f, "CELL_DATA %d\n", n_boxes);
  fprintf(f, "SCALARS box_gnum int\n LOOKUP_TABLE default\n");
  for (int i = 0; i < n_boxes; i++) {
    fprintf(f, ""PDM_FMT_G_NUM"\n", boxes_g_num[i]);
  }

  fclose(f);
}

/**
 *
 * Get the location of a point cloud
 *
 */

void
PDM_para_octree_points_inside_boxes
(
 const int           octree_id,
 const int           n_boxes,
 const double       *box_extents,
 const PDM_g_num_t  *box_g_num,
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord
 )
{
  const int DEBUG = 0;

  const int VISU = 0;
  const int VISU_POINTS_GNUM = 1;

  int n_recv_boxes = 0;
  PDM_morton_code_t *box_corners = NULL;
  double *recv_box_extents       = NULL;
  const PDM_g_num_t *recv_box_g_num;
  PDM_g_num_t *_recv_box_g_num   = NULL;
  PDM_part_to_block_t *ptb1      = NULL;
  PDM_g_num_t *block_distrib_idx = NULL;
  double s[3], d[3];


  _octree_t *octree = _get_from_id (octree_id);
  const int dim = octree->dim;
  const int two_dim = 2 * dim;

  const _l_octant_t *octants = octree->octants;

  int my_rank;
  PDM_MPI_Comm_rank (octree->comm, &my_rank);

  int n_ranks;
  PDM_MPI_Comm_size (octree->comm, &n_ranks);

  if (VISU) {
    char filename[999];

    sprintf(filename, "octants_%3.3d.vtk", my_rank);
    write_octree_octants (octree_id,
                          filename);

    sprintf(filename, "points_%3.3d.vtk", my_rank);
    write_octree_points (octree_id,
                         filename,
                         VISU_POINTS_GNUM);

    sprintf(filename, "boxes_%3.3d.vtk", my_rank);
    write_boxes (octree->comm,
                 n_boxes,
                 box_extents,
                 box_g_num,
                 filename);
  }

  /* Clip box extents (ensure Morton codes of box corners are properly computed) */
  double *_box_extents = malloc (sizeof(double) * two_dim * n_boxes);
  for (int ibox = 0; ibox < n_boxes; ibox++) {
    for (int idim = 0; idim < dim; idim++) {
      _box_extents[two_dim*ibox + idim]       = PDM_MAX (octree->global_extents[idim],
                                                         box_extents[two_dim*ibox + idim]);
      _box_extents[two_dim*ibox + dim + idim] = PDM_MIN (octree->global_extents[idim + dim],
                                                         box_extents[two_dim*ibox + dim + idim]);
    }
  }


  /* Multiple ranks */
  if (n_ranks > 1) {
    /* Part-to-block create (only to get block distribution) */
    ptb1 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                     PDM_PART_TO_BLOCK_POST_MERGE,
                                     1.,
                                     (PDM_g_num_t **) (&box_g_num),
                                     NULL,
                                     (int *) &n_boxes,
                                     1,
                                     octree->comm);

    block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb1);



    /***************************************
     * Redistribute bounding boxes
     ***************************************/
    int *send_count = malloc (sizeof(int) * n_ranks);
    for (int i = 0; i < n_ranks; i++) {
      send_count[i] = 0;
    }

    /* Encode box corners */
    box_corners = malloc (sizeof(PDM_morton_code_t) * 2 * n_boxes);
    PDM_morton_encode_coords (dim,
                              PDM_morton_max_level,
                              octree->global_extents,
                              2 * n_boxes,
                              _box_extents,
                              box_corners,
                              d,
                              s);

    size_t *box_rank = malloc (sizeof(int *) * 2 * n_boxes);

    /* Find which ranks possibly intersect each box */
    size_t start, end, tmp;
    for (int ibox = 0; ibox < n_boxes; ibox++) {
      PDM_morton_quantile_intersect (n_ranks,
                                     box_corners[2*ibox],
                                     octree->rank_octants_index,
                                     &start,
                                     &tmp);

      PDM_morton_quantile_intersect (n_ranks - start,
                                     box_corners[2*ibox+1],
                                     octree->rank_octants_index + start,
                                     &tmp,
                                     &end);
      end += start;

      box_rank[2*ibox]   = start;
      box_rank[2*ibox+1] = end;

      for (size_t irank = start; irank < end; irank++) {
        send_count[irank]++;
      }
    }
    free (box_corners);

    int *recv_count = malloc (sizeof(int) * n_ranks);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      octree->comm);

    int *send_shift = malloc (sizeof(int) * (n_ranks+1));
    int *recv_shift = malloc (sizeof(int) * (n_ranks+1));
    send_shift[0] = 0;
    recv_shift[0] = 0;
    for (int i = 0; i < n_ranks; i++) {
      send_shift[i+1] = send_shift[i] + send_count[i];
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
      send_count[i] = 0;
    }
    n_recv_boxes = recv_shift[n_ranks];

    /* Fill send buffers */
    PDM_g_num_t *send_box_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_ranks]);
    _recv_box_g_num = malloc (sizeof(PDM_g_num_t) * recv_shift[n_ranks]);
    recv_box_g_num = _recv_box_g_num;
    double *send_box_extents = malloc (sizeof(double) * send_shift[n_ranks] * two_dim);
    recv_box_extents = malloc (sizeof(double) * recv_shift[n_ranks] * two_dim);

    for (int ibox = 0; ibox < n_boxes; ibox++) {
      for (size_t irank = box_rank[2*ibox]; irank < box_rank[2*ibox+1]; irank++) {
        int idx = send_shift[irank] + send_count[irank];
        send_box_g_num[idx] = box_g_num[ibox];

        for (int k = 0; k < two_dim; k++) {
          send_box_extents[two_dim*idx + k] = _box_extents[two_dim*ibox + k];
        }

        send_count[irank]++;
      }
    }
    free (_box_extents);
    free (box_rank);

    /* Send boxes g_num buffer */
    PDM_MPI_Alltoallv (send_box_g_num,  send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       _recv_box_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       octree->comm);

    /* Send boxes extents buffer */
    for (int i = 0; i < n_ranks; i++) {
      send_shift[i+1] *= two_dim;
      recv_shift[i+1] *= two_dim;
      send_count[i]   *= two_dim;
      recv_count[i]   *= two_dim;
    }
    PDM_MPI_Alltoallv (send_box_extents, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_box_extents, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       octree->comm);

    free (send_count);
    free (recv_count);
    free (send_shift);
    free (recv_shift);
    free (send_box_g_num);
    free (send_box_extents);

  }

  /* Single rank */
  else {
    n_recv_boxes     = n_boxes;
    recv_box_extents = _box_extents;
    recv_box_g_num   = box_g_num;
  }



  if (VISU) {
    char filename[999];

    sprintf(filename, "recv_boxes_%3.3d.vtk", my_rank);
    write_boxes (octree->comm,
                 n_recv_boxes,
                 recv_box_extents,
                 recv_box_g_num,
                 filename);
  }

  /***************************************
   * Intersect redistributed boxes with local octree
   ***************************************/

  /* Encode corners of redistributed boxes */
  box_corners = malloc (sizeof(PDM_morton_code_t) * 2 * n_recv_boxes);
  PDM_morton_encode_coords (dim,
                            PDM_morton_max_level,
                            octree->global_extents,
                            2 * n_recv_boxes,
                            recv_box_extents,
                            box_corners,
                            d,
                            s);

  /* Root node of octree */
  PDM_morton_code_t root;
  PDM_morton_nearest_common_ancestor (octants->codes[0],
                                      octants->codes[octants->n_nodes - 1],
                                      &root);

  int *intersect_nodes = malloc (sizeof(int) * octants->n_nodes);
  size_t n_intersect_nodes;

  size_t s_box_pts = octree->n_points;
  int *box_pts = malloc (sizeof(int) * s_box_pts);
  int *box_pts_n = malloc (sizeof(int) * n_recv_boxes);
  int *box_pts_idx = malloc (sizeof(int) * (n_recv_boxes+1));
  box_pts_idx[0] = 0;

  /* Loop over redistributed boxes */
  for (int ibox = 0; ibox < n_recv_boxes; ibox++) {
    n_intersect_nodes = 0;
    box_pts_n[ibox] = 0;
    box_pts_idx[ibox+1] = box_pts_idx[ibox];

    /* Get list of all nodes (octants) that intersect the box */
    PDM_morton_intersect_box (dim,
                              root,
                              box_corners[2*ibox],
                              box_corners[2*ibox+1],
                              octants->codes,
                              0,
                              octants->n_nodes,
                              &n_intersect_nodes,
                              intersect_nodes);

    if (DEBUG) {
      printf("[%d]\tbox %d ("PDM_FMT_G_NUM") nodes:", my_rank, ibox, recv_box_g_num[ibox]);
      for (int j = 0; j < (int) n_intersect_nodes; j++) {
        printf(" %d", intersect_nodes[j]);
      }
      printf("\n");
    }

    size_t new_max_size = box_pts_idx[ibox] + n_intersect_nodes * octree->points_in_leaf_max;
    if (s_box_pts <= new_max_size) {
      s_box_pts = PDM_MAX (2*s_box_pts, new_max_size);
      box_pts = realloc (box_pts, sizeof(int) * s_box_pts);
    }

    /* Inspect nodes (octants) which intersect the box */
    double *box_min = recv_box_extents + two_dim*ibox;
    double *box_max = box_min + dim;

    for (size_t i = 0; i < n_intersect_nodes; i++) {
      int inode = intersect_nodes[i];

      /* Inspect points inside current node */
      for (int j = 0; j < octants->n_points[inode]; j++) {
        int ipt = octants->range[inode] + j;
        double *_pt = octree->points + ipt * dim;

        /* check whether current point lies inside the box */
        int inside = 1;
        for (int idim = 0; idim < dim; idim++) {
          if (_pt[idim] < box_min[idim] || _pt[idim] > box_max[idim]) {
            inside = 0;
            break;
          }
        }

        if (inside) {
          box_pts[box_pts_idx[ibox+1]] = ipt;
          box_pts_idx[ibox+1]++;
          box_pts_n[ibox]++;
        }
      }
    }

    if (DEBUG) {
      printf("[%d]\tbox %d ("PDM_FMT_G_NUM") pts:", my_rank, ibox, recv_box_g_num[ibox]);
      for (int j = box_pts_idx[ibox]; j < box_pts_idx[ibox+1]; j++) {
        printf(" %d", box_pts[j]);
      }
      printf("\n");
    }

  } // Loop over redistributed boxes
  free (recv_box_extents);
  free (box_corners);
  free (intersect_nodes);


  /* Get gnum and coords of points inside boxes */
  PDM_g_num_t *box_pts_g_num = malloc (sizeof(PDM_g_num_t) * box_pts_idx[n_recv_boxes]);
  double      *box_pts_coord = malloc (sizeof(double)      * box_pts_idx[n_recv_boxes] * dim);

  for (int i = 0; i < box_pts_idx[n_recv_boxes]; i++) {
    int ipt = box_pts[i];

    box_pts_g_num[i] = octree->points_gnum[ipt];

    for (int j = 0; j < dim; j++) {
      box_pts_coord[dim*i + j] = octree->points[dim*ipt + j];
    }
  }
  free (box_pts);

  /* Multiple ranks */
  if (n_ranks > 1) {
    free (box_pts_idx);

    /* Part#2 to Block */
    PDM_part_to_block_t *ptb2 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                          (PDM_g_num_t **) (&recv_box_g_num),
                                                          NULL,
                                                          &n_recv_boxes,
                                                          1,
                                                          octree->comm);
    free (_recv_box_g_num);

    int *block_pts_in_box_n = NULL;
    PDM_g_num_t *block_pts_in_box_g_num = NULL;
    PDM_part_to_block_exch (ptb2,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            1,
                            &box_pts_n,
                            (void **) &box_pts_g_num,
                            &block_pts_in_box_n,
                            (void **) &block_pts_in_box_g_num);
    free (box_pts_g_num);

    int *box_pts_n_coord = malloc (sizeof(int) * n_recv_boxes); // overwrite box_pts_n?
    for (int ibox = 0; ibox < n_recv_boxes; ibox++) {
      box_pts_n_coord[ibox] = box_pts_n[ibox] * dim;
    }

    int *block_pts_in_box_n_coord = NULL;
    double *block_pts_in_box_coord = NULL;
    PDM_part_to_block_exch (ptb2,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &box_pts_n_coord,
                            (void **) &box_pts_coord,
                            &block_pts_in_box_n_coord,
                            (void **) &block_pts_in_box_coord);
    free (box_pts_coord);
    free (box_pts_n);
    free (box_pts_n_coord);

    /* Block to Part#1 */
    PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                         (const PDM_g_num_t **) &box_g_num,
                                                         &n_boxes,
                                                         1,
                                                         octree->comm);

    int *pts_in_box_n = malloc (sizeof(int) * n_boxes);
    int one = 1;
    PDM_block_to_part_exch (btp,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            &one,
                            (void *) block_pts_in_box_n,
                            NULL,
                            (void **) &pts_in_box_n);

    *pts_in_box_idx = malloc (sizeof(int) * (n_boxes + 1));
    int *_pts_in_box_idx = *pts_in_box_idx;
    _pts_in_box_idx[0] = 0;
    for (int ibox = 0; ibox < n_boxes; ibox++) {
      _pts_in_box_idx[ibox+1] = _pts_in_box_idx[ibox] + pts_in_box_n[ibox];
    }

    *pts_in_box_g_num = malloc (sizeof(PDM_g_num_t) * _pts_in_box_idx[n_boxes]);

    PDM_block_to_part_exch (btp,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            block_pts_in_box_n,
                            (void *) block_pts_in_box_g_num,
                            &pts_in_box_n,
                            (void **) pts_in_box_g_num);

    int *pts_in_box_n_coord = malloc (sizeof(int) * n_boxes); // overwrite pts_in_box_n?
    for (int ibox = 0; ibox < n_boxes; ibox++) {
      pts_in_box_n_coord[ibox] = pts_in_box_n[ibox] * dim;
    }

    *pts_in_box_coord = malloc (sizeof(double) * _pts_in_box_idx[n_boxes] * dim);

    PDM_block_to_part_exch (btp,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            block_pts_in_box_n_coord,
                            (void *) block_pts_in_box_coord,
                            &pts_in_box_n_coord,
                            (void **) pts_in_box_coord);

    PDM_part_to_block_free (ptb1);
    PDM_part_to_block_free (ptb2);
    PDM_block_to_part_free (btp);

    free (pts_in_box_n);
    free (pts_in_box_n_coord);
    free (block_pts_in_box_n);
    free (block_pts_in_box_n_coord);
    free (block_pts_in_box_g_num);
    free (block_pts_in_box_coord);

  }

  /* Single rank */
  else {
    free (box_pts_n);

    *pts_in_box_idx   = box_pts_idx;
    *pts_in_box_g_num = box_pts_g_num;
    *pts_in_box_coord = box_pts_coord;
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
