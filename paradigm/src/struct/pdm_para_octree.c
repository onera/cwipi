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
#include "pdm_box.h"
#include "pdm_box_tree.h"
#include "pdm_box_priv.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define NTIMER 12
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
  BUILD_EXPLICIT_NODES          = 9,
  BUILD_TOTAL                   = 10,
  END                           = 11,

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


typedef struct {
  PDM_morton_code_t code;
  int               n_points;
  int               range;
  int               ancestor_id;
  int               children_id[8];
  int               leaf_id; //-1 if internal, >=0 if leaf
  double            pts_extents[6];
} _explicit_node_t;


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



  PDM_box_set_t  *rank_boxes;            /*!< Rank Boxes */
  int             n_used_rank;           /*!< Number of used ranks */
  int            *used_rank;             /*!< used ranks */
  double         *used_rank_extents;     /*!< Extents of processes */
  PDM_box_tree_t *bt_shared;             /*!< Shared Boundary box tree */
  PDM_MPI_Comm    rank_comm;             /*!< MPI communicator */



  int                 n_copied_ranks;      /*!< Number of copies from other ranks */
  int                *copied_ranks;        /*!< Copied ranks */
  _l_octant_t       **copied_octants;      /*!< Octants from copied ranks */
  int                *n_copied_points;     /*!< Number of points copied from other ranks */
  double            **copied_points;       /*!< Coordinates of copied points */
  PDM_g_num_t       **copied_points_gnum;  /*!< Global numbers of copied points  */
  PDM_morton_code_t **copied_points_code;  /*!< Morton codes of copied points */


  /*
   *  Shared 'coarse' octree
   */
  int               *shared_rank_idx;
  PDM_morton_code_t *shared_codes;
  int               *shared_pts_n;
  double            *shared_pts_extents;


  int explicit_nodes_to_build;
  _explicit_node_t  *explicit_nodes;
  int                n_explicit_nodes;

  _explicit_node_t **copied_explicit_nodes;
  int               *n_copied_explicit_nodes;

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

  int                size;
  int                count;
  PDM_morton_code_t *code;
  int               *start;
  int               *end;
  double            *dist2;

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
static _min_heap_t *
_min_heap_create
(
 const int size
 )
{
  _min_heap_t *h = malloc (sizeof(_min_heap_t));
  h->size        = size;
  h->count       = 0;
  h->code        = malloc (sizeof(PDM_morton_code_t) * size);
  h->start       = malloc (sizeof(int) * size);
  h->end         = malloc (sizeof(int) * size);
  h->dist2       = malloc (sizeof(double) * size);

  for (int i = 0; i < h->size; i++) {
    h->dist2[i] = HUGE_VAL;
  }

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
  free (h->code);
  free (h->start);
  free (h->end);
  free (h->dist2);
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
  PDM_morton_code_t c;
  PDM_morton_copy (h->code[a], &c);
  int s = h->start[a];
  int e = h->end[a];
  double d = h->dist2[a];

  PDM_morton_copy (h->code[b], h->code + a);
  h->start[a] = h->start[b];
  h->end[a]   = h->end[b];
  h->dist2[a] = h->dist2[b];

  PDM_morton_copy (c, h->code + b);
  h->start[b] = s;
  h->end[b]   = e;
  h->dist2[b] = d;
}


static void
_min_heap_push
(
 _min_heap_t       *h,
 PDM_morton_code_t  code,
 const int          start,
 const int          end,
 const double       dist2
 )
{
  int i = 0;
  /* make sure the heap is large enough to contain the new element */
  if (h->count >= h->size) {
    h->size  = (h->size) ? 2*h->size : 10;
    h->code  = realloc (h->code,  sizeof(PDM_morton_code_t) * h->size);
    h->start = realloc (h->start, sizeof(int)               * h->size);
    h->end   = realloc (h->end,   sizeof(int)               * h->size);
    h->dist2 = realloc (h->dist2, sizeof(double)            * h->size);

    for (i = h->count+1; i < h->size; i++) {
      h->dist2[i] = HUGE_VAL;
    }
  }

  i = h->count;
  h->count++;

  PDM_morton_copy (code, h->code + i);
  h->start[i] = start;
  h->end[i]   = end;
  h->dist2[i] = dist2;

  int parent = (i - 1)/2;
  while (i > 0 && h->dist2[parent] > dist2) {
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

  if (l < h->count && h->dist2[l] < h->dist2[i])
    s = l;

  if (r < h->count && h->dist2[r] < h->dist2[s])
    s = r;

  if (s != i) {
    _min_heap_swap (h, i, s);
    _min_heap_heapify (h, s);
  }
}

static int
_min_heap_pop
(
 _min_heap_t       *h,
 PDM_morton_code_t *code,
 int               *start,
 int               *end,
 double            *dist2
 )
{
  if (h->count < 1) {
    return 0;
  }

  PDM_morton_copy (h->code[0], code);
  *start = h->start[0];
  *end   = h->end[0];
  *dist2 = h->dist2[0];

  h->count--;
  PDM_morton_copy (h->code[h->count], h->code);
  h->start[0] = h->start[h->count];
  h->end[0]   = h->end[h->count];
  h->dist2[0] = h->dist2[h->count];

  _min_heap_heapify (h, 0);

  return 1;
}

/**
 *
 * \brief Compute squared distance between a pair of points
 *
 * \param [in]   dim        Dimension
 * \param [in]   a          Coordinates of point A
 * \param [in]   b          Coordinates of point B
 *
 * \return squared euclidean distance between points A and B
 *
 */

inline static double
_pt_to_pt_dist2
(
 const int    dim,
 const double a[],
 const double b[]
 )
{
  double dist2 = 0.;
  for (int i = 0; i < dim; i++) {
    double delta = a[i] - b[i];
    dist2 += delta * delta;
  }

  return dist2;
}

/**
 *
 * \brief Compute squared distance from a point to an octant
 *
 * \param [in]   dim        Dimension
 * \param [in]   code       Morton code of octant
 * \param [in]   d          Normalization vector (dilatation)
 * \param [in]   s          Normalization vector (translation)
 * \param [in]   coords     Point coordinates
 *
 * \return square of minimal distance from the point to the octant
 *
 */

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


/**
 *
 * \brief Binary search in a sorted array of doubles
 *
 * \param [in]   elem    Query element
 * \param [in]   array   Sorted array
 * \param [in]   n       Array length
 *
 * \return position of first element in array greater than or equal to query element
 *
 */

static int
_binary_search_double
(
 const double  elem,
 const double *array,
 const int     n
 )
{
  int l = 0;
  int r = n;

  if (n < 1)
    return 0;

  while (l + 1 < r) {
    int m = l + (r - l)/2;

    if (elem < array[m])
      r = m;
    else
      l = m;
  }

  if (array[l] < elem)
    return l + 1;
  else
    return l;
}


/**
 *
 * \brief Insertion sort (in order of ascending distance)
 *
 * \param [in]     dist2         Distance of element to insert
 * \param [in]     g_num         Global number of element to insert
 * \param [in]     n             Array length
 * \param [inout]  array_dist2   Distances of elements in array
 * \param [inout]  array_g_num   Global numbers  of elements in array
 *
 */

static void
_insertion_sort
(
 const double       dist2,
 const PDM_g_num_t  g_num,
 const int          n,
 double            *array_dist2,
 PDM_g_num_t       *array_g_num
 )
{
  int i = _binary_search_double (dist2,
                                 array_dist2,
                                 n);
  if (array_g_num[i] == g_num) {
    return;
  }

  for (int j = n-1; j > i; j--) {
    array_dist2[j] = array_dist2[j-1];
    array_g_num[j] = array_g_num[j-1];
  }

  array_dist2[i] = dist2;
  array_g_num[i] = g_num;
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
  octants->n_nodes_max = PDM_MAX(init_size,1); /* To avoid mem leaks if init_size == 0 */
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
  if (octants->n_nodes + n_free_node >= octants->n_nodes_max) {

    //octants->n_nodes_max *= 2;
    octants->n_nodes_max = PDM_MAX (2*octants->n_nodes_max, octants->n_nodes + n_free_node);

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

  int *send_count = PDM_array_zeros_int(n_ranks);
  size_t *send_shift = malloc(sizeof(size_t) * (n_ranks+1));

  int *recv_count = malloc(sizeof(int) * n_ranks);
  size_t *recv_shift = malloc(sizeof(size_t) * (n_ranks+1));


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

  PDM_array_reset_int(send_count, n_ranks, 0);


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

    /* Avoid mem leaks L2*/
    _octants_free (L2);

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

  PDM_array_reset_int(send_count, n_ranks, 0);

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

  PDM_array_reset_int(send_count, n_ranks, 0);

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

  PDM_array_reset_int(send_count, n_ranks, 0);

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





static void
_compress_octants
(
 const _l_octant_t  *octants,
 const double       *points,
 int                *n_nodes,
 int               **nodes,
 int               **n_pts,
 double            **pts_extents
 )
{
  int DEBUG = 0;

  *n_nodes = 0;
  *nodes       = malloc (sizeof(int)    * octants->n_nodes * 4);
  *n_pts       = malloc (sizeof(int)    * octants->n_nodes);
  *pts_extents = malloc (sizeof(double) * octants->n_nodes * 6);
  int *_nodes = *nodes;
  int *_n_pts = *n_pts;
  double *_pts_extents = *pts_extents;

  if (octants->n_nodes <= 0) {
    return;
  }

  /* Deepest common ancestor of all octants */
  PDM_morton_code_t root;
  PDM_morton_nearest_common_ancestor (octants->codes[0],
                                      octants->codes[octants->n_nodes - 1],
                                      &root);

  PDM_morton_int_t _X[3];
  PDM_morton_int_t _L;


  const int n_child = 8;
  const int depth_max = 31;
  int s_stack = ((n_child - 1) * (depth_max - 1) + n_child);

  int *start_stack = malloc ((sizeof(int)) * s_stack);
  int *end_stack   = malloc ((sizeof(int)) * s_stack);
  PDM_morton_code_t *code_stack = malloc (sizeof(PDM_morton_code_t) * s_stack);

  /* Push root in stack */
  int pos_stack = 0;
  PDM_morton_copy (root, code_stack + pos_stack);
  start_stack[pos_stack] = 0;
  end_stack[pos_stack] = octants->n_nodes;
  pos_stack++;

  PDM_morton_code_t node;
  while (pos_stack > 0) {
    pos_stack--;
    PDM_morton_copy (code_stack[pos_stack], &node);
    int start = start_stack[pos_stack];
    int end   = end_stack[pos_stack];

    if (DEBUG) {
      printf("\npopped {L=%d, X=%d %d %d}, start = %d, end = %d / %d\n",
             node.L, node.X[0], node.X[1], node.X[2], start, end, octants->n_nodes);
    }

    /* Leaf node */
    if (start == end-1) {
      double *_min = _pts_extents + 6*(*n_nodes);
      double *_max = _min + 3;
      _nodes[4*(*n_nodes)] = (int) octants->codes[start].L;
      for (int j = 0; j < 3; j++) {
        _nodes[4*(*n_nodes) + j + 1] = (int) octants->codes[start].X[j];
        _min[j] =  HUGE_VAL;
        _max[j] = -HUGE_VAL;
      }

      _n_pts[*n_nodes] = octants->n_points[start];
      for (int i = 0; i < octants->n_points[start]; i++) {
        const double *p = points + 3*(octants->range[start] + i);
        for (int j = 0; j < 3; j++) {
          _min[j] = PDM_MIN (_min[j], p[j]);
          _max[j] = PDM_MAX (_max[j], p[j]);
        }
      }

      (*n_nodes)++;
      if (DEBUG) printf("  --> add to nodes\n");
    }

    /* Internal node */
    else {
      // check if node is fully covered by leaves
      //   first, check if first descendant leaf is equal to current node
      if (DEBUG) {
        printf("  start = {L=%d, X=%d %d %d}, end = {L=%d, X=%d %d %d}\n", octants->codes[start].L, octants->codes[start].X[0], octants->codes[start].X[1], octants->codes[start].X[2], octants->codes[end-1].L, octants->codes[end-1].X[0], octants->codes[end-1].X[1], octants->codes[end-1].X[2]);
      }

      if (!PDM_morton_a_gtmin_b (octants->codes[start], node)) {
        // now check if last descendant leaf
        _L = octants->codes[end-1].L - node.L;
        for (int j = 0; j < 3; j++) {
          _X[j] = (node.X[j]+1) << _L;
        }

        //if (!PDM_morton_a_gt_b (node, octants->codes[end-1])) {
        if (_X[0] == octants->codes[end-1].X[0]+1 &&
            _X[1] == octants->codes[end-1].X[1]+1 &&
            _X[2] == octants->codes[end-1].X[2]+1) {

          double *_min = _pts_extents + 6*(*n_nodes);
          double *_max = _min + 3;
          _nodes[4*(*n_nodes)] = (int) node.L;
          for (int j = 0; j < 3; j++) {
            _nodes[4*(*n_nodes) + j + 1] = (int) node.X[j];
            _min[j] =  HUGE_VAL;
            _max[j] = -HUGE_VAL;
          }

          _n_pts[*n_nodes] = 0;
          for (int k = start; k < end; k++) {
            _n_pts[*n_nodes] += octants->n_points[k];

            for (int i = 0; i < octants->n_points[k]; i++) {
              const double *p = points + 3*(octants->range[k] + i);
              for (int j = 0; j < 3; j++) {
                _min[j] = PDM_MIN (_min[j], p[j]);
                _max[j] = PDM_MAX (_max[j], p[j]);
              }
            }
          }

          (*n_nodes)++;
          if (DEBUG) printf("  --> add to nodes\n");
          continue;
        }
        else {
          if (DEBUG) printf("  --> node gt end leaf {L=%d, X=%d %d %d}\n", octants->codes[end-1].L, octants->codes[end-1].X[0], octants->codes[end-1].X[1], octants->codes[end-1].X[2]);
        }
      }
      else {
        if (DEBUG) printf("  --> start leaf {L=%d, X=%d %d %d} gtmin node\n", octants->codes[start].L, octants->codes[start].X[0], octants->codes[start].X[1], octants->codes[start].X[2]);
      }

      /* Carry on with children */
      PDM_morton_code_t child_code[n_child];
      PDM_morton_get_children (3,
                               node,
                               child_code);

      int new_start, new_end;
      int prev_end = start;
      int n_select = 0;
      int select[n_child];
      int child_start[8];
      int child_end[8];
      for (int i = 0; i < n_child; i++) {
        /* get start and end of range in list of nodes covered by current child */
        /* new_start <-- first descendant of child in list */
        new_start = prev_end;
        while (new_start < end) {
          if (PDM_morton_ancestor_is (child_code[i], octants->codes[new_start])) {
            break;
          } else if (PDM_morton_a_gt_b (octants->codes[new_start], child_code[i])) {
            /* all the following nodes are clearly not descendants of current child */
            new_start = end+1;
            break;
          }
          new_start++;
        }

        if (new_start > end) {
          /* no need to go further for that child
             because it has no descendants in the node list */
          continue;
        }

        /* new_end <-- next of last descendant of child in list */
        int l = new_start;
        new_end = end;
        while (new_end > l + 1) {
          int m = l + (new_end - l) / 2;
          if (PDM_morton_ancestor_is (child_code[i], octants->codes[m])) {
            l = m;
          } else {
            new_end = m;
          }
        }
        prev_end = new_end;

        if (new_end > new_start) {
          child_start[i] = new_start;
          child_end[i]   = new_end;
          select[n_select++] = i;
        }
      } // End of loop on children

      /* Push selected children in stack */
      for (int i = n_select-1; i >= 0; i--) {
        int ichild = select[i];
        PDM_morton_copy (child_code[ichild], code_stack + pos_stack);
        start_stack[pos_stack] = child_start[ichild];
        end_stack[pos_stack]   = child_end[ichild];
        pos_stack++;
      }
    }
  }


  free (start_stack);
  free (end_stack);
  free (code_stack);

  if (*n_nodes < octants->n_nodes) {
    *nodes       = realloc (*nodes,       sizeof(int)    * (*n_nodes) * 4);
    *n_pts       = realloc (*n_pts,       sizeof(int)    * (*n_nodes));
    *pts_extents = realloc (*pts_extents, sizeof(double) * (*n_nodes) * 6);
  }

  /* Fix extents for empty nodes */
  for (int i = 0; i < *n_nodes; i++) {
    if ((*n_pts)[i] == 0) {
      for (int j = 0; j < 6; j++) {
        (*pts_extents)[6*i + j] = 0.;
      }
    }
  }
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

  int *send_count = PDM_array_zeros_int(n_ranks);
  int *send_shift = malloc(sizeof(int) * (n_ranks+1));

  int *recv_count = malloc(sizeof(int) * n_ranks);
  int *recv_shift = malloc(sizeof(int) * (n_ranks+1));

  int irank = 0;

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

  PDM_array_reset_int(send_count, n_ranks, 0);

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

  int *weight = PDM_array_zeros_int(G->n_nodes);

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


/**
 *
 * \brief Compute the connected parts of an octree
 *
 * \param [inout]   octree
 * \param [in]      neighbours
 *
 */

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

  int *visited = PDM_array_zeros_int(octree->octants->n_nodes);

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

    int *neighbour_rank_n = PDM_array_zeros_int(n_quantile);
    int *neighbour_rank_idx = malloc (sizeof(int) * (n_quantile + 1));


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

    int *recv_neighbour_rank_n = PDM_array_zeros_int(n_quantile);


    PDM_MPI_Request *recv_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);
    PDM_MPI_Request *send_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);

    int *used_ranks = malloc (sizeof(int) * n_ranks);

    PDM_MPI_Alltoall (neighbour_rank_n, n_direction, PDM_MPI_INT,
                      recv_neighbour_rank_n, n_direction, PDM_MPI_INT,
                      octree->comm);

    int *recv_neighbour_rank_idx = PDM_array_new_idx_from_sizes_int(recv_neighbour_rank_n, n_direction * n_ranks);


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

    int *neighbour_rank_n   = PDM_array_zeros_int(n_quantile);
    int *neighbour_rank_idx = malloc (sizeof(int) * (n_quantile + 1));


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

    int *recv_neighbour_rank_n = PDM_array_zeros_int(n_quantile);


    PDM_MPI_Request *recv_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);
    PDM_MPI_Request *send_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);

    int *used_ranks = malloc (sizeof(int) * n_ranks);

    PDM_MPI_Alltoall (neighbour_rank_n, n_direction, PDM_MPI_INT,
                      recv_neighbour_rank_n, n_direction, PDM_MPI_INT,
                      octree->comm);

    int *recv_neighbour_rank_idx = PDM_array_new_idx_from_sizes_int(recv_neighbour_rank_n, n_direction * n_ranks);


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

    recv_rank_ngb_idx = PDM_array_new_idx_from_sizes_int(recv_rank_ngb_n, n_ranks);


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
    assert (PDM_ABS(area[i] - exact_area) <= 1e-15 * exact_area);
  }

  free (area);
}



/**
 *
 * \brief Encode an array of coordinates (truncated to extents)
 *
 * \param [in]   dim       Dimension
 * \param [in]   level     Level in the grid
 * \param [in]   extents   Coordinate extents for normalization (size: dim*2)
 * \param [in]   coords    Coordinates in the grid (interlaced, not normalized)
 * \param [out]  m_code    Array of corresponding Morton codes
 * \param [out]  d         Normalization (dilatation component)
 * \param [out]  s         Normalization (translation component)
 *
 */

static void
_morton_encode_coords
(
 const int          dim,
 PDM_morton_int_t   level,
 const double       extents[],
 size_t             n_coords,
 const double       coords[],
 PDM_morton_code_t  m_code[],
 double             d[3],
 double             s[3]
 )
{
  size_t i, j;
  double n[3];
  double d_max = 0.0;

  PDM_morton_int_t  refinement = 1u << level;

  for (i = 0; i < (size_t)dim; i++) {
    s[i] = extents[i];
    d[i] = extents[i+dim] - extents[i];
    d_max = PDM_MAX(d_max, d[i]);
  }

  for (i = 0; i < (size_t)dim; i++) { /* Reduce effective dimension */
    if (d[i] < d_max * 1e-10)
      d[i] = d_max * 1e-10;
  }

  switch(dim) {

  case 3:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      for (j = 0; j < 3; j++) {
        n[j] = PDM_MAX (0., (coords[i*dim + j] - s[j]) / d[j]);//
        m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
      }
    }
    break;

  case 2:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      for (j = 0; j < 2; j++) {
        n[j] = PDM_MAX (0., (coords[i*dim + j] - s[j]) / d[j]);//
        m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
      }
      m_code[i].X[2] = 0;
    }
    break;

  case 1:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      n[0] = PDM_MAX (0., (coords[i] - s[0]) / d[0]);//
      m_code[i].X[0] = (PDM_morton_int_t) PDM_MIN(floor(n[0]*refinement), refinement - 1);
      m_code[i].X[1] = 0;
      m_code[i].X[2] = 0;
    }
    break;

  default:
    assert(dim > 0 && dim < 4);
    break;
  }
}


static inline int
_box_min_dist2
(
 const int              dim,
 const double *restrict extents,
 const double *restrict coords,
 double       *restrict min_dist2
 )
{
  int inbox = 0;
  double _min_dist2 = 0.;

  for (int i = 0; i < dim; i++) {
    double x = coords[i];
    double xmin = extents[i];
    double xmax = extents[i+dim];

    if (x > xmax) {
      double diff_max = x - xmax;
      _min_dist2 += diff_max * diff_max;
    }
    else if (x < xmin) {
      double diff_min = x - xmin;
      _min_dist2 += diff_min * diff_min;
    }
    else {
      inbox++;
    }
  }

  *min_dist2 = _min_dist2;
  return inbox == dim;
}



static void
_closest_points
(
 const int          n_closest_points,
 const int          dim,
 const double       d[],
 const double       s[],
 const _l_octant_t *octants,
 const PDM_g_num_t *src_g_num,
 const double      *src_coord,
 const int          n_tgt,
 const double      *tgt_coord,
 PDM_g_num_t       *closest_points_g_num,
 double            *closest_points_dist2
 )
{
  if (n_tgt < 1) {
    return;
  }

  /* Deepest common ancestor of all local octants */
  PDM_morton_code_t ancestor;
  PDM_morton_nearest_common_ancestor (octants->codes[0],
                                      octants->codes[octants->n_nodes-1],
                                      &ancestor);

  _min_heap_t *heap = _min_heap_create (octants->n_nodes);

  /* Loop over target points */
  double node_dist2;
  PDM_morton_code_t node_code;
  int node_start, node_end;

  const int n_child = 1 << dim;
  PDM_morton_code_t child_code[8];
  int new_start, new_end, prev_end;
  double child_dist2;

  for (int i_tgt = 0; i_tgt < n_tgt; i_tgt++) {

    const double *point = tgt_coord + dim*i_tgt;
    double *last_closest_pt_dist2 = closest_points_dist2 + (n_closest_points*(i_tgt+1) - 1);

    _min_heap_reset (heap);
    _min_heap_push (heap, ancestor, 0, octants->n_nodes, 0.);

    while (_min_heap_pop (heap, &node_code, &node_start, &node_end, &node_dist2)) {

      if (node_dist2 >= *last_closest_pt_dist2) {
        /* All the nodes in the heap are now farther than the provisional closest point */
        break;
      }

      /* Internal node */
      if (node_end > node_start + 1) {
        /* Carry on with children of popped node */
        PDM_morton_get_children (dim,
                                 node_code,
                                 child_code);
        prev_end = node_start;
        for (int i_child = 0; i_child < n_child; i_child++) {
          /* get start and end of range in list of nodes covered by current child */
          /* new_start <-- first descendant of child in list */
          new_start = prev_end; // end of previous child's range

          while (new_start < node_end) {
            if (PDM_morton_ancestor_is (child_code[i_child],
                                        octants->codes[new_start])) {
              break;
            } else if (PDM_morton_a_gt_b(octants->codes[new_start],
                                         child_code[i_child])) {
              /* all the following nodes are clearly not descendants of current child */
              new_start = node_end + 1;
              break;
            }
            new_start++;
          }

          if (new_start > node_end) {
            /* no need to go further for that child
               because it has no descendants in the node list */
            continue;
          }

          /* new_end <-- next of last descendant of child in list */
          int l = new_start;
          new_end = node_end;
          while (new_end > l + 1) {
            int m = l + (new_end - l) / 2;
            if (PDM_morton_ancestor_is (child_code[i_child],
                                        octants->codes[m])) {
              l = m;
            } else {
              new_end = m;
            }
          }

          prev_end = new_end;
          if (new_end <= new_start) {
            continue;
          }

          child_dist2 = _octant_min_dist2 (dim,
                                           child_code[i_child],
                                           d,
                                           s,
                                           point);
          /* Push child in heap */
          if (child_dist2 < *last_closest_pt_dist2) {
            _min_heap_push (heap, child_code[i_child], new_start, new_end, child_dist2);
          }

        } // End of loop on children
      }

      /* Leaf node */
      else {
        /* inspect source points inside popped leaf */
        for (int i = 0; i < octants->n_points[node_start]; i++) {
          int j = octants->range[node_start] + i;
          double point_dist2 = _pt_to_pt_dist2 (dim,
                                                point,
                                                src_coord + dim*j);
          if (point_dist2 < *last_closest_pt_dist2) {
            _insertion_sort (point_dist2,
                             src_g_num[j],
                             n_closest_points,
                             closest_points_dist2 + n_closest_points*i_tgt,
                             closest_points_g_num + n_closest_points*i_tgt);
          }
        } // End of loop on source points inside popped leaf
      } // End if internal/leaf node

    } // End while heap not empty

  } // End of loop on target points

  _min_heap_free (heap);
}




static void
_closest_points_explicit
(
 const int     n_closest_points,
 _octree_t    *octree,
 const int     i_copied_rank,
 const int     n_tgt,
 const double *tgt_coord,
 PDM_g_num_t  *closest_points_g_num,
 double       *closest_points_dist2
 )
{
  if (n_tgt < 1) {
    return;
  }

  int n_points;
  const _explicit_node_t *nodes;
  const _l_octant_t *leaves;
  const PDM_g_num_t *src_g_num;
  const double      *src_coord;
  if (i_copied_rank < 0) {
    n_points  = octree->n_points;
    nodes     = octree->explicit_nodes;
    leaves    = octree->octants;
    src_g_num = octree->points_gnum;
    src_coord = octree->points;
  } else {
    assert (i_copied_rank < octree->n_copied_ranks);
    n_points  = octree->n_copied_points[i_copied_rank];
    nodes     = octree->copied_explicit_nodes[i_copied_rank];
    leaves    = octree->copied_octants[i_copied_rank];
    src_g_num = octree->copied_points_gnum[i_copied_rank];
    src_coord = octree->copied_points[i_copied_rank];
  }

  if (n_points == 0) {
    return;
  }

  const int dim = octree->dim;
  const int n_child = 1 << dim;
  const int depth_max = 31;
  int s_stack = ((n_child - 1) * (depth_max - 1) + n_child);
  int    *stack_id     = malloc (sizeof(int)    * s_stack);
  double *stack_dist   = malloc (sizeof(double) * s_stack);
  int    *stack_inside = malloc (sizeof(int)    * s_stack);

  int    node_id;
  double node_dist;
  int    node_inside;
  double child_dist[n_child];
  int    child_inside[n_child];

  /* Loop on target points */
  for (int itgt = 0; itgt < n_tgt; itgt++) {

    const double *point = tgt_coord + 3*itgt;
    PDM_g_num_t *_closest_pt_g_num = closest_points_g_num + n_closest_points*itgt;
    double      *_closest_pt_dist2 = closest_points_dist2 + n_closest_points*itgt;
    double *last_closest_pt_dist2 = _closest_pt_dist2 + n_closest_points - 1;

    node_inside = _box_min_dist2 (dim,
                                  nodes[0].pts_extents,
                                  point,
                                  &node_dist);

    if (node_dist >= *last_closest_pt_dist2 && !node_inside) continue;

    /* Push root in stack */
    int pos_stack = 0;
    stack_id[pos_stack]     = 0;
    stack_dist[pos_stack]   = node_dist;
    stack_inside[pos_stack] = node_inside;
    pos_stack++;

    while (pos_stack > 0) {

      pos_stack--;
      node_dist   = stack_dist[pos_stack];
      node_inside = stack_inside[pos_stack];

      if (node_dist < *last_closest_pt_dist2 || node_inside) {
        node_id = stack_id[pos_stack];
        const _explicit_node_t *_node = nodes + node_id;

        /* Leaf node */
        if (_node->leaf_id >= 0) {
          int ileaf = _node->leaf_id;
          /*
           * optimization: if n_points == 1, distance has already been computed
           */
          if (leaves->n_points[ileaf] == 0) {
            _insertion_sort (node_dist,
                             src_g_num[leaves->range[ileaf]],
                             n_closest_points,
                             _closest_pt_dist2,
                             _closest_pt_g_num);
          }

          else {
            for (int i = 0; i < leaves->n_points[ileaf]; i++) {
              int j = leaves->range[ileaf] + i;
              double dist2 = _pt_to_pt_dist2 (dim,
                                              point,
                                              src_coord + dim*j);
              if (dist2 < *last_closest_pt_dist2) {
                _insertion_sort (dist2,
                                 src_g_num[j],
                                 n_closest_points,
                                 _closest_pt_dist2,
                                 _closest_pt_g_num);
              }
            }
          }
        }


        /* Internal node */
        else {

          int n_select = 0;
          for (int i = 0; i < 8; i++) {
            child_dist[i] = HUGE_VAL;
          }
          int selected_id[8];
          for (int i = 0; i < n_child; i++) {
            if (_node->children_id[i] >= 0) {
              const _explicit_node_t *_child = nodes + _node->children_id[i];

              double dist;
              int inside = _box_min_dist2 (dim,
                                           _child->pts_extents,
                                           point,
                                           &dist);

              if (dist < *last_closest_pt_dist2 || inside) {
                int j;
                for (j = n_select; (j > 0) && (child_dist[j-1] > dist); j--) {
                  selected_id[j]  = selected_id[j-1];
                  child_dist[j]   = child_dist[j-1];
                  child_inside[j] = child_inside[j-1];
                }

                selected_id[j]  = _node->children_id[i];
                child_dist[j]   = dist;
                child_inside[j] = inside;

                n_select++;
              }
            }
          }

          for (int i = n_select-1; i >= 0; i--) {
            stack_id[pos_stack]     = selected_id[i];
            stack_dist[pos_stack]   = child_dist[i];
            stack_inside[pos_stack] = child_inside[i];
            pos_stack++;
          }
        }
      }

    } // End of while loop

  } // End of loop on target points

  free (stack_id);
  free (stack_dist);
  free (stack_inside);
}



static void
_single_closest_point
(
 const int          dim,
 const double       d[],
 const double       s[],
 const _l_octant_t *octants,
 const PDM_g_num_t *src_g_num,
 const double      *src_coord,
 const int          n_tgt,
 const double      *tgt_coord,
 PDM_g_num_t       *closest_point_g_num,
 double            *closest_point_dist2
 )
{
  if (n_tgt < 1) {
    return;
  }

  PDM_morton_code_t ancestor;
  PDM_morton_nearest_common_ancestor (octants->codes[0],
                                      octants->codes[octants->n_nodes-1],
                                      &ancestor);

  const int n_child = 1 << dim;
  const int depth_max = 31;
  int s_stack = ((n_child - 1) * (depth_max - 1) + n_child);

  int *stack_start   = malloc (sizeof(int) * s_stack);
  int *stack_end     = malloc (sizeof(int) * s_stack);
  double *stack_dist = malloc (sizeof(double) * s_stack);
  PDM_morton_code_t *stack_code = malloc (sizeof(PDM_morton_code_t) * s_stack);

  int node_start;
  int node_end;
  double node_dist;
  PDM_morton_code_t node_code;
  int child_start[n_child];
  int child_end[n_child];
  double child_dist[n_child];
  PDM_morton_code_t child_code[n_child];

  //int count = 0;
  /* Loop on target points */
  for (int itgt = 0; itgt < n_tgt; itgt++) {
    const double *point = tgt_coord + 3*itgt;

    node_dist = _octant_min_dist2 (dim,
                                   ancestor,
                                   d,
                                   s,
                                   point);

    if (node_dist >= closest_point_dist2[itgt]) continue;

    /* Push ancestor in stack */
    int pos_stack = 0;
    PDM_morton_copy (ancestor, stack_code + pos_stack);
    stack_start[pos_stack] = 0;
    stack_end[pos_stack]   = octants->n_nodes;
    stack_dist[pos_stack]  = node_dist;
    pos_stack++;

    while (pos_stack > 0) {
      //count++;
      node_dist = stack_dist[--pos_stack];

      if (node_dist < closest_point_dist2[itgt]) {
        node_start = stack_start[pos_stack];
        node_end   = stack_end[pos_stack];
        PDM_morton_copy (stack_code[pos_stack], &node_code);

        /* Leaf node */
        if (node_start == node_end-1) {
          for (int i = 0; i < octants->n_points[node_start]; i++) {
            int j = octants->range[node_start] + i;
            double dist2 = _pt_to_pt_dist2 (dim,
                                            point,
                                            src_coord + dim*j);
            if (dist2 < closest_point_dist2[itgt]) {
              closest_point_dist2[itgt] = dist2;
              closest_point_g_num[itgt] = src_g_num[j];
            }
          }
        }

        /* Internal node */
        else {
          PDM_morton_get_children (dim,
                                   node_code,
                                   child_code);
          int prev_end = node_start;
          for (int i = 0; i < n_child; i++) {
            /* get start and end of range in list of nodes covered by current child */
            /* new_start <-- first descendant of child in list */
            int new_start = prev_end;
            while (new_start < node_end) {
              if (PDM_morton_ancestor_is (child_code[i], octants->codes[new_start])) {
                break;
              } else if (PDM_morton_a_gt_b (octants->codes[new_start], child_code[i])) {
                /* all the following nodes are clearly not descendants of current child */
                new_start = node_end+1;
                break;
              }
              new_start++;
            }

            if (new_start > node_end) {
              /* no need to go further for that child
                 because it has no descendants in the node list */
              child_dist[i] = HUGE_VAL;
              continue;
            }

            child_start[i] = new_start;

            /* new_end <-- next of last descendant of child in list */
            int new_end = node_end;
            while (new_end > new_start + 1) {
              int m = new_start + (new_end - new_start) / 2;
              if (PDM_morton_ancestor_is (child_code[i], octants->codes[m])) {
                new_start = m;
              } else {
                new_end = m;
              }
            }

            prev_end = new_end;
            child_end[i] = new_end;

            //if (child_end[i] > child_start[i] &&
            //    octants->range[child_start[i]] < octants->range[child_end[i]]) {
            if (child_end[i] > child_start[i]) {
              child_dist[i] = _octant_min_dist2 (dim,
                                                 child_code[i],
                                                 d,
                                                 s,
                                                 point);
            } else {
              child_dist[i] = HUGE_VAL;
            }

          } // End of loop on children

          int child_order[8] = {0, 1, 2, 3, 4, 5, 6, 7};
          PDM_sort_double (child_dist,
                           child_order,
                           n_child);

          for (int i = n_child-1; i >= 0; i--) {
            if (child_dist[i] < closest_point_dist2[itgt]) {
              int ichild = child_order[i];

              /* Empty leaf node */
              //if ((child_start[ichild] == child_end[ichild]-1 &&
              //     octants->n_points[child_start[ichild]] == 0) ||
              //    (octants->range[child_start[ichild]] == octants->range[child_end[ichild]-1])) {
              if ((child_start[ichild] == child_end[ichild]-1 &&
                   octants->n_points[child_start[ichild]] == 0)) {
                continue;
              }
              /*              if (octants->range[child_start[ichild]] == octants->range[child_end[ichild]]) {
continue;
}*/

              /* Push child in stack */
              PDM_morton_copy (child_code[ichild], stack_code + pos_stack);
              stack_start[pos_stack] = child_start[ichild];
              stack_end[pos_stack]   = child_end[ichild];
              stack_dist[pos_stack]  = child_dist[i];
              pos_stack++;
            }
          }
        }
      }

    } // End while stack not empty

  } // End of loop on target points

  //if (n_tgt != 0) printf("n nodes per point = %d\n", count / n_tgt);

  free (stack_start);
  free (stack_end);
  free (stack_dist);
  free (stack_code);
}



static void
_single_closest_point_explicit
(
 _octree_t    *octree,
 const int     i_copied_rank,
 const int     n_tgt,
 const double *tgt_coord,
 PDM_g_num_t  *closest_point_g_num,
 double       *closest_point_dist2
 )
{
  if (n_tgt < 1) {
    return;
  }

  int n_points;
  const _explicit_node_t *nodes;
  const _l_octant_t *leaves;
  const PDM_g_num_t *src_g_num;
  const double      *src_coord;
  if (i_copied_rank < 0) {
    n_points  = octree->n_points;
    nodes     = octree->explicit_nodes;
    leaves    = octree->octants;
    src_g_num = octree->points_gnum;
    src_coord = octree->points;
  } else {
    assert (i_copied_rank < octree->n_copied_ranks);
    n_points  = octree->n_copied_points[i_copied_rank];
    nodes     = octree->copied_explicit_nodes[i_copied_rank];
    leaves    = octree->copied_octants[i_copied_rank];
    src_g_num = octree->copied_points_gnum[i_copied_rank];
    src_coord = octree->copied_points[i_copied_rank];
  }

  if (n_points == 0) {
    return;
  }

  const int dim = octree->dim;
  const int n_child = 1 << dim;
  const int depth_max = 31;
  int s_stack = ((n_child - 1) * (depth_max - 1) + n_child);
  int    *stack_id     = malloc (sizeof(int)    * s_stack);
  double *stack_dist   = malloc (sizeof(double) * s_stack);
  int    *stack_inside = malloc (sizeof(int)    * s_stack);

  int    node_id;
  double node_dist;
  int    node_inside;
  double child_dist[n_child];
  int    child_inside[n_child];
  //int count = 0;
  /* Loop on target points */
  for (int itgt = 0; itgt < n_tgt; itgt++) {
    const double *point = tgt_coord + 3*itgt;

    node_inside = _box_min_dist2 (dim,
                                  nodes[0].pts_extents,
                                  point,
                                  &node_dist);

    if (node_dist >= closest_point_dist2[itgt] && !node_inside) continue;

    /* Push root in stack */
    int pos_stack = 0;
    stack_id[pos_stack]     = 0;
    stack_dist[pos_stack]   = node_dist;
    stack_inside[pos_stack] = node_inside;
    pos_stack++;

    while (pos_stack > 0) {
      //count++;
      pos_stack--;
      node_dist   = stack_dist[pos_stack];
      node_inside = stack_inside[pos_stack];

      if (node_dist < closest_point_dist2[itgt] || node_inside) {
        node_id = stack_id[pos_stack];
        const _explicit_node_t *_node = nodes + node_id;

        /* Leaf node */
        if (_node->leaf_id >= 0) {
          int ileaf = _node->leaf_id;
          /*
           * optimization: if n_points == 1, distance has already been computed
           */
          if (leaves->n_points[ileaf] == 0) {
            closest_point_dist2[itgt] = node_dist;
            closest_point_g_num[itgt] = src_g_num[leaves->range[ileaf]];
          }

          else {
            for (int i = 0; i < leaves->n_points[ileaf]; i++) {
              int j = leaves->range[ileaf] + i;
              double dist2 = _pt_to_pt_dist2 (dim,
                                              point,
                                              src_coord + dim*j);
              if (dist2 < closest_point_dist2[itgt]) {
                closest_point_dist2[itgt] = dist2;
                closest_point_g_num[itgt] = src_g_num[j];
              }
            }
          }
        }

        /* Internal node */
        else {

          int n_select = 0;
          for (int i = 0; i < 8; i++) {
            child_dist[i] = HUGE_VAL;
          }
          int selected_id[8];
          for (int i = 0; i < n_child; i++) {
            if (_node->children_id[i] >= 0) {
              const _explicit_node_t *_child = nodes + _node->children_id[i];

              double dist;
              int inside = _box_min_dist2 (dim,
                                           _child->pts_extents,
                                           point,
                                           &dist);

              if (dist < closest_point_dist2[itgt] || inside) {
                int j;
                for (j = n_select; (j > 0) && (child_dist[j-1] > dist); j--) {
                  selected_id[j]  = selected_id[j-1];
                  child_dist[j]   = child_dist[j-1];
                  child_inside[j] = child_inside[j-1];
                }

                selected_id[j]  = _node->children_id[i];
                child_dist[j]   = dist;
                child_inside[j] = inside;

                n_select++;
              }
            }
          }

          for (int i = n_select-1; i >= 0; i--) {
            stack_id[pos_stack]     = selected_id[i];
            stack_dist[pos_stack]   = child_dist[i];
            stack_inside[pos_stack] = child_inside[i];
            pos_stack++;
          }
        }
      }

    } // End of while loop

  } // End of loop on target points

  //if (n_tgt != 0) printf("n nodes per point = %d\n", count / n_tgt);

  free (stack_id);
  free (stack_dist);
  free (stack_inside);
}




static int
_intersect_node_box
(
 const int                dim,
 const PDM_morton_code_t  node,
 const PDM_morton_code_t  box_min,
 const PDM_morton_code_t  box_max,
 int                     *inside
 )
{
  const int DEBUG = 0;
  *inside = 1;

  assert (box_min.L >= node.L);

  const PDM_morton_int_t level_diff = box_min.L - node.L;

  const PDM_morton_int_t side = 1 << level_diff;

  for (int i = 0; i < dim; i++) {
    PDM_morton_int_t xmin = side * node.X[i];
    PDM_morton_int_t xmax = xmin + side;

    if (xmin > box_max.X[i]+1 || xmax < box_min.X[i]) {
      if (DEBUG) {
       printf("\t not intersecting\n");
      }
      return 0;
    } else if (xmin < box_min.X[i] || xmax > box_max.X[i]+1) {
      *inside = 0;
    };
  }

  if (DEBUG) {
    printf("\t intersecting\n");
  }

  return 1;
}


static void
_points_inside_boxes
(
 const _octree_t          *octree,
 const int                 i_copied_rank,
 const int                 n_box,
 const double              box_extents[],
 const PDM_morton_code_t   box_codes[],
 int                     **pts_idx,
 int                     **pts_l_num
 )
{
  *pts_idx = malloc (sizeof(int) * (n_box + 1));
  int *_pts_idx = *pts_idx;
  _pts_idx[0] = 0;

  if (n_box < 1) {
    *pts_l_num = malloc (sizeof(int) * _pts_idx[n_box]);
    return;
  }

  const _l_octant_t *octants;
  const double      *points;
  if (i_copied_rank < 0) {
    octants = octree->octants;
    points  = octree->points;
  } else {
    assert (i_copied_rank < octree->n_copied_ranks);
    octants = octree->copied_octants[i_copied_rank];
    points  = octree->copied_points[i_copied_rank];
  }


  /* Deepest common ancestor of all octants */
  PDM_morton_code_t root;
  if (octants->n_nodes > 0) {
    PDM_morton_nearest_common_ancestor (octants->codes[0],
                                        octants->codes[octants->n_nodes - 1],
                                        &root);
  }
  else {
    PDM_array_reset_int (_pts_idx+1, n_box, 0);
    *pts_l_num = malloc (sizeof(int) * _pts_idx[n_box]);
    return;
  }


  int tmp_size = 4 * n_box;
  *pts_l_num = malloc (sizeof(int) * tmp_size);
  int *_pts_l_num = *pts_l_num;


  const int n_child = 8;
  const int depth_max = 31;
  int s_stack = ((n_child - 1) * (depth_max - 1) + n_child);

  int *start_stack = malloc ((sizeof(int)) * s_stack);
  int *end_stack   = malloc ((sizeof(int)) * s_stack);
  PDM_morton_code_t *code_stack = malloc (sizeof(PDM_morton_code_t) * s_stack);

  PDM_morton_code_t node;
  int node_inside_box;
  int intersect;

  for (int ibox = 0; ibox < n_box; ibox++) {
    int DEBUG = 0;//(ibox == 201);
    _pts_idx[ibox+1] = _pts_idx[ibox];

    const double *box_min = box_extents + 6*ibox;
    const double *box_max = box_min + 3;
    const PDM_morton_code_t *box_code_min = box_codes + 2*ibox;
    const PDM_morton_code_t *box_code_max = box_code_min + 1;
    if (DEBUG) {
      printf("box %d: min = %f %f %f, max = %f %f %f\n",
             ibox,
             box_min[0], box_min[1], box_min[2],
             box_max[0], box_max[1], box_max[2]);
      printf("center: %f %f %f, length = %f %f %f\n",
             0.5*(box_max[0] + box_min[0]),
             0.5*(box_max[1] + box_min[1]),
             0.5*(box_max[2] + box_min[2]),
             box_max[0] - box_min[0],
             box_max[1] - box_min[1],
             box_max[2] - box_min[2]);
    }

    intersect = _intersect_node_box (3,
                                     root,
                                     *box_code_min,
                                     *box_code_max,
                                     &node_inside_box);

    if (!intersect) {
      if (DEBUG) {
        printf("box %d does not intersect root node\n", ibox);
      }
      continue;
    }

    if (node_inside_box) {
      /* The box must contain all points */
      int new_size = tmp_size;
      for (int i = 0; i < octants->n_nodes; i++) {
        new_size += octants->n_points[i];
      }

      if (tmp_size <= new_size) {
        tmp_size = PDM_MAX (2*tmp_size, new_size);
        *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
        _pts_l_num = *pts_l_num;

        for (int i = 0; i < octants->n_nodes; i++) {
          for (int j = 0; j < octants->n_points[i]; j++) {
            _pts_l_num[_pts_idx[ibox+1]++] = octants->range[i] + j;
          }
        }
      }

      continue;
    }

    /* Push root in stack */
    int pos_stack = 0;
    PDM_morton_copy (root, code_stack + pos_stack);
    start_stack[pos_stack] = 0;
    end_stack[pos_stack] = octants->n_nodes;
    pos_stack++;

    while (pos_stack > 0) {

      pos_stack--;
      PDM_morton_copy (code_stack[pos_stack], &node);
      int start = start_stack[pos_stack];
      int end   = end_stack[pos_stack];

      /* Leaf node */
      if (start == end-1) {
        for (int i = 0; i < octants->n_points[start]; i++) {
          int ipt = octants->range[start] + i;
          const double *_pt = points + ipt * 3;

          int pt_inside_box = 1;
          for (int idim = 0; idim < 3; idim++) {
            if (_pt[idim] < box_min[idim] || _pt[idim] > box_max[idim]) {
              pt_inside_box = 0;
              break;
            }
          }

          if (pt_inside_box) {
            if (_pts_idx[ibox+1] >= tmp_size) {
              tmp_size = PDM_MAX (2*tmp_size, _pts_idx[ibox+1] + 1);
              *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
              _pts_l_num = *pts_l_num;
            }

            _pts_l_num[_pts_idx[ibox+1]++] = ipt;
          }
        } // End of loop on points inside leaf
      }

      /* Internal node */
      else {
        PDM_morton_code_t child_code[n_child];
        PDM_morton_get_children (3,
                                 node,
                                 child_code);
        int new_start, new_end;
        int prev_end = start;
        for (int i = 0; i < n_child; i++) {
          /* get start and end of range in list of nodes covered by current child */
          /* new_start <-- first descendant of child in list */
          new_start = prev_end;
          while (new_start < end) {
            if (PDM_morton_ancestor_is (child_code[i], octants->codes[new_start])) {
              break;
            } else if (PDM_morton_a_gt_b (octants->codes[new_start], child_code[i])) {
              /* all the following nodes are clearly not descendants of current child */
              new_start = end+1;
              break;
            }
            new_start++;
          }

          if (new_start > end) {
            /* no need to go further for that child
               because it has no descendants in the node list */
            continue;
          }

          /* new_end <-- next of last descendant of child in list */
          int l = new_start;
          new_end = end;
          while (new_end > l + 1) {
            int m = l + (new_end - l) / 2;
            if (PDM_morton_ancestor_is (child_code[i], octants->codes[m])) {
              l = m;
            } else {
              new_end = m;
            }
          }
          prev_end = new_end;

          if (new_end > new_start) {
            if (DEBUG) {
              int n_points = octants->range[new_end-1] + octants->n_points[new_end-1] - octants->range[new_start];
              printf("child %d: L=%u, X=%u %u %u, start = %d, end = %d, range = %d, n_points = %d\n",
                     i,
                     child_code[i].L,
                     child_code[i].X[0],
                     child_code[i].X[1],
                     child_code[i].X[2],
                     new_start,
                     new_end,
                     octants->range[new_start],
                     n_points);
            }

            intersect = _intersect_node_box (3,
                                             child_code[i],
                                             *box_code_min,
                                             *box_code_max,
                                             &node_inside_box);
            if (DEBUG) {
              printf("intersect = %d\n", intersect);
            }

            if (intersect) {
              if (node_inside_box) {
                /* The box must contain all points contained in descendants of child */
                int new_size = _pts_idx[ibox+1];
                for (int j = new_start; j < new_end; j++) {
                  new_size += octants->n_points[j];
                }

                if (tmp_size <= new_size) {
                  tmp_size = PDM_MAX (2*tmp_size, new_size);
                  *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
                  _pts_l_num = *pts_l_num;
                }

                for (int j = new_start; j < new_end; j++) {
                  for (int k = 0; k < octants->n_points[j]; k++) {
                    _pts_l_num[_pts_idx[ibox+1]++] = octants->range[j] + k;
                  }
                }
              }

              else {
                /* Push child in stack */
                PDM_morton_copy (child_code[i], code_stack + pos_stack);
                start_stack[pos_stack] = new_start;
                end_stack[pos_stack]   = new_end;
                pos_stack++;
              }
            }
          }
        } // End of loop on child
      }

    }

  } // End loop on boxes

  free (start_stack);
  free (end_stack);
  free (code_stack);

  *pts_l_num = realloc (*pts_l_num, sizeof(int) * _pts_idx[n_box]);
}



static int
_intersect_node_box_explicit
(
 int           dim,
 const double *node_extents,
 const double *box_extents,
 int          *inside
 )
{
  *inside = 1;

  for (int i = 0; i < dim; i++) {
    if (node_extents[i]   > box_extents[i+3] ||
        node_extents[i+3] < box_extents[i]) {
      return 0;
    }
    else if (node_extents[i]  < box_extents[i] ||
             node_extents[i+3] > box_extents[i+3]) {
      *inside = 0;
    }
  }

  return 1;
}




static void
_points_inside_boxes_explicit
(
 const _octree_t          *octree,
 const int                 i_copied_rank,
 const int                 n_box,
 const double              box_extents[],
 // const PDM_morton_code_t   box_codes[],
 int                     **pts_idx,
 int                     **pts_l_num
 )
{
  *pts_idx = malloc (sizeof(int) * (n_box + 1));
  int *_pts_idx = *pts_idx;
  _pts_idx[0] = 0;

  if (n_box < 1) {
    *pts_l_num = malloc (sizeof(int) * _pts_idx[n_box]);
    return;
  }

  int n_nodes;
  const _explicit_node_t *nodes;
  const double *points;
  if (i_copied_rank < 0) {
    n_nodes = octree->n_explicit_nodes;
    nodes   = octree->explicit_nodes;
    points  = octree->points;
  } else {
    assert (i_copied_rank < octree->n_copied_ranks);
    n_nodes = octree->n_copied_explicit_nodes[i_copied_rank];
    nodes   = octree->copied_explicit_nodes[i_copied_rank];
    points  = octree->copied_points[i_copied_rank];
  }

  if (n_nodes == 0 || (n_nodes > 0 && nodes[0].n_points == 0)) {
    PDM_array_reset_int (_pts_idx+1, n_box, 0);
    *pts_l_num = malloc (sizeof(int) * _pts_idx[n_box]);
    return;
  }


  int tmp_size = 4 * n_box;
  *pts_l_num = malloc (sizeof(int) * tmp_size);
  int *_pts_l_num = *pts_l_num;


  const int dim = octree->dim;
  const int n_child = 1 << dim;
  const int depth_max = 31;
  int s_stack = ((n_child - 1) * (depth_max - 1) + n_child);
  int *stack_id = malloc (sizeof(int)    * s_stack);
  int node_inside_box;
  int intersect;

  for (int ibox = 0; ibox < n_box; ibox++) {
    int DEBUG = 0;//(ibox == 201);
    _pts_idx[ibox+1] = _pts_idx[ibox];

    const double *_box_extents = box_extents + 6*ibox;
    const double *box_min = box_extents + 6*ibox;
    const double *box_max = box_min + 3;
    //const PDM_morton_code_t *box_code_min = box_codes + 2*ibox;
    //const PDM_morton_code_t *box_code_max = box_code_min + 1;
    if (DEBUG) {
      printf("box %d: min = %f %f %f, max = %f %f %f\n",
             ibox,
             box_min[0], box_min[1], box_min[2],
             box_max[0], box_max[1], box_max[2]);
      printf("center: %f %f %f, length = %f %f %f\n",
             0.5*(box_max[0] + box_min[0]),
             0.5*(box_max[1] + box_min[1]),
             0.5*(box_max[2] + box_min[2]),
             box_max[0] - box_min[0],
             box_max[1] - box_min[1],
             box_max[2] - box_min[2]);
    }

    intersect = _intersect_node_box_explicit (3,
                                              nodes[0].pts_extents,
                                              _box_extents,
                                              &node_inside_box);

    if (!intersect) {
      if (DEBUG) {
        printf("box %d does not intersect root node\n", ibox);
      }
      continue;
    }

    if (node_inside_box) {
      /* The box must contain all points */
      int new_size = tmp_size + nodes[0].n_points;

      if (tmp_size <= new_size) {
        tmp_size = PDM_MAX (2*tmp_size, new_size);
        *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
        _pts_l_num = *pts_l_num;

        for (int j = 0; j < nodes[0].n_points; j++) {
          _pts_l_num[_pts_idx[ibox+1]++] = nodes[0].range + j;
        }
      }

      continue;
    }

    /* Push root in stack */
    int pos_stack = 0;
    stack_id[pos_stack++] = 0;

    while (pos_stack > 0) {
      int node_id = stack_id[--pos_stack];
      const _explicit_node_t *_node = nodes + node_id;
      if (DEBUG) {
        printf("node %d : L=%u, X=%u %u %u, range=%d, n_points=%d, leaf_id=%d\n",
               node_id,
               _node->code.L,
               _node->code.X[0],
               _node->code.X[1],
               _node->code.X[2],
               _node->range,
               _node->n_points,
               _node->leaf_id);
      }

      /* Leaf node */
      if (_node->leaf_id >= 0) {
        for (int i = 0; i < _node->n_points; i++) {
          int ipt = _node->range + i;
          const double *_pt = points + ipt * 3;

          int pt_inside_box = 1;
          for (int idim = 0; idim < 3; idim++) {
            if (_pt[idim] < box_min[idim] || _pt[idim] > box_max[idim]) {
              pt_inside_box = 0;
              break;
            }
          }

          if (pt_inside_box) {
            if (_pts_idx[ibox+1] >= tmp_size) {
              tmp_size = PDM_MAX (2*tmp_size, _pts_idx[ibox+1] + 1);
              *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
              _pts_l_num = *pts_l_num;
            }

            _pts_l_num[_pts_idx[ibox+1]++] = ipt;
          }
        } // End of loop on points inside leaf
      }

      /* Internal node */
      else {
        for (int i = 0; i < n_child; i++) {
          if (_node->children_id[i] < 0) {
            continue;
          }
          const _explicit_node_t *_child = nodes + _node->children_id[i];

          if (DEBUG) {
            printf("child %d: id=%d, L=%u, X=%u %u %u, range=%d, n_points=%d, leaf_id=%d\n",
                   i,
                   _node->children_id[i],
                   _child->code.L,
                   _child->code.X[0],
                   _child->code.X[1],
                   _child->code.X[2],
                   _child->range,
                   _child->n_points,
                   _child->leaf_id);
          }

          intersect = _intersect_node_box_explicit (3,
                                                    _child->pts_extents,
                                                    _box_extents,
                                                    &node_inside_box);

          if (DEBUG) {
            printf("intersect = %d\n", intersect);
          }

          if (intersect) {
            if (node_inside_box) {
              /* The box must contain all points */
              int new_size = tmp_size + _child->n_points;

              if (tmp_size <= new_size) {
                tmp_size = PDM_MAX (2*tmp_size, new_size);
                *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
                _pts_l_num = *pts_l_num;

                for (int j = 0; j < _child->n_points; j++) {
                  _pts_l_num[_pts_idx[ibox+1]++] = _child->range + j;
                }
              }
            }

            else {
              /* Push child in stack */
              stack_id[pos_stack++] = _node->children_id[i];
            }
          }
        } // End of loop on children
      }
    } // End of while loop

  } // End of loop on boxes

  free (stack_id);

  *pts_l_num = realloc (*pts_l_num, sizeof(int) * _pts_idx[n_box]);
}






static void
_build_explicit_nodes
(
 _octree_t *octree
 )
{
  int DEBUG = 0;

  int i_rank;
  PDM_MPI_Comm_rank (octree->comm, &i_rank);
  if (i_rank == 0) printf("BUILD OCTREE EXPLICIT NODES\n");

  _l_octant_t *octants = octree->octants;
  const int dim = octants->dim;



  int tmp_size = 2*octants->n_nodes;
  octree->explicit_nodes = malloc (sizeof(_explicit_node_t) * tmp_size);
  octree->n_explicit_nodes = 0;

  if (octants->n_nodes == 0 || octree->n_points == 0) return;

  _explicit_node_t *node = octree->explicit_nodes;

  PDM_morton_nearest_common_ancestor (octants->codes[0],
                                      octants->codes[octants->n_nodes-1],
                                      &(node->code));
  node->ancestor_id = -1;
  node->range = octants->range[0];
  node->n_points = octree->n_points;
  for (int j = 0; j < 3; j++) {
    node->pts_extents[j]   =  HUGE_VAL;
    node->pts_extents[j+3] = -HUGE_VAL;
  }
  if (octants->n_nodes == 1) {
    node->leaf_id = 0;
    octree->n_explicit_nodes = 1;
    for (int i = 0; i < octree->n_points; i++) {
      for (int j = 0; j < 3; j++) {
        node->pts_extents[j]   = PDM_MIN (node->pts_extents[j],   octree->points[3*i+j]);
        node->pts_extents[j+3] = PDM_MAX (node->pts_extents[j+3], octree->points[3*i+j]);
      }
    }
    return;
  } else {
    node->leaf_id = -1;
  }
  octree->n_explicit_nodes++;

  const int n_child = 1 << dim;
  const int depth_max = 31;
  int s_stack = ((n_child - 1) * (depth_max - 1) + n_child);

  int *stack_start = malloc (sizeof(int) * s_stack);
  int *stack_end   = malloc (sizeof(int) * s_stack);
  int *stack_id    = malloc (sizeof(int) * s_stack);

  /* Push ancestor in stack */
  int pos_stack = 0;
  if (octants->n_nodes > 1) {
    stack_id[pos_stack]    = 0;
    stack_start[pos_stack] = 0;
    stack_end[pos_stack]   = octants->n_nodes;
    pos_stack++;
  }


  int node_start;
  int node_end;
  int node_id;
  int child_id;
  PDM_morton_code_t child_code[n_child];

  while (pos_stack > 0) {
    pos_stack--;
    node_id    = stack_id[pos_stack];
    node_start = stack_start[pos_stack];
    node_end   = stack_end[pos_stack];

    if (octree->n_explicit_nodes + n_child >= tmp_size) {
      tmp_size = PDM_MAX (2*tmp_size, octree->n_explicit_nodes + n_child + 1);
      octree->explicit_nodes = realloc (octree->explicit_nodes,
                                        sizeof(_explicit_node_t) * tmp_size);
    }
    node = octree->explicit_nodes + node_id;

    node = octree->explicit_nodes + node_id;
    PDM_morton_get_children (dim,
                             node->code,
                             child_code);
    int new_start, new_end;
    int prev_end = node_start;
    for (int i = 0; i < n_child; i++) {
      /* get start and end of range in list of nodes covered by current child */
      /* new_start <-- first descendant of child in list */
      new_start = prev_end;
      while (new_start < node_end) {
        if (PDM_morton_ancestor_is (child_code[i], octants->codes[new_start])) {
          break;
        } else if (PDM_morton_a_gt_b (octants->codes[new_start], child_code[i])) {
          /* all the following nodes are clearly not descendants of current child */
          new_start = node_end+1;
          break;
        }
        new_start++;
      }

      if (new_start >= node_end) {
        // no need to go further for that child
        //   because it has no descendants in the node list
        // (descendant leaves are in another rank)
        node->children_id[i] = -1;
        continue;
      }

      /* new_end <-- next of last descendant of child in list */
      int l = new_start;
      new_end = node_end;
      while (new_end > l + 1) {
        int m = l + (new_end - l) / 2;
        if (PDM_morton_ancestor_is (child_code[i], octants->codes[m])) {
          l = m;
        } else {
          new_end = m;
        }
      }
      prev_end = new_end;

      int n_points = octants->range[new_end-1]
        + octants->n_points[new_end-1]
        - octants->range[new_start];

      if (n_points == 0) {
        node->children_id[i] = -1;
        continue;
      }

      child_id = octree->n_explicit_nodes++;
      _explicit_node_t *child = octree->explicit_nodes + child_id;
      node->children_id[i] = child_id;
      PDM_morton_copy (child_code[i], &(child->code));
      child->ancestor_id = node_id;
      child->range       = octants->range[new_start];
      child->n_points    = n_points;

      for (int j = 0; j < 3; j++) {
        child->pts_extents[j]   =  HUGE_VAL;
        child->pts_extents[j+3] = -HUGE_VAL;
      }

      /* Leaf node */
      if (new_start == new_end-1) {
        child->leaf_id = new_start;

        /* Compute leaf extents */
        for (int k = 0; k < n_points; k++) {
          double *point = octree->points + 3*(octants->range[new_start] + k);
          for (int j = 0; j < 3; j++) {
            child->pts_extents[j]   = PDM_MIN (child->pts_extents[j],   point[j]);
            child->pts_extents[j+3] = PDM_MAX (child->pts_extents[j+3], point[j]);
          }
        }

        /* Propagate extents to ancestors */
        int ancestor_id = node_id;
        while (ancestor_id >= 0) {
          _explicit_node_t *ancestor = octree->explicit_nodes + ancestor_id;
          for (int j = 0; j < 3; j++) {
            ancestor->pts_extents[j]   = PDM_MIN (child->pts_extents[j],
                                                  ancestor->pts_extents[j]);
            ancestor->pts_extents[j+3] = PDM_MAX (child->pts_extents[j+3],
                                                  ancestor->pts_extents[j+3]);
          }
          ancestor_id = ancestor->ancestor_id;
        }
      }

      /* Internal node */
      else {
        child->leaf_id = -1;

        /* Push child in stack */
        stack_id[pos_stack]    = child_id;
        stack_start[pos_stack] = new_start;
        stack_end[pos_stack]   = new_end;
        pos_stack++;
      }
      if (DEBUG) {
        printf("node %d : L=%d, X=%d %d %d, start=%d, end=%d, ancestor=%d, range=%d, n_points=%d, leaf_id=%d\n",
               child_id,
               child->code.L,
               child->code.X[0],
               child->code.X[1],
               child->code.X[2],
               new_start,
               new_end,
               child->ancestor_id,
               child->range,
               child->n_points,
               child->leaf_id);
      }

    } // End of loop on children
  }
  free (stack_start);
  free (stack_end);
  free (stack_id);

  octree->explicit_nodes = realloc (octree->explicit_nodes,
                                    sizeof(_explicit_node_t) * octree->n_explicit_nodes);
}





static void
_compute_rank_extents
(
 _octree_t *octree
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (octree->comm, &i_rank);
  PDM_MPI_Comm_size (octree->comm, &n_rank);

  int dim = octree->dim;
  int dim2 = 2*dim;
  double *local_extents = (double *) malloc (sizeof(double) * dim2);
  for (int j = 0; j < dim; j++) {
    local_extents[j]       =  HUGE_VAL;
    local_extents[dim + j] = -HUGE_VAL;
  }
  for (int i = 0; i < octree->n_points; i++) {
    for (int j = 0; j < dim; j++) {
      local_extents[j]       = PDM_MIN (local_extents[j],       octree->points[dim*i + j]);
      local_extents[dim + j] = PDM_MAX (local_extents[dim + j], octree->points[dim*i + j]);
    }
  }

  double *all_extents = (double *) malloc (sizeof(double) * dim2 * n_rank);
  PDM_MPI_Allgather (local_extents, 2*dim, PDM_MPI_DOUBLE,
                     all_extents,   2*dim, PDM_MPI_DOUBLE,
                     octree->comm);

  free (local_extents);


  int _n_pts = octree->n_points;
  int *n_pts_rank = (int *) malloc (sizeof(int) * n_rank);
  PDM_MPI_Allgather (&_n_pts,     1, PDM_MPI_INT,
                     n_pts_rank,  1, PDM_MPI_INT,
                     octree->comm);

  int n_used_rank = 0;
  for (int i = 0; i < n_rank; i++) {
    if (n_pts_rank[i] > 0) {
      n_used_rank++;
    }
  }

  octree->n_used_rank = n_used_rank;
  octree->used_rank = (int *) malloc (sizeof(int) * n_used_rank);
  octree->used_rank_extents = (double *) malloc (sizeof(double) * dim2 * n_used_rank);
  PDM_g_num_t *gnum_proc = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * octree->n_used_rank);

  n_used_rank = 0;
  for (int i = 0; i < n_rank; i++) {
    if (n_pts_rank[i] > 0) {
      octree->used_rank[n_used_rank] = i;
      gnum_proc[n_used_rank] = n_used_rank + 1;
      memcpy (octree->used_rank_extents + dim2*n_used_rank,
              all_extents + dim2*i,
              sizeof(double)*dim2);
      n_used_rank++;
    }
  }

  free (n_pts_rank);
  free (all_extents);


  // Box tree...
  const int n_info_location = 3;
  int *init_location_proc = PDM_array_zeros_int(n_info_location * octree->n_used_rank);

  PDM_MPI_Comm_split (octree->comm, i_rank, 0, &(octree->rank_comm));

  octree->rank_boxes = PDM_box_set_create (3,
                                           1,
                                           0,
                                           octree->n_used_rank,
                                           gnum_proc,
                                           octree->used_rank_extents,
                                           1,
                                           &n_used_rank,
                                           init_location_proc,
                                           octree->rank_comm);

  int   max_boxes_leaf_shared = 10; // Max number of boxes in a leaf for coarse shared BBTree

  int   max_tree_depth_shared = 6; // Max tree depth for coarse shared BBTree

  float max_box_ratio_shared = 5; // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)
  octree->bt_shared = PDM_box_tree_create (max_tree_depth_shared,
                                           max_boxes_leaf_shared,
                                           max_box_ratio_shared);

  PDM_box_tree_set_boxes (octree->bt_shared,
                          octree->rank_boxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);
  free (gnum_proc);
  free (init_location_proc);
}


/**
 *
 * \brief Assess load imbalance and identify ranks to copy
 *
 * \param [in]   comm                    MPI Communicator
 * \param [in]   f_threshold             Copy threshold (relative to mean nb of request)
 * \param [in]   f_max_copy              Maximum fraction of copied ranks
 * \param [out]  n_copied_ranks          Number of copied ranks
 * \param [out]  copied_ranks            Array of copied ranks
 * \param [out]  n_request_copied_ranks  Initial number or requests of copied ranks
 * \param [out]  mean_n_request          Mean number of requests
 *
 */

static void
_prepare_copies
(
 PDM_MPI_Comm   comm,
 const float    f_threshold,
 const float    f_max_copy,
 int            n_request,
 int           *n_copied_ranks,
 int          **copied_ranks,
 int          **n_request_copied_ranks,
 int           *mean_n_request
 )
{
  *n_copied_ranks = 0;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  int n_max_copy = (int) (f_max_copy * n_rank);
  if (n_max_copy < 1) {
    return;
  }

  int *all_n_request = malloc (sizeof(int) * n_rank);
  PDM_MPI_Allgather (&n_request,    1, PDM_MPI_INT,
                     all_n_request, 1, PDM_MPI_INT,
                     comm);

  // Mean number of requests
  long l_mean_n_request = 0;
  for (int i = 0; i < n_rank; i++) {
    l_mean_n_request += all_n_request[i];
  }
  *mean_n_request = (int) (l_mean_n_request / n_rank);

  float n_threshold = f_threshold * (*mean_n_request);

  // Sort ranks
  int *order = malloc (sizeof(int) * n_rank);
  for (int i = 0; i < n_rank; i++) {
    order[i] = i;
  }

  PDM_sort_int (all_n_request,
                order,
                n_rank);

  if (i_rank == 0) {
    printf("copy threshold = %d (%g), max copy = %d (%g)\n", (int) n_threshold, f_threshold, n_max_copy, f_max_copy);
    printf("avant: min = %d, max = %d (%g times mean), max-min = %d\n",
           all_n_request[0],
           all_n_request[n_rank-1],
           (float) all_n_request[n_rank-1] / (float) (*mean_n_request),
           all_n_request[n_rank-1] - all_n_request[0]);
  }

  // Identify ranks to copy
  *copied_ranks = malloc (sizeof(int) * n_max_copy);
  *n_request_copied_ranks = malloc (sizeof(int) * n_max_copy);
  for (int i = 0; i < n_max_copy; i++) {
    int j = n_rank - i - 1;

    if (all_n_request[j] > n_threshold) {
      (*copied_ranks)[*n_copied_ranks] = order[j];
      (*n_request_copied_ranks)[*n_copied_ranks] = all_n_request[j];
      (*n_copied_ranks)++;
    }
    else {
      break;
    }
  }
  free (all_n_request);
  free (order);

  if (*n_copied_ranks > 0) {
    *copied_ranks = realloc (*copied_ranks, sizeof(int) * (*n_copied_ranks));
    *n_request_copied_ranks = realloc (*n_request_copied_ranks,
                                       sizeof(int) * (*n_copied_ranks));
  }
}







static void _export_octree_points
(
 const char      *filename,
 const _octree_t *octree,
 const int        normalize
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\noctree points\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", octree->n_points);
  if (normalize) {
    double s[3] = {octree->global_extents[0],
                   octree->global_extents[1],
                   octree->global_extents[2]};
    double invd[3] = {1. / (octree->global_extents[3] - octree->global_extents[0]),
                      1. / (octree->global_extents[4] - octree->global_extents[1]),
                      1. / (octree->global_extents[5] - octree->global_extents[2])};
    for (int i = 0; i < octree->n_points; i++) {
      for (int j = 0; j < 3; j++) {
        fprintf(f, "%f ", (octree->points[3*i+j] - s[j]) * invd[j]);
      }
      fprintf(f, "\n");
    }
  } else {
    for (int i = 0; i < octree->n_points; i++) {
      for (int j = 0; j < 3; j++) {
        fprintf(f, "%f ", octree->points[3*i+j]);
      }
      fprintf(f, "\n");
    }
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
  fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
  for (int i = 0; i < octree->n_points; i++) {
    fprintf(f, ""PDM_FMT_G_NUM"\n", octree->points_gnum[i]);
  }

  fclose(f);
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
    if (e[3] < e[0] || e[4] < e[1] || e[5]  < e[2]) {
      fprintf(f, "0. 0. 0.\n0. 0. 0.\n0. 0. 0.\n0. 0. 0.\n0. 0. 0.\n0. 0. 0.\n0. 0. 0.\n0. 0. 0.\n");
    } else {
      fprintf(f, "%f %f %f\n", e[0], e[1], e[2]);
      fprintf(f, "%f %f %f\n", e[3], e[1], e[2]);
      fprintf(f, "%f %f %f\n", e[3], e[4], e[2]);
      fprintf(f, "%f %f %f\n", e[0], e[4], e[2]);
      fprintf(f, "%f %f %f\n", e[0], e[1], e[5]);
      fprintf(f, "%f %f %f\n", e[3], e[1], e[5]);
      fprintf(f, "%f %f %f\n", e[3], e[4], e[5]);
      fprintf(f, "%f %f %f\n", e[0], e[4], e[5]);
    }
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

  if (box_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_box);
    fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < n_box; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", box_g_num[i]);
    }
  }

  fclose(f);
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

  octree->rank_comm = PDM_MPI_COMM_NULL;
  octree->used_rank = NULL;
  octree->used_rank_extents = NULL;
  octree->rank_boxes = NULL;
  octree->bt_shared = NULL;

  octree->n_copied_ranks     = 0;
  octree->copied_ranks       = NULL;
  octree->copied_octants     = NULL;
  octree->n_copied_points    = NULL;
  octree->copied_points      = NULL;
  octree->copied_points_gnum = NULL;
  octree->copied_points_code = NULL;


  octree->shared_codes       = NULL;
  octree->shared_rank_idx    = NULL;
  octree->shared_pts_n       = NULL;
  octree->shared_pts_extents = NULL;


  octree->explicit_nodes_to_build = 1;
  char *env_var = NULL;
  env_var = getenv ("OCTREE_BUILD_EXPLICIT");
  if (env_var != NULL) {
    octree->explicit_nodes_to_build = atoi(env_var);
  }
  octree->explicit_nodes     = NULL;
  octree->n_explicit_nodes   = 0;

  octree->copied_explicit_nodes   = NULL;
  octree->n_copied_explicit_nodes = NULL;

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

  if (octree->used_rank != NULL) {
    free (octree->used_rank);
  }

  if (octree->used_rank_extents != NULL) {
    free (octree->used_rank_extents);
  }

  PDM_box_set_destroy (&octree->rank_boxes);

  if (octree->rank_comm != PDM_MPI_COMM_NULL) {
    PDM_MPI_Comm_free (&(octree->rank_comm));
  }
  PDM_box_tree_destroy (&octree->bt_shared);

  if (octree->shared_codes != NULL) {
    free (octree->shared_codes);
  }

  if (octree->shared_rank_idx != NULL) {
    free (octree->shared_rank_idx);
  }

  if (octree->shared_pts_n != NULL) {
    free (octree->shared_pts_n);
  }

  if (octree->shared_pts_extents != NULL) {
    free (octree->shared_pts_extents);
  }

  if (octree->explicit_nodes != NULL) {
    free (octree->explicit_nodes);
  }

  PDM_para_octree_free_copies (id);

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





static void
_export_nodes
(
 char              *filename,
 int                n_nodes,
 PDM_morton_code_t *nodes,
 double             s[],
 double             d[]
 )
{
  double s0[3] = {0., 0., 0.};
  double d0[3] = {1., 1., 1.};

  double *_s = s0;
  double *_d = d0;

  if (s != NULL) _s = s;
  if (d != NULL) _d = d;


  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nnodes\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_nodes);
  for (int inode = 0; inode < n_nodes; inode++) {
    int l = 1 << nodes[inode].L;
    double side = 1.0 / (double) l;

    for (int k = 0; k < 2; k++) {
      double z = _s[2] + _d[2] * (nodes[inode].X[2] + k) * side;
      for (int j = 0; j < 2; j++) {
        double y = _s[1] + _d[1] * (nodes[inode].X[1] + j) * side;
        for (int i = 0; i < 2; i++) {
          double x = _s[0] + _d[0] * (nodes[inode].X[0] + i) * side;
          fprintf(f, "%f %f %f\n", x, y, z);
        }
      }
    }
  }

  fprintf(f, "CELLS %d %d\n", n_nodes, 9*n_nodes);
  for (int inode = 0; inode < n_nodes; inode++) {
    fprintf(f, "8 ");
    for (int i = 0; i < 8; i++) {
      fprintf(f, "%d ", 8*inode + i);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELL_TYPES %d\n", n_nodes);
  for (int inode = 0; inode < n_nodes; inode++) {
    fprintf(f, "11\n");
  }

  fclose(f);
}


static void
_export_explicit_nodes
(
 const char             *filename,
 const int               n_nodes,
 const _explicit_node_t *nodes
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nexplicit nodes\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    const double *e = nodes[i].pts_extents;
    fprintf(f, "%f %f %f\n", e[0], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[5]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[5]);
  }

  fprintf(f, "CELLS %d %d\n", n_nodes, 9*n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    int j = 8*i;
    fprintf(f, "8 %d %d %d %d %d %d %d %d\n", j, j+1, j+2, j+3, j+4, j+5, j+6, j+7);
  }

  fprintf(f, "CELL_TYPES %d\n", n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    fprintf(f, "12\n");
  }

  fclose(f);
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

    int *weight = PDM_array_const_int(octree->n_points, 1);

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


    //-->>
    if (1) {
      /*
       *  Compress octants
       */
      int n_local_nodes;
      int *send_codes = NULL;
      int *send_n_pts = NULL;
      double *send_extents = NULL;
      _compress_octants (octree->octants,
                         octree->points,
                         &n_local_nodes,
                         &send_codes,
                         &send_n_pts,
                         &send_extents);

      /*// expand extents
      const double eps = 1e-3;
      const double min_delta = 1e-6;
      for (int i = 0; i < n_local_nodes; i++) {
        if (send_n_pts[i] > 0) {
          for (int j = 0; j < 3; j++) {
            double delta = PDM_MAX (min_delta,
                                    eps*(send_extents[6*i+j+3] - send_extents[6*i+j]));
            send_extents[6*i+j]   -= delta;
            send_extents[6*i+j+3] += delta;
          }
        }
        }*/

      int *recv_count = malloc (sizeof(int) * n_ranks);
      PDM_MPI_Allgather (&n_local_nodes, 1, PDM_MPI_INT,
                         recv_count,     1, PDM_MPI_INT,
                         octree->comm);

      octree->shared_rank_idx = PDM_array_new_idx_from_sizes_int (recv_count, n_ranks);
      int n_shared_nodes = octree->shared_rank_idx[n_ranks];
      int *recv_shift = PDM_array_new_idx_from_sizes_int (recv_count, n_ranks);

      octree->shared_pts_n = malloc (sizeof(int) * n_shared_nodes);
      PDM_MPI_Allgatherv (send_n_pts, n_local_nodes, PDM_MPI_INT,
                          octree->shared_pts_n, recv_count, recv_shift, PDM_MPI_INT,
                          octree->comm);
      free (send_n_pts);

      for (int i = 0; i < n_ranks; i++) {
        recv_count[i] *= 4;
        recv_shift[i+1] *= 4;
      }
      int *recv_codes = malloc (sizeof(int) * n_shared_nodes * 4);
      PDM_MPI_Allgatherv (send_codes, 4*n_local_nodes, PDM_MPI_INT,
                          recv_codes, recv_count, recv_shift, PDM_MPI_INT,
                          octree->comm);
      free (send_codes);

      for (int i = 0; i < n_ranks; i++) {
        recv_count[i] /= 4;
        recv_count[i] *= 6;
        recv_shift[i+1] = 6 * octree->shared_rank_idx[i+1];
      }
      octree->shared_pts_extents = malloc (sizeof(double) * n_shared_nodes * 6);
      PDM_MPI_Allgatherv (send_extents, 6*n_local_nodes, PDM_MPI_DOUBLE,
                          octree->shared_pts_extents, recv_count, recv_shift, PDM_MPI_DOUBLE,
                          octree->comm);
      free (send_extents);
      free (recv_count);
      free (recv_shift);

      octree->shared_codes = malloc (sizeof(PDM_morton_code_t) * n_shared_nodes);
      for (int i = 0; i < n_shared_nodes; i++) {
        octree->shared_codes[i].L = (PDM_morton_int_t) recv_codes[4*i];
        for (int j = 0; j < 3; j++) {
          octree->shared_codes[i].X[j] = (PDM_morton_int_t) recv_codes[4*i+j+1];
        }
      }
      free (recv_codes);
    }
    //<<--

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

      PDM_array_reset_int(n_points_children, n_child, 0);

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

      PDM_array_idx_from_sizes_int(n_points_children, n_child - 1, range_children);

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

  /*PDM_timer_hang_on(octree->timer);
  octree->times_elapsed[BUILD_TOTAL] = PDM_timer_elapsed(octree->timer);
  octree->times_cpu[BUILD_TOTAL]     = PDM_timer_cpu(octree->timer);
  octree->times_cpu_u[BUILD_TOTAL]   = PDM_timer_cpu_user(octree->timer);
  octree->times_cpu_s[BUILD_TOTAL]   = PDM_timer_cpu_sys(octree->timer);
  PDM_timer_resume(octree->timer);

  octree->times_elapsed[END] = octree->times_elapsed[BUILD_TOTAL];
  octree->times_cpu[END]     = octree->times_cpu[BUILD_TOTAL];
  octree->times_cpu_u[END]   = octree->times_cpu_u[BUILD_TOTAL];
  octree->times_cpu_s[END]   = octree->times_cpu_s[BUILD_TOTAL];*/

  //-->
  if (rank == 0 && octree->shared_rank_idx != NULL) {
    printf("shared_rank_idx = ");
    for (int i = 0; i <= n_ranks; i++) {
      printf("%d ", octree->shared_rank_idx[i]);
    }
    printf("\n");
  }
  if (0) {//octree->shared_rank_idx != NULL) {
    //const char *pref = "/stck/bandrieu/workspace/paradigma-dev/test/para_octree/shared_octree/";
    const char *pref = "";
    char filename[999];
    sprintf(filename, "%soctree_local_%4.4d.vtk", pref, rank);
    _export_nodes (filename,
                   octree->octants->n_nodes,
                   octree->octants->codes,
                   octree->s,
                   octree->d);

    sprintf(filename, "%soctree_pts_%4.4d.vtk", pref, rank);
    _export_octree_points (filename,
                           octree,
                           0);

    if (rank == 0) {
      for (int i = 0; i < n_ranks; i++) {
        sprintf(filename, "%soctree_shared_%4.4d.vtk", pref, i);
        _export_nodes (filename,
                       octree->shared_rank_idx[i+1] - octree->shared_rank_idx[i],
                       octree->shared_codes + octree->shared_rank_idx[i],
                       octree->s,
                       octree->d);

        printf("rank %d:\n", i);
        for (int j = octree->shared_rank_idx[i]; j < octree->shared_rank_idx[i+1]; j++) {
          double *e = octree->shared_pts_extents + 6*j;
          printf("  %d : L=%d, X=%d %d %d, n_pts = %d, extents = %f %f %f %f %f %f\n",
                 j,
                 octree->shared_codes[j].L,
                 octree->shared_codes[j].X[0],
                 octree->shared_codes[j].X[1],
                 octree->shared_codes[j].X[2],
                 octree->shared_pts_n[j],
                 e[0], e[1], e[2], e[3], e[4], e[5]);
        }

        sprintf(filename, "%soctree_shared_extents_%4.4d.vtk", pref, i);
        _export_boxes (filename,
                       octree->shared_rank_idx[i+1] - octree->shared_rank_idx[i],
                       octree->shared_pts_extents + 6* octree->shared_rank_idx[i],
                       NULL);
      }
    }
  }
  //<<--


  PDM_timer_hang_on(octree->timer);
  b_t_elapsed = PDM_timer_elapsed(octree->timer);
  b_t_cpu     = PDM_timer_cpu(octree->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);
  PDM_timer_resume(octree->timer);

  //-->>
  if (octree->explicit_nodes_to_build) {
    _build_explicit_nodes (octree);
    printf("[%d] %d explicit nodes\n", rank, octree->n_explicit_nodes);

    if (0) {
      for (int i = 0; i < octree->n_explicit_nodes; i++) {
        if (octree->explicit_nodes[i].leaf_id < 0) {
          for (int j = 0; j < 8; j++) {
            int k = octree->explicit_nodes[i].children_id[j];
            if (k < 0) continue;
            /*printf("node %d, n_pts = %d, extents = %f %f %f %f %f %f\n",
                   k,
                   octree->explicit_nodes[k].n_points,
                   octree->explicit_nodes[k].pts_extents[0],
                   octree->explicit_nodes[k].pts_extents[1],
                   octree->explicit_nodes[k].pts_extents[2],
                   octree->explicit_nodes[k].pts_extents[3],
                   octree->explicit_nodes[k].pts_extents[4],
                   octree->explicit_nodes[k].pts_extents[5]);*/
            assert (octree->explicit_nodes[k].ancestor_id == i);
            assert (octree->explicit_nodes[k].n_points > 0);
          }
        }
      }
    }

    if (0) {
      const char *pref = "";
      char filename[999];
      sprintf(filename, "%soctree_explicit_%4.4d.vtk", pref, rank);
      _export_explicit_nodes (filename,
                              octree->n_explicit_nodes,
                              octree->explicit_nodes);
    }
  }
  //<<--
  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_EXPLICIT_NODES] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_EXPLICIT_NODES]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_EXPLICIT_NODES]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_EXPLICIT_NODES]   += e_t_cpu_s - b_t_cpu_s;

  PDM_timer_resume(octree->timer);



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
PDM_para_octree_closest_points
(
 const int    id,
 const int    n_closest_points,
 const int    n_pts,
 double      *pts_coord,
 PDM_g_num_t *pts_g_num,
 PDM_g_num_t *closest_octree_pts_g_num,
 double      *closest_octree_pts_dist2
 )
{
  _octree_t *octree = _get_from_id (id);
  const int dim = octree->dim;

  float f_copy_threshold = 1.15;
  float f_max_copy = 0.1;

  char *env_var = NULL;
  env_var = getenv ("OCTREE_COPY_THRESHOLD");
  if (env_var != NULL) {
    f_copy_threshold = (float) atof(env_var);
  }

  env_var = getenv ("OCTREE_MAX_COPY");
  if (env_var != NULL) {
    f_max_copy = (float) atof(env_var);
  }

  int USE_SHARED_OCTREE = 1;
  env_var = getenv ("USE_SHARED_OCTREE");
  if (env_var != NULL) {
    USE_SHARED_OCTREE = (int) atoi(env_var);
  }


  /* Compute rank extents and build shared bounding-box tree */
  if (octree->used_rank_extents == NULL && !USE_SHARED_OCTREE) {
    _compute_rank_extents (octree);
  }



  int i_rank, n_rank;
  PDM_MPI_Comm_rank (octree->comm, &i_rank);
  PDM_MPI_Comm_size (octree->comm, &n_rank);


  if (i_rank == 0) {
    printf("USE_SHARED_OCTREE = %d\n", USE_SHARED_OCTREE);
  }


  PDM_MPI_Comm    bt_comm;
  PDM_box_set_t  *box_set   = NULL;
  PDM_box_tree_t *bt_shared = NULL;
  int   max_boxes_leaf_shared = 10; // Max number of boxes in a leaf for coarse shared BBTree
  int   max_tree_depth_shared = 6; // Max tree depth for coarse shared BBTree
  float max_box_ratio_shared  = 5; // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)
  if (USE_SHARED_OCTREE && n_rank > 1) {
    assert (octree->shared_rank_idx != NULL);

    PDM_MPI_Comm_split (octree->comm, i_rank, 0, &bt_comm);

    int n_boxes = octree->shared_rank_idx[n_rank];

    const int n_info_location = 3;
    int *init_location_proc = PDM_array_zeros_int (n_info_location * n_boxes);
    PDM_g_num_t *gnum_proc = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * n_boxes);
    for (int i = 0; i < n_boxes; i++) {
      gnum_proc[i] = i + 1;
    }

    box_set = PDM_box_set_create (3,
                                  1,
                                  0,
                                  n_boxes,
                                  gnum_proc,
                                  octree->shared_pts_extents,
                                  1,
                                  &n_boxes,
                                  init_location_proc,
                                  bt_comm);

    bt_shared = PDM_box_tree_create (max_tree_depth_shared,
                                     max_boxes_leaf_shared,
                                     max_box_ratio_shared);

    PDM_box_tree_set_boxes (bt_shared,
                            box_set,
                            PDM_BOX_TREE_ASYNC_LEVEL);

    free (gnum_proc);
    free (init_location_proc);
  }


  /********************************************
   * Distribute target points
   ********************************************/
  int *send_count = NULL;
  int *recv_count = NULL;
  int *send_shift = NULL;
  int *recv_shift = NULL;

  PDM_g_num_t *send_g_num = NULL;
  PDM_g_num_t *recv_g_num = NULL;
  double      *send_coord = NULL;
  double      *recv_coord = NULL;

  int n_recv_pts;
  double d[3], s[3];
  PDM_morton_code_t *pts_code = NULL;


  int n_copied_ranks1 = 0;
  int *copied_ranks1 = NULL;


  int *copied_shift1 = NULL;
  int n_pts_local1 = 0;
  int n_pts_copied1 = 0;
  int n_pts_recv1 = 0;
  int n_pts1;
  int idx_pts1[4];
  idx_pts1[0] = 0;

  PDM_g_num_t *pts_g_num1 = NULL;
  double      *pts_coord1 = NULL;

  if (n_rank > 1) {
    send_count = PDM_array_zeros_int (n_rank);
    int *rank_pt = malloc (sizeof(int) * n_pts);

    if (USE_SHARED_OCTREE) {
      double *node_min_dist = (double *) malloc (sizeof(double) * n_pts);
      PDM_box_tree_min_dist_max_box (bt_shared,
                                     n_pts,
                                     pts_coord,
                                     rank_pt,
                                     node_min_dist);
      /*PDM_box_tree_min_dist_max_box_disjoint (bt_shared,
        n_pts,
        pts_coord,
        rank_pt,
        node_min_dist);*/
      free (node_min_dist);

      for (int i = 0; i < n_pts; i++) {
        int inode = rank_pt[i];

        int l = 0;
        int r = n_rank;
        while (l + 1 < r) {
          int m = l + (r - l)/2;
          if (inode < octree->shared_rank_idx[m])
            r = m;
          else
            l = m;
        }
        int rank = l;
        rank_pt[i] = rank;
        send_count[rank]++;
      }
    }

    else {
      /*   1) Encode the coordinates of every target point */
      pts_code = malloc (sizeof(PDM_morton_code_t) * n_pts);
      _morton_encode_coords (dim,
                             PDM_morton_max_level,
                             octree->global_extents,
                             (size_t) n_pts,
                             pts_coord,
                             pts_code,
                             d,
                             s);

      /*   2) Use binary search to associate each target point to the appropriate process */
      for (int i = 0; i < n_pts; i++) {
        rank_pt[i] = PDM_morton_binary_search (n_rank,
                                               pts_code[i],
                                               octree->rank_octants_index);
        send_count[rank_pt[i]]++;
      }
      free (pts_code);
    }

    /*   3) Exchange send/recv counts */
    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      octree->comm);

    send_shift = malloc (sizeof(int) * (n_rank + 1));
    recv_shift = malloc (sizeof(int) * (n_rank + 1));
    send_shift[0] = 0;
    recv_shift[0] = 0;
    for (int i = 0; i < n_rank; i++) {
      send_shift[i+1] = send_shift[i] + send_count[i];
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }

    n_recv_pts = recv_shift[n_rank];

    int *n_recv_pts_copied_ranks = NULL;
    int mean_n_recv_pts;
    _prepare_copies (octree->comm,
                     f_copy_threshold,
                     f_max_copy,
                     n_recv_pts,
                     &n_copied_ranks1,
                     &copied_ranks1,
                     &n_recv_pts_copied_ranks,
                     &mean_n_recv_pts);
    if (n_recv_pts_copied_ranks != NULL) {
      free (n_recv_pts_copied_ranks);
    }

    if (n_copied_ranks1 > 0) {
      if (i_rank == 0) {
        if (n_copied_ranks1 == 1) {
          printf("phase 1: 1 copied rank: %d\n", copied_ranks1[0]);
        }
        else {
          printf("phase 1: %d copied ranks:", n_copied_ranks1);
          for (int i = 0; i < n_copied_ranks1; i++) {
            printf(" %d", copied_ranks1[i]);
          }
          printf("\n");
        }
      }

      PDM_para_octree_copy_ranks (id,
                                  n_copied_ranks1,
                                  copied_ranks1);
    } else {
      if (i_rank == 0) printf("phase 1: 0 copied ranks\n");
    }

    int *i_copied_rank1 = PDM_array_const_int(n_rank, -1);
    int *copied_count = malloc (sizeof(int) * octree->n_copied_ranks);

    for (int i = 0; i < octree->n_copied_ranks; i++) {
      i_copied_rank1[octree->copied_ranks[i]] = i;
      copied_count[i] = 0;
    }



    n_pts_local1 = 0;
    for (int i = 0; i < n_rank; i++) {
      send_shift[i+1] = send_shift[i];

      if (i == i_rank) {
        n_pts_local1 += send_count[i];
        send_count[i] = 0;
      } else if (i_copied_rank1[i] >= 0) {
        copied_count[i_copied_rank1[i]] = send_count[i];
        send_count[i] = 0;
      } else {
        send_shift[i+1] += send_count[i];
      }
    }

    copied_shift1 = malloc (sizeof(int) * (octree->n_copied_ranks + 1));
    copied_shift1[0] = 0;
    for (int i = 0; i < octree->n_copied_ranks; i++) {
      copied_shift1[i+1] = copied_shift1[i] + copied_count[i];
      copied_count[i] = 0;
    }
    n_pts_copied1 = copied_shift1[octree->n_copied_ranks];

    /*   3bis) Exchange new send/recv counts */
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      octree->comm);
    for (int i = 0; i < n_rank; i++) {
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
      send_count[i] = 0;
    }
    n_pts_recv1 = recv_shift[n_rank];


    idx_pts1[1] = n_pts_local1;
    idx_pts1[2] = idx_pts1[1] + n_pts_recv1;
    idx_pts1[3] = idx_pts1[2] + n_pts_copied1;
    n_pts1 = idx_pts1[3];

    pts_g_num1 = malloc (sizeof(PDM_g_num_t) * n_pts1);
    pts_coord1 = malloc (sizeof(double)      * n_pts1 * dim);

    /*   4) Fill send buffers */
    n_pts_local1 = 0;
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank]*dim);
    recv_g_num = pts_g_num1 + idx_pts1[1];
    recv_coord = pts_coord1 + idx_pts1[1] * dim;

    PDM_g_num_t *local_g_num = pts_g_num1 + idx_pts1[0];
    double      *local_coord = pts_coord1 + idx_pts1[0] * dim;

    PDM_g_num_t *copied_g_num = pts_g_num1 + idx_pts1[2];
    double      *copied_coord = pts_coord1 + idx_pts1[2] * dim;

    for (int i = 0; i < n_pts; i++) {
      int rank = rank_pt[i];

      if (rank == i_rank) {
        local_g_num[n_pts_local1] = pts_g_num[i];
        for (int j = 0; j < dim; j++) {
          local_coord[dim*n_pts_local1 + j] = pts_coord[dim*i + j];
        }
        n_pts_local1++;
      }

      else if (i_copied_rank1[rank] >= 0) {
        rank = i_copied_rank1[rank];
        int k = copied_shift1[rank] + copied_count[rank];
        copied_g_num[k] = pts_g_num[i];
        for (int j = 0; j < dim; j++) {
          copied_coord[dim*k + j] = pts_coord[dim*i + j];
        }
        copied_count[rank]++;
      }

      else {
        int k = send_shift[rank] + send_count[rank];
        send_g_num[k] = pts_g_num[i];
        for (int j = 0; j < dim; j++) {
          send_coord[dim*k + j] = pts_coord[dim*i + j];
        }
        send_count[rank]++;
      }
    }
    if (copied_count != NULL) {
      free (copied_count);
    }
    if (i_copied_rank1 != NULL) {
      free (i_copied_rank1);
    }
    free (rank_pt);

    /*   5) Send gnum buffer */
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       octree->comm);

    /*   6) Send coord buffer */
    for (int i = 0; i < n_rank; i++) {
      send_count[i] *= dim;
      recv_count[i] *= dim;
      send_shift[i+1] *= dim;
      recv_shift[i+1] *= dim;
    }
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       octree->comm);
  }

  /* Single proc */
  else {
    n_recv_pts = n_pts;
    n_pts_local1 = n_pts;
    n_pts_recv1 = 0;
    n_pts_copied1 = 0;

    idx_pts1[1] = n_pts_local1;
    idx_pts1[2] = idx_pts1[1] + n_pts_recv1;
    idx_pts1[3] = idx_pts1[2] + n_pts_copied1;
    n_pts1 = idx_pts1[3];

    pts_g_num1 = pts_g_num;
    pts_coord1 = pts_coord;
  }



  /********************************************
   * First guess : closest source points in the
   * sense of Morton codes
   ********************************************/
  int s_closest = n_closest_points * n_pts1;
  double      *_closest_pts_dist2 = NULL;
  PDM_g_num_t *_closest_pts_g_num = NULL;

  if (n_rank == 1) {
    _closest_pts_g_num = closest_octree_pts_g_num;
    _closest_pts_dist2 = closest_octree_pts_dist2;
  }
  else {
    _closest_pts_dist2 = malloc (sizeof(double)      * s_closest);
    _closest_pts_g_num = malloc (sizeof(PDM_g_num_t) * s_closest);
  }

  /* Encode the coordinates of the received target points */
  pts_code = malloc (sizeof(PDM_morton_code_t) * n_pts1);
  _morton_encode_coords (dim,
                         PDM_morton_max_level,
                         octree->global_extents,
                         (size_t) n_pts1,
                         pts_coord1,
                         pts_code,
                         d,
                         s);

  for (int i = 0; i < s_closest; i++) {
    _closest_pts_dist2[i] = HUGE_VAL;
    _closest_pts_g_num[i] = -1;
  }

  /* First guess in local octree */
  _l_octant_t       *_octants   = octree->octants;
  double            *_src_coord = octree->points;
  PDM_g_num_t       *_src_g_num = octree->points_gnum;
  //PDM_morton_code_t *_src_codes = octree->points_code;

  for (int i = 0; i < idx_pts1[2]; i++) {
    double *last_closest_pt_dist2 = _closest_pts_dist2 + (n_closest_points*(i+1) - 1);
    double *_coord = pts_coord1 + dim * i;

    /* Find base leaf */
    int base = PDM_morton_binary_search (_octants->n_nodes,
                                         pts_code[i],
                                         _octants->codes);

    /* Inspect source points inside base leaf */
    for (int k = 0; k < _octants->n_points[base]; k++) {
      int j = _octants->range[base] + k;
      double dist2 = _pt_to_pt_dist2 (dim,
                                      _coord,
                                      _src_coord + dim*j);
      if (dist2 < *last_closest_pt_dist2) {
        _insertion_sort (dist2,
                         _src_g_num[j],
                         n_closest_points,
                         _closest_pts_dist2 + n_closest_points*i,
                         _closest_pts_g_num + n_closest_points*i);
      }
    }
  }

  /* First guess in copied ranks octree */
  double *_pts_coord1 = pts_coord1 + idx_pts1[2] * dim;
  PDM_morton_code_t *_pts_code = pts_code + idx_pts1[2];
  double      *__closest_pts_dist2 = _closest_pts_dist2 + idx_pts1[2] * n_closest_points;
  PDM_g_num_t *__closest_pts_g_num = _closest_pts_g_num + idx_pts1[2] * n_closest_points;

  for (int rank = 0; rank < octree->n_copied_ranks; rank++) {
    _octants   = octree->copied_octants[rank];
    _src_coord = octree->copied_points[rank];
    _src_g_num = octree->copied_points_gnum[rank];
    //_src_codes = octree->copied_points_code[rank];

    for (int i = copied_shift1[rank]; i < copied_shift1[rank+1]; i++) {
      double *last_closest_pt_dist2 = __closest_pts_dist2 + (n_closest_points*(i+1) - 1);
      double *_coord = _pts_coord1 + dim * i;

      /* Find base leaf */
      int base = PDM_morton_binary_search (_octants->n_nodes,
                                           _pts_code[i],
                                           _octants->codes);

      /* Inspect source points inside base leaf */
      for (int k = 0; k < _octants->n_points[base]; k++) {
        int j = _octants->range[base] + k;
        double dist2 = _pt_to_pt_dist2 (dim,
                                        _coord,
                                        _src_coord + dim*j);
        if (dist2 < *last_closest_pt_dist2) {
          _insertion_sort (dist2,
                           _src_g_num[j],
                           n_closest_points,
                           __closest_pts_dist2 + n_closest_points*i,
                           __closest_pts_g_num + n_closest_points*i);
        }
      }
    }
  } // End of loops on copied ranks
  free (pts_code);


  /*
   *  Search closest points
   */
  if (octree->explicit_nodes_to_build) {
    _closest_points_explicit (n_closest_points,
                              octree,
                              -1,
                              idx_pts1[2],
                              pts_coord1,
                              _closest_pts_g_num,
                              _closest_pts_dist2);
  }
  else {
  _closest_points (n_closest_points,
                   dim,
                   octree->d,
                   octree->s,
                   octree->octants,
                   octree->points_gnum,
                   octree->points,
                   idx_pts1[2],
                   pts_coord1,
                   _closest_pts_g_num,
                   _closest_pts_dist2);
  }

  for (int i = 0; i < octree->n_copied_ranks; i++) {
    if (octree->explicit_nodes_to_build) {
      _closest_points_explicit (n_closest_points,
                                octree,
                                i,
                                copied_shift1[i+1] - copied_shift1[i],
                                _pts_coord1 + copied_shift1[i]*dim,
                                __closest_pts_g_num + n_closest_points*copied_shift1[i],
                                __closest_pts_dist2 + n_closest_points*copied_shift1[i]);
    }
    else {
      _closest_points (n_closest_points,
                       dim,
                       octree->d,
                       octree->s,
                       octree->copied_octants[i],
                       octree->copied_points_gnum[i],
                       octree->copied_points[i],
                       copied_shift1[i+1] - copied_shift1[i],
                       _pts_coord1 + copied_shift1[i]*dim,
                       __closest_pts_g_num + n_closest_points*copied_shift1[i],
                       __closest_pts_dist2 + n_closest_points*copied_shift1[i]);
    }
  }


  if (n_rank == 1) {
    return;
  }

  /*
   *  Fill block data
   */
  PDM_g_num_t *block_distrib_idx =
    PDM_compute_uniform_entity_distribution_from_partition (octree->comm,
                                                            1,
                                                            &n_pts,
                                                            (const PDM_g_num_t **) &pts_g_num);

  PDM_part_to_block_t *ptb1 = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                         PDM_PART_TO_BLOCK_POST_MERGE,
                                                         1.,
                                                         &pts_g_num1,
                                                         block_distrib_idx,
                                                         &n_pts1,
                                                         1,
                                                         octree->comm);

  int *part_stride = PDM_array_const_int(n_pts1, n_closest_points);

  int *block_stride = NULL;
  double *block_closest_pts_dist2 = NULL;
  PDM_part_to_block_exch (ptb1,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &_closest_pts_dist2,
                          &block_stride,
                          (void **) &block_closest_pts_dist2);
  free (block_stride);

  PDM_g_num_t *block_closest_pts_g_num = NULL;
  PDM_part_to_block_exch (ptb1,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &_closest_pts_g_num,
                          &block_stride,
                          (void **) &block_closest_pts_g_num);
  free (block_stride);
  free (part_stride);

  ptb1 = PDM_part_to_block_free (ptb1);



  /*
   *  Find ranks that may contain closer points
   */
  int *close_ranks_idx = NULL;
  int *close_ranks = NULL;

  double *upper_bound_dist2 = malloc (sizeof(double) * n_pts1);
  for (int i = 0; i < n_pts1; i++) {
    upper_bound_dist2[i] = _closest_pts_dist2[n_closest_points*(i+1) - 1];
  }

  if (USE_SHARED_OCTREE) {
    assert (octree->shared_rank_idx != NULL);

    int *close_nodes_idx = NULL;
    int *close_nodes     = NULL;
    PDM_box_tree_closest_upper_bound_dist_boxes_get (bt_shared,
                                                     n_pts1,
                                                     pts_coord1,
                                                     upper_bound_dist2,
                                                     &close_nodes_idx,
                                                     &close_nodes);
    int *tag_rank = PDM_array_zeros_int (n_rank);

    int tmp_size = 4 * n_pts1;
    close_ranks = malloc (sizeof(int) * tmp_size);
    close_ranks_idx = malloc (sizeof(int) * (n_pts1 + 1));
    close_ranks_idx[0] = 0;

    for (int i = 0; i < n_pts1; i++) {
      close_ranks_idx[i+1] = close_ranks_idx[i];

      for (int j = close_nodes_idx[i]; j < close_nodes_idx[i+1]; j++) {
        int inode = close_nodes[j];
        if (octree->shared_pts_n[inode] == 0) continue;

        int l = 0;
        int r = n_rank;
        while (l + 1 < r) {
          int m = l + (r - l)/2;
          if (inode < octree->shared_rank_idx[m])
            r = m;
          else
            l = m;
        }
        int rank = l;

        if (tag_rank[rank]) continue;

        if (tmp_size <= close_ranks_idx[i+1]) {
          tmp_size *= 2;
          close_ranks = realloc (close_ranks, sizeof(int) * tmp_size);
        }

        close_ranks[close_ranks_idx[i+1]++] = rank;
        tag_rank[rank] = 1;
      }

      for (int j = close_ranks_idx[i]; j < close_ranks_idx[i+1]; j++) {
        tag_rank[close_ranks[j]] = 0;
      }
    }

    free (close_nodes_idx);
    free (close_nodes);
    free (tag_rank);

  }

  else {
    PDM_box_tree_closest_upper_bound_dist_boxes_get (octree->bt_shared,
                                                     n_pts1,
                                                     pts_coord1,
                                                     upper_bound_dist2,
                                                     &close_ranks_idx,
                                                     &close_ranks);
    if (octree->n_used_rank < n_rank) {
      for (int i = 0; i < close_ranks_idx[n_pts1]; i++) {
        int j = close_ranks[i];
        close_ranks[i] = octree->used_rank[j];
      }
    }
  }

  PDM_array_reset_int(send_count, n_rank, 0);

  for (int i = 0; i < idx_pts1[2]; i++) {
    for (int j = close_ranks_idx[i]; j < close_ranks_idx[i+1]; j++) {
      int rank = close_ranks[j];
      if (rank != i_rank) {
        send_count[rank]++;
      }
    }
  }

  int *_close_ranks_idx = close_ranks_idx + idx_pts1[2];
  for (int c = 0; c < octree->n_copied_ranks; c++) {
    int _rank = octree->copied_ranks[c];

    for (int i = copied_shift1[c]; i < copied_shift1[c+1]; i++) {
      for (int j = _close_ranks_idx[i]; j < _close_ranks_idx[i+1]; j++) {
        int rank = close_ranks[j];
        if (rank != _rank) {
          send_count[rank]++;
        }
      }
    }
  }

  PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                    recv_count, 1, PDM_MPI_INT,
                    octree->comm);

  for (int i = 0; i < n_rank; i++) {
    send_shift[i+1] = send_shift[i] + send_count[i];
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
  }

  n_recv_pts = recv_shift[n_rank];

  PDM_para_octree_free_copies (id);

  int n_copied_ranks2;
  int *copied_ranks2 = NULL;
  int *n_recv_pts_copied_ranks = NULL;
  int mean_n_recv_pts;
  _prepare_copies (octree->comm,
                   f_copy_threshold,
                   f_max_copy,
                   n_recv_pts,
                   &n_copied_ranks2,
                   &copied_ranks2,
                   &n_recv_pts_copied_ranks,
                   &mean_n_recv_pts);

  if (n_copied_ranks2 > 0) {
    if (i_rank == 0) {
      if (n_copied_ranks2 == 1) {
        printf("phase 2: 1 copied rank: %d\n", copied_ranks2[0]);
      }
      else {
        printf("phase 2: %d copied ranks:", n_copied_ranks2);
        for (int i = 0; i < n_copied_ranks2; i++) {
          printf(" %d", copied_ranks2[i]);
        }
        printf("\n");
      }
    }

    PDM_para_octree_copy_ranks (id,
                                n_copied_ranks2,
                                copied_ranks2);
  } else {
    if (i_rank == 0) printf("phase 2: 0 copied ranks\n");
  }

  int *i_copied_rank2 = PDM_array_const_int(n_rank,-1);
  int *copied_count = malloc (sizeof(int) * octree->n_copied_ranks);

  for (int i = 0; i < octree->n_copied_ranks; i++) {
    i_copied_rank2[octree->copied_ranks[i]] = i;
    copied_count[i] = 0;
  }


  int n_pts_local2 = 0;
  int n_pts_recv2 = 0;
  int n_pts_copied2 = 0;
  for (int i = 0; i < octree->n_copied_ranks; i++) {
    int rank = octree->copied_ranks[i];
    if (rank != i_rank) {
      int si = send_count[rank];

      si = PDM_MIN (si, PDM_MAX (0, (n_recv_pts_copied_ranks[i] - n_recv_pts)/2));
      if (i_copied_rank2[i_rank] < 0) {
        si = PDM_MIN (si, PDM_MAX (0, mean_n_recv_pts - n_recv_pts));
      }

      copied_count[i] = si;
      n_recv_pts += si;
    }
  }

  for (int i = 0; i < n_rank; i++) {
    send_shift[i+1] = send_shift[i];

    if (i == i_rank) {
      n_pts_local2 += send_count[i];
      send_count[i] = 0;
    }
    else if (i_copied_rank2[i] >= 0) {
      send_count[i] -= copied_count[i_copied_rank2[i]];
    }

    send_shift[i+1] += send_count[i];
  }

  if (n_recv_pts_copied_ranks != NULL) {
    free (n_recv_pts_copied_ranks);
  }

  int *copied_shift2 = malloc (sizeof(int) * (octree->n_copied_ranks + 1));
  int *copied_count_tmp = malloc (sizeof(int) * octree->n_copied_ranks);
  copied_shift2[0] = 0;
  for (int i = 0; i < octree->n_copied_ranks; i++) {
    copied_shift2[i+1] = copied_shift2[i] + copied_count[i];
    copied_count_tmp[i] = 0;
  }
  n_pts_copied2 = copied_shift2[octree->n_copied_ranks];


  /* Exchange new send/recv counts */
  PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                    recv_count, 1, PDM_MPI_INT,
                    octree->comm);

  for (int i = 0; i < n_rank; i++) {
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
    send_count[i] = 0;
  }
  n_pts_recv2 = recv_shift[n_rank];

  int idx_pts2[4];
  idx_pts2[0] = 0;
  idx_pts2[1] = n_pts_local2;
  idx_pts2[2] = idx_pts2[1] + n_pts_recv2;
  idx_pts2[3] = idx_pts2[2] + n_pts_copied2;
  int n_pts2 = idx_pts2[3];

  PDM_g_num_t *pts_g_num2 = malloc (sizeof(PDM_g_num_t) * n_pts2);
  double      *pts_coord2 = malloc (sizeof(double)      * n_pts2 * dim);
  double      *_closest_pts_dist22 = malloc (sizeof(double) * n_pts2 * n_closest_points);

  send_g_num = realloc (send_g_num, sizeof(PDM_g_num_t) * send_shift[n_rank]);
  int s_data = dim + 1;
  send_coord = realloc (send_coord, sizeof(double) * send_shift[n_rank] * s_data);

  recv_g_num = pts_g_num2 + idx_pts2[1];
  recv_coord = pts_coord2 + idx_pts2[1] * dim;

  PDM_g_num_t *local_g_num = pts_g_num2 + idx_pts2[0];
  double      *local_coord = pts_coord2 + idx_pts2[0] * dim;
  double      *local_closest_pts_dist2 = _closest_pts_dist22 + idx_pts2[0] * n_closest_points;

  PDM_g_num_t *copied_g_num = pts_g_num2 + idx_pts2[2];
  double      *copied_coord = pts_coord2 + idx_pts2[2] * dim;
  double      *copied_closest_pts_dist2 = _closest_pts_dist22 + idx_pts2[2] * n_closest_points;

  n_pts_local2 = 0;
  for (int i = 0; i < idx_pts1[2]; i++) {
    for (int j = close_ranks_idx[i]; j < close_ranks_idx[i+1]; j++) {
      int rank = close_ranks[j];
      if (rank != i_rank) {

        if (i_copied_rank2[rank] >= 0) {
          int rank2 = i_copied_rank2[rank];

          if (copied_count_tmp[rank2] < copied_count[rank2]) {
            int k = copied_shift2[rank2] + copied_count_tmp[rank2];
            copied_g_num[k] = pts_g_num1[i];
            for (int l = 0; l < dim; l++) {
              copied_coord[dim*k+l] = pts_coord1[dim*i+l];
            }
            for (int l = 0; l < n_closest_points; l++) {
              copied_closest_pts_dist2[k*n_closest_points + l] = upper_bound_dist2[i];
            }
            copied_count_tmp[rank2]++;
          }
          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = pts_g_num1[i];
            for (int l = 0; l < dim; l++) {
              send_coord[s_data*k+l] = pts_coord1[dim*i+l];
            }
            send_coord[s_data*k+dim] = upper_bound_dist2[i];
            send_count[rank]++;
          }
        }

        else {
          int k = send_shift[rank] + send_count[rank];
          send_g_num[k] = pts_g_num1[i];
          for (int l = 0; l < dim; l++) {
            send_coord[s_data*k+l] = pts_coord1[dim*i+l];
          }
          send_coord[s_data*k+dim] = upper_bound_dist2[i];
          send_count[rank]++;
        }

      }
    }
  }

  PDM_g_num_t *_pts_g_num1 = pts_g_num1 + idx_pts1[2];
  double *_upper_bound_dist2 = upper_bound_dist2 + idx_pts1[2];
  for (int c = 0; c < n_copied_ranks1; c++) {
    int rank1 = copied_ranks1[c];

    for (int i = copied_shift1[c]; i < copied_shift1[c+1]; i++) {
      for (int j = _close_ranks_idx[i]; j < _close_ranks_idx[i+1]; j++) {
        int rank = close_ranks[j];
        if (rank != rank1) {

          if (rank == i_rank) {
            local_g_num[n_pts_local2] = _pts_g_num1[i];
            for (int l = 0; l < dim; l++) {
              local_coord[dim*n_pts_local2 + l] = _pts_coord1[dim*i + l];
            }
            for (int l = 0; l < n_closest_points; l++) {
              local_closest_pts_dist2[n_pts_local2*n_closest_points + l] = _upper_bound_dist2[i];
            }
            n_pts_local2++;
          }

          else if (i_copied_rank2[rank] >= 0) {
            int rank2 = i_copied_rank2[rank];

            if (copied_count_tmp[rank2] < copied_count[rank2]) {
              int k = copied_shift2[rank2] + copied_count_tmp[rank2];
              copied_g_num[k] = _pts_g_num1[i];
              for (int l = 0; l < dim; l++) {
                copied_coord[dim*k+l] = _pts_coord1[dim*i+l];
              }
              for (int l = 0; l < n_closest_points; l++) {
                copied_closest_pts_dist2[k*n_closest_points + l] = _upper_bound_dist2[i];
              }
              copied_count_tmp[rank2]++;
            }
            else {
              int k = send_shift[rank] + send_count[rank];
              send_g_num[k] = _pts_g_num1[i];
              for (int l = 0; l < dim; l++) {
                send_coord[s_data*k+l] = _pts_coord1[dim*i+l];
              }
              send_coord[s_data*k+dim] = _upper_bound_dist2[i];
              send_count[rank]++;
            }
          }

          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = _pts_g_num1[i];
            for (int l = 0; l < dim; l++) {
              send_coord[s_data*k+l] = _pts_coord1[dim*i+l];
            }
            send_coord[s_data*k+dim] = _upper_bound_dist2[i];
            send_count[rank]++;
          }
        }
      }
    }
  }
  free (upper_bound_dist2);


  if (copied_shift1 != NULL) {
    free (copied_shift1);
  }
  if (copied_count != NULL) {
    free (copied_count);
  }
  if (copied_count_tmp != NULL) {
    free (copied_count_tmp);
  }
  if (i_copied_rank2 != NULL) {
    free (i_copied_rank2);
  }
  free (close_ranks_idx);
  free (close_ranks);
  free (pts_coord1);
  free (pts_g_num1);
  free (_closest_pts_dist2);

  PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                     recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                     octree->comm);
  free (send_g_num);

  //printf ("[%4d] phase 2: n_recv_pts = %8d (wihtout copies: %8d)\n", i_rank, n_pts2, n_recv_pts_no_copies);

  for (int i = 0; i < n_rank; i++) {
    send_count[i] *= s_data;
    recv_count[i] *= s_data;
    send_shift[i+1] *= s_data;
    recv_shift[i+1] *= s_data;
  }
  double *recv_double = malloc (sizeof(double) * recv_shift[n_rank]);
  PDM_MPI_Alltoallv (send_coord,  send_count, send_shift, PDM_MPI_DOUBLE,
                     recv_double, recv_count, recv_shift, PDM_MPI_DOUBLE,
                     octree->comm);
  free (send_coord);

  double *recv_closest_pt_dist2 = _closest_pts_dist22 + idx_pts2[1] * n_closest_points;
  int idx1 = 0, idx2 = 0;
  for (int i = 0; i < n_pts_recv2; i++) {
    for (int j = 0; j < dim; j++) {
      recv_coord[idx1++] = recv_double[idx2++];
    }
    for (int j = 0; j < n_closest_points; j++) {
      recv_closest_pt_dist2[i*n_closest_points + j] = recv_double[idx2];
    }
    idx2++;
  }
  free (recv_double);
  free (send_count);
  free (send_shift);
  free (recv_count);
  free (recv_shift);


  _closest_pts_g_num = realloc (_closest_pts_g_num, sizeof(PDM_g_num_t) * n_pts2 * n_closest_points);
  PDM_array_reset_gnum(_closest_pts_g_num, n_pts2 * n_closest_points, -1);

  if (octree->explicit_nodes_to_build) {
    _closest_points_explicit (n_closest_points,
                              octree,
                              -1,
                              idx_pts2[2],
                              pts_coord2,
                              _closest_pts_g_num,
                              _closest_pts_dist22);
  }
  else {
      _closest_points (n_closest_points,
                       dim,
                       octree->d,
                       octree->s,
                       octree->octants,
                       octree->points_gnum,
                       octree->points,
                       idx_pts2[2],
                       pts_coord2,
                       _closest_pts_g_num,
                       _closest_pts_dist22);
  }

  double *_pts_coord2 = pts_coord2 + idx_pts2[2] * dim;
  __closest_pts_dist2 = _closest_pts_dist22 + idx_pts2[2] * n_closest_points;
  __closest_pts_g_num = _closest_pts_g_num  + idx_pts2[2] * n_closest_points;
  for (int i = 0; i < octree->n_copied_ranks; i++) {
    if (octree->explicit_nodes_to_build) {
      _closest_points_explicit (n_closest_points,
                                octree,
                                i,
                                copied_shift2[i+1] - copied_shift2[i],
                                _pts_coord2 + copied_shift2[i]*dim,
                                __closest_pts_g_num + copied_shift2[i] * n_closest_points,
                                __closest_pts_dist2 + copied_shift2[i] * n_closest_points);
    }
    else {
      _closest_points (n_closest_points,
                       dim,
                       octree->d,
                       octree->s,
                       octree->copied_octants[i],
                       octree->copied_points_gnum[i],
                       octree->copied_points[i],
                       copied_shift2[i+1] - copied_shift2[i],
                       _pts_coord2 + copied_shift2[i]*dim,
                       __closest_pts_g_num + copied_shift2[i] * n_closest_points,
                       __closest_pts_dist2 + copied_shift2[i] * n_closest_points);
    }
  }
  free (pts_coord2);

  if (copied_ranks2 != NULL) {
    free (copied_ranks2);
  }
  free (copied_shift2);


  /*
   *  Back to original partitioning
   */
  /* 1) Part-to-block */
  ptb1 = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                    PDM_PART_TO_BLOCK_POST_MERGE,
                                    1.,
                                    &pts_g_num2,
                                    block_distrib_idx,
                                    &n_pts2,
                                    1,
                                    octree->comm);
  part_stride = PDM_array_const_int(n_pts2, n_closest_points);

  double *tmp_block_closest_pts_dist2 = NULL;
  PDM_part_to_block_exch (ptb1,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &_closest_pts_dist22,
                          &block_stride,
                          (void **) &tmp_block_closest_pts_dist2);
  free (block_stride);
  free (_closest_pts_dist22);

  PDM_g_num_t *tmp_block_closest_pts_g_num = NULL;
  PDM_part_to_block_exch (ptb1,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &_closest_pts_g_num,
                          &block_stride,
                          (void **) &tmp_block_closest_pts_g_num);
  free (_closest_pts_g_num);
  free (part_stride);

  /* Merge block data */
  PDM_g_num_t *block_g_num1 = PDM_part_to_block_block_gnum_get (ptb1);
  const int n_pts_block1 = PDM_part_to_block_n_elt_block_get (ptb1);

  int idx = 0;
  for (int i1 = 0; i1 < n_pts_block1; i1++) {
    int i = (int) (block_g_num1[i1] - block_distrib_idx[i_rank] - 1);
    double *last_closest_pt_dist2 = block_closest_pts_dist2 + (n_closest_points*(i+1) - 1);

    // use more efficient method?
    for (int j = 0; j < block_stride[i1]; j++) {
      if (tmp_block_closest_pts_dist2[idx] < *last_closest_pt_dist2) {
        assert (tmp_block_closest_pts_g_num[idx] > 0);//
        _insertion_sort (tmp_block_closest_pts_dist2[idx],
                         tmp_block_closest_pts_g_num[idx],
                         n_closest_points,
                         block_closest_pts_dist2 + n_closest_points*i,
                         block_closest_pts_g_num + n_closest_points*i);
      }
      idx++;
    }
  }
  free (block_stride);
  free (tmp_block_closest_pts_dist2);
  free (tmp_block_closest_pts_g_num);


  /* 2) Block-to-part */
  PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                       (const PDM_g_num_t **) &pts_g_num,
                                                       &n_pts,
                                                       1,
                                                       octree->comm);
  int stride = n_closest_points;
  PDM_block_to_part_exch (btp,
                          sizeof(double),
                          PDM_STRIDE_CST,
                          &stride,
                          block_closest_pts_dist2,
                          NULL,
                          (void **) &closest_octree_pts_dist2);
  free (block_closest_pts_dist2);

  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          &stride,
                          block_closest_pts_g_num,
                          NULL,
                          (void **) &closest_octree_pts_g_num);
  free (block_closest_pts_g_num);
  free (pts_g_num2);

  ptb1 = PDM_part_to_block_free (ptb1);
  btp = PDM_block_to_part_free (btp);
  free (block_distrib_idx);

  if (copied_ranks1 != NULL) {
    free (copied_ranks1);
  }

  if (USE_SHARED_OCTREE) {
    PDM_box_set_destroy (&box_set);
    PDM_MPI_Comm_free (&bt_comm);
    PDM_box_tree_destroy (&bt_shared);
  }
}


/**
 *
 * Look for single closest point stored inside an octree
 *
 * \param [in]   id                     Identifier
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point coordinates
 * \param [in]   pts_g_num              Point global numbers
 * \param [out]  closest_octree_pt_id   Closest point in octree global number
 * \param [out]  closest_octree_pt_dist Closest point in octree distance
 *
 */

void
PDM_para_octree_single_closest_point
(
 const int    id,
 const int    n_pts,
 double      *pts_coord,
 PDM_g_num_t *pts_g_num,
 PDM_g_num_t *closest_octree_pt_g_num,
 double      *closest_octree_pt_dist2
 )
{
  _octree_t *octree = _get_from_id (id);
  const int dim = octree->dim;

  float f_copy_threshold = 1.1;
  float f_max_copy = 0.1;

  char *env_var = NULL;
  env_var = getenv ("OCTREE_COPY_THRESHOLD");
  if (env_var != NULL) {
    f_copy_threshold = (float) atof(env_var);
  }

  env_var = getenv ("OCTREE_MAX_COPY");
  if (env_var != NULL) {
    f_max_copy = (float) atof(env_var);
  }

  int USE_SHARED_OCTREE = 1;
  env_var = getenv ("USE_SHARED_OCTREE");
  if (env_var != NULL) {
    USE_SHARED_OCTREE = (int) atoi(env_var);
  }

  int DETAIL_TIMER = 0;
  env_var = getenv ("DETAIL_TIMER");
  if (env_var != NULL) {
    DETAIL_TIMER = (int) atoi(env_var);
  }

  const int ntimer=15;
  PDM_timer_t *timer; /*!< Timer */
  if (DETAIL_TIMER) {
    timer = PDM_timer_create ();
    PDM_timer_init (timer);
  }
  double times_elapsed[ntimer]; /*!< Elapsed time */

  double times_cpu[ntimer];     /*!< CPU time */

  double times_cpu_u[ntimer];  /*!< User CPU time */

  double times_cpu_s[ntimer];  /*!< System CPU time */

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  if (DETAIL_TIMER) {
    for (int i = 0; i < ntimer; i++) {
      times_elapsed[i] = 0;
      times_cpu[i]     = 0;
      times_cpu_u[i]   = 0;
      times_cpu_s[i]   = 0;
    }

    //PDM_timer_hang_on(timer);
    times_elapsed[0] = PDM_timer_elapsed(timer);
    times_cpu[0]     = PDM_timer_cpu(timer);
    times_cpu_u[0]   = PDM_timer_cpu_user(timer);
    times_cpu_s[0]   = PDM_timer_cpu_sys(timer);

    b_t_elapsed = times_elapsed[0];
    b_t_cpu     = times_cpu[0];
    b_t_cpu_u   = times_cpu_u[0];
    b_t_cpu_s   = times_cpu_s[0];

    PDM_timer_resume(timer);
  }


  /* Compute rank extents and build shared bounding-box tree */
  if (octree->used_rank_extents == NULL && !USE_SHARED_OCTREE) {
    _compute_rank_extents (octree);
  }


  int i_rank, n_rank;
  PDM_MPI_Comm_rank (octree->comm, &i_rank);
  PDM_MPI_Comm_size (octree->comm, &n_rank);

  if (i_rank == 0) {
    printf("USE_SHARED_OCTREE = %d\n", USE_SHARED_OCTREE);
  }


  PDM_MPI_Comm    bt_comm;
  PDM_box_set_t  *box_set   = NULL;
  PDM_box_tree_t *bt_shared = NULL;
  int   max_boxes_leaf_shared = 10; // Max number of boxes in a leaf for coarse shared BBTree
  int   max_tree_depth_shared = 6; // Max tree depth for coarse shared BBTree
  float max_box_ratio_shared  = 5; // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)
  if (USE_SHARED_OCTREE && n_rank > 1) {
    assert (octree->shared_rank_idx != NULL);

    PDM_MPI_Comm_split (octree->comm, i_rank, 0, &bt_comm);

    int n_boxes = octree->shared_rank_idx[n_rank];

    const int n_info_location = 3;
    int *init_location_proc = PDM_array_zeros_int (n_info_location * n_boxes);
    PDM_g_num_t *gnum_proc = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * n_boxes);
    for (int i = 0; i < n_boxes; i++) {
      gnum_proc[i] = i + 1;
    }

    box_set = PDM_box_set_create (3,
                                  1,
                                  0,
                                  n_boxes,
                                  gnum_proc,
                                  octree->shared_pts_extents,
                                  1,
                                  &n_boxes,
                                  init_location_proc,
                                  bt_comm);

    bt_shared = PDM_box_tree_create (max_tree_depth_shared,
                                     max_boxes_leaf_shared,
                                     max_box_ratio_shared);

    PDM_box_tree_set_boxes (bt_shared,
                            box_set,
                            PDM_BOX_TREE_ASYNC_LEVEL);

    free (gnum_proc);
    free (init_location_proc);
  }


  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[1] += e_t_elapsed - b_t_elapsed;
    times_cpu[1]     += e_t_cpu - b_t_cpu;
    times_cpu_u[1]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[1]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }


  /********************************************
   * Distribute target points
   ********************************************/
  int *send_count = NULL;
  int *recv_count = NULL;
  int *send_shift = NULL;
  int *recv_shift = NULL;

  PDM_g_num_t *send_g_num = NULL;
  PDM_g_num_t *recv_g_num = NULL;
  double      *send_coord = NULL;
  double      *recv_coord = NULL;

  int n_recv_pts;
  double d[3], s[3];
  PDM_morton_code_t *pts_code = NULL;


  int n_copied_ranks1 = 0;
  int *copied_ranks1 = NULL;


  int *copied_shift1 = NULL;
  int n_pts_local1 = 0;
  int n_pts_copied1 = 0;
  int n_pts_recv1 = 0;
  int n_pts1;
  int idx_pts1[4];
  idx_pts1[0] = 0;

  PDM_g_num_t *pts_g_num1 = NULL;
  double      *pts_coord1 = NULL;

  if (n_rank > 1) {
    send_count = PDM_array_zeros_int (n_rank);
    int *rank_pt = malloc (sizeof(int) * n_pts);

    if (USE_SHARED_OCTREE) {
      double *node_min_dist = (double *) malloc (sizeof(double) * n_pts);
      PDM_box_tree_min_dist_max_box (bt_shared,
                                     n_pts,
                                     pts_coord,
                                     rank_pt,
                                     node_min_dist);
      /*PDM_box_tree_min_dist_max_box_disjoint (bt_shared,
                                              n_pts,
                                              pts_coord,
                                              rank_pt,
                                              node_min_dist);*/
      free (node_min_dist);

      for (int i = 0; i < n_pts; i++) {
        int inode = rank_pt[i];

        int l = 0;
        int r = n_rank;
        while (l + 1 < r) {
          int m = l + (r - l)/2;
          if (inode < octree->shared_rank_idx[m])
            r = m;
          else
            l = m;
        }
        int rank = l;
        rank_pt[i] = rank;
        send_count[rank]++;
      }
    }

    else {
      /*   1) Encode the coordinates of every target point */
      pts_code = malloc (sizeof(PDM_morton_code_t) * n_pts);
      _morton_encode_coords (dim,
                             PDM_morton_max_level,
                             octree->global_extents,
                             (size_t) n_pts,
                             pts_coord,
                             pts_code,
                             d,
                             s);

      /*   2) Use binary search to associate each target point to the appropriate process */
      for (int i = 0; i < n_pts; i++) {
        rank_pt[i] = PDM_morton_binary_search (n_rank,
                                               pts_code[i],
                                               octree->rank_octants_index);
        send_count[rank_pt[i]]++;
      }
      free (pts_code);
    }

    /*   3) Exchange send/recv counts */
    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      octree->comm);

    send_shift = malloc (sizeof(int) * (n_rank + 1));
    recv_shift = malloc (sizeof(int) * (n_rank + 1));
    send_shift[0] = 0;
    recv_shift[0] = 0;
    for (int i = 0; i < n_rank; i++) {
      send_shift[i+1] = send_shift[i] + send_count[i];
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }

    n_recv_pts = recv_shift[n_rank];

    int *n_recv_pts_copied_ranks = NULL;
    int mean_n_recv_pts;
    _prepare_copies (octree->comm,
                     f_copy_threshold,
                     f_max_copy,
                     n_recv_pts,
                     &n_copied_ranks1,
                     &copied_ranks1,
                     &n_recv_pts_copied_ranks,
                     &mean_n_recv_pts);
    if (n_recv_pts_copied_ranks != NULL) {
      free (n_recv_pts_copied_ranks);
    }

    if (n_copied_ranks1 > 0) {
      if (i_rank == 0) {
        if (n_copied_ranks1 == 1) {
          printf("phase 1: 1 copied rank: %d\n", copied_ranks1[0]);
        }
        else {
          printf("phase 1: %d copied ranks:", n_copied_ranks1);
          for (int i = 0; i < n_copied_ranks1; i++) {
            printf(" %d", copied_ranks1[i]);
          }
          printf("\n");
        }
      }

      PDM_para_octree_copy_ranks (id,
                                  n_copied_ranks1,
                                  copied_ranks1);
    } else {
      if (i_rank == 0) printf("phase 1: 0 copied ranks\n");
    }

    int *i_copied_rank1 = PDM_array_const_int(n_rank, -1);
    int *copied_count = malloc (sizeof(int) * octree->n_copied_ranks);

    for (int i = 0; i < octree->n_copied_ranks; i++) {
      i_copied_rank1[octree->copied_ranks[i]] = i;
      copied_count[i] = 0;
    }



    n_pts_local1 = 0;
    for (int i = 0; i < n_rank; i++) {
      send_shift[i+1] = send_shift[i];

      if (i == i_rank) {
        n_pts_local1 += send_count[i];
        send_count[i] = 0;
      } else if (i_copied_rank1[i] >= 0) {
        copied_count[i_copied_rank1[i]] = send_count[i];
        send_count[i] = 0;
      } else {
        send_shift[i+1] += send_count[i];
      }
    }

    copied_shift1 = malloc (sizeof(int) * (octree->n_copied_ranks + 1));
    copied_shift1[0] = 0;
    for (int i = 0; i < octree->n_copied_ranks; i++) {
      copied_shift1[i+1] = copied_shift1[i] + copied_count[i];
      copied_count[i] = 0;
    }
    n_pts_copied1 = copied_shift1[octree->n_copied_ranks];

    /*   3bis) Exchange new send/recv counts */
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      octree->comm);
    for (int i = 0; i < n_rank; i++) {
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
      send_count[i] = 0;
    }
    n_pts_recv1 = recv_shift[n_rank];


    idx_pts1[1] = n_pts_local1;
    idx_pts1[2] = idx_pts1[1] + n_pts_recv1;
    idx_pts1[3] = idx_pts1[2] + n_pts_copied1;
    n_pts1 = idx_pts1[3];

    pts_g_num1 = malloc (sizeof(PDM_g_num_t) * n_pts1);
    pts_coord1 = malloc (sizeof(double)      * n_pts1 * dim);

    /*   4) Fill send buffers */
    n_pts_local1 = 0;
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank]*dim);
    recv_g_num = pts_g_num1 + idx_pts1[1];
    recv_coord = pts_coord1 + idx_pts1[1] * dim;

    PDM_g_num_t *local_g_num = pts_g_num1 + idx_pts1[0];
    double      *local_coord = pts_coord1 + idx_pts1[0] * dim;

    PDM_g_num_t *copied_g_num = pts_g_num1 + idx_pts1[2];
    double      *copied_coord = pts_coord1 + idx_pts1[2] * dim;

    for (int i = 0; i < n_pts; i++) {
      int rank = rank_pt[i];

      if (rank == i_rank) {
        local_g_num[n_pts_local1] = pts_g_num[i];
        for (int j = 0; j < dim; j++) {
          local_coord[dim*n_pts_local1 + j] = pts_coord[dim*i + j];
        }
        n_pts_local1++;
      }

      else if (i_copied_rank1[rank] >= 0) {
        rank = i_copied_rank1[rank];
        int k = copied_shift1[rank] + copied_count[rank];
        copied_g_num[k] = pts_g_num[i];
        for (int j = 0; j < dim; j++) {
          copied_coord[dim*k + j] = pts_coord[dim*i + j];
        }
        copied_count[rank]++;
      }

      else {
        int k = send_shift[rank] + send_count[rank];
        send_g_num[k] = pts_g_num[i];
        for (int j = 0; j < dim; j++) {
          send_coord[dim*k + j] = pts_coord[dim*i + j];
        }
        send_count[rank]++;
      }
    }
    if (copied_count != NULL) {
      free (copied_count);
    }
    if (i_copied_rank1 != NULL) {
      free (i_copied_rank1);
    }
    free (rank_pt);

    /*   5) Send gnum buffer */
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       octree->comm);

    /*   6) Send coord buffer */
    for (int i = 0; i < n_rank; i++) {
      send_count[i] *= dim;
      recv_count[i] *= dim;
      send_shift[i+1] *= dim;
      recv_shift[i+1] *= dim;
    }
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       octree->comm);
  }

  /* Single proc */
  else {
    n_recv_pts = n_pts;
    n_pts_local1 = n_pts;
    n_pts_recv1 = 0;
    n_pts_copied1 = 0;

    idx_pts1[1] = n_pts_local1;
    idx_pts1[2] = idx_pts1[1] + n_pts_recv1;
    idx_pts1[3] = idx_pts1[2] + n_pts_copied1;
    n_pts1 = idx_pts1[3];

    pts_g_num1 = pts_g_num;
    pts_coord1 = pts_coord;
  }
  printf ("[%4d] phase 1: n_recv_pts = %8d (wihtout copies: %8d)\n", i_rank, n_pts1, n_recv_pts);

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[2] += e_t_elapsed - b_t_elapsed;
    times_cpu[2]     += e_t_cpu - b_t_cpu;
    times_cpu_u[2]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[2]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  /********************************************
   * First guess : closest source point in the
   * sense of Morton code
   ********************************************/
  double      *_closest_pt_dist2 = NULL;
  PDM_g_num_t *_closest_pt_g_num = NULL;

  if (n_rank == 1) {
    _closest_pt_g_num = closest_octree_pt_g_num;
    _closest_pt_dist2 = closest_octree_pt_dist2;
  }
  else {
    _closest_pt_dist2 = malloc (sizeof(double)      * n_pts1);
    _closest_pt_g_num = malloc (sizeof(PDM_g_num_t) * n_pts1);
  }

  /* Encode the coordinates of the received target points */
  pts_code = malloc (sizeof(PDM_morton_code_t) * n_pts1);
  _morton_encode_coords (dim,
                         PDM_morton_max_level,
                         octree->global_extents,
                         (size_t) n_pts1,
                         pts_coord1,
                         pts_code,
                         d,
                         s);

  int _n_src = octree->n_points;

  if (_n_src == 0) {
    for (int i = 0; i < idx_pts1[2]; i++) {
      _closest_pt_dist2[i] = HUGE_VAL;
      _closest_pt_g_num[i] = -1;
    }
  }

  else {
    _l_octant_t       *_octants   = octree->octants;
    double            *_src_coord = octree->points;
    PDM_g_num_t       *_src_g_num = octree->points_gnum;
    PDM_morton_code_t *_src_codes = octree->points_code;

    for (int i = 0; i < idx_pts1[2]; i++) {

      double *_coord = pts_coord1 + dim * i;

      /* Find base leaf */
      int base = PDM_morton_binary_search (_octants->n_nodes,
                                           pts_code[i],
                                           _octants->codes);

      if (_octants->n_points[base] > 0) {
        _closest_pt_dist2[i] = HUGE_VAL;
        for (int k = 0; k < _octants->n_points[base]; k++) {
          int j = _octants->range[base] + k;
          double dist2 = _pt_to_pt_dist2 (dim,
                                          _coord,
                                          _src_coord + dim*j);
          if (dist2 < _closest_pt_dist2[i]) {
            _closest_pt_dist2[i] = dist2;
            _closest_pt_g_num[i] = _src_g_num[j];
          }
        }
      }

      else {
        /* Last source point with Morton code lower than current target point's code */
        int j = PDM_morton_binary_search (_n_src,
                                          pts_code[i],
                                          _src_codes);

        _closest_pt_dist2[i] = _pt_to_pt_dist2 (dim,
                                                _coord,
                                                _src_coord + dim*j);
        _closest_pt_g_num[i] = _src_g_num[j];
        /* Check next source point as well (if in same rank) */
        if (j+1 < _n_src) {
          double dist2 =  _pt_to_pt_dist2 (dim,
                                           _coord,
                                           _src_coord + dim*(j+1));
          if (dist2 < _closest_pt_dist2[i]) {
            _closest_pt_dist2[i] = dist2;
            _closest_pt_g_num[i] = _src_g_num[j+1];
          }
        }
      }

    }
  }



  double *_pts_coord1 = pts_coord1 + idx_pts1[2] * dim;
  PDM_morton_code_t *_pts_code = pts_code + idx_pts1[2];
  double      *__closest_pt_dist2 = _closest_pt_dist2 + idx_pts1[2];
  PDM_g_num_t *__closest_pt_g_num = _closest_pt_g_num + idx_pts1[2];

  for (int rank = 0; rank < octree->n_copied_ranks; rank++) {

    _n_src = octree->n_copied_points[rank];

    if (_n_src == 0) {
      for (int i = copied_shift1[rank]; i < copied_shift1[rank+1]; i++) {
        __closest_pt_dist2[i] = HUGE_VAL;
        __closest_pt_g_num[i] = -1;
      }
    }

    else {
      _l_octant_t       *_octants   = octree->copied_octants[rank];
      double            *_src_coord = octree->copied_points[rank];
      PDM_g_num_t       *_src_g_num = octree->copied_points_gnum[rank];
      PDM_morton_code_t *_src_codes = octree->copied_points_code[rank];

      for (int i = copied_shift1[rank]; i < copied_shift1[rank+1]; i++) {
        double *_coord = _pts_coord1 + dim * i;

        /* Find base leaf */
        int base = PDM_morton_binary_search (_octants->n_nodes,
                                             _pts_code[i],
                                             _octants->codes);

        if (_octants->n_points[base] > 0) {
          __closest_pt_dist2[i] = HUGE_VAL;
          for (int k = 0; k < _octants->n_points[base]; k++) {
            int j = _octants->range[base] + k;
            double dist2 = _pt_to_pt_dist2 (dim,
                                            _coord,
                                            _src_coord + dim*j);
            if (dist2 < _closest_pt_dist2[i]) {
              __closest_pt_dist2[i] = dist2;
              __closest_pt_g_num[i] = _src_g_num[j];
            }
          }
        }

        else {
          /* Last source point with Morton code lower than current target point's code */
          int j = PDM_morton_binary_search (_n_src,
                                            _pts_code[i],
                                            _src_codes);

          __closest_pt_dist2[i] = _pt_to_pt_dist2 (dim,
                                                   _coord,
                                                   _src_coord + dim*j);
          __closest_pt_g_num[i] = _src_g_num[j];
          /* Check next source point as well (if in same rank) */
          if (j+1 < _n_src) {
            double dist2 =  _pt_to_pt_dist2 (dim,
                                             _coord,
                                             _src_coord + dim*(j+1));
            if (dist2 < _closest_pt_dist2[i]) {
              __closest_pt_dist2[i] = dist2;
              __closest_pt_g_num[i] = _src_g_num[j+1];
            }
          }
        }
      }
    }
  }
  free (pts_code);

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[3] += e_t_elapsed - b_t_elapsed;
    times_cpu[3]     += e_t_cpu - b_t_cpu;
    times_cpu_u[3]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[3]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  /*
   *  Search closest point
   */
  if (octree->explicit_nodes_to_build) {
    _single_closest_point_explicit (octree,
                                    -1,
                                    idx_pts1[2],
                                    pts_coord1,
                                    _closest_pt_g_num,
                                    _closest_pt_dist2);
  }
  else {
    _single_closest_point (dim,
                           octree->d,
                           octree->s,
                           octree->octants,
                           octree->points_gnum,
                           octree->points,
                           idx_pts1[2],
                           pts_coord1,
                           _closest_pt_g_num,
                           _closest_pt_dist2);
  }

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[4] += e_t_elapsed - b_t_elapsed;
    times_cpu[4]     += e_t_cpu - b_t_cpu;
    times_cpu_u[4]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[4]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  for (int i = 0; i < octree->n_copied_ranks; i++) {
    if (octree->explicit_nodes_to_build) {
      _single_closest_point_explicit (octree,
                                      i,
                                      copied_shift1[i+1] - copied_shift1[i],
                                      _pts_coord1 + copied_shift1[i]*dim,
                                      __closest_pt_g_num + copied_shift1[i],
                                      __closest_pt_dist2 + copied_shift1[i]);
    }
    else {
      _single_closest_point (dim,
                             octree->d,
                             octree->s,
                             octree->copied_octants[i],
                             octree->copied_points_gnum[i],
                             octree->copied_points[i],
                             copied_shift1[i+1] - copied_shift1[i],
                             _pts_coord1 + copied_shift1[i]*dim,
                             __closest_pt_g_num + copied_shift1[i],
                             __closest_pt_dist2 + copied_shift1[i]);
    }
  }

  if (n_rank == 1) {
    if (DETAIL_TIMER)
      PDM_timer_free (timer);
    return;
  }

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[5] += e_t_elapsed - b_t_elapsed;
    times_cpu[5]     += e_t_cpu - b_t_cpu;
    times_cpu_u[5]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[5]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  /*
   *  Fill block data
   */
  PDM_g_num_t *block_distrib_idx =
    PDM_compute_uniform_entity_distribution_from_partition (octree->comm,
                                                            1,
                                                            &n_pts,
                                                            (const PDM_g_num_t **) &pts_g_num);

  PDM_part_to_block_t *ptb1 = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                         PDM_PART_TO_BLOCK_POST_MERGE,
                                                         1.,
                                                         &pts_g_num1,
                                                         block_distrib_idx,
                                                         &n_pts1,
                                                         1,
                                                         octree->comm);

  int *part_stride = PDM_array_const_int(n_pts1, 1);

  int *block_stride = NULL;
  double *block_closest_pt_dist2 = NULL;
  PDM_part_to_block_exch (ptb1,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &_closest_pt_dist2,
                          &block_stride,
                          (void **) &block_closest_pt_dist2);
  free (block_stride);

  PDM_g_num_t *block_closest_pt_g_num = NULL;
  PDM_part_to_block_exch (ptb1,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &_closest_pt_g_num,
                          &block_stride,
                          (void **) &block_closest_pt_g_num);
  free (block_stride);
  free (part_stride);

  ptb1 = PDM_part_to_block_free (ptb1);

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[6] += e_t_elapsed - b_t_elapsed;
    times_cpu[6]     += e_t_cpu - b_t_cpu;
    times_cpu_u[6]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[6]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  /*
   *  Find ranks that may contain closer points
   */
  int *close_ranks_idx = NULL;
  int *close_ranks = NULL;

  if (USE_SHARED_OCTREE) {
    assert (octree->shared_rank_idx != NULL);

    int *close_nodes_idx = NULL;
    int *close_nodes     = NULL;
    PDM_box_tree_closest_upper_bound_dist_boxes_get (bt_shared,
                                                     n_pts1,
                                                     pts_coord1,
                                                     _closest_pt_dist2,
                                                     &close_nodes_idx,
                                                     &close_nodes);
    int *tag_rank = PDM_array_zeros_int (n_rank);

    int tmp_size = 4 * n_pts1;
    close_ranks = malloc (sizeof(int) * tmp_size);
    close_ranks_idx = malloc (sizeof(int) * (n_pts1 + 1));
    close_ranks_idx[0] = 0;

    for (int i = 0; i < n_pts1; i++) {
      close_ranks_idx[i+1] = close_ranks_idx[i];

      for (int j = close_nodes_idx[i]; j < close_nodes_idx[i+1]; j++) {
        int inode = close_nodes[j];
        if (octree->shared_pts_n[inode] == 0) continue;

        int l = 0;
        int r = n_rank;
        while (l + 1 < r) {
          int m = l + (r - l)/2;
          if (inode < octree->shared_rank_idx[m])
            r = m;
          else
            l = m;
        }
        int rank = l;

        if (tag_rank[rank]) continue;

        if (tmp_size <= close_ranks_idx[i+1]) {
          tmp_size *= 2;
          close_ranks = realloc (close_ranks, sizeof(int) * tmp_size);
        }

        close_ranks[close_ranks_idx[i+1]++] = rank;
        tag_rank[rank] = 1;
      }

      for (int j = close_ranks_idx[i]; j < close_ranks_idx[i+1]; j++) {
        tag_rank[close_ranks[j]] = 0;
      }
    }

    free (close_nodes_idx);
    free (close_nodes);
    free (tag_rank);

  }

  else {
    PDM_box_tree_closest_upper_bound_dist_boxes_get (octree->bt_shared,
                                                     n_pts1,
                                                     pts_coord1,
                                                     _closest_pt_dist2,
                                                     &close_ranks_idx,
                                                     &close_ranks);
    if (octree->n_used_rank < n_rank) {
      for (int i = 0; i < close_ranks_idx[n_pts1]; i++) {
        int j = close_ranks[i];
        close_ranks[i] = octree->used_rank[j];
      }
    }
  }

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[7] += e_t_elapsed - b_t_elapsed;
    times_cpu[7]     += e_t_cpu - b_t_cpu;
    times_cpu_u[7]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[7]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  PDM_array_reset_int(send_count, n_rank, 0);

  for (int i = 0; i < idx_pts1[2]; i++) {
    for (int j = close_ranks_idx[i]; j < close_ranks_idx[i+1]; j++) {
      int rank = close_ranks[j];
      if (rank != i_rank) {
        send_count[rank]++;
      }
    }
  }

  int *_close_ranks_idx = close_ranks_idx + idx_pts1[2];
  for (int c = 0; c < octree->n_copied_ranks; c++) {
    int _rank = octree->copied_ranks[c];

    for (int i = copied_shift1[c]; i < copied_shift1[c+1]; i++) {
      for (int j = _close_ranks_idx[i]; j < _close_ranks_idx[i+1]; j++) {
        int rank = close_ranks[j];
        if (rank != _rank) {
          send_count[rank]++;
        }
      }
    }
  }

  PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                    recv_count, 1, PDM_MPI_INT,
                    octree->comm);

  for (int i = 0; i < n_rank; i++) {
    send_shift[i+1] = send_shift[i] + send_count[i];
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
  }

  n_recv_pts = recv_shift[n_rank];
  int n_recv_pts_no_copies = n_recv_pts;

  PDM_para_octree_free_copies (id);

  int n_copied_ranks2;
  int *copied_ranks2 = NULL;
  int *n_recv_pts_copied_ranks = NULL;
  int mean_n_recv_pts;
  _prepare_copies (octree->comm,
                   f_copy_threshold,
                   f_max_copy,
                   n_recv_pts,
                   &n_copied_ranks2,
                   &copied_ranks2,
                   &n_recv_pts_copied_ranks,
                   &mean_n_recv_pts);

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[8] += e_t_elapsed - b_t_elapsed;
    times_cpu[8]     += e_t_cpu - b_t_cpu;
    times_cpu_u[8]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[8]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  if (n_copied_ranks2 > 0) {
    if (i_rank == 0) {
      if (n_copied_ranks2 == 1) {
        printf("phase 2: 1 copied rank: %d\n", copied_ranks2[0]);
      }
      else {
        printf("phase 2: %d copied ranks:", n_copied_ranks2);
        for (int i = 0; i < n_copied_ranks2; i++) {
          printf(" %d", copied_ranks2[i]);
        }
        printf("\n");
      }
    }

    PDM_para_octree_copy_ranks (id,
                                n_copied_ranks2,
                                copied_ranks2);
  } else {
    if (i_rank == 0) printf("phase 2: 0 copied ranks\n");
  }

  int *i_copied_rank2 = PDM_array_const_int(n_rank, -1);
  int *copied_count = malloc (sizeof(int) * octree->n_copied_ranks);

  for (int i = 0; i < octree->n_copied_ranks; i++) {
    i_copied_rank2[octree->copied_ranks[i]] = i;
    copied_count[i] = 0;
  }

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[9] += e_t_elapsed - b_t_elapsed;
    times_cpu[9]     += e_t_cpu - b_t_cpu;
    times_cpu_u[9]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[9]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  int n_pts_local2 = 0;
  int n_pts_recv2 = 0;
  int n_pts_copied2 = 0;
  for (int i = 0; i < octree->n_copied_ranks; i++) {
    int rank = octree->copied_ranks[i];
    if (rank != i_rank) {
      int si = send_count[rank];

      si = PDM_MIN (si, PDM_MAX (0, (n_recv_pts_copied_ranks[i] - n_recv_pts)/2));
      if (i_copied_rank2[i_rank] < 0) {
        si = PDM_MIN (si, PDM_MAX (0, mean_n_recv_pts - n_recv_pts));
      }

      copied_count[i] = si;
      n_recv_pts += si;
    }
  }

  for (int i = 0; i < n_rank; i++) {
    send_shift[i+1] = send_shift[i];

    if (i == i_rank) {
      n_pts_local2 += send_count[i];
      send_count[i] = 0;
    }
    else if (i_copied_rank2[i] >= 0) {
      send_count[i] -= copied_count[i_copied_rank2[i]];
    }

    send_shift[i+1] += send_count[i];
  }

  if (n_recv_pts_copied_ranks != NULL) {
    free (n_recv_pts_copied_ranks);
  }

  int *copied_shift2 = malloc (sizeof(int) * (octree->n_copied_ranks + 1));
  int *copied_count_tmp = malloc (sizeof(int) * octree->n_copied_ranks);
  copied_shift2[0] = 0;
  for (int i = 0; i < octree->n_copied_ranks; i++) {
    copied_shift2[i+1] = copied_shift2[i] + copied_count[i];
    copied_count_tmp[i] = 0;
  }
  n_pts_copied2 = copied_shift2[octree->n_copied_ranks];


  /* Exchange new send/recv counts */
  PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                    recv_count, 1, PDM_MPI_INT,
                    octree->comm);

  for (int i = 0; i < n_rank; i++) {
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
    send_count[i] = 0;
  }
  n_pts_recv2 = recv_shift[n_rank];

  int idx_pts2[4];
  idx_pts2[0] = 0;
  idx_pts2[1] = n_pts_local2;
  idx_pts2[2] = idx_pts2[1] + n_pts_recv2;
  idx_pts2[3] = idx_pts2[2] + n_pts_copied2;
  int n_pts2 = idx_pts2[3];

  PDM_g_num_t *pts_g_num2 = malloc (sizeof(PDM_g_num_t) * n_pts2);
  double      *pts_coord2 = malloc (sizeof(double)      * n_pts2 * dim);
  double      *_closest_pt_dist22 = malloc (sizeof(double) * n_pts2);

  send_g_num = realloc (send_g_num, sizeof(PDM_g_num_t) * send_shift[n_rank]);
  int s_data = dim + 1;
  send_coord = realloc (send_coord, sizeof(double) * send_shift[n_rank] * s_data);

  recv_g_num = pts_g_num2 + idx_pts2[1];
  recv_coord = pts_coord2 + idx_pts2[1] * dim;

  PDM_g_num_t *local_g_num = pts_g_num2 + idx_pts2[0];
  double      *local_coord = pts_coord2 + idx_pts2[0] * dim;
  double      *local_closest_pt_dist2 = _closest_pt_dist22 + idx_pts2[0];

  PDM_g_num_t *copied_g_num = pts_g_num2 + idx_pts2[2];
  double      *copied_coord = pts_coord2 + idx_pts2[2] * dim;
  double      *copied_closest_pt_dist2 = _closest_pt_dist22 + idx_pts2[2];

  n_pts_local2 = 0;
  for (int i = 0; i < idx_pts1[2]; i++) {
    for (int j = close_ranks_idx[i]; j < close_ranks_idx[i+1]; j++) {
      int rank = close_ranks[j];
      if (rank != i_rank) {

        if (i_copied_rank2[rank] >= 0) {
          int rank2 = i_copied_rank2[rank];

          if (copied_count_tmp[rank2] < copied_count[rank2]) {
            int k = copied_shift2[rank2] + copied_count_tmp[rank2];
            copied_g_num[k] = pts_g_num1[i];
            for (int l = 0; l < dim; l++) {
              copied_coord[dim*k+l] = pts_coord1[dim*i+l];
            }
            copied_closest_pt_dist2[k] = _closest_pt_dist2[i];
            copied_count_tmp[rank2]++;
          }
          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = pts_g_num1[i];
            for (int l = 0; l < dim; l++) {
              send_coord[s_data*k+l] = pts_coord1[dim*i+l];
            }
            send_coord[s_data*k+dim] = _closest_pt_dist2[i];
            send_count[rank]++;
          }
        }

        else {
          int k = send_shift[rank] + send_count[rank];
          send_g_num[k] = pts_g_num1[i];
          for (int l = 0; l < dim; l++) {
            send_coord[s_data*k+l] = pts_coord1[dim*i+l];
          }
          send_coord[s_data*k+dim] = _closest_pt_dist2[i];
          send_count[rank]++;
        }

      }
    }
  }

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[10] += e_t_elapsed - b_t_elapsed;
    times_cpu[10]     += e_t_cpu - b_t_cpu;
    times_cpu_u[10]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[10]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  PDM_g_num_t *_pts_g_num1 = pts_g_num1 + idx_pts1[2];
  for (int c = 0; c < n_copied_ranks1; c++) {
    int rank1 = copied_ranks1[c];

    for (int i = copied_shift1[c]; i < copied_shift1[c+1]; i++) {
      for (int j = _close_ranks_idx[i]; j < _close_ranks_idx[i+1]; j++) {
        int rank = close_ranks[j];
        if (rank != rank1) {

          if (rank == i_rank) {
            local_g_num[n_pts_local2] = _pts_g_num1[i];
            for (int l = 0; l < dim; l++) {
              local_coord[dim*n_pts_local2 + l] = _pts_coord1[dim*i + l];
            }
            local_closest_pt_dist2[n_pts_local2] = __closest_pt_dist2[i];
            n_pts_local2++;
          }

          else if (i_copied_rank2[rank] >= 0) {
            int rank2 = i_copied_rank2[rank];

            if (copied_count_tmp[rank2] < copied_count[rank2]) {
              int k = copied_shift2[rank2] + copied_count_tmp[rank2];
              copied_g_num[k] = _pts_g_num1[i];
              for (int l = 0; l < dim; l++) {
                copied_coord[dim*k+l] = _pts_coord1[dim*i+l];
              }
              copied_closest_pt_dist2[k] = __closest_pt_dist2[i];
              copied_count_tmp[rank2]++;
            }
            else {
              int k = send_shift[rank] + send_count[rank];
              send_g_num[k] = _pts_g_num1[i];
              for (int l = 0; l < dim; l++) {
                send_coord[s_data*k+l] = _pts_coord1[dim*i+l];
              }
              send_coord[s_data*k+dim] = __closest_pt_dist2[i];
              send_count[rank]++;
            }
          }

          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = _pts_g_num1[i];
            for (int l = 0; l < dim; l++) {
              send_coord[s_data*k+l] = _pts_coord1[dim*i+l];
            }
            send_coord[s_data*k+dim] = __closest_pt_dist2[i];
            send_count[rank]++;
          }
        }
      }
    }
  }

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[11] += e_t_elapsed - b_t_elapsed;
    times_cpu[11]     += e_t_cpu - b_t_cpu;
    times_cpu_u[11]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[11]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  if (copied_shift1 != NULL) {
    free (copied_shift1);
  }
  if (copied_count != NULL) {
    free (copied_count);
  }
  if (copied_count_tmp != NULL) {
    free (copied_count_tmp);
  }
  if (i_copied_rank2 != NULL) {
    free (i_copied_rank2);
  }
  free (close_ranks_idx);
  free (close_ranks);
  free (pts_coord1);
  free (pts_g_num1);
  free (_closest_pt_dist2);

  PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                     recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                     octree->comm);
  free (send_g_num);

  printf ("[%4d] phase 2: n_recv_pts = %8d (wihtout copies: %8d)\n", i_rank, n_pts2, n_recv_pts_no_copies);

  for (int i = 0; i < n_rank; i++) {
    send_count[i] *= s_data;
    recv_count[i] *= s_data;
    send_shift[i+1] *= s_data;
    recv_shift[i+1] *= s_data;
  }
  double *recv_double = malloc (sizeof(double) * recv_shift[n_rank]);
  PDM_MPI_Alltoallv (send_coord,  send_count, send_shift, PDM_MPI_DOUBLE,
                     recv_double, recv_count, recv_shift, PDM_MPI_DOUBLE,
                     octree->comm);
  free (send_coord);

  double *recv_closest_pt_dist2 = _closest_pt_dist22 + idx_pts2[1];
  int idx1 = 0, idx2 = 0;
  for (int i = 0; i < n_pts_recv2; i++) {
    for (int j = 0; j < dim; j++) {
      recv_coord[idx1++] = recv_double[idx2++];
    }
    recv_closest_pt_dist2[i] = recv_double[idx2++];
  }
  free (recv_double);
  free (send_count);
  free (send_shift);
  free (recv_count);
  free (recv_shift);

  _closest_pt_g_num = realloc (_closest_pt_g_num, sizeof(PDM_g_num_t) * n_pts2);
  PDM_array_reset_gnum(_closest_pt_g_num, n_pts2, -1);

  if (octree->explicit_nodes_to_build) {
    _single_closest_point_explicit (octree,
                                    -1,
                                    idx_pts2[2],
                                    pts_coord2,
                                    _closest_pt_g_num,
                                    _closest_pt_dist22);
  }
  else {
    _single_closest_point (dim,
                           octree->d,
                           octree->s,
                           octree->octants,
                           octree->points_gnum,
                           octree->points,
                           idx_pts2[2],
                           pts_coord2,
                           _closest_pt_g_num,
                           _closest_pt_dist22);
  }

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[12] += e_t_elapsed - b_t_elapsed;
    times_cpu[12]     += e_t_cpu - b_t_cpu;
    times_cpu_u[12]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[12]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  double *_pts_coord2 = pts_coord2 + idx_pts2[2] * dim;
  __closest_pt_dist2 = _closest_pt_dist22 + idx_pts2[2];
  __closest_pt_g_num = _closest_pt_g_num + idx_pts2[2];
  for (int i = 0; i < octree->n_copied_ranks; i++) {
    if (octree->explicit_nodes_to_build) {
      _single_closest_point_explicit (octree,
                                      i,
                                      copied_shift2[i+1] - copied_shift2[i],
                                      _pts_coord2 + copied_shift2[i]*dim,
                                      __closest_pt_g_num + copied_shift2[i],
                                      __closest_pt_dist2 + copied_shift2[i]);
    }
    else {
      _single_closest_point (dim,
                             octree->d,
                             octree->s,
                             octree->copied_octants[i],
                             octree->copied_points_gnum[i],
                             octree->copied_points[i],
                             copied_shift2[i+1] - copied_shift2[i],
                             _pts_coord2 + copied_shift2[i]*dim,
                             __closest_pt_g_num + copied_shift2[i],
                             __closest_pt_dist2 + copied_shift2[i]);
    }
  }
  free (pts_coord2);

  if (copied_ranks2 != NULL) {
    free (copied_ranks2);
  }
  free (copied_shift2);

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);

    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[13] += e_t_elapsed - b_t_elapsed;
    times_cpu[13]     += e_t_cpu - b_t_cpu;
    times_cpu_u[13]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[13]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(timer);
  }

  /*
   *  End of phase 2 -- Back to original partitioning
   */
  /* 1) Part-to-block */
  ptb1 = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                    PDM_PART_TO_BLOCK_POST_MERGE,
                                    1.,
                                    &pts_g_num2,
                                    block_distrib_idx,
                                    &n_pts2,
                                    1,
                                    octree->comm);
  part_stride = PDM_array_const_int(n_pts2, 1);

  double *tmp_block_closest_pt_dist2 = NULL;
  PDM_part_to_block_exch (ptb1,
                          sizeof(double),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &_closest_pt_dist22,
                          &block_stride,
                          (void **) &tmp_block_closest_pt_dist2);
  free (block_stride);
  free (_closest_pt_dist22);

  PDM_g_num_t *tmp_block_closest_pt_g_num = NULL;
  PDM_part_to_block_exch (ptb1,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stride,
                          (void **) &_closest_pt_g_num,
                          &block_stride,
                          (void **) &tmp_block_closest_pt_g_num);
  free (_closest_pt_g_num);
  free (part_stride);

  /* Merge block data (keep closest point if multiple candidates) */
  PDM_g_num_t *block_g_num1 = PDM_part_to_block_block_gnum_get (ptb1);
  const int n_pts_block1 = PDM_part_to_block_n_elt_block_get (ptb1);

  int k = 0;
  for (int i = 0; i < n_pts_block1; i++) {
    int l = (int) (block_g_num1[i] - block_distrib_idx[i_rank] - 1);

    for (int j = 0; j < block_stride[i]; j++) {
      if (tmp_block_closest_pt_g_num[k] > 0 &&
          tmp_block_closest_pt_dist2[k] < block_closest_pt_dist2[l]) {
        block_closest_pt_dist2[l] = tmp_block_closest_pt_dist2[k];
        block_closest_pt_g_num[l] = tmp_block_closest_pt_g_num[k];
      }
      k++;
    }
  }

  free (block_stride);
  free (tmp_block_closest_pt_dist2);
  free (tmp_block_closest_pt_g_num);


  /* 2) Block-to-part */
  PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                       (const PDM_g_num_t **) &pts_g_num,
                                                       &n_pts,
                                                       1,
                                                       octree->comm);
  int stride = 1;
  PDM_block_to_part_exch (btp,
                          sizeof(double),
                          PDM_STRIDE_CST,
                          &stride,
                          block_closest_pt_dist2,
                          NULL,
                          (void **) &closest_octree_pt_dist2);
  free (block_closest_pt_dist2);

  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          &stride,
                          block_closest_pt_g_num,
                          NULL,
                          (void **) &closest_octree_pt_g_num);
  free (block_closest_pt_g_num);
  free (pts_g_num2);

  ptb1 = PDM_part_to_block_free (ptb1);
  btp = PDM_block_to_part_free (btp);
  free (block_distrib_idx);

  if (copied_ranks1 != NULL) {
    free (copied_ranks1);
  }

  if (USE_SHARED_OCTREE) {
    PDM_box_set_destroy (&box_set);
    PDM_MPI_Comm_free (&bt_comm);
    PDM_box_tree_destroy (&bt_shared);
  }

  if (DETAIL_TIMER) {
    PDM_MPI_Barrier(octree->comm);
    PDM_timer_hang_on(timer);
    e_t_elapsed = PDM_timer_elapsed(timer);
    e_t_cpu     = PDM_timer_cpu(timer);
    e_t_cpu_u   = PDM_timer_cpu_user(timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(timer);

    times_elapsed[14] += e_t_elapsed - b_t_elapsed;
    times_cpu[14]     += e_t_cpu - b_t_cpu;
    times_cpu_u[14]   += e_t_cpu_u - b_t_cpu_u;
    times_cpu_s[14]   += e_t_cpu_s - b_t_cpu_s;

    if (i_rank == 0) {
      for (int i = 0; i < ntimer; i++) {
        printf("timer single point step %d : %12.5e\n",i,times_elapsed[i]);

      }
    }

    PDM_timer_free(timer);
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
PDM_para_octree_dump_times
(
 const int id
 )
{
  _octree_t *octree = _get_from_id (id);

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

    PDM_printf( "PDM_para_octree timer : build octree : build explicit nodes (elapsed and cpu)           :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_EXPLICIT_NODES],
                t_cpu_max[BUILD_EXPLICIT_NODES]);

  }

}


/**
 *
 * Get points located inside a set of boxes
 *
 * \param [in]   id                     Identifier
 * \param [in]   n_boxes                Number of boxes
 * \param [in]   box_extents            Extents of boxes
 * \param [in]   box_g_num              Global numbers of boxes
 * \param [out]  pts_in_box_idx         Index of points located in boxes
 * \param [out]  pts_in_box_g_num       Global numbers of points located in boxes
 * \param [out]  pts_in_box_coord       Coordinates of points located in boxes
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

  int n_recv_boxes = 0;
  PDM_morton_code_t *box_corners = NULL;
  double *recv_box_extents       = NULL;
  const PDM_g_num_t *recv_box_g_num;
  PDM_g_num_t *_recv_box_g_num   = NULL;
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
    sprintf(filename, "octree_points_%3.3d.vtk", my_rank);
    _export_octree_points (filename, octree, 0);
  }

  if (DEBUG) {
    if (my_rank == 0) {
      printf("octree extents = %f %f %f %f %f %f\n",
             octree->global_extents[0], octree->global_extents[1],
             octree->global_extents[2], octree->global_extents[3],
             octree->global_extents[4], octree->global_extents[5]);
    }
  }

  int NEW_PTS_IN_BOXES = 0;
  char *env_var = NULL;
  env_var = getenv ("NEW_PTS_IN_BOXES");
  if (env_var != NULL) {
    NEW_PTS_IN_BOXES = atoi(env_var);
  }
  if (my_rank == 0) {
    printf("NEW_PTS_IN_BOXES = %d\n", NEW_PTS_IN_BOXES);
  }

  /* Multiple ranks */
  if (n_ranks > 1) {
    block_distrib_idx = PDM_compute_uniform_entity_distribution_from_partition(octree->comm,
                                                                               1,
                                                                               &n_boxes,
                                                                               &box_g_num);


    /***************************************
     * Redistribute bounding boxes
     ***************************************/
    int *send_count = PDM_array_zeros_int(n_ranks);

    /* Encode box corners */
    box_corners = malloc (sizeof(PDM_morton_code_t) * 2 * n_boxes);
    _morton_encode_coords (dim,
                           PDM_morton_max_level,
                           octree->global_extents,
                           2 * n_boxes,
                           box_extents,
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

    printf("[%d] octree->n_pts = %d, n_recv_boxes = %d\n", my_rank, octree->n_points, n_recv_boxes);

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
          send_box_extents[two_dim*idx + k] = box_extents[two_dim*ibox + k];
        }

        send_count[irank]++;
      }
    }
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
    recv_box_extents = (double *) box_extents;
    recv_box_g_num   = box_g_num;
  }


  /***************************************
   * Intersect redistributed boxes with local octree
   ***************************************/
  if (VISU) {
    char filename[999];
    sprintf(filename, "recv_boxes0_%3.3d.vtk", my_rank);
    _export_boxes (filename, n_recv_boxes, recv_box_extents, recv_box_g_num);
  }

  /* Encode corners of redistributed boxes */
  box_corners = malloc (sizeof(PDM_morton_code_t) * 2 * n_recv_boxes);
  _morton_encode_coords (dim,
                         PDM_morton_max_level,
                         octree->global_extents,
                         2 * n_recv_boxes,
                         recv_box_extents,
                         box_corners,
                         d,
                         s);

  int *box_pts_n = NULL;
  int *box_pts_idx = NULL;
  int *box_pts = NULL;
  if (NEW_PTS_IN_BOXES) {
    _points_inside_boxes (octree,
                          -1,
                          n_recv_boxes,
                          recv_box_extents,
                          box_corners,
                          &box_pts_idx,
                          &box_pts);
    box_pts_n = malloc (sizeof(int) * n_recv_boxes);
    for (int i = 0; i < n_recv_boxes; i++) {
      box_pts_n[i] = box_pts_idx[i+1] - box_pts_idx[i];
    }

  } else {
    /* Root node of octree */
    PDM_morton_code_t root;
    root.L = 0;
    root.X[0] = 0;
    root.X[1] = 0;
    root.X[2] = 0;

    if (octants->n_nodes > 0) {
      PDM_morton_nearest_common_ancestor (octants->codes[0],
                                          octants->codes[octants->n_nodes - 1],
                                          &root);
    }

    int *intersect_nodes = malloc (sizeof(int) * octants->n_nodes);
    size_t n_intersect_nodes;

    size_t s_box_pts = octree->n_points;
    box_pts = malloc (sizeof(int) * s_box_pts);
    box_pts_n = malloc (sizeof(int) * n_recv_boxes);

    box_pts_idx = malloc (sizeof(int) * (n_recv_boxes+1));
    box_pts_idx[0] = 0;

    size_t min_intersect = 100000000;
    size_t max_intersect = 0;
    size_t mean_intersect = 0;
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
                                octants->n_points,
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

      // printf("n_intersect_nodes =%li \n", n_intersect_nodes);
      min_intersect = PDM_MIN(min_intersect, n_intersect_nodes);
      max_intersect = PDM_MAX(max_intersect, n_intersect_nodes);
      mean_intersect += n_intersect_nodes;
      //size_t new_max_size = box_pts_idx[ibox] + n_intersect_nodes * octree->points_in_leaf_max;
      size_t new_max_size = box_pts_idx[ibox];
      for (size_t i = 0; i < n_intersect_nodes; i++) {
        new_max_size += octants->n_points[intersect_nodes[i]];
      }

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
    free (intersect_nodes);
  }
  free (box_corners);
  if (recv_box_extents != box_extents) free (recv_box_extents);


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
    // PDM_MPI_Barrier(octree->comm);
    double t1 = PDM_MPI_Wtime();
    PDM_part_to_block_t *ptb2 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                          (PDM_g_num_t **) (&recv_box_g_num),
                                                          NULL,
                                                          &n_recv_boxes,
                                                          1,
                                                          octree->comm);
    double dt = PDM_MPI_Wtime() - t1;
    log_debug(" dt = %12.5e \n",  dt);
    if (_recv_box_g_num != box_g_num) free (_recv_box_g_num);

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
    //free (box_pts_n_coord);

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

    *pts_in_box_idx = PDM_array_new_idx_from_sizes_int(pts_in_box_n, n_boxes);
    *pts_in_box_g_num = malloc (sizeof(PDM_g_num_t) * (*pts_in_box_idx)[n_boxes]);

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

    *pts_in_box_coord = malloc (sizeof(double) * (*pts_in_box_idx)[n_boxes] * dim);

    PDM_block_to_part_exch (btp,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            block_pts_in_box_n_coord,
                            (void *) block_pts_in_box_coord,
                            &pts_in_box_n_coord,
                            (void **) pts_in_box_coord);

    free (block_distrib_idx);
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



#define NTIMER_PIB 8

typedef enum {
  PIB_BEGIN,
  PIB_REDISTRIBUTE,
  PIB_COPIES,
  PIB_EXCHANGE,
  PIB_LOCAL,
  PIB_PTB,
  PIB_BTP,
  PIB_TOTAL
} _pib_step_t;

#define PIB_TIME_FMT "%f" //"12.5e"



void
PDM_para_octree_points_inside_boxes_with_copies
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
  int DEBUG = 0;
  float f_copy_threshold = 1.05;
  float f_max_copy = 0.05;

  char *env_var = NULL;
  env_var = getenv ("OCTREE_COPY_THRESHOLD");
  if (env_var != NULL) {
    f_copy_threshold = (float) atof(env_var);
  }

  env_var = getenv ("OCTREE_MAX_COPY");
  if (env_var != NULL) {
    f_max_copy = (float) atof(env_var);
  }


  int USE_SHARED_OCTREE = 1;
  env_var = getenv ("USE_SHARED_OCTREE");
  if (env_var != NULL) {
    USE_SHARED_OCTREE = (float) atof(env_var);
  }

  _octree_t *octree = _get_from_id (octree_id);
  const int dim = octree->dim;
  const int two_dim = 2 * dim;


  int i_rank, n_rank;
  PDM_MPI_Comm_rank (octree->comm, &i_rank);
  PDM_MPI_Comm_size (octree->comm, &n_rank);

  if (i_rank == 0) {
    printf("USE_SHARED_OCTREE = %d\n", USE_SHARED_OCTREE);
  }

  printf("[%d] n_boxes = %d\n", i_rank, n_boxes);

  double times_elapsed[NTIMER_PIB], b_t_elapsed, e_t_elapsed;
  for (_pib_step_t step = PIB_BEGIN; step <= PIB_TOTAL; step++) {
    times_elapsed[step] = 0.;
  }

  PDM_timer_hang_on (octree->timer);
  times_elapsed[PIB_BEGIN] = PDM_timer_elapsed (octree->timer);
  b_t_elapsed = times_elapsed[PIB_BEGIN];
  PDM_timer_resume (octree->timer);


  PDM_morton_code_t *box_corners = NULL;
  double d[3], s[3];

  int n_copied_ranks = 0;
  int *copied_ranks = NULL;

  int n_box_local;
  int n_box_recv;
  int n_box_copied;
  int n_box1;

  int *copied_shift = NULL;

  PDM_g_num_t *box_g_num1   = NULL;
  double      *box_extents1 = NULL;

  if (n_rank > 1) {
    /*
     *  Redistribute boxes
     */
    /* Encode box corners */
    box_corners = malloc (sizeof(PDM_morton_code_t) * 2 * n_boxes);
    _morton_encode_coords (dim,
                           PDM_morton_max_level,
                           octree->global_extents,
                           2 * n_boxes,
                           box_extents,
                           box_corners,
                           d,
                           s);


    /* Find which ranks possibly intersect each box */
    int *send_count = PDM_array_zeros_int (n_rank);

    int tmp_size = 4 * n_boxes;
    int *box_rank = malloc (sizeof(int) * tmp_size);
    int *box_rank_idx = malloc (sizeof(int) * (n_boxes + 1));
    box_rank_idx[0] = 0;

    if (USE_SHARED_OCTREE) {
      assert (octree->shared_rank_idx != NULL);

      size_t n_intersect_nodes;
      int *intersect_nodes = malloc (sizeof(int) * octree->shared_rank_idx[n_rank]);

      PDM_morton_code_t root;
      root.L    = 0;
      root.X[0] = 0;
      root.X[1] = 0;
      root.X[2] = 0;

      int *tag_rank = PDM_array_zeros_int (n_rank);

      for (int ibox = 0; ibox < n_boxes; ibox++) {
        box_rank_idx[ibox+1] = box_rank_idx[ibox];

        const double *box_min = box_extents + two_dim*ibox;
        const double *box_max = box_min + dim;
        if (DEBUG) {
          printf("\n[%d] box %d ("PDM_FMT_G_NUM") extents %f %f %f %f %f %f\n",
                 i_rank, ibox, box_g_num[ibox],
                 box_min[0], box_min[1], box_min[2],
                 box_max[0], box_max[1], box_max[2]);
          printf("[%d] code min/max: L = %u, X = %u %u %u / L = %u, X = %u %u %u\n",
                 i_rank,
                 box_corners[2*ibox].L, box_corners[2*ibox].X[0], box_corners[2*ibox].X[1], box_corners[2*ibox].X[2],
                 box_corners[2*ibox+1].L, box_corners[2*ibox+1].X[0], box_corners[2*ibox+1].X[1], box_corners[2*ibox+1].X[2]);

        }

        n_intersect_nodes = 0;
        PDM_morton_intersect_box (dim,
                                  root,
                                  box_corners[2*ibox],
                                  box_corners[2*ibox+1],
                                  octree->shared_codes,
                                  octree->shared_pts_n,
                                  0,
                                  octree->shared_rank_idx[n_rank],
                                  &n_intersect_nodes,
                                  intersect_nodes);

        for (size_t i = 0; i < n_intersect_nodes; i++) {
          int inode = intersect_nodes[i];

          int l = 0;
          int r = n_rank;
          while (l + 1 < r) {
            int m = l + (r - l)/2;
            if (inode < octree->shared_rank_idx[m])
              r = m;
            else
              l = m;
          }
          int rank = l;
          if (DEBUG) {
            printf("[%d]  intersects shared node %d (rank %d)\n", i_rank, inode, rank);
          }

          if (tag_rank[rank]) continue;

          int intersect = 1;
          double *node_min = octree->shared_pts_extents + 6*inode;
          double *node_max = node_min + 3;

          for (int j = 0; j < dim; j++) {
            if (box_min[j] > node_max[j] || box_max[j] < node_min[j]) {
              intersect = 0;
              break;
            }
          }

          if (DEBUG) {
            printf("[%d]    intersects node pts extents? %d\n", i_rank, intersect);
          }

          if (intersect) {
            if (tmp_size <= box_rank_idx[ibox+1]) {
              tmp_size *= 2;
              box_rank = realloc (box_rank, sizeof(int) * tmp_size);
            }
            box_rank[box_rank_idx[ibox+1]++] = rank;
            tag_rank[rank] = 1;
            send_count[rank]++;
          }
        }

        for (int i = box_rank_idx[ibox]; i < box_rank_idx[ibox+1]; i++) {
          tag_rank[box_rank[i]] = 0;
        }
      }
      free (tag_rank);
      free (intersect_nodes);
    }

    else {
      size_t start, end, tmp;
      for (int ibox = 0; ibox < n_boxes; ibox++) {
        PDM_morton_quantile_intersect (n_rank,
                                       box_corners[2*ibox],
                                       octree->rank_octants_index,
                                       &start,
                                       &tmp);

        PDM_morton_quantile_intersect (n_rank - start,
                                       box_corners[2*ibox+1],
                                       octree->rank_octants_index + start,
                                       &tmp,
                                       &end);
        end += start;

        int new_size = box_rank_idx[ibox] + (int) (end - start);
        if (tmp_size <= new_size) {
          tmp_size = PDM_MAX (2*tmp_size, new_size);
          box_rank = realloc (box_rank, sizeof(int) * tmp_size);
        }

        box_rank_idx[ibox+1] = box_rank_idx[ibox];
        for (size_t irank = start; irank < end; irank++) {
          box_rank[box_rank_idx[ibox+1]++] = (int) irank;
          send_count[irank]++;
        }
      }

    }
    free (box_corners);


    PDM_timer_hang_on (octree->timer);
    e_t_elapsed = PDM_timer_elapsed (octree->timer);
    times_elapsed[PIB_REDISTRIBUTE] = e_t_elapsed - b_t_elapsed;
    b_t_elapsed = e_t_elapsed;
    PDM_timer_resume (octree->timer);

    //-->>
    if (0) {
      printf("[%d] --- Box Rank ---\n", i_rank);
      for (int i = 0; i < n_boxes; i++) {
        printf("[%d] %d ("PDM_FMT_G_NUM") : ", i_rank, i, box_g_num[i]);
        for (int j = box_rank_idx[i]; j < box_rank_idx[i+1]; j++) {
          printf("%d ", box_rank[j]);
        }
        printf("\n");
      }
      printf("[%d] ------------------\n", i_rank);
    }
    //<<--

    /* Exchange provisional send/recv counts */
    int *recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      octree->comm);

    int n_recv_box = 0;
    for (int i = 0; i < n_rank; i++) {
      n_recv_box += recv_count[i];
    }
    int n_recv_box_no_copies = n_recv_box;

    int *n_recv_box_copied_ranks = NULL;
    int avg_n_recv_box;
    _prepare_copies (octree->comm,
                     f_copy_threshold,
                     f_max_copy,
                     n_recv_box,
                     &n_copied_ranks,
                     &copied_ranks,
                     &n_recv_box_copied_ranks,
                     &avg_n_recv_box);

    if (n_copied_ranks > 0) {
      if (i_rank == 0) {
        if (n_copied_ranks == 1) {
          printf("1 copied rank: %d\n", copied_ranks[0]);
        }
        else {
          printf("%d copied ranks:", n_copied_ranks);
          for (int i = 0; i < n_copied_ranks; i++) {
            printf(" %d", copied_ranks[i]);
          }
          printf("\n");
        }
      }

      PDM_para_octree_copy_ranks (octree_id,
                                  n_copied_ranks,
                                  copied_ranks);
      free (copied_ranks);
    } else {
      if (i_rank == 0) printf("0 copied ranks\n");
    }

    int *i_copied_rank = PDM_array_const_int (n_rank, -1);
    int *copied_count = malloc (sizeof(int) * octree->n_copied_ranks);

    for (int i = 0; i < octree->n_copied_ranks; i++) {
      i_copied_rank[octree->copied_ranks[i]] = i;
      copied_count[i] = 0;
    }

    int *send_shift = PDM_array_new_idx_from_sizes_int (send_count, n_rank);

    n_box_local = 0;
    for (int i = 0; i < octree->n_copied_ranks; i++) {
      int rank = octree->copied_ranks[i];
      if (rank != i_rank) {
        int si = send_count[rank];

        si = PDM_MIN (si, PDM_MAX (0, (n_recv_box_copied_ranks[i] - n_recv_box)/2));
        if (i_copied_rank[i_rank] < 0) {
          si = PDM_MIN (si, PDM_MAX (0, avg_n_recv_box - n_recv_box));
        }

        copied_count[i] = si;
        n_recv_box += si;
      }
    }

    if (n_recv_box_copied_ranks != NULL) {
      free (n_recv_box_copied_ranks);
    }

    for (int i = 0; i < n_rank; i++) {

      send_shift[i+1] = send_shift[i];

      if (i == i_rank) {
        n_box_local += send_count[i];
        send_count[i] = 0;
      }
      else if (i_copied_rank[i] >= 0) {
        send_count[i] -= copied_count[i_copied_rank[i]];
      }

      send_shift[i+1] += send_count[i];
    }

    copied_shift = PDM_array_new_idx_from_sizes_int (copied_count, octree->n_copied_ranks);
    int *copied_count_tmp = PDM_array_zeros_int (octree->n_copied_ranks);
    n_box_copied = copied_shift[octree->n_copied_ranks];

    PDM_timer_hang_on (octree->timer);
    e_t_elapsed = PDM_timer_elapsed (octree->timer);
    times_elapsed[PIB_COPIES] = e_t_elapsed - b_t_elapsed;
    b_t_elapsed = e_t_elapsed;
    PDM_timer_resume (octree->timer);

    /* Exchange new send/recv counts */
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      octree->comm);

    int *recv_shift = PDM_array_new_idx_from_sizes_int (recv_count, n_rank);
    PDM_array_reset_int (send_count, n_rank, 0);

    n_box_recv = recv_shift[n_rank];

    n_box1 = n_box_local + n_box_recv + n_box_copied;
    printf("[%d] octree->n_pts = %d, n_recv_boxes = %d (without copies : %d)\n", i_rank, octree->n_points, n_box1, n_recv_box_no_copies);

    box_g_num1   = malloc (sizeof(PDM_g_num_t) * n_box1);
    box_extents1 = malloc (sizeof(double)      * n_box1 * two_dim);

    /* Fill send buffers */
    PDM_g_num_t *send_g_num   = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    double      *send_extents = malloc (sizeof(double)      * send_shift[n_rank] * two_dim);
    PDM_g_num_t *recv_g_num   = box_g_num1 + n_box_local;
    double      *recv_extents = box_extents1 + n_box_local * two_dim;

    int idx_copied = n_box_local + n_box_recv;
    PDM_g_num_t *copied_g_num   = box_g_num1   + idx_copied;
    double      *copied_extents = box_extents1 + idx_copied * two_dim;

    n_box_local = 0;
    for (int ibox = 0; ibox < n_boxes; ibox++) {
      for (int i = box_rank_idx[ibox]; i < box_rank_idx[ibox+1]; i++) {
        int rank = box_rank[i];

        if ((int) rank == i_rank) {
          box_g_num1[n_box_local] = box_g_num[ibox];
          for (int j = 0; j < two_dim; j++) {
            box_extents1[two_dim*n_box_local + j] = box_extents[two_dim*ibox + j];
          }
          n_box_local++;
        }

        else if (i_copied_rank[rank] >= 0) {
          int _rank = i_copied_rank[rank];

          if (copied_count_tmp[_rank] < copied_count[_rank]) {
            int k = copied_shift[_rank] + copied_count_tmp[_rank];
            copied_g_num[k] = box_g_num[ibox];
            for (int j = 0; j < two_dim; j++) {
              copied_extents[two_dim*k + j] = box_extents[two_dim*ibox + j];
            }
            copied_count_tmp[_rank]++;
          }
          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = box_g_num[ibox];
            for (int j = 0; j < two_dim; j++) {
              send_extents[two_dim*k + j] = box_extents[two_dim*ibox + j];
            }
            send_count[rank]++;
          }
        }

        else {
          int k = send_shift[rank] + send_count[rank];
          send_g_num[k] = box_g_num[ibox];
          for (int j = 0; j < two_dim; j++) {
            send_extents[two_dim*k + j] = box_extents[two_dim*ibox + j];
          }
          send_count[rank]++;
        }
      }
    }
    if (copied_count != NULL) {
      free (copied_count);
    }
    if (copied_count_tmp != NULL) {
      free (copied_count_tmp);
    }
    if (i_copied_rank != NULL) {
      free (i_copied_rank);
    }
    free (box_rank);
    free (box_rank_idx);

    /* Send boxes g_num buffer */
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       octree->comm);

    /* Send boxes extents buffer */
    for (int i = 0; i < n_rank; i++) {
      send_shift[i+1] *= two_dim;
      recv_shift[i+1] *= two_dim;
      send_count[i]   *= two_dim;
      recv_count[i]   *= two_dim;
    }
    PDM_MPI_Alltoallv (send_extents, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_extents, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       octree->comm);

    free (send_count);
    free (recv_count);
    free (send_shift);
    free (recv_shift);
    free (send_g_num);
    free (send_extents);

    PDM_timer_hang_on (octree->timer);
    e_t_elapsed = PDM_timer_elapsed (octree->timer);
    times_elapsed[PIB_EXCHANGE] = e_t_elapsed - b_t_elapsed;
    b_t_elapsed = e_t_elapsed;
    PDM_timer_resume (octree->timer);
  }

  /* Single proc */
  else {
    n_box_local  = n_boxes;
    n_box_recv   = 0;
    n_box_copied = 0;

    n_box1 = n_boxes;

    box_g_num1 = (PDM_g_num_t *) box_g_num;
    box_extents1 = (double *) box_extents;
  }


  /***************************************
   * Intersect redistributed boxes with local octree
   ***************************************/
  /* Encode corners of redistributed boxes */
  box_corners = malloc (sizeof(PDM_morton_code_t) * 2 * n_box1);
  _morton_encode_coords (dim,
                         PDM_morton_max_level,
                         octree->global_extents,
                         2 * n_box1,
                         box_extents1,
                         box_corners,
                         d,
                         s);

  /*
   *  Get points inside boxes
   */
  int n_part = 1 + octree->n_copied_ranks;
  int *part_n_box = malloc (sizeof(int) * n_part);
  part_n_box[0] = n_box_local + n_box_recv;

  int **box_pts_idx = malloc (sizeof(int *) * n_part);
  int **box_pts_l_num = malloc (sizeof(int *) * n_part);

  int size_box_pts = 0;

  /*
   *  Search in local tree
   */
  if (octree->explicit_nodes_to_build) {
    _points_inside_boxes_explicit (octree,
                                   -1,
                                   part_n_box[0],
                                   box_extents1,
                                   //box_corners,
                                   &(box_pts_idx[0]),
                                   &(box_pts_l_num[0]));
  }
  else {
    _points_inside_boxes (octree,
                          -1,
                          part_n_box[0],
                          box_extents1,
                          box_corners,
                          &(box_pts_idx[0]),
                          &(box_pts_l_num[0]));
  }
  size_box_pts += box_pts_idx[0][part_n_box[0]];

  /*
   *  Search in copied trees
   */
  if (octree->n_copied_ranks > 0) {
    double            *box_extents_copied = box_extents1 + part_n_box[0] * two_dim;
    PDM_morton_code_t *box_corners_copied = box_corners  + part_n_box[0] * 2;
    for (int i = 0; i < octree->n_copied_ranks; i++) {
      part_n_box[i+1] = copied_shift[i+1] - copied_shift[i];

      if (octree->explicit_nodes_to_build) {
        _points_inside_boxes_explicit (octree,
                                       i,
                                       part_n_box[i+1],
                                       box_extents_copied + copied_shift[i] * two_dim,
                                       //box_corners_copied + copied_shift[i] * 2,
                                       &(box_pts_idx[i+1]),
                                       &(box_pts_l_num[i+1]));
      }
      else {
        _points_inside_boxes (octree,
                              i,
                              part_n_box[i+1],
                              box_extents_copied + copied_shift[i] * two_dim,
                              box_corners_copied + copied_shift[i] * 2,
                              &(box_pts_idx[i+1]),
                              &(box_pts_l_num[i+1]));
      }

      size_box_pts += box_pts_idx[i+1][part_n_box[i+1]];
    }
  }
  if (copied_shift != NULL) free (copied_shift);
  if (box_extents1 != box_extents) free (box_extents1);
  free (box_corners);

  PDM_g_num_t *box_pts_g_num = malloc (sizeof(PDM_g_num_t) * size_box_pts);
  double      *box_pts_coord = malloc (sizeof(double)      * size_box_pts * dim);
  int idx = 0;
  for (int j = 0; j < box_pts_idx[0][part_n_box[0]]; j++) {
    box_pts_g_num[idx] = octree->points_gnum[box_pts_l_num[0][j]];
    for (int k = 0; k < 3; k++) {
      box_pts_coord[3*idx + k] = octree->points[3*box_pts_l_num[0][j] + k];
    }
    idx++;
  }

  for (int i = 0; i < octree->n_copied_ranks; i++) {
    for (int j = 0; j < box_pts_idx[i+1][part_n_box[i+1]]; j++) {
      box_pts_g_num[idx] = octree->copied_points_gnum[i][box_pts_l_num[i+1][j]];
      for (int k = 0; k < 3; k++) {
        box_pts_coord[3*idx + k] = octree->copied_points[i][3*box_pts_l_num[i+1][j] + k];
      }
      idx++;
    }
  }

  PDM_para_octree_free_copies (octree_id);

  PDM_timer_hang_on (octree->timer);
  e_t_elapsed = PDM_timer_elapsed (octree->timer);
  times_elapsed[PIB_LOCAL] = e_t_elapsed - b_t_elapsed;
  b_t_elapsed = e_t_elapsed;
  PDM_timer_resume (octree->timer);


  if (n_rank == 1) {
    *pts_in_box_g_num = box_pts_g_num;
    *pts_in_box_coord = box_pts_coord;

    *pts_in_box_idx   = malloc (sizeof(int) * (n_boxes + 1));
    memcpy (*pts_in_box_idx, box_pts_idx[0], sizeof(int) * (n_boxes + 1));

    free (box_pts_idx[0]);
    free (box_pts_l_num[0]);
    free (box_pts_idx);
    free (box_pts_l_num);
    free (part_n_box);
  }

  else {
    int *part_stride = malloc (sizeof(int) * n_box1);
    double *weight = malloc (sizeof(double) * n_box1);
    idx = 0;
    for (int i = 0; i < n_part; i++) {
      for (int j = 0; j < part_n_box[i]; j++) {
        part_stride[idx] = box_pts_idx[i][j+1] - box_pts_idx[i][j];
        weight[idx] = (double) part_stride[idx];
        idx++;
      }
      free (box_pts_idx[i]);
      free (box_pts_l_num[i]);
    }
    free (part_n_box);
    free (box_pts_idx);
    free (box_pts_l_num);

    /* Part#2 to Block */
    PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                         PDM_PART_TO_BLOCK_POST_MERGE,
                                                         1.,
                                                         (PDM_g_num_t **) &box_g_num1,
                                                         &weight,
                                                         &n_box1,
                                                         1,
                                                         octree->comm);

    PDM_g_num_t l_max_box_g_num = 0;
    for (int i = 0; i < n_boxes; i++) {
      l_max_box_g_num = PDM_MAX (l_max_box_g_num, box_g_num[i]);
    }
    PDM_g_num_t g_max_box_g_num;
    PDM_MPI_Allreduce (&l_max_box_g_num, &g_max_box_g_num, 1,
                       PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, octree->comm);
    //printf("g_max_box_g_num = "PDM_FMT_G_NUM"\n", g_max_box_g_num);
    free (weight);

    int *block_pts_in_box_n = NULL;
    PDM_g_num_t *block_pts_in_box_g_num = NULL;
    PDM_part_to_block_exch (ptb,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            1,
                            &part_stride,
                            (void **) &box_pts_g_num,
                            &block_pts_in_box_n,
                            (void **) &block_pts_in_box_g_num);
    free (box_pts_g_num);

    /*for (int i = 0; i < n_box1; i++) {
      part_stride[i] *= dim;
      }*/

    int *block_stride = NULL;
    double *block_pts_in_box_coord = NULL;
    PDM_part_to_block_exch (ptb,
                            dim*sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &part_stride,
                            (void **) &box_pts_coord,
                            &block_stride,
                            (void **) &block_pts_in_box_coord);
    free (box_pts_coord);
    free (part_stride);
    free (block_stride);

    PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);
    PDM_g_num_t *_block_distrib_idx = block_distrib_idx;
    if (block_distrib_idx[n_rank] < g_max_box_g_num) {
      _block_distrib_idx = malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
      for (int i = 0; i < n_rank; i++) {
        _block_distrib_idx[i] = block_distrib_idx[i];
      }
      _block_distrib_idx[n_rank] = g_max_box_g_num;
    }

    int n_elt_block = PDM_part_to_block_n_elt_block_get (ptb);
    int n_elt_block_full = (int) (_block_distrib_idx[i_rank+1] - _block_distrib_idx[i_rank]);

    if (n_elt_block < n_elt_block_full) {

      PDM_g_num_t *block_g_num = PDM_part_to_block_block_gnum_get (ptb);
      int *block_pts_in_box_n_full = PDM_array_zeros_int (n_elt_block_full);

      int i1 = 0;
      for (int i = 0; i < n_elt_block; i++) {
        while (_block_distrib_idx[i_rank] + 1 + i1 < block_g_num[i]) {
          i1++;
        }

        block_pts_in_box_n_full[i1] = block_pts_in_box_n[i];
      }

      free (block_pts_in_box_n);
      block_pts_in_box_n = block_pts_in_box_n_full;
    }
    free (box_g_num1);

    PDM_timer_hang_on (octree->timer);
    e_t_elapsed = PDM_timer_elapsed (octree->timer);
    times_elapsed[PIB_PTB] = e_t_elapsed - b_t_elapsed;
    b_t_elapsed = e_t_elapsed;
    PDM_timer_resume (octree->timer);

    /*
     *  Block to part
     */
    PDM_block_to_part_t *btp = PDM_block_to_part_create (_block_distrib_idx,
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

    *pts_in_box_idx = PDM_array_new_idx_from_sizes_int(pts_in_box_n, n_boxes);
    *pts_in_box_g_num = malloc (sizeof(PDM_g_num_t) * (*pts_in_box_idx)[n_boxes]);

    PDM_block_to_part_exch (btp,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            block_pts_in_box_n,
                            (void *) block_pts_in_box_g_num,
                            &pts_in_box_n,
                            (void **) pts_in_box_g_num);
    free (block_pts_in_box_g_num);


    /*for (int i = 0; i < n_elt_block_full; i++) {
      block_pts_in_box_n[i] *= dim;
    }

    for (int i = 0; i < n_boxes; i++) {
      pts_in_box_n[i] *= dim;
      }*/

    *pts_in_box_coord = malloc (sizeof(double) * (*pts_in_box_idx)[n_boxes] * dim);

    PDM_block_to_part_exch (btp,
                            dim*sizeof(double),
                            PDM_STRIDE_VAR,
                            block_pts_in_box_n,
                            (void *) block_pts_in_box_coord,
                            &pts_in_box_n,
                            (void **) pts_in_box_coord);
    free (block_pts_in_box_n);
    free (block_pts_in_box_coord);
    free (pts_in_box_n);

    if (_block_distrib_idx != block_distrib_idx) free (_block_distrib_idx);
    PDM_part_to_block_free (ptb);
    PDM_block_to_part_free (btp);

    PDM_timer_hang_on (octree->timer);
    e_t_elapsed = PDM_timer_elapsed (octree->timer);
    times_elapsed[PIB_BTP] = e_t_elapsed - b_t_elapsed;
    times_elapsed[PIB_TOTAL] = e_t_elapsed - times_elapsed[PIB_BEGIN];
    PDM_timer_resume (octree->timer);
  }

  printf ("[%d] PiB timers: "PIB_TIME_FMT" "PIB_TIME_FMT" "PIB_TIME_FMT" "PIB_TIME_FMT" "PIB_TIME_FMT" "PIB_TIME_FMT" "PIB_TIME_FMT"\n",
          i_rank,
          times_elapsed[PIB_TOTAL],
          times_elapsed[PIB_REDISTRIBUTE],
          times_elapsed[PIB_COPIES],
          times_elapsed[PIB_EXCHANGE],
          times_elapsed[PIB_LOCAL],
          times_elapsed[PIB_PTB],
          times_elapsed[PIB_BTP]);
}













void
PDM_para_octree_points_inside_boxes2
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

  int n_recv_boxes = 0;
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


  int i_rank, n_rank;
  PDM_MPI_Comm_rank (octree->comm, &i_rank);
  PDM_MPI_Comm_size (octree->comm, &n_rank);

  if (VISU) {
    char filename[999];
    sprintf(filename, "octree_points_%3.3d.vtk", i_rank);
    _export_octree_points (filename, octree, 0);

    sprintf(filename, "octree_points_normalized_%3.3d.vtk", i_rank);
    _export_octree_points (filename, octree, 1);

    sprintf(filename, "boxes_%3.3d.vtk", i_rank);
    _export_boxes (filename, n_boxes, box_extents, box_g_num);
  }


  /* Multiple ranks */
  if (n_rank > 1) {
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
    /* Compute rank extents and build shared bounding-box tree */
    if (octree->used_rank_extents == NULL) {
      _compute_rank_extents (octree);
    }

    int *init_location = malloc (sizeof(int) * n_boxes * 3);
    for (int i = 0; i < n_boxes; i++) {
      init_location[3*i]     = i_rank;
      init_location[3*i + 1] = 0;
      init_location[3*i + 2] = i;
    }

    // Normalize using octree extents
    double _s[3] = {octree->global_extents[0],
                    octree->global_extents[1],
                    octree->global_extents[2]};
    double _d[3] = {octree->global_extents[3] - octree->global_extents[0],
                    octree->global_extents[4] - octree->global_extents[1],
                    octree->global_extents[5] - octree->global_extents[2]};
    double invd[3] = {1. / _d[0],
                      1. / _d[1],
                      1. / _d[2]};
    double *_box_extents = malloc (sizeof(double) * two_dim * n_boxes);
    for (int i = 0; i < 2 * n_boxes; i++) {
      for (int j = 0; j < dim; j++) {
        _box_extents[dim*i + j] = (box_extents[dim*i + j] - _s[j]) * invd[j];
      }
    }
    if (VISU) {
      char filename[999];
      sprintf(filename, "boxes_normalized_%3.3d.vtk", i_rank);
      _export_boxes (filename, n_boxes, _box_extents, box_g_num);
    }

    PDM_box_set_t *box_set = PDM_box_set_create (dim,
                                                 0,
                                                 0,
                                                 n_boxes,
                                                 box_g_num,
                                                 _box_extents,
                                                 1,
                                                 &n_boxes,
                                                 init_location,
                                                 octree->comm);

    /* Find which ranks intersect each box */
    int *box_rank_index = NULL;
    int *box_rank = NULL;
    PDM_box_tree_get_boxes_intersects (octree->bt_shared,
                                       box_set,
                                       &box_rank_index,
                                       &box_rank);

    int *send_count = PDM_array_zeros_int (n_rank);
    for (int i = 0; i < octree->n_used_rank; i++) {
      int rank = octree->used_rank[i];
      send_count[rank] = box_rank_index[rank+1] - box_rank_index[rank];
    }

    int *recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      octree->comm);

    int *send_shift = PDM_array_new_idx_from_sizes_int (send_count, n_rank);
    int *recv_shift = PDM_array_new_idx_from_sizes_int (recv_count, n_rank);

    n_recv_boxes = recv_shift[n_rank];

    const double      *box_set_extents = PDM_box_set_get_extents (box_set);
    const PDM_g_num_t *box_set_g_num   = PDM_box_set_get_g_num (box_set);

    /* Fill send buffers */
    PDM_g_num_t *send_box_g_num   = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    double      *send_box_extents = malloc (sizeof(double)      * send_shift[n_rank] * two_dim);
    for (int i = 0; i < octree->n_used_rank; i++) {
      int rank = octree->used_rank[i];
      int idx = send_shift[rank];
      for (int j = box_rank_index[i]; j < box_rank_index[i+1]; j++) {
        int ibox = box_rank[j];
        send_box_g_num[idx] = box_set_g_num[ibox];
        /*for (int k = 0; k < two_dim; k++) {
          send_box_extents[two_dim*idx + k] = _s[k] + _d[k]*box_set_extents[two_dim*ibox + k];
          }*/
        for (int k = 0; k < dim; k++) {
          send_box_extents[two_dim*idx + k] = _s[k] + _d[k]*box_set_extents[two_dim*ibox + k];
          send_box_extents[two_dim*idx + k + dim] = _s[k] + _d[k]*box_set_extents[two_dim*ibox + k + dim];
        }
        idx++;
      }
    }

    PDM_box_set_destroy (&box_set);
    free (init_location);
    free (box_rank_index);
    free (box_rank);


    /* Send boxes g_num buffer */
    _recv_box_g_num = malloc (sizeof(PDM_g_num_t) * recv_shift[n_rank]);
    recv_box_g_num = _recv_box_g_num;
    PDM_MPI_Alltoallv (send_box_g_num,  send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       _recv_box_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       octree->comm);

    /* Send boxes extents buffer */
    for (int i = 0; i < n_rank; i++) {
      send_shift[i+1] *= two_dim;
      recv_shift[i+1] *= two_dim;
      send_count[i]   *= two_dim;
      recv_count[i]   *= two_dim;
    }
    recv_box_extents = malloc (sizeof(double) * recv_shift[n_rank]);
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
    recv_box_extents = (double *) box_extents;
    recv_box_g_num   = box_g_num;
  }


  /***************************************
   * Intersect redistributed boxes with local octree
   ***************************************/
  if (VISU) {
    char filename[999];
    sprintf(filename, "recv_boxes1_%3.3d.vtk", i_rank);
    _export_boxes (filename, n_recv_boxes, recv_box_extents, recv_box_g_num);
  }

  /* Encode corners of redistributed boxes */
  PDM_morton_code_t *box_corners = malloc (sizeof(PDM_morton_code_t) * 2 * n_recv_boxes);
  _morton_encode_coords (dim,
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
  int *box_pts     = malloc (sizeof(int) * s_box_pts);
  int *box_pts_n   = malloc (sizeof(int) * n_recv_boxes);
  int *box_pts_idx = malloc (sizeof(int) * (n_recv_boxes + 1));
  box_pts_idx[0] = 0;

  size_t min_intersect = 100000000;
  size_t max_intersect = 0;
  size_t mean_intersect = 0;
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
                              octants->n_points,
                              0,
                              octants->n_nodes,
                              &n_intersect_nodes,
                              intersect_nodes);

    if (DEBUG) {
      printf("[%d]\tbox %d ("PDM_FMT_G_NUM") nodes:", i_rank, ibox, recv_box_g_num[ibox]);
      for (int j = 0; j < (int) n_intersect_nodes; j++) {
        printf(" %d", intersect_nodes[j]);
      }
      printf("\n");
    }

    // printf("n_intersect_nodes =%li \n", n_intersect_nodes);
    min_intersect = PDM_MIN(min_intersect, n_intersect_nodes);
    max_intersect = PDM_MAX(max_intersect, n_intersect_nodes);
    mean_intersect += n_intersect_nodes;
    //size_t new_max_size = box_pts_idx[ibox] + n_intersect_nodes * octree->points_in_leaf_max;
    size_t new_max_size = box_pts_idx[ibox];
    for (size_t i = 0; i < n_intersect_nodes; i++) {
      new_max_size += octants->n_points[intersect_nodes[i]];
    }

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
      printf("[%d]\tbox %d ("PDM_FMT_G_NUM") pts:", i_rank, ibox, recv_box_g_num[ibox]);
      for (int j = box_pts_idx[ibox]; j < box_pts_idx[ibox+1]; j++) {
        printf(" %d", box_pts[j]);
      }
      printf("\n");
    }

  } // Loop over redistributed boxes
  if (recv_box_extents != box_extents) free (recv_box_extents);
  free (box_corners);
  free (intersect_nodes);
  //log_debug("approx_size = %i | TRUE s_box_pts = %i \n", approx_size, s_box_pts);

  box_pts = realloc (box_pts, sizeof(int) * box_pts_idx[n_recv_boxes]);

  if (n_recv_boxes != 0) mean_intersect /= n_recv_boxes;
  //log_debug("min_intersect = %li | max_intersect = %li | mean_intersect = %li \n", min_intersect, max_intersect, mean_intersect);
  //log_debug("box_pts_idx[%i] = %li \n", n_recv_boxes, box_pts_idx[n_recv_boxes]);


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
  if (n_rank > 1) {
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
    if (_recv_box_g_num != box_g_num) free (_recv_box_g_num);

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

    for (int ibox = 0; ibox < n_recv_boxes; ibox++) {
      box_pts_n[ibox] *=  dim;
    }

    int *block_pts_in_box_n_coord = NULL;
    double *block_pts_in_box_coord = NULL;
    PDM_part_to_block_exch (ptb2,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &box_pts_n,
                            (void **) &box_pts_coord,
                            &block_pts_in_box_n_coord,
                            (void **) &block_pts_in_box_coord);
    free (box_pts_coord);
    free (box_pts_n);
    free (block_pts_in_box_n_coord);

    /* Fix partial block --> full block */
    int n_elt_block_full = PDM_part_to_block_n_elt_block_get (ptb1);
    PDM_g_num_t *block_elt_g_num_full = PDM_part_to_block_block_gnum_get (ptb1);

    int n_elt_block = PDM_part_to_block_n_elt_block_get (ptb2);
    PDM_g_num_t *block_elt_g_num = PDM_part_to_block_block_gnum_get (ptb2);

    int *block_pts_in_box_n_full = PDM_array_zeros_int (n_elt_block_full);
    int i1 = 0;
    for (int ielt = 0; ielt < n_elt_block; ielt++) {

      while (block_elt_g_num_full[i1] < block_elt_g_num[ielt]) {
        i1++;
      }

      block_pts_in_box_n_full[i1] = block_pts_in_box_n[ielt];
    }
    free (block_pts_in_box_n);


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
                            (void *) block_pts_in_box_n_full,
                            NULL,
                            (void **) &pts_in_box_n);

    *pts_in_box_idx = PDM_array_new_idx_from_sizes_int (pts_in_box_n, n_boxes);
    *pts_in_box_g_num = malloc (sizeof(PDM_g_num_t) * (*pts_in_box_idx)[n_boxes]);

    PDM_block_to_part_exch (btp,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            block_pts_in_box_n_full,
                            (void *) block_pts_in_box_g_num,
                            &pts_in_box_n,
                            (void **) pts_in_box_g_num);

    for (int ibox = 0; ibox < n_boxes; ibox++) {
      pts_in_box_n[ibox] *= dim;
    }
    for (int i = 0; i < n_elt_block_full; i++) {
      block_pts_in_box_n_full[i] *= dim;
    }

    *pts_in_box_coord = malloc (sizeof(double) * (*pts_in_box_idx)[n_boxes] * dim);

    PDM_block_to_part_exch (btp,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            block_pts_in_box_n_full,
                            (void *) block_pts_in_box_coord,
                            &pts_in_box_n,
                            (void **) pts_in_box_coord);

    PDM_part_to_block_free (ptb1);
    PDM_part_to_block_free (ptb2);
    PDM_block_to_part_free (btp);

    free (pts_in_box_n);
    free (block_pts_in_box_n_full);
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


/**
 *
 * \brief Copy octree data of some ranks into all ranks
 *
 * \param [in]   id                 Identifier
 * \param [in]   n_copied_ranks     Number of ranks to copy
 * \param [in]   copied_ranks       Array of ranks to copy
 *
 */

void
PDM_para_octree_copy_ranks
(
 const int  id,
 const int  n_copied_ranks,
 const int *copied_ranks
 )
{
  _octree_t *octree = _get_from_id (id);
  int dim = octree->dim;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (octree->comm, &i_rank);
  PDM_MPI_Comm_size (octree->comm, &n_rank);

  PDM_timer_hang_on (octree->timer);
  double b_t_elapsed = PDM_timer_elapsed (octree->timer);
  PDM_timer_resume (octree->timer);

  octree->n_copied_ranks = n_copied_ranks;

  octree->copied_ranks = malloc (sizeof(int) * n_copied_ranks);
  memcpy (octree->copied_ranks, copied_ranks, sizeof(int) * n_copied_ranks);

  octree->copied_octants     = malloc (sizeof(_l_octant_t *)       * n_copied_ranks);
  octree->n_copied_points    = malloc (sizeof(int)                 * n_copied_ranks);
  octree->copied_points      = malloc (sizeof(double *)            * n_copied_ranks);
  octree->copied_points_gnum = malloc (sizeof(PDM_g_num_t *)       * n_copied_ranks);
  octree->copied_points_code = malloc (sizeof(PDM_morton_code_t *) * n_copied_ranks);

  if (octree->explicit_nodes_to_build) {
    octree->copied_explicit_nodes = malloc (sizeof(_explicit_node_t *) * n_copied_ranks);
    octree->n_copied_explicit_nodes = malloc (sizeof(int) * n_copied_ranks);
  }

  int n[3];
  double      *pts_coord = NULL;
  PDM_g_num_t *pts_g_num = NULL;
  int         *codes     = NULL;
  int         *oct_n_pts = NULL;
  int         *ibuf      = NULL;
  double      *dbuf      = NULL;

  const int s_explicit_data = 16;

  for (int i = 0; i < n_copied_ranks; i++) {

    int rank = copied_ranks[i];

    if (rank == i_rank) {
      n[0] = octree->octants->n_nodes;
      n[1] = octree->n_points;
      n[2] = octree->n_explicit_nodes;
    }

    PDM_MPI_Bcast (n, 3, PDM_MPI_INT, rank, octree->comm);
    int n_copied_octants  = n[0];
    int n_copied_points   = n[1];
    int n_copied_explicit = n[2];

    int s_codes = (n_copied_octants + n_copied_points) * 4;
    codes     = malloc (sizeof(int) * s_codes);
    oct_n_pts = malloc (sizeof(int) * n_copied_octants * 2);
    if (octree->explicit_nodes_to_build) {
      ibuf = malloc (sizeof(int)    * n_copied_explicit * s_explicit_data);
      dbuf = malloc (sizeof(double) * n_copied_explicit * 6);
    }

    if (rank == i_rank) {
      octree->copied_octants[i]     = NULL;
      octree->n_copied_points[i]    = 0;
      octree->copied_points[i]      = NULL;
      octree->copied_points_gnum[i] = NULL;
      octree->copied_points_code[i] = NULL;

      for (int j = 0; j < n_copied_octants; j++) {
        codes[4*j] = (int) octree->octants->codes[j].L;
        for (int k = 0; k < 3; k++) {
          codes[4*j + k + 1] = (int) octree->octants->codes[j].X[k];
        }

        oct_n_pts[2*j]     = octree->octants->n_points[j];
        oct_n_pts[2*j + 1] = octree->octants->range[j];
      }

      int *_codes = codes + 4*n_copied_octants;
      for (int j = 0; j < n_copied_points; j++) {
        _codes[4*j] = octree->points_code[j].L;
        for (int k = 0; k < 3; k++) {
          _codes[4*j + k + 1] = octree->points_code[j].X[k];
        }
      }

      pts_coord = octree->points;
      pts_g_num = octree->points_gnum;

      if (octree->explicit_nodes_to_build) {
        octree->copied_explicit_nodes[i] = NULL;

        for (int j = 0; j < n_copied_explicit; j++) {
          ibuf[s_explicit_data*j] = (int) octree->explicit_nodes[j].code.L;
          for (int k = 0; k < 3; k++) {
            ibuf[s_explicit_data*j + 1 + k] = (int) octree->explicit_nodes[j].code.X[k];
          }
          ibuf[s_explicit_data*j + 4] = octree->explicit_nodes[j].n_points;
          ibuf[s_explicit_data*j + 5] = octree->explicit_nodes[j].range;
          ibuf[s_explicit_data*j + 6] = octree->explicit_nodes[j].ancestor_id;
          for (int k = 0; k < 8; k++) {
            ibuf[s_explicit_data*j + 7 + k] = octree->explicit_nodes[j].children_id[k];
          }
          ibuf[s_explicit_data*j + 15] = octree->explicit_nodes[j].leaf_id;

          for (int k = 0; k < 6; k++) {
            dbuf[6*j + k] = octree->explicit_nodes[j].pts_extents[k];
          }
        }
      }
    }

    else {
      octree->n_copied_points[i]    = n_copied_points;
      octree->copied_points[i]      = malloc (sizeof(double)      * dim * n_copied_points);
      octree->copied_points_gnum[i] = malloc (sizeof(PDM_g_num_t)       * n_copied_points);
      octree->copied_points_code[i] = malloc (sizeof(PDM_morton_code_t) * n_copied_points);

      pts_coord = octree->copied_points[i];
      pts_g_num = octree->copied_points_gnum[i];
    }

    PDM_MPI_Bcast (codes,     s_codes,               PDM_MPI_INT,        rank, octree->comm);
    PDM_MPI_Bcast (oct_n_pts, n_copied_octants * 2,  PDM_MPI_INT,        rank, octree->comm);
    PDM_MPI_Bcast (pts_coord, n_copied_points * dim, PDM_MPI_DOUBLE,     rank, octree->comm);
    PDM_MPI_Bcast (pts_g_num, n_copied_points,       PDM__PDM_MPI_G_NUM, rank, octree->comm);
    if (octree->explicit_nodes_to_build) {
      PDM_MPI_Bcast (ibuf, n_copied_explicit * s_explicit_data,
                     PDM_MPI_INT, rank, octree->comm);
      PDM_MPI_Bcast (dbuf, n_copied_explicit * 6,
                     PDM_MPI_DOUBLE, rank, octree->comm);
    }


    if (rank != i_rank) {

      octree->copied_octants[i] = malloc (sizeof(_l_octant_t));
      octree->copied_octants[i]->dim = dim;
      octree->copied_octants[i]->n_nodes = n_copied_octants;
      octree->copied_octants[i]->n_nodes_max = n_copied_octants;//?
      octree->copied_octants[i]->codes = malloc (sizeof(PDM_morton_code_t) * n_copied_octants);
      octree->copied_octants[i]->n_points = malloc (sizeof(int) * n_copied_octants);
      octree->copied_octants[i]->range    = malloc (sizeof(int) * n_copied_octants);
      octree->copied_octants[i]->neighbours    = NULL;
      octree->copied_octants[i]->neighbour_idx = NULL;

      for (int j = 0; j < n_copied_octants; j++) {
        octree->copied_octants[i]->codes[j].L = (PDM_morton_int_t) codes[4*j];
        for (int k = 0; k < 3; k++) {
          octree->copied_octants[i]->codes[j].X[k] = (PDM_morton_int_t) codes[4*j + k + 1];
        }

        octree->copied_octants[i]->n_points[j] = oct_n_pts[2*j];
        octree->copied_octants[i]->range[j]    = oct_n_pts[2*j + 1];
      }

      PDM_morton_code_t *_pts_codes = octree->copied_points_code[i];
      int *_codes = codes + 4*n_copied_octants;
      for (int j = 0; j < n_copied_points; j++) {
        _pts_codes[j].L = (PDM_morton_int_t) _codes[4*j];
        for (int k = 0; k < 3; k++) {
          _pts_codes[j].X[k] = (PDM_morton_int_t) _codes[4*j + k + 1];
        }
      }

      // TO DO: copy neighbours...


      /* Copy explicit nodes */
      if (octree->explicit_nodes_to_build) {
        octree->n_copied_explicit_nodes[i] = n_copied_explicit;
        octree->copied_explicit_nodes[i] = malloc (sizeof(_explicit_node_t) * n_copied_explicit);
        _explicit_node_t *cen = octree->copied_explicit_nodes[i];
        for (int j = 0; j < n_copied_explicit; j++) {
          cen[j].code.L = (PDM_morton_int_t) ibuf[s_explicit_data*j];
          for (int k = 0; k < 3; k++) {
            cen[j].code.X[k] = (PDM_morton_int_t) ibuf[s_explicit_data*j + 1 + k];
          }
          cen[j].n_points    = ibuf[s_explicit_data*j + 4];
          cen[j].range       = ibuf[s_explicit_data*j + 5];
          cen[j].ancestor_id = ibuf[s_explicit_data*j + 6];
          for (int k = 0; k < 8; k++) {
            cen[j].children_id[k] = ibuf[s_explicit_data*j + 7 + k];
          }
          cen[j].leaf_id = ibuf[s_explicit_data*j + 15];

          for (int k = 0; k < 6; k++) {
            cen[j].pts_extents[k] = dbuf[6*j + k];
          }
        }
      }

    }

    free (codes);
    free (oct_n_pts);
    if (ibuf != NULL) free (ibuf);
    if (dbuf != NULL) free (dbuf);
  }


  PDM_timer_hang_on (octree->timer);
  double e_t_elapsed = PDM_timer_elapsed (octree->timer);
  if (i_rank == 0) {
    printf("PDM_para_octree_copy_ranks: elapsed = %12.5es\n", e_t_elapsed - b_t_elapsed);
  }
  PDM_timer_resume (octree->timer);
}



/**
 *
 * \brief Free copied data in an octree structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_free_copies
(
 const int          id
 )
{
  _octree_t *octree = _get_from_id (id);

  if (octree->copied_ranks != NULL) {
    free (octree->copied_ranks);
    octree->copied_ranks = NULL;
  }

  if (octree->copied_octants != NULL) {
    for (int i = 0; i < octree->n_copied_ranks; i++) {
      if (octree->copied_octants[i] != NULL) {
        _octants_free (octree->copied_octants[i]);
      }
    }
    free (octree->copied_octants);
    octree->copied_octants = NULL;
  }

  if (octree->n_copied_points != NULL) {
    free (octree->n_copied_points);
    octree->n_copied_points = NULL;
  }

  if (octree->copied_points != NULL) {
    for (int i = 0; i < octree->n_copied_ranks; i++) {
      if (octree->copied_points[i] != NULL) {
        free (octree->copied_points[i]);
        octree->copied_points[i] = NULL;
      }
    }
    free (octree->copied_points);
    octree->copied_points = NULL;
  }

  if (octree->copied_points_gnum != NULL) {
    for (int i = 0; i < octree->n_copied_ranks; i++) {
      if (octree->copied_points_gnum[i] != NULL) {
        free (octree->copied_points_gnum[i]);
        octree->copied_points_gnum[i] = NULL;
      }
    }
    free (octree->copied_points_gnum);
    octree->copied_points_gnum = NULL;
  }

  if (octree->copied_points_code != NULL) {
    for (int i = 0; i < octree->n_copied_ranks; i++) {
      if (octree->copied_points_code[i] != NULL) {
        free (octree->copied_points_code[i]);
        octree->copied_points_code[i] = NULL;
      }
    }
    free (octree->copied_points_code);
    octree->copied_points_code = NULL;
  }

  if (octree->copied_explicit_nodes != NULL) {
    for (int i = 0; i < octree->n_copied_ranks; i++) {
      if (octree->copied_explicit_nodes[i] != NULL) {
        free (octree->copied_explicit_nodes[i]);
        octree->copied_explicit_nodes[i] = NULL;
      }
    }
    free (octree->copied_explicit_nodes);
    octree->copied_explicit_nodes = NULL;
  }

  if (octree->n_copied_explicit_nodes != NULL) {
    free (octree->n_copied_explicit_nodes);
    octree->n_copied_explicit_nodes = NULL;
  }

  octree->n_copied_ranks = 0;
}



#ifdef __cplusplus
}
#endif /* __cplusplus */
