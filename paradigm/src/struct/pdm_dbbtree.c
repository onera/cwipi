/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_morton.h"
#include "pdm_mpi.h"
#include "pdm_box.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_hash_tab.h"
#include "pdm_box_priv.h"
#include "pdm_box_tree.h"
#include "pdm_dbbtree.h"
#include "pdm_dbbtree_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_timer.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/**
 * \brief  Normalize
 */

static void
_normalize
(
 _PDM_dbbtree_t *dbbt,
 const double *pt_origin,
 double *pt_nomalized
 )
{
  for (int j = 0; j < dbbt->dim; j++) {
    pt_nomalized[j] = (pt_origin[j] - dbbt->s[j]) / dbbt->d[j];
  }
}

/**
 * \brief  Initialize box_tree statistics
 *
 * \param [inout]  bts  pointer to box tree statistics structure
 *
 */

static void
_init_bt_statistics
(
 _box_tree_stats_t  *bts
 )
{
  size_t i;

  assert(bts != NULL);

  bts->dim = 0;

  for (i = 0; i < 3; i++) {
    bts->depth[i] = 0;
    bts->n_leaves[i] = 0;
    bts->n_boxes[i] = 0;
    bts->n_threshold_leaves[i] = 0;
    bts->n_leaf_boxes[i] = 0;
    bts->mem_used[i] = 0;
    bts->mem_required[i] = 0;
  }
}


/**
 * \brief Update box-tree statistics.
 *
 * For most fields, we replace previous values with the current ones.
 *
 * For memory required, we are interested in the maximum values over time
 * (i.e. algorthm steps); this is the case even for the minimal memory
 * required, we is thus the time maximum of the rank minimum.
 *
 * \param [inout]   bts   Pointer to box tree statistics structure
 * \param [inout]   bt    Pointer to box tree structure
 *
 */

static void
_update_bt_statistics
(
 _box_tree_stats_t     *bts,
 const PDM_box_tree_t  *bt
 )
{
  int dim;
  size_t i;
  size_t mem_required[3];

  assert(bts != NULL);

  dim = PDM_box_tree_get_stats (bt,
                                bts->depth,
                                bts->n_leaves,
                                bts->n_boxes,
                                bts->n_threshold_leaves,
                                bts->n_leaf_boxes,
                                bts->mem_used,
                                mem_required);

  bts->dim = dim;

  for (i = 0; i < 3; i++)
    bts->mem_required[i] = _MAX(bts->mem_required[i], mem_required[i]);
}


/**
 * \brief Distribute bounding boxes over the ranks according to a Morton encoding
 * index. Try to get a well-balanced distribution and spatially coherent.
 *
 * \param [inout]  n     <-> pointer to neighborhood management structure
 * \param [inout]  boxes <-> box set to redistribute
 *
 */

static void
_redistribute_boxes
(
 _PDM_dbbtree_t      *dbbt
 )
{

  /* Sanity checks */

  assert (dbbt != NULL);

  PDM_box_tree_t  *coarse_tree = PDM_box_tree_create (dbbt->maxTreeDepthCoarse,
                                                      dbbt->maxBoxesLeafCoarse,
                                                      dbbt->maxBoxRatioCoarse);

  /* Build a tree and associate boxes */

  PDM_box_tree_set_boxes (coarse_tree,
                          dbbt->boxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics(&(dbbt->btsCoarse), coarse_tree);

  if (1 == 0) {
    PDM_printf ("-- dump stats\n");

    PDM_box_tree_dump_statistics(coarse_tree);

    PDM_printf ("-- fin dump stats\n");

    PDM_printf ("-- dump \n");

    PDM_box_tree_dump(coarse_tree);

    PDM_printf ("-- fin dump\n");
  }

  /*
   * Compute an index based on Morton encoding to ensure a good distribution
   * of bounding boxes among the ranks.
   */

  PDM_box_distrib_t  *distrib = PDM_box_tree_get_distrib (coarse_tree, dbbt->boxes);

  PDM_box_tree_destroy (&coarse_tree);

  if (1 == 0) {
    PDM_box_distrib_dump_statistics (distrib, dbbt->comm);
  }

  /* Define a new distribution of boxes according to the Morton
     encoding index */

  if (1 == 0) {
    PDM_printf("affichage 1\n");
    PDM_box_set_dump( dbbt->boxes,1);
    PDM_printf("fin affichage 1\n");
  }

  PDM_box_set_redistribute (distrib, dbbt->boxes);

  if (1 == 0) {
    PDM_printf("affichage 2\n");
    PDM_box_set_dump( dbbt->boxes,1);
    PDM_printf("fin affichage 2\n");
  }

  /* Delete intermediate structures */

  PDM_box_distrib_destroy (&distrib);

}



static void _export_point_cloud
(
 char         *filename,
 int           n_part,
 const int    *n_pts,
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


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Return an intialized \ref PDM_dbbtree_t structure
 *
 * This function returns an initialized \ref PDM_dbbtree_t structure
 *
 * \param [in]  comm             Associated communicator
 * \param [in]  dim              boxes dimension
 * \param [in]  global_extents   Globals of elements to storage into the tree
 *                               (automatic computation if NULL)
 *
 * \return      A new initialized \ref PDM_dbbtree_t structure
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_create
(
 PDM_MPI_Comm  comm,
 int           dim,
 double       *global_extents
 )
{
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) malloc(sizeof(_PDM_dbbtree_t));

  _dbbt->comm                 = comm;
  _dbbt->rankComm             = PDM_MPI_COMM_NULL;
  _dbbt->dim                  = dim;

  _dbbt->maxTreeDepth         = 30;
  _dbbt->maxBoxRatio          = 10.;
  _dbbt->maxBoxesLeaf         = 30;

  _dbbt->maxTreeDepthShared   = 10;
  _dbbt->maxBoxRatioShared    =  6;
  _dbbt->maxBoxesLeafShared   =  5;

  _dbbt->maxTreeDepthCoarse   = 10;
  _dbbt->maxBoxRatioCoarse    =  4;
  _dbbt->maxBoxesLeafCoarse   = 30;

  _dbbt->rankBoxes            = NULL;
  _dbbt->btShared             = NULL;

  _dbbt->nUsedRank            = 0;
  _dbbt->usedRank             = NULL;

  _dbbt->boxes                = NULL;
  _dbbt->btLoc                = NULL;

  _dbbt->global_extents       = NULL;

  if (global_extents != NULL) {
    _dbbt->global_extents = malloc (sizeof(double) * dim * 2);
    memcpy(_dbbt->global_extents, global_extents, sizeof(double) * dim * 2);
    for (int j = 0; j < dim; j++) {
      _dbbt->s[j] = _dbbt->global_extents[j];
      _dbbt->d[j] = _dbbt->global_extents[j+dim] - _dbbt->global_extents[j];
    }
  }

  else {
    for (int j = 0; j < dim; j++) {
      _dbbt->s[j] = 0.;
      _dbbt->d[j] = 1.;
    }
  }

  _init_bt_statistics (&(_dbbt->btsShared));
  _init_bt_statistics (&(_dbbt->btsLoc));
  _init_bt_statistics (&(_dbbt->btsCoarse));

  return (PDM_dbbtree_t *) _dbbt;
}



/**
 * \brief Free a \ref PDM_dbbtree_t structure
 *
 * \param [in]  dbbt   Pointer to a distributed bounding box tree
 *
 * \return      NULL
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_free
(
 PDM_dbbtree_t     *dbbt
 )
{
  if (dbbt != NULL) {
    _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

    PDM_box_set_destroy (&_dbbt->rankBoxes);
    //PDM_box_set_destroy (&_dbbt->boxes);

    if (_dbbt->global_extents != NULL) {
      free (_dbbt->global_extents);
    }

    free (_dbbt->usedRank);
    if (_dbbt->rankComm != PDM_MPI_COMM_NULL) {
      PDM_MPI_Comm_free (&(_dbbt->rankComm));
    }

    PDM_box_tree_destroy (&_dbbt->btShared);
    PDM_box_tree_destroy (&_dbbt->btLoc);

    free (_dbbt);
  }

  return NULL;
}



/**
 * \brief Assign a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * This function assigns a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  n_part    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 *
 * \return associated \ref PDM_box_set_t structure distributed according to
 * the tree location
 *
 */

PDM_box_set_t *
PDM_dbbtree_boxes_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const double      **extents,
 const PDM_g_num_t **gNum
 )
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  const int nInfoLocation = 3;
  const int sExtents = _dbbt->dim * 2;

  int nEltsProc = 0;
  for (int i = 0; i < n_part; i++) {
    nEltsProc += nElts[i];
  }

  PDM_g_num_t *_boxGnum = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nEltsProc);
  double *_extents      = (double *) malloc (sizeof(double) * nEltsProc * sExtents);
  int *_initLocation   = (int *) malloc (sizeof(int) * nEltsProc * nInfoLocation);

  int idx = 0;
  int idx1 = 0;
  int idx2 = 0;

  for (int i = 0; i < n_part; i++) {

    for (int j = 0; j < nElts[i]; j++) {
      _boxGnum[idx++] = gNum[i][j];

      for (int k = 0; k < sExtents; k++) {
        _extents[idx1++] = extents[i][sExtents*j+k];
      }

      _initLocation[idx2++] = myRank;
      _initLocation[idx2++] = i;
      _initLocation[idx2++] = j;
    }
  }

  /*
   * Redistribute boxes of mesh A
   */

  if (0 == 1) {

    PDM_printf ("nEltsProc : %d\n", nEltsProc);

    PDM_printf ("_boxGnum :");
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, _boxGnum[i]);
    }
    PDM_printf ("\n");

    PDM_printf ("_extents m2 :\n");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM":", _boxGnum[i]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");

    PDM_printf ("_initLocation :");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }

  int ind_norm = 1;
  if (_dbbt->global_extents != NULL) {
    ind_norm = 0;
    for (int i = 0; i < nEltsProc; i++) {
      _normalize (_dbbt,
                  _extents+2*i*_dbbt->dim,
                  _extents+2*i*_dbbt->dim);
      _normalize (_dbbt,
                  _extents+(2*i+1)*_dbbt->dim,
                  _extents+(2*i+1)*_dbbt->dim);
    }
  }

  _dbbt->boxes = PDM_box_set_create(3,
                                    ind_norm,  // No normalization to preserve initial extents
                                    0,  // No projection to preserve initial extents
                                    nEltsProc,
                                    _boxGnum,
                                    _extents,
                                    n_part,
                                    nElts,
                                    _initLocation,
                                    _dbbt->comm);

  if (_dbbt->global_extents == NULL) {
    memcpy (_dbbt->d, _dbbt->boxes->d, sizeof(double) * 3);
    memcpy (_dbbt->s, _dbbt->boxes->s, sizeof(double) * 3);
  }

  free (_boxGnum);
  free (_extents);
  free (_initLocation);

  if (lComm > 1) {
    _redistribute_boxes(_dbbt);

    /*
     * Compute processus extents
     */

    int nBoxes = PDM_box_set_get_size (_dbbt->boxes);
    const double *extents2 = PDM_box_set_get_extents (_dbbt->boxes);

    double gExtents[sExtents];
    for (int i = 0; i < _dbbt->dim; i++) {
      gExtents[i]   =  DBL_MAX;
      gExtents[_dbbt->dim+i] = -DBL_MAX;
    }


    for (int i = 0; i < nBoxes; i++) {
      for (int k1 = 0; k1 < _dbbt->dim; k1++) {
        gExtents[k1]   = _MIN (gExtents[k1], extents2[sExtents * i + k1]);
        gExtents[_dbbt->dim+k1] = _MAX (gExtents[_dbbt->dim+k1], extents2[sExtents * i
                                                                          + _dbbt->dim + k1]);
      }
    }

    /*
     * Exchange extents to build shared boxes
     */

    int *allNBoxes = (int *) malloc (sizeof(int) * lComm);
    PDM_MPI_Allgather (&nBoxes, 1, PDM_MPI_INT,
                       allNBoxes, 1, PDM_MPI_INT,
                       _dbbt->comm);

    int nUsedRank = 0;
    for (int i = 0; i < lComm; i++) {
      if (allNBoxes[i] > 0) {
        nUsedRank += 1;
      }
    }


    double *allGExtents = (double *) malloc (sizeof(double) * sExtents * lComm);
    PDM_MPI_Allgather (gExtents, sExtents, PDM__PDM_MPI_REAL,
                       allGExtents, sExtents, PDM__PDM_MPI_REAL,
                       _dbbt->comm);

    int *numProc = (int *) malloc (sizeof(int *) * nUsedRank);

    _dbbt->usedRank = numProc;
    _dbbt->nUsedRank = nUsedRank;

    PDM_g_num_t *gNumProc = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nUsedRank);

    idx = 0;
    for (int i = 0; i < lComm; i++) {
      if (allNBoxes[i] > 0) {
        gNumProc[idx] = idx;
        numProc[idx] = i;
        for (int j = 0; j < sExtents; j++) {
          allGExtents[idx*sExtents + j] = allGExtents[i*sExtents + j];
        }
        idx += 1;
      }
    }

    free (allNBoxes);

    allGExtents = (double *) realloc (allGExtents, sizeof(double) * sExtents * nUsedRank);

    int *initLocation_proc = PDM_array_zeros_int(3 * nUsedRank);

    //TODO: Faire un PDM_box_set et PDM_box_tree_create sequentiel ! Le comm split a u n 1 proc ici : pas terrible

    PDM_MPI_Comm_split(_dbbt->comm, myRank, 0, &(_dbbt->rankComm));

    _dbbt->rankBoxes = PDM_box_set_create(3,
                                          0,  // No normalization to preserve initial extents
                                          0,  // No projection to preserve initial extents
                                          nUsedRank,
                                          gNumProc,
                                          allGExtents,
                                          1,
                                          &nUsedRank,
                                          initLocation_proc,
                                          _dbbt->rankComm);

    memcpy (_dbbt->rankBoxes->d, _dbbt->d, sizeof(double) * 3);
    memcpy (_dbbt->rankBoxes->s, _dbbt->s, sizeof(double) * 3);

    _dbbt->btShared = PDM_box_tree_create (_dbbt->maxTreeDepthShared,
                                           _dbbt->maxBoxesLeafShared,
                                           _dbbt->maxBoxRatioShared);

    /* Build a tree and associate boxes */

    PDM_box_tree_set_boxes (_dbbt->btShared,
                            _dbbt->rankBoxes,
                            PDM_BOX_TREE_ASYNC_LEVEL);

    _update_bt_statistics(&(_dbbt->btsShared), _dbbt->btShared);

    free (allGExtents);
    free (gNumProc);
    free (initLocation_proc);

  }

  /*
   * Build local bt
   */

  _dbbt->btLoc = PDM_box_tree_create (_dbbt->maxTreeDepth,
                                      _dbbt->maxBoxesLeaf,
                                      _dbbt->maxBoxRatio);

  /* Build a tree and associate boxes */

  PDM_box_tree_set_boxes (_dbbt->btLoc,
                          _dbbt->boxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics(&(_dbbt->btsLoc), _dbbt->btLoc);

  return _dbbt->boxes;

}



/**
 * \brief Assign boxes to intersect to the tree.
 *
 * This function  assigns boxes to intersect to the tree.
 *
 * \param [in]  dbbt       Pointer to a distributed bounding box tree
 * \param [in]  n_part     Number of partitions
 * \param [in]  nElts      Number of elements of each partition
 * \param [in]  extents    Extents of each element of each partition
 * \param [in]  gNum       Global number of each element of each partition
 * \param [out] box_index  Pointer to the index array on associated tree bounding boxeq
 * \param [out] box_g_num  Pointer to the list of intersecting bounding boxes
 *
 * \return associated \ref PDM_box_set_t structure distributed according
 * to the tree intersection
 *
 */

PDM_box_set_t *
PDM_dbbtree_intersect_boxes_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const double      **extents,
 const PDM_g_num_t **gNum,
 int                *box_index[],
 int                *box_l_num[]
 )
{

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  const int nInfoLocation = 3;
  const int sExtents = _dbbt->dim * 2;

  int nEltsProc = 0;
  for (int i = 0; i < n_part; i++) {
    nEltsProc += nElts[i];
  }

  PDM_g_num_t *_boxGnum = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nEltsProc);
  double *_extents     = (double *) malloc (sizeof(double) * nEltsProc * sExtents);
  int *_initLocation   = (int *) malloc (sizeof(int) * nEltsProc * nInfoLocation);

  int idx = 0;
  int idx1 = 0;
  int idx2 = 0;

  for (int i = 0; i < n_part; i++) {

    for (int j = 0; j < nElts[i]; j++) {
      _boxGnum[idx++] = gNum[i][j];

      for (int k = 0; k < sExtents; k++) {
        _extents[idx1++] = extents[i][sExtents*j+k];
      }

      _initLocation[idx2++] = myRank;
      _initLocation[idx2++] = i;
      _initLocation[idx2++] = j;
    }

  }
  if (1 == 0) {
    PDM_printf ("_extents m1 :\n");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM":", _boxGnum[i]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");

    PDM_printf ("_initLocation :");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }

  for (int i = 0; i < nEltsProc; i++) {
    _normalize (_dbbt,
                _extents+2*i*_dbbt->boxes->dim,
                _extents+2*i*_dbbt->boxes->dim);
    _normalize (_dbbt,
                _extents+(2*i+1)*_dbbt->boxes->dim,
                _extents+(2*i+1)*_dbbt->boxes->dim);
  }

  if (1 == 0) {

    PDM_printf ("nEltsProc : %d\n", nEltsProc);

    PDM_printf ("_boxGnum :");
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, _boxGnum[i]);
    }
    PDM_printf ("\n");

    PDM_printf ("_extents m2 :\n");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM":", _boxGnum[i]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");

    PDM_printf ("_initLocation :");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }

  PDM_box_set_t  *boxes = PDM_box_set_create (3,
                                              0,  // No normalization to preserve initial extents
                                              0,  // No projection to preserve initial extents
                                              nEltsProc,
                                              _boxGnum,
                                              _extents,
                                              n_part,
                                              nElts,
                                              _initLocation,
                                              _dbbt->comm);


  /*
   * Intersection boxes whith shared tree
   */

  if (_dbbt->btShared != NULL) {

    PDM_box_tree_get_boxes_intersects (_dbbt->btShared,
                                       boxes,
                                       box_index,
                                       box_l_num);

    int nUsedRank = PDM_box_set_get_size (_dbbt->rankBoxes);
    const int *usedRanks = _dbbt->usedRank;

    /*
     * Distribute boxes on intersection ranks
     */

    if (1 == 0){
      PDM_printf ("box_l_num_shared : \n");
      for (int i = 0; i < nUsedRank; i++) {
        printf("[%d] : %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
               i, _dbbt->rankBoxes->local_boxes->extents[6*i]
               , _dbbt->rankBoxes->local_boxes->extents[6*i+1]
               , _dbbt->rankBoxes->local_boxes->extents[6*i+2]
               , _dbbt->rankBoxes->local_boxes->extents[6*i+3]
               , _dbbt->rankBoxes->local_boxes->extents[6*i+4]
               , _dbbt->rankBoxes->local_boxes->extents[6*i+5]);
        printf("[%d] : ",i);
        for (int j = (*box_index)[i]; j < (*box_index)[i+1]; j++) {
          PDM_printf (" "PDM_FMT_G_NUM"",  boxes->local_boxes->g_num[(*box_l_num)[j]]);
        }
        PDM_printf ("\n");
      }
    }

    PDM_box_distrib_t  *distrib = NULL;
    distrib = PDM_box_distrib_create (boxes->local_boxes->n_boxes,//*
                                      boxes->n_g_boxes,
                                      1, // Don't use in this case
                                      boxes->comm);

    PDM_array_reset_int(distrib->index, lComm + 1, 0);

    for (int i = 0; i < nUsedRank; i++) {
      distrib->index[usedRanks[i]+1] = (*box_index)[i+1] - (*box_index)[i];
    }

    for (int i = 0; i < lComm; i++) {
      distrib->index[i+1] = distrib->index[i+1] + distrib->index[i];
    }

    distrib->list = (int *) malloc (sizeof(int) * distrib->index[lComm]);

    for (int i = 0; i < distrib->index[lComm]; i++) {
      distrib->list[i] = (*box_l_num)[i];
    }

    /*
     * Redistribute boxes on intersecting ranks
     */

    PDM_box_set_redistribute (distrib, boxes);

    if (1 == 0) {
      printf ("Boxes B apres redistribution : %d\n", boxes->local_boxes->n_boxes);
      for (int i = 0; i < boxes->local_boxes->n_boxes; i++) {
        printf (" "PDM_FMT_G_NUM"", boxes->local_boxes->g_num[i]);
      }
      printf("\n");
    }

    /*
     * Free
     */

    PDM_box_distrib_destroy (&distrib);

    free (*box_l_num);
    free (*box_index);

  }

  free (_boxGnum);
  free (_extents);
  free (_initLocation);

  /*
   * Intersection boxes whith local tree
   */

  //PDM_box_tree_dump_statistics(_dbbt->btLoc);

  PDM_box_tree_get_boxes_intersects (_dbbt->btLoc,
                                     boxes,
                                     box_index,
                                     box_l_num);

  /*
   * Sort boxes and remove double boxes
   */

  int nBoxesA = _dbbt->boxes->local_boxes->n_boxes;//*

  int *newIndex = PDM_array_zeros_int(nBoxesA + 1);

  int *_box_index = *box_index;
  int *_box_l_num = *box_l_num;

  idx = 0;
  for (int i = 0; i < nBoxesA; i++) {

    int *ideb = _box_l_num + _box_index[i];
    int length = _box_index[i+1] - _box_index[i];

    PDM_sort_int (ideb, NULL, length);

    int pre = -1;
    for (int k = 0; k < length; k++) {
      if (pre != ideb[k]) {
        (*box_l_num)[idx++] = ideb[k];
        newIndex[i+1] += 1;
        pre = ideb[k];
      }
    }
  }

  for (int i = 0; i < nBoxesA; i++) {
    newIndex[i+1] += newIndex[i];
  }

  free (*box_index);
  *box_index = newIndex;

  *box_l_num = (int *) realloc (*box_l_num, sizeof (int) * newIndex[nBoxesA]);

  if (1 == 0) {
    printf ("Intersections : %d\n", boxes->local_boxes->n_boxes);
    for (int i = 0; i < nBoxesA; i++) {
      printf ("A elt "PDM_FMT_G_NUM" :", _dbbt->boxes->local_boxes->g_num[i]);
      for (int j = newIndex[i]; j < newIndex[i+1]; j++) {
        printf (" "PDM_FMT_G_NUM"", boxes->local_boxes->g_num[(*box_l_num)[j]]);
      }
      printf("\n");
    }
  }



  return boxes;

}



/**
 *
 * \brief Get the boxes closer than the upper bound distance
 *
 *   \param [in]  bt                 Pointer to box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts                Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  pts_g_num          Point global ids
 *   \param [in]  upper_bound_dist2  Upper bound of the square of the distance (size = \ref n_pts)
 *   \param [out] box_index          Index of boxes (size = \ref n_pts + 1)
 *   \param [out] box_g_num          Global ids of boxes (size = \ref i_boxes[\ref n_pts])
 *
 */

void
PDM_dbbtree_closest_upper_bound_dist_boxes_get
(
 PDM_dbbtree_t   *dbbt,
 const int        n_pts,
 double           pts[],
 PDM_g_num_t      pts_g_num[],
 double           upper_bound_dist2[],
 int             *box_index[],
 PDM_g_num_t     *box_g_num[]
 )
{
  PDM_UNUSED(pts_g_num);
  /*
   * RANK DATA COPY PARAMETERS
   */
  const double RANK_COPY_threshold  = 1.2;  // factor of the mean nb of requests
  const double RANK_COPY_max_copies = 0.15; // factor of the total nb of processes

  /*
   * Initialization
   */
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  double *_pts = pts;
  if (_dbbt->global_extents != NULL) {
    _pts = malloc (sizeof(double) * n_pts * 3);

    for (int i = 0; i < n_pts; i++) {
      _normalize (_dbbt,
                  pts + 3*i,
                  _pts + 3*i);
    }
  }

  int     n_pts_local            = n_pts;
  double *pts_local              = _pts;
  double *upper_bound_dist_local = upper_bound_dist2;


  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  /*
   * Determine for each point the list of involved processes
   */
  int *n_send_pts = NULL;
  int *n_recv_pts = NULL;

  int *box_index_tmp = NULL;
  int *box_l_num_tmp = NULL;

  int *n_pts_rank = NULL;
  int *n_pts_send = NULL;
  int *n_pts_recv = NULL;
  int n_pts_recv_total = 0;

  int *i_pts_rank = NULL;
  int *i_pts_send = NULL;
  int *i_pts_recv = NULL;

  double *pts_rank              = NULL;
  double *upper_bound_dist_rank = NULL;
  double *pts_recv              = NULL;
  double *upper_bound_dist_recv = NULL;

  int n_copied_ranks = 0;
  int *copied_ranks  = NULL;
  int *rank_copy_num = NULL;

  const int *usedRanks = _dbbt->usedRank;

  const int idebug = 0;

  int i1 = 0, i2 = 0, i3 = 0;
  int i_rank = 0;
  if (_dbbt->btShared != NULL) {

    if (idebug) {
      printf ("  **** deb PDM_box_tree_closest_upper_bound_dist_boxes_get shared _pts : %d\n", n_pts);
    }

    PDM_box_tree_closest_upper_bound_dist_boxes_get (_dbbt->btShared,
                                                     n_pts,
                                                     pts,
                                                     upper_bound_dist2,
                                                     &box_index_tmp,
                                                     &box_l_num_tmp);

    if (idebug) {
      printf ("  **** fin PDM_box_tree_closest_upper_bound_dist_boxes_get shared n_pts : %d\n", n_pts);
      for (int i = 0; i < n_pts; i++) {
        printf ("%d : (%12.5e %12.5e %12.5e) %12.5e\n", i,
                pts[3*i], pts[3*i+1], pts[3*i+2],
                upper_bound_dist2[i]);
        printf ("  boxes %d :" , box_index_tmp[i+1] - box_index_tmp[i]);
        for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
          printf (" %d", box_l_num_tmp[j]);
        }
        printf ("\n");
      }
    }

    /*
     * Count (provisional) nb of points to send to each process
     */
    n_send_pts = PDM_array_zeros_int(lComm);
    n_recv_pts = malloc (sizeof(int) * lComm);


    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        n_send_pts[usedRanks[box_l_num_tmp[j]]]++;
      }
    }

    PDM_MPI_Alltoall (n_send_pts, 1, PDM_MPI_INT,
                      n_recv_pts, 1, PDM_MPI_INT,
                      _dbbt->comm);
    free(n_send_pts);

    /*
     * Prepare copies
     */
    // total nb of requests received by current process
    int local_sum_nrecv = 0;
    for (int i = 0; i < lComm; i++) {
      local_sum_nrecv += n_recv_pts[i];
    }
    free(n_recv_pts);

    int *n_requests = malloc (lComm * sizeof(int));
    PDM_MPI_Allgather (&local_sum_nrecv, 1, PDM_MPI_INT,
                       n_requests,       1, PDM_MPI_INT,
                       _dbbt->comm);

    // mean nb of requests
    int mean_n_requests = 0;
    for (int i = 0; i < lComm; i++) {
      mean_n_requests += n_requests[i];
    }
    mean_n_requests /= lComm;

    /* sort the ranks in ascending order of
     * the total nb of points they are supposed to receive */
    int *order = malloc (lComm * sizeof(int));
    for (int i = 0; i < lComm; i ++) {
      order[i] = i;
    }

    PDM_sort_int (n_requests, order, lComm);

    // identify ranks to be copied
    double threshold_n_req = RANK_COPY_threshold*mean_n_requests;
    int max_copied_ranks   = (int) _MAX (1, RANK_COPY_max_copies*lComm);

    n_copied_ranks = 0;
    copied_ranks = malloc (max_copied_ranks * sizeof(int));

    for (int i = 0; i < max_copied_ranks; i++) {
      i_rank = lComm - 1 - i;

      if ( n_requests[i_rank] > threshold_n_req ) {
        copied_ranks[n_copied_ranks++] = order[i_rank];
      } else {
        break;
      }
    }
    free(order);
    free(n_requests);

    //------------->>>
    if ( idebug && myRank == 0 ) {
      if ( n_copied_ranks == 0 ) {
        printf("n_copied_ranks = 0\n");
      } else {
        printf("copied rank(s) = ");
        for (int i = 0; i < n_copied_ranks; i++) {
          printf("%d ", copied_ranks[i]);
        }
        printf("\n");
      }
    }
    //<<<-------------

    /*
     * Copy the data of selected ranks
     */
    rank_copy_num = (int *) malloc (sizeof(int) * lComm);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                &n_copied_ranks,
                                copied_ranks,
                                rank_copy_num);
    /* rank_copy_num[_dbbt->btLoc->copied_ranks[i]] (def)= i*/
    free(copied_ranks);


    /*
     * Distribution of points...
     *    ..._local --> search in local box tree
     *    ..._rank  --> search in copied box trees
     *    ..._send  --> search in distant box trees (send to other processes)
     */
    n_pts_local = 0;

    n_pts_rank = PDM_array_zeros_int(n_copied_ranks);

    n_pts_send = PDM_array_zeros_int(lComm);
    n_pts_recv = malloc (sizeof(int) * lComm);


    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        i_rank = usedRanks[box_l_num_tmp[j]];
        if ( i_rank == myRank ) {
          // ---> search in btLoc->local_data of current process
          n_pts_local++;
        } else if ( rank_copy_num[i_rank] >= 0 ) {
          // ---> search in btLoc->rank_data[rank_copy_num[i_rank]] of current process
          n_pts_rank[rank_copy_num[i_rank]]++;
        } else {
          // ---> search in btLoc->local_data of process with rank i_rank
          n_pts_send[i_rank]++;
        }
      }
    }

    PDM_MPI_Alltoall (n_pts_send, 1, PDM_MPI_INT,
                      n_pts_recv, 1, PDM_MPI_INT,
                      _dbbt->comm);

    i_pts_rank = PDM_array_new_idx_from_sizes_int(n_pts_rank, n_copied_ranks);

    i_pts_send = malloc (sizeof(int) * (lComm+1));
    i_pts_recv = malloc (sizeof(int) * (lComm+1));
    i_pts_send[0] = 0;
    i_pts_recv[0] = 0;
    for (int i = 0; i < lComm; i++) {
      i_pts_send[i+1] = i_pts_send[i] + 4 * n_pts_send[i];
      i_pts_recv[i+1] = i_pts_recv[i] + 4 * n_pts_recv[i];
      n_pts_recv[i] *= 4;
    }

    pts_local              = malloc (sizeof(double) * n_pts_local*3);
    upper_bound_dist_local = malloc (sizeof(double) * n_pts_local);

    pts_rank              = malloc (sizeof(double) * i_pts_rank[n_copied_ranks]*3);
    upper_bound_dist_rank = malloc (sizeof(double) * i_pts_rank[n_copied_ranks]);

    double *data_send = malloc (sizeof(double) * i_pts_send[lComm]);
    double *data_recv = malloc (sizeof(double) * i_pts_recv[lComm]);

    PDM_array_reset_int(n_pts_send, lComm, 0);
    PDM_array_reset_int(n_pts_rank, n_copied_ranks, 0);

    i1 = 0; i2 = 0; i3 = 0;
    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        i_rank = usedRanks[box_l_num_tmp[j]];
        if ( i_rank == myRank ) {
          // pts_local, upper_bound_dist_local (local points, local box tree)
          upper_bound_dist_local[i1] = upper_bound_dist2[i];
          pts_local[3*i1]            = _pts[3*i];
          pts_local[3*i1+1]          = _pts[3*i+1];
          pts_local[3*i1+2]          = _pts[3*i+2];
          i1++;
        } else if ( rank_copy_num[i_rank] >= 0 ) {
          // pts_rank, upper_bound_dist_rank (local points, distant (copied) box trees)
          int j_rank = rank_copy_num[i_rank];
          i2 = i_pts_rank[j_rank] + n_pts_rank[j_rank];
          upper_bound_dist_rank[i2] = upper_bound_dist2[i];
          pts_rank[3*i2]            = _pts[3*i];
          pts_rank[3*i2+1]          = _pts[3*i+1];
          pts_rank[3*i2+2]          = _pts[3*i+2];
          n_pts_rank[j_rank]++;
        } else {
          // data_send (local points, distant (not copied) box trees)
          i3 = i_pts_send[i_rank] + 4*n_pts_send[i_rank];
          data_send[i3++] = _pts[3*i];
          data_send[i3++] = _pts[3*i+1];
          data_send[i3++] = _pts[3*i+2];
          data_send[i3++] = upper_bound_dist2[i];

          n_pts_send[i_rank]++;
        }
      }
    }

    for (int i = 0; i < lComm; i++) {
      n_pts_send[i] *= 4;
    }


    // Send points to search in distant (not copied) box trees
    PDM_MPI_Alltoallv (data_send, n_pts_send, i_pts_send, PDM_MPI_DOUBLE,
                       data_recv, n_pts_recv, i_pts_recv, PDM_MPI_DOUBLE,
                       _dbbt->comm);
    free(data_send);

    n_pts_recv_total = i_pts_recv[lComm] / 4;

    pts_recv              = (double *) malloc (sizeof(double) * n_pts_recv_total * 3);
    upper_bound_dist_recv = (double *) malloc (sizeof(double) * n_pts_recv_total);

    for (int i = 0; i < n_pts_recv_total; i++) {
      for (int j = 0; j < 3; j++) {
        pts_recv[3*i+j] = data_recv[4*i+j];
      }
      upper_bound_dist_recv[i] = data_recv[4*i+3];
    }
    free(data_recv);
  }
  if (_pts != pts && _pts != pts_local) {
    free (_pts);
  }


  // Determine candidate boxes in local box tree (local points)
  int *box_index_local;
  int *box_l_num_local;
  PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                      -1, // search in local box tree
                                                      n_pts_local,
                                                      pts_local,
                                                      upper_bound_dist_local,
                                                      &box_index_local,
                                                      &box_l_num_local,
                                                      _dbbt->d);
  if (pts_local != pts) {
    free(pts_local);
  }
  if (upper_bound_dist_local != upper_bound_dist2) {
    free(upper_bound_dist_local);
  }

  // conversion local --> global numbering
  const PDM_g_num_t *gnum_boxes_local = PDM_box_set_get_g_num (_dbbt->boxes);
  PDM_g_num_t *box_g_num_local = malloc(sizeof(PDM_g_num_t) * box_index_local[n_pts_local]);

  for (int i = 0; i < box_index_local[n_pts_local]; i++) {
    box_g_num_local[i] = gnum_boxes_local[box_l_num_local[i]];
  }
  free(box_l_num_local);


  if (_dbbt->btShared == NULL) {

    *box_index = box_index_local;
    *box_g_num = box_g_num_local;

  } else {
    // Determine candidate boxes in copied box trees (local points)
    int **box_index_rank;
    int **box_l_num_rank;

    box_index_rank = (int **) malloc (sizeof(int *) * n_copied_ranks);
    box_l_num_rank = (int **) malloc (sizeof(int *) * n_copied_ranks);

    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      int n_pts_copied_rank = i_pts_rank[i_copied_rank+1] - i_pts_rank[i_copied_rank];
      double *pts_copied_rank = pts_rank + 3*i_pts_rank[i_copied_rank];
      double *upper_bound_dist_copied_rank = upper_bound_dist_rank + i_pts_rank[i_copied_rank];
      PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                          i_copied_rank,
                                                          n_pts_copied_rank,
                                                          pts_copied_rank,
                                                          upper_bound_dist_copied_rank,
                                                          &(box_index_rank[i_copied_rank]),
                                                          &(box_l_num_rank[i_copied_rank]),
                                                          _dbbt->d);
    }

    free(pts_rank);
    free(upper_bound_dist_rank);

    // conversion local --> global numbering for each copied rank
    PDM_g_num_t **box_g_num_rank = malloc(sizeof(PDM_g_num_t *) * n_copied_ranks);

    PDM_g_num_t *gnum_boxes_rank = NULL;
    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      gnum_boxes_rank = PDM_box_set_get_rank_boxes_g_num (_dbbt->boxes,
                                                          i_copied_rank);
      box_g_num_rank[i_copied_rank] = malloc(sizeof(PDM_g_num_t) * box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]]);
      for (int i = 0; i < box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]]; i++) {
        box_g_num_rank[i_copied_rank][i] = gnum_boxes_rank[box_l_num_rank[i_copied_rank][i]];
      }

      free(box_l_num_rank[i_copied_rank]);
    }
    free(box_l_num_rank);



    // Determine candidate boxes in local box tree (received points)
    int *n_pts_send2 = NULL;
    int *i_pts_send2 = NULL;
    int *n_box_l_num_per_pts = NULL;
    PDM_g_num_t *box_g_num_per_pts = NULL;


    int *box_index_recv = NULL;
    int *box_l_num_recv = NULL;

    PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                        -1, // search in local box tree
                                                        n_pts_recv_total,
                                                        pts_recv,
                                                        upper_bound_dist_recv,
                                                        &box_index_recv,
                                                        &box_l_num_recv,
                                                        _dbbt->d);
    free(pts_recv);
    free(upper_bound_dist_recv);

    /*
     * Send back results for distant points to original processes:
     *     - nb of boxes for each point
     *     - global numbering of these boxes
     */

    int *n_box_l_num_recv = malloc (sizeof(int) * n_pts_recv_total);

    for (int i = 0; i < n_pts_recv_total; i++) {
      n_box_l_num_recv[i] = box_index_recv[i+1] - box_index_recv[i];
    }

    for (int i = 0; i < lComm; i++) {
      i_pts_send[i+1] = i_pts_send[i+1]/4;
      i_pts_recv[i+1] = i_pts_recv[i+1]/4;
      n_pts_send[i]   = n_pts_send[i]/4;
      n_pts_recv[i]   = n_pts_recv[i]/4;
    }

    n_box_l_num_per_pts = malloc (sizeof(int) * i_pts_send[lComm]);

    PDM_MPI_Alltoallv (n_box_l_num_recv,    n_pts_recv, i_pts_recv, PDM_MPI_INT,
                       n_box_l_num_per_pts, n_pts_send, i_pts_send, PDM_MPI_INT,
                       _dbbt->comm);

    n_pts_send2 = malloc (sizeof(int) * lComm);
    i_pts_send2 = malloc (sizeof(int) * (lComm+1));

    int *n_pts_recv2 = malloc (sizeof(int) * lComm);
    int *i_pts_recv2 = malloc (sizeof(int) * (lComm+1));

    for (int i = 0; i < lComm; i++) {
      n_pts_send2[i] = 0;
      n_pts_recv2[i] = 0;
    }

    for (int i = 0; i < lComm; i++) {
      for (int j = i_pts_recv[i]; j < i_pts_recv[i+1]; j++) {
        n_pts_recv2[i] += n_box_l_num_recv[j];
      }
      for (int j = i_pts_send[i]; j < i_pts_send[i+1]; j++) {
        n_pts_send2[i] += n_box_l_num_per_pts[j];
      }
    }

    free(n_box_l_num_recv);

    i_pts_send2[0] = 0;
    i_pts_recv2[0] = 0;
    for (int i = 0; i < lComm; i++) {
      i_pts_send2[i+1] = i_pts_send2[i] + n_pts_send2[i];
      i_pts_recv2[i+1] = i_pts_recv2[i] + n_pts_recv2[i];
    }


    // Conversion local --> global numbering
    PDM_g_num_t *box_g_num_recv = malloc(sizeof(PDM_g_num_t) * box_index_recv[n_pts_recv_total]);

    for (int i = 0; i < box_index_recv[n_pts_recv_total]; i++) {
      box_g_num_recv[i] = gnum_boxes_local[box_l_num_recv[i]];
    }

    box_g_num_per_pts = malloc(sizeof(PDM_g_num_t) * i_pts_send2[lComm]);
    PDM_MPI_Alltoallv (box_g_num_recv,    n_pts_recv2, i_pts_recv2, PDM__PDM_MPI_G_NUM,
                       box_g_num_per_pts, n_pts_send2, i_pts_send2, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);

    free(box_index_recv);
    free(box_l_num_recv);
    free(box_g_num_recv);

    free(n_pts_recv);
    free(i_pts_recv);
    free(n_pts_recv2);
    free(i_pts_recv2);




    /*
     * Merge all results and resolve duplicates
     */
    *box_index = malloc(sizeof(int) * (n_pts + 1));

    int max_n_box_g_num = box_index_local[n_pts_local];
    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      max_n_box_g_num += box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]];
    }
    max_n_box_g_num += i_pts_send2[lComm];

    *box_g_num = malloc(sizeof(PDM_g_num_t) * max_n_box_g_num);


    int *rank_index = (int *) malloc(sizeof(int) * (n_pts+1));
    memcpy(rank_index, box_index_tmp, sizeof(int) * (n_pts+1));

    for (int i = 0; i < lComm; i++) {
      n_pts_send[i] = 0;
      n_pts_send2[i] = 0;
    }

    PDM_array_reset_int(n_pts_rank, n_copied_ranks, 0);



    int keyMax = 3 * n_pts;
    int key = 0;
    int found = 0;
    PDM_hash_tab_t *ht = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                              &keyMax);

    PDM_g_num_t i_box = 0;

    box_index_tmp[0] = 0;
    int idx = 0;
    i1 = 0; i2 = 0; i3 = 0;
    if ( 1 ) {
      for (int i = 0; i < n_pts; i++) { // loop over local points
        box_index_tmp[i+1] = box_index_tmp[i];
        for (int j = rank_index[i]; j < rank_index[i+1]; j++) { // loop over procs to which the current point was sent
          i_rank = usedRanks[box_l_num_tmp[j]]; // i_rank = rank j-th proc to which the current point was sent

          if ( i_rank == myRank ) { // boxes local to current process
            for (int k = box_index_local[i1]; k < box_index_local[i1+1]; k++) {
              i_box = box_g_num_local[k];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            i1++;

          } else if ( rank_copy_num[i_rank] >= 0 ) { // distant boxes copied in current process
            int j_rank = rank_copy_num[i_rank];
            i2 = n_pts_rank[j_rank];

            for (int k = box_index_rank[j_rank][i2]; k < box_index_rank[j_rank][i2+1]; k++) {
              i_box = box_g_num_rank[j_rank][k];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            n_pts_rank[j_rank]++;

          } else { // distant boxes (not copied)
            i3 = n_pts_send[i_rank];
            int i4 = i_pts_send2[i_rank] + n_pts_send2[i_rank];
            int i5 = i_pts_send[i_rank] + i3;
            for (int k = 0; k < n_box_l_num_per_pts[i5]; k++) {
              i_box = box_g_num_per_pts[i4++];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            n_pts_send2[i_rank] += n_box_l_num_per_pts[i5];
            n_pts_send[i_rank]++;

          }
        }
        PDM_hash_tab_purge (ht, PDM_FALSE);
      }
    }
    PDM_hash_tab_free (ht);

    free (n_box_l_num_per_pts);


    if ( *box_index != NULL ) {
      free (*box_index);
    }
    *box_index = box_index_tmp;


    *box_g_num = realloc (*box_g_num, sizeof(PDM_g_num_t) * box_index_tmp[n_pts]);

    /*
     * Deallocate stuff
     */
    free(box_l_num_tmp);
    free(box_g_num_per_pts);

    free(box_index_local);
    free(box_g_num_local);

    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      free(box_index_rank[i_copied_rank]);
      free(box_g_num_rank[i_copied_rank]);
    }
    free(box_index_rank);
    free(box_g_num_rank);


    free(n_pts_send);
    free(i_pts_send);
    free(n_pts_send2);
    free(i_pts_send2);

    free(i_pts_rank);
    free(n_pts_rank);

    free(rank_index);
    free(rank_copy_num);
  }

}



/**
 *
 * \brief Get the boxes closer than the upper bound distance (Asynchronous)
 *
 *   \param [in]  bt                 Pointer to box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts                Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  pts_g_num          Point global ids
 *   \param [in]  upper_bound_dist2  Upper bound of the square of the distance (size = \ref n_pts)
 *   \param [out] box_index          Index of boxes (size = \ref n_pts + 1)
 *   \param [out] box_g_num          Global ids of boxes (size = \ref i_boxes[\ref n_pts])
 *
 */

void
PDM_dbbtree_closest_upper_bound_dist_boxes_get_async
(
 PDM_dbbtree_t   *dbbt,
 const int        n_pts,
 double           pts[],
 PDM_g_num_t      pts_g_num[],
 double           upper_bound_dist2[],
 int             *box_index[],
 PDM_g_num_t     *box_g_num[]
 )
{
  PDM_UNUSED(pts_g_num);
  /*
   * RANK DATA COPY PARAMETERS
   */
  const double RANK_COPY_threshold  = 1.2;  // factor of the mean nb of requests
  const double RANK_COPY_max_copies = 0.15; // factor of the total nb of processes

  /*
   * Initialization
   */
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  double *_pts = pts;
  if (_dbbt->global_extents != NULL) {
    _pts = malloc (sizeof(double) * n_pts * 3);

    for (int i = 0; i < n_pts; i++) {
      _normalize (_dbbt,
                  pts + 3*i,
                  _pts + 3*i);
    }
  }

  int     n_pts_local            = n_pts;
  double *pts_local              = _pts;
  double *upper_bound_dist_local = upper_bound_dist2;


  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  /*
   * Determine for each point the list of involved processes
   */
  int *n_send_pts = NULL;
  int *n_recv_pts = NULL;

  int *box_index_tmp = NULL;
  int *box_l_num_tmp = NULL;

  int *n_pts_rank = NULL;
  int *n_pts_send = NULL;
  int *n_pts_recv = NULL;
  int n_pts_recv_total = 0;

  int *i_pts_rank = NULL;
  int *i_pts_send = NULL;
  int *i_pts_recv = NULL;

  double *pts_rank              = NULL;
  double *upper_bound_dist_rank = NULL;
  double *pts_recv              = NULL;
  double *upper_bound_dist_recv = NULL;

  int n_copied_ranks = 0;
  int *copied_ranks  = NULL;
  int *rank_copy_num = NULL;

  const int *usedRanks = _dbbt->usedRank;

  const int idebug = 0;

  int i1 = 0, i2 = 0, i3 = 0;
  int i_rank = 0;
  PDM_MPI_Request request;
  double *data_send = NULL;
  double *data_recv = NULL;
  if (_dbbt->btShared != NULL) {

    if (idebug) {
      printf ("  **** deb PDM_box_tree_closest_upper_bound_dist_boxes_get shared _pts : %d\n", n_pts);
    }

    PDM_box_tree_closest_upper_bound_dist_boxes_get (_dbbt->btShared,
                                                     n_pts,
                                                     pts,
                                                     upper_bound_dist2,
                                                     &box_index_tmp,
                                                     &box_l_num_tmp);

    if (idebug) {
      printf ("  **** fin PDM_box_tree_closest_upper_bound_dist_boxes_get shared n_pts : %d\n", n_pts);
      for (int i = 0; i < n_pts; i++) {
        printf ("%d : (%12.5e %12.5e %12.5e) %12.5e\n", i,
                pts[3*i], pts[3*i+1], pts[3*i+2],
                upper_bound_dist2[i]);
        printf ("  boxes %d :" , box_index_tmp[i+1] - box_index_tmp[i]);
        for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
          printf (" %d", box_l_num_tmp[j]);
        }
        printf ("\n");
      }
    }

    /*
     * Count (provisional) nb of points to send to each process
     */
    n_send_pts = PDM_array_zeros_int(lComm);
    n_recv_pts = malloc (sizeof(int) * lComm);


    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        n_send_pts[usedRanks[box_l_num_tmp[j]]]++;
      }
    }

    PDM_MPI_Alltoall (n_send_pts, 1, PDM_MPI_INT,
                      n_recv_pts, 1, PDM_MPI_INT,
                      _dbbt->comm);
    free(n_send_pts);

    /*
     * Prepare copies
     */
    // total nb of requests received by current process
    int local_sum_nrecv = 0;
    for (int i = 0; i < lComm; i++) {
      local_sum_nrecv += n_recv_pts[i];
    }
    free(n_recv_pts);

    int *n_requests = malloc (lComm * sizeof(int));
    PDM_MPI_Allgather (&local_sum_nrecv, 1, PDM_MPI_INT,
                       n_requests,       1, PDM_MPI_INT,
                       _dbbt->comm);

    // mean nb of requests
    int mean_n_requests = 0;
    for (int i = 0; i < lComm; i++) {
      mean_n_requests += n_requests[i];
    }
    mean_n_requests /= lComm;

    /* sort the ranks in ascending order of
     * the total nb of points they are supposed to receive */
    int *order = malloc (lComm * sizeof(int));
    for (int i = 0; i < lComm; i ++) {
      order[i] = i;
    }

    PDM_sort_int (n_requests, order, lComm);

    // identify ranks to be copied
    double threshold_n_req = RANK_COPY_threshold*mean_n_requests;
    int max_copied_ranks   = (int) _MAX (1, RANK_COPY_max_copies*lComm);

    n_copied_ranks = 0;
    copied_ranks = malloc (max_copied_ranks * sizeof(int));

    for (int i = 0; i < max_copied_ranks; i++) {
      i_rank = lComm - 1 - i;

      if ( n_requests[i_rank] > threshold_n_req ) {
        copied_ranks[n_copied_ranks++] = order[i_rank];
      } else {
        break;
      }
    }
    free(order);
    free(n_requests);

    //------------->>>
    if ( myRank == 0 ) {
      if ( n_copied_ranks == 0 ) {
        printf("n_copied_ranks = 0\n");
      } else {
        printf("copied rank(s) = ");
        for (int i = 0; i < n_copied_ranks; i++) {
          printf("%d ", copied_ranks[i]);
        }
        printf("\n");
      }
    }
    //<<<-------------

    /*
     * Copy the data of selected ranks
     */
    rank_copy_num = (int *) malloc (sizeof(int) * lComm);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                &n_copied_ranks,
                                copied_ranks,
                                rank_copy_num);
    /* rank_copy_num[_dbbt->btLoc->copied_ranks[i]] (def)= i*/
    free(copied_ranks);


    /*
     * Distribution of points...
     *    ..._local --> search in local box tree
     *    ..._rank  --> search in copied box trees
     *    ..._send  --> search in distant box trees (send to other processes)
     */
    n_pts_local = 0;

    n_pts_rank = PDM_array_zeros_int(n_copied_ranks);

    n_pts_send = PDM_array_zeros_int(lComm);
    n_pts_recv = malloc (sizeof(int) * lComm);


    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        i_rank = usedRanks[box_l_num_tmp[j]];
        if ( i_rank == myRank ) {
          // ---> search in btLoc->local_data of current process
          n_pts_local++;
        } else if ( rank_copy_num[i_rank] >= 0 ) {
          // ---> search in btLoc->rank_data[rank_copy_num[i_rank]] of current process
          n_pts_rank[rank_copy_num[i_rank]]++;
        } else {
          // ---> search in btLoc->local_data of process with rank i_rank
          n_pts_send[i_rank]++;
        }
      }
    }

    PDM_MPI_Alltoall (n_pts_send, 1, PDM_MPI_INT,
                      n_pts_recv, 1, PDM_MPI_INT,
                      _dbbt->comm);

    i_pts_rank = PDM_array_new_idx_from_sizes_int(n_pts_rank, n_copied_ranks);

    i_pts_send = malloc (sizeof(int) * (lComm+1));
    i_pts_recv = malloc (sizeof(int) * (lComm+1));
    i_pts_send[0] = 0;
    i_pts_recv[0] = 0;
    for (int i = 0; i < lComm; i++) {
      i_pts_send[i+1] = i_pts_send[i] + 4 * n_pts_send[i];
      i_pts_recv[i+1] = i_pts_recv[i] + 4 * n_pts_recv[i];
      n_pts_recv[i] *= 4;
    }

    pts_local              = malloc (sizeof(double) * n_pts_local*3);
    upper_bound_dist_local = malloc (sizeof(double) * n_pts_local);

    pts_rank              = malloc (sizeof(double) * i_pts_rank[n_copied_ranks]*3);
    upper_bound_dist_rank = malloc (sizeof(double) * i_pts_rank[n_copied_ranks]);

    data_send = malloc (sizeof(double) * i_pts_send[lComm]);
    data_recv = malloc (sizeof(double) * i_pts_recv[lComm]);

    PDM_array_reset_int(n_pts_send, lComm, 0);
    PDM_array_reset_int(n_pts_rank, n_copied_ranks, 0);

    i1 = 0; i2 = 0; i3 = 0;
    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        i_rank = usedRanks[box_l_num_tmp[j]];
        if ( i_rank == myRank ) {
          // pts_local, upper_bound_dist_local (local points, local box tree)
          upper_bound_dist_local[i1] = upper_bound_dist2[i];
          pts_local[3*i1]            = _pts[3*i];
          pts_local[3*i1+1]          = _pts[3*i+1];
          pts_local[3*i1+2]          = _pts[3*i+2];
          i1++;
        } else if ( rank_copy_num[i_rank] >= 0 ) {
          // pts_rank, upper_bound_dist_rank (local points, distant (copied) box trees)
          int j_rank = rank_copy_num[i_rank];
          i2 = i_pts_rank[j_rank] + n_pts_rank[j_rank];
          upper_bound_dist_rank[i2] = upper_bound_dist2[i];
          pts_rank[3*i2]            = _pts[3*i];
          pts_rank[3*i2+1]          = _pts[3*i+1];
          pts_rank[3*i2+2]          = _pts[3*i+2];
          n_pts_rank[j_rank]++;
        } else {
          // data_send (local points, distant (not copied) box trees)
          i3 = i_pts_send[i_rank] + 4*n_pts_send[i_rank];
          data_send[i3++] = _pts[3*i];
          data_send[i3++] = _pts[3*i+1];
          data_send[i3++] = _pts[3*i+2];
          data_send[i3++] = upper_bound_dist2[i];

          n_pts_send[i_rank]++;
        }
      }
    }

    for (int i = 0; i < lComm; i++) {
      n_pts_send[i] *= 4;
    }


    // Send points to search in distant (not copied) box trees
    // PDM_MPI_Alltoallv (data_send, n_pts_send, i_pts_send, PDM_MPI_DOUBLE,
    //                    data_recv, n_pts_recv, i_pts_recv, PDM_MPI_DOUBLE,
    //                    _dbbt->comm);
    // free(data_send);

    // n_pts_recv_total = i_pts_recv[lComm] / 4;

    // pts_recv              = (double *) malloc (sizeof(double) * n_pts_recv_total * 3);
    // upper_bound_dist_recv = (double *) malloc (sizeof(double) * n_pts_recv_total);

    // for (int i = 0; i < n_pts_recv_total; i++) {
    //   for (int j = 0; j < 3; j++) {
    //     pts_recv[3*i+j] = data_recv[4*i+j];
    //   }
    //   upper_bound_dist_recv[i] = data_recv[4*i+3];
    // }
    // free(data_recv);

    /* Asynchrone */
    PDM_MPI_Ialltoallv (data_send, n_pts_send, i_pts_send, PDM_MPI_DOUBLE,
                        data_recv, n_pts_recv, i_pts_recv, PDM_MPI_DOUBLE,
                        _dbbt->comm, &request);



  }
  if (_pts != pts && _pts != pts_local) {
    free (_pts);
  }


  // Determine candidate boxes in local box tree (local points)
  int *box_index_local;
  int *box_l_num_local;

  // Determine candidate boxes in local box tree (received points)
  int *n_pts_send2 = NULL;
  int *i_pts_send2 = NULL;
  int *n_box_l_num_per_pts = NULL;
  PDM_g_num_t *box_g_num_per_pts = NULL;
  int *box_index_recv = NULL;
  int *box_l_num_recv = NULL;
  int *n_box_l_num_recv = NULL;

  if (_dbbt->btShared != NULL) {
    PDM_MPI_Wait(&request);
    free(data_send);

    n_pts_recv_total = i_pts_recv[lComm] / 4;

    pts_recv              = (double *) malloc (sizeof(double) * n_pts_recv_total * 3);
    upper_bound_dist_recv = (double *) malloc (sizeof(double) * n_pts_recv_total);

    for (int i = 0; i < n_pts_recv_total; i++) {
      for (int j = 0; j < 3; j++) {
        pts_recv[3*i+j] = data_recv[4*i+j];
      }
      upper_bound_dist_recv[i] = data_recv[4*i+3];
    }
    free(data_recv);

    PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                        -1, // search in local box tree
                                                        n_pts_recv_total,
                                                        pts_recv,
                                                        upper_bound_dist_recv,
                                                        &box_index_recv,
                                                        &box_l_num_recv,
                                                        _dbbt->d);
    free(pts_recv);
    free(upper_bound_dist_recv);

    /*
     * Send back results for distant points to original processes:
     *     - nb of boxes for each point
     *     - global numbering of these boxes
     */

    n_box_l_num_recv = malloc (sizeof(int) * n_pts_recv_total);

    // log_debug("n_pts_recv_total %i \n", n_pts_recv_total);
    for (int i = 0; i < n_pts_recv_total; i++) {
      n_box_l_num_recv[i] = box_index_recv[i+1] - box_index_recv[i];
    }

    for (int i = 0; i < lComm; i++) {
      i_pts_send[i+1] = i_pts_send[i+1]/4;
      i_pts_recv[i+1] = i_pts_recv[i+1]/4;
      n_pts_send[i]   = n_pts_send[i]/4;
      n_pts_recv[i]   = n_pts_recv[i]/4;
    }

    n_box_l_num_per_pts = malloc (sizeof(int) * i_pts_send[lComm]);

    // double t1i = PDM_MPI_Wtime();
    // PDM_MPI_Alltoallv (n_box_l_num_recv,    n_pts_recv, i_pts_recv, PDM_MPI_INT,
    //                    n_box_l_num_per_pts, n_pts_send, i_pts_send, PDM_MPI_INT,
    //                    _dbbt->comm);
    PDM_MPI_Ialltoallv (n_box_l_num_recv,    n_pts_recv, i_pts_recv, PDM_MPI_INT,
                        n_box_l_num_per_pts, n_pts_send, i_pts_send, PDM_MPI_INT,
                        _dbbt->comm, &request);
    // double dti = PDM_MPI_Wtime() - t1i;
    // log_debug("dt first all_to_all %12.5e \n", dti);
    // PDM_log_trace_array_int(n_pts_recv, lComm, "n_pts_recv :: ");
    // PDM_log_trace_array_int(n_pts_send, lComm, "n_pts_send :: ");
  }

  // double t1 = PDM_MPI_Wtime();
  // log_debug("n_pts_local %i \n", n_pts_local);
  PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                      -1, // search in local box tree
                                                      n_pts_local,
                                                      pts_local,
                                                      upper_bound_dist_local,
                                                      &box_index_local,
                                                      &box_l_num_local,
                                                      _dbbt->d);
  if (pts_local != pts) {
    free(pts_local);
  }
  if (upper_bound_dist_local != upper_bound_dist2) {
    free(upper_bound_dist_local);
  }

  // conversion local --> global numbering
  const PDM_g_num_t *gnum_boxes_local = PDM_box_set_get_g_num (_dbbt->boxes);
  PDM_g_num_t *box_g_num_local = malloc(sizeof(PDM_g_num_t) * box_index_local[n_pts_local]);

  for (int i = 0; i < box_index_local[n_pts_local]; i++) {
    box_g_num_local[i] = gnum_boxes_local[box_l_num_local[i]];
  }
  free(box_l_num_local);

  if (_dbbt->btShared == NULL) {

    *box_index = box_index_local;
    *box_g_num = box_g_num_local;

  } else {
    // Determine candidate boxes in copied box trees (local points)
    int **box_index_rank;
    int **box_l_num_rank;

    box_index_rank = (int **) malloc (sizeof(int *) * n_copied_ranks);
    box_l_num_rank = (int **) malloc (sizeof(int *) * n_copied_ranks);

    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      int n_pts_copied_rank = i_pts_rank[i_copied_rank+1] - i_pts_rank[i_copied_rank];
      double *pts_copied_rank = pts_rank + 3*i_pts_rank[i_copied_rank];
      double *upper_bound_dist_copied_rank = upper_bound_dist_rank + i_pts_rank[i_copied_rank];
      PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                          i_copied_rank,
                                                          n_pts_copied_rank,
                                                          pts_copied_rank,
                                                          upper_bound_dist_copied_rank,
                                                          &(box_index_rank[i_copied_rank]),
                                                          &(box_l_num_rank[i_copied_rank]),
                                                          _dbbt->d);
    }

    free(pts_rank);
    free(upper_bound_dist_rank);

    // conversion local --> global numbering for each copied rank
    PDM_g_num_t **box_g_num_rank = malloc(sizeof(PDM_g_num_t *) * n_copied_ranks);

    PDM_g_num_t *gnum_boxes_rank = NULL;
    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      gnum_boxes_rank = PDM_box_set_get_rank_boxes_g_num (_dbbt->boxes,
                                                          i_copied_rank);
      box_g_num_rank[i_copied_rank] = malloc(sizeof(PDM_g_num_t) * box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]]);
      for (int i = 0; i < box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]]; i++) {
        box_g_num_rank[i_copied_rank][i] = gnum_boxes_rank[box_l_num_rank[i_copied_rank][i]];
      }

      free(box_l_num_rank[i_copied_rank]);
    }
    free(box_l_num_rank);
    // double dt = PDM_MPI_Wtime() - t1;
    // log_debug("dt first management copie %12.5e \n", dt);

    if (_dbbt->btShared != NULL) {
      PDM_MPI_Wait(&request);
    }



    n_pts_send2 = malloc (sizeof(int) * lComm);
    i_pts_send2 = malloc (sizeof(int) * (lComm+1));

    int *n_pts_recv2 = malloc (sizeof(int) * lComm);
    int *i_pts_recv2 = malloc (sizeof(int) * (lComm+1));

    for (int i = 0; i < lComm; i++) {
      n_pts_send2[i] = 0;
      n_pts_recv2[i] = 0;
    }

    for (int i = 0; i < lComm; i++) {
      for (int j = i_pts_recv[i]; j < i_pts_recv[i+1]; j++) {
        n_pts_recv2[i] += n_box_l_num_recv[j];
      }
      for (int j = i_pts_send[i]; j < i_pts_send[i+1]; j++) {
        n_pts_send2[i] += n_box_l_num_per_pts[j];
      }
    }

    free(n_box_l_num_recv);

    i_pts_send2[0] = 0;
    i_pts_recv2[0] = 0;
    for (int i = 0; i < lComm; i++) {
      i_pts_send2[i+1] = i_pts_send2[i] + n_pts_send2[i];
      i_pts_recv2[i+1] = i_pts_recv2[i] + n_pts_recv2[i];
    }


    // Conversion local --> global numbering
    PDM_g_num_t *box_g_num_recv = malloc(sizeof(PDM_g_num_t) * box_index_recv[n_pts_recv_total]);

    for (int i = 0; i < box_index_recv[n_pts_recv_total]; i++) {
      box_g_num_recv[i] = gnum_boxes_local[box_l_num_recv[i]];
    }

    box_g_num_per_pts = malloc(sizeof(PDM_g_num_t) * i_pts_send2[lComm]);
    // t1 = PDM_MPI_Wtime();
    PDM_MPI_Alltoallv (box_g_num_recv,    n_pts_recv2, i_pts_recv2, PDM__PDM_MPI_G_NUM,
                       box_g_num_per_pts, n_pts_send2, i_pts_send2, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    // dt = PDM_MPI_Wtime() - t1;
    // log_debug("dt second all_to_all %12.5e \n", dt);
    // PDM_log_trace_array_int(n_pts_recv2, lComm, "n_pts_recv2 :: ");
    // PDM_log_trace_array_int(n_pts_send2, lComm, "n_pts_send2 :: ");

    free(box_index_recv);
    free(box_l_num_recv);
    free(box_g_num_recv);

    free(n_pts_recv);
    free(i_pts_recv);
    free(n_pts_recv2);
    free(i_pts_recv2);




    /*
     * Merge all results and resolve duplicates
     */
    *box_index = malloc(sizeof(int) * (n_pts + 1));

    int max_n_box_g_num = box_index_local[n_pts_local];
    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      max_n_box_g_num += box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]];
    }
    max_n_box_g_num += i_pts_send2[lComm];

    *box_g_num = malloc(sizeof(PDM_g_num_t) * max_n_box_g_num);


    int *rank_index = (int *) malloc(sizeof(int) * (n_pts+1));
    memcpy(rank_index, box_index_tmp, sizeof(int) * (n_pts+1));

    for (int i = 0; i < lComm; i++) {
      n_pts_send[i] = 0;
      n_pts_send2[i] = 0;
    }

    PDM_array_reset_int(n_pts_rank, n_copied_ranks, 0);



    int keyMax = 3 * n_pts;
    int key = 0;
    int found = 0;
    PDM_hash_tab_t *ht = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                              &keyMax);

    PDM_g_num_t i_box = 0;

    box_index_tmp[0] = 0;
    int idx = 0;
    i1 = 0; i2 = 0; i3 = 0;
    if ( 1 ) {
      for (int i = 0; i < n_pts; i++) { // loop over local points
        box_index_tmp[i+1] = box_index_tmp[i];
        for (int j = rank_index[i]; j < rank_index[i+1]; j++) { // loop over procs to which the current point was sent
          i_rank = usedRanks[box_l_num_tmp[j]]; // i_rank = rank j-th proc to which the current point was sent

          if ( i_rank == myRank ) { // boxes local to current process
            for (int k = box_index_local[i1]; k < box_index_local[i1+1]; k++) {
              i_box = box_g_num_local[k];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            i1++;

          } else if ( rank_copy_num[i_rank] >= 0 ) { // distant boxes copied in current process
            int j_rank = rank_copy_num[i_rank];
            i2 = n_pts_rank[j_rank];

            for (int k = box_index_rank[j_rank][i2]; k < box_index_rank[j_rank][i2+1]; k++) {
              i_box = box_g_num_rank[j_rank][k];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            n_pts_rank[j_rank]++;

          } else { // distant boxes (not copied)
            i3 = n_pts_send[i_rank];
            int i4 = i_pts_send2[i_rank] + n_pts_send2[i_rank];
            int i5 = i_pts_send[i_rank] + i3;
            for (int k = 0; k < n_box_l_num_per_pts[i5]; k++) {
              i_box = box_g_num_per_pts[i4++];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            n_pts_send2[i_rank] += n_box_l_num_per_pts[i5];
            n_pts_send[i_rank]++;

          }
        }
        PDM_hash_tab_purge (ht, PDM_FALSE);
      }
    }
    PDM_hash_tab_free (ht);

    free (n_box_l_num_per_pts);


    if ( *box_index != NULL ) {
      free (*box_index);
    }
    *box_index = box_index_tmp;


    *box_g_num = realloc (*box_g_num, sizeof(PDM_g_num_t) * box_index_tmp[n_pts]);

    /*
     * Deallocate stuff
     */
    free(box_l_num_tmp);
    free(box_g_num_per_pts);

    free(box_index_local);
    free(box_g_num_local);

    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      free(box_index_rank[i_copied_rank]);
      free(box_g_num_rank[i_copied_rank]);
    }
    free(box_index_rank);
    free(box_g_num_rank);


    free(n_pts_send);
    free(i_pts_send);
    free(n_pts_send2);
    free(i_pts_send2);

    free(i_pts_rank);
    free(n_pts_rank);

    free(rank_index);
    free(rank_copy_num);
  }

}



/**
 *
 * \brief Export boxes and points to ASCII VTK format for debugging purposes
 *
 * A set of files (normalized/unnormalized points and boxes) is written
 * for each MPI rank.
 *
 *   \param [in]  filename_pattern  File name pattern (prefix)
 *   \param [in]  dbbt              Pointer to distributed box tree structure
 *   \param [in]  n_pts             Number of points
 *   \param [in]  pts_g_num         Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord         Point coordinates (size = 3 * \ref n_pts)
 *
 */

void
PDM_dbbtree_points_debug
(
 char*               filename_pattern,
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 double             *pts_coord,
 PDM_g_num_t        *pts_g_num
)
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int my_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &my_rank);
  char filename[999];
  //-->>
  double *_pts_coord = malloc (sizeof(double) * n_pts * 3);
  for (int i = 0; i < n_pts; i++) {
    _normalize (_dbbt,
                pts_coord + 3*i,
                _pts_coord + 3*i);
  }

  sprintf(filename, "%s_dbbt_pts_n_%3.3d.vtk", filename_pattern, my_rank);
  _export_point_cloud (filename,
                       1,
                       &n_pts,
                       &_pts_coord,
                       &pts_g_num,
                       NULL);

  sprintf(filename, "%s_dbbt_boxes_n_%3.3d.vtk", filename_pattern, my_rank);
  PDM_vtk_write_boxes (filename,
                       _dbbt->boxes->local_boxes->n_boxes,
                       _dbbt->boxes->local_boxes->extents,
                       _dbbt->boxes->local_boxes->g_num);


  double *_extents = malloc (sizeof(double) * _dbbt->boxes->local_boxes->n_boxes * 6);
  for (int i = 0; i < 2*_dbbt->boxes->local_boxes->n_boxes; i++) {
    for (int j = 0; j < 3; j++) {
      _extents[3*i+j] = _dbbt->s[j] + _dbbt->d[j] * _dbbt->boxes->local_boxes->extents[3*i+j];
    }
  }

  sprintf(filename, "%s_dbbt_pts_%3.3d.vtk", filename_pattern, my_rank);
  _export_point_cloud (filename,
                       1,
                       &n_pts,
                       &pts_coord,
                       &pts_g_num,
                       NULL);

  sprintf(filename, "%s_dbbt_boxes_%3.3d.vtk", filename_pattern, my_rank);
  PDM_vtk_write_boxes (filename,
                       _dbbt->boxes->local_boxes->n_boxes,
                       _extents,
                       _dbbt->boxes->local_boxes->g_num);

  free(_extents);
  free(_pts_coord);
}



/**
 *
 * \brief Get an indexed list of all points inside the boxes of a distributed box tree
 *
 *   \param [in]  dbbt               Pointer to distributed box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts_g_num          Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord          Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  n_boxes            Number of boxes
 *   \param [in]  box_g_num          Global ids of boxes (size = \ref n_boxes)
 *   \param [out] pts_in_box_idx     Index of points in boxes (size = \ref n_boxes + 1, allocated inside function)
 *   \param [out] pts_in_box_g_num   Global ids of points in boxes (size = \ref pts_in_box_idx[\ref n_boxes], allocated inside function)
 *   \param [out] pts_in_box_coord   Coordinates of points in boxes (size = 3 *\ref pts_in_box_idx[\ref n_boxes], allocated inside function)
 *   \param [in]  ellipsoids         Consider boxes as axis-aligned ellipsoids (1 or 0)
 *
 */

void
PDM_dbbtree_points_inside_boxes
(
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 PDM_g_num_t         pts_g_num[],
 double              pts_coord[],
 const int           n_boxes,
 const PDM_g_num_t   box_g_num[],
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord,
 const int           ellipsoids
 )
{
  const float f_threshold = 1.1;  // factor of the mean nb of requests
  const float f_max_copy  = 0.1;  // factor of the total nb of processes

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  /*
   *  Normalize target points coords
   */
  double *_pts_coord = malloc (sizeof(double) * n_pts * 3);
  for (int i = 0; i < n_pts; i++) {
    _normalize (_dbbt,
                pts_coord  + 3*i,
                _pts_coord + 3*i);
  }

  /*
   *  For each point, find all ranks that might have boxes containing that point
   */
  int *send_count = NULL;
  int *send_shift = NULL;
  int *recv_count = NULL;
  int *recv_shift = NULL;
  PDM_g_num_t *send_g_num = NULL;
  PDM_g_num_t *recv_g_num = NULL;
  double      *send_coord = NULL;
  double      *recv_coord = NULL;

  int n_copied_ranks = 0;
  int *copied_ranks = NULL;
  int n_pts_local  = 0;
  int n_pts_recv   = 0;
  int n_pts_copied = 0;
  int n_pts1;

  int *copied_shift = NULL;

  PDM_g_num_t *pts_g_num1 = NULL;
  double      *pts_coord1 = NULL;

  if (_dbbt->btShared != NULL) {
    int *pts_rank_idx = NULL;
    int *pts_rank = NULL;
    PDM_box_tree_boxes_containing_points (_dbbt->btShared,
                                       -1,
                                       n_pts,
                                       pts_coord,
                                       &pts_rank_idx,
                                       &pts_rank);

    /* Count points to send to each rank */
    send_count = PDM_array_zeros_int (n_rank);

    for (int i = 0; i < pts_rank_idx[n_pts]; i++) {
      int rank = _dbbt->usedRank[pts_rank[i]];
      pts_rank[i] = rank;
      send_count[rank]++;
    }

    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    n_pts_recv = 0;
    for (int i = 0; i < n_rank; i++) {
      n_pts_recv += recv_count[i];
    }
    int n_pts_recv_no_copies = n_pts_recv;

    /* Prepare copies */
    int n_max_copy = (int) (f_max_copy * n_rank);
    int *i_copied_rank = NULL;
    int *copied_count = NULL;
    int mean_n_pts_recv = 0;
    int *n_pts_recv_copied_ranks = NULL;

    if (n_max_copy > 0) {
      int *all_n_pts_recv = malloc (sizeof(int) * n_rank);
      PDM_MPI_Allgather (&n_pts_recv,    1, PDM_MPI_INT,
                         all_n_pts_recv, 1, PDM_MPI_INT,
                         _dbbt->comm);

      // Mean number of pts_recvs
      long l_mean_n_pts_recv = 0;
      for (int i = 0; i < n_rank; i++) {
        l_mean_n_pts_recv += all_n_pts_recv[i];
      }
      mean_n_pts_recv = (int) (l_mean_n_pts_recv / n_rank);

      float n_threshold = f_threshold * mean_n_pts_recv;

      // Sort ranks
      int *order = malloc (sizeof(int) * n_rank);
      for (int i = 0; i < n_rank; i++) {
        order[i] = i;
      }

      PDM_sort_int (all_n_pts_recv,
                    order,
                    n_rank);

      // Identify ranks to copy
      copied_ranks = malloc (sizeof(int) * n_max_copy);
      n_pts_recv_copied_ranks = malloc (sizeof(int) * n_max_copy);
      for (int i = 0; i < n_max_copy; i++) {
        int j = n_rank - i - 1;

        if (all_n_pts_recv[j] > n_threshold) {
          copied_ranks[n_copied_ranks] = order[j];
          n_pts_recv_copied_ranks[n_copied_ranks] = all_n_pts_recv[j];
          n_copied_ranks++;
        }
        else {
          break;
        }
      }
      free (all_n_pts_recv);
      free (order);

      if (n_copied_ranks > 0) {
        copied_ranks = realloc (copied_ranks, sizeof(int) * n_copied_ranks);
        n_pts_recv_copied_ranks = realloc (n_pts_recv_copied_ranks,
                                           sizeof(int) * n_copied_ranks);
      }
    }

    if (i_rank == 0) {
      if (n_copied_ranks > 0) {
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
      else {
        printf("0 copied ranks\n");
      }
    }

    /*
     * Copy the data of selected ranks
     */
    i_copied_rank = malloc (sizeof(int) * n_rank);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                &n_copied_ranks,
                                copied_ranks,
                                i_copied_rank);

    /* Re-compute send/recv counts */
    copied_count = PDM_array_zeros_int (n_copied_ranks);
    n_pts_local = 0;
    for (int i = 0; i < n_copied_ranks; i++) {
      int rank = copied_ranks[i];
      if (rank != i_rank) {
        int si = send_count[rank];

        si = PDM_MIN (si, PDM_MAX (0, (n_pts_recv_copied_ranks[i] - n_pts_recv)/2));
        if (i_copied_rank[i_rank] < 0) {
          si = PDM_MIN (si, PDM_MAX (0, mean_n_pts_recv - n_pts_recv));
        }

        copied_count[i] = si;
        n_pts_recv += si;
      }
    }
    if (copied_ranks != NULL) {
      free (copied_ranks);
    }

    if (n_pts_recv_copied_ranks != NULL) {
      free (n_pts_recv_copied_ranks);
    }

    for (int i = 0; i < n_rank; i++) {
      if (i == i_rank) {
        n_pts_local += send_count[i];
        send_count[i] = 0;
      }
      else if (i_copied_rank[i] >= 0) {
        send_count[i] -= copied_count[i_copied_rank[i]];
      }
    }

    copied_shift = PDM_array_new_idx_from_sizes_int (copied_count, n_copied_ranks);
    int *copied_count_tmp = PDM_array_zeros_int (n_copied_ranks);
    n_pts_copied = copied_shift[n_copied_ranks];

    /* Exchange new send/recv counts */
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    send_shift = PDM_array_new_idx_from_sizes_int (send_count, n_rank);
    recv_shift = PDM_array_new_idx_from_sizes_int (recv_count, n_rank);
    PDM_array_reset_int (send_count, n_rank, 0);

    n_pts_recv = recv_shift[n_rank];
    n_pts1 = n_pts_local + n_pts_recv + n_pts_copied;
    printf("[%d] n_pts1 = %d (without copies : %d)\n", i_rank, n_pts1, n_pts_recv_no_copies);

    pts_g_num1 = malloc (sizeof(PDM_g_num_t) * n_pts1);
    pts_coord1 = malloc (sizeof(double)      * n_pts1 * 3);


    /* Fill send buffers */
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank] * 3);

    int idx_copied = n_pts_local + n_pts_recv;
    PDM_g_num_t *copied_g_num = pts_g_num1 + idx_copied;
    double      *copied_coord = pts_coord1 + idx_copied * 3;
    n_pts_local = 0;
    for (int ipt = 0; ipt < n_pts; ipt++) {
      for (int i = pts_rank_idx[ipt]; i < pts_rank_idx[ipt+1]; i++) {
        int rank = pts_rank[i];

        if (rank == i_rank) {
          pts_g_num1[n_pts_local] = pts_g_num[ipt];
          for (int j = 0; j < 3; j++) {
            pts_coord1[3*n_pts_local + j] = _pts_coord[3*ipt + j];
          }
          n_pts_local++;
        }

        else if (i_copied_rank[rank] >= 0) {
          int _rank = i_copied_rank[rank];

          if (copied_count_tmp[_rank] < copied_count[_rank]) {
            int k = copied_shift[_rank] + copied_count_tmp[_rank];
            copied_g_num[k] = pts_g_num[ipt];
            for (int j = 0; j < 3; j++) {
              copied_coord[3*k + j] = _pts_coord[3*ipt + j];
            }
            copied_count_tmp[_rank]++;
          }
          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = pts_g_num[ipt];
            for (int j = 0; j < 3; j++) {
              send_coord[3*k + j] = _pts_coord[3*ipt + j];
            }
            send_count[rank]++;
          }
        }

        else {
          int k = send_shift[rank] + send_count[rank];
          send_g_num[k] = pts_g_num[ipt];
          for (int j = 0; j < 3; j++) {
            send_coord[3*k + j] = _pts_coord[3*ipt + j];
          }
          send_count[rank]++;
        }
      }
    }
    free (_pts_coord);
    if (copied_count != NULL) {
      free (copied_count);
    }
    if (copied_count_tmp != NULL) {
      free (copied_count_tmp);
    }
    if (i_copied_rank != NULL) {
      free (i_copied_rank);
    }
    free (pts_rank);
    free (pts_rank_idx);


    /* Exchange points */
    recv_g_num = pts_g_num1 + n_pts_local;
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    free (send_g_num);

    for (int i = 0; i < n_rank; i++) {
      send_count[i] *= 3;
      recv_count[i] *= 3;
      send_shift[i+1] *= 3;
      recv_shift[i+1] *= 3;
    }

    recv_coord = pts_coord1 + 3*n_pts_local;
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free (send_coord);
    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);
  }

  else {
    n_pts_local  = n_pts;
    n_pts_recv   = 0;
    n_pts_copied = 0;

    n_pts1 = n_pts;

    pts_g_num1 = (PDM_g_num_t *) pts_g_num;
    pts_coord1 = _pts_coord;
  }



  /*
   *  Get boxes containing points
   */
  int n_part = 1 + n_copied_ranks;
  int *part_n_pts = malloc (sizeof(int) * n_part);
  part_n_pts[0] = n_pts_local + n_pts_recv;

  int **pts_box_idx   = malloc (sizeof(int *) * n_part);
  int **pts_box_l_num = malloc (sizeof(int *) * n_part);

  int size_pts_box = 0;

  void (*_points_inside_volumes) (PDM_box_tree_t *bt,
                                  const int,
                                  const int,
                                  const double *,
                                  int **,
                                  int **);
  if (ellipsoids) {
    _points_inside_volumes = &PDM_box_tree_ellipsoids_containing_points;
  } else {
    _points_inside_volumes = &PDM_box_tree_boxes_containing_points;
  }

  /*
   *  Search in local tree
   */
  _points_inside_volumes (_dbbt->btLoc,
                          -1,
                          part_n_pts[0],
                          pts_coord1,
                          &(pts_box_idx[0]),
                          &(pts_box_l_num[0]));
  size_pts_box += pts_box_idx[0][part_n_pts[0]];


  /*
   *  Search in copied trees
   */
  double      *pts_coord_copied = NULL;
  PDM_g_num_t *pts_g_num_copied = NULL;
  if (n_copied_ranks > 0) {
    pts_coord_copied = pts_coord1 + part_n_pts[0] * 3;
    pts_g_num_copied = pts_g_num1 + part_n_pts[0];
    for (int i = 0; i < n_copied_ranks; i++) {
      part_n_pts[i+1] = copied_shift[i+1] - copied_shift[i];

      _points_inside_volumes (_dbbt->btLoc,
                              i,
                              part_n_pts[i+1],
                              pts_coord_copied + copied_shift[i] * 3,
                              &(pts_box_idx[i+1]),
                              &(pts_box_l_num[i+1]));

      size_pts_box += pts_box_idx[i+1][part_n_pts[i+1]];
    }
  }
  if (copied_shift != NULL) free (copied_shift);


  PDM_g_num_t *part_box_g_num = malloc (sizeof(PDM_g_num_t) * size_pts_box);
  PDM_g_num_t *part_pts_g_num = malloc (sizeof(PDM_g_num_t) * size_pts_box);
  double      *part_pts_coord = malloc (sizeof(double)      * size_pts_box * 3);
  int idx = 0;
  for (int j = 0; j < part_n_pts[0]; j++) {
    for (int k = pts_box_idx[0][j]; k < pts_box_idx[0][j+1]; k++) {
      part_box_g_num[idx] = _dbbt->boxes->local_boxes->g_num[pts_box_l_num[0][k]];
      part_pts_g_num[idx] = pts_g_num1[j];
      for (int l = 0; l < 3; l++) {
        part_pts_coord[3*idx + l] = pts_coord1[3*j + l];
      }
      idx++;
    }
  }
  free (pts_box_idx[0]);
  free (pts_box_l_num[0]);

  for (int i = 0; i < n_copied_ranks; i++) {
    for (int j = 0; j < part_n_pts[i+1]; j++) {
      for (int k = pts_box_idx[i+1][j]; k < pts_box_idx[i+1][j+1]; k++) {
        part_box_g_num[idx] = _dbbt->boxes->rank_boxes[i].g_num[pts_box_l_num[i+1][k]];
        part_pts_g_num[idx] = pts_g_num_copied[j];
        for (int l = 0; l < 3; l++) {
          part_pts_coord[3*idx + l] = pts_coord_copied[3*j + l];
        }
        idx++;
      }
    }
    pts_coord_copied += part_n_pts[i+1] * 3;
    pts_g_num_copied += part_n_pts[i+1];
    free (pts_box_idx[i+1]);
    free (pts_box_l_num[i+1]);
  }
  free (part_n_pts);
  free (pts_box_idx);
  free (pts_box_l_num);
  free (pts_coord1);
  if (pts_g_num1 != pts_g_num) free (pts_g_num1);

  /*
   *  Part-to-block
   */
  int    *part_stride = malloc (sizeof(int)    * size_pts_box);
  double *weight      = malloc (sizeof(double) * size_pts_box);
  for (int i = 0; i < size_pts_box; i++) {
    part_stride[i] = 1;
    weight[i]      = 1.;
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &part_box_g_num,
                                                       &weight,
                                                       &size_pts_box,
                                                       1,
                                                       _dbbt->comm);
  free (weight);

  int *block_box_pts_n = NULL;
  PDM_g_num_t *block_box_pts_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &part_stride,
                          (void **) &part_pts_g_num,
                          &block_box_pts_n,
                          (void **) &block_box_pts_g_num);
  free (part_pts_g_num);

  /* Fix partial block stride */
  PDM_g_num_t l_max_box_g_num = 0;
  for (int i = 0; i < n_boxes; i++) {
    l_max_box_g_num = PDM_MAX (l_max_box_g_num, box_g_num[i]);
  }

  PDM_g_num_t g_max_box_g_num;
  PDM_MPI_Allreduce (&l_max_box_g_num, &g_max_box_g_num, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, _dbbt->comm);

  PDM_g_num_t *block_distrib_idx =
    PDM_part_to_block_adapt_partial_block_to_block (ptb,
                                                    &block_box_pts_n,
                                                    g_max_box_g_num);

  int *block_stride = NULL;
  double *block_box_pts_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          3*sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &part_stride,
                          (void **) &part_pts_coord,
                          &block_stride,
                          (void **) &block_box_pts_coord);
  free (block_stride);
  free (part_stride);
  free (part_pts_coord);

  /*
   *  Block-to-part
   */
  PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                       (const PDM_g_num_t **) &box_g_num,
                                                       &n_boxes,
                                                       1,
                                                       _dbbt->comm);

  int *pts_in_box_n = malloc (sizeof(int) * n_boxes);
  int one = 1;
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          &one,
                          block_box_pts_n,
                          NULL,
                          (void **) &pts_in_box_n);

  *pts_in_box_idx = PDM_array_new_idx_from_sizes_int (pts_in_box_n, n_boxes);

  *pts_in_box_g_num = malloc (sizeof(PDM_g_num_t) * (*pts_in_box_idx)[n_boxes]);
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_box_pts_n,
                          block_box_pts_g_num,
                          &pts_in_box_n,
                          (void **) pts_in_box_g_num);
  free (block_box_pts_g_num);

  *pts_in_box_coord = malloc (sizeof(double) * (*pts_in_box_idx)[n_boxes] * 3);
  PDM_block_to_part_exch_in_place (btp,
                          3*sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_box_pts_n,
                          block_box_pts_coord,
                          &pts_in_box_n,
                          (void **) pts_in_box_coord);
  free (block_box_pts_coord);
  free (block_box_pts_n);
  free (pts_in_box_n);

  btp = PDM_block_to_part_free (btp);
  ptb = PDM_part_to_block_free (ptb);
  free (block_distrib_idx);
  free (part_box_g_num);

  //-->>
  for (int i = 0; i < (*pts_in_box_idx)[n_boxes]; i++) {
    for (int j = 0; j < 3; j++) {
      (*pts_in_box_coord)[3*i+j] = _dbbt->s[j] + _dbbt->d[j] * ((*pts_in_box_coord)[3*i+j]);
    }
  }
  //<<--
}



/**
 *
 * \brief Get an indexed list of all boxes containing points
 *
 *   \param [in]  dbbt        Pointer to distributed box tree structure
 *   \param [in]  n_pts       Number of points
 *   \param [in]  pts_g_num   Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord   Point coordinates (size = 3 * \ref n_pts)
 *   \param [out] box_idx     Index of boxes (size = \ref n_pts + 1, allocated inside function)
 *   \param [out] box_g_num   Global ids of boxes (size = \ref box_idx[\ref n_pts], allocated inside function)
 *   \param [in]  ellipsoids  Consider boxes as axis-aligned ellipsoids (1 or 0)
 *
 */

void
PDM_dbbtree_boxes_containing_points
(
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 PDM_g_num_t         pts_g_num[],
 double              pts_coord[],
 int               **box_idx,
 PDM_g_num_t       **box_g_num,
 const int           ellipsoids
 )
{
  const float f_threshold = 1.1;  // factor of the mean nb of requests
  const float f_max_copy  = 0.1;  // factor of the total nb of processes

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  /*
   *  Normalize target points coords
   */
  double *_pts_coord = malloc (sizeof(double) * n_pts * 3);
  for (int i = 0; i < n_pts; i++) {
    _normalize (_dbbt,
                pts_coord  + 3*i,
                _pts_coord + 3*i);
  }

  /*
   *  For each point, find all ranks that might have boxes containing that point
   */
  int *send_count = NULL;
  int *send_shift = NULL;
  int *recv_count = NULL;
  int *recv_shift = NULL;
  PDM_g_num_t *send_g_num = NULL;
  PDM_g_num_t *recv_g_num = NULL;
  double      *send_coord = NULL;
  double      *recv_coord = NULL;

  int n_copied_ranks = 0;
  int *copied_ranks = NULL;
  int n_pts_local  = 0;
  int n_pts_recv   = 0;
  int n_pts_copied = 0;
  int n_pts1;

  int *copied_shift = NULL;

  PDM_g_num_t *pts_g_num1 = NULL;
  double      *pts_coord1 = NULL;

  if (_dbbt->btShared != NULL) {
    int *pts_rank_idx = NULL;
    int *pts_rank = NULL;
    PDM_box_tree_boxes_containing_points (_dbbt->btShared,
                                       -1,
                                       n_pts,
                                       pts_coord,
                                       &pts_rank_idx,
                                       &pts_rank);

    /* Count points to send to each rank */
    send_count = PDM_array_zeros_int (n_rank);

    for (int i = 0; i < pts_rank_idx[n_pts]; i++) {
      int rank = _dbbt->usedRank[pts_rank[i]];
      pts_rank[i] = rank;
      send_count[rank]++;
    }

    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    n_pts_recv = 0;
    for (int i = 0; i < n_rank; i++) {
      n_pts_recv += recv_count[i];
    }
    int n_pts_recv_no_copies = n_pts_recv;

    /* Prepare copies */
    int n_max_copy = (int) (f_max_copy * n_rank);
    int *i_copied_rank = NULL;
    int *copied_count = NULL;
    int mean_n_pts_recv = 0;
    int *n_pts_recv_copied_ranks = NULL;

    if (n_max_copy > 0) {
      int *all_n_pts_recv = malloc (sizeof(int) * n_rank);
      PDM_MPI_Allgather (&n_pts_recv,    1, PDM_MPI_INT,
                         all_n_pts_recv, 1, PDM_MPI_INT,
                         _dbbt->comm);

      // Mean number of pts_recvs
      long l_mean_n_pts_recv = 0;
      for (int i = 0; i < n_rank; i++) {
        l_mean_n_pts_recv += all_n_pts_recv[i];
      }
      mean_n_pts_recv = (int) (l_mean_n_pts_recv / n_rank);

      float n_threshold = f_threshold * mean_n_pts_recv;

      // Sort ranks
      int *order = malloc (sizeof(int) * n_rank);
      for (int i = 0; i < n_rank; i++) {
        order[i] = i;
      }

      PDM_sort_int (all_n_pts_recv,
                    order,
                    n_rank);

      // Identify ranks to copy
      copied_ranks = malloc (sizeof(int) * n_max_copy);
      n_pts_recv_copied_ranks = malloc (sizeof(int) * n_max_copy);
      for (int i = 0; i < n_max_copy; i++) {
        int j = n_rank - i - 1;

        if (all_n_pts_recv[j] > n_threshold) {
          copied_ranks[n_copied_ranks] = order[j];
          n_pts_recv_copied_ranks[n_copied_ranks] = all_n_pts_recv[j];
          n_copied_ranks++;
        }
        else {
          break;
        }
      }
      free (all_n_pts_recv);
      free (order);

      if (n_copied_ranks > 0) {
        copied_ranks = realloc (copied_ranks, sizeof(int) * n_copied_ranks);
        n_pts_recv_copied_ranks = realloc (n_pts_recv_copied_ranks,
                                           sizeof(int) * n_copied_ranks);
      }
    }

    if (i_rank == 0) {
      if (n_copied_ranks > 0) {
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
      else {
        printf("0 copied ranks\n");
      }
    }

    /*
     * Copy the data of selected ranks
     */
    i_copied_rank = malloc (sizeof(int) * n_rank);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                &n_copied_ranks,
                                copied_ranks,
                                i_copied_rank);

    /* Re-compute send/recv counts */
    copied_count = PDM_array_zeros_int (n_copied_ranks);
    n_pts_local = 0;
    for (int i = 0; i < n_copied_ranks; i++) {
      int rank = copied_ranks[i];
      if (rank != i_rank) {
        int si = send_count[rank];

        si = PDM_MIN (si, PDM_MAX (0, (n_pts_recv_copied_ranks[i] - n_pts_recv)/2));
        if (i_copied_rank[i_rank] < 0) {
          si = PDM_MIN (si, PDM_MAX (0, mean_n_pts_recv - n_pts_recv));
        }

        copied_count[i] = si;
        n_pts_recv += si;
      }
    }
    if (copied_ranks != NULL) {
      free (copied_ranks);
    }

    if (n_pts_recv_copied_ranks != NULL) {
      free (n_pts_recv_copied_ranks);
    }

    for (int i = 0; i < n_rank; i++) {
      if (i == i_rank) {
        n_pts_local += send_count[i];
        send_count[i] = 0;
      }
      else if (i_copied_rank[i] >= 0) {
        send_count[i] -= copied_count[i_copied_rank[i]];
      }
    }

    copied_shift = PDM_array_new_idx_from_sizes_int (copied_count, n_copied_ranks);
    int *copied_count_tmp = PDM_array_zeros_int (n_copied_ranks);
    n_pts_copied = copied_shift[n_copied_ranks];

    /* Exchange new send/recv counts */
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    send_shift = PDM_array_new_idx_from_sizes_int (send_count, n_rank);
    recv_shift = PDM_array_new_idx_from_sizes_int (recv_count, n_rank);
    PDM_array_reset_int (send_count, n_rank, 0);

    n_pts_recv = recv_shift[n_rank];
    n_pts1 = n_pts_local + n_pts_recv + n_pts_copied;
    printf("[%d] n_pts1 = %d (without copies : %d)\n", i_rank, n_pts1, n_pts_recv_no_copies);

    pts_g_num1 = malloc (sizeof(PDM_g_num_t) * n_pts1);
    pts_coord1 = malloc (sizeof(double)      * n_pts1 * 3);


    /* Fill send buffers */
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank] * 3);

    int idx_copied = n_pts_local + n_pts_recv;
    PDM_g_num_t *copied_g_num = pts_g_num1 + idx_copied;
    double      *copied_coord = pts_coord1 + idx_copied * 3;
    n_pts_local = 0;
    for (int ipt = 0; ipt < n_pts; ipt++) {
      for (int i = pts_rank_idx[ipt]; i < pts_rank_idx[ipt+1]; i++) {
        int rank = pts_rank[i];

        if (rank == i_rank) {
          pts_g_num1[n_pts_local] = pts_g_num[ipt];
          for (int j = 0; j < 3; j++) {
            pts_coord1[3*n_pts_local + j] = _pts_coord[3*ipt + j];
          }
          n_pts_local++;
        }

        else if (i_copied_rank[rank] >= 0) {
          int _rank = i_copied_rank[rank];

          if (copied_count_tmp[_rank] < copied_count[_rank]) {
            int k = copied_shift[_rank] + copied_count_tmp[_rank];
            copied_g_num[k] = pts_g_num[ipt];
            for (int j = 0; j < 3; j++) {
              copied_coord[3*k + j] = _pts_coord[3*ipt + j];
            }
            copied_count_tmp[_rank]++;
          }
          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = pts_g_num[ipt];
            for (int j = 0; j < 3; j++) {
              send_coord[3*k + j] = _pts_coord[3*ipt + j];
            }
            send_count[rank]++;
          }
        }

        else {
          int k = send_shift[rank] + send_count[rank];
          send_g_num[k] = pts_g_num[ipt];
          for (int j = 0; j < 3; j++) {
            send_coord[3*k + j] = _pts_coord[3*ipt + j];
          }
          send_count[rank]++;
        }
      }
    }
    free (_pts_coord);
    if (copied_count != NULL) {
      free (copied_count);
    }
    if (copied_count_tmp != NULL) {
      free (copied_count_tmp);
    }
    if (i_copied_rank != NULL) {
      free (i_copied_rank);
    }
    free (pts_rank);
    free (pts_rank_idx);


    /* Exchange points */
    recv_g_num = pts_g_num1 + n_pts_local;
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    free (send_g_num);

    for (int i = 0; i < n_rank; i++) {
      send_count[i] *= 3;
      recv_count[i] *= 3;
      send_shift[i+1] *= 3;
      recv_shift[i+1] *= 3;
    }

    recv_coord = pts_coord1 + 3*n_pts_local;
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free (send_coord);
    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);
  }

  else {
    n_pts_local  = n_pts;
    n_pts_recv   = 0;
    n_pts_copied = 0;

    n_pts1 = n_pts;

    pts_g_num1 = (PDM_g_num_t *) pts_g_num;
    pts_coord1 = _pts_coord;
  }



  /*
   *  Get boxes containing points
   */
  int n_part = 1 + n_copied_ranks;
  int *part_n_pts = malloc (sizeof(int) * n_part);
  part_n_pts[0] = n_pts_local + n_pts_recv;

  int **pts_box_idx   = malloc (sizeof(int *) * n_part);
  int **pts_box_l_num = malloc (sizeof(int *) * n_part);

  void (*_points_inside_volumes) (PDM_box_tree_t *bt,
                                  const int,
                                  const int,
                                  const double *,
                                  int **,
                                  int **);
  if (ellipsoids) {
    _points_inside_volumes = &PDM_box_tree_ellipsoids_containing_points;
  } else {
    _points_inside_volumes = &PDM_box_tree_boxes_containing_points;
  }

  /*
   *  Search in local tree
   */
  _points_inside_volumes (_dbbt->btLoc,
                          -1,
                          part_n_pts[0],
                          pts_coord1,
                          &(pts_box_idx[0]),
                          &(pts_box_l_num[0]));

  /*
   *  Search in copied trees
   */
  double      *pts_coord_copied = NULL;
  //PDM_g_num_t *pts_g_num_copied = NULL;
  if (n_copied_ranks > 0) {
    pts_coord_copied = pts_coord1 + part_n_pts[0] * 3;
    //pts_g_num_copied = pts_g_num1 + part_n_pts[0];
    for (int i = 0; i < n_copied_ranks; i++) {
      part_n_pts[i+1] = copied_shift[i+1] - copied_shift[i];

      _points_inside_volumes (_dbbt->btLoc,
                              i,
                              part_n_pts[i+1],
                              pts_coord_copied + copied_shift[i] * 3,
                              &(pts_box_idx[i+1]),
                              &(pts_box_l_num[i+1]));
    }
  }
  if (copied_shift != NULL) free (copied_shift);
  free (pts_coord1);

  PDM_g_num_t **pts_box_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  pts_box_g_num[0] = malloc (sizeof(PDM_g_num_t) * pts_box_idx[0][part_n_pts[0]]);
  for (int j = 0; j < part_n_pts[0]; j++) {
    for (int k = pts_box_idx[0][j]; k < pts_box_idx[0][j+1]; k++) {
      pts_box_g_num[0][k] = _dbbt->boxes->local_boxes->g_num[pts_box_l_num[0][k]];
    }
  }

  for (int i = 0; i < n_copied_ranks; i++) {
    pts_box_g_num[i+1] = malloc (sizeof(PDM_g_num_t) * pts_box_idx[i+1][part_n_pts[i+1]]);
    for (int j = 0; j < part_n_pts[i+1]; j++) {
      for (int k = pts_box_idx[i+1][j]; k < pts_box_idx[i+1][j+1]; k++) {
        pts_box_g_num[i+1][k] = _dbbt->boxes->rank_boxes[i].g_num[pts_box_l_num[i+1][k]];
      }
    }
  }


  PDM_g_num_t **part_pts_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  int    **part_stride = malloc (sizeof(int *)    * n_part);
  double **part_weight = malloc (sizeof(double *) * n_part);
  int idx = 0;
  for (int ipart = 0; ipart < n_part; ipart++) {
    part_pts_g_num[ipart] = pts_g_num1 + idx;
    idx += part_n_pts[ipart];

    part_stride[ipart] = malloc (sizeof(int)    * part_n_pts[ipart]);
    part_weight[ipart] = malloc (sizeof(double) * part_n_pts[ipart]);
    for (int i = 0; i < part_n_pts[ipart]; i++) {
      part_stride[ipart][i] = pts_box_idx[ipart][i+1] - pts_box_idx[ipart][i];
      part_weight[ipart][i] = (double) part_stride[ipart][i];
    }
    free (pts_box_idx[ipart]);
    free (pts_box_l_num[ipart]);
  }
  free (pts_box_idx);
  free (pts_box_l_num);


  if (0) {
    printf("[%d] n_pts1 = %d\n", i_rank, n_pts1);
    printf("[%d] before merge:\n", i_rank);
    for (int ipart = 0; ipart < n_part; ipart++) {
      for (int i = 0; i< part_n_pts[ipart]; i++) {
        printf("[%d] pt "PDM_FMT_G_NUM" inside boxes", i_rank, part_pts_g_num[ipart][i]);
        for (int j = pts_box_idx[ipart][i]; j < pts_box_idx[ipart][i+1]; j++) {
          printf(" "PDM_FMT_G_NUM, pts_box_g_num[ipart][j]);
        }
        printf("\n");
      }
    }
  }


  /*
   *  Merge results
   */
  /* 1) Part-to-Block */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       part_pts_g_num,
                                                       part_weight,
                                                       part_n_pts,
                                                       n_part,
                                                       _dbbt->comm);

  int *block_box_n = NULL;
  PDM_g_num_t *block_box_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          part_stride,
                          (void **) pts_box_g_num,
                          &block_box_n,
                          (void **) &block_box_g_num);

  for (int ipart = 0; ipart < n_part; ipart++) {
    free (part_stride[ipart]);
    free (part_weight[ipart]);
    free (pts_box_g_num[ipart]);
  }
  free (part_stride);
  free (part_weight);
  free (pts_box_g_num);


  /* Remove doubles */
  int idx1 = 0, idx2 = 0;
  int n_pts_block = PDM_part_to_block_n_elt_block_get (ptb);
  int max_n = 0;
  for (int i = 0; i < n_pts_block; i++) {
    max_n = PDM_MAX (max_n, block_box_n[i]);
  }

  int *order = malloc (sizeof(int) * max_n);
  idx1 = 0;
  idx2 = 0;
  for (int i = 0; i < n_pts_block; i++) {
    if (block_box_n[i] == 0) continue;

    PDM_g_num_t *_g_num1 = block_box_g_num + idx1;
    PDM_g_num_t *_g_num2 = block_box_g_num + idx2;

    for (int j = 0; j < block_box_n[i]; j++) {
      order[j] = j;
    }
    PDM_sort_long (_g_num1,
                   order,
                   block_box_n[i]);

    _g_num2[0] = _g_num1[0];
    int tmp_n = 1;
    for (int j = 1; j < block_box_n[i]; j++) {
      if (_g_num1[j] != _g_num2[tmp_n-1]) {
        _g_num2[tmp_n++] = _g_num1[j];
      }
    }

    idx1 += block_box_n[i];
    idx2 += tmp_n;
    block_box_n[i] = tmp_n;
  }
  free (order);

  /* Fix partial block stride */
  PDM_g_num_t l_max_g_num = 0;
  for (int i = 0; i < n_pts; i++) {
    l_max_g_num = PDM_MAX (l_max_g_num, pts_g_num[i]);
  }

  PDM_g_num_t g_max_g_num;
  PDM_MPI_Allreduce (&l_max_g_num, &g_max_g_num, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, _dbbt->comm);

  PDM_g_num_t *block_distrib_idx =
    PDM_part_to_block_adapt_partial_block_to_block (ptb,
                                                    &block_box_n,
                                                    g_max_g_num);

  /* 2) Block-to-Part */
  PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                       (const PDM_g_num_t **) &pts_g_num,
                                                       &n_pts,
                                                       1,
                                                       _dbbt->comm);

  int *box_n = malloc (sizeof(int) * n_pts);
  int one = 1;
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          &one,
                          block_box_n,
                          NULL,
                          (void **) &box_n);

  *box_idx = PDM_array_new_idx_from_sizes_int (box_n, n_pts);

  *box_g_num = malloc (sizeof(PDM_g_num_t) * (*box_idx)[n_pts]);
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_box_n,
                          block_box_g_num,
                          &box_n,
                          (void **) box_g_num);
  free (block_box_g_num);
  free (block_box_n);
  free (box_n);

  btp = PDM_block_to_part_free (btp);
  ptb = PDM_part_to_block_free (ptb);
  free (block_distrib_idx);
  free (part_n_pts);

  free (part_pts_g_num);
  if (pts_g_num1 != pts_g_num) free (pts_g_num1);
}



/**
 *
 * \brief Get an indexed list of all boxes intersecting lines
 *
 *   \param [in]  dbbt        Pointer to distributed box tree structure
 *   \param [in]  n_line      Number of points
 *   \param [in]  line_g_num  Line global ids (size = \ref n_line)
 *   \param [in]  line_coord  Line coordinates (size = 6 * \ref n_line)
 *   \param [out] box_idx     Index of boxes (size = \ref n_line + 1, allocated inside function)
 *   \param [out] box_g_num   Global ids of boxes (size = \ref box_idx[\ref n_line], allocated inside function)
 *
 */

void
PDM_dbbtree_lines_intersect_boxes
(
 PDM_dbbtree_t  *dbbt,
 const int       n_line,
 PDM_g_num_t    *line_g_num,
 double         *line_coord,
 int           **box_idx,
 PDM_g_num_t   **box_g_num
 )
{
  const float f_threshold = 1.1;  // factor of the mean nb of requests
  const float f_max_copy  = 0.1;  // factor of the total nb of processes

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  /*
   *  Normalize coordinates
   */
  double *_line_coord = malloc (sizeof(double) * n_line * 6);
  for (int i = 0; i < 2*n_line; i++) {
    _normalize (_dbbt,
                line_coord  + 3*i,
                _line_coord + 3*i);
  }



  /*
   *  For each line, find all ranks that might have boxes intersecting that line
   */
  int *send_count = NULL;
  int *send_shift = NULL;
  int *recv_count = NULL;
  int *recv_shift = NULL;
  PDM_g_num_t *send_g_num = NULL;
  PDM_g_num_t *recv_g_num = NULL;
  double      *send_coord = NULL;
  double      *recv_coord = NULL;

  int n_copied_ranks = 0;
  int *copied_ranks = NULL;
  int n_line_local  = 0;
  int n_line_recv   = 0;
  int n_line_copied = 0;
  int n_line1;

  int *copied_shift = NULL;

  PDM_g_num_t *line_g_num1 = NULL;
  double      *line_coord1 = NULL;

  if (_dbbt->btShared != NULL) {
    int *line_rank_idx = NULL;
    int *line_rank     = NULL;
    PDM_box_tree_intersect_lines_boxes (_dbbt->btShared,
                                        -1,
                                        n_line,
                                        line_coord,
                                        &line_rank_idx,
                                        &line_rank);

    /* Count points to send to each rank */
    send_count = PDM_array_zeros_int (n_rank);

    for (int i = 0; i < line_rank_idx[n_line]; i++) {
      int rank = _dbbt->usedRank[line_rank[i]];
      line_rank[i] = rank;
      send_count[rank]++;
    }

    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    n_line_recv = 0;
    for (int i = 0; i < n_rank; i++) {
      n_line_recv += recv_count[i];
    }
    int n_line_recv_no_copies = n_line_recv;

    /* Prepare copies */
    int n_max_copy = (int) (f_max_copy * n_rank);
    int *i_copied_rank = NULL;
    int *copied_count = NULL;
    int mean_n_line_recv = 0;
    int *n_line_recv_copied_ranks = NULL;

    if (n_max_copy > 0) {
      int *all_n_line_recv = malloc (sizeof(int) * n_rank);
      PDM_MPI_Allgather (&n_line_recv,    1, PDM_MPI_INT,
                         all_n_line_recv, 1, PDM_MPI_INT,
                         _dbbt->comm);

      // Mean number of line_recvs
      long l_mean_n_line_recv = 0;
      for (int i = 0; i < n_rank; i++) {
        l_mean_n_line_recv += all_n_line_recv[i];
      }
      mean_n_line_recv = (int) (l_mean_n_line_recv / n_rank);

      float n_threshold = f_threshold * mean_n_line_recv;

      // Sort ranks
      int *order = malloc (sizeof(int) * n_rank);
      for (int i = 0; i < n_rank; i++) {
        order[i] = i;
      }

      PDM_sort_int (all_n_line_recv,
                    order,
                    n_rank);

      // Identify ranks to copy
      copied_ranks = malloc (sizeof(int) * n_max_copy);
      n_line_recv_copied_ranks = malloc (sizeof(int) * n_max_copy);
      for (int i = 0; i < n_max_copy; i++) {
        int j = n_rank - i - 1;

        if (all_n_line_recv[j] > n_threshold) {
          copied_ranks[n_copied_ranks] = order[j];
          n_line_recv_copied_ranks[n_copied_ranks] = all_n_line_recv[j];
          n_copied_ranks++;
        }
        else {
          break;
        }
      }
      free (all_n_line_recv);
      free (order);

      if (n_copied_ranks > 0) {
        copied_ranks = realloc (copied_ranks, sizeof(int) * n_copied_ranks);
        n_line_recv_copied_ranks = realloc (n_line_recv_copied_ranks,
                                           sizeof(int) * n_copied_ranks);
      }
    }

    if (i_rank == 0) {
      if (n_copied_ranks > 0) {
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
      else {
        printf("0 copied ranks\n");
      }
    }

    /*
     * Copy the data of selected ranks
     */
    i_copied_rank = malloc (sizeof(int) * n_rank);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                &n_copied_ranks,
                                copied_ranks,
                                i_copied_rank);

    /* Re-compute send/recv counts */
    copied_count = PDM_array_zeros_int (n_copied_ranks);
    n_line_local = 0;
    for (int i = 0; i < n_copied_ranks; i++) {
      int rank = copied_ranks[i];
      if (rank != i_rank) {
        int si = send_count[rank];

        si = PDM_MIN (si, PDM_MAX (0, (n_line_recv_copied_ranks[i] - n_line_recv)/2));
        if (i_copied_rank[i_rank] < 0) {
          si = PDM_MIN (si, PDM_MAX (0, mean_n_line_recv - n_line_recv));
        }

        copied_count[i] = si;
        n_line_recv += si;
      }
    }
    if (copied_ranks != NULL) {
      free (copied_ranks);
    }

    if (n_line_recv_copied_ranks != NULL) {
      free (n_line_recv_copied_ranks);
    }

    for (int i = 0; i < n_rank; i++) {
      if (i == i_rank) {
        n_line_local += send_count[i];
        send_count[i] = 0;
      }
      else if (i_copied_rank[i] >= 0) {
        send_count[i] -= copied_count[i_copied_rank[i]];
      }
    }

    copied_shift = PDM_array_new_idx_from_sizes_int (copied_count, n_copied_ranks);
    int *copied_count_tmp = PDM_array_zeros_int (n_copied_ranks);
    n_line_copied = copied_shift[n_copied_ranks];

    /* Exchange new send/recv counts */
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    send_shift = PDM_array_new_idx_from_sizes_int (send_count, n_rank);
    recv_shift = PDM_array_new_idx_from_sizes_int (recv_count, n_rank);
    PDM_array_reset_int (send_count, n_rank, 0);

    n_line_recv = recv_shift[n_rank];
    n_line1 = n_line_local + n_line_recv + n_line_copied;
    printf("[%d] n_line1 = %d (without copies : %d)\n", i_rank, n_line1, n_line_recv_no_copies);

    line_g_num1 = malloc (sizeof(PDM_g_num_t) * n_line1);
    line_coord1 = malloc (sizeof(double)      * n_line1 * 6);


    /* Fill send buffers */
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank] * 6);

    int idx_copied = n_line_local + n_line_recv;
    PDM_g_num_t *copied_g_num = line_g_num1 + idx_copied;
    double      *copied_coord = line_coord1 + idx_copied * 6;
    n_line_local = 0;
    for (int ipt = 0; ipt < n_line; ipt++) {
      for (int i = line_rank_idx[ipt]; i < line_rank_idx[ipt+1]; i++) {
        int rank = line_rank[i];

        if (rank == i_rank) {
          line_g_num1[n_line_local] = line_g_num[ipt];
          for (int j = 0; j < 6; j++) {
            line_coord1[6*n_line_local + j] = _line_coord[6*ipt + j];
          }
          n_line_local++;
        }

        else if (i_copied_rank[rank] >= 0) {
          int _rank = i_copied_rank[rank];

          if (copied_count_tmp[_rank] < copied_count[_rank]) {
            int k = copied_shift[_rank] + copied_count_tmp[_rank];
            copied_g_num[k] = line_g_num[ipt];
            for (int j = 0; j < 6; j++) {
              copied_coord[6*k + j] = _line_coord[6*ipt + j];
            }
            copied_count_tmp[_rank]++;
          }
          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = line_g_num[ipt];
            for (int j = 0; j < 6; j++) {
              send_coord[6*k + j] = _line_coord[6*ipt + j];
            }
            send_count[rank]++;
          }
        }

        else {
          int k = send_shift[rank] + send_count[rank];
          send_g_num[k] = line_g_num[ipt];
          for (int j = 0; j < 6; j++) {
            send_coord[6*k + j] = _line_coord[6*ipt + j];
          }
          send_count[rank]++;
        }
      }
    }
    free (_line_coord);
    if (copied_count != NULL) {
      free (copied_count);
    }
    if (copied_count_tmp != NULL) {
      free (copied_count_tmp);
    }
    if (i_copied_rank != NULL) {
      free (i_copied_rank);
    }
    free (line_rank);
    free (line_rank_idx);


    /* Exchange points */
    recv_g_num = line_g_num1 + n_line_local;
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    free (send_g_num);

    for (int i = 0; i < n_rank; i++) {
      send_count[i]   *= 6;
      recv_count[i]   *= 6;
      send_shift[i+1] *= 6;
      recv_shift[i+1] *= 6;
    }

    recv_coord = line_coord1 + 6*n_line_local;
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free (send_coord);
    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);
  }

  else {
    n_line_local  = n_line;
    n_line_recv   = 0;
    n_line_copied = 0;

    n_line1 = n_line;

    line_g_num1 = (PDM_g_num_t *) line_g_num;
    line_coord1 = _line_coord;
  }



  /*
   *  Get boxes intersecting lines
   */
  int n_part = 1 + n_copied_ranks;
  int *part_n_line = malloc (sizeof(int) * n_part);
  part_n_line[0] = n_line_local + n_line_recv;

  int **line_box_idx   = malloc (sizeof(int *) * n_part);
  int **line_box_l_num = malloc (sizeof(int *) * n_part);

  /*
   *  Search in local tree
   */
  PDM_box_tree_intersect_lines_boxes (_dbbt->btLoc,
                                      -1,
                                      part_n_line[0],
                                      line_coord1,
                                      &(line_box_idx[0]),
                                      &(line_box_l_num[0]));

  /*
   *  Search in copied trees
   */
  double *line_coord_copied = NULL;
  if (n_copied_ranks > 0) {
    line_coord_copied = line_coord1 + part_n_line[0] * 6;
    for (int i = 0; i < n_copied_ranks; i++) {
      part_n_line[i+1] = copied_shift[i+1] - copied_shift[i];

      PDM_box_tree_intersect_lines_boxes (_dbbt->btLoc,
                                          i,
                                          part_n_line[i+1],
                                          line_coord_copied + copied_shift[i] * 6,
                                          &(line_box_idx[i+1]),
                                          &(line_box_l_num[i+1]));
    }
  }
  if (copied_shift != NULL) free (copied_shift);
  free (line_coord1);

  PDM_g_num_t **line_box_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  line_box_g_num[0] = malloc (sizeof(PDM_g_num_t) * line_box_idx[0][part_n_line[0]]);
  for (int j = 0; j < part_n_line[0]; j++) {
    for (int k = line_box_idx[0][j]; k < line_box_idx[0][j+1]; k++) {
      line_box_g_num[0][k] = _dbbt->boxes->local_boxes->g_num[line_box_l_num[0][k]];
    }
  }

  for (int i = 0; i < n_copied_ranks; i++) {
    line_box_g_num[i+1] = malloc (sizeof(PDM_g_num_t) * line_box_idx[i+1][part_n_line[i+1]]);
    for (int j = 0; j < part_n_line[i+1]; j++) {
      for (int k = line_box_idx[i+1][j]; k < line_box_idx[i+1][j+1]; k++) {
        line_box_g_num[i+1][k] = _dbbt->boxes->rank_boxes[i].g_num[line_box_l_num[i+1][k]];
      }
    }
  }


  PDM_g_num_t **part_line_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  int    **part_stride = malloc (sizeof(int *)    * n_part);
  double **part_weight = malloc (sizeof(double *) * n_part);
  int idx = 0;
  for (int ipart = 0; ipart < n_part; ipart++) {
    part_line_g_num[ipart] = line_g_num1 + idx;
    idx += part_n_line[ipart];

    part_stride[ipart] = malloc (sizeof(int)    * part_n_line[ipart]);
    part_weight[ipart] = malloc (sizeof(double) * part_n_line[ipart]);
    for (int i = 0; i < part_n_line[ipart]; i++) {
      part_stride[ipart][i] = line_box_idx[ipart][i+1] - line_box_idx[ipart][i];
      part_weight[ipart][i] = (double) part_stride[ipart][i];
    }
    free (line_box_idx[ipart]);
    free (line_box_l_num[ipart]);
  }
  free (line_box_idx);
  free (line_box_l_num);


  /*
   *  Merge results
   */
  /* 1) Part-to-Block */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       part_line_g_num,
                                                       part_weight,
                                                       part_n_line,
                                                       n_part,
                                                       _dbbt->comm);

  int *block_box_n = NULL;
  PDM_g_num_t *block_box_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          part_stride,
                          (void **) line_box_g_num,
                          &block_box_n,
                          (void **) &block_box_g_num);

  for (int ipart = 0; ipart < n_part; ipart++) {
    free (part_stride[ipart]);
    free (part_weight[ipart]);
    free (line_box_g_num[ipart]);
  }
  free (part_stride);
  free (part_weight);
  free (line_box_g_num);


  /* Remove doubles */
  int idx1 = 0, idx2 = 0;
  int n_line_block = PDM_part_to_block_n_elt_block_get (ptb);
  int max_n = 0;
  for (int i = 0; i < n_line_block; i++) {
    max_n = PDM_MAX (max_n, block_box_n[i]);
  }

  int *order = malloc (sizeof(int) * max_n);
  idx1 = 0;
  idx2 = 0;
  for (int i = 0; i < n_line_block; i++) {
    if (block_box_n[i] == 0) continue;

    PDM_g_num_t *_g_num1 = block_box_g_num + idx1;
    PDM_g_num_t *_g_num2 = block_box_g_num + idx2;

    for (int j = 0; j < block_box_n[i]; j++) {
      order[j] = j;
    }
    PDM_sort_long (_g_num1,
                   order,
                   block_box_n[i]);

    _g_num2[0] = _g_num1[0];
    int tmp_n = 1;
    for (int j = 1; j < block_box_n[i]; j++) {
      if (_g_num1[j] != _g_num2[tmp_n-1]) {
        _g_num2[tmp_n++] = _g_num1[j];
      }
    }

    idx1 += block_box_n[i];
    idx2 += tmp_n;
    block_box_n[i] = tmp_n;
  }
  free (order);

  /* Fix partial block stride */
  PDM_g_num_t l_max_g_num = 0;
  for (int i = 0; i < n_line; i++) {
    l_max_g_num = PDM_MAX (l_max_g_num, line_g_num[i]);
  }

  PDM_g_num_t g_max_g_num;
  PDM_MPI_Allreduce (&l_max_g_num, &g_max_g_num, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, _dbbt->comm);

  PDM_g_num_t *block_distrib_idx =
    PDM_part_to_block_adapt_partial_block_to_block (ptb,
                                                    &block_box_n,
                                                    g_max_g_num);

  /* 2) Block-to-Part */
  PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                       (const PDM_g_num_t **) &line_g_num,
                                                       &n_line,
                                                       1,
                                                       _dbbt->comm);

  int *box_n = malloc (sizeof(int) * n_line);
  int one = 1;
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          &one,
                          block_box_n,
                          NULL,
                          (void **) &box_n);

  *box_idx = PDM_array_new_idx_from_sizes_int (box_n, n_line);

  *box_g_num = malloc (sizeof(PDM_g_num_t) * (*box_idx)[n_line]);
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_box_n,
                          block_box_g_num,
                          &box_n,
                          (void **) box_g_num);
  free (block_box_g_num);
  free (block_box_n);
  free (box_n);

  btp = PDM_block_to_part_free (btp);
  ptb = PDM_part_to_block_free (ptb);
  free (block_distrib_idx);
  free (part_n_line);

  free (part_line_g_num);
  if (line_g_num1 != line_g_num) free (line_g_num1);
}






#undef _MIN
#undef _MAX

#ifdef __cplusplus
}
#endif /* __cplusplus */
