
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
#include "pdm_hash_tab.h"
#include "pdm_box_priv.h"
#include "pdm_box_tree.h"
#include "pdm_dbbtree.h"
#include "pdm_dbbtree_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"

#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

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
 PDM_MPI_Comm          comm,
 const int         dim,
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


PDM_box_set_t  *
PDM_dbbtree_boxes_set
(
 PDM_dbbtree_t     *dbbt,
 const int          n_part,
 const int         *nElts,
 const double     **extents,
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

    int *initLocation_proc = (int *) malloc (sizeof(int) * 3 * nUsedRank);
    for (int i = 0; i < 3 * nUsedRank; i++) {
      initLocation_proc[i] = 0;
    }

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
 * This function assigns boxes to intersect to the tree.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  n_part    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 *
 * \return associated \ref PDM_box_set_t structure distributed according
 * to the tree intersection
 *
 */

PDM_box_set_t  *
PDM_dbbtree_intersect_boxes_set
(
 PDM_dbbtree_t    *dbbt,
 const int         n_part,
 const int        *nElts,
 const double     **extents,
 const PDM_g_num_t **gNum,
 int              *box_index[],
 int              *box_l_num[]
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

    for (int i = 0; i < lComm + 1; i++) {
      distrib->index[i] = 0;
    }

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

  int *newIndex = (int *) malloc (sizeof(int) * (nBoxesA + 1));
  for (int i = 0; i < nBoxesA + 1; i++) {
    newIndex[i] = 0;
  }

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
 * Get the boxes closer than the upper bound distance
 *
 *   \param [in] bt               Pointer to box tree structure
 *   \param [in] n_pts            Number of points
 *   \param [in] pts              Point coordinates (size = 3 * n_pts)
 *   \param [in] pts_g_num        Point global numbers
 *   \param [in] upper_bound_dist2 Upper bound of the squer of the distance (size = n_pts)
 *   \param [out] box_index       Index of boxes (size = n_pts + 1)
 *   \param [out] box_g_num       Global num of boxes (size = i_boxes[n_pts])
 *
 */

void
PDM_dbbtree_closest_upper_bound_dist_boxes_get_OLD
(
 PDM_dbbtree_t    *dbbt,
 const int        n_pts,
 double           pts[],
 PDM_g_num_t      pts_g_num[],
 double           upper_bound_dist2[],
 int             *box_index[],
 PDM_g_num_t     *box_g_num[]
 )
{
  /*
   * Initialization
   */
  PDM_UNUSED(pts_g_num);

  int npts_in_rank = n_pts;

  double *pts_in_rank =  pts;
  double *upper_bound_dist_in_rank = upper_bound_dist2;

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  /*
   * Determination de liste des procs concernes pour chaque sommet
   */

  int *n_send_pts = NULL;
  int *i_send_pts = NULL;

  int *n_recv_pts = NULL;
  int *i_recv_pts = NULL;

  int *box_index_tmp = NULL;
  int *box_l_num_tmp = NULL;


  const int *usedRanks = _dbbt->usedRank;

  const int idebug = 0;

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
     * Envoi des points a chaque proc concerne
     */

    n_send_pts = malloc (sizeof(int) * lComm);
    i_send_pts = malloc (sizeof(int) * (lComm+1));

    n_recv_pts = malloc (sizeof(int) * lComm);
    i_recv_pts = malloc (sizeof(int) * (lComm+1));

    for (int i = 0; i < lComm; i++) {
      n_send_pts[i] = 0;
    }

    i_send_pts[0] = 0;
    i_recv_pts[0] = 0;

    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        n_send_pts[usedRanks[box_l_num_tmp[j]]]++;
      }
    }

    PDM_MPI_Alltoall (n_send_pts, 1, PDM_MPI_INT,
                      n_recv_pts, 1, PDM_MPI_INT,
                      _dbbt->comm);

    for (int i = 0; i < lComm; i++) {
      i_send_pts[i+1] = i_send_pts[i] + 4 * n_send_pts[i];
      i_recv_pts[i+1] = i_recv_pts[i] + 4 * n_recv_pts[i];
      n_recv_pts[i] *= 4;
    }

    double *send_pts = malloc (sizeof(double) * i_send_pts[lComm]);
    double *recv_pts = malloc (sizeof(double) * i_recv_pts[lComm]);

    for (int i = 0; i < lComm; i++) {
      n_send_pts[i] = 0;
    }

    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        const int i_rank = usedRanks[box_l_num_tmp[j]];
        int idx            = i_send_pts[i_rank] + n_send_pts[i_rank];
        send_pts[idx++]    = pts[3*i];
        send_pts[idx++]    = pts[3*i+1];
        send_pts[idx++]    = pts[3*i+2];
        send_pts[idx++]    = upper_bound_dist2[i];
        n_send_pts[i_rank] += 4;
      }
    }

    PDM_MPI_Alltoallv (send_pts, n_send_pts, i_send_pts, PDM_MPI_DOUBLE,
                       recv_pts, n_recv_pts, i_recv_pts, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free(send_pts);

    npts_in_rank = i_recv_pts[lComm] / 4;

    pts_in_rank =  malloc (sizeof(double) * 3 * npts_in_rank);
    upper_bound_dist_in_rank =  malloc (sizeof(double) * npts_in_rank);

    for (int i = 0; i < npts_in_rank; i++) {
      for (int j = 0; j < 3; j++) {
        pts_in_rank[3*i+j] = recv_pts[4*i+j];
      }
      upper_bound_dist_in_rank[i] = recv_pts[4*i+3];
    }

    free(recv_pts);

  }

  /*
   * Determination des candidats localement
   */

  int *box_index_in_rank;
  int *box_l_num_in_rank;

  PDM_box_tree_closest_upper_bound_dist_boxes_get (_dbbt->btLoc,
                                                   npts_in_rank,
                                                   pts_in_rank,
                                                   upper_bound_dist_in_rank,
                                                   &box_index_in_rank,
                                                   &box_l_num_in_rank);

  if (idebug) {
    printf (" **** PDM_dbbtree_closest_upper_bound_dist_boxes_get n_pts_rank : %d\n", npts_in_rank);
    for (int i = 0; i < npts_in_rank; i++) {
      printf ("%d : (%12.5e %12.5e %12.5e) %12.5e\n", i,
              pts_in_rank[3*i], pts_in_rank[3*i+1], pts_in_rank[3*i+2],
              upper_bound_dist_in_rank[i]);
      printf ("  boxes %d :" , box_index_in_rank[i+1] - box_index_in_rank[i]);
      for (int j = box_index_in_rank[i]; j < box_index_in_rank[i+1]; j++) {
        printf (" %d", box_l_num_in_rank[j]);
      }
      printf ("\n");
    }
  }

  if (_dbbt->btShared == NULL) {
    *box_index = box_index_in_rank;
    *box_g_num = malloc (sizeof(PDM_g_num_t) * box_index_in_rank[npts_in_rank]);
    const PDM_g_num_t *gnum_boxes = PDM_box_set_get_g_num (_dbbt->boxes);
    for (int i = 0; i < box_index_in_rank[npts_in_rank]; i++) {
      (*box_g_num)[i] = gnum_boxes[box_l_num_in_rank[i]];
    }
    free (box_l_num_in_rank);
  }

  else {
    /*
     * Retour des resultats (AlltoAll inverse)
     *     - Envoi du nombre de boites trouvees pour chaque point
     *     - Envoi du numero des boites en numabs
     */

    free (pts_in_rank);
    free (upper_bound_dist_in_rank);

    int *n_box_l_num_in_rank = malloc (sizeof(int) * npts_in_rank);

    for (int i = 0; i < npts_in_rank; i++) {
      n_box_l_num_in_rank[i] = box_index_in_rank[i+1] - box_index_in_rank[i];
    }

    for (int i = 0; i < lComm; i++) {
      i_send_pts[i+1] = i_send_pts[i+1]/4;
      i_recv_pts[i+1] = i_recv_pts[i+1]/4;
      n_send_pts[i]   = n_send_pts[i]/4;
      n_recv_pts[i]   = n_recv_pts[i]/4;
    }

    int *n_box_l_num_per_pts = malloc (sizeof(int) * i_send_pts[lComm]);

    PDM_MPI_Alltoallv (n_box_l_num_in_rank, n_recv_pts, i_recv_pts, PDM_MPI_INT,
                       n_box_l_num_per_pts, n_send_pts, i_send_pts, PDM_MPI_INT,
                       _dbbt->comm);

    int *n_send_pts2 = malloc (sizeof(int) * lComm);
    int *i_send_pts2 = malloc (sizeof(int) * (lComm+1));

    int *n_recv_pts2 = malloc (sizeof(int) * lComm);
    int *i_recv_pts2 = malloc (sizeof(int) * (lComm+1));

    for (int i = 0; i < lComm; i++) {
      n_send_pts2[i] = 0;
      n_recv_pts2[i] = 0;
    }

    for (int i = 0; i < lComm; i++) {
      for (int j = i_recv_pts[i]; j < i_recv_pts[i+1]; j++) {
        n_recv_pts2[i] += n_box_l_num_in_rank[j];
      }
      for (int j = i_send_pts[i]; j < i_send_pts[i+1]; j++) {
        n_send_pts2[i] += n_box_l_num_per_pts[j];
      }
    }

    free (n_box_l_num_in_rank);

    i_send_pts2[0] = 0;
    i_recv_pts2[0] = 0;
    for (int i = 0; i < lComm; i++) {
      i_send_pts2[i+1] = i_send_pts2[i] + n_send_pts2[i];
      i_recv_pts2[i+1] = i_recv_pts2[i] + n_recv_pts2[i];
    }

    PDM_g_num_t *box_g_num_in_rank =
      malloc(sizeof(PDM_g_num_t) * box_index_in_rank[npts_in_rank]);

    const PDM_g_num_t *gnum_boxes = PDM_box_set_get_g_num (_dbbt->boxes);

    for (int i = 0; i < box_index_in_rank[npts_in_rank]; i++) {
      box_g_num_in_rank[i] = gnum_boxes[box_l_num_in_rank[i]];
    }

    PDM_g_num_t *box_g_num_per_pts = malloc(sizeof(PDM_g_num_t) * i_send_pts2[lComm]);

    PDM_MPI_Alltoallv (box_g_num_in_rank, n_recv_pts2, i_recv_pts2, PDM__PDM_MPI_G_NUM,
                       box_g_num_per_pts, n_send_pts2, i_send_pts2, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);

    free (box_index_in_rank);
    free (box_l_num_in_rank);
    free (box_g_num_in_rank);

    free (n_recv_pts);
    free (i_recv_pts);
    free (n_recv_pts2);
    free (i_recv_pts2);

    /*
     * Tri du tableau de retour
     */

    *box_index = malloc(sizeof(int) * (n_pts + 1));
    int* box_n = malloc(sizeof(int) * n_pts);
    *box_g_num = malloc(sizeof(PDM_g_num_t) * i_send_pts2[lComm]);
    (*box_index)[0] = 0;

    for (int i = 0; i < n_pts; i++) {
      box_n[i] = 0;
    }

    for (int i = 0; i < lComm; i++) {
      n_send_pts[i] = 0;
    }

    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        const int i_rank = usedRanks[box_l_num_tmp[j]];
        const int idx   = i_send_pts[i_rank] + n_send_pts[i_rank];
        box_n[i]        += n_box_l_num_per_pts[idx];
        n_send_pts[i_rank] += 1;
      }
    }

    for (int i = 0; i < n_pts; i++) {
      (*box_index)[i+1] = (*box_index)[i] +  box_n[i];
    }

    for (int i = 0; i < lComm; i++) {
      n_send_pts[i] = 0;
      n_send_pts2[i] = 0;
    }

    int k1 = 0;

    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        const int i_rank = usedRanks[box_l_num_tmp[j]];
        int idx  = i_send_pts[i_rank] + n_send_pts[i_rank];
        int idx2 = i_send_pts2[i_rank] + n_send_pts2[i_rank];
        for (int k = 0; k < n_box_l_num_per_pts[idx]; k++) {
          (*box_g_num)[k1++] = box_g_num_per_pts[idx2++];
        }
        n_send_pts2[i_rank] += n_box_l_num_per_pts[idx];
        n_send_pts[i_rank] += 1;
      }
    }

    int keyMax = 3 * n_pts;
    PDM_hash_tab_t * ht = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                               &keyMax);

    box_index_tmp[0] = 0;
    int idx = 0;
    for (int i = 0; i < n_pts; i++) {
      box_index_tmp[i+1] = box_index_tmp[i];
      for (int j = (*box_index)[i]; j < (*box_index)[i+1]; j++) {
        PDM_g_num_t curr_box = (*box_g_num)[j];
        int key = curr_box % keyMax;

        int n_data = PDM_hash_tab_n_data_get (ht, &key);

        int found = 0;

        PDM_g_num_t **data = (PDM_g_num_t **) PDM_hash_tab_data_get (ht, &key);
        for (int k = 0; k < n_data; k++) {
          if (*(data[k]) == curr_box) {
            found = 1;
            break;
          }
        }

        if (!found) {
          PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
          (*box_g_num)[idx++] = curr_box;
          box_index_tmp[i+1] += 1;
        }
      }
      PDM_hash_tab_purge (ht, PDM_FALSE);
    }

    PDM_hash_tab_free (ht);

    free (n_box_l_num_per_pts);

    free (*box_index);
    *box_index = box_index_tmp;

    *box_g_num = realloc (*box_g_num, sizeof(PDM_g_num_t) * box_index_tmp[n_pts]);

    free (box_l_num_tmp);
    free (box_g_num_per_pts);
    free (box_n);

    free (n_send_pts);
    free (i_send_pts);
    free (n_send_pts2);
    free (i_send_pts2);
  }

}


void
PDM_dbbtree_closest_upper_bound_dist_boxes_get
(
 PDM_dbbtree_t    *dbbt,
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
  int     n_pts_local            = n_pts;
  double *pts_local              = pts;
  double *upper_bound_dist_local = upper_bound_dist2;

    assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

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
    n_send_pts = malloc (sizeof(int) * lComm);
    n_recv_pts = malloc (sizeof(int) * lComm);

    for (int i = 0; i < lComm; i++) {
      n_send_pts[i] = 0;
    }

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

    n_pts_rank = malloc (sizeof(int) * n_copied_ranks);
    for (int i = 0; i < n_copied_ranks; i++) {
      n_pts_rank[i] = 0;
    }

    n_pts_send = malloc (sizeof(int) * lComm);
    n_pts_recv = malloc (sizeof(int) * lComm);

    for (int i = 0; i < lComm; i++) {
      n_pts_send[i] = 0;
    }

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

    i_pts_rank = malloc (sizeof(int) * (n_copied_ranks+1));
    i_pts_rank[0] = 0;
    for (int i = 0; i < n_copied_ranks; i++) {
      i_pts_rank[i+1] = i_pts_rank[i] + n_pts_rank[i];
    }

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

    for (int i = 0; i < lComm; i++) {
      n_pts_send[i] = 0;
    }
    for (int i = 0; i < n_copied_ranks; i++) {
      n_pts_rank[i] = 0;
    }

    i1 = 0; i2 = 0; i3 = 0;
    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        i_rank = usedRanks[box_l_num_tmp[j]];
        if ( i_rank == myRank ) {
          // pts_local, upper_bound_dist_local (local points, local box tree)
          upper_bound_dist_local[i1] = upper_bound_dist2[i];
          pts_local[3*i1]            = pts[3*i];
          pts_local[3*i1+1]          = pts[3*i+1];
          pts_local[3*i1+2]          = pts[3*i+2];
          i1++;
        } else if ( rank_copy_num[i_rank] >= 0 ) {
          // pts_rank, upper_bound_dist_rank (local points, distant (copied) box trees)
          int j_rank = rank_copy_num[i_rank];
          i2 = i_pts_rank[j_rank] + n_pts_rank[j_rank];
          upper_bound_dist_rank[i2] = upper_bound_dist2[i];
          pts_rank[3*i2]            = pts[3*i];
          pts_rank[3*i2+1]          = pts[3*i+1];
          pts_rank[3*i2+2]          = pts[3*i+2];
          n_pts_rank[j_rank]++;
        } else {
          // data_send (local points, distant (not copied) box trees)
          i3 = i_pts_send[i_rank] + 4*n_pts_send[i_rank];
          data_send[i3++] = pts[3*i];
          data_send[i3++] = pts[3*i+1];
          data_send[i3++] = pts[3*i+2];
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



  // Determine candidate boxes in local box tree (local points)
  int *box_index_local;
  int *box_l_num_local;
  #if 0
  PDM_box_tree_closest_upper_bound_dist_boxes_get (_dbbt->btLoc,
						   n_pts_local,
						   pts_local,
						   upper_bound_dist_local,
						   &box_index_local,
						   &box_l_num_local);
  #else
  PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
						      -1, // search in local box tree
						      n_pts_local,
						      pts_local,
						      upper_bound_dist_local,
						      &box_index_local,
						      &box_l_num_local);
  #endif
  free(pts_local);
  free(upper_bound_dist_local);

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

    #if 0
    PDM_box_tree_closest_upper_bound_dist_boxes_get_from_copied_ranks (_dbbt->btLoc,
								       i_pts_rank,
								       pts_rank,
								       upper_bound_dist_rank,
								       &box_index_rank,
								       &box_l_num_rank);
    #else
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
							  &(box_l_num_rank[i_copied_rank]));
    }
    #endif

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

    #if 0
    PDM_box_tree_closest_upper_bound_dist_boxes_get (_dbbt->btLoc,
                                                     n_pts_recv_total,
                                                     pts_recv,
                                                     upper_bound_dist_recv,
                                                     &box_index_recv,
                                                     &box_l_num_recv);
    #else
    PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
							-1, // search in local box tree
							n_pts_recv_total,
							pts_recv,
							upper_bound_dist_recv,
							&box_index_recv,
							&box_l_num_recv);
    #endif
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

    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      n_pts_rank[i_copied_rank] = 0;
    }



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
 * Get the location of a point cloud
 *
 *   \param [in] bt                    Pointer to box tree structure
 *   \param [in] n_pts                 Number of points
 *   \param [in] pts_coord             Point coordinates (size = 3 * n_pts)
 *   \param [in] pts_g_num             Point global numbers (size = n_pts)
 *   \param [in] n_boxes               Number of boxes
 *   \param [in] box_g_num             Global num of boxes (size = n_boxes)
 *   \param [inout] pts_in_box_idx     Index of points in boxes (size = n_boxes+1, allocated inside function)
 *   \param [inout] pts_in_box_g_num   Global num of points in boxes (size = pts_in_box_idx[n_boxes], allocated inside function)
 *   \param [inout] pts_in_box_coord   Coordinates of points in boxes (size = 3*pts_in_box_idx[n_boxes], allocated inside function)
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
 double            **pts_in_box_coord
 )
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;



  int my_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &my_rank);

  int n_ranks;
  PDM_MPI_Comm_size (_dbbt->comm, &n_ranks);

  //const int *used_ranks = _dbbt->usedRank;

  /***************************************
   * Redistribute points
   ***************************************/
  int *send_count = NULL;
  int *send_shift = NULL;

  int *recv_count = NULL;
  int *recv_shift = NULL;

  int n_recv_pts = 0;
  PDM_g_num_t *send_g_num = NULL;
  PDM_g_num_t *recv_g_num = NULL;
  double      *send_coord = NULL;
  double      *recv_coord = NULL;

  if (_dbbt->btShared != NULL) {
    PDM_box_tree_points_inside_boxes (_dbbt->btShared,
                                      n_pts,
                                      pts_g_num,
                                      pts_coord,
                                      &send_shift,
                                      &send_g_num,
                                      &send_coord);

    /* Count points to send to each rank */
    send_count = malloc (sizeof(int) * n_ranks);
    for (int i = 0; i < n_ranks; i++) {
      send_count[i] = send_shift[i+1] - send_shift[i];
    }

    recv_count = malloc (sizeof(int) * n_ranks);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    recv_shift = malloc (sizeof(int) * (n_ranks + 1));
    recv_shift[0] = 0;
    for (int i = 0; i < n_ranks; i++) {
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }
    n_recv_pts = recv_shift[n_ranks];

    /* Exchange points */
    recv_g_num = malloc (sizeof(PDM_g_num_t) * recv_shift[n_ranks]);
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    free (send_g_num);


    for (int i = 0; i < n_ranks; i++) {
      send_count[i] *= 3;
      recv_count[i] *= 3;
      send_shift[i+1] *= 3;
      recv_shift[i+1] *= 3;
    }

    recv_coord = malloc (sizeof(double) * recv_shift[n_ranks]);
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free (send_coord);
    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);
  }

  /* dbbt == NULL (Single rank) */
  else {
    n_recv_pts = n_pts;
    recv_g_num = (PDM_g_num_t *) pts_g_num;
    recv_coord = (double *) pts_coord;
  }



  /***************************************
   * Get list of boxes that contain each redistributed point
   ***************************************/
  int         *pts_in_box_idx2   = NULL;
  PDM_g_num_t *pts_in_box_g_num2 = NULL;
  double      *pts_in_box_coord2 = NULL;
  PDM_box_tree_points_inside_boxes (_dbbt->btLoc,
                                    n_recv_pts,
                                    recv_g_num,
                                    recv_coord,
                                    &pts_in_box_idx2,
                                    &pts_in_box_g_num2,
                                    &pts_in_box_coord2);


  /* Conform to original partitioning */
  if (box_g_num != NULL) {

    /* Part2 to Block */
    int n_boxes2 = _dbbt->boxes->local_boxes->n_boxes;
    PDM_part_to_block_t *ptb2 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                          &(_dbbt->boxes->local_boxes->g_num),
                                                          NULL,
                                                          &n_boxes2,
                                                          1,
                                                          _dbbt->comm);

    int *pts_in_box_count = malloc (sizeof(int) * n_boxes2);
    for (int ibox = 0; ibox < n_boxes2; ibox++) {
      pts_in_box_count[ibox] = pts_in_box_idx2[ibox+1] - pts_in_box_idx2[ibox];
    }
    free (pts_in_box_idx2);

    int *block_pts_in_box_count = NULL;
    PDM_g_num_t *block_pts_in_box_g_num = NULL;
    PDM_part_to_block_exch (ptb2,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            1,
                            &pts_in_box_count,
                            (void **) &pts_in_box_g_num2,
                            &block_pts_in_box_count,
                            (void **) &block_pts_in_box_g_num);
    free (pts_in_box_g_num2);


    for (int ibox = 0; ibox < n_boxes2; ibox++) {
      pts_in_box_count[ibox] *= 3;
    }
    double *block_pts_in_box_coord = NULL;
    int *block_pts_in_box_3count = NULL;
    PDM_part_to_block_exch (ptb2,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &pts_in_box_count,
                            (void **) &pts_in_box_coord2,
                            &block_pts_in_box_3count,
                            (void **) &block_pts_in_box_coord);
    free (pts_in_box_coord2);

    PDM_part_to_block_free (ptb2);


    /* Block to Part1 */
    int _n_boxes = n_boxes;
    PDM_part_to_block_t *ptb1 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                          (PDM_g_num_t **) (&box_g_num),
                                                          NULL,
                                                          &_n_boxes,
                                                          1,
                                                          _dbbt->comm);
    PDM_g_num_t *block_distrib_idx1 = PDM_part_to_block_distrib_index_get (ptb1);

    PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx1,
                                                         (const PDM_g_num_t **) &box_g_num,
                                                         &n_boxes,
                                                         1,
                                                         _dbbt->comm);

    if (n_boxes > n_boxes2) {
      pts_in_box_count = realloc (pts_in_box_count, sizeof(int) * n_boxes);
    }
    int one = 1;
    PDM_block_to_part_exch (btp,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            &one,
                            block_pts_in_box_count,
                            NULL,
                            (void **) &pts_in_box_count);

    *pts_in_box_idx = malloc (sizeof(int) * (n_boxes + 1));
    (*pts_in_box_idx)[0] = 0;
    for (int ibox = 0; ibox < n_boxes; ibox++) {
      (*pts_in_box_idx)[ibox+1] = (*pts_in_box_idx)[ibox] + pts_in_box_count[ibox];
    }

    *pts_in_box_g_num = malloc (sizeof(PDM_g_num_t) * (*pts_in_box_idx)[n_boxes]);
    PDM_block_to_part_exch (btp,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            block_pts_in_box_count,
                            block_pts_in_box_g_num,
                            &pts_in_box_count,
                            (void **) pts_in_box_g_num);
    free (block_pts_in_box_g_num);
    free (block_pts_in_box_count);

    *pts_in_box_coord = malloc (sizeof(double) * (*pts_in_box_idx)[n_boxes] * 3);
    PDM_block_to_part_exch (btp,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            block_pts_in_box_3count,
                            block_pts_in_box_coord,
                            &pts_in_box_count,
                            (void **) pts_in_box_coord);
    free (block_pts_in_box_coord);
    free (block_pts_in_box_3count);
    free (pts_in_box_count);


    PDM_part_to_block_free (ptb1);
    PDM_block_to_part_free (btp);
  }

  else {
    *pts_in_box_idx   = pts_in_box_idx2;
    *pts_in_box_g_num = pts_in_box_g_num2;
    *pts_in_box_coord = pts_in_box_coord2;
  }

  if (_dbbt->btShared != NULL) {
    free (recv_g_num);
    free (recv_coord);
  }
}




#undef _MIN
#undef _MAX

#ifdef __cplusplus
}
#endif /* __cplusplus */
