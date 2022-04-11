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
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_part.h"
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_extract_part.h"
#include "pdm_extract_part_priv.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_para_graph_dual.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_plane.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_face_center_2d
(
  int       n_part_in,
  int      *n_extract,
  int     **extract_lnum,
  int     **pface_edge_idx,
  int     **pface_edge,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***face_center
)
{
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    assert(pface_edge    [i_part] != NULL);
    assert(pface_edge_idx[i_part] != NULL);
    assert(pedge_vtx     [i_part] != NULL);
    assert(pvtx_coord    [i_part] != NULL);
  }

  double** entity_center = malloc(n_part_in * sizeof(double * ));
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

    int    *_pface_edge     = pface_edge    [i_part];
    int    *_pface_edge_idx = pface_edge_idx[i_part];
    int    *_pedge_vtx      = pedge_vtx     [i_part];
    double *_pvtx_coord     = pvtx_coord    [i_part];

    for(int idx_face = 0; idx_face < n_extract[i_part]; ++idx_face) {

      int i_face = extract_lnum[i_part][idx_face];
      entity_center[i_part][3*idx_face  ] = 0.;
      entity_center[i_part][3*idx_face+1] = 0.;
      entity_center[i_part][3*idx_face+2] = 0.;

      double inv = 1./((double) _pface_edge_idx[i_face+1] - _pface_edge_idx[i_face]);

      for(int idx_edge = _pface_edge_idx[i_face]; idx_edge < _pface_edge_idx[i_face+1]; ++idx_edge) {
        int i_edge = PDM_ABS(_pface_edge[idx_edge])-1;
        int i_vtx1 = _pedge_vtx[2*i_edge  ] - 1;
        int i_vtx2 = _pedge_vtx[2*i_edge+1] - 1;

        entity_center[i_part][3*idx_face  ] += 0.5 * (_pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
        entity_center[i_part][3*idx_face+1] += 0.5 * (_pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
        entity_center[i_part][3*idx_face+2] += 0.5 * (_pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);

      }
      entity_center[i_part][3*idx_face  ] = entity_center[i_part][3*idx_face  ] * inv;
      entity_center[i_part][3*idx_face+1] = entity_center[i_part][3*idx_face+1] * inv;
      entity_center[i_part][3*idx_face+2] = entity_center[i_part][3*idx_face+2] * inv;
    }
  }

  *face_center = entity_center;
}


static
void
_cell_center_3d
(
  int       n_part_in,
  int      *n_extract,
  int     **extract_lnum,
  int     **pcell_face_idx,
  int     **pcell_face,
  int     **pface_edge_idx,
  int     **pface_edge,
  int     **pface_vtx_idx,
  int     **pface_vtx,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***cell_center
)
{
  int from_edge = 0;
  int from_face = 0;
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    if(pface_edge    [i_part] != NULL) {
      from_edge = 1;
    }
    if(pface_vtx    [i_part] != NULL) {
      from_face = 1;
    }
    assert(pvtx_coord    [i_part] != NULL);
  }

  double** entity_center = malloc(n_part_in * sizeof(double * ));

  if(from_face == 1) {
    for(int i_part = 0; i_part < n_part_in; ++i_part) {
      entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_vtx      = pface_vtx     [i_part];
      int    *_pface_vtx_idx  = pface_vtx_idx [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      // PDM_log_trace_array_int(extract_lnum[i_part], n_extract[i_part], "extract_lnum ::");
      for(int idx_cell = 0; idx_cell < n_extract[i_part]; ++idx_cell) {
        int i_cell = extract_lnum[i_part][idx_cell];

        entity_center[i_part][3*idx_cell  ] = 0.;
        entity_center[i_part][3*idx_cell+1] = 0.;
        entity_center[i_part][3*idx_cell+2] = 0.;

        double inv = 1./((double)  _pcell_face_idx[idx_cell+1] - _pcell_face_idx[idx_cell]);

        double fcx = 0;
        double fcy = 0;
        double fcz = 0;
        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double inv2 = 1./((double)  _pface_vtx_idx[i_face+1] - _pface_vtx_idx[i_face]);

          for(int idx_vtx = _pface_vtx_idx[i_face]; idx_vtx < _pface_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = _pface_vtx[idx_vtx]-1;
            fcx += _pvtx_coord[3*i_vtx  ];
            fcy += _pvtx_coord[3*i_vtx+1];
            fcz += _pvtx_coord[3*i_vtx+2];
          }
          fcx = fcx * inv2;
          fcy = fcy * inv2;
          fcz = fcz * inv2;

          entity_center[i_part][3*idx_cell  ] += fcx;
          entity_center[i_part][3*idx_cell+1] += fcy;
          entity_center[i_part][3*idx_cell+2] += fcz;
        }

        entity_center[i_part][3*idx_cell  ] = entity_center[i_part][3*idx_cell  ] * inv;
        entity_center[i_part][3*idx_cell+1] = entity_center[i_part][3*idx_cell+1] * inv;
        entity_center[i_part][3*idx_cell+2] = entity_center[i_part][3*idx_cell+2] * inv;
      } /* End cell */
    }
  } else if( from_edge == 1) {
    for(int i_part = 0; i_part < n_part_in; ++i_part) {
      entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_edge     = pface_edge    [i_part];
      int    *_pface_edge_idx = pface_edge_idx[i_part];
      int    *_pedge_vtx      = pedge_vtx     [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      for(int idx_cell = 0; idx_cell < n_extract[i_part]; ++idx_cell) {
        int i_cell = extract_lnum[i_part][idx_cell];

        entity_center[i_part][3*idx_cell  ] = 0.;
        entity_center[i_part][3*idx_cell+1] = 0.;
        entity_center[i_part][3*idx_cell+2] = 0.;

        double inv = 1./((double)  _pcell_face_idx[idx_cell+1] - _pcell_face_idx[idx_cell]);

        double fcx = 0;
        double fcy = 0;
        double fcz = 0;
        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double inv2 = 1./((double)  _pface_edge_idx[i_face+1] - _pface_edge_idx[i_face]);

          for(int idx_edge = _pface_edge_idx[i_face]; idx_edge < _pface_edge_idx[i_face+1]; ++idx_edge) {
            int i_edge = PDM_ABS(_pface_edge[idx_edge])-1;
            int i_vtx1 = _pedge_vtx[2*i_edge  ] - 1;
            int i_vtx2 = _pedge_vtx[2*i_edge+1] - 1;
            fcx += 0.5 * (_pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
            fcy += 0.5 * (_pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
            fcz += 0.5 * (_pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);
          }
          fcx = fcx * inv2;
          fcy = fcy * inv2;
          fcz = fcz * inv2;

          entity_center[i_part][3*idx_cell  ] += fcx;
          entity_center[i_part][3*idx_cell+1] += fcy;
          entity_center[i_part][3*idx_cell+2] += fcz;
        }

        entity_center[i_part][3*idx_cell  ] = entity_center[i_part][3*idx_cell  ] * inv;
        entity_center[i_part][3*idx_cell+1] = entity_center[i_part][3*idx_cell+1] * inv;
        entity_center[i_part][3*idx_cell+2] = entity_center[i_part][3*idx_cell+2] * inv;
      } /* End cell */
    }
  }

  *cell_center = entity_center;
}



static
void
extract_entity1_entity2
(
int            n_part,
int           *n_entity1,
int           *n_entity2,
int           *n_extract_entity1,
int          **extract_entity1_lnum,
int          **entity1_entity2_idx,
int          **entity1_entity2,
PDM_g_num_t  **entity2_ln_to_gn,
int         ***extract_entity1_entity2_n,
PDM_g_num_t ***extract_entity1_entity2,
int          **n_extract_entity2,
int         ***extract_entity2_lnum,
int         ***idx_visited
)
{
  int         **_extract_entity1_entity2_n = malloc(n_part * sizeof(int         *));
  int         **_extract_entity2_lnum      = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **_extract_entity1_entity2   = malloc(n_part * sizeof(PDM_g_num_t *));
  int          *_n_extract_entity2         = malloc(n_part * sizeof(int          ));
  int         **_idx_visited               = malloc(n_part * sizeof(int         *));

  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_extract = n_extract_entity1[i_part];
    int         *_pentity1_entity2     = entity1_entity2    [i_part];
    int         *_pentity1_entity2_idx = entity1_entity2_idx[i_part];
    PDM_g_num_t *_pentity2_ln_to_gn    = entity2_ln_to_gn   [i_part];
    int          _pn_entity1           = n_entity1          [i_part];
    int          _pn_entity2           = n_entity2          [i_part];

    _extract_entity1_entity2  [i_part] = (PDM_g_num_t  *) malloc( _pentity1_entity2_idx[_pn_entity1] * sizeof(PDM_g_num_t));
    _extract_entity1_entity2_n[i_part] = (int          *) malloc( (n_extract)     * sizeof(int        ));
    _extract_entity2_lnum     [i_part] = (int          *) malloc( (_pn_entity2+1) * sizeof(int        ));
    _idx_visited              [i_part] = (int          *) malloc( _pn_entity2     * sizeof(int        ));

    int         *is_visited = (int *) malloc( _pn_entity2 * sizeof(int));
    for(int i = 0; i < _pn_entity2; ++i) {
      is_visited[i] = 0;
    }

    int idx_write = 0;
    _n_extract_entity2[i_part] = 0;
    for(int idx_entity = 0; idx_entity < n_extract; ++idx_entity) {
      int i_entity = extract_entity1_lnum[i_part][idx_entity];

      int n_tmp = _pentity1_entity2_idx[i_entity+1] - _pentity1_entity2_idx[i_entity];
      _extract_entity1_entity2_n[i_part][idx_entity] = n_tmp;

      for(int idx_entity2 = _pentity1_entity2_idx[i_entity]; idx_entity2 < _pentity1_entity2_idx[i_entity+1]; ++idx_entity2) {
        int i_entity2 = PDM_ABS(_pentity1_entity2[idx_entity2])-1;
        int sgn       = PDM_SIGN(_pentity1_entity2[idx_entity2]);
        _extract_entity1_entity2[i_part][idx_write++] = sgn*_pentity2_ln_to_gn[i_entity2];

        if(is_visited[i_entity2] == 0) {
          int idx = _n_extract_entity2[i_part]++;
          _extract_entity2_lnum[i_part][idx] = i_entity2;
          is_visited[i_entity2] = 1;
          _idx_visited[i_part][idx] = idx_write-1;
        }
      }
    }
    _extract_entity1_entity2  [i_part] = realloc(_extract_entity1_entity2[i_part], idx_write                  * sizeof(PDM_g_num_t));
    _extract_entity2_lnum     [i_part] = realloc(_extract_entity2_lnum   [i_part], _n_extract_entity2[i_part] * sizeof(int        ));
    _idx_visited              [i_part] = realloc(_idx_visited            [i_part], _n_extract_entity2[i_part] * sizeof(int        ));

    free(is_visited);
  }

  *n_extract_entity2           = _n_extract_entity2;
  *extract_entity2_lnum        = _extract_entity2_lnum;
  *extract_entity1_entity2_n   = _extract_entity1_entity2_n;
  *extract_entity1_entity2     = _extract_entity1_entity2;
  *idx_visited                 = _idx_visited;
}


static
PDM_part_to_block_t*
extract_and_renum_entity1_entity2
(
PDM_MPI_Comm           comm,
PDM_part_to_block_t   *ptb_entity1,
int                    n_part,
int                   *n_entity1,
int                   *n_entity2,
int                   *n_extract_entity1,
int                  **extract_entity1_lnum,
int                  **entity1_entity2_idx,
int                  **entity1_entity2,
PDM_g_num_t          **entity2_ln_to_gn,
int                  **dequi_entity1_entity2_idx,
PDM_g_num_t          **dequi_entity1_entity2,
PDM_g_num_t          **dequi_extract_parent_entity2_ln_to_gn,
int                  **n_extract_entity2,
int                 ***extract_entity2_lnum
)
{
  int         **extract_entity1_entity2_idx  = NULL;
  int         **_selected_entity1_entity2_n  = NULL;
  PDM_g_num_t **_selected_entity1_entity2    = NULL;
  extract_entity1_entity2(n_part,
                          n_entity1,
                          n_entity2,
                          n_extract_entity1,
                          extract_entity1_lnum,
                          entity1_entity2_idx,
                          entity1_entity2,
                          entity2_ln_to_gn,
                          &_selected_entity1_entity2_n,
                          &_selected_entity1_entity2,
                          n_extract_entity2,
                          extract_entity2_lnum,
                          &extract_entity1_entity2_idx);


  int *_n_extract_entity2 = *n_extract_entity2;

  PDM_gen_gnum_t* gnum_extract_entity2 = PDM_gnum_create(3,
                                                         n_part,
                                                         PDM_FALSE,
                                                         1.e-6,
                                                         comm,
                                                         PDM_OWNERSHIP_USER);
  //
  int *n_connect_tot = (int *) malloc( n_part * sizeof(int));
  int **sgn_entity1_entity2 = (int **) malloc( n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    n_connect_tot[i_part] = 0;
    for(int i = 0; i < n_extract_entity1[i_part]; ++i) {
      n_connect_tot[i_part] += _selected_entity1_entity2_n[i_part][i];
    }
    // TODO --> compute idx

    sgn_entity1_entity2[i_part] = (int *) malloc(sizeof(int) * n_connect_tot[i_part]);
    for(int i = 0; i < n_connect_tot[i_part]; ++i) {
      sgn_entity1_entity2[i_part][i]       = PDM_SIGN(_selected_entity1_entity2[i_part][i]);
      _selected_entity1_entity2[i_part][i] = PDM_ABS(_selected_entity1_entity2[i_part][i]);
    }

    PDM_gnum_set_from_parents(gnum_extract_entity2, i_part, n_connect_tot[i_part], _selected_entity1_entity2[i_part]);
  }

  PDM_gnum_compute(gnum_extract_entity2);

  PDM_g_num_t **_child_entity1_entity2           = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **_child_entity2_ln_to_gn          = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **_extract_parent_entity2_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    _child_entity1_entity2[i_part] = PDM_gnum_get(gnum_extract_entity2, i_part);

    if(0 == 1) {
      PDM_log_trace_array_long(_child_entity1_entity2[i_part], n_connect_tot[i_part], "_child_entity1_entity2 :: ");
      PDM_log_trace_array_int(extract_entity1_entity2_idx[i_part], _n_extract_entity2[i_part], "extract_entity1_entity2_idx :: ");
    }

    // Reforme a temporary ln_to_gn to setup recurence
    _child_entity2_ln_to_gn         [i_part] = (PDM_g_num_t *) malloc( _n_extract_entity2[i_part] * sizeof(PDM_g_num_t));
    _extract_parent_entity2_ln_to_gn[i_part] = (PDM_g_num_t *) malloc( _n_extract_entity2[i_part] * sizeof(PDM_g_num_t));
    for(int i = 0; i < _n_extract_entity2[i_part]; ++i) {
      int idx = extract_entity1_entity2_idx[i_part][i];
      _child_entity2_ln_to_gn         [i_part][i] = _child_entity1_entity2   [i_part][idx];
      _extract_parent_entity2_ln_to_gn[i_part][i] = _selected_entity1_entity2[i_part][idx]; // To exchange to keep link with original part
    }
    for(int i = 0; i < n_connect_tot[i_part]; ++i) {
      _child_entity1_entity2[i_part][i] *= sgn_entity1_entity2[i_part][i];
    }

    free(extract_entity1_entity2_idx[i_part]);
    free(sgn_entity1_entity2[i_part]);
  }
  free(sgn_entity1_entity2);
  free(extract_entity1_entity2_idx);
  free(n_connect_tot);
  PDM_gnum_free(gnum_extract_entity2);

  int         *_dequi_entity1_entity2_n = NULL;
  PDM_g_num_t *_dequi_entity1_entity2   = NULL;
  int dn_entity1_entity2 = PDM_part_to_block_exch(ptb_entity1,
                                                  sizeof(PDM_g_num_t),
                                                  PDM_STRIDE_VAR_INTERLACED,
                                                  -1,
                                                  _selected_entity1_entity2_n,
                                   (void **)      _child_entity1_entity2,
                                                  &_dequi_entity1_entity2_n,
                                   (void **)      &_dequi_entity1_entity2);
  PDM_UNUSED(dn_entity1_entity2);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(_child_entity1_entity2 [i_part]);
  }
  free(_child_entity1_entity2 );


  PDM_part_to_block_t *ptb_entity2_equi = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                   PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                   1.,
                                                                   _child_entity2_ln_to_gn,
                                                                   NULL,
                                                                   _n_extract_entity2,
                                                                   n_part,
                                                                   comm);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(_child_entity2_ln_to_gn    [i_part]);
    free(_selected_entity1_entity2_n[i_part]);
    free(_selected_entity1_entity2  [i_part]);
  }
  free(_child_entity2_ln_to_gn);
  free(_selected_entity1_entity2_n);
  free(_selected_entity1_entity2  );

  int          dn_entity2_equi        = PDM_part_to_block_n_elt_block_get(ptb_entity2_equi);
  PDM_g_num_t *dextract_entity2_gnum  = PDM_part_to_block_block_gnum_get (ptb_entity2_equi);

  if(0 == 1) {
    PDM_log_trace_array_long(dextract_entity2_gnum, dn_entity2_equi, "dextract_entity2_gnum :: ");
  }


  /*
   *  Exchange parent to keep link
   */
  PDM_g_num_t *_dequi_extract_parent_entity2_ln_to_gn = NULL;
  PDM_part_to_block_exch(ptb_entity2_equi,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
          (void **)      _extract_parent_entity2_ln_to_gn,
                         NULL,
          (void **)      &_dequi_extract_parent_entity2_ln_to_gn);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(_extract_parent_entity2_ln_to_gn[i_part]);
  }
  free(_extract_parent_entity2_ln_to_gn);

  // PDM_part_to_block_free(ptb_entity2_equi);

  int dn_equi_entity1 = PDM_part_to_block_n_elt_block_get(ptb_entity1);
  int *_dequi_entity1_entity2_idx = (int * ) malloc((dn_equi_entity1 + 1)  * sizeof(int));
  _dequi_entity1_entity2_idx[0] = 0;
  for(int i = 0; i < dn_equi_entity1; ++i) {
    _dequi_entity1_entity2_idx[i+1] = _dequi_entity1_entity2_idx[i] + _dequi_entity1_entity2_n[i];
  }

  free(_dequi_entity1_entity2_n);

  *dequi_entity1_entity2_idx             = _dequi_entity1_entity2_idx;
  *dequi_entity1_entity2                 = _dequi_entity1_entity2;
  *dequi_extract_parent_entity2_ln_to_gn = _dequi_extract_parent_entity2_ln_to_gn;


  return ptb_entity2_equi;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_extract_part_t*
PDM_extract_part_create
(
 const int                    dim,
 const int                    n_part_in,
 const int                    n_part_out,
       PDM_split_dual_t       split_dual_method,
       PDM_ownership_t        ownership,
       PDM_MPI_Comm           comm
)
{
  PDM_extract_part_t *extrp = (PDM_extract_part_t *) malloc(sizeof(PDM_extract_part_t));

  extrp->dim               = dim;
  extrp->n_part_in         = n_part_in;
  extrp->n_part_out        = n_part_out;
  extrp->split_dual_method = split_dual_method;
  extrp->ownership         = ownership;
  extrp->comm              = comm;

  extrp->n_cell         = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_face         = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_edge         = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_vtx          = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_extract      = (int          *) malloc(n_part_in * sizeof(int          ));

  extrp->pcell_face     = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pcell_face_idx = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pface_edge     = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pface_edge_idx = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pedge_vtx      = (int         **) malloc(n_part_in * sizeof(int         *));

  extrp->extract_lnum   = (int         **) malloc(n_part_in * sizeof(int         *));

  extrp->cell_ln_to_gn  = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));
  extrp->face_ln_to_gn  = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));
  extrp->edge_ln_to_gn  = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));
  extrp->vtx_ln_to_gn   = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));

  extrp->pface_vtx_idx  = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pface_vtx      = (int         **) malloc(n_part_in * sizeof(int         *));

  extrp->pvtx_coord     = (double      **) malloc(n_part_in * sizeof(double      *));

  extrp->dn_equi_cell               = 0;
  extrp->dn_equi_face               = 0;
  extrp->dn_equi_edge               = 0;
  extrp->dn_equi_vtx                = 0;

  extrp->dequi_cell_face            = NULL;
  extrp->dequi_cell_face_idx        = NULL;
  extrp->dequi_face_edge            = NULL;
  extrp->dequi_face_edge_idx        = NULL;
  extrp->dequi_edge_vtx             = NULL;
  extrp->dequi_edge_vtx_idx         = NULL;
  extrp->dequi_face_vtx             = NULL;
  extrp->dequi_face_vtx_idx         = NULL;

  extrp->dequi_parent_cell_ln_to_gn = NULL;
  extrp->dequi_parent_face_ln_to_gn = NULL;
  extrp->dequi_parent_edge_ln_to_gn = NULL;
  extrp->dequi_parent_vtx_ln_to_gn  = NULL;
  extrp->dequi_vtx_coord            = NULL;

  extrp->ptb_equi_cell              = NULL;
  extrp->ptb_equi_face              = NULL;
  extrp->ptb_equi_edge              = NULL;
  extrp->ptb_equi_vtx               = NULL;

  extrp->is_owner_connectivity    = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t   ) );
  extrp->is_owner_ln_to_gn        = malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_bool_t   ) );
  extrp->is_owner_parent_ln_to_gn = malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_bool_t   ) );
  extrp->is_owner_vtx_coord       = PDM_TRUE;

  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    extrp->is_owner_connectivity    [i] = PDM_TRUE;
    extrp->pextract_connectivity    [i] = NULL;
    extrp->pextract_connectivity_idx[i] = NULL;
  }
  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    extrp->is_owner_ln_to_gn              [i] = PDM_TRUE;
    extrp->is_owner_parent_ln_to_gn       [i] = PDM_TRUE;
    extrp->pextract_entity_ln_to_gn       [i] = NULL;
    extrp->pextract_entity_parent_ln_to_gn[i] = NULL;
  }

  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    extrp->pextract_n_entity[i] = NULL;
  }

  return extrp;
}



void
PDM_extract_part_compute
(
  PDM_extract_part_t        *extrp
)
{
  assert(extrp->dim >= 2);

  int          *pn_entity    = 0;
  PDM_g_num_t **entity_g_num = NULL;
  if(extrp->dim == 3) {
    pn_entity    = extrp->n_cell;
    entity_g_num = extrp->cell_ln_to_gn;
  } else {
    pn_entity    = extrp->n_face;
    entity_g_num = extrp->face_ln_to_gn;
  }

  PDM_UNUSED(pn_entity);

  /*
   *  Create array selected in gnum
   */
  PDM_g_num_t** entity_extract_g_num = (PDM_g_num_t **) malloc( extrp->n_part_in * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    entity_extract_g_num[i_part] = (PDM_g_num_t *) malloc( extrp->n_extract[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity = 0; i_entity < extrp->n_extract[i_part]; ++i_entity) {
      entity_extract_g_num[i_part][i_entity] = entity_g_num[i_part][extrp->extract_lnum[i_part][i_entity]];
    }
    if(0 == 1) {
      PDM_log_trace_array_long(entity_extract_g_num[i_part], extrp->n_extract[i_part], "entity_extract_g_num ::" );
    }
  }

  PDM_gen_gnum_t* gnum_extract = PDM_gnum_create(3, extrp->n_part_in, PDM_FALSE,
                                                 1.e-6,
                                                 extrp->comm,
                                                 PDM_OWNERSHIP_KEEP);

  /*
   * Calcul des coordonnées to setup hilbert ordering (independant of parallelism )
   */
  double **entity_center = NULL;
  if(extrp->split_dual_method == PDM_SPLIT_DUAL_WITH_HILBERT) {
    if(extrp->dim == 2) {
      // Compute entity_center with face_edge + edge_vtx
      _face_center_2d(extrp->n_part_in,
                      extrp->n_extract,
                      extrp->extract_lnum,
                      extrp->pface_edge_idx,
                      extrp->pface_edge,
                      extrp->pedge_vtx,
                      extrp->pvtx_coord,
                      &entity_center);
    } else {  // dim == 3
      _cell_center_3d(extrp->n_part_in,
                      extrp->n_extract,
                      extrp->extract_lnum,
                      extrp->pcell_face_idx,
                      extrp->pcell_face,
                      extrp->pface_edge_idx,
                      extrp->pface_edge,
                      extrp->pface_vtx_idx,
                      extrp->pface_vtx,
                      extrp->pedge_vtx,
                      extrp->pvtx_coord,
                      &entity_center);
    }

    for (int i_part = 0; i_part < extrp->n_part_in; i_part++){
      PDM_gnum_set_from_coords(gnum_extract, i_part, extrp->n_extract[i_part], entity_center[i_part], NULL);
    }
  } else {

    for (int i_part = 0; i_part < extrp->n_part_in; i_part++){
      PDM_gnum_set_from_parents(gnum_extract, i_part, extrp->n_extract[i_part], entity_extract_g_num[i_part]);
    }

  }

  /*
   * Global numering computation
   */
  PDM_g_num_t **child_selected_g_num = (PDM_g_num_t **) malloc( extrp->n_part_in * sizeof(PDM_g_num_t *));
  PDM_gnum_compute(gnum_extract);

  for (int i_part = 0; i_part < extrp->n_part_in; i_part++){
    child_selected_g_num[i_part] = PDM_gnum_get(gnum_extract, i_part);
    // PDM_log_trace_array_long(child_selected_g_num[i_part], extrp->n_extract[i_part], "child_selected_g_num : ");
  }

  /*
   *  Remake equilibrate block -> Block is not partial because we use child_gnum
   */
  PDM_part_to_block_t *ptb_equi = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                           1.,
                                                           child_selected_g_num,
                                                           NULL,
                                                           extrp->n_extract,
                                                           extrp->n_part_in,
                                                           extrp->comm);

  free(child_selected_g_num);
  PDM_gnum_free(gnum_extract);

  PDM_g_num_t *dequi_parent_entity_ln_to_gn = NULL;
  PDM_part_to_block_exch(ptb_equi,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
          (void **)      entity_extract_g_num,
                         NULL,
          (void **)      &dequi_parent_entity_ln_to_gn);


  int          dn_entity_equi = PDM_part_to_block_n_elt_block_get(ptb_equi);
  // PDM_g_num_t *dextract_gnum  = PDM_part_to_block_block_gnum_get (ptb_equi);

  if(extrp->split_dual_method == PDM_SPLIT_DUAL_WITH_HILBERT) {
    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      free(entity_center[i_part]);
    }
    free(entity_center);
  } else {
    /*
     * Si scotch ou metis, on doit calculer le dcell_face + pdm_gnum (pour avoir les faces child)
     *   Puis dconnectiviy_transpose puis combine
     */
    // abort();
  }

  if(extrp->dim == 3) {
    extrp->dequi_parent_cell_ln_to_gn = dequi_parent_entity_ln_to_gn;
  } else {
    extrp->dequi_parent_face_ln_to_gn = dequi_parent_entity_ln_to_gn;
  }

  /*
   * Extraction des connectivités
   */
  int from_face_edge = 0;
  int from_face_vtx  = 0;
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    if(extrp->pface_edge    [i_part] != NULL) {
      from_face_edge = 1;
    }
    if(extrp->pface_vtx    [i_part] != NULL) {
      from_face_vtx = 1;
    }
  }
  int                  *n_extract_vtx    = NULL;
  int                 **extract_vtx_lnum = NULL;
  if(extrp->dim == 3) {

    int          *n_extract_face              = NULL;
    int         **extract_face_lnum           = NULL;
    extrp->ptb_equi_face = extract_and_renum_entity1_entity2(extrp->comm,
                                                             ptb_equi,
                                                             extrp->n_part_in,
                                                             extrp->n_cell,
                                                             extrp->n_face,
                                                             extrp->n_extract,
                                                             extrp->extract_lnum,
                                                             extrp->pcell_face_idx,
                                                             extrp->pcell_face,
                                                             extrp->face_ln_to_gn,
                                                             &extrp->dequi_cell_face_idx,
                                                             &extrp->dequi_cell_face,
                                                             &extrp->dequi_parent_face_ln_to_gn,
                                                             &n_extract_face,
                                                             &extract_face_lnum);


    extrp->dn_equi_face = PDM_part_to_block_n_elt_block_get(extrp->ptb_equi_face);
    if(0 == 1) {
      PDM_log_trace_array_int (extrp->dequi_cell_face_idx, dn_entity_equi+1                          , "dequi_cell_face_idx :: ");
      PDM_log_trace_array_long(extrp->dequi_cell_face    , extrp->dequi_cell_face_idx[dn_entity_equi], "dequi_cell_face     :: ");
    }

    /*
     *  If face_vtx available descent by it
     */
    if(from_face_edge == 1) {
      int                  *n_extract_edge    = NULL;
      int                 **extract_edge_lnum = NULL;
      extrp->ptb_equi_edge = extract_and_renum_entity1_entity2(extrp->comm,
                                                              extrp->ptb_equi_face,
                                                              extrp->n_part_in,
                                                              extrp->n_face,
                                                              extrp->n_edge,
                                                              n_extract_face,
                                                              extract_face_lnum,
                                                              extrp->pface_edge_idx,
                                                              extrp->pface_edge,
                                                              extrp->edge_ln_to_gn,
                                                              &extrp->dequi_face_edge_idx,
                                                              &extrp->dequi_face_edge,
                                                              &extrp->dequi_parent_edge_ln_to_gn,
                                                              &n_extract_edge,
                                                              &extract_edge_lnum);


      extrp->dn_equi_edge = PDM_part_to_block_n_elt_block_get(extrp->ptb_equi_edge);

      int **pedge_vtx_idx = (int **) malloc(extrp->n_part_in * sizeof(int *));
      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        pedge_vtx_idx[i_part] = (int *) malloc((extrp->n_edge[i_part]+1) * sizeof(int));
        pedge_vtx_idx[i_part][0] = 0;
        for(int i_edge = 0; i_edge < extrp->n_edge[i_part]; ++i_edge) {
          pedge_vtx_idx[i_part][i_edge+1] = pedge_vtx_idx[i_part][i_edge] + 2;
        }
      }
      extrp->ptb_equi_vtx = extract_and_renum_entity1_entity2(extrp->comm,
                                                              extrp->ptb_equi_edge,
                                                              extrp->n_part_in,
                                                              extrp->n_edge,
                                                              extrp->n_vtx,
                                                              n_extract_edge,
                                                              extract_edge_lnum,
                                                              pedge_vtx_idx,
                                                              extrp->pedge_vtx,
                                                              extrp->vtx_ln_to_gn,
                                                              &extrp->dequi_edge_vtx_idx,
                                                              &extrp->dequi_edge_vtx,
                                                              &extrp->dequi_parent_vtx_ln_to_gn,
                                                              &n_extract_vtx,
                                                              &extract_vtx_lnum);

      extrp->dn_equi_vtx = PDM_part_to_block_n_elt_block_get(extrp->ptb_equi_vtx);
      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        free(extract_edge_lnum   [i_part]);
        free(pedge_vtx_idx   [i_part]);
      }
      free(extract_edge_lnum   );
      free(n_extract_edge      );
      free(pedge_vtx_idx      );

    } else if(from_face_vtx == 1){
      extrp->ptb_equi_vtx = extract_and_renum_entity1_entity2(extrp->comm,
                                                              extrp->ptb_equi_face,
                                                              extrp->n_part_in,
                                                              extrp->n_face,
                                                              extrp->n_vtx,
                                                              n_extract_face,
                                                              extract_face_lnum,
                                                              extrp->pface_vtx_idx,
                                                              extrp->pface_vtx,
                                                              extrp->vtx_ln_to_gn,
                                                              &extrp->dequi_face_vtx_idx,
                                                              &extrp->dequi_face_vtx,
                                                              &extrp->dequi_parent_vtx_ln_to_gn,
                                                              &n_extract_vtx,
                                                              &extract_vtx_lnum);


      extrp->dn_equi_vtx = PDM_part_to_block_n_elt_block_get(extrp->ptb_equi_vtx);
      if(0 == 1) {
        PDM_log_trace_array_int (extrp->dequi_face_vtx_idx, extrp->dn_equi_face+1                         , "dequi_face_vtx_idx :: ");
        PDM_log_trace_array_long(extrp->dequi_face_vtx    , extrp->dequi_face_vtx_idx[extrp->dn_equi_face], "dequi_face_vtx     :: ");
      }

      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        free(extract_face_lnum[i_part]);
      }
      free(extract_face_lnum);
      free(n_extract_face   );
    }
  } else { // dim == 2
    if(from_face_edge == 1) {
       int                  *n_extract_edge    = NULL;
       int                 **extract_edge_lnum = NULL;
       extrp->ptb_equi_edge = extract_and_renum_entity1_entity2(extrp->comm,
                                                                ptb_equi,
                                                                extrp->n_part_in,
                                                                extrp->n_face,
                                                                extrp->n_edge,
                                                                extrp->n_extract,
                                                                extrp->extract_lnum,
                                                                extrp->pface_edge_idx,
                                                                extrp->pface_edge,
                                                                extrp->edge_ln_to_gn,
                                                                &extrp->dequi_face_edge_idx,
                                                                &extrp->dequi_face_edge,
                                                                &extrp->dequi_parent_edge_ln_to_gn,
                                                                &n_extract_edge,
                                                                &extract_edge_lnum);

       extrp->dn_equi_edge = PDM_part_to_block_n_elt_block_get(extrp->ptb_equi_edge);

       int **pedge_vtx_idx = (int **) malloc(extrp->n_part_in * sizeof(int *));
       for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
         pedge_vtx_idx[i_part] = (int *) malloc((extrp->n_edge[i_part]+1) * sizeof(int));
         pedge_vtx_idx[i_part][0] = 0;
         for(int i_edge = 0; i_edge < extrp->n_edge[i_part]; ++i_edge) {
           pedge_vtx_idx[i_part][i_edge+1] = pedge_vtx_idx[i_part][i_edge] + 2;
         }
       }
       extrp->ptb_equi_vtx = extract_and_renum_entity1_entity2(extrp->comm,
                                                               extrp->ptb_equi_edge,
                                                               extrp->n_part_in,
                                                               extrp->n_edge,
                                                               extrp->n_vtx,
                                                               n_extract_edge,
                                                               extract_edge_lnum,
                                                               pedge_vtx_idx,
                                                               extrp->pedge_vtx,
                                                               extrp->vtx_ln_to_gn,
                                                               &extrp->dequi_edge_vtx_idx,
                                                               &extrp->dequi_edge_vtx,
                                                               &extrp->dequi_parent_vtx_ln_to_gn,
                                                               &n_extract_vtx,
                                                               &extract_vtx_lnum);

       extrp->dn_equi_vtx = PDM_part_to_block_n_elt_block_get(extrp->ptb_equi_vtx);
       for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
         free(extract_edge_lnum   [i_part]);
         free(pedge_vtx_idx   [i_part]);
       }
       free(extract_edge_lnum   );
       free(n_extract_edge      );
       free(pedge_vtx_idx      );

     } else if(from_face_vtx == 1){
       extrp->ptb_equi_vtx = extract_and_renum_entity1_entity2(extrp->comm,
                                                               ptb_equi,
                                                               extrp->n_part_in,
                                                               extrp->n_face,
                                                               extrp->n_vtx,
                                                               extrp->n_extract,
                                                               extrp->extract_lnum,
                                                               extrp->pface_vtx_idx,
                                                               extrp->pface_vtx,
                                                               extrp->vtx_ln_to_gn,
                                                               &extrp->dequi_face_vtx_idx,
                                                               &extrp->dequi_face_vtx,
                                                               &extrp->dequi_parent_vtx_ln_to_gn,
                                                               &n_extract_vtx,
                                                               &extract_vtx_lnum);


       extrp->dn_equi_vtx = PDM_part_to_block_n_elt_block_get(extrp->ptb_equi_vtx);
       if(0 == 1) {
         PDM_log_trace_array_int (extrp->dequi_face_vtx_idx, extrp->dn_equi_face+1                         , "(2D) dequi_face_vtx_idx :: ");
         PDM_log_trace_array_long(extrp->dequi_face_vtx    , extrp->dequi_face_vtx_idx[extrp->dn_equi_face], "(2D) dequi_face_vtx     :: ");
       }

     }
  }


  /*
   * Exchange coordinates
   */
  double **pextract_vtx_coord = (double **) malloc( extrp->n_part_in * sizeof(double *));
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    pextract_vtx_coord[i_part] = (double *) malloc( 3 * n_extract_vtx[i_part] * sizeof(double));

    for(int idx_vtx = 0; idx_vtx < n_extract_vtx[i_part]; ++idx_vtx) {
      int i_vtx = extract_vtx_lnum[i_part][idx_vtx];
      pextract_vtx_coord[i_part][3*idx_vtx  ] = extrp->pvtx_coord[i_part][3*i_vtx  ];
      pextract_vtx_coord[i_part][3*idx_vtx+1] = extrp->pvtx_coord[i_part][3*i_vtx+1];
      pextract_vtx_coord[i_part][3*idx_vtx+2] = extrp->pvtx_coord[i_part][3*i_vtx+2];
    }
  }

  PDM_part_to_block_exch(extrp->ptb_equi_vtx,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
          (void **)      pextract_vtx_coord,
                         NULL,
          (void **)      &extrp->dequi_vtx_coord);

  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    free(extract_vtx_lnum  [i_part]);
    free(pextract_vtx_coord[i_part]);
  }
  free(extract_vtx_lnum  );
  free(n_extract_vtx     );
  free(pextract_vtx_coord);

  int i_rank;
  PDM_MPI_Comm_rank(extrp->comm, &i_rank);

  /*
   * At this stage we have all information in block frame but we need to repart them
   *   - if hilbert -> block and part is the same (implicite block numbering is hilbert)
   *   - if scotch  -> Need to setup graph and split it before doing partitioning
   */
  if(extrp->split_dual_method == PDM_SPLIT_DUAL_WITH_HILBERT) {
    assert(extrp->n_part_out == 1);


    PDM_g_num_t* cell_distri = PDM_part_to_block_distrib_index_get(ptb_equi);
    int dn_cell_equi = cell_distri[i_rank+1] - cell_distri[i_rank];

    extrp->pextract_n_entity[PDM_MESH_ENTITY_CELL] = (int  *) malloc(extrp->n_part_out * sizeof(int  ));
    extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_CELL] = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));

    int i_part0 = 0;
    extrp->pextract_n_entity[PDM_MESH_ENTITY_CELL][i_part0] = dn_cell_equi;
    extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_CELL][i_part0] = (PDM_g_num_t *) malloc( extrp->pextract_n_entity[PDM_MESH_ENTITY_CELL][i_part0] * sizeof(PDM_g_num_t));

    for(int i = 0; i < dn_cell_equi; ++i) {
      extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_CELL][i_part0][i] = cell_distri[i_rank] + i + 1;
    }


  } else {

    int         *delmt_to_arc_idx = NULL; // Donc cell_face OU face_edge
    PDM_g_num_t *delmt_to_arc     = NULL;
    PDM_g_num_t *distrib_elmt     = PDM_part_to_block_distrib_index_get(ptb_equi);;
    PDM_g_num_t *distrib_arc      = NULL;

    if(extrp->dim == 3) {
      delmt_to_arc_idx = extrp->dequi_cell_face_idx;
      delmt_to_arc     = extrp->dequi_cell_face;
      distrib_arc      = PDM_part_to_block_distrib_index_get(extrp->ptb_equi_face);
    } else {
      if(extrp->ptb_equi_edge != NULL) {
        distrib_arc      = PDM_part_to_block_distrib_index_get(extrp->ptb_equi_edge);
        delmt_to_arc_idx = extrp->dequi_face_edge_idx;
        delmt_to_arc     = extrp->dequi_face_edge;
      } else {
        distrib_arc      = PDM_part_to_block_distrib_index_get(extrp->ptb_equi_vtx);
        delmt_to_arc_idx = extrp->dequi_face_vtx_idx;
        delmt_to_arc     = extrp->dequi_face_vtx;
      }
    }

    int         *darc_to_elmt_idx = NULL;
    PDM_g_num_t *darc_to_elmt     = NULL;
    PDM_dconnectivity_transpose(extrp->comm,
                                distrib_elmt,
                                distrib_arc,
                                delmt_to_arc_idx,
                                delmt_to_arc,
                                1,
                                &darc_to_elmt_idx,
                                &darc_to_elmt);

    PDM_g_num_t *dual_graph_idx = NULL;
    PDM_g_num_t *dual_graph     = NULL;
    PDM_deduce_combine_connectivity_dual(extrp->comm,
                                         distrib_elmt,
                                         distrib_arc,
                                         delmt_to_arc_idx,
                                         delmt_to_arc,
                                         darc_to_elmt_idx,
                                         darc_to_elmt,
                                         1,
                                         &dual_graph_idx,
                                         &dual_graph);
    free(darc_to_elmt_idx);
    free(darc_to_elmt);

    /* Shift to 0 dual */
    int dn_elmt = distrib_elmt[i_rank+1] - distrib_elmt[i_rank];
    int *_elmt_part = malloc(dn_elmt * sizeof(int));
    for(int i = 0; i < dual_graph_idx[dn_elmt]; ++i) {
      dual_graph[i] = dual_graph[i] - 1;
    }

    int tn_part;
    PDM_MPI_Allreduce(& extrp->n_part_out, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, extrp->comm);

    PDM_para_graph_split(extrp->split_dual_method,
                         distrib_elmt,
                         dual_graph_idx,
                         dual_graph,
                         NULL,
                         NULL,
                         tn_part,
                         NULL,
                         _elmt_part,
                         extrp->comm);

    PDM_g_num_t *distrib_partition = PDM_compute_entity_distribution(extrp->comm, extrp->n_part_out);

    if(extrp->dim == 3) {
      PDM_part_assemble_partitions(extrp->comm,
                                   distrib_partition,
                                   distrib_elmt,
                                   _elmt_part,
                                   &extrp->pextract_n_entity[PDM_MESH_ENTITY_CELL],
                                   &extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_CELL]);
    } else {
      PDM_part_assemble_partitions(extrp->comm,
                                   distrib_partition,
                                   distrib_elmt,
                                   _elmt_part,
                                   &extrp->pextract_n_entity[PDM_MESH_ENTITY_FACE],
                                   &extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_FACE]);
    }

    free(distrib_partition);
    free(dual_graph_idx);
    free(dual_graph);
    free(_elmt_part);

  }


  if(extrp->dim >= 2) {

    if(extrp->dim == 3) {
      PDM_g_num_t* cell_distri = PDM_part_to_block_distrib_index_get(ptb_equi);
      PDM_part_dconnectivity_to_pconnectivity_sort(extrp->comm,
                                                   cell_distri,
                                                   extrp->dequi_cell_face_idx,
                                                   extrp->dequi_cell_face,
                                                   extrp->n_part_out,
                                                   extrp->pextract_n_entity        [PDM_MESH_ENTITY_CELL],
                            (const PDM_g_num_t **) extrp->pextract_entity_ln_to_gn [PDM_MESH_ENTITY_CELL],
                                                  &extrp->pextract_n_entity        [PDM_MESH_ENTITY_FACE],
                                                  &extrp->pextract_entity_ln_to_gn [PDM_MESH_ENTITY_FACE],
                                                  &extrp->pextract_connectivity_idx[PDM_CONNECTIVITY_TYPE_CELL_FACE],
                                                  &extrp->pextract_connectivity    [PDM_CONNECTIVITY_TYPE_CELL_FACE]);

      PDM_part_dfield_to_pfield(extrp->comm,
                                extrp->n_part_out,
                                sizeof(PDM_g_num_t),
                                cell_distri,
         (unsigned char *)      extrp->dequi_parent_cell_ln_to_gn,
                                extrp->pextract_n_entity       [PDM_MESH_ENTITY_CELL],
         (const PDM_g_num_t **) extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_CELL],
        (unsigned char     ***) &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_CELL]);
    }

    if(extrp->dequi_face_edge_idx != NULL) {
      PDM_g_num_t* face_distri = PDM_part_to_block_distrib_index_get(extrp->ptb_equi_face);
      PDM_part_dconnectivity_to_pconnectivity_sort(extrp->comm,
                                                   face_distri,
                                                   extrp->dequi_face_edge_idx,
                                                   extrp->dequi_face_edge,
                                                   extrp->n_part_out,
                                                   extrp->pextract_n_entity       [PDM_MESH_ENTITY_FACE],
                            (const PDM_g_num_t **) extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_FACE],
                                                  &extrp->pextract_n_entity       [PDM_MESH_ENTITY_EDGE],
                                                  &extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_EDGE],
                                                  &extrp->pextract_connectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                                  &extrp->pextract_connectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE]);

      PDM_part_dfield_to_pfield(extrp->comm,
                                extrp->n_part_out,
                                sizeof(PDM_g_num_t),
                                face_distri,
         (unsigned char *)      extrp->dequi_parent_face_ln_to_gn,
                                extrp->pextract_n_entity       [PDM_MESH_ENTITY_FACE],
         (const PDM_g_num_t **) extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_FACE],
        (unsigned char     ***) &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_FACE]);

      PDM_g_num_t* edge_distri = PDM_part_to_block_distrib_index_get(extrp->ptb_equi_edge);
      int **pextract_edge_vtx_idx = NULL;
      PDM_part_dconnectivity_to_pconnectivity_sort(extrp->comm,
                                                   edge_distri,
                                                   extrp->dequi_edge_vtx_idx,
                                                   extrp->dequi_edge_vtx,
                                                   extrp->n_part_out,
                                                   extrp->pextract_n_entity       [PDM_MESH_ENTITY_EDGE],
                            (const PDM_g_num_t **) extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_EDGE],
                                                  &extrp->pextract_n_entity       [PDM_MESH_ENTITY_VERTEX],
                                                  &extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_VERTEX],
                                                  &pextract_edge_vtx_idx,
                                                  &extrp->pextract_connectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX]);
      for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
        free(pextract_edge_vtx_idx[i_part]);
      }
      free(pextract_edge_vtx_idx);

      PDM_part_dfield_to_pfield(extrp->comm,
                                extrp->n_part_out,
                                sizeof(PDM_g_num_t),
                                edge_distri,
         (unsigned char *)      extrp->dequi_parent_edge_ln_to_gn,
                                extrp->pextract_n_entity       [PDM_MESH_ENTITY_EDGE],
         (const PDM_g_num_t **) extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_EDGE],
        (unsigned char     ***) &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_EDGE]);

    } else {
      PDM_g_num_t* face_distri = NULL;
      if(extrp->ptb_equi_face != NULL) {
        face_distri = PDM_part_to_block_distrib_index_get(extrp->ptb_equi_face);
      } else {
        face_distri = PDM_part_to_block_distrib_index_get(ptb_equi);
      }
      PDM_part_dconnectivity_to_pconnectivity_sort(extrp->comm,
                                                   face_distri,
                                                   extrp->dequi_face_vtx_idx,
                                                   extrp->dequi_face_vtx,
                                                   extrp->n_part_out,
                                                   extrp->pextract_n_entity        [PDM_MESH_ENTITY_FACE],
                            (const PDM_g_num_t **) extrp->pextract_entity_ln_to_gn [PDM_MESH_ENTITY_FACE],
                                                  &extrp->pextract_n_entity        [PDM_MESH_ENTITY_VERTEX],
                                                  &extrp->pextract_entity_ln_to_gn [PDM_MESH_ENTITY_VERTEX],
                                                  &extrp->pextract_connectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                                  &extrp->pextract_connectivity    [PDM_CONNECTIVITY_TYPE_FACE_VTX]);

      PDM_part_dfield_to_pfield(extrp->comm,
                                extrp->n_part_out,
                                sizeof(PDM_g_num_t),
                                face_distri,
         (unsigned char *)      extrp->dequi_parent_face_ln_to_gn,
                                extrp->pextract_n_entity       [PDM_MESH_ENTITY_FACE],
         (const PDM_g_num_t **) extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_FACE],
        (unsigned char     ***) &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_FACE]);
    }

    /*
     * Finalize with vtx
     */
    PDM_g_num_t* vtx_distri = PDM_part_to_block_distrib_index_get(extrp->ptb_equi_vtx);
    PDM_part_dcoordinates_to_pcoordinates(extrp->comm,
                                          extrp->n_part_out,
                                          vtx_distri,
                                          extrp->dequi_vtx_coord,
                                          extrp->pextract_n_entity       [PDM_MESH_ENTITY_VERTEX],
                   (const PDM_g_num_t **) extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_VERTEX],
                                         &extrp->pextract_vtx_coord);


    PDM_part_dfield_to_pfield(extrp->comm,
                              extrp->n_part_out,
                              sizeof(PDM_g_num_t),
                              vtx_distri,
       (unsigned char *)      extrp->dequi_parent_vtx_ln_to_gn,
                              extrp->pextract_n_entity       [PDM_MESH_ENTITY_VERTEX],
       (const PDM_g_num_t **) extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_VERTEX],
      (unsigned char     ***) &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VERTEX]);

  }

  /*
   * Cleaning
   */
  if(extrp->dequi_cell_face != NULL ) {
    free(extrp->dequi_cell_face);
  }
  if(extrp->dequi_cell_face_idx != NULL ) {
    free(extrp->dequi_cell_face_idx);
  }
  if(extrp->dequi_face_edge != NULL ) {
    free(extrp->dequi_face_edge);
  }
  if(extrp->dequi_face_edge_idx != NULL ) {
    free(extrp->dequi_face_edge_idx);
  }
  if(extrp->dequi_edge_vtx != NULL ) {
    free(extrp->dequi_edge_vtx);
  }
  if(extrp->dequi_edge_vtx_idx != NULL ) {
    free(extrp->dequi_edge_vtx_idx);
  }
  if(extrp->dequi_face_vtx != NULL ) {
    free(extrp->dequi_face_vtx);
  }
  if(extrp->dequi_face_vtx_idx != NULL ) {
    free(extrp->dequi_face_vtx_idx);
  }
  if(extrp->dequi_parent_cell_ln_to_gn != NULL ) {
    free(extrp->dequi_parent_cell_ln_to_gn);
  }
  if(extrp->dequi_parent_face_ln_to_gn != NULL ) {
    free(extrp->dequi_parent_face_ln_to_gn);
  }
  if(extrp->dequi_parent_edge_ln_to_gn != NULL ) {
    free(extrp->dequi_parent_edge_ln_to_gn);
  }
  if(extrp->dequi_parent_vtx_ln_to_gn != NULL ) {
    free(extrp->dequi_parent_vtx_ln_to_gn);
  }
  if(extrp->dequi_vtx_coord != NULL ) {
    free(extrp->dequi_vtx_coord);
  }

  PDM_part_to_block_free(ptb_equi);
  if(extrp->ptb_equi_face != NULL) {
    PDM_part_to_block_free(extrp->ptb_equi_face);
  }
  if(extrp->ptb_equi_edge != NULL) {
    PDM_part_to_block_free(extrp->ptb_equi_edge);
  }
  if(extrp->ptb_equi_vtx != NULL) {
    PDM_part_to_block_free(extrp->ptb_equi_vtx);
  }

  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    free(entity_extract_g_num[i_part]);
  }
  free(entity_extract_g_num);

  /*
   * Exchange parent (dommage on refais les block_to_part du dconnectivity_to_pconnectivity)
   */

}


void
PDM_extract_part_part_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
)
{
  extrp->n_cell        [i_part] = n_cell;
  extrp->n_face        [i_part] = n_face;
  extrp->n_edge        [i_part] = n_edge;
  extrp->n_vtx         [i_part] = n_vtx;
  extrp->pcell_face    [i_part] = cell_face;
  extrp->pcell_face_idx[i_part] = cell_face_idx;
  extrp->pface_edge    [i_part] = face_edge;
  extrp->pface_edge_idx[i_part] = face_edge_idx;
  extrp->pedge_vtx     [i_part] = edge_vtx;
  extrp->cell_ln_to_gn [i_part] = cell_ln_to_gn;
  extrp->face_ln_to_gn [i_part] = face_ln_to_gn;
  extrp->edge_ln_to_gn [i_part] = edge_ln_to_gn;
  extrp->vtx_ln_to_gn  [i_part] = vtx_ln_to_gn;
  extrp->pface_vtx_idx [i_part] = face_vtx_idx;
  extrp->pface_vtx     [i_part] = face_vtx;
  extrp->pvtx_coord    [i_part] = vtx_coord;
}

void
PDM_extract_part_selected_lnum_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_extract,
  int                      *extract_lnum
)
{
  extrp->n_extract   [i_part] = n_extract;
  extrp->extract_lnum[i_part] = extract_lnum;
}


void
PDM_extract_part_n_entity_get
(
 PDM_extract_part_t       *extrp,
 PDM_mesh_entities_t       entity_type,
 int                     **pn_entity
)
{
  *pn_entity = malloc( extrp->n_part_out * sizeof(int));
  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    (*pn_entity)[i_part] = extrp->pextract_n_entity[entity_type][i_part];
  }
}


void
PDM_extract_part_connectivity_get
(
 PDM_extract_part_t        *extrp,
 PDM_connectivity_type_t    connectivity_type,
 int                     ***connect,
 int                     ***connect_idx,
 PDM_ownership_t           ownership
)
{
  *connect     = extrp->pextract_connectivity    [connectivity_type];
  *connect_idx = extrp->pextract_connectivity_idx[connectivity_type];

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    extrp->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else {
    extrp->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }
}




void
PDM_extract_part_ln_to_gn_get
(
 PDM_extract_part_t        *extrp,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t            ***pentity_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  *pentity_ln_to_gn = extrp->pextract_entity_ln_to_gn[entity_type];
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    extrp->is_owner_ln_to_gn[entity_type] = PDM_FALSE;
  } else {
    extrp->is_owner_ln_to_gn[entity_type] = PDM_TRUE;
  }
}


void
PDM_extract_part_parent_ln_to_gn_get
(
 PDM_extract_part_t        *extrp,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t            ***parent_entity_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  *parent_entity_ln_to_gn = extrp->pextract_entity_parent_ln_to_gn[entity_type];
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    extrp->is_owner_parent_ln_to_gn[entity_type] = PDM_FALSE;
  } else {
    extrp->is_owner_parent_ln_to_gn[entity_type] = PDM_TRUE;
  }
}

void
PDM_extract_part_vtx_coord_get
(
 PDM_extract_part_t         *extrp,
 double                   ***pvtx_coord,
 PDM_ownership_t           ownership
)
{
  *pvtx_coord = extrp->pextract_vtx_coord;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    extrp->is_owner_vtx_coord = PDM_FALSE;
  } else {
    extrp->is_owner_vtx_coord = PDM_TRUE;
  }
}


void
PDM_extract_part_free
(
  PDM_extract_part_t  *extrp
)
{
  free(extrp->n_cell        );
  free(extrp->n_face        );
  free(extrp->n_edge        );
  free(extrp->n_vtx         );
  free(extrp->n_extract     );

  free(extrp->pcell_face    );
  free(extrp->pcell_face_idx);

  free(extrp->pface_edge    );
  free(extrp->pface_edge_idx);

  free(extrp->pface_vtx     );
  free(extrp->pface_vtx_idx );
  free(extrp->pedge_vtx     );
  free(extrp->extract_lnum  );

  free(extrp->cell_ln_to_gn );
  free(extrp->face_ln_to_gn );
  free(extrp->edge_ln_to_gn );
  free(extrp->vtx_ln_to_gn  );

  free(extrp->pvtx_coord);

  /*
   * Free extracted partition if owner
   */
  /* Free connectivity */
  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    if(extrp->is_owner_connectivity[i] == PDM_TRUE) {
      if(extrp->pextract_connectivity[i] != NULL) {
        for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
          if(extrp->pextract_connectivity[i][i_part] != NULL) {
            free(extrp->pextract_connectivity[i][i_part]);
          }
          if(extrp->pextract_connectivity_idx[i] != NULL) {
            if(extrp->pextract_connectivity_idx[i][i_part] != NULL) {
              free(extrp->pextract_connectivity_idx[i][i_part]);
            }
          }
        }
      }

      if(extrp->pextract_connectivity[i] != NULL) {
        free(extrp->pextract_connectivity[i]);
        extrp->pextract_connectivity[i] = NULL;
      }

      if(extrp->pextract_connectivity_idx[i] != NULL) {
        free(extrp->pextract_connectivity_idx[i]);
        extrp->pextract_connectivity_idx[i] = NULL;
      }
    }
  }

  /* Free ln_to_gn */
  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    if(extrp->is_owner_ln_to_gn[i] == PDM_TRUE) {
      if(extrp->pextract_entity_ln_to_gn[i] != NULL) {
        for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
          if(extrp->pextract_entity_ln_to_gn[i][i_part] != NULL) {
            free(extrp->pextract_entity_ln_to_gn[i][i_part]);
          }
        }

        free(extrp->pextract_entity_ln_to_gn[i]);
        extrp->pextract_entity_ln_to_gn[i] = NULL;
      }
    }
  }


  /* Free parent_ln_to_gn */
  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    if(extrp->is_owner_parent_ln_to_gn[i] == PDM_TRUE) {
      if(extrp->pextract_entity_parent_ln_to_gn[i] != NULL) {
        for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
          if(extrp->pextract_entity_parent_ln_to_gn[i][i_part] != NULL) {
            free(extrp->pextract_entity_parent_ln_to_gn[i][i_part]);
          }
        }

        free(extrp->pextract_entity_parent_ln_to_gn[i]);
        extrp->pextract_entity_parent_ln_to_gn[i] = NULL;
      }
    }
  }

  if(extrp->is_owner_vtx_coord == PDM_TRUE) {
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      free(extrp->pextract_vtx_coord[i_part]);
    }
    free(extrp->pextract_vtx_coord);
  }



  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    free(extrp->pextract_n_entity[i]);
  }


  free(extrp->is_owner_connectivity   );
  free(extrp->is_owner_ln_to_gn       );
  free(extrp->is_owner_parent_ln_to_gn);


  free(extrp);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
