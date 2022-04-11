
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_para_graph_dual.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_unique.h"
#include "pdm_array.h"
#include "pdm_logging.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
void
PDM_compress_connectivity
(
 int  n_entity1,
 int *entity1_entity2_idx,
 int *entity1_entity2
)
{

  int idx_comp  = 0; /* Compressed index use to fill the buffer */
  int idx_block = 0; /* Index in the block to post-treat        */
  entity1_entity2_idx[0] = 0;
  int need_shift = 0;

  int *entity1_entity2_n = (int * ) malloc( n_entity1 * sizeof(int));
  for(int i = 0; i < n_entity1; ++i) {
    entity1_entity2_n[i] = entity1_entity2_idx[i+1] - entity1_entity2_idx[i];
  }

  for(int i = 0; i < n_entity1; ++i){
    int n_cell_connect = entity1_entity2_n[i];

    /* Reshift next value in compressed block to avoid create a new shift */
    if(need_shift) {
      int idx_new = idx_comp;
      for(int j = idx_block; j < idx_block+n_cell_connect; ++j){
        entity1_entity2[idx_new++] = entity1_entity2[j];
      }
    }

    int end_connect         = idx_comp + n_cell_connect - 1;
    // printf(" idx_comp:: %d | end_connect:: %d \n", idx_comp, end_connect);

    int n_cell_connect_comp = PDM_inplace_unique(entity1_entity2, idx_comp, end_connect);
    // printf(" n_cell_connect:: %d | n_cell_connect_comp:: %d \n", n_cell_connect, n_cell_connect_comp);

    //Once a shift is needed, need_shift must stay at one
    if(n_cell_connect_comp < n_cell_connect) {
      need_shift = 1;
    }

    entity1_entity2_idx[i+1] = entity1_entity2_idx[i] + n_cell_connect_comp;
    idx_comp  += n_cell_connect_comp;
    idx_block += n_cell_connect;

    // printf("idx_comp : %d | idx_block : %d \n", idx_comp, idx_block);
  }

  free(entity1_entity2_n);

}



void
PDM_combine_connectivity
(
 int   n_entity1,
 int  *entity1_entity2_idx,
 int  *entity1_entity2,
 int  *entity2_entity3_idx,
 int  *entity2_entity3,
 int **entity1_entity3_idx,
 int **entity1_entity3
)
{
  int* _entity1_entity3_idx = (int * ) malloc( (n_entity1 + 1) * sizeof(int));
  _entity1_entity3_idx[0] = 0;
  for(int i_entity1 = 0; i_entity1 < n_entity1; ++i_entity1) {
    _entity1_entity3_idx[i_entity1+1] = _entity1_entity3_idx[i_entity1];
    for(int idx_1 = entity1_entity2_idx[i_entity1]; idx_1 < entity1_entity2_idx[i_entity1+1]; ++idx_1 ) {
      int i_entity2 = PDM_ABS(entity1_entity2[idx_1])-1;
      _entity1_entity3_idx[i_entity1+1] += entity2_entity3_idx[i_entity2+1] - entity2_entity3_idx[i_entity2];
    }
  }

  // printf("_entity1_entity3_idx[%i] = %i \n", n_entity1, _entity1_entity3_idx[n_entity1]);
  int* _entity1_entity3 = (int *) malloc( _entity1_entity3_idx[n_entity1] * sizeof(int));

  int idx = 0;
  for(int i_entity1 = 0; i_entity1 < n_entity1; ++i_entity1) {
    for(int idx_1 = entity1_entity2_idx[i_entity1]; idx_1 < entity1_entity2_idx[i_entity1+1]; ++idx_1 ) {
      int i_entity2 = PDM_ABS(entity1_entity2[idx_1])-1;
      for(int idx_2 = entity2_entity3_idx[i_entity2]; idx_2 < entity2_entity3_idx[i_entity2+1]; ++idx_2 ) {
        _entity1_entity3[idx++] = entity2_entity3[idx_2];
      }
    }
  }

  PDM_compress_connectivity(n_entity1, _entity1_entity3_idx, _entity1_entity3);

  // PDM_log_trace_array_int(_entity1_entity3_idx, n_entity1+1                     , "_entity1_entity3_idx::");
  // PDM_log_trace_array_int(_entity1_entity3    , _entity1_entity3_idx[n_entity1], "_entity1_entity3::");

  _entity1_entity3 = realloc(_entity1_entity3, _entity1_entity3_idx[n_entity1] * sizeof(int));

  *entity1_entity3_idx = _entity1_entity3_idx;
  *entity1_entity3     = _entity1_entity3;
}



void
PDM_connectivity_transpose
(
const int   n_entity1,
const int   n_entity2,
      int  *entity1_entity2_idx,
      int  *entity1_entity2,
      int **entity2_entity1_idx,
      int **entity2_entity1
)
{
  int* _entity2_entity1_idx = (int * ) malloc( (n_entity2 + 1) * sizeof(int));
  int* entity2_entity1_n    = PDM_array_zeros_int(n_entity2);


  for(int i_entity1 = 0; i_entity1 < n_entity1; ++i_entity1) {
    for(int idx_1 = entity1_entity2_idx[i_entity1]; idx_1 < entity1_entity2_idx[i_entity1+1]; ++idx_1 ) {
      int i_entity2 = PDM_ABS(entity1_entity2[idx_1])-1;
      entity2_entity1_n[i_entity2] += 1;
    }
  }

  PDM_array_idx_from_sizes_int(entity2_entity1_n, n_entity2, _entity2_entity1_idx);
  PDM_array_reset_int(entity2_entity1_n, n_entity2, 0);

  int* _entity2_entity1 = (int * ) malloc( _entity2_entity1_idx[n_entity2] * sizeof(int));

  for(int i_entity1 = 0; i_entity1 < n_entity1; ++i_entity1) {
    for(int idx_1 = entity1_entity2_idx[i_entity1]; idx_1 < entity1_entity2_idx[i_entity1+1]; ++idx_1 ) {
      int i_entity2 = PDM_ABS(entity1_entity2[idx_1])-1;
      int idx = _entity2_entity1_idx[i_entity2] + entity2_entity1_n[i_entity2]++;
      int sign = PDM_SIGN(entity1_entity2[idx_1]);
      _entity2_entity1[idx] = sign * (i_entity1+1);
    }
  }

  // PDM_log_trace_array_int(_entity2_entity1_idx, n_entity2+1, "_entity2_entity1_idx::");
  // PDM_log_trace_array_int(entity2_entity1_n   , n_entity2  , "entity2_entity1_n::");
  // PDM_log_trace_array_int(_entity2_entity1   , _entity2_entity1_idx[n_entity2]  , "_entity2_entity1::");

  free(entity2_entity1_n);

  *entity2_entity1_idx = _entity2_entity1_idx;
  *entity2_entity1     = _entity2_entity1;


}

void
PDM_part_connectivity_transpose
(
const int    n_part,
const int   *n_entity1,
const int   *n_entity2,
      int  **entity1_entity2_idx,
      int  **entity1_entity2,
      int ***entity2_entity1_idx,
      int ***entity2_entity1
)
{
  *entity2_entity1_idx = (int ** ) malloc( n_part * sizeof(int *));
  *entity2_entity1     = (int ** ) malloc( n_part * sizeof(int *));

  int** _entity2_entity1_idx = *entity2_entity1_idx;
  int** _entity2_entity1     = *entity2_entity1;

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int* _pentity2_entity1_idx;
    int* _pentity2_entity1;

    PDM_connectivity_transpose(n_entity1[i_part],
                               n_entity2[i_part],
                               entity1_entity2_idx[i_part],
                               entity1_entity2[i_part],
                              &_pentity2_entity1_idx,
                              &_pentity2_entity1);

    _entity2_entity1_idx[i_part] = _pentity2_entity1_idx;
    _entity2_entity1[i_part]     = _pentity2_entity1;
  }

}



void
PDM_part_connectivity_to_connectity_idx
(
const int    n_part,
const int   *n_entity1,
      int  **entity1_entity2_in,
      int ***entity1_entity2_idx,
      int ***entity1_entity2
)
{

  *entity1_entity2_idx = (int ** ) malloc( n_part * sizeof(int *));
  *entity1_entity2     = (int ** ) malloc( n_part * sizeof(int *));

  int** _entity1_entity2_idx = *entity1_entity2_idx;
  int** _entity1_entity2     = *entity1_entity2;

  for(int i_part = 0; i_part < n_part; ++i_part) {

    _entity1_entity2_idx[i_part] = (int * ) malloc( (    n_entity1[i_part] + 1) * sizeof(int));
    _entity1_entity2    [i_part] = (int * ) malloc( (2 * n_entity1[i_part]    ) * sizeof(int));

    int idx = 0;
    _entity1_entity2_idx[i_part][0] = 0;
    for(int i_entity = 0; i_entity < n_entity1[i_part]; ++i_entity) {
      _entity1_entity2_idx[i_part][i_entity+1] = _entity1_entity2_idx[i_part][i_entity];
      if(entity1_entity2_in[i_part][2*i_entity + 1 ] == 0) {
        _entity1_entity2_idx[i_part][i_entity+1]++;
        _entity1_entity2[i_part][idx++] = entity1_entity2_in[i_part][2*i_entity];
      } else if(entity1_entity2_in[i_part][2*i_entity] == 0) {
        _entity1_entity2_idx[i_part][i_entity+1]++;
        _entity1_entity2[i_part][idx++] = - entity1_entity2_in[i_part][2*i_entity+1];
      } else {
        _entity1_entity2_idx[i_part][i_entity+1] += 2;
        _entity1_entity2[i_part][idx++] =   entity1_entity2_in[i_part][2*i_entity  ];
        _entity1_entity2[i_part][idx++] = - entity1_entity2_in[i_part][2*i_entity+1];
      }
    }

    _entity1_entity2    [i_part] = realloc(_entity1_entity2    [i_part], idx * sizeof(int));


  }
  // int* face_cell_idx = (int *) malloc( (pn_faces[0] + 1 ) * sizeof(int));
  // int* face_cell     = (int *) malloc( (2 * pn_faces[0] ) * sizeof(int));
  // int idx = 0;
  // face_cell_idx[0] = 0;
  // for(int i_face = 0; i_face < pn_faces[0]; ++i_face) {
  //   face_cell_idx[i_face+1] = face_cell_idx[i_face];
  //   if(pface_cell[0][2*i_face + 1 ] == 0) {
  //     face_cell_idx[i_face+1]++;
  //     face_cell[idx++] = pface_cell[0][2*i_face];
  //   } else {
  //     face_cell_idx[i_face+1] += 2;
  //     face_cell[idx++] = pface_cell[0][2*i_face  ];
  //     face_cell[idx++] = pface_cell[0][2*i_face+1];
  //   }
  // }

}


void
PDM_compute_face_vtx_from_face_and_edge
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
)
{
  int dbg = 0;

  *face_vtx = malloc (sizeof(int) * face_edge_idx[n_face]);

  int n_edge = 0;
  for (int i = 0; i < face_edge_idx[n_face]; i++) {
    n_edge = PDM_MAX(n_edge, PDM_ABS(face_edge[i]));
  }


  int *edge_tag = PDM_array_zeros_int(n_edge);

  for (int iface = 0; iface < n_face; iface++) {
    int *_face_vtx  = *face_vtx  + face_edge_idx[iface];
    int *_face_edge =  face_edge + face_edge_idx[iface];

    if (dbg) {
      log_trace("\nFace %d\n", iface);
      for (int idx_edge = face_edge_idx[iface]; idx_edge < face_edge_idx[iface+1]; idx_edge++) {
        int iedge = PDM_ABS(face_edge[idx_edge]) - 1;
        log_trace("  edge %d: %d %d\n",
                  face_edge[idx_edge],
                  edge_vtx[2*iedge], edge_vtx[2*iedge+1]);
      }
    }

    int _n_edge = face_edge_idx[iface+1] - face_edge_idx[iface];
    // first edge
    int iedge = PDM_ABS(_face_edge[0]) - 1;
    edge_tag[iedge] = 1;
    _face_vtx[0] = edge_vtx[2*iedge  ];
    _face_vtx[1] = edge_vtx[2*iedge+1];

    for (int i = 2; i < _n_edge; i++) {

      for (int j = 1; j < _n_edge; j++) {
        iedge = PDM_ABS(_face_edge[j]) - 1;

        if (edge_tag[iedge]) {
          continue;
        }

        if (edge_vtx[2*iedge] == _face_vtx[i-1]) {
          _face_vtx[i] = edge_vtx[2*iedge+1];
          edge_tag[iedge] = 1;
          break;
        }
        else if (edge_vtx[2*iedge+1] == _face_vtx[i-1]) {
          _face_vtx[i] = edge_vtx[2*iedge];
          edge_tag[iedge] = 1;
          break;
        }
      }
    }

    if (dbg) {
      log_trace("  face_vtx = ");
      for (int ivtx = 0; ivtx < face_edge_idx[iface+1] - face_edge_idx[iface]; ivtx++) {
        log_trace("%d ", _face_vtx[ivtx]);
      }
      log_trace("\n");
    }

    // reset tags
    for (int i = 0; i < _n_edge; i++) {
      iedge = PDM_ABS(_face_edge[i]) - 1;
      edge_tag[iedge] = 0;
    }
  }
  free(edge_tag);
}
