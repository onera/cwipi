
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
#include "pdm_dconnectivity_transform.h"
#include "pdm_para_graph_dual.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_unique.h"
#include "pdm_logging.h"
#include "pdm_array.h"

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
_deduce_combine_connectivity_impl
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
 const int             *dentity2_entity3_idx,
 const PDM_g_num_t     *dentity2_entity3,
 const int              is_signed,
       int            **dentity1_entity3_n,
        PDM_g_num_t   **dentity1_entity3
)
{
  int i_rank, n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity1 = entity1_distrib[i_rank+1] - entity1_distrib[i_rank];
  int dn_entity2 = entity2_distrib[i_rank+1] - entity2_distrib[i_rank];

  /*
   * Some print to create unit_test
   */
  // PDM_log_trace_array_long(entity1_distrib  , n_rank+1               , "entity1_distrib::");
  // PDM_log_trace_array_long(entity2_distrib  , n_rank+1               , "entity2_distrib::");
  // PDM_log_trace_array_int (dentity1_entity2_idx, dn_entity1+1              , "dentity1_entity2_idx::");
  // PDM_log_trace_array_long(dentity1_entity2    , dentity1_entity2_idx[dn_entity1], "dentity1_entity2::");
  // PDM_log_trace_array_int (dentity2_entity3_idx , dn_entity2+1              , "dentity2_entity3_idx::");
  // PDM_log_trace_array_long(dentity2_entity3     , dentity2_entity3_idx[dn_entity2] , "dentity2_entity3::");

  int* dentity2_entity3_n = (int * ) malloc( dn_entity2 * sizeof(int));
  for(int i = 0; i < dn_entity2; ++i) {
    dentity2_entity3_n[i] = dentity2_entity3_idx[i+1] - dentity2_entity3_idx[i];
  }

  PDM_g_num_t* dentity1_entity2_cur = NULL;

  if(is_signed) {
    dentity1_entity2_cur = (PDM_g_num_t* ) malloc( dentity1_entity2_idx[dn_entity1] * sizeof(PDM_g_num_t));
    for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
      dentity1_entity2_cur[i] = PDM_ABS(dentity1_entity2[i]);
    }
  } else {
    dentity1_entity2_cur = (PDM_g_num_t *) dentity1_entity2;
  }

  // PDM_log_trace_array_int(dentity1_entity2_idx, dn_entity1+1, "dentity1_entity2_idx::");
  // PDM_log_trace_array_long(dentity1_entity2_cur, dentity1_entity2_idx[dn_entity1], "dentity1_entity2::");

  /*
   * First compute the dentity1_entity3 connectivity
   *      -  We use ln_to_gn = dentity1_entity2 in order to have for each partition the entity2s
   *      -  We exhange the entity2_vtx, entity2_vtx_idx
   * So, for each entity2 describe by dentity1_entity2, we receive all vtx
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity2_distrib,
                               (const PDM_g_num_t **) &dentity1_entity2_cur,
                                                      &dentity1_entity2_idx[dn_entity1],
                                                      1,
                                                      comm);

  /*
   * Free
   */
  if(is_signed) {
    free(dentity1_entity2_cur);
  }

  /*
   * Exchange
   */
  int**         pentity1_entity3_n;
  PDM_g_num_t** pentity1_entity3;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          dentity2_entity3_n,
             (void *  )   dentity2_entity3,
             (int  ***)  &pentity1_entity3_n,
             (void ***)  &pentity1_entity3);
  free(dentity2_entity3_n);

  if(is_signed) {
    int idx_read = 0;
    for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
      int sign = PDM_SIGN(dentity1_entity2[i]);
      for(int j = 0; j < pentity1_entity3_n[0][i]; ++j){
        pentity1_entity3[0][idx_read] = sign * pentity1_entity3[0][idx_read];
        idx_read++;
      }
    }
  }


  /*
   * Panic Verbose
   */
  // int* pentity1_entity3_idx = (int*) malloc( (dentity1_entity2_idx[dn_entity1] + 1) * sizeof(int));
  // pentity1_entity3_idx[0] = 0;
  // for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
  //   pentity1_entity3_idx[i+1] = pentity1_entity3_idx[i] + pentity1_entity3_n[0][i];
  // }

  // PDM_log_trace_array_int(pentity1_entity3_n[0], dentity1_entity2_idx[dn_entity1], "pentity1_entity3_n::");
  // PDM_log_trace_array_int(pentity1_entity3_idx, dentity1_entity2_idx[dn_entity1]+1, "pentity1_entity3_idx::");
  // PDM_log_trace_array_long(pentity1_entity3[0], pentity1_entity3_idx[dentity1_entity2_idx[dn_entity1]], "pentity1_entity3::");

  /*
   * Free
   */
  PDM_block_to_part_free(btp);

  /*
   * Assign pointer
   */
  *dentity1_entity3   = pentity1_entity3[0];
  *dentity1_entity3_n = pentity1_entity3_n[0];

  /*
   * Free first level of pointer - the second level is hold by dentity1_entity3/dentity1_entity3_n
   */
  free(pentity1_entity3  );
  free(pentity1_entity3_n);

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Compute the combine connectivty of entity1 with entity2 to entity3
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   entity1_distrib       distribution of entity1 over the procs (size=n_rank+1)
 * \param [in]   entity2_distrib       distribution of entity2 over the procs (size=n_rank+1)
 * \param [in]   dentity1_entity2_idx
 * \param [in]   dentity1_entity2
 * \param [in]   dentity1_entity2      is array is signed
 * \param [in]   dentity2_entity1_idx
 * \param [in]   dentity2_entity1
 * \param [in]   dentity1_entity3_idx
 * \param [in]   dentity1_entity3
 */
void
PDM_deduce_combine_connectivity
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
 const int             *dentity2_entity3_idx,
 const PDM_g_num_t     *dentity2_entity3,
 const int              is_signed,
       int            **dentity1_entity3_idx,
       PDM_g_num_t    **dentity1_entity3
)
{
  int* pentity1_entity3_n;
  _deduce_combine_connectivity_impl(comm,
                                    entity1_distrib,
                                    entity2_distrib,
                                    dentity1_entity2_idx,
                                    dentity1_entity2,
                                    dentity2_entity3_idx,
                                    dentity2_entity3,
                                    is_signed,
                                    &pentity1_entity3_n,
                                    dentity1_entity3);

  PDM_g_num_t* _dentity1_entity3 = *dentity1_entity3;

  int i_rank, n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity1 = entity1_distrib[i_rank+1] - entity1_distrib[i_rank];

  /*
   * Post-treat
   */
  *dentity1_entity3_idx = (int*) malloc( (dn_entity1 + 1) * sizeof(int));
  int* _dentity1_entity3_idx = *dentity1_entity3_idx;
  int*  dentity1_entity3_n   = (int*) malloc( (dn_entity1    ) * sizeof(int));

  int idx = 0;
  _dentity1_entity3_idx[0] = 0;
  for(int i_entity1 = 0; i_entity1 < dn_entity1; ++i_entity1) {

    dentity1_entity3_n   [i_entity1  ] = 0;
    _dentity1_entity3_idx[i_entity1+1] = _dentity1_entity3_idx[i_entity1];

    int n_entity2_per_entity1 = dentity1_entity2_idx[i_entity1+1] - dentity1_entity2_idx[i_entity1];
    // printf("n_entity2_per_entity1::%i\n", n_entity2_per_entity1);

    for(int i_entity2 = 0; i_entity2 < n_entity2_per_entity1; ++i_entity2) {
      _dentity1_entity3_idx[i_entity1+1] += pentity1_entity3_n[idx];
      // if(is_signed) {
      //   int beg = dentity1_entity2_idx[i_entity1];
      //   int sgn = PDM_SIGN(dentity1_entity2[beg+i_entity2]);
      //   _dentity1_entity3[idx] = _dentity1_entity3[idx] * sgn;
      //   dentity1_entity3_n[i_entity1]      += pentity1_entity3_n[idx++];
      // } else {
        dentity1_entity3_n[i_entity1]      += pentity1_entity3_n[idx++];
      // }
    }
  }
  // printf("idx::%i\n", idx);
  // printf("dentity1_entity2_idx[dn_entity1]::%i\n", dentity1_entity2_idx[dn_entity1]);

  assert(idx == dentity1_entity2_idx[dn_entity1]);
  free(pentity1_entity3_n);

  // PDM_log_trace_array_int(_dentity1_entity3_idx, dn_entity1+1, "_dentity1_entity3_idx::");
  // PDM_log_trace_array_int(dentity1_entity3_n  , dn_entity1  , "dentity1_entity3_n::");

  PDM_para_graph_compress_connectivity(dn_entity1, _dentity1_entity3_idx, dentity1_entity3_n, _dentity1_entity3);

  /*
   * Realloc
   */
  *dentity1_entity3 = (PDM_g_num_t *) realloc(*dentity1_entity3, sizeof(PDM_g_num_t) * _dentity1_entity3_idx[dn_entity1] );
  _dentity1_entity3 = *dentity1_entity3;

  // PDM_log_trace_array_int (_dentity1_entity3_idx, dn_entity1+1              , "after -> _dentity1_entity3_idx::");
  // PDM_log_trace_array_int (dentity1_entity3_n   , dn_entity1                , "after -> dentity1_entity3_n::");
  // PDM_log_trace_array_long(_dentity1_entity3    , _dentity1_entity3_idx[dn_entity1], "after -> dentity1_entity3::");

  /*
   * Free
   */
  free(dentity1_entity3_n);
}


/**
 *
 * \brief Compute the combine connectivty of entity1 with entity2 to entity3
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   entity1_distrib       distribution of entity1 over the procs (size=n_rank+1)
 * \param [in]   entity2_distrib       distribution of entity2 over the procs (size=n_rank+1)
 * \param [in]   dentity1_entity2_idx
 * \param [in]   dentity1_entity2
 * \param [in]   dentity1_entity2      is array is signed
 * \param [in]   dentity2_entity1_idx
 * \param [in]   dentity2_entity1
 * \param [out]   dentity1_entity3_idx
 * \param [out]   dentity1_entity3
 */
void
PDM_deduce_combine_connectivity_dual
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
 const int             *dentity2_entity3_idx,
 const PDM_g_num_t     *dentity2_entity3,
 const int              is_signed,
       PDM_g_num_t    **dentity1_entity3_idx,
       PDM_g_num_t    **dentity1_entity3
)
{
  int* pentity1_entity3_n;
  _deduce_combine_connectivity_impl(comm,
                                    entity1_distrib,
                                    entity2_distrib,
                                    dentity1_entity2_idx,
                                    dentity1_entity2,
                                    dentity2_entity3_idx,
                                    dentity2_entity3,
                                    is_signed,
                                    &pentity1_entity3_n,
                                    dentity1_entity3);

  PDM_g_num_t* _dentity1_entity3 = *dentity1_entity3;

  int i_rank, n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity1 = entity1_distrib[i_rank+1] - entity1_distrib[i_rank];

  /*
   * Post-treat
   */
  *dentity1_entity3_idx = (PDM_g_num_t*) malloc( (dn_entity1 + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t* _dentity1_entity3_idx = *dentity1_entity3_idx;
  int*  dentity1_entity3_n   = (int*) malloc( (dn_entity1    ) * sizeof(int));

  int idx = 0;
  _dentity1_entity3_idx[0] = 0;
  for(int i_entity1 = 0; i_entity1 < dn_entity1; ++i_entity1) {

    dentity1_entity3_n   [i_entity1  ] = 0;
    _dentity1_entity3_idx[i_entity1+1] = _dentity1_entity3_idx[i_entity1];

    int n_entity2_per_entity1 = dentity1_entity2_idx[i_entity1+1] - dentity1_entity2_idx[i_entity1];
    // printf("n_entity2_per_entity1::%i\n", n_entity2_per_entity1);

    for(int i_entity2 = 0; i_entity2 < n_entity2_per_entity1; ++i_entity2) {
      _dentity1_entity3_idx[i_entity1+1] += (PDM_g_num_t) pentity1_entity3_n[idx];
      dentity1_entity3_n[i_entity1]      += (PDM_g_num_t) pentity1_entity3_n[idx++];
    }
  }
  // printf("idx::%i\n", idx);
  // printf("dentity1_entity2_idx[dn_entity1]::%i\n", dentity1_entity2_idx[dn_entity1]);

  for(int i = 0; i < _dentity1_entity3_idx[dn_entity1]; ++i) {
    _dentity1_entity3[i] = PDM_ABS(_dentity1_entity3[i]);
  }

  assert(idx == dentity1_entity2_idx[dn_entity1]);
  free(pentity1_entity3_n);

  // PDM_log_trace_array_long(_dentity1_entity3_idx, dn_entity1+1, "_dentity1_entity3_idx::");
  // PDM_log_trace_array_int (dentity1_entity3_n   , dn_entity1  , "dentity1_entity3_n::");

  PDM_para_graph_compress_connectivity_dual(dn_entity1,
                                            entity1_distrib[i_rank],
                                            _dentity1_entity3_idx,
                                            dentity1_entity3_n,
                                            _dentity1_entity3);

  /*
   * Realloc
   */
  *dentity1_entity3 = (PDM_g_num_t *) realloc(*dentity1_entity3, sizeof(PDM_g_num_t) * _dentity1_entity3_idx[dn_entity1] );

  // PDM_log_trace_array_long(_dentity1_entity3_idx, dn_entity1+1              , "after -> _dentity1_entity3_idx::");
  // PDM_log_trace_array_int(dentity1_entity3_n   , dn_entity1                , "after -> dentity1_entity3_n::");
  // PDM_log_trace_array_long(dentity1_entity3[0]    , _dentity1_entity3_idx[dn_entity1], "after -> dentity1_entity3::");

  /*
   * Free
   */
  free(dentity1_entity3_n);
}

/**
 *
 * \brief Compute the dual connectivty of entity1
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   entity1_distrib       distribution of entity1 over the procs (size=n_rank+1)
 * \param [in]   entity2_distrib       distribution of entity2 over the procs (size=n_rank+1)
 * \param [in]   dentity1_entity2_idx
 * \param [in]   dentity1_entity2
 * \param [in]   dentity1_entity2      is array is signed
 * \param [in]   dentity2_entity1_idx
 * \param [in]   dentity2_entity1
 */
void
PDM_dconnectivity_transpose
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
       PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
       int              is_signed,
       int            **dentity2_entity1_idx,
       PDM_g_num_t    **dentity2_entity1
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity1  = entity1_distrib[i_rank+1] - entity1_distrib[i_rank];
  // int dn_entity2  = entity2_distrib[i_rank+1] - entity2_distrib[i_rank];

  PDM_g_num_t* ln_to_gn = NULL;
  PDM_g_num_t* gnum = (PDM_g_num_t * ) malloc( dentity1_entity2_idx[dn_entity1] * sizeof(PDM_g_num_t));

  PDM_g_num_t shift_g = 1 + entity1_distrib[i_rank]; // Entre 1 et N

  if(is_signed) {
    ln_to_gn = (PDM_g_num_t * ) malloc( dentity1_entity2_idx[dn_entity1] * sizeof(PDM_g_num_t));

    for(int i_entity = 0; i_entity < dn_entity1; ++i_entity) {
      for(int j = dentity1_entity2_idx[i_entity]; j < dentity1_entity2_idx[i_entity+1]; ++j) {
        int g_sign = PDM_SIGN(dentity1_entity2[j]);
        gnum[j] = g_sign * (i_entity + shift_g);
        ln_to_gn[j] = PDM_ABS(dentity1_entity2[j]);
      }
    }
  } else {
    ln_to_gn = (PDM_g_num_t * ) dentity1_entity2;
    for(int i_entity = 0; i_entity < dn_entity1; ++i_entity) {
      for(int j = dentity1_entity2_idx[i_entity]; j < dentity1_entity2_idx[i_entity+1]; ++j) {
        gnum[j] = i_entity + shift_g;
      }
    }
  }

  // PDM_log_trace_array_long(ln_to_gn, dentity1_entity2_idx[dn_entity1], "ln_to_gn::");
  // PDM_log_trace_array_long(gnum,  dentity1_entity2_idx[dn_entity1], "gnum::");
  /*
   * In order to revert the conncectivty we use the global numbering property
   */
  PDM_part_to_block_t *ptb = NULL;

  int save_entity_distrib = 0;
  if(entity2_distrib != NULL && entity2_distrib[0] != -1) {
    ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                    PDM_PART_TO_BLOCK_POST_MERGE,
                                    1.,
                                    &ln_to_gn,
                                    entity2_distrib,
                         (int *)    &dentity1_entity2_idx[dn_entity1],
                                     1,
                                     comm);
  } else {
    save_entity_distrib = 1;
    ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                   PDM_PART_TO_BLOCK_POST_MERGE,
                                   1.,
                                   &ln_to_gn,
                                   NULL,
                         (int *)   &dentity1_entity2_idx[dn_entity1],
                                    1,
                                    comm);
  }

  int* send_stri = PDM_array_const_int(dentity1_entity2_idx[dn_entity1], 1);

  if(is_signed) {
    free(ln_to_gn);
  }

  int         *dentity2_entity1_n = NULL;
  PDM_g_num_t *recv_data          = NULL;

  int blk_size = PDM_part_to_block_exch (ptb,
                                         sizeof(PDM_g_num_t),
                                         PDM_STRIDE_VAR,
                                         1,
                                         &send_stri,
                               (void **) &gnum,
                                         &dentity2_entity1_n,
                               (void **) &recv_data);
  PDM_UNUSED(blk_size);

  int dn_entity2_recv = PDM_part_to_block_n_elt_block_get(ptb);

  if(entity2_distrib != NULL && save_entity_distrib != 1) {
    PDM_g_num_t *distrib2_idx_full =
      PDM_part_to_block_adapt_partial_block_to_block (ptb,
                                                      &dentity2_entity1_n,
                                                      entity2_distrib[n_rank]);
    dn_entity2_recv = distrib2_idx_full[i_rank+1] - distrib2_idx_full[i_rank];
    free(distrib2_idx_full);
  }

  if(save_entity_distrib == 1) {
    // Update distrib
    PDM_g_num_t* ptb_distrib = PDM_part_to_block_distrib_index_get(ptb);
    for(int i = 0; i < n_rank+1; ++i) {
      entity2_distrib[i] = ptb_distrib[i];
    }
  }

  /*
   * Free
   */
  PDM_part_to_block_free(ptb);
  free(gnum);
  free(send_stri);

  /*
   * Allocate
   */
  *dentity2_entity1_idx = (int *) malloc( (dn_entity2_recv + 1) * sizeof(int));
  int* _dentity2_entity1_idx = *dentity2_entity1_idx;

  // printf("blk_size       ::%i\n", blk_size       );
  // printf("dn_entity2_recv::%i\n", dn_entity2_recv);
  // log_trace("blk_size = %i | dn_entity2_recv = %i \n", blk_size, dn_entity2_recv);

  // PDM_log_trace_array_long(dentity2_entity1_n, dn_entity2_recv, "Before : dentity2_entity1_n::");
  // PDM_log_trace_array_long(recv_data, blk_size, "Before : recv_data::");

  PDM_para_graph_compress_connectivity(dn_entity2_recv,
                                        _dentity2_entity1_idx,
                                        dentity2_entity1_n,
                                        recv_data);
  // printf("*dentity2_entity1_idx[dn_entity2_recv]       ::%i\n", _dentity2_entity1_idx[dn_entity2_recv]       );

  // *dentity2_entity1 = recv_data;

  /*
   * Realloc
   */
  *dentity2_entity1 = realloc(recv_data, _dentity2_entity1_idx[dn_entity2_recv] * sizeof(PDM_g_num_t));
  // PDM_g_num_t* _dentity2_entity1 = *dentity2_entity1;

  // PDM_log_trace_array_int (_dentity2_entity1_idx, dn_entity2_recv+1         , "_dentity2_entity1_idx::");
  // PDM_log_trace_array_long(*dentity2_entity1, _dentity2_entity1_idx[dn_entity2_recv], "recv_data::");

  free(dentity2_entity1_n);
}


/**
 *
 * \brief Compute the dual connectivty of entity1
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   entity1_distrib       distribution of entity1 over the procs (size=n_rank+1)
 * \param [in]   entity2_distrib       distribution of entity2 over the procs (size=n_rank+1)
 * \param [in]   dentity1_entity2_idx
 * \param [in]   dentity1_entity2
 * \param [in]   dentity1_entity2      is array is signed
 * \param [in]   dentity2_entity1_idx
 * \param [in]   dentity2_entity1
 */
void
PDM_dorder_reverse
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity_distrib,
 const PDM_g_num_t     *dentity1_entity2,
       PDM_g_num_t    **dentity2_entity1
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int dn_entity1  = entity_distrib[i_rank+1] - entity_distrib[i_rank];
  PDM_g_num_t* gnum = (PDM_g_num_t * ) malloc( dn_entity1 * sizeof(PDM_g_num_t));

  for(int i_entity = 0; i_entity < dn_entity1; ++i_entity) {
    gnum[i_entity] = entity_distrib[i_rank] + i_entity + 1;
  }

  /*
   * In order to revert the conncectivty we use the global numbering property
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                       1.,
                                      (PDM_g_num_t **) &dentity1_entity2,
                                                       entity_distrib,
                                            (int *)    &dn_entity1,
                                                       1,
                                                       comm);

  PDM_g_num_t *recv_data          = NULL;

  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          1,
                          NULL,
                (void **) &gnum,
                           NULL,
                 (void **) &recv_data);
  free(gnum);

  *dentity2_entity1 = recv_data;
  PDM_part_to_block_free(ptb);
}

// void
// PDM_dconnectivity_update_child_connectivity
// (
//  const PDM_MPI_Comm     comm,
//  const PDM_g_num_t     *entity1_distrib,
//        PDM_g_num_t     *entity2_distrib,
//        int             *dentity1_entity2_idx,
//        PDM_g_num_t     *dentity1_entity2,
//        int              is_signed,
//  const PDM_g_num_t     *dold_to_new_entity2
// )
// {

//   // On peut également updater sans echanger les doublons je pense :
//   //  PDM_unique (avec unique order ) Puis on replace aprés l'échange !!!
//   //  Maybe il faut order + unique order ?
// }

// void
// PDM_dconnectivity_reorder
// (
//  const PDM_MPI_Comm     comm,
//  const PDM_g_num_t     *entity1_distrib,
//        int             *dentity1_entity2_idx,
//        PDM_g_num_t     *dentity1_entity2,
//        int              is_signed,
//  const PDM_g_num_t     *dold_to_new_entity1
// )
// {

// }


void
PDM_dgroup_entity_transpose
(
 int            n_group,
 int           *dgroup_entity_idx,
 PDM_g_num_t   *dgroup_entity,
 PDM_g_num_t   *distrib_entity,
 int          **dentity_group_idx,
 int          **dentity_group,
 PDM_MPI_Comm   comm
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  //printf(" dgroup_entity_idx[%i] = %i \n", n_group, dgroup_entity_idx[n_group]);
  PDM_part_to_block_t *ptb = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                        1.,
                                                        &dgroup_entity,
                                                        distrib_entity,
                                                        &dgroup_entity_idx[n_group],
                                                        1,
                                                        comm);


  int* pgroup_id_n = (int *) malloc(dgroup_entity_idx[n_group] * sizeof(int));
  int* pgroup_id   = (int *) malloc(dgroup_entity_idx[n_group] * sizeof(int));

  for(int i_group = 0; i_group < n_group; ++i_group) {
    for(int i = dgroup_entity_idx[i_group]; i < dgroup_entity_idx[i_group+1]; ++i){
      pgroup_id_n[i] = 1;
      pgroup_id  [i] = i_group;
    }
  }

  /*
   *  Exchange group id
   */

  int *tmp_dentity_group_n = NULL;
  int *tmp_dentity_group   = NULL;
  int s_block = PDM_part_to_block_exch (ptb,
                                        sizeof(int),
                                        PDM_STRIDE_VAR,
                                        1,
                                        &pgroup_id_n,
                              (void **) &pgroup_id,
                                        &tmp_dentity_group_n,
                              (void **) &tmp_dentity_group);
  free(pgroup_id_n);
  free(pgroup_id);


  int dn_entity = distrib_entity[i_rank+1] - distrib_entity[i_rank];

  PDM_g_num_t* tmp_distrib = PDM_part_to_block_adapt_partial_block_to_block(ptb, &tmp_dentity_group_n, distrib_entity[n_rank]);
  free(tmp_distrib);
  PDM_part_to_block_free (ptb);

  /*
   * Post-treatment
   */
  *dentity_group     = malloc(s_block       * sizeof(int));
  *dentity_group_idx = malloc((dn_entity+1) * sizeof(int));
  int *_dentity_group     = *dentity_group;
  int *_dentity_group_idx = *dentity_group_idx;

  int idx_read  = 0;
  int idx_write = 0;
  _dentity_group_idx[0] = 0;
  for(int i = 0; i < dn_entity; ++i) {
    int n_id = tmp_dentity_group_n[i];
    int n_unique_id = 0;
    if(n_id > 0) {
      n_unique_id = PDM_inplace_unique(&tmp_dentity_group[idx_read], 0, n_id-1);
    }

    for(int j = 0; j < n_unique_id; ++j) {
      _dentity_group[idx_write++] = tmp_dentity_group[idx_read+j];
    }
    _dentity_group_idx[i+1] = _dentity_group_idx[i] + n_unique_id;

    idx_read += n_id;
  }

  free(tmp_dentity_group_n);
  free(tmp_dentity_group);

  if(0 == 1) {
    PDM_log_trace_connectivity_int(_dentity_group_idx, _dentity_group, dn_entity, "_dentity_group ::");
  }

  *dentity_group = realloc(*dentity_group, _dentity_group_idx[dn_entity] * sizeof(int));
}


#ifdef  __cplusplus
}
#endif
