
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

  PDM_g_num_t* entity_distribution_ptb = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * (n_rank+1) );
  for(int i = 0; i < n_rank+1; ++i){
    entity_distribution_ptb[i] = entity2_distrib[i] - 1;
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
  // PDM_log_trace_array_int(dentity1_entity2_cur, dentity1_entity2_idx[dn_entity1], "dentity1_entity2::");

  /*
   * First compute the dentity1_entity3 connectivity
   *      -  We use ln_to_gn = dentity1_entity2 in order to have for each partition the entity2s
   *      -  We exhange the entity2_vtx, entity2_vtx_idx
   * So, for each entity2 describe by dentity1_entity2, we receive all vtx
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distribution_ptb,
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
  int**         pentity2_vtx_n;
  PDM_g_num_t** pentity2_vtx;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          dentity2_entity3_n,
             (void *  )   dentity2_entity3,
             (int  ***)  &pentity2_vtx_n,
             (void ***)  &pentity2_vtx);
  free(dentity2_entity3_n);

  /*
   * Panic Verbose
   */
  // int* pentity2_vtx_idx = (int*) malloc( (dentity1_entity2_idx[dn_entity1] + 1) * sizeof(int));
  // pentity2_vtx_idx[0] = 0;
  // for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
  //   pentity2_vtx_idx[i+1] = pentity2_vtx_idx[i] + pentity2_vtx_n[0][i];
  // }

  // PDM_log_trace_array_int(pentity2_vtx_n[0], dentity1_entity2_idx[dn_entity1], "pentity2_vtx_n::");
  // PDM_log_trace_array_int(pentity2_vtx_idx, dentity1_entity2_idx[dn_entity1]+1, "pentity2_vtx_idx::");
  // PDM_log_trace_array_long(pentity2_vtx[0], pentity2_vtx_idx[dentity1_entity2_idx[dn_entity1]], "pentity2_vtx::");

  /*
   * Free
   */
  free(entity_distribution_ptb);
  PDM_block_to_part_free(btp);

  /*
   * Assign pointer
   */
  *dentity1_entity3   = pentity2_vtx[0];
  *dentity1_entity3_n = pentity2_vtx_n[0];

  /*
   * Free first level of pointer - the second level is hold by dentity1_entity3/dentity1_entity3_n
   */
  free(pentity2_vtx  );
  free(pentity2_vtx_n);

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
      dentity1_entity3_n[i_entity1]      += pentity1_entity3_n[idx++];
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
 * \param [in]   dentity1_entity3_idx
 * \param [in]   dentity1_entity3
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

  assert(idx == dentity1_entity2_idx[dn_entity1]);
  free(pentity1_entity3_n);

  // PDM_log_trace_array_int(_dentity1_entity3_idx, dn_entity1+1, "_dentity1_entity3_idx::");
  // PDM_log_trace_array_int(dentity1_entity3_n  , dn_entity1  , "dentity1_entity3_n::");

  PDM_para_graph_compress_connectivity_dual(dn_entity1,
                                            entity1_distrib[i_rank],
                                            _dentity1_entity3_idx,
                                            dentity1_entity3_n,
                                            _dentity1_entity3);

  /*
   * Realloc
   */
  *dentity1_entity3 = (PDM_g_num_t *) realloc(*dentity1_entity3, sizeof(PDM_g_num_t) * _dentity1_entity3_idx[dn_entity1] );

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
 const PDM_g_num_t     *entity2_distrib,
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

  PDM_g_num_t shift_g = entity1_distrib[i_rank]; // Entre 1 et N

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

  PDM_g_num_t* entity2_distrib_ptb = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_rank+1));
  for(int i = 0; i < n_rank+1; ++i){
    entity2_distrib_ptb[i] = entity2_distrib[i] - 1;
  }

  // PDM_log_trace_array_long(ln_to_gn, dentity1_entity2_idx[dn_entity1], "ln_to_gn::");
  // PDM_log_trace_array_long(gnum,  dentity1_entity2_idx[dn_entity1], "gnum::");
  /*
   * In order to revert the conncectivty we use the global numbering property
   */
  PDM_part_to_block_t *ptb =
   PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE,
                             1.,
                            &ln_to_gn,
                             entity2_distrib_ptb,
                 (int *)    &dentity1_entity2_idx[dn_entity1],
                             1,
                             comm);

  int* send_stri = (int *) malloc(dentity1_entity2_idx[dn_entity1] * sizeof(int));
  for (int i = 0; i < dentity1_entity2_idx[dn_entity1]; i++) {
    send_stri[i] = 1;
  }

  // if(is_signed) {
  //   free(ln_to_gn);
  // }

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

  int dn_entity2_recv = PDM_part_to_block_n_elt_block_get(ptb);

  /*
   * Free
   */
  PDM_part_to_block_free(ptb);
  free(entity2_distrib_ptb);
  free(gnum);
  free(send_stri);

  /*
   * Allocate
   */
  *dentity2_entity1_idx = (int *) malloc( (dn_entity2_recv + 1) * sizeof(int));
  int* _dentity2_entity1_idx = *dentity2_entity1_idx;

  printf("blk_size       ::%i\n", blk_size       );
  printf("dn_entity2_recv::%i\n", dn_entity2_recv);

  // PDM_log_trace_array_long(dentity2_entity1_n, dn_entity2_recv, "Before : dentity2_entity1_n::");
  // PDM_log_trace_array_long(recv_data, blk_size, "Before : recv_data::");

  PDM_para_graph_compress_connectivity(dn_entity2_recv,
                                        _dentity2_entity1_idx,
                                        dentity2_entity1_n,
                                        recv_data);
  // printf("*dentity2_entity1_idx[dn_entity2_recv]       ::%i\n", _dentity2_entity1_idx[dn_entity2_recv]       );

  *dentity2_entity1 = recv_data;

  /*
   * Realloc
   */
  *dentity2_entity1 = realloc(*dentity2_entity1, _dentity2_entity1_idx[dn_entity2_recv] * sizeof(PDM_g_num_t));
  // PDM_g_num_t* _dentity2_entity1 = *dentity2_entity1;

  // PDM_log_trace_array_int (_dentity2_entity1_idx, dn_entity2_recv+1         , "_dentity2_entity1_idx::");
  // PDM_log_trace_array_long(recv_data, _dentity2_entity1_idx[dn_entity2_recv], "recv_data::");

  free(dentity2_entity1_n);
}

#ifdef  __cplusplus
}
#endif
