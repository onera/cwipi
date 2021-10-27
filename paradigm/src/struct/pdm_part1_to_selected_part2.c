/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part1_to_selected_part2.h"
#include "pdm_part1_to_selected_part2_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_gnum_location.h"
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/


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
 *
 * \brief Free the asynchronous properties of an exchange 
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   request       Request
 * 
 */

static void
_free_async_send_exch
(
 PDM_part1_to_selected_part2_t *ptp,
 const int                     request
)
{
  if (ptp->async_send_request[request] != NULL) {
    free (ptp->async_send_request[request]);
  }
  ptp->async_send_request[request] = NULL;

  ptp->async_send_s_data[request]     = -1;
  ptp->async_send_cst_stride[request] = -1;
  ptp->async_send_tag[request]        = -1;
  if (ptp->async_send_buffer[request] != NULL) {
    free (ptp->async_send_buffer[request]);
  }
  ptp->async_send_buffer[request] = NULL;

  if ((ptp->async_n_send_buffer[request] != NULL) && 
      (ptp->async_n_send_buffer[request] != ptp->default_n_send_buffer)) {
    free (ptp->async_n_send_buffer[request]);
  }
  ptp->async_n_send_buffer[request] = NULL;

  if ((ptp->async_i_send_buffer[request] != NULL) && 
      (ptp->async_i_send_buffer[request] != ptp->default_i_send_buffer)) {
    free (ptp->async_i_send_buffer[request]);
  }
  ptp->async_i_send_buffer[request] = NULL;

  ptp->async_send_open[ptp->async_send_n_open++] = request;     
}


/**
 *
 * \brief Free the asynchronous properties of an exchange 
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   request       Request
 * 
 */

static void
_free_async_recv_exch
(
 PDM_part1_to_selected_part2_t *ptp,
 const int                     request
)
{
  if (ptp->async_recv_request[request] != NULL) {
    free (ptp->async_recv_request[request]);
  }
  ptp->async_recv_request[request] = NULL;

  ptp->async_recv_s_data[request]     = -1;
  ptp->async_recv_cst_stride[request] = -1;
  ptp->async_recv_tag[request]        = -1;

  if (ptp->async_recv_buffer[request]  != NULL) {
    free (ptp->async_recv_buffer[request]);
  }
  ptp->async_recv_buffer[request]   = NULL;

  if ((ptp->async_n_recv_buffer[request] != NULL) && 
      (ptp->async_n_recv_buffer[request] != ptp->default_n_recv_buffer)) {
    free (ptp->async_n_recv_buffer[request]);
  }
  ptp->async_n_recv_buffer[request] = NULL;

  if ((ptp->async_i_recv_buffer[request] != NULL) && 
      (ptp->async_i_recv_buffer[request] != ptp->default_i_recv_buffer)) {
    free (ptp->async_i_recv_buffer[request]);
  }
  ptp->async_i_recv_buffer[request]   = NULL;

  ptp->async_recv_part2_data[request] = NULL;   

  ptp->async_recv_open[ptp->async_recv_n_open++] = request;     
}


/**
 *
 * \brief Initialize an exchange
 *
 * \param [in]   ptp           Block to part structure
 *
 */

static void
_check_async_send_alloc
(
 PDM_part1_to_selected_part2_t *ptp
)
{
  if (ptp->async_send_l_array == 0) {
    ptp->async_send_l_array    = 10;
    ptp->async_send_s_data     = malloc (sizeof(size_t) * ptp->async_send_l_array);
    ptp->async_send_cst_stride = malloc (sizeof(int) * ptp->async_send_l_array);
    ptp->async_send_tag        = malloc (sizeof(int) * ptp->async_send_l_array);
    ptp->async_send_request    = malloc (sizeof(PDM_MPI_Request *) * ptp->async_send_l_array);
    ptp->async_send_buffer     = malloc (sizeof(unsigned char *) * ptp->async_send_l_array);
    ptp->async_n_send_buffer   = malloc (sizeof(int *) * ptp->async_send_l_array);
    ptp->async_i_send_buffer   = malloc (sizeof(int *) * ptp->async_send_l_array);
    ptp->async_send_open       = malloc (sizeof(int) * ptp->async_send_l_array);

    for (int i = 0; i < ptp->async_send_l_array; i++) {
      ptp->async_send_open[ptp->async_send_n_open++] = ptp->async_send_l_array -1 - i; 
      ptp->async_send_s_data[i]     = -1;
      ptp->async_send_cst_stride[i] = -1;
      ptp->async_send_tag[i]        = -1;
      ptp->async_send_request[i]    = NULL;
      ptp->async_send_buffer[i]     = NULL;
      ptp->async_n_send_buffer[i]   = NULL;
      ptp->async_i_send_buffer[i]   = NULL;
    }
  }

  if (ptp->async_send_n_open == 0) {
    const int pre_val = ptp->async_send_l_array; 
    ptp->async_send_l_array *= 2;
    ptp->async_send_open       = realloc (ptp->async_send_open       , sizeof(int) * ptp->async_send_l_array);
    ptp->async_send_s_data     = realloc (ptp->async_send_s_data     , sizeof(size_t) * ptp->async_send_l_array);
    ptp->async_send_cst_stride = realloc (ptp->async_send_cst_stride , sizeof(int) * ptp->async_send_l_array);
    ptp->async_send_tag        = realloc (ptp->async_send_tag        , sizeof(int) * ptp->async_send_l_array);
    ptp->async_send_request    = realloc (ptp->async_send_request    , sizeof(PDM_MPI_Request *) * ptp->async_send_l_array);
    ptp->async_send_buffer     = realloc (ptp->async_send_buffer     , sizeof(unsigned char *) * ptp->async_send_l_array);
    ptp->async_n_send_buffer   = realloc (ptp->async_n_send_buffer   , sizeof(int *) * ptp->async_send_l_array);
    ptp->async_i_send_buffer   = realloc (ptp->async_i_send_buffer   , sizeof(int *) * ptp->async_send_l_array);

    for (int i = pre_val; i < ptp->async_send_l_array; i++) {
      ptp->async_send_open[ptp->async_send_n_open++] = i;
      ptp->async_send_s_data[i]     = -1;
      ptp->async_send_cst_stride[i] = -1;
      ptp->async_send_tag[i]        = -1;
      ptp->async_send_request[i]    = NULL;
      ptp->async_send_buffer[i]     = NULL;
      ptp->async_n_send_buffer[i]   = NULL;
      ptp->async_i_send_buffer[i]   = NULL;
    }
  }
}


/**
 *
 * \brief Initialize an exchange
 *
 * \param [in]   ptp           Block to part structure
 *
 */

static void
_check_async_recv_alloc
(
 PDM_part1_to_selected_part2_t *ptp
)
{
  if (ptp->async_recv_l_array == 0) {
    ptp->async_recv_l_array    = 10;
    ptp->async_recv_s_data     = malloc (sizeof(size_t) * ptp->async_recv_l_array);
    ptp->async_recv_cst_stride = malloc (sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_tag        = malloc (sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_request    = malloc (sizeof(PDM_MPI_Request *) * ptp->async_recv_l_array);
    ptp->async_recv_buffer     = malloc (sizeof(unsigned char *) * ptp->async_recv_l_array);
    ptp->async_n_recv_buffer   = malloc (sizeof(int *) * ptp->async_recv_l_array);
    ptp->async_i_recv_buffer   = malloc (sizeof(int *) * ptp->async_recv_l_array);
    ptp->async_recv_open       = malloc (sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_part2_data = malloc (sizeof(void *) * ptp->async_recv_l_array);   

    for (int i = 0; i < ptp->async_recv_l_array; i++) {
      ptp->async_recv_open[ptp->async_recv_n_open++] = ptp->async_recv_l_array -1 - i; 
      ptp->async_recv_s_data[i]     = -1;
      ptp->async_recv_cst_stride[i] = -1;
      ptp->async_recv_tag[i]        = -1;
      ptp->async_recv_request[i]    = NULL;
      ptp->async_recv_buffer[i]     = NULL;
      ptp->async_n_recv_buffer[i]   = NULL;
      ptp->async_i_recv_buffer[i]   = NULL;
      ptp->async_recv_part2_data[i] = NULL;
    }
  }

  if (ptp->async_recv_n_open == 0) {
    const int pre_val = ptp->async_recv_l_array; 
    ptp->async_recv_l_array *= 2;
    ptp->async_recv_open       = realloc (ptp->async_recv_open       , sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_s_data     = realloc (ptp->async_recv_s_data     , sizeof(size_t) * ptp->async_recv_l_array);
    ptp->async_recv_cst_stride = realloc (ptp->async_recv_cst_stride , sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_tag        = realloc (ptp->async_recv_tag        , sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_request    = realloc (ptp->async_recv_request    , sizeof(PDM_MPI_Request *) * ptp->async_recv_l_array);
    ptp->async_recv_buffer     = realloc (ptp->async_recv_buffer     , sizeof(unsigned char *) * ptp->async_recv_l_array);
    ptp->async_n_recv_buffer   = realloc (ptp->async_n_recv_buffer   , sizeof(int *) * ptp->async_recv_l_array);
    ptp->async_i_recv_buffer   = realloc (ptp->async_i_recv_buffer   , sizeof(int *) * ptp->async_recv_l_array);
    ptp->async_recv_part2_data = realloc (ptp->async_recv_part2_data , sizeof(void *) * ptp->async_recv_l_array);   

    for (int i = pre_val; i < ptp->async_recv_l_array; i++) {
      ptp->async_recv_open[ptp->async_recv_n_open++] = i;
      ptp->async_recv_s_data[i]     = -1;
      ptp->async_recv_cst_stride[i] = -1;
      ptp->async_recv_tag[i]        = -1;
      ptp->async_recv_request[i]    = NULL;
      ptp->async_recv_buffer[i]     = NULL;
      ptp->async_n_recv_buffer[i]   = NULL;
      ptp->async_i_recv_buffer[i]   = NULL;
      ptp->async_recv_part2_data[i] = NULL;
    }
  }
}

/**
 *
 * \brief Init an exchange 
 *
 * \param [in]   ptp           Block to part structure
 * 
 */

static int
_find_open_async_send_exch
(
 PDM_part1_to_selected_part2_t *ptp
)
{
  _check_async_send_alloc (ptp);

  return ptp->async_send_open[--ptp->async_send_n_open];
}

/**
 *
 * \brief Init an exchange 
 *
 * \param [in]   ptp           Block to part structure
 * 
 */

static int
_find_open_async_recv_exch
(
 PDM_part1_to_selected_part2_t *ptp
)
{
  _check_async_recv_alloc (ptp);

  return ptp->async_recv_open[--ptp->async_recv_n_open];
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a partitions to partitions redistribution
 *
 * \param [in]   gnum_elt1          Element global number (size : \ref n_part1)
 * \param [in]   n_elt1             Local number of elements (size : \ref n_part1)
 * \param [in]   n_part1            Number of partition
 * \param [in]   gnum_elt2          Element global number (size : \ref n_part2)
 * \param [in]   n_elt2             Local number of elements (size : \ref n_part2)
 * \param [in]   n_part2            Number of partition
 * \param [in]   send_to_gnum2_idx  Index of data to send to gnum2 from gnum1 
 *                                  (for each part size : \ref n_elt1+1) 
 * \param [in]   send_to_gnum2      Data to send to gnum2 from gnum1 
 * \param [in]   comm               MPI communicator
 *
 * \return   Initialized \ref PDM_part1_to_selected_part2 instance
 *
 */

PDM_part1_to_selected_part2_t *
PDM_part1_to_selected_part2_create
(
 const PDM_g_num_t   **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t   **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **selected_part2_idx,
 const PDM_g_num_t   **selected_part2,
 const PDM_MPI_Comm    comm
)
{
  PDM_part1_to_selected_part2_t *ptp =
    (PDM_part1_to_selected_part2_t *) malloc (sizeof(PDM_part1_to_selected_part2_t));

  /* Init */

  ptp->n_part1                  = n_part1;    
  ptp->gnum_elt1                = gnum_elt1;
  ptp->n_elt1                   = n_elt1;

  ptp->n_part2                  = n_part2;                
  ptp->gnum_elt2                = gnum_elt2;              
  ptp->n_elt2                   = n_elt2;                

  ptp->selected_part2_idx       = selected_part2_idx;    
  ptp->selected_part2           = selected_part2;        
  ptp->comm                     = comm;

  PDM_MPI_Comm_size (comm, &(ptp->n_rank));
  PDM_MPI_Comm_rank (comm, &(ptp->my_rank));
  int my_rank = ptp->my_rank;
  int n_rank = ptp->n_rank;


  ptp->n_ref_gnum2                = NULL;         
  ptp->ref_gnum2                  = NULL;          
  ptp->n_unref_gnum2              = NULL;      
  ptp->unref_gnum2                = NULL;        
  ptp->gnum1_come_from_idx        = NULL;
  ptp->gnum1_come_from            = NULL;    

  ptp->gnum1_to_send_buffer       = NULL;  
  ptp->recv_buffer_to_ref_gnum2   = NULL;  

  ptp->default_n_send_buffer      = malloc (sizeof(int) * n_rank);
  ptp->default_i_send_buffer      = malloc (sizeof(int) * (n_rank + 1));
  ptp->default_i_send_buffer[0]   = 0;
  ptp->default_n_recv_buffer      = malloc (sizeof(int) * n_rank);
  ptp->default_i_recv_buffer      = malloc (sizeof(int) * (n_rank + 1));
  ptp->default_i_recv_buffer[0]   = 0;

  for (int i = 0; i < n_rank; i++) {
    ptp->default_n_send_buffer[i] = 0;
    ptp->default_n_recv_buffer[i] = 0;
  }

  ptp->n_active_rank_send         = 1;
  ptp->active_rank_send           = NULL;
  
  ptp->n_active_rank_recv         = 1;
  ptp->active_rank_recv           = NULL;

  ptp->async_send_n_open          = 0;           
  ptp->async_send_open            = NULL;           
  ptp->async_send_l_array         = 0;          
  ptp->async_send_s_data          = NULL;          
  ptp->async_send_cst_stride      = NULL;      
  ptp->async_send_tag             = NULL;      
  ptp->async_send_request         = NULL;    
  ptp->async_send_buffer          = NULL;     
  ptp->async_n_send_buffer        = NULL;   
  ptp->async_i_send_buffer        = NULL;   

  ptp->async_recv_n_open          = 0;           
  ptp->async_recv_open            = NULL;           
  ptp->async_recv_l_array         = 0;          
  ptp->async_recv_s_data          = NULL;          
  ptp->async_recv_cst_stride      = NULL;      
  ptp->async_recv_tag             = NULL;      
  ptp->async_recv_request         = NULL;    
  ptp->async_recv_buffer          = NULL;     
  ptp->async_n_recv_buffer        = NULL;   
  ptp->async_i_recv_buffer        = NULL;
  ptp->async_recv_part2_data      = NULL;   

  /* 1 - gnum_location in 2 1D array selected_part2_rank selected_part2_part   selected_part2_part elt*/

  PDM_gnum_location_t *gl = PDM_gnum_location_create (n_part2, n_part1, comm);

  for (int i = 0; i < n_part2; i++) {
    PDM_gnum_location_elements_set (gl, i, n_elt2[i], gnum_elt2[i]);
  }

  for (int i = 0; i < n_part1; i++) {
    PDM_gnum_location_requested_elements_set (gl, i, selected_part2_idx[i][n_elt1[i]], selected_part2[i]);
  }

  PDM_gnum_location_compute(gl);

  int n_total_elt = 0;
  for (int i = 0; i < n_part1; i++) {
    int *location_selected_part2_idx;
    int *location_selected_part2;
    PDM_gnum_location_get (gl,
                           i,
                           &location_selected_part2_idx,
                           &location_selected_part2);
    n_total_elt += location_selected_part2_idx[selected_part2_idx[i][n_elt1[i]]];
  }

  n_total_elt /= 3;

  int *merge_selected_part2_rank2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_selected_part2_part2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_selected_part2_lnum2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_selected_part2_rank1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_selected_part2_part1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_selected_part2_lnum1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_selected_part2_addr1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *order                      = (int *) malloc (sizeof(int) * n_total_elt);
  int *i_send_buffer              = (int *) malloc (sizeof(int) * n_total_elt);
  int *n_selected_part2_rank      = (int *) malloc (sizeof(int) * n_rank);
  int *idx_selected_part2_rank    = (int *) malloc (sizeof(int) * (n_rank + 1));

  for (int i = 0; i < n_rank; i++) {
    n_selected_part2_rank[i] = 0;
  }
  idx_selected_part2_rank[0] = 0;

  n_total_elt = 0;
  for (int i = 0; i < n_part1; i++) {
    int *location_selected_part2_idx;
    int *location_selected_part2;
    PDM_gnum_location_get (gl,
                           i,
                           &location_selected_part2_idx,
                           &location_selected_part2);
  
    for (int j = 0; j < n_elt1[i]; j++) {

      for (int k = location_selected_part2_idx[selected_part2_idx[i][j]]/3; 
               k < location_selected_part2_idx[selected_part2_idx[i][j+1]]/3; k++) {
        int i_rank2 = location_selected_part2[3*k];
        n_selected_part2_rank[i_rank2]++;
        merge_selected_part2_rank2[n_total_elt] = i_rank2;
        merge_selected_part2_part2[n_total_elt] = location_selected_part2[3*k+1];
        merge_selected_part2_lnum2[n_total_elt] = location_selected_part2[3*k+2]-1;
        merge_selected_part2_rank1[n_total_elt] = my_rank;
        merge_selected_part2_part1[n_total_elt] = i;
        merge_selected_part2_lnum1[n_total_elt] = j; 
        merge_selected_part2_addr1[n_total_elt] = k; 
        order[n_total_elt]                      = n_total_elt;
        n_total_elt++;
      }
    }
  }

  PDM_gnum_location_free (gl, 0);

  for (int i = 0; i < n_rank; i++) {
    idx_selected_part2_rank[i+1] = n_selected_part2_rank[i] + 
                                   idx_selected_part2_rank[i];
  }

  /* 2 - Sort selected_part2_rank selected_part2_part according successively rank, part and elemt  */

  PDM_sort_int (merge_selected_part2_rank2, order, n_total_elt);

  int *_merge_selected_part2_rank2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_selected_part2_part2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_selected_part2_lnum2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_selected_part2_rank1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_selected_part2_part1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_selected_part2_lnum1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_selected_part2_addr1 = (int *) malloc (sizeof(int) * n_total_elt);

  for (int i = 0; i < n_total_elt; i++) {
    i_send_buffer[i] = -2;
    _merge_selected_part2_part2[i] = merge_selected_part2_part2[order[i]];
    _merge_selected_part2_lnum2[i] = merge_selected_part2_lnum2[order[i]];
    _merge_selected_part2_rank1[i] = merge_selected_part2_rank1[order[i]];
    _merge_selected_part2_part1[i] = merge_selected_part2_part1[order[i]];
    _merge_selected_part2_lnum1[i] = merge_selected_part2_lnum1[order[i]]; 
    _merge_selected_part2_addr1[i] = merge_selected_part2_addr1[order[i]]; 
  }

  int *_tmp_merge_selected_part2_rank2 = merge_selected_part2_rank2;
  int *_tmp_merge_selected_part2_part2 = merge_selected_part2_part2;
  int *_tmp_merge_selected_part2_lnum2 = merge_selected_part2_lnum2;
  int *_tmp_merge_selected_part2_rank1 = merge_selected_part2_rank1;
  int *_tmp_merge_selected_part2_part1 = merge_selected_part2_part1;
  int *_tmp_merge_selected_part2_lnum1 = merge_selected_part2_lnum1;
  int *_tmp_merge_selected_part2_addr1 = merge_selected_part2_addr1;

  merge_selected_part2_rank2 = _merge_selected_part2_rank2;
  merge_selected_part2_part2 = _merge_selected_part2_part2;
  merge_selected_part2_lnum2 = _merge_selected_part2_lnum2;
  merge_selected_part2_rank1 = _merge_selected_part2_rank1;
  merge_selected_part2_part1 = _merge_selected_part2_part1;
  merge_selected_part2_lnum1 = _merge_selected_part2_lnum1;
  merge_selected_part2_addr1 = _merge_selected_part2_addr1;

  _merge_selected_part2_rank2 = _tmp_merge_selected_part2_rank2;
  _merge_selected_part2_part2 = _tmp_merge_selected_part2_part2;
  _merge_selected_part2_lnum2 = _tmp_merge_selected_part2_lnum2;
  _merge_selected_part2_rank1 = _tmp_merge_selected_part2_rank1;
  _merge_selected_part2_part1 = _tmp_merge_selected_part2_part1;
  _merge_selected_part2_lnum1 = _tmp_merge_selected_part2_lnum1;
  _merge_selected_part2_addr1 = _tmp_merge_selected_part2_addr1;

  int n_part2_max = 0;
  PDM_MPI_Allreduce ((int* )&n_part2, &n_part2_max, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  int *n_elt_part = malloc (sizeof(int) * n_part2_max);
  int *idx_elt_part = malloc (sizeof(int) * (n_part2_max + 1));
  idx_elt_part[0] = 0;

  int cpt_buff = -1;
  for (int i = 0; i < n_rank; i++) {
    int n_elt_rank = n_selected_part2_rank[i];
    for (int j = 0; j < n_elt_rank; j++) {
      order[j] = j;
    }

    for (int j = 0; j < n_part2_max; j++) {
      n_elt_part[j] = 0;
    }

    int *rank_i_send_buffer              = i_send_buffer              + idx_selected_part2_rank[i];  
    int *rank_merge_selected_part2_part2 = merge_selected_part2_part2 + idx_selected_part2_rank[i];
    int *rank_merge_selected_part2_lnum2 = merge_selected_part2_lnum2 + idx_selected_part2_rank[i];
    int *rank_merge_selected_part2_rank1 = merge_selected_part2_rank1 + idx_selected_part2_rank[i];
    int *rank_merge_selected_part2_part1 = merge_selected_part2_part1 + idx_selected_part2_rank[i];
    int *rank_merge_selected_part2_lnum1 = merge_selected_part2_lnum1 + idx_selected_part2_rank[i];
    int *rank_merge_selected_part2_addr1 = merge_selected_part2_addr1 + idx_selected_part2_rank[i];

    int *_rank_merge_selected_part2_part2 = _merge_selected_part2_part2 + idx_selected_part2_rank[i];
    int *_rank_merge_selected_part2_lnum2 = _merge_selected_part2_lnum2 + idx_selected_part2_rank[i];
    int *_rank_merge_selected_part2_rank1 = _merge_selected_part2_rank1 + idx_selected_part2_rank[i];
    int *_rank_merge_selected_part2_part1 = _merge_selected_part2_part1 + idx_selected_part2_rank[i];
    int *_rank_merge_selected_part2_lnum1 = _merge_selected_part2_lnum1 + idx_selected_part2_rank[i];
    int *_rank_merge_selected_part2_addr1 = _merge_selected_part2_addr1 + idx_selected_part2_rank[i];

    PDM_sort_int (_rank_merge_selected_part2_part2, order, n_elt_rank);

    int _max_part = 0;
    for (int k = 0; k < n_elt_rank; k++) {
      int i_part = rank_merge_selected_part2_part2[order[k]];
      n_elt_part[i_part]++;
      _max_part = PDM_MAX (_max_part, i_part+1);
      _rank_merge_selected_part2_part2[k] = i_part;
      _rank_merge_selected_part2_lnum2[k] = rank_merge_selected_part2_lnum2[order[k]];
      _rank_merge_selected_part2_rank1[k] = rank_merge_selected_part2_rank1[order[k]];
      _rank_merge_selected_part2_part1[k] = rank_merge_selected_part2_part1[order[k]];
      _rank_merge_selected_part2_lnum1[k] = rank_merge_selected_part2_lnum1[order[k]]; 
      _rank_merge_selected_part2_addr1[k] = rank_merge_selected_part2_addr1[order[k]]; 
    }

    for (int k = 0; k < _max_part; k++) {
      idx_elt_part[k+1] = idx_elt_part[k] + n_elt_part[k]; 
    }

    for (int k1 = 0; k1 < _max_part; k1++) {

      int _n_elt_part = n_elt_part[k1];

      for (int j = 0; j < _n_elt_part; j++) {
        order[j] = j;
      }

      int *_part_rank_merge_selected_part2_part2 = _rank_merge_selected_part2_part2 + idx_elt_part[k1];
      int *_part_rank_merge_selected_part2_lnum2 = _rank_merge_selected_part2_lnum2 + idx_elt_part[k1];
      int *_part_rank_merge_selected_part2_rank1 = _rank_merge_selected_part2_rank1 + idx_elt_part[k1];
      int *_part_rank_merge_selected_part2_part1 = _rank_merge_selected_part2_part1 + idx_elt_part[k1];
      int *_part_rank_merge_selected_part2_lnum1 = _rank_merge_selected_part2_lnum1 + idx_elt_part[k1];
      int *_part_rank_merge_selected_part2_addr1 = _rank_merge_selected_part2_addr1 + idx_elt_part[k1];

      int *part_rank_merge_selected_part2_part2 = rank_merge_selected_part2_part2 + idx_elt_part[k1];
      int *part_rank_merge_selected_part2_lnum2 = rank_merge_selected_part2_lnum2 + idx_elt_part[k1];
      int *part_rank_merge_selected_part2_rank1 = rank_merge_selected_part2_rank1 + idx_elt_part[k1];
      int *part_rank_merge_selected_part2_part1 = rank_merge_selected_part2_part1 + idx_elt_part[k1];
      int *part_rank_merge_selected_part2_lnum1 = rank_merge_selected_part2_lnum1 + idx_elt_part[k1];
      int *part_rank_merge_selected_part2_addr1 = rank_merge_selected_part2_addr1 + idx_elt_part[k1];
      int *part_rank_i_send_buffer              = rank_i_send_buffer              + idx_elt_part[k1];  

      PDM_sort_int (_part_rank_merge_selected_part2_lnum2, order, _n_elt_part);

      int pre_val = -1;
      for (int k2 = 0; k2 < _n_elt_part; k2++) {
        part_rank_merge_selected_part2_part2[k2] = _part_rank_merge_selected_part2_part2[k2];
        part_rank_merge_selected_part2_lnum2[k2] = _part_rank_merge_selected_part2_lnum2[k2];
        part_rank_merge_selected_part2_rank1[k2] = _part_rank_merge_selected_part2_rank1[order[k2]];
        part_rank_merge_selected_part2_part1[k2] = _part_rank_merge_selected_part2_part1[order[k2]];
        part_rank_merge_selected_part2_lnum1[k2] = _part_rank_merge_selected_part2_lnum1[order[k2]]; 
        part_rank_merge_selected_part2_addr1[k2] = _part_rank_merge_selected_part2_addr1[order[k2]]; 

        if (pre_val != part_rank_merge_selected_part2_lnum2[k2]) {
          cpt_buff++;   
          ptp->default_n_send_buffer[i]++;
        }
        part_rank_i_send_buffer[k2] = cpt_buff;

      }

    }

  }

  free (idx_elt_part);

  /* 3 - Define Default_n_send_buffer and  Default_i_send_buffer */

  for (int i = 0; i < n_rank; i++) {
    ptp->default_i_send_buffer[i+1] = ptp->default_i_send_buffer[i] + ptp->default_n_send_buffer[i];
  }

  if (1 == 0) {
    printf ("ptp->default_i_send_buffer :");
    for (int i = 0; i < n_rank + 1; i++) {
      printf (" %d", ptp->default_i_send_buffer[i]);
    }
    printf("\n");

    printf ("ptp->default_n_send_buffer :");
    for (int i = 0; i < n_rank; i++) {
      printf (" %d", ptp->default_n_send_buffer[i]);
    }
    printf("\n");
  }

  /* 4 - Define gnum1_to_send_buffer */

  ptp->gnum1_to_send_buffer_idx = malloc (sizeof (int*) * n_part1);
  ptp->gnum1_to_send_buffer     = malloc (sizeof (int*) * n_part1);
  int **gnum1_to_send_buffer_n    = malloc (sizeof (int*) * n_part1);

  for (int i = 0; i < n_part1; i++) {
    ptp->gnum1_to_send_buffer_idx[i] = malloc (sizeof (int) * (selected_part2_idx[i][n_elt1[i]]+1));    
    gnum1_to_send_buffer_n[i] = malloc (sizeof (int) * selected_part2_idx[i][n_elt1[i]]);    
    ptp->gnum1_to_send_buffer_idx[i][0] = 0;
    for (int k = 0; k < selected_part2_idx[i][n_elt1[i]]; k++) {
      ptp->gnum1_to_send_buffer_idx[i][k+1] = 0;
      gnum1_to_send_buffer_n[i][k] = 0;
    }
  }

  for (int i = 0; i < n_total_elt; i++) {
    int ipart1 = merge_selected_part2_part1[i];
    int addr1  = merge_selected_part2_addr1[i];
    ptp->gnum1_to_send_buffer_idx[ipart1][addr1+1]++;
  }


  for (int i = 0; i < n_part1; i++) {
    for (int k = 0; k < selected_part2_idx[i][n_elt1[i]]; k++) {
      ptp->gnum1_to_send_buffer_idx[i][k+1] += ptp->gnum1_to_send_buffer_idx[i][k] ;
    }
  }

  if (1 == 0) {
    for (int i = 0; i < n_part1; i++) {
      printf ("gnum1_to_send_buffer_idx 2 %d :", i);
      for (int k = 0; k < n_elt1[i] + 1; k++) {
        printf (" %d",ptp->gnum1_to_send_buffer_idx[i][k]);
      }
      printf("\n");
    }
  }

  for (int i = 0; i < n_part1; i++) {
    int size = ptp->gnum1_to_send_buffer_idx[i][selected_part2_idx[i][n_elt1[i]]];
    ptp->gnum1_to_send_buffer[i] = malloc (sizeof (int) * size);
    for (int k = 0; k < size; k++) {
      ptp->gnum1_to_send_buffer[i][k] = -1;
    }
  }

  for (int i = 0; i < n_total_elt; i++) {
    int ipart1 = merge_selected_part2_part1[i];
    int addr1  = merge_selected_part2_addr1[i];
    int idx    = i_send_buffer[i];
    int idx2   = ptp->gnum1_to_send_buffer_idx[ipart1][addr1] + 
                 gnum1_to_send_buffer_n[ipart1][addr1]++; 
    ptp->gnum1_to_send_buffer[ipart1][idx2] = idx;
  }

  if (1 == 0) {
    for (int i = 0; i < n_part1; i++) {
      for (int j = 0; j < selected_part2_idx[i][n_elt1[i]]; j++) { 
        for (int k = ptp->gnum1_to_send_buffer_idx[i][j]; 
                 k < ptp->gnum1_to_send_buffer_idx[i][j+1];
                 k++) {
          printf ("2 - %d %d : %d\n", i, j, ptp->gnum1_to_send_buffer[i][k]);
        }
      }
    }
  }

  for (int i = 0; i < n_part1; i++) {
    free (gnum1_to_send_buffer_n[i]);
  }
  free(gnum1_to_send_buffer_n); 



  fflush(stdout);

  free (n_elt_part);
  free (order);

  free (_merge_selected_part2_rank2);
  free (_merge_selected_part2_part2);
  free (_merge_selected_part2_lnum2);
  free (_merge_selected_part2_rank1);
  free (_merge_selected_part2_part1);
  free (_merge_selected_part2_lnum1);
  free (_merge_selected_part2_addr1);

  free (idx_selected_part2_rank);
  free (n_selected_part2_rank);


  /* 5 - Define Default_n_recv_buffer and  Default_i_recv_buffer */

  PDM_MPI_Alltoall (ptp->default_n_send_buffer, 1, PDM_MPI_INT,
                    ptp->default_n_recv_buffer, 1, PDM_MPI_INT,
                    comm);
  
  for (int i = 0; i < n_rank; i++) {
    ptp->default_i_recv_buffer[i+1] = ptp->default_i_recv_buffer[i] + ptp->default_n_recv_buffer[i];
  }

  for (int i = 0; i < n_rank; i++) {
    ptp->default_n_send_buffer[i]   *= 4;
    ptp->default_i_send_buffer[i+1] *= 4;
    ptp->default_n_recv_buffer[i]   *= 4;
    ptp->default_i_recv_buffer[i+1] *= 4;
  }

  int *int_s_buff = malloc (sizeof(int) * ptp->default_i_send_buffer[n_rank]);  
  int *int_r_buff = malloc (sizeof(int) * ptp->default_i_recv_buffer[n_rank]);  

  PDM_g_num_t *gnum_s_buff = malloc (sizeof(PDM_g_num_t) * ptp->default_i_send_buffer[n_rank]);  
  PDM_g_num_t *gnum_r_buff = malloc (sizeof(PDM_g_num_t) * ptp->default_i_recv_buffer[n_rank]);  

  for (int i = 0; i < n_total_elt; i++) {
    int ipart1      = merge_selected_part2_part1[i];
    int ielt1       = merge_selected_part2_lnum1[i];
    int ipart2      = merge_selected_part2_part2[i];
    int ielt2       = merge_selected_part2_lnum2[i];
    int idx         = i_send_buffer[i];
    if (idx >= 0) {
      int_s_buff[4 * idx    ] = ipart1;
      int_s_buff[4 * idx + 1] = ielt1;
      int_s_buff[4 * idx + 2] = ipart2;
      int_s_buff[4 * idx + 3] = ielt2;
      gnum_s_buff[idx]        = gnum_elt1[ipart1][ielt1];
    }
  }

  free (merge_selected_part2_rank2);
  free (merge_selected_part2_part2);
  free (merge_selected_part2_lnum2);
  free (merge_selected_part2_rank1);
  free (merge_selected_part2_part1);
  free (merge_selected_part2_lnum1);
  free (merge_selected_part2_addr1);
  free (i_send_buffer);

  PDM_MPI_Alltoallv (int_s_buff, ptp->default_n_send_buffer, ptp->default_i_send_buffer, PDM_MPI_INT,
                     int_r_buff, ptp->default_n_recv_buffer, ptp->default_i_recv_buffer, PDM_MPI_INT,
                     comm);

  for (int i = 0; i < n_rank; i++) {
    ptp->default_n_send_buffer[i]   /= 4;
    ptp->default_i_send_buffer[i+1] /= 4;
    ptp->default_n_recv_buffer[i]   /= 4;
    ptp->default_i_recv_buffer[i+1] /= 4;
  }

  PDM_MPI_Alltoallv (gnum_s_buff, ptp->default_n_send_buffer, ptp->default_i_send_buffer, PDM__PDM_MPI_G_NUM,
                     gnum_r_buff, ptp->default_n_recv_buffer, ptp->default_i_recv_buffer, PDM__PDM_MPI_G_NUM,
                     comm);

  /* 6 - Build the arrays for the reveived view */

  ptp->n_ref_gnum2              = malloc (sizeof (int) * n_part2);         
  ptp->ref_gnum2                = malloc (sizeof (int *) * n_part2);          
  ptp->n_unref_gnum2            = malloc (sizeof (int) * n_part2);      
  ptp->unref_gnum2              = malloc (sizeof (int *) * n_part2);        
  ptp->gnum1_come_from_idx      = malloc (sizeof (int *) * n_part2);
  ptp->gnum1_come_from          = malloc (sizeof (PDM_g_num_t *) * n_part2);
  ptp->recv_buffer_to_ref_gnum2 = malloc (sizeof (int *) * n_part2);   

  //ptp->gnum1_to_send_buffer     = NULL;  
  //ptp->recv_buffer_to_ref_gnum2 = NULL;  

  int **tag_elt2 = malloc (sizeof (int *) * n_part2);

  for (int i = 0; i < n_part2; i++) {
    tag_elt2[i] = malloc (sizeof (int) * n_elt2[i]);
    for (int j = 0; j < n_elt2[i]; j++) {
      tag_elt2[i][j] = 0;
    }
  }

  for (int i = 0; i < n_rank; i++) {
    //int iproc1 = i;
    for (int j = ptp->default_i_recv_buffer[i]; j < ptp->default_i_recv_buffer[i+1]; j++) {
      //int recv_ipart1 = int_r_buff[4 * j + 0];
      //int recv_ielt1  = int_r_buff[4 * j + 1];
      //int recv_gnum1  = gnum_r_buff[j];
      int recv_ipart2 = int_r_buff[4 * j + 2];
      int recv_ielt2  = int_r_buff[4 * j + 3];
      tag_elt2[recv_ipart2][recv_ielt2]++;
    }
  }

  int **ielt_to_ref = malloc (sizeof (int *) * n_part2);
  for (int i = 0; i < n_part2; i++) {
    ielt_to_ref[i] = malloc (sizeof (int) * n_elt2[i]);
    ptp->n_ref_gnum2[i]   = 0;
    ptp->n_unref_gnum2[i] = 0;
    ptp->ref_gnum2[i] = malloc (sizeof (int) * n_elt2[i]);
    ptp->unref_gnum2[i] = malloc (sizeof (int) * n_elt2[i]);
    ptp->gnum1_come_from_idx[i] = malloc (sizeof (int) * (n_elt2[i] + 1));
    ptp->gnum1_come_from_idx[i][0] = 0;

    for (int j = 0; j < n_elt2[i]; j++) {
      int _tag = tag_elt2[i][j];
      tag_elt2[i][j] = 0;
      if (_tag > 0) {
        ptp->gnum1_come_from_idx[i][ptp->n_ref_gnum2[i]+1] = _tag;
        ielt_to_ref[i][j] = ptp->n_ref_gnum2[i];
        ptp->ref_gnum2[i][ptp->n_ref_gnum2[i]++] = j+1;
      }
      else {
        ptp->unref_gnum2[i][ptp->n_unref_gnum2[i]++] = j+1;
      }
    }
  
    ptp->ref_gnum2[i]           = realloc (ptp->ref_gnum2[i], sizeof (int) * ptp->n_ref_gnum2[i]);
    ptp->unref_gnum2[i]         = realloc (ptp->unref_gnum2[i], sizeof (int) * ptp->n_unref_gnum2[i]); 
    ptp->gnum1_come_from_idx[i] = realloc (ptp->gnum1_come_from_idx[i], sizeof (int) * (ptp->n_ref_gnum2[i] + 1));
    
    for (int j = 0; j < ptp->n_ref_gnum2[i]; j++) {
      ptp->gnum1_come_from_idx[i][j+1] += ptp->gnum1_come_from_idx[i][j];   
    }

    ptp->gnum1_come_from[i]          = malloc (sizeof (PDM_g_num_t) * ptp->gnum1_come_from_idx[i][ptp->n_ref_gnum2[i]]);
    ptp->recv_buffer_to_ref_gnum2[i] = malloc (sizeof (int) * ptp->gnum1_come_from_idx[i][ptp->n_ref_gnum2[i]]);
  }

  cpt_buff = 0;
  for (int i = 0; i < n_rank; i++) {
    //int iproc1 = i;
    for (int j = ptp->default_i_recv_buffer[i]; j < ptp->default_i_recv_buffer[i+1]; j++) {
      //int recv_ipart1 = int_r_buff[4 * j + 0];
      //int recv_ielt1  = int_r_buff[4 * j + 1];
      int recv_gnum1  = gnum_r_buff[j];
      int recv_ipart2 = int_r_buff[4 * j + 2];
      int recv_ielt2  = int_r_buff[4 * j + 3];
      int iref = ielt_to_ref[recv_ipart2][recv_ielt2];
      int idx = ptp->gnum1_come_from_idx[recv_ipart2][iref] + tag_elt2[recv_ipart2][iref];
      ptp->gnum1_come_from[recv_ipart2][idx] = recv_gnum1;
      ptp->recv_buffer_to_ref_gnum2[recv_ipart2][idx] = cpt_buff++;   
      tag_elt2[recv_ipart2][iref]++;
    }
  }

  /* 8 - Sort and remove duplicate */

  int max_n_gnum1 = 0;
  for (int i = 0; i < n_part2; i++) {
    for (int j = 0; j < ptp->n_ref_gnum2[i]; j++) {
      int idx = ptp->gnum1_come_from_idx[i][j];
      int n_gnum1 = ptp->gnum1_come_from_idx[i][j+1] - idx;
      max_n_gnum1 = PDM_MAX(max_n_gnum1, n_gnum1); 
    }
  }

  if (1 == 0) {
    printf ("ptp->gnum1_come_from 1\n");
    for (int i = 0; i < n_part2; i++) {
      printf(" - ipart : %d\n", i);
      for (int j = 0; j < ptp->n_ref_gnum2[i]; j++) { 
        printf ("%d :", j);
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) { 
          printf (" "PDM_FMT_G_NUM"", ptp->gnum1_come_from[i][k]);
        }
        printf ("\n");
      }
    }

    printf ("ptp->recv_buffer_to_ref_gnum2 1\n");
    for (int i = 0; i < n_part2; i++) {
      printf(" - ipart : %d\n", i);
      for (int j = 0; j < ptp->n_ref_gnum2[i]; j++) { 
        printf ("%d :", j);
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) { 
          printf (" %d", ptp->recv_buffer_to_ref_gnum2[i][k]);
        }
        printf ("\n");
      }
    }
  }
  order = (int *) malloc (sizeof(int) * max_n_gnum1);

  for (int i = 0; i < n_part2; i++) {
    int *_recv_buffer_to_ref_gnum2 = malloc(sizeof(int) * ptp->gnum1_come_from_idx[i][ptp->n_ref_gnum2[i]]);

    for (int j = 0; j < ptp->n_ref_gnum2[i]; j++) { 
      int idx = ptp->gnum1_come_from_idx[i][j];
      int n_gnum1 = ptp->gnum1_come_from_idx[i][j+1] - idx;
      PDM_g_num_t *_gnum1_come_from = ptp->gnum1_come_from[i] + idx;
      int *__recv_buffer_to_ref_gnum2 = _recv_buffer_to_ref_gnum2 + idx;
      int *___recv_buffer_to_ref_gnum2 = ptp->recv_buffer_to_ref_gnum2[i] + idx;
      
      for (int k = 0; k < n_gnum1; k++) { 
        order[k] = k;
      }

      PDM_sort_long (_gnum1_come_from, order, n_gnum1);

      for (int k = 0; k < n_gnum1; k++) { 
        __recv_buffer_to_ref_gnum2[k] = ___recv_buffer_to_ref_gnum2[order[k]];
      }
    }

    free (ptp->recv_buffer_to_ref_gnum2[i]);
    ptp->recv_buffer_to_ref_gnum2[i] = _recv_buffer_to_ref_gnum2;

    int *_old_gnum1_come_from_idx = malloc(sizeof(int) * (ptp->n_ref_gnum2[i] + 1));

    for (int j = 0; j < ptp->n_ref_gnum2[i] + 1; j++) {
      _old_gnum1_come_from_idx[j] = ptp->gnum1_come_from_idx[i][j];

      ptp->gnum1_come_from_idx[i][j] = 0;
    }

    int cpt = 0; 
    for (int j = 0; j < ptp->n_ref_gnum2[i]; j++) {
      int current_val = -1;
      for (int k = _old_gnum1_come_from_idx[j]; k < _old_gnum1_come_from_idx[j+1]; k++) { 
        if (ptp->gnum1_come_from[i][k] != current_val) {
          current_val = ptp->gnum1_come_from[i][k];
          ptp->gnum1_come_from[i][cpt] = ptp->gnum1_come_from[i][k]; 
          ptp->recv_buffer_to_ref_gnum2[i][cpt] = ptp->recv_buffer_to_ref_gnum2[i][k]; 
          cpt++; 
        }
      }
      ptp->gnum1_come_from_idx[i][j+1] = cpt;
    }


    // realloc

    free (_old_gnum1_come_from_idx);

  }

  if (1 == 0) {
    printf ("ptp->gnum1_come_from 2\n");
    for (int i = 0; i < n_part2; i++) {
      printf(" - ipart : %d\n", i);
      for (int j = 0; j < ptp->n_ref_gnum2[i]; j++) { 
        printf ("%d :", j);
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) { 
          printf (" "PDM_FMT_G_NUM"", ptp->gnum1_come_from[i][k]);
        }
        printf ("\n");
      }
    }

    printf ("ptp->recv_buffer_to_ref_gnum2 2\n");
    for (int i = 0; i < n_part2; i++) {
      printf(" - ipart : %d\n", i);
      for (int j = 0; j < ptp->n_ref_gnum2[i]; j++) { 
        printf ("%d :", j);
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) { 
          printf (" %d", ptp->recv_buffer_to_ref_gnum2[i][k]);
        }
        printf ("\n");
      }
    }
  }

  /* 7 - Look for the active ranks */

  ptp->n_active_rank_send = 0;
  ptp->active_rank_send = malloc (sizeof(int) * n_rank);
  for (int i = 0; i < n_rank; i++) {
    if (ptp->default_n_send_buffer[i] > 0) {
      ptp->active_rank_send[ptp->n_active_rank_send++] = i;
    }
  }

  ptp->n_active_rank_recv = 0;
  ptp->active_rank_recv = malloc (sizeof(int) * n_rank);
  for (int i = 0; i < n_rank; i++) {
    if (ptp->default_n_recv_buffer[i] > 0) {
      ptp->active_rank_recv[ptp->n_active_rank_recv++] = i;
    }
  }

  free (order);

  free (int_s_buff);
  free (int_r_buff);

  free (gnum_s_buff);
  free (gnum_r_buff);

  for (int i = 0; i < n_part2; i++) {
    free (ielt_to_ref[i]);
    free (tag_elt2[i]);
  }
  free (tag_elt2);
  free (ielt_to_ref);

  return ptp;

}


PDM_part1_to_selected_part2_t *
PDM_part1_to_selected_part2_create_cf
(
 const PDM_g_num_t    **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t    **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **selected_part2_idx,
 const PDM_g_num_t   **selected_part2,
 const PDM_MPI_Fint    fcomm
)
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(fcomm);
  return PDM_part1_to_selected_part2_create (gnum_elt1, n_elt1, n_part1,
                                            gnum_elt2, n_elt2, n_part2,
                                            selected_part2_idx, selected_part2, _comm);
}


/**
 *
 * \brief Get selected numbers of part2
 *
 * \param [in]   ptp                 Block to part structure
 * \param [out]  n_elt1              Number of gnum1 element  
 * \param [out]  selected_part2_idx  Index of data to send to gnum2 from gnum1 
 *                                  (for each part size : \ref n_elt1+1) 
 * \param [out]  selected_part2      Data to send to gnum2 from gnum1 for each part
 *
 */

void
PDM_part1_to_selected_part2_selected_part2_get
(
 PDM_part1_to_selected_part2_t *ptp,
 int                           *n_elt1,
 int                         ***selected_part2_idx,
 PDM_g_num_t                 ***selected_part2
)
{
  *n_elt1             = ptp->n_elt1;
  *selected_part2_idx = ptp->selected_part2_idx;
  *selected_part2     = ptp->selected_part2;
}


/**
 *
 * \brief Initialize an exchange based on MPI_ialltoall
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   selected_part2_data Data in same order than selected_part2 array
 * \param [out]  ref_part2_data      Data to referenced part2 elements
 * \param [out]  request             Request
 *
 */

void
PDM_part1_to_selected_part2_ialltoall
(
PDM_part1_to_selected_part2_t *ptp,
 const size_t                  s_data,
 const int                     cst_stride,
 void                        **selected_part2_data,
 void                        **ref_part2_data,
 int                          *request
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (s_data);
  PDM_UNUSED (cst_stride);
  PDM_UNUSED (selected_part2_data);
  PDM_UNUSED (ref_part2_data);
  PDM_UNUSED (request);
  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_part1_to_selected_part2_ialltoall not yet implemented\n");
}


/**
 *
 * \brief Wait a asynchronus issend
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part1_to_selected_part2_ialltoall_wait
(
 PDM_part1_to_selected_part2_t *ptp,
 int                           request
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (request);
  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_part1_to_selected_part2_ialltoall_wait not yet implemented\n");
}


/**
 *
 * \brief Initialize an exchange based on MPI_ineighbor_alltoall
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   selected_part2_data Data in same order than selected_part2 array
 * \param [out]  ref_part2_data          Data to referenced part2 elements
 * \param [out]  request             Request
 *
 */

void
PDM_part1_to_selected_part2_ineighbor_alltoall
(
PDM_part1_to_selected_part2_t *ptp,
 const size_t                  s_data,
 const int                     cst_stride,
 void                        **selected_part2_data,
 void                        **ref_part2_data,
 int                          *request
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (s_data);
  PDM_UNUSED (cst_stride);
  PDM_UNUSED (selected_part2_data);
  PDM_UNUSED (ref_part2_data);
  PDM_UNUSED (request);
  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_part1_to_selected_part2_ineighbor_alltoall not yet implemented\n");
}

/**
 *
 * \brief Wait a asynchronus issend
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part1_to_selected_part2_ineighbor_alltoall_wait
(
 PDM_part1_to_selected_part2_t *ptp,
 int                           request
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (request);
  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_part1_to_selected_part2_ineighbor_alltoall_wait not yet implemented\n");
}


/**
 *
 * \brief Get referenced gnum2 elements
 *
 * \param [in]   ptp           Block to part structure
 * \param [out]  n_ref_gnum2   Number of referenced gnum2
 * \param [out]  ref_gnum2     Referenced gnum2
 *
 */

void
PDM_part1_to_selected_part2_ref_gnum2_get
(
 PDM_part1_to_selected_part2_t *ptp,
 int                         **n_ref_gnum2,
 int                        ***ref_gnum2
)
{
  *n_ref_gnum2 = ptp->n_ref_gnum2;
  *ref_gnum2   = ptp->ref_gnum2;
}


/**
 *
 * \brief Get unreferenced gnum2 elements
 *
 * \param [in]   ptp           Block to part structure
 * \param [out]  n_unref_gnum2   Number of referenced gnum2
 * \param [out]  unref_gnum2     Referenced gnum2
 *
 */

void
PDM_part1_to_selected_part2_unref_gnum2_get
(
 PDM_part1_to_selected_part2_t *ptp,
 int                         **n_unref_gnum2,
 int                        ***unref_gnum2
)
{
  *n_unref_gnum2 = ptp->n_unref_gnum2;
  *unref_gnum2   = ptp->unref_gnum2;  
}


/**
 *
 * \brief Get gnum come from gnum1 for each referenced gnum2
 *
 * \param [in]   ptp                 Block to part structure
 * \param [out]  gnum1_come_from_idx Index for gnum1_come_from array (size = \ref n_part2)  
 * \param [out]  gnum1_come_from     Gnum come from gnum1 for each referenced gnum2
 *
 */

void
PDM_part1_to_selected_part2_gnum1_come_from_get
(
 PDM_part1_to_selected_part2_t *ptp,
 int                        ***gnum1_come_from_idx,
 PDM_g_num_t                ***gnum1_come_from
)
{
  *gnum1_come_from_idx = ptp->gnum1_come_from_idx;
  *gnum1_come_from     = ptp->gnum1_come_from;
}


/**
 *
 * \brief Initialize a asynchronus issend
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   selected_part2_data Data in same order than selected_part2 array
 * \param [in]   tag                 Tag of the exchange 
 * \param [out]  request             Request
 *
 */

void
PDM_part1_to_selected_part2_issend
(
 PDM_part1_to_selected_part2_t *ptp,
 const size_t                  s_data,
 const int                     cst_stride,
 void                        **selected_part2_data,
 int                           tag,
 int                          *request
)
{
  unsigned char ** _part1_data = (unsigned char **) selected_part2_data;

  *request = _find_open_async_send_exch (ptp);
  int _request = *request;

  ptp->async_send_s_data[_request]      = s_data;      
  ptp->async_send_cst_stride[_request]  = cst_stride;      
  ptp->async_send_tag[_request]         = tag;      
  ptp->async_send_request[_request]     = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_send);      
  ptp->async_n_send_buffer[_request]  = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_send_buffer[_request]  = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_send_buffer[_request][0] = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_send_buffer[_request][i]   = cst_stride * ptp->default_n_send_buffer[i] * (int) s_data; 
    ptp->async_i_send_buffer[_request][i+1] = cst_stride * ptp->default_i_send_buffer[i+1] * (int) s_data; 
  }
  ptp->async_send_buffer[_request]      = malloc (sizeof (unsigned char) * ptp->async_i_send_buffer[_request][ptp->n_rank]);      

  int delta = (int) s_data * cst_stride;
  for (int i = 0; i < ptp->n_part1; i++) {
    for (int j = 0; j < ptp->selected_part2_idx[i][ptp->n_elt1[i]]; j++) {
      for (int k = ptp->gnum1_to_send_buffer_idx[i][j]; 
               k < ptp->gnum1_to_send_buffer_idx[i][j+1]; 
               k++) {

        if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
          int idx = ptp->gnum1_to_send_buffer[i][k] * delta;
          int idx1 = j* delta;
          for (int k1 = 0; k1 < delta; k1++) {
            ptp->async_send_buffer[_request][idx+k1] = _part1_data[i][idx1+k1]; 
          }
        }
      } 
    }
  }

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    int dest = ptp->active_rank_send[i];
    unsigned char *buf =  ptp->async_send_buffer[_request] + ptp->async_i_send_buffer[_request][dest];
    int count = ptp->async_n_send_buffer[_request][dest];
    PDM_MPI_Issend (buf, count, PDM_MPI_UNSIGNED_CHAR, dest, 
                    tag, ptp->comm, &(ptp->async_send_request[_request][i])); 
  }
}


/**
 *
 * \brief Wait a asynchronus issend
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part1_to_selected_part2_issend_wait
(
 PDM_part1_to_selected_part2_t *ptp,
 int                        request
)
{

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    PDM_MPI_Wait (&(ptp->async_send_request[request][i]));
  }

  _free_async_send_exch (ptp, request);
  
}


/**
 *
 * \brief Initialize a asynchronus irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  s_data        Data size
 * \param [in]  cst_stride    Constant stride
 * \param [in]  part1_data    Partition 2 data
 * \param [in]  tag           Tag of the exchange 
 * \param [out] request       Request
 *
 */

void
PDM_part1_to_selected_part2_irecv
(
 PDM_part1_to_selected_part2_t *ptp,
 const size_t                  s_data,
 const int                     cst_stride,
 void                        **part2_data,
 int                           tag,
 int                          *request
)
{

  *request = _find_open_async_recv_exch (ptp);
  int _request = *request;

  ptp->async_recv_s_data[_request]      = s_data;      
  ptp->async_recv_cst_stride[_request]  = cst_stride;      
  ptp->async_recv_tag[_request]         = tag;
  ptp->async_recv_part2_data[_request]  = part2_data;      
  ptp->async_recv_request[_request]     = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_recv);      
  ptp->async_n_recv_buffer[_request]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_recv_buffer[_request]    = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_recv_buffer[_request][0] = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_recv_buffer[_request][i]   = cst_stride * ptp->default_n_recv_buffer[i] * (int) s_data; 
    ptp->async_i_recv_buffer[_request][i+1] = cst_stride * ptp->default_i_recv_buffer[i+1] * (int) s_data; 
  }
  ptp->async_recv_buffer[_request]      = malloc (sizeof (unsigned char) * ptp->async_i_recv_buffer[_request][ptp->n_rank]);      

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    int source = ptp->active_rank_recv[i];
    unsigned char *buf =  ptp->async_recv_buffer[_request] + ptp->async_i_recv_buffer[_request][source];
    int count = ptp->async_n_recv_buffer[_request][source];
    PDM_MPI_Irecv (buf, count, PDM_MPI_UNSIGNED_CHAR, source, 
                    tag, ptp->comm, &(ptp->async_recv_request[_request][i])); 
  }
  
}


/**
 *
 * \brief Initialize a asynchronus irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part1_to_selected_part2_irecv_wait
(
 PDM_part1_to_selected_part2_t *ptp,
 int                           request
)
{

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    PDM_MPI_Wait (&(ptp->async_recv_request[request][i]));
  }

  size_t s_data  = ptp->async_recv_s_data[request];      
  int cst_stride = ptp->async_recv_cst_stride[request];      

  unsigned char ** _part2_data = (unsigned char **) ptp->async_recv_part2_data[request];

  int delta = (int) s_data * cst_stride;
  for (int i = 0; i < ptp->n_part2; i++) {
    for (int j = 0; j < ptp->n_ref_gnum2[i]; j++) {
      for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
        int idx = ptp->recv_buffer_to_ref_gnum2[i][k] * delta;
        int idx1 = k* delta;
        for (int k1 = 0; k1 < delta; k1++) {
          _part2_data[i][idx1+k1] = ptp->async_recv_buffer[request][idx+k1]; 
        }
      }
    }
  }

  _free_async_recv_exch (ptp, request);

  
}


/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] ptp  Block to part structure
 *
 * \return       NULL
 */

PDM_part1_to_selected_part2_t *
PDM_part1_to_selected_part2_free
(
 PDM_part1_to_selected_part2_t *ptp
)
{
  if (ptp == NULL) {
    return NULL;
  }

  if (ptp->gnum1_to_send_buffer != NULL) {
    for (int i = 0; i < ptp->n_part1; i++) {    
      free (ptp->gnum1_to_send_buffer[i]);
    }
    free (ptp->gnum1_to_send_buffer);
  }

  if (ptp->gnum1_to_send_buffer_idx != NULL) {
    for (int i = 0; i < ptp->n_part1; i++) {    
      free (ptp->gnum1_to_send_buffer_idx[i]);
    }
    free (ptp->gnum1_to_send_buffer_idx);
  }

  if (ptp->recv_buffer_to_ref_gnum2 != NULL) {  
    for (int i = 0; i < ptp->n_part2; i++) {
      if (ptp->ref_gnum2[i] != NULL) {
        free (ptp->ref_gnum2[i]);
      }    
      if (ptp->unref_gnum2[i] != NULL) {
        free (ptp->unref_gnum2[i]);
      }    
      free (ptp->gnum1_come_from_idx[i]);    
      free (ptp->gnum1_come_from[i]);    
      free (ptp->recv_buffer_to_ref_gnum2[i]);
    }
    free (ptp->recv_buffer_to_ref_gnum2);
    free (ptp->ref_gnum2);    
    free (ptp->unref_gnum2);    
    free (ptp->n_ref_gnum2);    
    free (ptp->n_unref_gnum2);    
    free (ptp->gnum1_come_from_idx);    
    free (ptp->gnum1_come_from);    
  }  
  
  free (ptp->active_rank_send);
  free (ptp->active_rank_recv);

  if (ptp->async_send_l_array != 0) {
    for (int i = 0; i < ptp->async_send_l_array; i++) {    
      if (ptp->async_send_buffer[i] != NULL) {
        free (ptp->async_send_buffer[i]);
      }      
      if (ptp->async_n_send_buffer[i] != NULL) {
        free (ptp->async_n_send_buffer[i]);
      }      
      if (ptp->async_i_send_buffer[i] != NULL) {
        free (ptp->async_i_send_buffer[i]);
      }      
      if (ptp->async_send_request[i] != NULL) {
        free (ptp->async_send_request[i]);
      }      
    }
    free (ptp->async_send_open); 
    free (ptp->async_send_s_data); 
    free (ptp->async_send_cst_stride);        
    free (ptp->async_send_tag);        
    free (ptp->async_send_request);    
    free (ptp->async_send_buffer);      
    free (ptp->async_n_send_buffer);  
    free (ptp->async_i_send_buffer);  
  }

  if (ptp->async_recv_l_array != 0) {
    for (int i = 0; i < ptp->async_recv_l_array; i++) {    
      if (ptp->async_recv_buffer[i] != NULL) {
        free (ptp->async_recv_buffer[i]);
      }      
      if (ptp->async_n_recv_buffer[i] != NULL) {
        free (ptp->async_n_recv_buffer[i]);
      }      
      if (ptp->async_i_recv_buffer[i] != NULL) {
        free (ptp->async_i_recv_buffer[i]);
      }      
      if (ptp->async_recv_request[i] != NULL) {
        free (ptp->async_recv_request[i]);
      }      
    }
    free (ptp->async_recv_open); 
    free (ptp->async_recv_s_data); 
    free (ptp->async_recv_cst_stride);        
    free (ptp->async_recv_tag);        
    free (ptp->async_recv_request);    
    free (ptp->async_recv_buffer);      
    free (ptp->async_n_recv_buffer);  
    free (ptp->async_i_recv_buffer);  
    free (ptp->async_recv_part2_data);   
  }

  free (ptp->default_n_send_buffer);
  free (ptp->default_i_send_buffer);
  free (ptp->default_n_recv_buffer);
  free (ptp->default_i_recv_buffer);


  ptp->async_send_n_open  = 0;           
  ptp->async_send_l_array = 0;          

  ptp->async_recv_n_open  = 0;           
  ptp->async_recv_l_array = 0;

  free(ptp);
  return NULL;
  
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
