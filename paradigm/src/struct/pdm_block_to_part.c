/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_block_to_part.h"
#include "pdm_block_to_part_priv.h"
#include "pdm_binary_search.h"
#include "pdm_array.h"
#include "pdm_priv.h"

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

static
void
_comm_graph_statistics
(
 PDM_block_to_part_t* btp
)
{
  /*
   *  Statistic of send --> requested_data_idx
   */
  int min_n_rank_connected  = btp->n_rank+1;
  int max_n_rank_connected  = -1;

  int n_connect_rank = 0;
  for(int i = 0; i < btp->n_rank; ++i) {
    if(btp->requested_data_n[i] > 0) {
      n_connect_rank++;
    }
  }

  double d_n_rank_connected = n_connect_rank;
  double mean_n_rank_connected = 0;
  PDM_MPI_Allreduce(&n_connect_rank    , &max_n_rank_connected , 1, PDM_MPI_INT   , PDM_MPI_MAX, btp->comm);
  PDM_MPI_Allreduce(&n_connect_rank    , &min_n_rank_connected , 1, PDM_MPI_INT   , PDM_MPI_MIN, btp->comm);
  PDM_MPI_Allreduce(&d_n_rank_connected, &mean_n_rank_connected, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, btp->comm);

  mean_n_rank_connected = mean_n_rank_connected/btp->n_rank;

  if(btp->i_rank == 0) {
    printf("PDM_block_to_part requested statistics : [min/max/mean] = %i / %i / %12.5e \n", min_n_rank_connected, max_n_rank_connected, mean_n_rank_connected);
  }

  n_connect_rank = 0;
  for(int i = 0; i < btp->n_rank; ++i) {
    if(btp->distributed_data_n[i] > 0) {
      n_connect_rank++;
    }
  }

  d_n_rank_connected = n_connect_rank;
  mean_n_rank_connected = 0;
  PDM_MPI_Allreduce(&n_connect_rank    , &max_n_rank_connected , 1, PDM_MPI_INT   , PDM_MPI_MAX, btp->comm);
  PDM_MPI_Allreduce(&n_connect_rank    , &min_n_rank_connected , 1, PDM_MPI_INT   , PDM_MPI_MIN, btp->comm);
  PDM_MPI_Allreduce(&d_n_rank_connected, &mean_n_rank_connected, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, btp->comm);

  mean_n_rank_connected = mean_n_rank_connected/btp->n_rank;

  if(btp->i_rank == 0) {
    printf("PDM_block_to_part to send statistics : [min/max/mean] = %i / %i / %12.5e \n", min_n_rank_connected, max_n_rank_connected, mean_n_rank_connected);
  }

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a block to partitions redistribution
 *
 * \param [in]   block_distrib_idx Block distribution (size : \ref size of \ref comm + 1)
 *                               C numbering (block_distrib_idx[0] = 0)
 * \param [in]   gnum_elt        Element global number (size : \ref n_part)
 * \param [in]   n_elt           Local number of elements (size : \ref n_part)
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized \ref PDM_block_to_part instance
 *
 */


PDM_block_to_part_t *
PDM_block_to_part_create_cf
(
 const PDM_g_num_t     *block_distrib_idx,
 const PDM_g_num_t    **gnum_elt,
 const int            *n_elt,
 const int             n_part,
 const PDM_MPI_Fint    fcomm
 )
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(fcomm);
  return PDM_block_to_part_create (block_distrib_idx, gnum_elt, n_elt, n_part, _comm);
}



PDM_block_to_part_t *
PDM_block_to_part_create
(
 const PDM_g_num_t     *block_distrib_idx,
 const PDM_g_num_t    **gnum_elt,
 const int             *n_elt,
 const int              n_part,
 const PDM_MPI_Comm     comm
)
{

  PDM_block_to_part_t *btp =
    (PDM_block_to_part_t *) malloc (sizeof(PDM_block_to_part_t));

  btp->comm = comm;

  btp->pttopt_comm = 0;

  PDM_MPI_Comm_size (comm, &btp->n_rank);
  PDM_MPI_Comm_rank (comm, &btp->i_rank);

  /*
   * Define requested data for each process
   */

  btp->block_distrib_idx = malloc (sizeof(PDM_g_num_t) * (btp->n_rank + 1));
  int max_data_block = -1;
  for (int i = 0; i < btp->n_rank + 1; i++) {
    btp->block_distrib_idx[i] = block_distrib_idx[i];
  }
  for (int i = 0; i < btp->n_rank; i++) {
    max_data_block = PDM_MAX(max_data_block, block_distrib_idx[i+1] - block_distrib_idx[i]) ;
  }

  btp->n_part = n_part;

  btp->requested_data_idx = malloc (sizeof(int) * (btp->n_rank + 1));
  btp->requested_data_n = malloc (sizeof(int) * btp->n_rank);
  for (int i = 0; i < btp->n_rank; i++) {
    btp->requested_data_idx[i] = 0;
    btp->requested_data_n[i] = 0;
  }

  btp->n_elt = malloc (sizeof(int  ) * n_part);
  btp->ind   = malloc (sizeof(int *) * n_part);

  for (int i = 0; i < n_part; i++) {

    btp->n_elt[i] = n_elt[i];
    btp->ind[i] = malloc (sizeof(int) * n_elt[i]);

    const PDM_g_num_t *_gnum_elt = gnum_elt[i];

    for (int j = 0; j < n_elt[i]; j++) {

      int ind = PDM_binary_search_gap_long (_gnum_elt[j] - 1,
                                            block_distrib_idx,
                                            btp->n_rank + 1);
      // printf(" [%i][%i] --> ind = %i (g_num = %i )\n", i, j, ind, (int) _gnum_elt[j]);
      btp->requested_data_n[ind]++;

    }
  }

  for (int i = 0; i < btp->n_rank; i++) {
    btp->requested_data_idx[i+1] = btp->requested_data_idx[i] + btp->requested_data_n  [i];
  }

  int s_requested_data = btp->requested_data_idx[btp->n_rank - 1]
                       + btp->requested_data_n  [btp->n_rank - 1];

  int *requested_data = malloc (sizeof(int) *  s_requested_data);

  for (int i = 0; i < btp->n_rank; i++) {
    btp->requested_data_n[i] = 0;
  }

  for (int i = 0; i < n_part; i++) {

    const PDM_g_num_t *_gnum_elt = gnum_elt[i];

    // printf("n_elt[%i] = %i \n", i, (int) n_elt[i]);
    for (int j = 0; j < n_elt[i]; j++) {

      int ind = PDM_binary_search_gap_long (_gnum_elt[j] - 1,
                                            block_distrib_idx,
                                            btp->n_rank + 1);
      int idx = btp->requested_data_idx[ind] + btp->requested_data_n[ind]++;

      btp->ind[i][j] = idx;

      PDM_g_num_t _requested_data = _gnum_elt[j] - 1 - block_distrib_idx[ind];
      // printf("requested_data[%i] = %i / size_max = %i and gn_m = %i \n", idx, (int) _requested_data, s_requested_data, (int)_gnum_elt[j]);
      requested_data[idx] = (int) _requested_data;
    }
  }

  btp->distributed_data_n = malloc (sizeof(int) * btp->n_rank);

  PDM_MPI_Alltoall (btp->requested_data_n,   1, PDM_MPI_INT,
                    btp->distributed_data_n, 1, PDM_MPI_INT,
                    comm);

  btp->distributed_data_idx = PDM_array_new_idx_from_sizes_int(btp->distributed_data_n, btp->n_rank);

  btp->distributed_data = malloc (sizeof(int) *
                                  btp->distributed_data_idx[btp->n_rank]);

  PDM_MPI_Alltoallv (requested_data,
                     btp->requested_data_n,
                     btp->requested_data_idx,
                     PDM_MPI_INT,
                     btp->distributed_data,
                     btp->distributed_data_n,
                     btp->distributed_data_idx,
                     PDM_MPI_INT,
                     comm);

  int coeff = 10;
  if (btp->distributed_data_idx[btp->n_rank] >= coeff * max_data_block) {
    btp->pttopt_comm = 1;
  }

  if(0 == 1) {
    _comm_graph_statistics(btp);
  }

  //PDM_log_trace_array_long(btp->distributed_data_idx, btp->n_rank+1, "block_distrib");

  int tmp;
  PDM_MPI_Allreduce (&(btp->pttopt_comm), &tmp, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  btp->pttopt_comm = tmp;

  free (requested_data);

  return (PDM_block_to_part_t *) btp;

}


/**
 *
 * \brief Initialize an exchange
 *
 * \param [in]   btp          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
 *                            Constant stride for \ref PDM_STRIDE_VAR
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Partition stride or NULL
 * \param [out]  part_data    Partition data
 *
 */

void
PDM_block_to_part_exch
(
 PDM_block_to_part_t *btp,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                 *block_stride,
 void                *block_data,
 int                **part_stride,
 void               **part_data
)
{
  unsigned char *_block_data = (unsigned char *) block_data;
  unsigned char **_part_data = (unsigned char **) part_data;

  int n_elt_block = btp->block_distrib_idx[btp->i_rank+1] - btp->block_distrib_idx[btp->i_rank];

  size_t *i_send_buffer = (size_t *) malloc (sizeof(size_t) * btp->n_rank);
  size_t *i_recv_buffer = (size_t *) malloc (sizeof(size_t) * btp->n_rank);
  int *n_send_buffer = (int *) malloc (sizeof(int) * btp->n_rank);
  int *n_recv_buffer = (int *) malloc (sizeof(int) * btp->n_rank);
  int max_n_send_buffer = -1;
  int max_n_recv_buffer = -1;
  int *block_stride_idx = NULL;

  for (int i = 0; i < btp->n_rank; i++) {
    n_send_buffer[i] = 0;
    n_recv_buffer[i] = 0;
    i_send_buffer[i] = 0;
    i_recv_buffer[i] = 0;
  }

  unsigned char **send_buffer = NULL;
  unsigned char *recv_buffer = NULL;

  size_t s_send_buffer = 0;
  size_t s_recv_buffer = 0;

  int n_rank1 = btp->n_rank - 1;

  int s_distributed_data = btp->distributed_data_idx[btp->n_rank];

  /* int step; */

  int rank;
  PDM_MPI_Comm_rank(btp->comm, &rank);

  /*
   * Exchange Stride and build buffer properties
   */

  int *recv_stride = NULL;
  if (t_stride == PDM_STRIDE_VAR) {

    int s_send_stride = btp->distributed_data_idx[btp->n_rank];

    int s_recv_stride = btp->requested_data_idx[btp->n_rank];

    int *send_stride = (int *) malloc (sizeof(int) * s_send_stride);
    recv_stride = (int *) malloc (sizeof(int) * s_recv_stride);

    for (int i = 0; i < s_send_stride; i++) {
      send_stride[i] = block_stride[btp->distributed_data[i]];
    }

    PDM_MPI_Alltoallv (send_stride,
                       btp->distributed_data_n,
                       btp->distributed_data_idx,
                       PDM_MPI_INT,
                       recv_stride,
                       btp->requested_data_n,
                       btp->requested_data_idx,
                       PDM_MPI_INT,
                       btp->comm);

    for (int i = 0; i < btp->n_part; i++) {
      for (int j = 0; j < btp->n_elt[i]; j++) {
        int ielt = btp->ind[i][j];
        part_stride[i][j] = recv_stride[ielt];
      }
    }

    /*
     * Build buffers
     */

    for (int i = 0; i < btp->n_rank; i++) {
      int ibeg = btp->distributed_data_idx[i];
      int iend = btp->distributed_data_idx[i] +
                 btp->distributed_data_n[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)  {
        n_send_buffer[i] += send_stride[k];
      }

      n_send_buffer[i] *= (int) s_data;
      max_n_send_buffer = PDM_MAX(max_n_send_buffer, n_send_buffer[i]);

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      } else {
        i_send_buffer[i] = 0;
      }

      ibeg = btp->requested_data_idx[i];
      iend = btp->requested_data_idx[i] +
        btp->requested_data_n[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++) {
        n_recv_buffer[i] += recv_stride[k];
      }

      n_recv_buffer[i] *= (int) s_data;
      max_n_recv_buffer = PDM_MAX(max_n_recv_buffer, n_recv_buffer[i]);

      if (i > 0) {
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      } else {
        i_recv_buffer[i] = 0;
      }
    }

    block_stride_idx = PDM_array_new_idx_from_sizes_int(block_stride, n_elt_block);
    free(send_stride);
  }

  else {

    int cst_stride = *block_stride;
    max_n_send_buffer = 0;
    max_n_recv_buffer = 0;

    for (int i = 0; i < btp->n_rank; i++) {

      i_send_buffer[i] = btp->distributed_data_idx[i] * cst_stride * (int) s_data;
      i_recv_buffer[i] = btp->requested_data_idx[i] * cst_stride * (int) s_data;

      n_send_buffer[i] = btp->distributed_data_n[i] * cst_stride * (int) s_data;
      n_recv_buffer[i] = btp->requested_data_n[i] * cst_stride * (int) s_data;
      max_n_send_buffer = PDM_MAX(max_n_send_buffer, n_send_buffer[i]);
      max_n_recv_buffer = PDM_MAX(max_n_recv_buffer, n_recv_buffer[i]);

    }

    s_send_buffer = i_send_buffer[n_rank1] + n_send_buffer[n_rank1];
    s_recv_buffer = i_recv_buffer[n_rank1] + n_recv_buffer[n_rank1];

  }

  s_recv_buffer = i_recv_buffer[n_rank1] + n_recv_buffer[n_rank1];

  int n_active_buffer;

  if (btp->pttopt_comm) {
    n_active_buffer = 5;
  }
  else {
    n_active_buffer = 1;
  }

  send_buffer = (unsigned char **) malloc(sizeof(unsigned char *) * n_active_buffer);

  if (btp->pttopt_comm) {
    for (int i = 0; i < n_active_buffer; i++) {
      send_buffer[i] = (unsigned char *) malloc(sizeof(unsigned char) *  max_n_send_buffer);
    }
  }
  else {
    s_send_buffer = i_send_buffer[n_rank1] + n_send_buffer[n_rank1];
    send_buffer[0] = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
  }

  recv_buffer = (unsigned char *) malloc(sizeof(unsigned char ) * s_recv_buffer);

  if (btp->pttopt_comm) {

    PDM_MPI_Request *s_request =  malloc (sizeof(PDM_MPI_Request) * n_active_buffer);
    PDM_MPI_Request *r_request = malloc (sizeof(PDM_MPI_Request) * btp->n_rank);

    for (int i = 0; i < btp->n_rank; i++) {
      if (n_recv_buffer[i] > 0) {
        PDM_MPI_Irecv(recv_buffer + i_recv_buffer[i],
                      n_recv_buffer[i],
                      PDM_MPI_BYTE,
                      i,
                      0,
                      btp->comm,
                      r_request + i);
      }
    }

    int *active_rank = malloc(sizeof(int) * n_active_buffer);
    for (int i = 0; i < n_active_buffer; i++) {
      active_rank[i] = i;
    }

    while (1) {
      int _n_active_buffer = 0;
      for (int i = 0; i < n_active_buffer; i++) {
        if (active_rank[i] < btp->n_rank) {
          _n_active_buffer += 1;
        }
      }

      if (_n_active_buffer == 0) {
        break;
      }

      for (int i = 0; i < _n_active_buffer; i++) {
        if (n_send_buffer[active_rank[i]] > 0) {

          int s_distributed_active_rank = btp->distributed_data_idx[active_rank[i]] +
                                          btp->distributed_data_n [active_rank[i]];

          if (t_stride == PDM_STRIDE_VAR) {
            int idx1 = 0;
            for (int j = btp->distributed_data_idx[active_rank[i]];
                 j < s_distributed_active_rank; j++) {

              int ind =  block_stride_idx[btp->distributed_data[j]] * (int) s_data;

              int s_block_unit =  block_stride[btp->distributed_data[j]] * (int) s_data;

              unsigned char *_block_data_deb = _block_data + ind;

              for (int k = 0; k < s_block_unit; k++) {
                send_buffer[i][idx1++] = _block_data_deb[k];
              }
            }
          }
          else {
            int cst_stride = *block_stride;
            int s_block_unit = cst_stride * (int) s_data;

            int idx1 = 0;
            for (int j = btp->distributed_data_idx[active_rank[i]];
                 j < s_distributed_active_rank; j++) {
              int ind = btp->distributed_data[j];
              unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;
              for (int k = 0; k < s_block_unit; k++) {
                send_buffer[i][idx1++] = _block_data_deb[k];
              }
            }
          }

          PDM_MPI_Issend(send_buffer[i],
                         n_send_buffer[active_rank[i]],
                         PDM_MPI_BYTE,
                         active_rank[i],
                         0,
                         btp->comm,
                         s_request + i);
        }
      }

      for (int i = 0; i < _n_active_buffer; i++) {
        if (n_send_buffer[active_rank[i]] > 0) {
          PDM_MPI_Wait (s_request + i);
        }
      }

      for (int i = 0; i < n_active_buffer; i++) {
        active_rank[i] += n_active_buffer;
      }

    }

    for (int i = 0; i < btp->n_rank; i++) {
      if (n_recv_buffer[i] > 0) {
        PDM_MPI_Wait (r_request + i);
      }
    }

    free (s_request);
    free (r_request);
    free (active_rank);

  }

  else {

    if (t_stride == PDM_STRIDE_VAR) {
      int idx1 = 0;
      for (int i = 0; i < s_distributed_data; i++) {
        int ind =  block_stride_idx[btp->distributed_data[i]] * (int) s_data;
        int s_block_unit =  block_stride[btp->distributed_data[i]] * (int) s_data;
        unsigned char *_block_data_deb = _block_data + ind;
        for (int k = 0; k < s_block_unit; k++) {
          send_buffer[0][idx1++] = _block_data_deb[k];
        }
      }
    }
    else {
      int idx1 = 0;
      int cst_stride = *block_stride;
      int s_block_unit = cst_stride * (int) s_data;
      for (int i = 0; i < s_distributed_data; i++) {
        int ind = btp->distributed_data[i];
        unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;
        for (int k = 0; k < s_block_unit; k++) {
          send_buffer[0][idx1++] = _block_data_deb[k];
        }
      }
    }

    PDM_MPI_Alltoallv_l(send_buffer[0],
                        n_send_buffer,
                        i_send_buffer,
                        PDM_MPI_BYTE,
                        recv_buffer,
                        n_recv_buffer,
                        i_recv_buffer,
                        PDM_MPI_BYTE,
                        btp->comm);
  }

  for (int i = 0; i < n_active_buffer; i++) {
    free(send_buffer[i]);
  }
  free(send_buffer);
  free(n_send_buffer);
  free(i_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  if (block_stride_idx != NULL) {
    free (block_stride_idx);
  }

  /*
   * Partitions filling
   */

  if (t_stride == PDM_STRIDE_VAR) {

    int s_recv_elt = btp->requested_data_idx[n_rank1] +
      btp->requested_data_n[n_rank1];

    int **part_idx = malloc (sizeof(int *) * btp->n_part);
    int  *recv_idx = PDM_array_new_idx_from_sizes_int(recv_stride, s_recv_elt);

    for (int i = 0; i < btp->n_part; i++) {
      part_idx[i] = PDM_array_new_idx_from_sizes_int(part_stride[i], btp->n_elt[i]);
    }

    for (int i = 0; i < btp->n_part; i++) {

      for (int j = 0; j < btp->n_elt[i]; j++) {

        int idx1  = part_idx[i][j] * (int) s_data;
        int n_elt = part_stride[i][j] * (int) s_data;

        int idx2 = recv_idx[btp->ind[i][j]] * (int) s_data;

        for (int k = 0; k < n_elt; k++) {
          _part_data[i][idx1+k] = recv_buffer[idx2+k];
        }
      }
    }

    for (int i = 0; i < btp->n_part; i++) {
      free (part_idx[i]);
    }

    free(recv_idx);
    free(part_idx);
    free (recv_stride);
  }

  else if (t_stride == PDM_STRIDE_CST) {

    const int cst_stride = *block_stride;
    const int s_block_unit = cst_stride * (int) s_data;

    for (int i = 0; i < btp->n_part; i++) {

      for (int j = 0; j < btp->n_elt[i]; j++) {

        int idx1  = j * s_block_unit;
        int idx2 = btp->ind[i][j] * s_block_unit;

        for (int k = 0; k < s_block_unit; k++) {
          _part_data[i][idx1+k] = recv_buffer[idx2+k];
        }
      }
    }
  }

  free(recv_buffer);

}



/**
 *
 * \brief Initialize an exchange
 * (part_stride and part_data are allocated in function)
 *
 * \param [in]   btp          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
 *                            Constant stride for \ref PDM_STRIDE_VAR
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Partition stride or NULL
 * \param [out]  part_data    Partition data
 *
 */

void
PDM_block_to_part_exch2
(
 PDM_block_to_part_t *btp,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                 *block_stride,
 void                *block_data,
 int               ***part_stride,
 void              ***part_data
)
{

  int n_elt_block = btp->block_distrib_idx[btp->i_rank+1] - btp->block_distrib_idx[btp->i_rank];

  unsigned char *_block_data = (unsigned char *) block_data;
  unsigned char **_part_data;

  size_t *i_send_buffer = (size_t *) malloc (sizeof(size_t) * btp->n_rank);
  size_t *i_recv_buffer = (size_t *) malloc (sizeof(size_t) * btp->n_rank);
  int *n_send_buffer = (int *) malloc (sizeof(int) * btp->n_rank);
  int *n_recv_buffer = (int *) malloc (sizeof(int) * btp->n_rank);

  for (int i = 0; i < btp->n_rank; i++) {
    n_send_buffer[i] = 0;
    n_recv_buffer[i] = 0;
    i_send_buffer[i] = 0;
    i_recv_buffer[i] = 0;
  }

  unsigned char *send_buffer = NULL;
  unsigned char *recv_buffer = NULL;

  size_t s_send_buffer = 0;
  size_t s_recv_buffer = 0;

  int n_rank1 = btp->n_rank - 1;

  int s_distributed_data = btp->distributed_data_idx[btp->n_rank];

  /*
   * Exchange Stride and build buffer properties
   */

  int *recv_stride = NULL;
  int **_part_stride = NULL;

  if (t_stride == PDM_STRIDE_VAR) {

    int s_send_stride = btp->distributed_data_idx[btp->n_rank];

    int s_recv_stride = btp->requested_data_idx[btp->n_rank];

    int *send_stride = (int *) malloc (sizeof(int) * s_send_stride);
    recv_stride = (int *) malloc (sizeof(int) * s_recv_stride);

    for (int i = 0; i < s_send_stride; i++) {
      send_stride[i] = block_stride[btp->distributed_data[i]];
    }

    PDM_MPI_Alltoallv (send_stride,
                       btp->distributed_data_n,
                       btp->distributed_data_idx,
                       PDM_MPI_INT,
                       recv_stride,
                       btp->requested_data_n,
                       btp->requested_data_idx,
                       PDM_MPI_INT,
                       btp->comm);

    *part_stride = (int **) malloc(sizeof(int *) * btp->n_part);
    _part_stride = *part_stride;

    for (int i = 0; i < btp->n_part; i++) {

      _part_stride[i] = malloc (sizeof(int) * btp->n_elt[i]);

      for (int j = 0; j < btp->n_elt[i]; j++) {

        int ielt = btp->ind[i][j];
        _part_stride[i][j] = recv_stride[ielt];

      }
    }

    /*
     * Build buffers
     */

    for (int i = 0; i < btp->n_rank; i++) {
      int ibeg = btp->distributed_data_idx[i];
      int iend = btp->distributed_data_idx[i] +
                 btp->distributed_data_n[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)  {
        n_send_buffer[i] += send_stride[k];
      }

      n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      }
      else {
        i_send_buffer[i] = 0;
      }

      ibeg = btp->requested_data_idx[i];
      iend = btp->requested_data_idx[i] +
             btp->requested_data_n[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++) {
        n_recv_buffer[i] += recv_stride[k];
      }

      n_recv_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      }
      else {
        i_recv_buffer[i] = 0;
      }

    }

    s_send_buffer = i_send_buffer[n_rank1] + n_send_buffer[n_rank1];
    s_recv_buffer = i_recv_buffer[n_rank1] + n_recv_buffer[n_rank1];

    send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

    // int *send_stride_idx = (int *) malloc(sizeof(int) * (s_distributed_data+1));
    // send_stride_idx[0] = 0;
    // for (int i = 0; i < s_distributed_data; i++) {
    //   send_stride_idx[i+1] = send_stride_idx[i] + send_stride[i];
    // }

    int idx1 = 0;

    int *block_stride_idx = PDM_array_new_idx_from_sizes_int(block_stride, n_elt_block);

    for (int i = 0; i < s_distributed_data; i++) {

      int ind =  block_stride_idx[btp->distributed_data[i]] * (int) s_data;

      int s_block_unit =  block_stride[btp->distributed_data[i]] * (int) s_data;

      unsigned char *_block_data_deb = _block_data + ind;
      for (int k = 0; k < s_block_unit; k++) {
        send_buffer[idx1++] = _block_data_deb[k];
      }
    }
    free(send_stride);
    // free(send_stride_idx);
    free(block_stride_idx);

  }

  else if (t_stride == PDM_STRIDE_CST) {

    int cst_stride = *block_stride;
    int s_block_unit = cst_stride * (int) s_data;

    for (int i = 0; i < btp->n_rank; i++) {

      i_send_buffer[i] = btp->distributed_data_idx[i] * cst_stride * (int) s_data;
      i_recv_buffer[i] = btp->requested_data_idx[i] * cst_stride * (int) s_data;

      n_send_buffer[i] = btp->distributed_data_n[i] * cst_stride * (int) s_data;
      n_recv_buffer[i] = btp->requested_data_n[i] * cst_stride * (int) s_data;

    }

    s_send_buffer = i_send_buffer[n_rank1] + n_send_buffer[n_rank1];
    s_recv_buffer = i_recv_buffer[n_rank1] + n_recv_buffer[n_rank1];

    send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

    int idx1 = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int ind = btp->distributed_data[i];
      unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;
      for (int k = 0; k < s_block_unit; k++) {
        send_buffer[idx1++] = _block_data_deb[k];
      }
    }
  }

  /*
   * Data exchange
   */

  PDM_MPI_Alltoallv_l(send_buffer,
                      n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_BYTE,
                      recv_buffer,
                      n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_BYTE,
                      btp->comm);

  free(send_buffer);
  free(n_send_buffer);
  free(i_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  /*
   * Partitions filling
   */

  *part_data = malloc(sizeof(unsigned char *) * btp->n_part);
  _part_data = (*(unsigned char ***) part_data);

  if (t_stride == PDM_STRIDE_VAR) {

    int s_recv_elt = btp->requested_data_idx[n_rank1] +
                     btp->requested_data_n[n_rank1];

    int **part_idx = malloc (sizeof(int *) * btp->n_part);
    int *recv_idx = PDM_array_new_idx_from_sizes_int(recv_stride, s_recv_elt);

    for (int i = 0; i < btp->n_part; i++) {
      part_idx[i] = PDM_array_new_idx_from_sizes_int(_part_stride[i], btp->n_elt[i]);
    }

    for (int i = 0; i < btp->n_part; i++) {

      int s_part =  part_idx[i][btp->n_elt[i]] * (int) s_data;

      _part_data[i] = malloc(sizeof(unsigned char) * s_part);

      for (int j = 0; j < btp->n_elt[i]; j++) {

        int idx1  = part_idx[i][j] * (int) s_data;
        int n_elt = _part_stride[i][j] * (int) s_data;

        int idx2 = recv_idx[btp->ind[i][j]] * (int) s_data;

        for (int k = 0; k < n_elt; k++) {
           _part_data[i][idx1+k] = recv_buffer[idx2+k];
        }
      }
    }

    for (int i = 0; i < btp->n_part; i++) {
      free (part_idx[i]);
    }

    free(recv_idx);
    free(part_idx);
    free (recv_stride);
  }

  else if (t_stride == PDM_STRIDE_CST) {

    const int cst_stride = *block_stride;
    const int s_block_unit = cst_stride * (int) s_data;

    for (int i = 0; i < btp->n_part; i++) {

      _part_data[i] = malloc(sizeof(unsigned char) * s_block_unit * btp->n_elt[i]);

      for (int j = 0; j < btp->n_elt[i]; j++) {

        int idx1  = j * s_block_unit;
        int idx2 = btp->ind[i][j] * s_block_unit;

        for (int k = 0; k < s_block_unit; k++) {
           _part_data[i][idx1+k] = recv_buffer[idx2+k];
        }
      }
    }
  }

  free(recv_buffer);

}


/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] btp  Block to part structure
 *
 * \return       NULL
 */

PDM_block_to_part_t *
PDM_block_to_part_free
(
 PDM_block_to_part_t *btp
)
{

  for (int i = 0; i < btp->n_part; i++) {
    free (btp->ind[i]);
  }

  free (btp->ind);
  free (btp->n_elt);
  free (btp->block_distrib_idx);
  free (btp->distributed_data);
  free (btp->distributed_data_idx);
  free (btp->distributed_data_n);
  free (btp->requested_data_idx);
  free (btp->requested_data_n);

  free (btp);

  return NULL;
}


/**
 *
 * \brief Return index in the block for a gnum
 *
 * \param [in] ptb         Part to block structure
 * \param [in] gNum        Global number
 *
 * \return  Index
 */

PDM_l_num_t
PDM_block_to_part_gnum_idx_get
(
 PDM_block_to_part_t *btp,
 PDM_g_num_t gNum
)
{
  return (PDM_l_num_t) (gNum - 1 - btp->block_distrib_idx[btp->i_rank]);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
