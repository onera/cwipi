#ifndef __PDM_PARTGNUM1_PARTGNUM2_PRIV_H__
#define __PDM_PARTGNUM1_PARTGNUM2_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/


/**
 * \struct _pdm_partgnum1_partgnum2_t
 *
 * \brief  Data transfer from partitions to blocks
 *
 */

struct _pdm_partgnum1_to_partgnum2_t {

        int           n_part1;                  /*!< Number of parts for gnum1 */    
  const PDM_g_num_t **gnum_elt1;                /*!< gnum of elements in the partition for gnum1 */
  const int          *n_elt1;                   /*!< Number of elements in the partition for gnum1 */

        int           n_part2;                  /*!< Number of parts for gnum2 */
  const PDM_g_num_t **gnum_elt2;                /*!< gnum of elements in the partition for gnum2 */
  const int          *n_elt2;                   /*!< Number of elements in the partition for gnum2 */

  const int         **gnum1_to_gnum2_idx;       /*!< gnum1 to send to gnum2 index */
  const PDM_g_num_t **gnum1_to_gnum2;           /*!< gnum1 to send to gnum2 */
        PDM_MPI_Comm  comm;                     /*!< MPI communicator */  

  int                 n_rank;                   /*!< Number of MPI ranks */
  int                 my_rank;                  /*!< Current rank */

  int                *n_ref_gnum2;              /*!< Numbers of referenced gnum2 (size = \ref n_part2) */
  int               **ref_gnum2;                /*!< Lists of referenced gnum2 (size = \ref n_part2) */

  int                *n_unref_gnum2;            /*!< Numbers of unreferenced gnum2 (size = \ref n_part2) */
  int               **unref_gnum2;              /*!< Lists of unreferenced gnum2 (size = \ref n_part2) */

  int               **gnum1_come_from_idx;      /*!< Index for gnum1_come_from array (size = \ref n_part2) */
  PDM_g_num_t       **gnum1_come_from;          /*!< Gnum come from gnum1 for each referenced gnum2 */

  int               **gnum1_to_send_buffer_idx; /*!< Indirection to store send buffer */
  int               **gnum1_to_send_buffer;     /*!< Indirection to store send buffer */
  int               **recv_buffer_to_ref_gnum2; /*!< Indirection to store gnum2 data from receive buffer */

  int                *default_n_send_buffer;    /*!< Default number of points in the sent buffer */
  int                *default_i_send_buffer;    /*!< Default index in the sent buffer */
  int                *default_n_recv_buffer;    /*!< Default number of points in the received buffer */
  int                *default_i_recv_buffer;    /*!< Default index in the received buffer */

  int                 n_active_rank_send;       /*!< Number of ranks to which data is sent */
  int                *active_rank_send;         /*!< Ranks to which data is sent */

  int                 n_active_rank_recv;       /*!< Number of ranks from which data is received */
  int                *active_rank_recv;         /*!< Ranks from which data is received */

  int                 async_send_n_open;        /*!< Number of open asynchronous exchange */ 
  int                *async_send_open;          /*!< Open asynchronous exchanges */ 
  int                 async_send_l_array;       /*!< Size of arrays to store asynchonous exchanges */ 
  size_t             *async_send_s_data;        /*!< Size of datas of asynchonous exchanges */
  int                *async_send_cst_stride;    /*!< Constant strides of asynchonous exchanges */
  int                *async_send_tag;           /*!< Tag of asynchonous exchanges */
  PDM_MPI_Request   **async_send_request;       /*!< Send requests of asynchonous exchanges */
  unsigned char     **async_send_buffer;        /*!< Send buffers of asynchonous exchanges */
  int               **async_n_send_buffer;      /*!< Number of data in the buffer to send to each rank of asynchonous exchanges */
  int               **async_i_send_buffer;      /*!< Index in the send buffer of each rank of asynchonous exchanges */

  int                 async_recv_n_open;        /*!< Number of open asynchronous exchange */ 
  int                *async_recv_open;          /*!< Open asynchronous exchanges */ 
  int                 async_recv_l_array;       /*!< Size of arrays to store asynchonous exchanges */ 
  size_t             *async_recv_s_data;        /*!< Size of datas of asynchonous exchanges */
  int                *async_recv_cst_stride;    /*!< Constant strides of asynchonous exchanges */
  int                *async_recv_tag;           /*!< Tag of asynchonous exchanges */
  PDM_MPI_Request   **async_recv_request;       /*!< Receive requests of asynchonous exchanges */
  unsigned char     **async_recv_buffer;        /*!< Receive buffers of asynchonous exchanges */
  int               **async_n_recv_buffer;      /*!< Number of data in the buffer received from each rank of asynchonous exchanges */
  int               **async_i_recv_buffer;      /*!< Index in the receive buffer of each rank of asynchonous exchanges */
  void              **async_recv_part2_data;    /*!< Store adress to store data after the receipt of data*/

};


/*=============================================================================
 * Static global variables
 *============================================================================*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PARTGNUM1_PARTGNUM2_PRIV_H__ */
