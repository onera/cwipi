#ifndef __PDM_PART1_TO_SELECTED_PART2_H__
#define	__PDM_PART1_TO_SELECTED_PART2_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct PDM_part1_to_selected_part2_t
 * \brief  Block to partition redistribution
 *
 */

typedef struct _pdm_part1_to_selected_part2_t PDM_part1_to_selected_part2_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
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
 * \param [in]   gnum1_to_gnum2_idx  Index of data to send to gnum2 from gnum1 
 *                                  (for each part size : \ref n_elt1+1) 
 * \param [in]   gnum1_to_gnum2      Data to send to gnum2 from gnum1 
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
 const int           **gnum1_to_gnum2_idx,
 const PDM_g_num_t   **gnum1_to_gnum2,
 const PDM_MPI_Comm    comm
);


PDM_part1_to_selected_part2_t *
PDM_part1_to_selected_part2_create_cf
(
 const PDM_g_num_t    **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t    **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **gnum1_to_gnum2_idx,
 const PDM_g_num_t   **gnum1_to_gnum2,
 const PDM_MPI_Fint    fcomm
);


/**
 *
 * \brief Initialize an exchange
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   s_data        Data size
 * \param [in]   t_stride      Stride type
 * \param [in]   part1_stride  Partition 1 stride
 * \param [in]   part1_data    Partition 1 data
 * \param [out]  part2_stride  Partition 2 stride
 * \param [out]  part2_data    Partition 2 data
 *
 */

void
PDM_part1_to_selected_part2_exch
(
 PDM_part1_to_selected_part2_t    *ptp,
 const size_t                  s_data,
 const PDM_stride_t            t_stride,
 int                         **part1_stride,
 void                        **part1_data,
 int                         **part2_stride,
 void                        **part2_data
);


/**
 *
 * \brief Initialize an exchange with allocation of result arrays
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   s_data        Data size
 * \param [in]   t_stride      Stride type
 * \param [in]   part1_stride  Partition 1 stride
 * \param [in]   part1_data    Partition 1 data
 * \param [out]  part2_stride  Partition 2 stride
 * \param [out]  part2_data    Partition 2 data
 *
 */

void
PDM_part1_to_selected_part2_exch_with_alloc
(
 PDM_part1_to_selected_part2_t *ptp,
 const size_t               s_data,
 const PDM_stride_t         t_stride,
 int                      **part1_stride,
 void                     **part1_data,
 int                     ***part2_stride,
 void                    ***part2_data
);


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
);


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
);


/**
 *
 * \brief Get gnum come from gnum1 for each referenced gnum2
 *
 * \param [in]   ptp                 Block to part structure
 * \param [out]  gnum1_come_from_idx Index for gnum1_come_from array (size = \ref n_part2)  
 * \param [out]  gnum1_come_from     gnum come from gnum1 for each referenced gnum2
 *
 */

void
PDM_part1_to_selected_part2_gnum1_come_from_get
(
 PDM_part1_to_selected_part2_t *ptp,
 int                        ***gnum1_come_from_idx,
 PDM_g_num_t                ***gnum1_come_from
);


/**
 *
 * \brief Initialize a asynchronus issend
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   s_data        Data size
 * \param [in]   cst_stride    Constant stride
 * \param [in]   part1_data    Partition 1 data
 * \param [in]   tag           Tag of the exchange 
 * \param [out]  request       Request
 *
 */

void
PDM_part1_to_selected_part2_issend
(
 PDM_part1_to_selected_part2_t *ptp,
 const size_t                  s_data,
 const int                     cst_stride,
 void                        **part1_data,
 int                           tag,
 int                          *request
);


/**
 *
 * \brief Wait a asynchronus issend
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  tag           Tag of the exchange 
 * \param [in]  request       Request
 *
 */

void
PDM_part1_to_selected_part2_issend_wait
(
 PDM_part1_to_selected_part2_t *ptp,
 int                           request
);


/**
 *
 * \brief Initialize a asynchronus irecv from reference gnum2
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
);


/**
 *
 * \brief Initialize a asynchronus irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  tag           Tag of the exchange 
 * \param [in]  request       Request
 *
 */

void
PDM_part1_to_selected_part2_irecv_wait
(
 PDM_part1_to_selected_part2_t *ptp,
 int                           request
);


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
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_tO_pART_H */
