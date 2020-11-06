/*
 * File:   pdm_block_to_part.h
 * Author: equemera
 *
 * Created on April 14, 2016, 7:56 AM
 */

#ifndef PDM_BLOCK_TO_PART_H
#define	PDM_BLOCK_TO_PART_H

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

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct PDM_block_to_part_t
 * \brief  Block to partition redistribution
 *
 */

typedef struct _pdm_block_to_part_t PDM_block_to_part_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a block to partitions redistribution
 *
 * \param [in]   block_distrib_idx Block distribution (size : \ref size of \ref comm + 1)
 * \param [in]   gnum_elt          Element global number (size : \ref n_part)
 * \param [in]   n_elt             Local number of elements (size : \ref n_part)
 * \param [in]   n_part            Number of partition
 * \param [in]   comm              MPI communicator
 *
 * \return   Initialized \ref PDM_block_to_part instance
 *
 */

PDM_block_to_part_t *
PDM_block_to_part_create
(
 const PDM_g_num_t    *block_distrib_idx,
 const PDM_g_num_t    **gnum_elt,
 const int            *n_elt,
 const int             n_part,
 const PDM_MPI_Comm        comm
);


PDM_block_to_part_t *
PDM_block_to_part_create_cf
(
 const PDM_g_num_t    *block_distrib_idx,
 const PDM_g_num_t    **gnum_elt,
 const int            *n_elt,
 const int             n_part,
 const PDM_MPI_Fint    fcomm
);


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
 * \param [out]  part_stride  Partition stride
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
);


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
 * \param [out]  part_stride  Partition stride
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
);


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
);


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
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_BLOCK_TO_PART_H */
