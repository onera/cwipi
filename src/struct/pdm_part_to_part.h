#ifndef __PDM_PART_TO_PART_H__
#define	__PDM_PART_TO_PART_H__

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
 * \struct PDM_part_to_part_t
 * \brief  Block to partition redistribution
 *
 */

typedef struct _pdm_part_to_part_t PDM_part_to_part_t;


/**
 * \enum PDM_part_to_part_data_def_t
 * \brief Kind of data definition
 *
 */

typedef enum {

  PDM_PART_TO_PART_DATA_DEF_ORDER_PART1           = 0, /*!< Data defined according 
                                                         to the part1 arrays order*/
  PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2  = 1, /*!< Data defined according 
                                                         to the part1_to_part2 arrays order */
  PDM_PART_TO_PART_DATA_DEF_ORDER_PART2           = 2, /*!< Data defined according 
                                                         to the part2 arrays order*/
  PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM = 3, /*!< Data defined according 
                                                         to the gnum1_come_from arrays order */

} PDM_part_to_part_data_def_t;

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
 * \param [in]   part1_to_part2_idx Index of data to send to gnum2 from gnum1 
 *                                  (for each part size : \ref n_elt1+1) 
 * \param [in]   part1_to_part2     Data to send to gnum2 from gnum1 
 * \param [in]   comm               MPI communicator
 *
 * \return   Initialized \ref PDM_part_to_part instance
 *
 */

PDM_part_to_part_t *
PDM_part_to_part_create
(
 const PDM_g_num_t   **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t   **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **part1_to_part2_idx,
 const PDM_g_num_t   **part1_to_part2,
 const PDM_MPI_Comm    comm
);


PDM_part_to_part_t *
PDM_part_to_part_create_cf
(
 const PDM_g_num_t    **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t    **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **part1_to_part2_idx,
 const PDM_g_num_t   **part1_to_part2,
 const PDM_MPI_Fint    fcomm
);


/**
 *
 * \brief Initialize an exchange based on MPI_ialltoall
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Data in same order than part1_to_part2 array
 * \param [out]  ref_part2_data      Data to referenced part2 elements
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_ialltoall
(
PDM_part_to_part_t *ptp,
 const size_t       s_data,
 const int          cst_stride,
 void             **part1_to_part2_data,
 void             **ref_part2_data,
 int               *request
);


/**
 *
 * \brief Initialize an exchange based on MPI_ialltoall
 *
 * \param [in]   ptp                     Block to part structure
 * \param [in]   s_data                  Data size
 * \param [in]   cst_stride              Constant stride
 * \param [in]   ref_part2_to_part1_data Data in ref_gnum2 order
 * \param [out]  part1_part2_data        Data in part1_to_part2 order
 * \param [out]  request                 Request
 *
 */

void
PDM_part_to_part_reserve_ialltoall
(
PDM_part_to_part_t *ptp,
 const size_t       s_data,
 const int          cst_stride,
 void             **ref_part2_to_part1_data,
 void             **part1_part2_data,
 int               *request
);


/**
 *
 * \brief Wait a asynchronus alltoall wait
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_ialltoall_wait
(
 PDM_part_to_part_t *ptp,
 int                 request
);


/**
 *
 * \brief Wait a asynchronus reverse alltoall wait
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_reverse_ialltoall_wait
(
 PDM_part_to_part_t *ptp,
 int                 request
);


/**
 *
 * \brief Initialize an exchange based on MPI_ineighbor_alltoall
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Data in same order than part1_to_part2 array
 * \param [out]  ref_part2_data      Data to referenced part2 elements
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_ineighbor_alltoall
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void              **part1_to_part2_data,
 void              **ref_part2_data,
 int                *request
);


/**
 *
 * \brief Wait a asynchronus issend
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_ineighbor_alltoall_wait
(
 PDM_part_to_part_t *ptp,
 int                 request
);


/**
 *
 * \brief Get selected numbers of part2
 *
 * \param [in]   ptp                 Block to part structure
 * \param [out]  n_elt1              Number of gnum1 element  
 * \param [out]  part1_to_part2_idx  Index of data to send to gnum2 from gnum1 
 *                                  (for each part size : \ref n_elt1+1) 
 * \param [out]  part1_to_part2      Data to send to gnum2 from gnum1 for each part
 *
 */

void
PDM_part_to_part_part1_to_part2_get
(
 PDM_part_to_part_t *ptp,
 int               **n_elt1,
 int              ***part1_to_part2_idx,
 PDM_g_num_t      ***part1_to_part2
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
PDM_part_to_part_ref_gnum2_get
(
 PDM_part_to_part_t *ptp,
 int               **n_ref_gnum2,
 int              ***ref_gnum2
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
PDM_part_to_part_unref_gnum2_get
(
 PDM_part_to_part_t *ptp,
 int               **n_unref_gnum2,
 int               ***unref_gnum2
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
PDM_part_to_part_gnum1_come_from_get
(
 PDM_part_to_part_t *ptp,
 int              ***gnum1_come_from_idx,
 PDM_g_num_t      ***gnum1_come_from
);


/**
 *
 * \brief Initialize an asynchronus issend
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Data (order given by part1_to_part2 array)
 * \param [in]   tag                 Tag of the exchange 
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_issend
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 const void        **part1_to_part2_data,
 const int           tag,
 int                *request
);


/**
 *
 * \brief Wait an asynchronus issend (part1 to part2)
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  tag           Tag of the exchange 
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_issend_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
);


/**
 *
 * \brief Initialize an asynchronus reverse issend (part2 to part1)
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part2_to_part1_data Data (order given by gnum1_come_from and ref_gnum2 arrays)  
 * \param [in]   tag                 Tag of the exchange 
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_reverse_issend
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 const void        **part2_to_part1_data,
 const int           tag,
 int                *request
);


/**
 *
 * \brief Wait an asynchronus reverse issend (part2 to part1)
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  tag           Tag of the exchange 
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_reverse_issend_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
);


/**
 *
 * \brief Initialize a asynchronus irecv (from part1)
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  s_data        Data size
 * \param [in]  cst_stride    Constant stride
 * \param [in]  part2_data    Partition 2 data (order given by gnum1_come_from and ref_gnum2 arrays) 
 * \param [in]  tag           Tag of the exchange 
 * \param [out] request       Request
 *
 */

void
PDM_part_to_part_irecv
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void              **part2_data,
 int                 tag,
 int                *request
);


/**
 *
 * \brief Initialize a asynchronus irecv (from part1)
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  tag           Tag of the exchange 
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_irecv_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
);


/**
 *
 * \brief Initialize a asynchronus reverse irecv (from part2)
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  s_data        Data size
 * \param [in]  cst_stride    Constant stride
 * \param [in]  part1_data    Partition 1 data (order given by part1_to_part2 array)
 * \param [in]  tag           Tag of the exchange 
 * \param [out] request       Request
 *
 */

void
PDM_part_to_part_reverse_irecv
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void              **part2_data,
 const int           tag,
 int                *request
);


/**
 *
 * \brief Initialize a asynchronus reverse irecv (from part2)
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  tag           Tag of the exchange 
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_reverse_irecv_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
);


/**
 *
 * \brief Initialize a partial asynchronus exchange
 *
 * \param [in]   ptp              Part to part structure
 * \param [in]   k_commm          Kind of MPI communication
 * \param [in]   t_stride         Kind of stride
 * \param [in]   t_part1_data_def Kind of part1 data definition
 * \param [in]   cst_stride       Constant stride
 * \param [in]   s_data           Data size
 * \param [in]   part1_stride     Stride of partition 1 data 
 * \param [in]   part1_data       Partition 1 data 
 * \param [out]  part2_stride     Stride of partition 2 data (order given by gnum1_come_from and ref_gnum2 arrays)
 * \param [out]  part2_data       Partition 2 data (order given by gnum1_come_from and ref_gnum2 arrays)
 * \param [out]  request          Request
 *
 */

void
PDM_part_to_part_iexch
(
 PDM_part_to_part_t                *ptp,
 const PDM_mpi_comm_kind_t          k_comm, 
 const PDM_stride_t                 t_stride,
 const PDM_part_to_part_data_def_t  t_part1_data_def,
 const int                          cst_stride,
 const size_t                       s_data,
 const int                        **part1_stride,
 const void                       **part1_data,
 int                             ***part2_stride,
 void                            ***part2_data,
 int                               *request
);


/**
 *
 * \brief Wait a partial asynchronus exchange
 *
 * \param [in]  ptp      Part to part structure
 * \param [in]  request  Request
 *
 */

void
PDM_part_to_part_iexch_wait
(
 PDM_part_to_part_t                *ptp,
 const int                          request
);


/**
 *
 * \brief Initialize a partial reverse asynchronus exchange
 *
 * \param [in]   ptp              Part to part structure
 * \param [in]   k_commm          Kind of MPI communication
 * \param [in]   t_stride         Kind of stride
 * \param [in]   t_part2_data_def Kind of part2 data definition
 * \param [in]   cst_stride       Constant stride
 * \param [in]   s_data           Data size
 * \param [in]   part2_stride     Stride of partition 1 data (Accordding to t_part2_data_def) 
 * \param [in]   part2_data       Partition 1 data (Accordding to t_part2_data_def)
 * \param [out]  part1_stride     Stride of partition 2 data (order given by part1_to_part2)
 * \param [out]  part1_data       Partition 2 data (order given by part1_to_part2)
 * \param [out]  request          Request
 *
 */

void
PDM_part_to_part_reverse_iexch
(
 PDM_part_to_part_t                *ptp,
 const PDM_mpi_comm_kind_t          k_comm, 
 const PDM_stride_t                 t_stride,
 const PDM_part_to_part_data_def_t  t_part2_data_def,
 const int                          cst_stride,
 const size_t                       s_data,
 const int                        **part2_stride,
 const void                       **part2_data,
 int                             ***part1_stride,
 void                            ***part1_data,
 int                               *request
);


/**
 *
 * \brief Wait a partial asynchronus exchange
 *
 * \param [in]  ptp      Part to part structure
 * \param [in]  request  Request
 *
 */

void
PDM_part_to_part_reverse_iexch_wait
(
 PDM_part_to_part_t                *ptp,
 const int                          request
);


/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] ptp  Block to part structure
 *
 * \return       NULL
 */

PDM_part_to_part_t *
PDM_part_to_part_free
(
 PDM_part_to_part_t *ptp
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_tO_pART_H */
