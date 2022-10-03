#ifndef __PDM_ABSTRACT_DISTRIB_H__
#define __PDM_ABSTRACT_DISTRIB_H__
/*----------------------------------------------------------------------------*/

#include <stdbool.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_morton.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro and type definitions
 *============================================================================*/


/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Create a \ref PDM_abstract_distrib_t structure.
 *
 */

PDM_abstract_distrib_t *
PDM_abstract_distrib_create
(
	const int             n_part,
	const int             n_track_data,
  const int            *n_elt,
  const PDM_g_num_t   **gnum_elt,
        PDM_MPI_Comm    comm
);


/**
 * \brief Free a \ref PDM_box_distrib_t structure.
 *
 * \param [in]  ad Pointer to the structure to destroy
 */

void
PDM_abstract_distrib_free
(
	PDM_abstract_distrib_t  ad
);


/**
 * \brief Set a tracked data
 *
 * \param [in]  ad Pointer to the structure 
 */

void
PDM_abstract_distrib_track_data_set 
(
	      PDM_box_distrib_t  ad,
  const int                i_data,
  const PDM_stride_t       t_stride,
  const int                stride_cst,
  const size_t             data_size,
	      int              **data_stride,
	      void             **data
);


/**
 * \brief Get a tracked data
 *
 * \param [in]  ad Pointer to the structure 
 */

void
PDM_abstract_distrib_track_data_get 
(
	      PDM_box_distrib_t  ad,
  const int                i_data,
        PDM_stride_t      *t_stride,
        int               *stride_cst,
        size_t            *data_size,
	      int              **data_stride,
	      void             **data
);

/**
 * \brief Redistribute the distribution and its tracked data
 *
 * \param [in]  ad Pointer to the structure to destroy
 */

void
PDM_abstract_distrib_redistribute_start
(
	PDM_abstract_distrib_t  ad,
  const int              *n_part,
  const int              *n_elt,
  const PDM_g_num_t     **gnum_elt,
);

/**
 * \brief Redistribute the distribution and its tracked data
 *
 * \param [in]  ad Pointer to the structure to destroy
 */

void
PDM_abstract_distrib_redistribute_wait
(
	PDM_abstract_distrib_t  ad
);


/**
 * \brief Receive data from origin distribution 
 *
 * \param[in]     ad                     Pointer to the structure 
 * \param[in]     t_stride               Type of stride
 * \param[in]     stride_cst             Constant stride
 * \param[in]     data_size              Size of data
 * \param[in]     origin_distrib_stride  Origin stride distribution
 * \param[in]     origin_distrib_data    Origin data distribution
 * \param[in,out] current_distrib_stride Current stride distribution (Allocate and compute if input is NULL,
 *                                       otherwise nothing)
 * \param[out]    current_distrib_data   Current data distribution
 *
 * \return  Size of current_distrib_data
 */

void
PDM_abstract_distrib_recv_data_from_origin_start
(
	      PDM_abstract_distrib_t  ad,
  const PDM_stride_t            t_stride,
  const int                     stride_cst,
  const size_t                  data_size,
        int                   **current_distrib_stride,
        void                  **current_distrib_data,
        int                   **origin_distrib_stride,
        void                  **origin_distrib_data,
);

/**
 * \brief Receive data from origin distribution 
 *
 * \param[in]     ad                     Pointer to the structure 
 */

void
PDM_abstract_distrib_recv_data_from_origin_wait
(
	      PDM_abstract_distrib_t  ad
);


/**
 * \brief Send data to origin for any box
 *
 * \param [in]  box_set                 pointer to the \ref PDM_box_t structure
 * \param [in]  t_stride                Type of stride
 * \param [in]  stride_cst              Constant stride
 * \param [in]  data_size               Size of data
 * \param [in]  current_distrib_stride  Current stride distribution
 * \param [in]  current_distrib_data    Current data distribution
 * \param [out] origin_distrib_stride   Origin stride distribution
 * \param [out] origin_distrib_data     Origin data distribution
 *
 * \return   Size of origin_distrib_data
 */

void
PDM_abstract_distrib_send_data_to_origin_start
(
	      PDM_abstract_distrib_t  ad,
 const  PDM_stride_t            t_stride,
 const  int                     stride_cst,
 const  size_t                  data_size,
        int                    *current_distrib_stride,
        void                   *current_distrib_data,
        int                   **origin_distrib_stride,
        void                  **origin_distrib_data
);


void
PDM_abstract_distrib_send_data_to_origin_wait
(
	      PDM_abstract_distrib_t  ad,
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_ABSTRACT_DISTRIB_H__ */
