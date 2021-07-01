#ifndef PDM_DISTRIB_H
#define PDM_DISTRIB_H

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

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
 * \brief Compute distribution from dNelmt
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     dnelt                Number of element on current process
 * \param [in]     comm                 MPI Communicator
 */
void
PDM_distrib_compute
(
 const int           dnelt,
       PDM_g_num_t  *elt_distrib,
       int           offset,
 const PDM_MPI_Comm  comm
);

/**
 * \brief Compute distribution from dNelmt
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     dnelt                Number of element on current process
 * \param [in]     comm                 MPI Communicator
 */
PDM_g_num_t*
PDM_compute_entity_distribution
(
 const PDM_MPI_Comm     comm,
 const int              dn_entity
);


int
PDM_compute_uniform_dn_entity
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t      n_g_entity
);


PDM_g_num_t*
PDM_compute_uniform_entity_distribution
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t      n_g_entity
);


PDM_g_num_t*
PDM_compute_uniform_entity_distribution_from_partition
(
 const PDM_MPI_Comm     comm,
 const int              n_part,
 const int             *n_elmts,
 const PDM_g_num_t    **ln_to_gn
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif // PDM_DISTRIB_H
