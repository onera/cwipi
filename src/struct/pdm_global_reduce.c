
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_morton.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_global_reduce_priv.h"
#include "pdm_global_reduce.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/


/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure that computes a global reduction
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Pointer to global reduction object
 */

PDM_global_reduce_t *
PDM_global_reduce_create
(
 const int          n_part,
 const PDM_MPI_Comm comm
)
{
  PDM_global_reduce_t *gre = (PDM_global_reduce_t *) malloc (sizeof(PDM_global_reduce_t));

  gre->n_part  = n_part;
  gre->comm    = comm;
  gre->g_nums  = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part);
  gre->n_elts  = (int          *) malloc (sizeof(int          ) * n_part);
  gre->strides = (int         **) malloc (sizeof(int         *) * n_part);
  gre->ptb     = NULL;
  gre->btp     = NULL;

  gre->operation = PDM_REDUCE_OP_MEAN;

  gre->s_weight = NULL;

  for (int i = 0; i < n_part; i++) {
    gre->g_nums [i] = NULL;
    gre->strides[i] = NULL;
  }

  gre->local_field          = (double **) malloc (sizeof(double * ) * n_part);
  gre->local_weight         = (double **) malloc (sizeof(double * ) * n_part);
  gre->global_reduced_field = (double **) malloc (sizeof(double * ) * n_part);

  for (int i = 0; i < n_part; i++) {
    gre->local_field         [i] = NULL;
    gre->local_weight        [i] = NULL;
    gre->global_reduced_field[i] = NULL;
  }

  return gre;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
