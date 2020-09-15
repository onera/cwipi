#ifndef __PDM_ORDER_H__
#define __PDM_ORDER_H__

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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Order an array
 *
 * \param [in]      size_array       Number of elements
 * \param [in]      new_to_old_order New order (size = \ref nElt
 * \param [in, out] Array            Array to renumber
 *
 */

void
PDM_order_array
(
const int     size_array,
const size_t  elt_size,
const int    *new_to_old_order,
void         *array
);

/**
 * This function is part of Code_Saturne, a general-purpose CFD tool.
 *  Copyright (C) 1998-2014 EDF S.A.
 *
 * \brief Order a strided array of global numbers lexicographically.
 *
 * \param [in]     number array of entity numbers (if NULL, a default 1 to n numbering is considered)
 * \param [in]     stride stride of array (number of values to compare)
 * \param [in,out] order  pre-allocated ordering table
 * \param [in]     nb_ent number of entities considered
 */

void
PDM_order_lnum_s
(
const int    number[],
size_t       stride,
int          order[],
const size_t nb_ent
);



void
PDM_order_gnum_s
(
const PDM_g_num_t number[],
size_t            stride,
int               order[],
const size_t      nb_ent
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_BINARY_SEARCH_H__ */
