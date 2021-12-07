#ifndef __PDM_HO_ORDERING_H_
#define __PDM_HO_ORDERING_H_


/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_hash_tab.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/


/*=============================================================================
 * Global variables
 *============================================================================*/


/*=============================================================================
 * Public function interfaces
 *============================================================================*/

void
PDM_ho_ordering_init
(
void
);

void
PDM_ho_ordering_free
(
void
);

int
PDM_ho_ordering_user_to_ijk_add
(
 const char                 *name,
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const int                   n_nodes,
 const int                  *user_to_ijk
);


int *
PDM_ho_ordering_user_to_ijk_get
(
 const char                 *name,
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order
);

int *
PDM_ho_ordering_ijk_to_user_get
(
 const char                 *name,
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order
);


int *
PDM_ho_ordering_compute_ijk_to_user
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const int                  *user_to_ijk
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HO_ORDERING_H__ */
