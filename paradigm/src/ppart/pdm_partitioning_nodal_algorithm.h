#ifndef __PDM_PARTITIONING_NODAL_ALGORITHM_H__
#define __PDM_PARTITIONING_NODAL_ALGORITHM_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"

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

void
PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts
(
 PDM_dmesh_nodal_elmts_t      *dmne,
 PDM_part_mesh_nodal_elmts_t  *pmne,
 int                           n_part,
 int                          *pn_elmt,
 PDM_g_num_t                 **elmt_ln_to_gn
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PARTITIONING_NODAL_ALGORITHM_H__ */
