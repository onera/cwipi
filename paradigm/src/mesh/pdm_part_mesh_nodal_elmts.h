#ifndef __PDM_PART_MESH_NODAL_ELMTS_H__
#define __PDM_PART_MESH_NODAL_ELMTS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"
#include "pdm_mesh_nodal.h"

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

typedef struct _pdm_part_mesh_nodal_elmts_t PDM_part_mesh_nodal_elmts_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */
PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_elmts_create
(
const int          mesh_dimension,
const int          n_part,
const PDM_MPI_Comm comm
);

int
PDM_part_mesh_nodal_elmts_add
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const PDM_Mesh_nodal_elt_t         t_elt
);

void
PDM_part_mesh_nodal_elmts_free
(
 PDM_part_mesh_nodal_elmts_t *pmne
);

void
PDM_part_mesh_nodal_elmts_std_set
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part,
const int                          n_elt,
const int                         *connec,
const PDM_g_num_t                 *numabs,
const int                         *parent_num,
      PDM_ownership_t              owner
);

void
PDM_part_mesh_nodal_elmts_block_std_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_block,
const int                           id_part,
      int                         **connec
);

int
PDM_part_mesh_nodal_elmts_block_n_elt_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ELMTS_H__ */
