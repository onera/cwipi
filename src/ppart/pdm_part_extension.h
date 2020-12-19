#ifndef __PDM_PART_EXTENSION_H__
#define __PDM_PART_EXTENSION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
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

typedef struct _pdm_part_extension_t PDM_part_extension_t;


typedef enum {

  PDM_EXTEND_FROM_FACE = 0,
  PDM_EXTEND_FROM_EDGE = 1,
  PDM_EXTEND_FROM_VTX  = 2

} PDM_extend_type_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


PDM_part_extension_t*
PDM_part_extension_create
(
 const int              n_domain,
 const int             *n_part,
 const PDM_MPI_Comm     comm,
 const PDM_ownership_t  owner
);

/**
 *
 * \brief Compute a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_compute
(
  PDM_part_extension_t *part_ext,
  int                   depth,
  PDM_extend_type_t     extend_type
);

void
PDM_part_extension_set_part
(
  PDM_part_extension_t *part_ext,
  int                   i_domain,
  int                   i_part,
  int                  *cell_face_idx,
  int                  *cell_face,
  int                  *face_cell,
  int                  *face_edge_idx,
  int                  *face_edge,
  int                  *face_vtx_idx,
  int                  *face_vtx,
  int                  *edge_vtx,
  int                  *face_bound_idx,
  int                  *face_bound,
  int                  *face_join_idx,
  int                  *face_join,
  int                  *face_part_bound_proc_idx,
  int                  *face_part_bound_part_idx,
  int                  *face_part_bound,
  PDM_g_num_t          *cell_ln_to_gn,
  PDM_g_num_t          *face_ln_to_gn,
  PDM_g_num_t          *edge_ln_to_gn,
  PDM_g_num_t          *vtx_ln_to_gn,
  PDM_g_num_t          *face_group_ln_to_gn
);

void
PDM_part_extension_free
(
 PDM_part_extension_t *part_ext
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_H__ */
