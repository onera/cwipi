#ifndef __PDM_PART_RENUM_H__
#define	__PDM_PART_RENUM_H__

/*============================================================================
 * Mesh entities renumbering 
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part.h"
#include "pdm_part_priv.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro and type definitions
 *============================================================================*/


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Perform mesh entities renumbering
 *
 * \param [in,out]  ppart       ppart structure
 * \param [in]      entity      Mesh entity to renumber
 * \param [in]      method      Renumbering method
 *
 */

void 
PDM_part_renum_cell
(
 _PDM_part_t           *ppart,
 PDM_part_renum_cell_t method                 
);        


/**
 *
 * \brief Perform mesh entities renumbering
 *
 * \param [in,out]  ppart       ppart structure
 * \param [in]      entity      Mesh entity to renumber
 * \param [in]      method      Renumbering method
 *
 */

void 
PDM_part_renum_face
(
 _PDM_part_t              *ppart,
 PDM_part_renum_face_t  method                 
);        


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_RENUM_H */

