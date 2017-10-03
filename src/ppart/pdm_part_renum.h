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

/**
 * \struct PDM_part_renum_fct_t
 *
 * \brief  Function pointer used to define a renumbering method 
 *
 */

typedef void (*PDM_part_renum_fct_t) (_PDM_part_t  *ppart);  

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Add a new method for cell renumbering 
 *
 * \param [in]      name           Mesh entity to renumber
 * \param [in]      renum_fct      Renumbering function
 *
 */

void 
PDM_part_renum_cell_add
(
 const char                 *name,     /*!< Name          */ 
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_cell function for the format */             
);        

/**
 *
 * \brief Add a new method for face renumbering 
 *
 * \param [in]      name           Mesh entity to renumber
 * \param [in]      renum_fct      Renumbering function
 *
 */

void 
PDM_part_renum_face_add
(
 const char                 *name,     /*!< Name          */ 
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_face function for the format */             
);        

/**
 *
 * \brief Purge renumbering methods 
 *
 */

void 
PDM_part_renum_purge
(
);        

/**
 *
 * \brief Load local renumbering methods 
 *
 */

void 
PDM_part_renum_load_local
(
);        

  
/**
 *
 * \brief Perform cell renumbering
 *
 * \param [in]      method      Renumbering method
 *
 */

void 
PDM_part_renum_cell
(
const char *method
);        


/**
 *
 * \brief Perform face renumbering
 *
 * \param [in]      method      Renumbering method
 *
 */

void 
PDM_part_renum_face
(
const char *method
);        

/**
 *
 * \brief Perform cells renumbering from a new order 
 *        Actualise all cells array according to the new numbering 
 *        Connectivities/cellTag/cellColor/cellLNToGN
 *
 * \param [in,out]  part        Current partition
 * \param [in]      newToOldOrder    NewOrder
 *
 */

void 
PDM_part_reorder_cell
(
 _part_t *part, 
 int     *newToOldOrder               
);        


/**
 *
 * \brief Perform faces renumbering from a new order 
 *        Actualise all cells array according to the new numbering 
 *        Connectivities/faceTag/faceColor/faceLNToGN
 *
 * \param [in,out]  part        Current partition
 * \param [in]      newToOldOrder    NewOrder
 *
 */

void 
PDM_part_reorder_face
(
 _part_t *part, 
 int     *newToOldOrder               
);        

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_RENUM_H */

