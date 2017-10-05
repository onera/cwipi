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

int 
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

int 
PDM_part_renum_face_add
(
 const char                 *name,     /*!< Name          */ 
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_face function for the format */             
);        


/**
 *
 * \brief Get index of a renumbering cell method
 * 
 * \param [in]  name   Name of the method
 * 
 * \return Index (-1 if not found)
 */

int 
PDM_part_renum_cell_method_idx_get
(
const char *name
);


/**
 *
 * \brief Get index of a renumbering face method
 * 
 * \param [in]  name   Name of the method
 * 
 * \return Index (-1 if not found)
 */

int 
PDM_part_renum_face_method_idx_get
(
const char *name
);        


/**
 *
 * \brief Get name of the cell renumbering method 
 * 
 * \param [in]  idx     Index of the method
 * 
 * \return Name of the method
 *
 */

const char * 
PDM_part_renum_cell_method_name_get
(
const int idx
);        


/**
 *
 * \brief Get name of the face renumbering method 
 * 
 * \param [in]  idx     Index of the method
 * 
 * \return Name of the method
 *
 */

const char * 
PDM_part_renum_face_method_name_get
(
const int idx
);        


/**
 *
 * \brief Get the number of renumbering face methods 
 * 
 * \return Name of the method
 *
 */

int  
PDM_part_n_renum_cell_method_get
(
void 
);        


/**
 *
 * \brief Get the number of renumbering face methods 
 * 
 * \return Name of the method
 *
 */

int  
PDM_part_n_renum_face_method_get
(
void 
);        


/**
 *
 * \brief Purge renumbering methods 
 *
 */

void 
PDM_part_renum_purge
(
void
);        

/**
 *
 * \brief Load local renumbering methods 
 *
 */

void 
PDM_part_renum_load_local
(
void
);        

  
/**
 *
 * \brief Perform cell renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void 
PDM_part_renum_cell
(
 _PDM_part_t   *ppart                
);       


/**
 *
 * \brief Perform face renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void 
PDM_part_renum_face
(
  _PDM_part_t   *ppart                
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

