#ifndef __PDM_DOMAIN_INTERFACE_H__
#define __PDM_DOMAIN_INTERFACE_H__

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

typedef enum {

  PDM_DOMAIN_INTERFACE_MULT_NO  = 0,  /*!< Each interface involves only 2 zones */
  PDM_DOMAIN_INTERFACE_MULT_YES = 1,  /*!< Each interface involves several zones */

} PDM_domain_interface_mult_t;

typedef struct _pdm_domain_interface_t PDM_domain_interface_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

PDM_domain_interface_t *
PDM_domain_interface_create
(
 const int                   n_interface,
 const int                   n_zone,
 PDM_domain_interface_mult_t multizone_interface,
 PDM_ownership_t             ownership,
 PDM_MPI_Comm                comm
);

void
PDM_domain_interface_set
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                    *interface_dn,
 PDM_g_num_t           **interface_ids,
 int                   **interface_dom
);

void
PDM_domain_interface_get
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                   **interface_dn,
 PDM_g_num_t          ***interface_ids,
 int                  ***interface_dom
);

void
PDM_domain_interface_translate_face2vtx
(
 PDM_domain_interface_t  *dom_intrf,
 int                     *dn_vtx,
 int                     *dn_face,
 int                    **dface_vtx_idx,
 PDM_g_num_t            **dface_vtx
);

void
PDM_domain_interface_free
(
 PDM_domain_interface_t *dom_intrf
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DOMAIN_INTERFACE_H__ */
