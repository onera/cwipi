#ifndef __pdm_dmesh_nodal_tO_DMESH_PRIV_H__
#define __pdm_dmesh_nodal_tO_DMESH_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct _pdm_dmesh_nodal_to_dmesh_t
 * \brief  link between dmesh_nodal and dmesh
 *
 */

typedef struct _pdm_link_dmesh_nodal_to_dmesh_t {

  PDM_dmesh_nodal_t *dmesh_nodal;
  PDM_dmesh_t       *dmesh;
  PDM_ownership_t    owner;                   /*!< Which have the responsabilities of results */
  PDM_bool_t         results_is_getted;       /*!< Flags to indicate if result is getted      */

  /* We keep the link between the dmesh_nodal and dmesh */
  int           dn_elmt;
  PDM_g_num_t  *elmt_distrib;
  PDM_g_num_t  *_delmt_face;
  int          *_delmt_face_idx;

  PDM_g_num_t  *_dface_elmt;
  int          *_dface_elmt_idx;

  PDM_g_num_t  *_delmt_edge;
  int          *_delmt_edge_idx;

  PDM_g_num_t  *_dedge_elmt;
  int          *_dedge_elmt_idx;

} _pdm_link_dmesh_nodal_to_dmesh_t;


struct _pdm_dmesh_nodal_to_dmesh_t {

  PDM_MPI_Comm         comm;                    /*!< MPI communicator */
  PDM_ownership_t      owner;                   /*!< Which have the responsabilities of results */
  PDM_bool_t           results_is_getted;       /*!< Flags to indicate if result is getted      */

  int                  n_mesh;                  /*!< Number of meshes to manages                */

  _pdm_link_dmesh_nodal_to_dmesh_t **link;

};

#ifdef  __cplusplus
}
#endif

#endif  /* __pdm_dmesh_nodal_tO_DMESH_PRIV_H__ */
