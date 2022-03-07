#ifndef __PDM_DCUBE_NODAL_GEN_H__
#define __PDM_DCUBE_NODAL_GEN_H__

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_ho_ordering.h"
#include "pdm_domain_interface.h"

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
 * Type definitions
 *============================================================================*/

typedef struct _pdm_dcube_nodal_t PDM_dcube_nodal_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Initialize a \ref PDM_dcube_nodal_t structure
 *
 * \param [in]   comm           MPI communicator
 * \param [in]   n_vtx_x        Number of vertices on segments in x-direction
 * \param [in]   n_vtx_y        Number of vertices on segments in y-direction
 * \param [in]   n_vtx_z        Number of vertices on segments in z-direction
 * \param [in]   length         Segment length
 * \param [in]   zero_x         X-coordinate of the origin
 * \param [in]   zero_y         Y-coordinate of the origin
 * \param [in]   zero_z         Z-coordinate of the origin
 * \param [in]   t_elt          Element type
 * \param [in]   order          Element order
 * \param [in]   owner          Ownership
 *
 * \return   Pointer to new \ref PDM_dcube_nodal_t object
 *
 */

PDM_dcube_nodal_t *
PDM_dcube_nodal_gen_create
(
 PDM_MPI_Comm          comm,
 const PDM_g_num_t     n_vtx_x,
 const PDM_g_num_t     n_vtx_y,
 const PDM_g_num_t     n_vtx_z,
 const double          length,
 const double          zero_x,
 const double          zero_y,
 const double          zero_z,
 PDM_Mesh_nodal_elt_t  t_elt,
 const int             order,
 PDM_ownership_t       owner
);


/**
 *
 * \brief Free a \ref PDM_dcube_nodal_t structure
 *
 * \param [in]  dcube      Pointer to \ref PDM_dcube_nodal_t object
 *
 */

void
PDM_dcube_nodal_gen_free
(
 PDM_dcube_nodal_t *dcube
);


/**
 *
 * \brief Set the HO-ordering for a \ref PDM_dcube_nodal_t structure
 *
 * \param [in]  dcube      Pointer to \ref PDM_dcube_nodal_t object
 * \param [in]  ordering   Name of the HO-ordering
 *
 */

void PDM_dcube_nodal_gen_ordering_set
(
       PDM_dcube_nodal_t *dcube,
 const char              *ordering
);


/**
 *
 * \brief Build a \ref PDM_dcube_nodal_t structure
 *
 * \param [in]  dcube      Pointer to \ref PDM_dcube_nodal_t object
 *
 * \return   Pointer to the associated \ref PDM_dmesh_nodal_t object
 *
 */

PDM_dmesh_nodal_t *
PDM_dcube_nodal_gen_build
(
 PDM_dcube_nodal_t *dcube
);


/**
 *
 * \brief Get the \ref PDM_dmesh_nodal_t associated to a \ref PDM_dcube_nodal_t
 *
 * \param [in]  dcube      Pointer to \ref PDM_dcube_nodal_t object
 *
 * \return   Pointer to the associated \ref PDM_dmesh_nodal_t object
 *
 */

PDM_dmesh_nodal_t *
PDM_dcube_nodal_gen_dmesh_nodal_get
(
 PDM_dcube_nodal_t *dcube
);


void
PDM_dcube_nodal_cart_topo
(
       PDM_MPI_Comm              comm,
       int                       n_dom_i,
       int                       n_dom_j,
       int                       n_dom_k,
       int                       periodic_i,
       int                       periodic_j,
       int                       periodic_k,
 const PDM_g_num_t               n_vtx_x_in,
 const PDM_g_num_t               n_vtx_y_in,
 const PDM_g_num_t               n_vtx_z_in,
 const double                    length,
 const double                    zero_x,
 const double                    zero_y,
 const double                    zero_z,
       PDM_Mesh_nodal_elt_t      t_elt,
 const int                       order,
       PDM_dcube_nodal_t      ***dcube,
       PDM_domain_interface_t  **dom_intrf,
       PDM_ownership_t           owner
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_DCUBE_NODAL_GEN_H__ */
