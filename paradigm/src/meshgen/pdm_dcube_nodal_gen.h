#ifndef __PDM_DCUBE_NODAL_GEN_H__
#define __PDM_DCUBE_NODAL_GEN_H__

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_ho_ordering.h"

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

PDM_dcube_nodal_t *
PDM_dcube_nodal_gen_create
(
 PDM_MPI_Comm          comm,
 const PDM_g_num_t     nx,
 const PDM_g_num_t     ny,
 const PDM_g_num_t     nz,
 const double          length,
 const double          zero_x,
 const double          zero_y,
 const double          zero_z,
 PDM_Mesh_nodal_elt_t  t_elt,
 const int             order,
 PDM_ownership_t       owner
);


void
PDM_dcube_nodal_gen_free
(
 PDM_dcube_nodal_t *dcube
 );

void PDM_dcube_nodal_gen_ordering_set
(
 PDM_dcube_nodal_t *dcube,
 char              *ordering
 );


PDM_dmesh_nodal_t *
PDM_dcube_nodal_gen_build
(
 PDM_dcube_nodal_t *dcube
 );

/**
 *
 * \brief Create a distributed cube
 *
 * \param [out]  id             dcube_nodal identifier
 * \param [in]   comm           Communicator
 * \param [in]   n_vtx_seg        Number of vertices in segments
 * \param [in]   length         Segment length
 * \param [in]   zero_x         Coordinates of the origin
 * \param [in]   zero_y         Coordinates of the origin
 * \param [in]   zero_z         Coordinates of the origin
 *
 */

PDM_dmesh_nodal_t *
PDM_dcube_nodal_gen_dmesh_nodal_get
(
 PDM_dcube_nodal_t *dcube
 );



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_DCUBE_NODAL_GEN_H__ */
