#ifndef __PDM_DCUBE_NODAL_GEN_H__
#define __PDM_DCUBE_NODAL_GEN_H__

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"

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

PDM_dcube_nodal_t*
PDM_dcube_nodal_gen_init
(
      PDM_MPI_Comm          comm,
const PDM_g_num_t           n_vtx_seg,
const double                length,
const double                zero_x,
const double                zero_y,
const double                zero_z,
      PDM_Mesh_nodal_elt_t  t_elt,
      PDM_ownership_t       owner
);

/**
 *
 * \brief Return distributed cube size
 *
 * \param [in]   id          dcube_nodal identifier
 * \param [out]  n_face_group  Number of faces groups
 * \param [out]  dn_cell       Number of cells stored in this process
 * \param [out]  dn_face       Number of faces stored in this process
 * \param [out]  dn_vtx        Number of vertices stored in this process
 * \param [out]  sface_vtx     Length of dface_vtx array
 * \param [out]  sface_group   Length of dface_group array
 *
 */

void
PDM_dcube_nodal_gen_dim_get
(
 PDM_dcube_nodal_t  *pdm_dcube_nodal,
 int                *n_face_group,
 int                *dn_cell,
 int                *dn_face,
 int                *dn_vtx,
 int                *sface_vtx,
 int                *sface_group
);

/**
 *
 * \brief Return distributed cube data
 *
 * \param [in]  id              dcube_nodal identifier
 * \param [out] dface_cell      Faces from cells connectivity (size = 2 * dn_face)
 * \param [out] dface_vtx_idx   Faces from vertices connectivity index (size = dn_face + 1)
 * \param [out] dface_vtx       Faces from vertices connectivity (size = dface_vtxL)
 * \param [out] dvtx_coord      Vertices coordinates (size = 3 * dn_vtx)
 * \param [out] dface_group_idx Faces groups index (size = n_face_group + 1)
 * \param [out] dface_group     Faces groups (size = dFacegroupL)
 *
 */

void
PDM_dcube_nodal_gen_data_get
(
 PDM_dcube_nodal_t  *pdm_dcube_nodal,
 PDM_g_num_t       **delmt_vtx,
 double            **dvtx_coord,
 int               **dface_group_idx,
 PDM_g_num_t       **dface_group
);



PDM_dmesh_nodal_t*
PDM_dcube_nodal_gen_dmesh_nodal_get
(
 PDM_dcube_nodal_t  *pdm_dcube_nodal
);

/**
 *
 * \brief Free a distributed cube
 *
 * \param [in]  id            dcube_nodal identifier
 *
 */

void
PDM_dcube_nodal_gen_free
(
 PDM_dcube_nodal_t       *pdm_dcube_nodal
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PART_DCUBE_NODAL_H__ */
