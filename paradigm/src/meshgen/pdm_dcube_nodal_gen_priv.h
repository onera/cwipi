#ifndef __PDM_DCUBE_NODAL_GEN_PRIV_H__
#define __PDM_DCUBE_NODAL_GEN_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_dcube_t
 * \brief  Distributed cube
 *
 * _dcube_t define a distributed mesh of a cube
 *
 */

struct _pdm_dcube_nodal_t {
  PDM_MPI_Comm          comm;                   /*!< MPI communicator                          */
  PDM_ownership_t       owner;                  /*!< Which have the responsabilities of results*/
  PDM_ownership_t       owner_for_dmesh_nodal;  /*!< Which have the responsabilities of results*/
  PDM_g_num_t           n_vtx_seg;              /*!< Number of vertices in segments            */
  double                length;                 /*!< Segment length                            */
  double                zero_x;                 /*!< Coordinates of the origin                 */
  double                zero_y;                 /*!< Coordinates of the origin                 */
  double                zero_z;                 /*!< Coordinates of the origin                 */

  int                   dn_hexa_cell;           /*!< Number of hexa stored in this process     */
  int                   dn_quad_lim;            /*!< Number of quad stored in this process     */
  int                   dn_quad_seq_lim;        /*!< Number of quad stored in this process     */

  int                   n_face_group;           /*!< Number of faces groups                    */
  int                   dn_cell;                /*!< Number of cells stored in this process    */
  int                   dn_vtx;                 /*!< Number of vertices stored in this process */

  PDM_g_num_t           n_g_hexa_cell;
  PDM_g_num_t           n_g_hexa_cell_seg;
  PDM_g_num_t*          distrib_hexa;
  PDM_g_num_t*          distrib_quad_lim;
  PDM_g_num_t*          distrib_quad_seg_lim;

  PDM_g_num_t          *delmt_vtx;              /*!< Faces from vertices connectivity          */
  PDM_g_num_t          *delmt_lim_vtx;          /*!< Faces from vertices connectivity          */
  double               *dvtx_coord;             /*!< Vertices coordinates                      */
  int                  *dface_group_idx;        /*!< Faces groups index                        */
  PDM_g_num_t          *dface_group;            /*!< Faces groups                              */
  PDM_Mesh_nodal_elt_t  t_elt;                  /*!< Type of elements to generate              */

  PDM_dmesh_nodal_t    *dmesh_nodal;            /*!< Results                                   */

  int                   api_type;

} ;


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_DCUBE_NODAL_GEN_PRIV_H__ */
