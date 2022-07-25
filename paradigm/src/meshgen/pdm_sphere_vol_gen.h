#ifndef __PDM_SPHERE_VOL_GEN_H__
#define __PDM_SPHERE_VOL_GEN_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_mesh_nodal.h"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/**
 *
 * \brief Create a volume mesh bounded by a sphere (deformed cube)
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n_vtx_x         Number of vertices on segments in x-direction
 * \param[in]  n_vtx_y         Number of vertices on segments in y-direction
 * \param[in]  n_vtx_z         Number of vertices on segments in z-direction
 * \param[in]  radius          Radius of the sphere
 * \param[in]  center_x        x coordinate of the center of the sphere
 * \param[in]  center_y        y coordinate of the center of the sphere
 * \param[in]  center_z        z coordinate of the center of the sphere
 * \param[in]  t_elt           Element type
 * \param[in]  order           Element order
 * \param[out] dmn             Pointer to a \ref PDM_dmesh_nodal object
 *
 */

void
PDM_sphere_vol_gen_nodal
(
 PDM_MPI_Comm           comm,
 const PDM_g_num_t      n_vtx_x,
 const PDM_g_num_t      n_vtx_y,
 const PDM_g_num_t      n_vtx_z,
 const double           radius,
 const double           center_x,
 const double           center_y,
 const double           center_z,
 PDM_Mesh_nodal_elt_t   t_elt,
 const int              order,
 PDM_dmesh_nodal_t    **dmn
 );

#ifdef __cplusplus
}
#endif
#endif
