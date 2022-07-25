#ifndef __PDM_SPHERE_SURF_GEN_H__
#define __PDM_SPHERE_SURF_GEN_H__
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
 * \brief Create a surface mesh of a sphere
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  nu              Number of points in longitude
 * \param[in]  nv              Number of points in latitude
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius          Radius of the sphere
 * \param[out] dvtx_coord      Connectivity of distributed vertex to coordinates
 * \param[out] dface_vtx_idx   Index of distributed face to vertex
 * \param[out] dface_vtx       Connectivity of distributed face to vertex
 * \param[out] distrib_vtx     Distribution of vertices
 * \param[out] distrib_face    Distribution of faces
 *
 */

void
PDM_sphere_surf_gen
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         nu,
 const PDM_g_num_t         nv,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       double            **dvtx_coord,
       int               **dface_vtx_idx,
       PDM_g_num_t       **dface_vtx,
       PDM_g_num_t       **distrib_vtx,
       PDM_g_num_t       **distrib_face
);


/**
 *
 * \brief Create a surface mesh of a sphere
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  nu              Number of points in longitude
 * \param[in]  nv              Number of points in latitude
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius          Radius of the sphere
 * \param[out] _dmn            Sphere mesh in the form of a distributed nodal mesh
 *
 */

void
PDM_sphere_surf_gen_nodal
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         nu,
 const PDM_g_num_t         nv,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       PDM_dmesh_nodal_t **_dmn
);


/**
 *
 * \brief Create a surface mesh of a sphere (icosphere)
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Number of icosphere subdivisions
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius          Radius of the sphere
 * \param[out] dvtx_coord      Connectivity of distributed vertex to coordinates
 * \param[out] dface_vtx_idx   Index of distributed face to vertex
 * \param[out] dface_vtx       Connectivity of distributed face to vertex
 * \param[out] distrib_vtx     Distribution of vertices
 * \param[out] distrib_face    Distribution of faces
 *
 */

void
PDM_sphere_surf_icosphere_gen
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       double            **dvtx_coord,
       int               **dface_vtx_idx,
       PDM_g_num_t       **dface_vtx,
       PDM_g_num_t       **distrib_vtx,
       PDM_g_num_t       **distrib_face
);


/**
 *
 * \brief Create a surface mesh of a sphere (icosphere)
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Number of icosphere subdivisions
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius          Radius of the sphere
 * \param[out] _dmn            Sphere mesh in the form of a distributed nodal mesh
 *
 */

void
PDM_sphere_surf_icosphere_gen_nodal
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       PDM_dmesh_nodal_t **_dmn
);

#ifdef __cplusplus
}
#endif
#endif
