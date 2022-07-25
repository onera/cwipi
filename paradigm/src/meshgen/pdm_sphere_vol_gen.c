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

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_config.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_array.h"
#include "pdm_distrib.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"

#include "pdm_sphere_vol_gen.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Local macro definitions
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/


static inline double
_smooth_step
(
 const double x
 )
{
  return PDM_MIN(1., PDM_MAX(0., 3*x*x*(2 - x)));
}

/**
 *
 * -1 <= xc, yc, zc <= 1
 *
 */

static void
_cube_to_sphere
(
 const double  xc,
 const double  yc,
 const double  zc,
 double       *xs,
 double       *ys,
 double       *zs
 )
{
  double r2 = xc*xc + yc*yc + zc*zc;
  double r  = sqrt(r2);
  double rinf = PDM_MAX(PDM_MAX(PDM_ABS(xc), PDM_ABS(yc)), PDM_ABS(zc));

  double scale1 = rinf;
  if (r > 1e-15) {
    scale1 /= r;
  }
  double x1 = xc * scale1;
  double y1 = yc * scale1;
  double z1 = zc * scale1;

  double scale2 = sqrt(r2 -
                       xc*xc*yc*yc - yc*yc*zc*zc - zc*zc*xc*xc +
                       xc*xc*yc*yc*zc*zc);
  if (r > 1e-15) {
    scale2 /= r;
  }
  double x2 = xc * scale2;
  double y2 = yc * scale2;
  double z2 = zc * scale2;

  double t = 0;//0.7*_smooth_step(r);

  *xs = (1 - t)*x2 + t*x1;
  *ys = (1 - t)*y2 + t*y1;
  *zs = (1 - t)*z2 + t*z1;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

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
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int mesh_dimension = PDM_Mesh_nodal_elt_dim_get(t_elt);
  if (mesh_dimension != 3) {
    PDM_error(__FILE__, __LINE__, 0,
              "Not implemented yes for dimension %d\n", mesh_dimension);
  }

  /* First: generate a dcube nodal */
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_x,
                                                        n_vtx_y,
                                                        n_vtx_z,
                                                        2*radius,
                                                        center_x - radius,
                                                        center_y - radius,
                                                        center_z - radius,
                                                        t_elt,
                                                        order,
                                                        PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t *_dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(_dmn);

  /* Second: "spherify" */
  PDM_g_num_t *distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(_dmn);
  int dn_vtx = (int) (distrib_vtx[i_rank+1] - distrib_vtx[i_rank]);
  double *dvtx_coord  = PDM_DMesh_nodal_vtx_get(_dmn);

  double center[3] = {center_x, center_y, center_z};

  double *_dvtx_coord = malloc(sizeof(double) * dn_vtx * 3);
  for (int i = 0; i < dn_vtx; i++) {

    // double r2 = 0., rinf = 0.;
    // for (int j = 0; j < 3; j++) {
    //   double x = dvtx_coord[3*i + j] - center[j];
    //   r2 += x*x;
    //   rinf = PDM_MAX(rinf, PDM_ABS(x));
    // }

    // r2 = sqrt(r2);

    // double scale = rinf/r2;
    // for (int j = 0; j < 3; j++) {
    //   double x = dvtx_coord[3*i + j] - center[j];

    //   _dvtx_coord[3*i + j] = center[j] + scale*x;
    // }
    double xyzc[3];
    for (int j = 0; j < 3; j++) {
      xyzc[j] = (dvtx_coord[3*i + j] - center[j]) / radius;
    }

    _cube_to_sphere(xyzc[0], xyzc[1], xyzc[2],
                    _dvtx_coord + 3*i,
                    _dvtx_coord + 3*i+1,
                    _dvtx_coord + 3*i+2);

    for (int j = 0; j < 3; j++) {
      _dvtx_coord[3*i + j] = center[j] + radius*_dvtx_coord[3*i + j];
    }

  }


  /* Third: get rid of ridges and unify surfaces */

  *dmn = PDM_DMesh_nodal_create(comm,
                                mesh_dimension,
                                distrib_vtx[n_rank],
                                1,
                                0,
                                0);

  PDM_DMesh_nodal_coord_set(*dmn,
                            dn_vtx,
                            _dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  // Volume
  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
  int n_section    = PDM_DMesh_nodal_n_section_get  (_dmn, geom_kind);
  int *sections_id = PDM_DMesh_nodal_sections_id_get(_dmn, geom_kind);

  PDM_g_num_t gn_elt_vol = 0;
  for (int i_section = 0; i_section < n_section; i_section++) {

    int id_section = sections_id[i_section];
    int _order;
    const char *ho_ordering = NULL;
    const PDM_g_num_t    *distrib_elt = PDM_DMesh_nodal_distrib_section_get(_dmn, geom_kind, id_section);
    int                   dn_elt      = PDM_DMesh_nodal_section_n_elt_get  (_dmn, geom_kind, id_section);
    PDM_g_num_t          *delt_vtx    = PDM_DMesh_nodal_section_std_ho_get (_dmn, geom_kind, id_section, &_order, &ho_ordering);
    PDM_Mesh_nodal_elt_t  elt_type    = PDM_DMesh_nodal_section_type_get   (_dmn, geom_kind, id_section);

    int elt_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);
    PDM_g_num_t *_delt_vtx = malloc(sizeof(PDM_g_num_t) * dn_elt * elt_vtx_n);
    memcpy(_delt_vtx, delt_vtx, sizeof(PDM_g_num_t) * dn_elt * elt_vtx_n);

    int _id_section = PDM_DMesh_nodal_elmts_section_ho_add((*dmn)->volumic,
                                                           elt_type,
                                                           order,
                                                           ho_ordering);
    PDM_DMesh_nodal_elmts_section_std_set((*dmn)->volumic,
                                          _id_section,
                                          dn_elt,
                                          _delt_vtx,
                                          PDM_OWNERSHIP_KEEP);

    gn_elt_vol += distrib_elt[n_rank];
  }
  (*dmn)->volumic->n_g_elmts = gn_elt_vol;


  // Surface
  geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
  n_section   = PDM_DMesh_nodal_n_section_get  (_dmn, geom_kind);
  sections_id = PDM_DMesh_nodal_sections_id_get(_dmn, geom_kind);

  PDM_g_num_t gn_elt_surf = 0;
  for (int i_section = 0; i_section < n_section; i_section++) {

    int id_section = sections_id[i_section];
    int _order;
    const char *ho_ordering = NULL;
    const PDM_g_num_t    *distrib_elt = PDM_DMesh_nodal_distrib_section_get(_dmn, geom_kind, id_section);
    int                   dn_elt      = PDM_DMesh_nodal_section_n_elt_get  (_dmn, geom_kind, id_section);
    PDM_g_num_t          *delt_vtx    = PDM_DMesh_nodal_section_std_ho_get (_dmn, geom_kind, id_section, &_order, &ho_ordering);
    PDM_Mesh_nodal_elt_t  elt_type    = PDM_DMesh_nodal_section_type_get   (_dmn, geom_kind, id_section);

    int elt_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);
    PDM_g_num_t *_delt_vtx = malloc(sizeof(PDM_g_num_t) * dn_elt * elt_vtx_n);
    memcpy(_delt_vtx, delt_vtx, sizeof(PDM_g_num_t) * dn_elt * elt_vtx_n);

    int _id_section = PDM_DMesh_nodal_elmts_section_ho_add((*dmn)->surfacic,
                                                           elt_type,
                                                           order,
                                                           ho_ordering);
    PDM_DMesh_nodal_elmts_section_std_set((*dmn)->surfacic,
                                          _id_section,
                                          dn_elt,
                                          _delt_vtx,
                                          PDM_OWNERSHIP_KEEP);

    gn_elt_surf += distrib_elt[n_rank];
  }
  (*dmn)->surfacic->n_g_elmts = gn_elt_surf;


  // Groups
  PDM_g_num_t *distrib_face = PDM_compute_uniform_entity_distribution(comm,
                                                                      gn_elt_surf);
  int n_group = 1;
  int *dgroup_elt_idx = malloc(sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;
  dgroup_elt_idx[n_group] = (int) (distrib_face[i_rank+1] - distrib_face[i_rank]);

  PDM_g_num_t *dgroup_elt = malloc(sizeof(PDM_g_num_t) * dgroup_elt_idx[n_group]);
  for (int i = 0; i < dgroup_elt_idx[n_group]; i++) {
    dgroup_elt[i] = distrib_face[i_rank] + i + 1;
  }

  free(distrib_face);

  PDM_DMesh_nodal_elmts_group_set((*dmn)->surfacic,
                                  n_group,
                                  dgroup_elt_idx,
                                  dgroup_elt,
                                  PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_generate_distribution(*dmn);


  PDM_dcube_nodal_gen_free(dcube);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
