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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_sphere_surf_gen.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_config.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_array.h"
#include "pdm_distrib.h"
#include "pdm_printf.h"
#include "pdm_error.h"

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

/*=============================================================================
 * Public function definitions
 *============================================================================*/

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
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  assert(nu > 1);
  assert(nv > 1);

  /*
   *  Vertices
   */
  PDM_g_num_t gn_vtx = (nu - 1) * nv + 2;

  PDM_g_num_t *_distrib_vtx = PDM_compute_uniform_entity_distribution(comm, gn_vtx);

  const double step_u = 2*PDM_PI / (double) (nu - 2);
  const double step_v =   PDM_PI / (double) (nv + 1);

  int dn_vtx = (int) (_distrib_vtx[i_rank+1] - _distrib_vtx[i_rank]);
  double *_dvtx_coord = (double *) malloc(sizeof(double) * dn_vtx * 3);

  for (int i = 0; i < dn_vtx; i++) {

    PDM_g_num_t g = _distrib_vtx[i_rank] + i;

    if (g == gn_vtx-1) {
      // north pole
      _dvtx_coord[3*i  ] = x_center;
      _dvtx_coord[3*i+1] = y_center;
      _dvtx_coord[3*i+2] = z_center + radius;

    } else if (g == gn_vtx-2) {
      // south pole
      _dvtx_coord[3*i  ] = x_center;
      _dvtx_coord[3*i+1] = y_center;
      _dvtx_coord[3*i+2] = z_center - radius;

    } else {

      PDM_g_num_t jj = g / (nu - 1);
      PDM_g_num_t ii = g % (nu - 1);

      double v = -0.5*PDM_PI + (jj+1)*step_v;
      double u = ii*step_u;
      double c = cos(v);

      _dvtx_coord[3*i  ] = x_center + radius * c * cos(u);
      _dvtx_coord[3*i+1] = y_center + radius * c * sin(u);
      _dvtx_coord[3*i+2] = z_center + radius * sin(v);

    }

  }

  /*
   *  Faces
   */
  PDM_g_num_t gn_face = 2*((nu - 1)*(nv - 1) + nu - 1);

  PDM_g_num_t *_distrib_face = PDM_compute_uniform_entity_distribution(comm, gn_face);

  int           dn_face       = (int) (_distrib_face[i_rank+1] - _distrib_face[i_rank]);
  int         *_dface_vtx_idx = (int *)         malloc(sizeof(int)         * (dn_face + 1));
  PDM_g_num_t *_dface_vtx     = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dn_face * 3);

  _dface_vtx_idx[0] = 0;
  for (int i = 0; i < dn_face; i++) {

    PDM_g_num_t g = _distrib_face[i_rank] + i;

    _dface_vtx_idx[i+1] = _dface_vtx_idx[i] + 3;
    PDM_g_num_t *_fv = _dface_vtx + _dface_vtx_idx[i];


    if (g >= 2*(nu - 1)*(nv - 1) + nu - 1) {

      // north pole cap
      PDM_g_num_t ii = g - (2*(nu - 1)*(nv - 1) + nu - 1);

      _fv[0] = 1 + ii            + (nu-1)*(nv-1);
      _fv[1] = 1 + (ii+1)%(nu-1) + (nu-1)*(nv-1);
      _fv[2] = gn_vtx;
    }

    else if (g >= 2*(nu - 1)*(nv - 1)) {

      // south pole cap
      PDM_g_num_t ii = g - 2*(nu - 1)*(nv - 1);

      _fv[0] = 1 + ii;
      _fv[1] = gn_vtx - 1;
      _fv[2] = 1 + (ii+1)%(nu-1);
    }

    else {

      if (g >= (nu - 1)*(nv - 1)) {
        g -= (nu - 1)*(nv - 1);

        PDM_g_num_t jj = g / (nu - 1);
        PDM_g_num_t ii = g % (nu - 1);

        _fv[0] = 1 + (ii+1)%(nu-1) + (nu-1)*jj;
        _fv[1] = 1 + (ii+1)%(nu-1) + (nu-1)*(jj+1);
        _fv[2] = 1 + ii            + (nu-1)*(jj+1);
      }

      else {

        PDM_g_num_t jj = g / (nu - 1);
        PDM_g_num_t ii = g % (nu - 1);

        _fv[0] = 1 + ii +            (nu-1)*jj;
        _fv[1] = 1 + (ii+1)%(nu-1) + (nu-1)*jj;
        _fv[2] = 1 + ii            + (nu-1)*(jj+1);

      }
    }
  }

  *distrib_vtx  = _distrib_vtx;
  *distrib_face = _distrib_face;

  *dvtx_coord    = _dvtx_coord;
  *dface_vtx_idx = _dface_vtx_idx;
  *dface_vtx     = _dface_vtx;


}



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
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  double      *dvtx_coord    = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *distrib_vtx   = NULL;
  PDM_g_num_t *distrib_face  = NULL;
  PDM_sphere_surf_gen(comm,
                      nu,
                      nv,
                      x_center,
                      y_center,
                      z_center,
                      radius,
                      &dvtx_coord,
                      &dface_vtx_idx,
                      &dface_vtx,
                      &distrib_vtx,
                      &distrib_face);

  /*
   *  Create dmesh nodal
   */
  PDM_g_num_t gn_vtx  = distrib_vtx [n_rank];
  PDM_g_num_t gn_face = distrib_face[n_rank];
  int dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank];
  int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];

  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  2,
                                                  gn_vtx,
                                                  0,
                                                  gn_face,
                                                  0);

  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  dmn->surfacic->n_g_elmts = gn_face;
  int id_section = PDM_DMesh_nodal_elmts_section_add(dmn->surfacic,
                                                     PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmn->surfacic,
                                        id_section,
                                        dn_face,
                                        dface_vtx,
                                        PDM_OWNERSHIP_KEEP);

  int n_group = 0;
  int *dgroup_elt_idx = (int *) malloc(sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;
  PDM_g_num_t *dgroup_elt = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dgroup_elt_idx[n_group]);
  PDM_DMesh_nodal_elmts_group_set(dmn->surfacic,
                                  n_group,
                                  dgroup_elt_idx,
                                  dgroup_elt,
                                  PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_generate_distribution(dmn);

  free(distrib_vtx);
  free(distrib_face);
  free(dface_vtx_idx);

  *_dmn = dmn;

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
