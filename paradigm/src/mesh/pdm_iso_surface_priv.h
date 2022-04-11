#ifndef __PDM_ISO_SURFACE_PRIV_H__
#define __PDM_ISO_SURFACE_PRIV_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define N_TIMER_ISO_SURFACE 5

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {

  SCAN_EDGES,
  EXTRACTION,
  BUILD_ISOSURF,
  TOTAL,
  WRITE_ISOSURF

} _iso_surface_timer_step_t;


/**
 * \struct _pdm_iso_surface
 * \brief  Define a partition mesh. Arrays are shared
 *
 */

struct _pdm_iso_surface_t
{
  int                    dim;
  int                    n_part;
  PDM_iso_surface_kind_t iso_kind;
  PDM_ownership_t        ownership;
  PDM_MPI_Comm           comm;

  int                  is_dist; // Ins are distributed

  /* Plane equation */
  double plane_equation[4];

  void (*eval_field_and_gradient) (const double, const double, const double,
                                   double *,
                                   double *, double *, double *);


  /* Distributed view */
  PDM_g_num_t         *dcell_face;
  int                 *dcell_face_idx;
  PDM_g_num_t         *dface_edge;
  int                 *dface_edge_idx;
  PDM_g_num_t         *dedge_vtx;

  PDM_g_num_t         *distrib_cell;
  PDM_g_num_t         *distrib_face;
  PDM_g_num_t         *distrib_edge;
  PDM_g_num_t         *distrib_vtx;

  /* Shortcut possible but we need to compute edge */
  PDM_g_num_t         *dface_vtx;
  int                 *dface_vtx_idx;

  double              *dvtx_coord;

  double              *dfield;
  double              *dgradient_field;

  /* Partitioned view - To do with extract for selected gnum + part_to_part */
  int                 *n_cell;
  int                 *n_face;
  int                 *n_edge;
  int                 *n_vtx;
  int                **pcell_face;
  int                **pcell_face_idx;
  int                **pface_edge;
  int                **pface_edge_idx;
  int                **pedge_vtx;

  PDM_g_num_t        **cell_ln_to_gn;
  PDM_g_num_t        **face_ln_to_gn;
  PDM_g_num_t        **edge_ln_to_gn;
  PDM_g_num_t        **vtx_ln_to_gn;

  /* Shortcut possible but we need to compute edge */
  int                **pface_vtx;
  int                **pface_vtx_idx;

  double             **pvtx_coord;
  double             **pfield;
  double             **pgradient_field;

  /* Result */
  int                  isosurf_n_vtx;
  int                  isosurf_n_edge; // only in 2d
  int                  isosurf_n_face;
  int                 *isosurf_face_vtx_idx;
  int                 *isosurf_face_vtx;
  int                 *isosurf_edge_vtx_idx; // only in 2d
  int                 *isosurf_edge_vtx; // only in 2d
  double              *isosurf_vtx_coord;
  PDM_g_num_t         *isosurf_vtx_ln_to_gn;
  PDM_g_num_t         *isosurf_edge_ln_to_gn; // only in 2d
  PDM_g_num_t         *isosurf_face_ln_to_gn;
  PDM_g_num_t         *isosurf_edge_parent_ln_to_gn; // parent face, only in 2d
  PDM_g_num_t         *isosurf_face_parent_ln_to_gn; // parent cell

  double times_elapsed[N_TIMER_ISO_SURFACE];

  int debug;
};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_H__ */
